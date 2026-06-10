/*---------------------------------------------------------------------------*\
Simplified Actuator Disk Model (ADM) for Wind Turbines
Based on SOWFA-6 horizontalAxisWindTurbinesADM
Adapted for OpenFOAM v2512

Force calculation based on Cp/Ct curves instead of BEM theory.
\*---------------------------------------------------------------------------*/

#include "horizontalAxisWindTurbinesADM_Simple.H"
#include "fvc.H"

namespace Foam
{
namespace turbineModels
{

horizontalAxisWindTurbinesADM_Simple::horizontalAxisWindTurbinesADM_Simple
(
    const volVectorField& U
)
:
    runTime_(U.time()),
    mesh_(U.mesh()),
    U_(U),
    degRad(Foam::constant::mathematical::pi / 180.0),
    dt(runTime_.deltaT().value()),
    time(runTime_.timeName()),
    t(runTime_.value()),
    pastFirstTimeStep(false),

    gradU
    (
        IOobject
        (
            "gradU",
            time,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor("gradU", dimVelocity/dimLength, tensor::zero)
    ),

    bodyForce
    (
        IOobject
        (
            "bodyForce",
            time,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("bodyForce", dimForce/dimVolume/dimDensity, vector::zero)
    )
{
    // -----------------------------------------------------------------------
    // Step 1: Read turbineArrayProperties
    // -----------------------------------------------------------------------
    IOdictionary turbineArrayProperties
    (
        IOobject
        (
            "turbineArrayProperties",
            runTime_.constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Extract turbine names (skip "globalProperties")
    {
        List<word> listTemp = turbineArrayProperties.toc();
        for (int i = 0; i < listTemp.size(); i++)
        {
            if (listTemp[i] != "globalProperties")
            {
                turbineName.append(listTemp[i]);
            }
        }
    }

    numTurbines = turbineName.size();

    outputControl = turbineArrayProperties.subDict("globalProperties")
                    .lookupOrDefault<word>("outputControl", "timeStep");
    outputInterval = turbineArrayProperties.subDict("globalProperties")
                     .lookupOrDefault<scalar>("outputInterval", 1);
    lastOutputTime = runTime_.startTime().value();
    outputIndex = 0;

    // Step 1a: Read basic parameters (defer nRadial/azimuthMaxDis/epsilon)
    forAll(turbineName, i)
    {
        const dictionary& turbDict = turbineArrayProperties.subDict(turbineName[i]);

        turbineType.append(turbDict.get<word>("turbineType"));
        baseLocation.append(turbDict.get<vector>("baseLocation"));
        pointDistType.append(turbDict.lookupOrDefault<word>("pointDistType", "uniform"));
        pointInterpType.append(turbDict.lookupOrDefault<word>("pointInterpType", "cellCenter"));
        inflowVelocityScalar.append(turbDict.lookupOrDefault<scalar>("inflowVelocityScalar", 1.0));
        nacYaw.append(turbDict.get<scalar>("NacYaw"));
        fluidDensity.append(turbDict.get<scalar>("fluidDensity"));

        // Temporarily store user-provided values (or -1 if not provided)
        nRadial.append(turbDict.lookupOrDefault<label>("nRadial", -1));
        azimuthMaxDis.append(turbDict.lookupOrDefault<scalar>("azimuthMaxDis", -1.0));
        epsilon.append(turbDict.lookupOrDefault<scalar>("epsilon", -1.0));
    }

    // -----------------------------------------------------------------------
    // Step 2: Catalog distinct turbine types
    // -----------------------------------------------------------------------
    numTurbinesDistinct = 1;
    {
        turbineTypeDistinct.append(turbineType[0]);
        forAll(turbineType, i)
        {
            bool flag = false;
            for (int j = 0; j < numTurbinesDistinct; j++)
            {
                if (turbineType[i] == turbineTypeDistinct[j])
                {
                    flag = true;
                }
            }
            if (!flag)
            {
                numTurbinesDistinct++;
                turbineTypeDistinct.append(turbineType[i]);
            }
        }
    }
    forAll(turbineType, i)
    {
        for (int j = 0; j < numTurbinesDistinct; j++)
        {
            if (turbineType[i] == turbineTypeDistinct[j])
            {
                turbineTypeID.append(j);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Step 3: Read turbineProperties for each distinct type
    // -----------------------------------------------------------------------
    for (int i = 0; i < numTurbinesDistinct; i++)
    {
        IOdictionary turbineProperties
        (
            IOobject
            (
                turbineTypeDistinct[i],
                runTime_.constant(), "turbineProperties",
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        TipRad.append(turbineProperties.get<scalar>("TipRad"));
        HubRad.append(turbineProperties.get<scalar>("HubRad"));
        TowerHt.append(turbineProperties.get<scalar>("TowerHt"));
        Twr2Shft.append(turbineProperties.lookupOrDefault<scalar>("Twr2Shft", 0.0));
        UndSling.append(turbineProperties.lookupOrDefault<scalar>("UndSling", 0.0));
        OverHang.append(turbineProperties.get<scalar>("OverHang"));
        ShftTilt.append(turbineProperties.lookupOrDefault<scalar>("ShftTilt", 0.0));

        // Read performance curves: support both PowerCtData and CpCtData formats
        List<scalar> wsTable, powerVals, ctVals;

        if (turbineProperties.found("PowerCtData"))
        {
            // Industry-standard format: (windSpeed[m/s], Power[kW], Ct)
            List<List<scalar>> PowerCtData(turbineProperties.lookup("PowerCtData"));

            scalar rho = 1.225;  // default density
            forAll(turbineType, t)
            {
                if (turbineType[t] == turbineTypeDistinct[i])
                {
                    rho = fluidDensity[t];
                    break;
                }
            }

            forAll(PowerCtData, j)
            {
                scalar V   = PowerCtData[j][0];   // wind speed [m/s]
                scalar P_W = PowerCtData[j][1] * 1000.0;  // kW → W
                scalar Ct  = PowerCtData[j][2];

                wsTable.append(V);
                ctVals.append(Ct);
                // Store power divided by density (for incompressible solver)
                powerVals.append(P_W / rho);
            }
        }
        else
        {
            // Legacy format: (windSpeed[m/s], Cp, Ct)
            List<List<scalar>> CpCtData(turbineProperties.lookup("CpCtData"));

            scalar A = Foam::constant::mathematical::pi * Foam::sqr(TipRad[i]);
            scalar rho = 1.225;
            forAll(turbineType, t)
            {
                if (turbineType[t] == turbineTypeDistinct[i])
                {
                    rho = fluidDensity[t];
                    break;
                }
            }

            forAll(CpCtData, j)
            {
                scalar V  = CpCtData[j][0];
                scalar Cp = CpCtData[j][1];
                scalar Ct = CpCtData[j][2];

                wsTable.append(V);
                ctVals.append(Ct);
                // Convert Cp to power: P = 0.5 * rho * A * V^3 * Cp
                scalar P_W = (V > SMALL)
                    ? 0.5 * rho * A * Foam::pow(V, 3) * Cp
                    : 0.0;
                powerVals.append(P_W / rho);
            }
        }

        windSpeedTable.append(wsTable);
        PowerTable.append(powerVals);
        CtTable.append(ctVals);
    }

    // -----------------------------------------------------------------------
    // Step 4: Unit conversions
    // -----------------------------------------------------------------------
    // Convert nacYaw from compass degrees to standard radians
    forAll(nacYaw, i)
    {
        nacYaw[i] = degRad * compassToStandard(nacYaw[i]);
    }
    // Convert ShftTilt from degrees to radians
    forAll(ShftTilt, i)
    {
        ShftTilt[i] = degRad * ShftTilt[i];
    }

    // -----------------------------------------------------------------------
    // Step 5: Compute key geometry (tower-shaft intersection, rotor apex)
    // -----------------------------------------------------------------------
    for (int i = 0; i < numTurbines; i++)
    {
        int j = turbineTypeID[i];

        // Tower-shaft intersection
        towerShaftIntersect.append(baseLocation[i]);
        towerShaftIntersect[i].z() += TowerHt[j] + Twr2Shft[j];

        // Rotor apex (accounting for overhang and shaft tilt, without nacYaw)
        // nacYaw rotation will be applied later in Step 8
        rotorApex.append(towerShaftIntersect[i]);
        scalar overhangTotal = OverHang[j] + UndSling[j];
        rotorApex[i].x() += overhangTotal * Foam::cos(ShftTilt[j]);
        rotorApex[i].z() += overhangTotal * Foam::sin(ShftTilt[j]);
    }

    // -----------------------------------------------------------------------
    // Step 5b: Compute automatic parameters if not provided by user
    // -----------------------------------------------------------------------
    // CRITICAL: In parallel runs, each processor has a different mesh partition.
    // Auto-calculated epsilon/azimuthMaxDis/nRadial must be synchronized across
    // all processors to ensure totDiskPointsArray is identical everywhere.
    for (int i = 0; i < numTurbines; i++)
    {
        int j = turbineTypeID[i];

        // If epsilon not provided, estimate from local mesh at rotor center
        if (epsilon[i] < 0.0)
        {
            label cellID = mesh_.findNearestCell(rotorApex[i]);
            scalar localMeshSize = GREAT;

            if (cellID != -1)
            {
                localMeshSize = Foam::cbrt(mesh_.V()[cellID]);
            }

            // In parallel: find the processor that owns the cell closest to rotorApex
            // and use its mesh size. Use min to find the processor with valid data.
            if (Pstream::parRun())
            {
                reduce(localMeshSize, minOp<scalar>());
            }

            if (localMeshSize < GREAT)
            {
                epsilon[i] = 2.5 * localMeshSize;
                if (Pstream::master())
                {
                    Info << "Turbine " << i << ": Auto epsilon = " << epsilon[i]
                         << " m (2.5 × mesh size)" << endl;
                }
            }
            else
            {
                FatalErrorInFunction
                    << "Cannot find mesh cell at turbine " << i
                    << " rotor apex " << rotorApex[i]
                    << exit(FatalError);
            }
        }

        // If azimuthMaxDis not provided, match epsilon
        if (azimuthMaxDis[i] < 0.0)
        {
            azimuthMaxDis[i] = epsilon[i];
            if (Pstream::master())
            {
                Info << "Turbine " << i << ": Auto azimuthMaxDis = "
                     << azimuthMaxDis[i] << " m" << endl;
            }
        }

        // If nRadial not provided, ensure ~epsilon spacing
        if (nRadial[i] < 0)
        {
            scalar bladeSpan = TipRad[j] - HubRad[j];
            nRadial[i] = max(10, label(bladeSpan / epsilon[i]));
            if (Pstream::master())
            {
                Info << "Turbine " << i << ": Auto nRadial = " << nRadial[i] << endl;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Step 6: Define sphere of influence for each turbine
    // -----------------------------------------------------------------------
    for (int i = 0; i < numTurbines; i++)
    {
        // Projection radius: where Gaussian decays to 0.1% of peak
        projectionRadius.append(epsilon[i] * Foam::sqrt(Foam::log(1.0/0.001)));

        int j = turbineTypeID[i];
        scalar sphereRadius = Foam::sqrt(
            Foam::sqr(OverHang[j] + UndSling[j])
            + Foam::sqr(TipRad[j])
        );
        sphereRadius += projectionRadius[i];

        DynamicList<label> sphereCellsI;
        forAll(U_.mesh().cells(), cellI)
        {
            if (mag(U_.mesh().C()[cellI] - towerShaftIntersect[i]) <= sphereRadius)
            {
                sphereCellsI.append(cellI);
            }
        }
        sphereCells.append(sphereCellsI);
        sphereCellsI.clear();

        if (sphereCells[i].size() > 0)
        {
            turbinesControlled.append(i);
        }
    }

    // -----------------------------------------------------------------------
    // Step 7: Generate actuator disk points
    // -----------------------------------------------------------------------
    totDiskPointsArray = 0;

    for (int i = 0; i < numTurbines; i++)
    {
        totDiskPoints.append(0);
        int j = turbineTypeID[i];

        // Unit vector along shaft (in yaw=0 reference frame, will be rotated in Step 8)
        // x-axis aligned (pointing downwind for zero yaw), accounting only for shaft tilt
        vector uvShaft_temp;
        uvShaft_temp.x() = Foam::cos(ShftTilt[j]);
        uvShaft_temp.y() = 0.0;
        uvShaft_temp.z() = Foam::sin(ShftTilt[j]);
        uvShaft.append(uvShaft_temp);

        // Unit vector along tower
        uvTower.append(towerShaftIntersect[i] - baseLocation[i]);
        uvTower[i] = uvTower[i] / mag(uvTower[i]);

        // Radial width of each actuator section
        dr.append(List<scalar>(nRadial[i], 0.0));
        if (pointDistType[i] == "uniform")
        {
            scalar actuatorRadialWidth = (TipRad[j] - HubRad[j]) / nRadial[i];
            for (int m = 0; m < nRadial[i]; m++)
            {
                dr[i][m] = actuatorRadialWidth;
            }
        }

        // Azimuthal discretization at each radial station
        bladeRadius.append(List<scalar>(nRadial[i], 0.0));
        solidity.append(List<scalar>(nRadial[i], 0.0));
        nAzimuth.append(List<int>(nRadial[i], 0));

        scalar dist = 0.0;
        for (int m = 0; m < nRadial[i]; m++)
        {
            dist += 0.5 * dr[i][m];
            bladeRadius[i][m] = HubRad[j] + dist;
            dist += 0.5 * dr[i][m];

            scalar circ = 2.0 * Foam::constant::mathematical::pi * bladeRadius[i][m];
            nAzimuth[i][m] = max(1, int(circ / azimuthMaxDis[i]));

            // Solidity: 1 per azimuthal point (uniform disk, no blade count)
            solidity[i][m] = 1.0 / scalar(nAzimuth[i][m]);
        }

        // Allocate point arrays
        bladePoints.append(List<List<vector>>(nRadial[i]));
        elementAzimuth.append(List<List<scalar>>(nRadial[i]));
        bladeForce.append(List<List<vector>>(nRadial[i]));
        windVectors.append(List<List<vector>>(nRadial[i]));
        minDisCellID.append(List<List<label>>(nRadial[i]));

        // Generate point coordinates
        // ADM: disk is perpendicular to shaft axis
        scalar beta = -ShftTilt[j];
        vector root = rotorApex[i];
        root.x() += HubRad[j] * Foam::sin(beta);
        root.z() += HubRad[j] * Foam::cos(beta);

        dist = 0.0;
        for (int m = 0; m < nRadial[i]; m++)
        {
            dist += 0.5 * dr[i][m];
            scalar dAzimuth = 2.0 * Foam::constant::mathematical::pi / nAzimuth[i][m];

            bladePoints[i][m].setSize(nAzimuth[i][m], vector::zero);
            elementAzimuth[i][m].setSize(nAzimuth[i][m], 0.0);
            bladeForce[i][m].setSize(nAzimuth[i][m], vector::zero);
            windVectors[i][m].setSize(nAzimuth[i][m], vector::zero);
            minDisCellID[i][m].setSize(nAzimuth[i][m], -1);

            for (int k = 0; k < nAzimuth[i][m]; k++)
            {
                totDiskPointsArray++;
                totDiskPoints[i]++;
                elementAzimuth[i][m][k] = k * dAzimuth;

                bladePoints[i][m][k].x() = root.x() + dist * Foam::sin(beta);
                bladePoints[i][m][k].y() = root.y();
                bladePoints[i][m][k].z() = root.z() + dist * Foam::cos(beta);

                if (k > 0)
                {
                    bladePoints[i][m][k] = rotatePoint(
                        bladePoints[i][m][k],
                        rotorApex[i],
                        uvShaft[i],
                        elementAzimuth[i][m][k]
                    );
                }
            }
            dist += 0.5 * dr[i][m];
        }

        // Initialize output scalars
        inflowVelocity.append(0.0);
        V_inf_.append(0.0);
        thrust.append(0.0);
        powerRotor.append(0.0);
    }

    // -----------------------------------------------------------------------
    // Step 8: Yaw nacelle to initial position
    // -----------------------------------------------------------------------
    // Step 7 generated all geometry in the yaw=0 frame. Here we rotate
    // uvShaft, rotorApex and bladePoints together by nacYaw around the tower
    // axis. Rotating a direction vector uses rotationPoint = zero.
    forAll(uvTower, i)
    {
        scalar cosYaw = Foam::cos(nacYaw[i]);
        scalar sinYaw = Foam::sin(nacYaw[i]);

        rotorApex[i] = rotatePointCached(rotorApex[i], towerShaftIntersect[i], uvTower[i], cosYaw, sinYaw);
        uvShaft[i]   = rotatePointCached(uvShaft[i],   vector::zero,           uvTower[i], cosYaw, sinYaw);

        forAll(bladePoints[i], j)
        {
            forAll(bladePoints[i][j], k)
            {
                bladePoints[i][j][k] = rotatePointCached(
                    bladePoints[i][j][k],
                    towerShaftIntersect[i],
                    uvTower[i],
                    cosYaw,
                    sinYaw
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Step 9: Find control processors and do initial calculation
    // -----------------------------------------------------------------------
    findControlProcNo();

    // 预分配 windVectorsLocal_ 内存（避免每个时间步重新分配）
    windVectorsLocal_.setSize(totDiskPointsArray, vector::zero);

    // 预计算 Gaussian 核常量和插值类型标志（避免每个时间步重新计算）
    cachedEps3Pi15_.setSize(numTurbines);
    cachedInvEps2_.setSize(numTurbines);
    cachedProjRadSq_.setSize(numTurbines);
    useLinearInterp_.setSize(numTurbines);
    for (int i = 0; i < numTurbines; i++)
    {
        cachedEps3Pi15_[i] = Foam::pow(epsilon[i], 3)
                           * Foam::pow(Foam::constant::mathematical::pi, 1.5);
        cachedInvEps2_[i]  = 1.0 / Foam::sqr(epsilon[i]);
        cachedProjRadSq_[i] = Foam::sqr(projectionRadius[i]);
        useLinearInterp_[i] = (pointInterpType[i] == "linear");
    }

    computeWindVectors();
    computeAverageInflow();
    computeBladeForce();
    computeBodyForce();

    openOutputFiles();
    printOutputFiles();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

horizontalAxisWindTurbinesADM_Simple::~horizontalAxisWindTurbinesADM_Simple()
{
    if (Pstream::master())
    {
        delete thrustFile_;
        delete powerRotorFile_;
        delete inflowVelocityFile_;
        delete CtFile_;
    }
}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

void horizontalAxisWindTurbinesADM_Simple::update()
{
    dt   = runTime_.deltaT().value();
    time = runTime_.timeName();
    t    = runTime_.value();

    computeWindVectors();
    computeAverageInflow();
    computeBladeForce();
    computeBodyForce();

    outputIndex++;

    if (outputControl == "timeStep")
    {
        if (outputIndex >= outputInterval)
        {
            outputIndex = 0;
            printOutputFiles();
        }
    }
    else if (outputControl == "runTime")
    {
        if ((runTime_.value() - lastOutputTime) >= outputInterval)
        {
            lastOutputTime += outputInterval;
            printOutputFiles();
        }
    }

    pastFirstTimeStep = true;
}


void horizontalAxisWindTurbinesADM_Simple::findControlProcNo()
{
    // Phase 1: Find local nearest cell and distance for each blade point
    // Initialize all distances to GREAT (ensures consistent data across all processors)
    List<scalar> flatLocalDist(totDiskPointsArray, GREAT);

    // Only fill distances for turbines controlled by this processor
    forAll(turbinesControlled, p)
    {
        int i = turbinesControlled[p];

        // Calculate starting index for this turbine in the flat array
        int startIter = 0;
        for (int n = 0; n < i; n++)
        {
            startIter += totDiskPoints[n];
        }
        int iter = startIter;

        forAll(bladePoints[i], j)
        {
            forAll(bladePoints[i][j], k)
            {
                scalar minDis = GREAT;
                label  minCell = -1;

                forAll(sphereCells[i], m)
                {
                    scalar dis = mag(
                        mesh_.C()[sphereCells[i][m]]
                        - bladePoints[i][j][k]
                    );
                    if (dis < minDis)
                    {
                        minDis  = dis;
                        minCell = sphereCells[i][m];
                    }
                }
                minDisCellID[i][j][k] = minCell;
                flatLocalDist[iter] = minDis;
                iter++;
            }
        }
    }

    // Phase 2: Communicate global minimum distance across all processors
    // listCombineReduce ensures result is consistent on all ranks (no broadcast needed)
    List<scalar> flatGlobalDist = flatLocalDist;
    if (Pstream::parRun())
    {
        Pstream::listCombineReduce(flatGlobalDist, minEqOp<scalar>());
    }
    // flatGlobalDist[idx] now holds the global minimum distance for blade point idx

    // Phase 3: Only the processor owning the global minimum distance keeps valid minDisCellID
    // Other processors clear to -1 to avoid duplicate velocity contributions in computeWindVectors
    // CRITICAL: Must iterate over ALL turbines (not just turbinesControlled) to maintain
    // consistent indexing into flatLocalDist/flatGlobalDist across all processors
    int iter = 0;
    for (int i = 0; i < numTurbines; i++)
    {
        forAll(bladePoints[i], j)
        {
            forAll(bladePoints[i][j], k)
            {
                // If this processor's distance is NOT the global minimum, clear the cell ID
                if (flatLocalDist[iter] > flatGlobalDist[iter] + SMALL)
                {
                    minDisCellID[i][j][k] = -1;
                }
                iter++;
            }
        }
    }
}


void horizontalAxisWindTurbinesADM_Simple::computeWindVectors()
{
    // 使用预分配的成员变量，清零而非重新分配
    windVectorsLocal_ = vector::zero;

    bool anyLinearInterp = false;
    forAll(useLinearInterp_, ii)
    {
        if (useLinearInterp_[ii]) { anyLinearInterp = true; break; }
    }
    if (anyLinearInterp)
    {
        gradU = fvc::grad(U_);
    }

    forAll(turbinesControlled, p)
    {
        int i = turbinesControlled[p];
        int iter = 0;
        for (int n = 0; n < i; n++)
        {
            iter += totDiskPoints[n];
        }

        forAll(bladePoints[i], j)
        {
            forAll(bladePoints[i][j], k)
            {
                if (minDisCellID[i][j][k] != -1)
                {
                    windVectorsLocal_[iter] = U_[minDisCellID[i][j][k]];

                    if (useLinearInterp_[i])
                    {
                        vector dx = bladePoints[i][j][k]
                                  - mesh_.C()[minDisCellID[i][j][k]];
                        windVectorsLocal_[iter] += dx & gradU[minDisCellID[i][j][k]];
                    }
                }
                iter++;
            }
        }
    }

    // Combine across all processors (single collective call instead of per-element reduce)
    // listCombineReduce ensures result is consistent on all ranks (no broadcast needed)
    if (Pstream::parRun())
    {
        Pstream::listCombineReduce(windVectorsLocal_, plusEqOp<vector>());
    }

    int iter = 0;
    forAll(windVectors, i)
    {
        forAll(windVectors[i], j)
        {
            forAll(windVectors[i][j], k)
            {
                windVectors[i][j][k] = windVectorsLocal_[iter] * inflowVelocityScalar[i];
                iter++;
            }
        }
    }
}


void horizontalAxisWindTurbinesADM_Simple::computeAverageInflow()
{
    forAll(windVectors, i)
    {
        scalar vSum  = 0.0;
        label  nPts  = 0;

        forAll(windVectors[i], j)
        {
            forAll(windVectors[i][j], k)
            {
                // Axial component: uvShaft points downwind, wind also downwind → positive
                vSum += (windVectors[i][j][k] & uvShaft[i]);
                nPts++;
            }
        }

        inflowVelocity[i] = (nPts > 0) ? (vSum / nPts) : 0.0;
    }
}


void horizontalAxisWindTurbinesADM_Simple::computeBladeForce()
{
    forAll(bladeForce, i)
    {
        int n = turbineTypeID[i];

        // Disk-averaged wind speed measured in CFD (reduced by turbine induction)
        scalar V_disk = max(inflowVelocity[i], SMALL);

        // Recover V_inf (local undisturbed wind speed) from V_disk via momentum theory:
        //   V_disk = V_inf * (1 - a)
        //   a = 0.5 * (1 - sqrt(1 - Ct(V_inf)))
        //
        // Ct table is defined based on V_inf (IEC 61400-12 standard).
        // This fixed-point iteration is valid for both flat and complex terrain:
        // it removes only the turbine's own induction effect, preserving
        // terrain-induced flow variations already captured in V_disk.
        //
        // Convergence: f'(V_inf) < 0 for typical Ct curves (Ct decreases with V),
        // so iteration oscillates; under-relaxation (omega=0.5) damps it.
        scalar V_inf = V_disk;  // initial guess
        const scalar tol   = 1e-4;
        const int maxIter  = 30;
        const scalar omega = 0.5;  // under-relaxation factor

        for (int iter = 0; iter < maxIter; iter++)
        {
            // Ct must be in [0, 0.99]: Ct >= 1 makes sqrt(1-Ct) imaginary
            scalar Ct_k = interpolate(V_inf, windSpeedTable[n], CtTable[n]);
            Ct_k = max(0.0, min(Ct_k, 0.99));

            scalar a_k = 0.5 * (1.0 - Foam::sqrt(1.0 - Ct_k));

            // Raw update from momentum theory
            scalar V_inf_raw = V_disk / max(1.0 - a_k, 0.01);

            // Under-relaxation to damp oscillations
            scalar V_inf_new = (1.0 - omega) * V_inf + omega * V_inf_raw;

            if (mag(V_inf_new - V_inf) < tol)
            {
                V_inf = V_inf_new;
                break;
            }
            V_inf = V_inf_new;
        }

        // Store V_inf for output file consistency
        V_inf_[i] = V_inf;

        // Final lookup at converged V_inf
        scalar Ct    = interpolate(V_inf, windSpeedTable[n], CtTable[n]);
        scalar Power = interpolate(V_inf, windSpeedTable[n], PowerTable[n]);
        Ct = max(0.0, min(Ct, 0.99));

        // Rotor swept area
        scalar A = Foam::constant::mathematical::pi * Foam::sqr(TipRad[n]);

        // Thrust and power based on V_inf (consistent with table definition)
        thrust[i]     = 0.5 * Ct * Foam::sqr(V_inf) * A;  // N / rho
        powerRotor[i] = Power;                              // W / rho

        // Distribute thrust uniformly over all actuator points.
        // bladeForce is stored as force-per-point; solidity is NOT applied here
        // (the solidity multiplication in the old computeBodyForce loop cancelled
        // out exactly, so both were redundant).
        scalar totalPts = scalar(totDiskPoints[i]);
        vector fPerPoint = (totalPts > 0)
            ? (-thrust[i] / totalPts) * uvShaft[i]
            : vector::zero;

        forAll(bladeForce[i], j)
        {
            forAll(bladeForce[i][j], k)
            {
                bladeForce[i][j][k] = fPerPoint;
            }
        }

        // Induction factor for diagnostic output
        scalar a_final = 0.5 * (1.0 - Foam::sqrt(1.0 - Ct));

        Info << "Turbine " << i
             << ": V_disk = " << V_disk << " m/s"
             << "  V_inf = "  << V_inf  << " m/s"
             << "  a = "      << a_final
             << "  Ct = "     << Ct
             << "  Thrust = " << thrust[i] * fluidDensity[i] << " N"
             << "  Power = "  << powerRotor[i] * fluidDensity[i] / 1e6 << " MW"
             << endl;
    }
}


void horizontalAxisWindTurbinesADM_Simple::computeBodyForce()
{
    bodyForce *= 0.0;

    scalar thrustSum          = 0.0;
    scalar thrustBodyForceSum = 0.0;

    forAll(bladeForce, i)
    {
        if (sphereCells[i].size() > 0)
        {
            const scalar invEps2    = cachedInvEps2_[i];
            const scalar eps3_pi15  = cachedEps3Pi15_[i];
            const scalar projRadSq  = cachedProjRadSq_[i];

            // 交换循环顺序：以球内单元为外层，执行器点为内层
            // 优势：每个网格单元坐标只读取一次（原来是读 N_bladePoints 次），
            //       bodyForce[cell] 只写入一次（原来是读-改-写 N_bladePoints 次）
            forAll(sphereCells[i], m)
            {
                const label  cellID = sphereCells[i][m];
                const vector cellC  = mesh_.C()[cellID];
                const scalar cellV  = mesh_.V()[cellID];

                vector sumForce(vector::zero);
                scalar sumThrust(0.0);

                forAll(bladeForce[i], j)
                {
                    forAll(bladeForce[i][j], k)
                    {
                        scalar disSq = magSqr(cellC - bladePoints[i][j][k]);

                        if (disSq <= projRadSq)
                        {
                            scalar gaussKernel =
                                Foam::exp(-disSq * invEps2) / eps3_pi15;

                            // bladeForce already stores force-per-point without solidity scaling
                            vector contrib = bladeForce[i][j][k] * gaussKernel;

                            sumForce  += contrib;
                            sumThrust += (-contrib * cellV) & uvShaft[i];
                        }
                    }
                }

                bodyForce[cellID]   += sumForce;
                thrustBodyForceSum  += sumThrust;
            }
        }
        thrustSum += thrust[i];
    }

    reduce(thrustBodyForceSum, sumOp<scalar>());

    if (mag(thrustSum) > SMALL)
    {
        Info << "Thrust from Body Force = " << thrustBodyForceSum
             << "  Thrust from Act. Disk = " << thrustSum
             << "  Ratio = " << thrustBodyForceSum / thrustSum
             << endl;
    }
}


vector horizontalAxisWindTurbinesADM_Simple::rotatePoint
(
    vector point,
    vector rotationPoint,
    vector axis,
    scalar angle
)
{
    tensor RM;
    RM.xx() = Foam::sqr(axis.x()) + (1.0 - Foam::sqr(axis.x())) * Foam::cos(angle);
    RM.xy() = axis.x()*axis.y()*(1.0 - Foam::cos(angle)) - axis.z()*Foam::sin(angle);
    RM.xz() = axis.x()*axis.z()*(1.0 - Foam::cos(angle)) + axis.y()*Foam::sin(angle);
    RM.yx() = axis.x()*axis.y()*(1.0 - Foam::cos(angle)) + axis.z()*Foam::sin(angle);
    RM.yy() = Foam::sqr(axis.y()) + (1.0 - Foam::sqr(axis.y())) * Foam::cos(angle);
    RM.yz() = axis.y()*axis.z()*(1.0 - Foam::cos(angle)) - axis.x()*Foam::sin(angle);
    RM.zx() = axis.x()*axis.z()*(1.0 - Foam::cos(angle)) - axis.y()*Foam::sin(angle);
    RM.zy() = axis.y()*axis.z()*(1.0 - Foam::cos(angle)) + axis.x()*Foam::sin(angle);
    RM.zz() = Foam::sqr(axis.z()) + (1.0 - Foam::sqr(axis.z())) * Foam::cos(angle);

    point = point - rotationPoint;
    point = RM & point;
    point = point + rotationPoint;
    return point;
}


vector horizontalAxisWindTurbinesADM_Simple::rotatePointCached
(
    vector point,
    vector rotationPoint,
    vector axis,
    scalar cosAngle,
    scalar sinAngle
)
{
    tensor RM;
    RM.xx() = Foam::sqr(axis.x()) + (1.0 - Foam::sqr(axis.x())) * cosAngle;
    RM.xy() = axis.x()*axis.y()*(1.0 - cosAngle) - axis.z()*sinAngle;
    RM.xz() = axis.x()*axis.z()*(1.0 - cosAngle) + axis.y()*sinAngle;
    RM.yx() = axis.x()*axis.y()*(1.0 - cosAngle) + axis.z()*sinAngle;
    RM.yy() = Foam::sqr(axis.y()) + (1.0 - Foam::sqr(axis.y())) * cosAngle;
    RM.yz() = axis.y()*axis.z()*(1.0 - cosAngle) - axis.x()*sinAngle;
    RM.zx() = axis.x()*axis.z()*(1.0 - cosAngle) - axis.y()*sinAngle;
    RM.zy() = axis.y()*axis.z()*(1.0 - cosAngle) + axis.x()*sinAngle;
    RM.zz() = Foam::sqr(axis.z()) + (1.0 - Foam::sqr(axis.z())) * cosAngle;

    point = point - rotationPoint;
    point = RM & point;
    point = point + rotationPoint;
    return point;
}


scalar horizontalAxisWindTurbinesADM_Simple::interpolate
(
    scalar xNew,
    List<scalar>& xOld,
    List<scalar>& yOld
)
{
    if (xOld.size() < 2) return (xOld.size() == 1 ? yOld[0] : 0.0);

    // Clamp to table range
    if (xNew <= xOld[0])          return yOld[0];
    if (xNew >= xOld.last())      return yOld.last();

    // 二分查找：O(log n) 替代线性搜索 O(n)
    label lo = 0;
    label hi = xOld.size() - 2;

    while (lo < hi)
    {
        label mid = (lo + hi) / 2;
        if (xNew > xOld[mid + 1])
        {
            lo = mid + 1;
        }
        else
        {
            hi = mid;
        }
    }

    // 线性插值
    scalar t = (xNew - xOld[lo]) / (xOld[lo+1] - xOld[lo]);
    return yOld[lo] + t * (yOld[lo+1] - yOld[lo]);
}


scalar horizontalAxisWindTurbinesADM_Simple::compassToStandard(scalar dir)
{
    dir += 180.0;
    if (dir >= 360.0) dir -= 360.0;
    dir = 90.0 - dir;
    if (dir < 0.0)    dir += 360.0;
    return dir;
}


void horizontalAxisWindTurbinesADM_Simple::openOutputFiles()
{
    if (Pstream::master())
    {
        fileName rootDir;
        if (Pstream::parRun())
            rootDir = runTime_.path() / "../turbineOutput";
        else
            rootDir = runTime_.path() / "turbineOutput";

        if (!isDir(rootDir))  mkDir(rootDir);
        if (!isDir(rootDir/time)) mkDir(rootDir/time);

        thrustFile_         = new OFstream(rootDir/time/"thrust");
        powerRotorFile_     = new OFstream(rootDir/time/"powerRotor");
        inflowVelocityFile_ = new OFstream(rootDir/time/"inflowVelocity");
        CtFile_             = new OFstream(rootDir/time/"Ct");

        *thrustFile_         << "#Turbine  Time(s)  dt(s)  thrust(N)" << endl;
        *powerRotorFile_     << "#Turbine  Time(s)  dt(s)  powerRotor(W)" << endl;
        *inflowVelocityFile_ << "#Turbine  Time(s)  dt(s)  inflowVelocity(m/s)" << endl;
        *CtFile_             << "#Turbine  Time(s)  dt(s)  Ct" << endl;
    }
}


void horizontalAxisWindTurbinesADM_Simple::printOutputFiles()
{
    if (Pstream::master())
    {
        forAll(thrust, i)
        {
            int n = turbineTypeID[i];
            // Ct must be looked up at V_inf (the table reference velocity per IEC 61400-12),
            // not at inflowVelocity (= V_disk, which is already reduced by induction).
            scalar Ct = interpolate(V_inf_[i], windSpeedTable[n], CtTable[n]);

            *thrustFile_         << i << " " << time << " " << dt << " "
                                 << thrust[i] * fluidDensity[i] << endl;
            *powerRotorFile_     << i << " " << time << " " << dt << " "
                                 << powerRotor[i] * fluidDensity[i] << endl;
            *inflowVelocityFile_ << i << " " << time << " " << dt << " "
                                 << inflowVelocity[i] << endl;
            *CtFile_             << i << " " << time << " " << dt << " "
                                 << Ct << endl;
        }
    }
}

} // End namespace turbineModels
} // End namespace Foam

// ************************************************************************* //
