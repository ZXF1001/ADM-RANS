/*---------------------------------------------------------------------------*\
Simplified Actuator Disk Model (ADM) for Wind Turbines
Based on SOWFA-6 horizontalAxisWindTurbinesADM
Adapted for OpenFOAM v2512

Simplifications:
- Removed: Controllers (torque, pitch, yaw)
- Removed: Detailed blade aerodynamics (BEM, airfoil Cl/Cd)
- Removed: Blade geometry (chord, twist)
- Added: Cp/Ct curve-based force calculation
- Retained: Actuator point distribution, Gaussian projection, CFD interface

\*---------------------------------------------------------------------------*/

#ifndef horizontalAxisWindTurbinesADM_Simple_H
#define horizontalAxisWindTurbinesADM_Simple_H

#include "HashPtrTable.H"
#include "IOdictionary.H"
#include "IFstream.H"
#include "OFstream.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbineModels
{

/*---------------------------------------------------------------------------*\
              Class horizontalAxisWindTurbinesADM_Simple declaration
\*---------------------------------------------------------------------------*/

class horizontalAxisWindTurbinesADM_Simple
{

private:
    // Private Data

        //- Constants
            //- Runtime pointer
            const Time& runTime_;

            //- Mesh pointer
            const fvMesh& mesh_;

            //- Velocity field pointer
            const volVectorField& U_;

            //- Degrees to radians conversion factor
            const scalar degRad;


        //- Current time step size
        scalar dt;

        //- Current simulation time
        word time;
        scalar t;

        //- Boolean that is true past first time step
        bool pastFirstTimeStep;

        //- Velocity gradient for linear interpolation
        volTensorField gradU;

        //- Body force field applied to fluid by turbine
        volVectorField bodyForce;

        //- Output control
        word outputControl;
        scalar outputInterval;
        scalar lastOutputTime;
        label outputIndex;


        //- Turbine Array Level Data
            //- List of turbine names
            List<word> turbineName;

            //- Number of turbines in array
            int numTurbines;

            //- List of turbine types
            DynamicList<word> turbineType;

            //- Base locations (m)
            DynamicList<vector> baseLocation;

            //- Number of radial actuator points
            DynamicList<int> nRadial;

            //- Maximum azimuthal element length (m)
            DynamicList<scalar> azimuthMaxDis;

            //- Point distribution type ("uniform")
            DynamicList<word> pointDistType;

            //- Point interpolation type ("cellCenter" or "linear")
            DynamicList<word> pointInterpType;

            //- Gaussian kernel width (m)
            DynamicList<scalar> epsilon;

            //- Force scaling factor
            DynamicList<scalar> forceScalar;

            //- Inflow velocity scaling factor
            DynamicList<scalar> inflowVelocityScalar;

            //- Projection radius (m)
            DynamicList<scalar> projectionRadius;

            //- Rotation direction ("cw" or "ccw") — unused in ADM, removed

            //- Fixed nacelle yaw angle (radians, standard convention)
            DynamicList<scalar> nacYaw;

            //- Fluid density (kg/m^3)
            DynamicList<scalar> fluidDensity;

            //- Number of distinct turbine types
            int numTurbinesDistinct;

            //- List of distinct turbine type names
            DynamicList<word> turbineTypeDistinct;

            //- Turbine type ID for each turbine
            DynamicList<label> turbineTypeID;


        //- Turbine Type Level Data (Geometry)
            //- Tip radius (m)
            DynamicList<scalar> TipRad;

            //- Hub radius (m)
            DynamicList<scalar> HubRad;

            //- Undersling distance (m)
            DynamicList<scalar> UndSling;

            //- Overhang distance (m)
            DynamicList<scalar> OverHang;

            //- Tower height (m)
            DynamicList<scalar> TowerHt;

            //- Tower to shaft vertical distance (m)
            DynamicList<scalar> Twr2Shft;

            //- Shaft tilt angle (radians)
            DynamicList<scalar> ShftTilt;


        //- Wind speed table for performance curves (m/s)
            DynamicList<List<scalar>> windSpeedTable;

            //- Thrust coefficient curve
            DynamicList<List<scalar>> CtTable;

            //- Power curve (W, divided by density)
            DynamicList<List<scalar>> PowerTable;


        //- Turbine Instance Level Data (Computed)
            //- Tower-shaft intersection point
            DynamicList<vector> towerShaftIntersect;

            //- Rotor apex location
            DynamicList<vector> rotorApex;

            //- Unit vector along shaft (pointing downwind)
            DynamicList<vector> uvShaft;

            //- Unit vector along tower (pointing up)
            DynamicList<vector> uvTower;

            //- Cells within sphere of influence for each turbine
            DynamicList<DynamicList<label>> sphereCells;

            //- List of turbines controlled by this processor
            DynamicList<label> turbinesControlled;


        //- Actuator Point Level Data
            //- Radial width of each actuator section (m)
            DynamicList<List<scalar>> dr;

            //- Blade radius at each radial station (m)
            DynamicList<List<scalar>> bladeRadius;

            //- Solidity factor at each radial station
            DynamicList<List<scalar>> solidity;

            //- Number of azimuthal points at each radial station
            DynamicList<List<int>> nAzimuth;

            //- Total number of disk points per turbine
            DynamicList<label> totDiskPoints;

            //- Total number of disk points in array
            label totDiskPointsArray;

            //- Actuator point locations [turbine][radial][azimuth]
            DynamicList<List<List<vector>>> bladePoints;

            //- Element azimuth angles (radians)
            DynamicList<List<List<scalar>>> elementAzimuth;

            //- Wind velocity vectors at actuator points (blade coordinates)
            DynamicList<List<List<vector>>> windVectors;

            //- Blade force vectors at actuator points (N/rho)
            DynamicList<List<List<vector>>> bladeForce;

            //- Closest cell ID for each actuator point
            DynamicList<List<List<label>>> minDisCellID;


        //- Turbine Output Data
            //- Disk-averaged inflow velocity (m/s)
            DynamicList<scalar> inflowVelocity;

            //- Total thrust force (N/rho)
            DynamicList<scalar> thrust;

            //- Rotor power (W/rho)
            DynamicList<scalar> powerRotor;


        //- Output file pointers
            OFstream* thrustFile_;
            OFstream* powerRotorFile_;
            OFstream* inflowVelocityFile_;
            OFstream* CtFile_;


    // Private Member Functions

        //- Rotate a point about an axis
        vector rotatePoint(vector point, vector rotationPoint, vector axis, scalar angle);

        //- Linear interpolation
        scalar interpolate(scalar xNew, List<scalar>& xOld, List<scalar>& yOld);

        //- Convert compass direction to standard mathematical convention
        scalar compassToStandard(scalar dir);

        //- Convert standard mathematical convention to compass direction
        scalar standardToCompass(scalar dir);

        //- Find which processor controls each actuator point
        void findControlProcNo();

        //- Compute wind velocity vectors at actuator points
        void computeWindVectors();

        //- Compute disk-averaged inflow velocity
        void computeAverageInflow();

        //- Compute blade forces using Cp/Ct curves
        void computeBladeForce();

        //- Project blade forces to body force field using Gaussian kernel
        void computeBodyForce();

        //- Open output files
        void openOutputFiles();

        //- Write output to files
        void printOutputFiles();


public:

    // Constructors

        //- Construct from components
        horizontalAxisWindTurbinesADM_Simple
        (
            const volVectorField& U
        );


    //- Destructor
    ~horizontalAxisWindTurbinesADM_Simple();


    // Member Functions

        //- Update the turbine state and forces
        void update();

        //- Return the body force field
        volVectorField& force()
        {
            return bodyForce;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbineModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
