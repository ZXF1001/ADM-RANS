# ADM_Simple 独立库

简化的执行器盘模型（Actuator Disk Model），从 SOWFA-6 提取并改进。

## 目录结构

```
ADM_Simple_Standalone/
├── Allwmake                    <- 一键编译（库 + solver）
├── Allwclean                   <- 一键清理
├── README.md
├── src/                        <- ADM 库源码
│   ├── horizontalAxisWindTurbinesADM_Simple.H
│   ├── horizontalAxisWindTurbinesADM_Simple.C
│   └── Make/
│       ├── files
│       └── options
└── applications/
    └── solvers/
        └── simpleFoamADM/      <- 稳态 RANS solver
            ├── simpleFoamADM.C
            ├── createFields.H
            ├── UEqn.H
            ├── pEqn.H
            └── Make/
```

## 特性

- 基于 Cp/Ct 曲线的力计算（无需 BEM 理论）
- 支持 PowerCtData（功率曲线）和 CpCtData（功率系数曲线）格式
- 诱导因子迭代修正（从盘面风速反推上游风速）
- 高斯核投影（正确归一化，积分=1.0）
- 移除预锥角（ADM 是平面盘，PreCone 无意义）
- 完全脱离 SOWFA，只依赖 OpenFOAM 标准库

## 编译

```bash
cd /home/zxf/OpenFOAM/zxf-v2512/run/ADM/ADM_Simple_Standalone
./Allwmake
```

编译后生成：
- 库：`$FOAM_USER_LIBBIN/libturbineModelsSimple.so`
- Solver：`$FOAM_USER_APPBIN/simpleFoamADM`

## 使用

### 1. 配置风机阵列

`constant/turbineArrayProperties`：
```
globalProperties
{
    outputControl   timeStep;
    outputInterval  1;
}

turbine0
{
    turbineType         "NREL5MW_Simple";
    baseLocation        (500 500 0);
    numBladePoints      50;
    pointDistType       "uniform";
    pointInterpType     "cellCenter";
    epsilon             10.0;
    tipRad              63.0;
    smearRadius         1.0;
    fluidDensity        1.225;
    azimuthMaxDis       10.0;
    nacYaw              0.0;
    inflowVelocityScalar 1.0;
}
```

### 2. 配置风机性能

`constant/turbineProperties/NREL5MW_Simple`：
```
TipRad    63.0;
HubRad    1.5;
TowerHt   90.0;
Twr2Shft  2.4;
OverHang  5.0;
UndSling  0.0;
ShftTilt  5.0;
PreCone   (2.5 2.5 2.5);  // ADM 中被忽略

PowerCtData
(
    (  3.0      0.0  0.000)    // (风速[m/s], 功率[kW], Ct)
    (  4.0    156.4  0.720)
    ...
    ( 11.4   5000.0  0.800)    // 额定风速，额定功率
    ...
    ( 25.0   5000.0  0.300)
);
```

### 3. 运行

```bash
simpleFoamADM
```

不需要在 `system/controlDict` 中手动加载库，solver 已经链接了 `libturbineModelsSimple.so`。

## 改进内容（相比原始 SOWFA-6）

1. 数据格式：支持行业标准的功率曲线格式（PowerCtData）
2. 来流修正：诱导因子迭代修正，从盘面风速反推上游风速
3. 几何修正：移除 ADM 中无意义的预锥角
4. 高斯核：确认原始公式正确（积分=1.0）
5. 独立编译：完全脱离 SOWFA，只依赖 OpenFOAM 标准库

## 依赖

- OpenFOAM v2512（或兼容版本）
- 标准库：finiteVolume, meshTools, sampling, turbulenceModels

## 参考

- 原始代码：SOWFA-6 (https://github.com/NatLabRockies/SOWFA-6)
