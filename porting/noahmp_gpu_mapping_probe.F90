! Actual model datatypes and runtime, without HRLDAS input or physical timesteps.
program noahmp_gpu_mapping_probe
  use Machine, only: kind_noahmp
  use NoahmpColumnFinalizeMod, only: NoahmpUpdateOutputAlbedo
  use NoahmpColumnResetMod, only: NoahmpResetStepTemporaries
  use NoahmpVarType, only: noahmp_type
  use NoahmpWorkspaceMod, only: NoahmpWorkspace_type,NoahmpAllocateWorkspacePool
  use NoahmpGpuDataMod, only: NoahmpGpuStart,NoahmpGpuUploadRefresh,NoahmpGpuBeforeRefresh, &
                            NoahmpGpuSyncHost,NoahmpGpuFinish,NoahmpGpuTraceStage
  use ConfigVarInitMod, only: ConfigVarInitDefault
  use ForcingVarInitMod, only: ForcingVarInitDefault
  use EnergyVarInitMod, only: EnergyVarInitDefault
  use WaterVarInitMod, only: WaterVarInitDefault
  use BiochemVarInitMod, only: BiochemVarInitDefault
  use openacc, only: acc_on_device,acc_device_nvidia
  implicit none
  integer, parameter :: ncolumns=205,capacity=64,nsoil=4,nsnow=3,nrad=2
  type(noahmp_type), allocatable :: columns(:)
  type(NoahmpWorkspace_type), allocatable :: workspaces(:)
  logical, allocatable :: active(:)
  integer :: p,device_flag,lane,errors
  allocate(columns(ncolumns),active(ncolumns))
  active=.true.
  ! Include skipped cells as well as a final partial workspace block.
  do p=1,ncolumns
    if (mod(p,7)==0) active(p)=.false.
    call ConfigVarInitDefault(columns(p))
    columns(p)%config%domain%NumSoilLayer=nsoil
    columns(p)%config%domain%NumSnowLayerMax=nsnow
    columns(p)%config%domain%NumSnowLayerNeg=0
    columns(p)%config%domain%NumSwRadBand=nrad
    columns(p)%config%domain%NumCropGrowStage=8
    columns(p)%config%domain%LandUseDataName='MODIFIED_IGBP_MODIS_NOAH'
    call ForcingVarInitDefault(columns(p))
    call EnergyVarInitDefault(columns(p))
    columns(p)%energy%state%AlbedoSfcOutput=0.125_kind_noahmp
    call WaterVarInitDefault(columns(p))
    call BiochemVarInitDefault(columns(p))
    columns(p)%water%state%SoilMoisture=real(p,kind_noahmp)
    columns(p)%energy%state%TemperatureSoilSnow=250.0_kind_noahmp
    columns(p)%forcing%TemperatureAirRefHeight=270.0_kind_noahmp
    columns(p)%initialized=.true.
  enddo
  ! Exercise sparse/empty groups and optional matrices without physical steps.
  deallocate(columns(1)%water%flux%TranspWatLossSoilMean)
  allocate(columns(1)%water%flux%TranspWatLossSoilMean(1:0))
  deallocate(columns(1)%biochem%param%LeafDeathTempCoeffCrop)
  allocate(columns(1)%energy%flux%FracRadSwAbsSnowDir(-2:0,2:3))
  allocate(columns(1)%energy%flux%FracRadSwAbsSnowDif(-2:0,2:3))
  columns(1)%energy%flux%FracRadSwAbsSnowDir=17.0_kind_noahmp
  columns(1)%energy%flux%FracRadSwAbsSnowDif=23.0_kind_noahmp
  call NoahmpAllocateWorkspacePool(workspaces,capacity,nsoil,nsnow,nrad)
  ! startup scratch coverage: every component starts with a distinct lane value.
  do lane=1,capacity
    workspaces(1)%GlacierPhaseChange%EnergyRes(:,lane)=real(1000*lane+1,kind_noahmp)
    workspaces(1)%GlacierPhaseChange%GlacierPhaseChg(:,lane)=real(1000*lane+2,kind_noahmp)
    workspaces(1)%GlacierPhaseChange%MassWatTotInit(:,lane)=real(1000*lane+3,kind_noahmp)
    workspaces(1)%GlacierPhaseChange%MassWatIceInit(:,lane)=real(1000*lane+4,kind_noahmp)
    workspaces(1)%GlacierPhaseChange%MassWatLiqInit(:,lane)=real(1000*lane+5,kind_noahmp)
    workspaces(1)%GlacierPhaseChange%MassWatIceTmp(:,lane)=real(1000*lane+6,kind_noahmp)
    workspaces(1)%GlacierPhaseChange%MassWatLiqTmp(:,lane)=real(1000*lane+7,kind_noahmp)
    workspaces(1)%GlacierPhaseChange%EnergyResLeft(:,lane)=real(1000*lane+8,kind_noahmp)
    workspaces(1)%GlacierTemperatureMain%MatRight(:,lane)=real(1000*lane+9,kind_noahmp)
    workspaces(1)%GlacierTemperatureMain%MatLeft1(:,lane)=real(1000*lane+10,kind_noahmp)
    workspaces(1)%GlacierTemperatureMain%MatLeft2(:,lane)=real(1000*lane+11,kind_noahmp)
    workspaces(1)%GlacierTemperatureMain%MatLeft3(:,lane)=real(1000*lane+12,kind_noahmp)
    workspaces(1)%GlacierTemperatureSolver%MatRightTmp(:,lane)=real(1000*lane+13,kind_noahmp)
    workspaces(1)%GlacierTemperatureSolver%MatLeft3Tmp(:,lane)=real(1000*lane+14,kind_noahmp)
    workspaces(1)%GlacierThermalDiffusion%DepthSnowSoilInv(:,lane)=real(1000*lane+15,kind_noahmp)
    workspaces(1)%GlacierThermalDiffusion%HeatCapacPerArea(:,lane)=real(1000*lane+16,kind_noahmp)
    workspaces(1)%GlacierThermalDiffusion%TempGradDepth(:,lane)=real(1000*lane+17,kind_noahmp)
    workspaces(1)%GlacierThermalDiffusion%EnergyExcess(:,lane)=real(1000*lane+18,kind_noahmp)
    workspaces(1)%GroundWaterTopModel%DepthSoilMid(:,lane)=real(1000*lane+19,kind_noahmp)
    workspaces(1)%GroundWaterTopModel%ThicknessSoil(:,lane)=real(1000*lane+20,kind_noahmp)
    workspaces(1)%GroundWaterTopModel%SoilLiqTmp(:,lane)=real(1000*lane+21,kind_noahmp)
    workspaces(1)%GroundWaterTopModel%SoilEffPorosity(:,lane)=real(1000*lane+22,kind_noahmp)
    workspaces(1)%GroundWaterTopModel%SoilWatConductTmp(:,lane)=real(1000*lane+23,kind_noahmp)
    workspaces(1)%GroundWaterTopModel%SoilMoisture(:,lane)=real(1000*lane+24,kind_noahmp)
    workspaces(1)%RunoffSurfaceFreeDrain%SoilWatMaxHold(:,lane)=real(1000*lane+25,kind_noahmp)
    workspaces(1)%SnowLayerDivide%SnowIceTmp(:,lane)=real(1000*lane+26,kind_noahmp)
    workspaces(1)%SnowLayerDivide%SnowLiqTmp(:,lane)=real(1000*lane+27,kind_noahmp)
    workspaces(1)%SnowLayerDivide%TemperatureSnowTmp(:,lane)=real(1000*lane+28,kind_noahmp)
    workspaces(1)%SnowLayerDivide%SnowThickTmp(:,lane)=real(1000*lane+29,kind_noahmp)
    workspaces(1)%SnowThermalProperty%SnowDensBulk(:,lane)=real(1000*lane+30,kind_noahmp)
    workspaces(1)%SoilMoistureSolver%MatRightTmp(:,lane)=real(1000*lane+31,kind_noahmp)
    workspaces(1)%SoilMoistureSolver%MatLeft3Tmp(:,lane)=real(1000*lane+32,kind_noahmp)
    workspaces(1)%SoilSnowTemperatureMain%MatRight(:,lane)=real(1000*lane+33,kind_noahmp)
    workspaces(1)%SoilSnowTemperatureMain%MatLeft1(:,lane)=real(1000*lane+34,kind_noahmp)
    workspaces(1)%SoilSnowTemperatureMain%MatLeft2(:,lane)=real(1000*lane+35,kind_noahmp)
    workspaces(1)%SoilSnowTemperatureMain%MatLeft3(:,lane)=real(1000*lane+36,kind_noahmp)
    workspaces(1)%SoilSnowTemperatureSolver%MatRightTmp(:,lane)=real(1000*lane+37,kind_noahmp)
    workspaces(1)%SoilSnowTemperatureSolver%MatLeft3Tmp(:,lane)=real(1000*lane+38,kind_noahmp)
    workspaces(1)%SoilSnowThermalDiffusion%DepthSnowSoilInv(:,lane)=real(1000*lane+39,kind_noahmp)
    workspaces(1)%SoilSnowThermalDiffusion%HeatCapacPerArea(:,lane)=real(1000*lane+40,kind_noahmp)
    workspaces(1)%SoilSnowThermalDiffusion%TempGradDepth(:,lane)=real(1000*lane+41,kind_noahmp)
    workspaces(1)%SoilSnowThermalDiffusion%EnergyExcess(:,lane)=real(1000*lane+42,kind_noahmp)
    workspaces(1)%SoilSnowWaterPhaseChange%EnergyRes(:,lane)=real(1000*lane+43,kind_noahmp)
    workspaces(1)%SoilSnowWaterPhaseChange%WaterPhaseChg(:,lane)=real(1000*lane+44,kind_noahmp)
    workspaces(1)%SoilSnowWaterPhaseChange%MassWatTotInit(:,lane)=real(1000*lane+45,kind_noahmp)
    workspaces(1)%SoilSnowWaterPhaseChange%MassWatIceInit(:,lane)=real(1000*lane+46,kind_noahmp)
    workspaces(1)%SoilSnowWaterPhaseChange%MassWatLiqInit(:,lane)=real(1000*lane+47,kind_noahmp)
    workspaces(1)%SoilSnowWaterPhaseChange%MassWatIceTmp(:,lane)=real(1000*lane+48,kind_noahmp)
    workspaces(1)%SoilSnowWaterPhaseChange%MassWatLiqTmp(:,lane)=real(1000*lane+49,kind_noahmp)
    workspaces(1)%SoilThermalProperty%SoilIceTmp(:,lane)=real(1000*lane+50,kind_noahmp)
    workspaces(1)%SoilWaterDiffusionRichards%DepthSnowSoilInv(:,lane)=real(1000*lane+51,kind_noahmp)
    workspaces(1)%SoilWaterDiffusionRichards%SoilThickTmp(:,lane)=real(1000*lane+52,kind_noahmp)
    workspaces(1)%SoilWaterDiffusionRichards%SoilWaterGrad(:,lane)=real(1000*lane+53,kind_noahmp)
    workspaces(1)%SoilWaterDiffusionRichards%WaterExcess(:,lane)=real(1000*lane+54,kind_noahmp)
    workspaces(1)%SoilWaterDiffusionRichards%SoilMoistureTmp(:,lane)=real(1000*lane+55,kind_noahmp)
    workspaces(1)%SoilWaterMain%MatRight(:,lane)=real(1000*lane+56,kind_noahmp)
    workspaces(1)%SoilWaterMain%MatLeft1(:,lane)=real(1000*lane+57,kind_noahmp)
    workspaces(1)%SoilWaterMain%MatLeft2(:,lane)=real(1000*lane+58,kind_noahmp)
    workspaces(1)%SoilWaterMain%MatLeft3(:,lane)=real(1000*lane+59,kind_noahmp)
    workspaces(1)%SoilWaterMain%SoilLiqTmp(:,lane)=real(1000*lane+60,kind_noahmp)
    workspaces(1)%SurfaceEnergyFluxGlacier%SoilIceTmp(:,lane)=real(1000*lane+61,kind_noahmp)
    workspaces(1)%SurfaceRadiation%RadSwAbsCanDir(:,lane)=real(1000*lane+62,kind_noahmp)
    workspaces(1)%SurfaceRadiation%RadSwAbsCanDif(:,lane)=real(1000*lane+63,kind_noahmp)
    workspaces(1)%TileDrainageHooghoudt%ThickSatZone(:,lane)=real(1000*lane+64,kind_noahmp)
    workspaces(1)%TileDrainageHooghoudt%LateralWatCondTmp(:,lane)=real(1000*lane+65,kind_noahmp)
    workspaces(1)%TileDrainageHooghoudt%WatExcFieldCapTmp(:,lane)=real(1000*lane+66,kind_noahmp)
    workspaces(1)%TileDrainageHooghoudt%SoilLiqWaterAftDrain(:,lane)=real(1000*lane+67,kind_noahmp)
    workspaces(1)%TileDrainageSimple%DrainFracTmp(:,lane)=real(1000*lane+68,kind_noahmp)
    workspaces(1)%TileDrainageSimple%SoilFieldCapLiq(:,lane)=real(1000*lane+69,kind_noahmp)
    workspaces(1)%TileDrainageSimple%WatExcFieldCap(:,lane)=real(1000*lane+70,kind_noahmp)
    workspaces(1)%WaterMainGlacier%SoilIceTmp(:,lane)=real(1000*lane+71,kind_noahmp)
    workspaces(1)%WaterMainGlacier%SoilLiqWaterTmp(:,lane)=real(1000*lane+72,kind_noahmp)
    workspaces(1)%WaterTableEquilibrium%DepthSoilFineLy(:,lane)=real(1000*lane+73,kind_noahmp)
  enddo
  ! sparse workspace subgroups: missing payloads must not unbatch siblings.
  deallocate(workspaces(1)%SoilMoistureSolver%MatRightTmp)
  allocate(workspaces(1)%SoilMoistureSolver%MatRightTmp(1:0,1:capacity))
  deallocate(workspaces(1)%WaterMainGlacier%SoilLiqWaterTmp)
  allocate(workspaces(1)%WaterMainGlacier%SoilLiqWaterTmp(1:0,1:capacity))
  call NoahmpGpuTraceStage('probe startup mapping: begin')
  call NoahmpGpuStart(columns,workspaces,active)
  call NoahmpGpuTraceStage('probe startup mapping: passed')
  errors=0
  !$acc parallel loop present(workspaces) reduction(+:errors)
  do lane=1,capacity
    if (lbound(workspaces(1)%GlacierPhaseChange%EnergyRes,2)/=1 .or. ubound(workspaces(1)%GlacierPhaseChange%EnergyRes,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierPhaseChange%GlacierPhaseChg,2)/=1 .or. ubound(workspaces(1)%GlacierPhaseChange%GlacierPhaseChg,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierPhaseChange%MassWatTotInit,2)/=1 .or. ubound(workspaces(1)%GlacierPhaseChange%MassWatTotInit,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierPhaseChange%MassWatIceInit,2)/=1 .or. ubound(workspaces(1)%GlacierPhaseChange%MassWatIceInit,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierPhaseChange%MassWatLiqInit,2)/=1 .or. ubound(workspaces(1)%GlacierPhaseChange%MassWatLiqInit,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierPhaseChange%MassWatIceTmp,2)/=1 .or. ubound(workspaces(1)%GlacierPhaseChange%MassWatIceTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierPhaseChange%MassWatLiqTmp,2)/=1 .or. ubound(workspaces(1)%GlacierPhaseChange%MassWatLiqTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierPhaseChange%EnergyResLeft,2)/=1 .or. ubound(workspaces(1)%GlacierPhaseChange%EnergyResLeft,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierTemperatureMain%MatRight,2)/=1 .or. ubound(workspaces(1)%GlacierTemperatureMain%MatRight,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierTemperatureMain%MatLeft1,2)/=1 .or. ubound(workspaces(1)%GlacierTemperatureMain%MatLeft1,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierTemperatureMain%MatLeft2,2)/=1 .or. ubound(workspaces(1)%GlacierTemperatureMain%MatLeft2,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierTemperatureMain%MatLeft3,2)/=1 .or. ubound(workspaces(1)%GlacierTemperatureMain%MatLeft3,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierTemperatureSolver%MatRightTmp,2)/=1 .or. ubound(workspaces(1)%GlacierTemperatureSolver%MatRightTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierTemperatureSolver%MatLeft3Tmp,2)/=1 .or. ubound(workspaces(1)%GlacierTemperatureSolver%MatLeft3Tmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierThermalDiffusion%DepthSnowSoilInv,2)/=1 .or. ubound(workspaces(1)%GlacierThermalDiffusion%DepthSnowSoilInv,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierThermalDiffusion%HeatCapacPerArea,2)/=1 .or. ubound(workspaces(1)%GlacierThermalDiffusion%HeatCapacPerArea,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierThermalDiffusion%TempGradDepth,2)/=1 .or. ubound(workspaces(1)%GlacierThermalDiffusion%TempGradDepth,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GlacierThermalDiffusion%EnergyExcess,2)/=1 .or. ubound(workspaces(1)%GlacierThermalDiffusion%EnergyExcess,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GroundWaterTopModel%DepthSoilMid,2)/=1 .or. ubound(workspaces(1)%GroundWaterTopModel%DepthSoilMid,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GroundWaterTopModel%ThicknessSoil,2)/=1 .or. ubound(workspaces(1)%GroundWaterTopModel%ThicknessSoil,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GroundWaterTopModel%SoilLiqTmp,2)/=1 .or. ubound(workspaces(1)%GroundWaterTopModel%SoilLiqTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GroundWaterTopModel%SoilEffPorosity,2)/=1 .or. ubound(workspaces(1)%GroundWaterTopModel%SoilEffPorosity,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GroundWaterTopModel%SoilWatConductTmp,2)/=1 .or. ubound(workspaces(1)%GroundWaterTopModel%SoilWatConductTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%GroundWaterTopModel%SoilMoisture,2)/=1 .or. ubound(workspaces(1)%GroundWaterTopModel%SoilMoisture,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%RunoffSurfaceFreeDrain%SoilWatMaxHold,2)/=1 .or. ubound(workspaces(1)%RunoffSurfaceFreeDrain%SoilWatMaxHold,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SnowLayerDivide%SnowIceTmp,2)/=1 .or. ubound(workspaces(1)%SnowLayerDivide%SnowIceTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SnowLayerDivide%SnowLiqTmp,2)/=1 .or. ubound(workspaces(1)%SnowLayerDivide%SnowLiqTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SnowLayerDivide%TemperatureSnowTmp,2)/=1 .or. ubound(workspaces(1)%SnowLayerDivide%TemperatureSnowTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SnowLayerDivide%SnowThickTmp,2)/=1 .or. ubound(workspaces(1)%SnowLayerDivide%SnowThickTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SnowThermalProperty%SnowDensBulk,2)/=1 .or. ubound(workspaces(1)%SnowThermalProperty%SnowDensBulk,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilMoistureSolver%MatRightTmp,2)/=1 .or. ubound(workspaces(1)%SoilMoistureSolver%MatRightTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilMoistureSolver%MatLeft3Tmp,2)/=1 .or. ubound(workspaces(1)%SoilMoistureSolver%MatLeft3Tmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowTemperatureMain%MatRight,2)/=1 .or. ubound(workspaces(1)%SoilSnowTemperatureMain%MatRight,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowTemperatureMain%MatLeft1,2)/=1 .or. ubound(workspaces(1)%SoilSnowTemperatureMain%MatLeft1,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowTemperatureMain%MatLeft2,2)/=1 .or. ubound(workspaces(1)%SoilSnowTemperatureMain%MatLeft2,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowTemperatureMain%MatLeft3,2)/=1 .or. ubound(workspaces(1)%SoilSnowTemperatureMain%MatLeft3,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowTemperatureSolver%MatRightTmp,2)/=1 .or. ubound(workspaces(1)%SoilSnowTemperatureSolver%MatRightTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowTemperatureSolver%MatLeft3Tmp,2)/=1 .or. ubound(workspaces(1)%SoilSnowTemperatureSolver%MatLeft3Tmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowThermalDiffusion%DepthSnowSoilInv,2)/=1 .or. ubound(workspaces(1)%SoilSnowThermalDiffusion%DepthSnowSoilInv,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowThermalDiffusion%HeatCapacPerArea,2)/=1 .or. ubound(workspaces(1)%SoilSnowThermalDiffusion%HeatCapacPerArea,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowThermalDiffusion%TempGradDepth,2)/=1 .or. ubound(workspaces(1)%SoilSnowThermalDiffusion%TempGradDepth,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowThermalDiffusion%EnergyExcess,2)/=1 .or. ubound(workspaces(1)%SoilSnowThermalDiffusion%EnergyExcess,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowWaterPhaseChange%EnergyRes,2)/=1 .or. ubound(workspaces(1)%SoilSnowWaterPhaseChange%EnergyRes,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowWaterPhaseChange%WaterPhaseChg,2)/=1 .or. ubound(workspaces(1)%SoilSnowWaterPhaseChange%WaterPhaseChg,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowWaterPhaseChange%MassWatTotInit,2)/=1 .or. ubound(workspaces(1)%SoilSnowWaterPhaseChange%MassWatTotInit,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowWaterPhaseChange%MassWatIceInit,2)/=1 .or. ubound(workspaces(1)%SoilSnowWaterPhaseChange%MassWatIceInit,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowWaterPhaseChange%MassWatLiqInit,2)/=1 .or. ubound(workspaces(1)%SoilSnowWaterPhaseChange%MassWatLiqInit,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowWaterPhaseChange%MassWatIceTmp,2)/=1 .or. ubound(workspaces(1)%SoilSnowWaterPhaseChange%MassWatIceTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilSnowWaterPhaseChange%MassWatLiqTmp,2)/=1 .or. ubound(workspaces(1)%SoilSnowWaterPhaseChange%MassWatLiqTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilThermalProperty%SoilIceTmp,2)/=1 .or. ubound(workspaces(1)%SoilThermalProperty%SoilIceTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilWaterDiffusionRichards%DepthSnowSoilInv,2)/=1 .or. ubound(workspaces(1)%SoilWaterDiffusionRichards%DepthSnowSoilInv,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilWaterDiffusionRichards%SoilThickTmp,2)/=1 .or. ubound(workspaces(1)%SoilWaterDiffusionRichards%SoilThickTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilWaterDiffusionRichards%SoilWaterGrad,2)/=1 .or. ubound(workspaces(1)%SoilWaterDiffusionRichards%SoilWaterGrad,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilWaterDiffusionRichards%WaterExcess,2)/=1 .or. ubound(workspaces(1)%SoilWaterDiffusionRichards%WaterExcess,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilWaterDiffusionRichards%SoilMoistureTmp,2)/=1 .or. ubound(workspaces(1)%SoilWaterDiffusionRichards%SoilMoistureTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilWaterMain%MatRight,2)/=1 .or. ubound(workspaces(1)%SoilWaterMain%MatRight,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilWaterMain%MatLeft1,2)/=1 .or. ubound(workspaces(1)%SoilWaterMain%MatLeft1,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilWaterMain%MatLeft2,2)/=1 .or. ubound(workspaces(1)%SoilWaterMain%MatLeft2,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilWaterMain%MatLeft3,2)/=1 .or. ubound(workspaces(1)%SoilWaterMain%MatLeft3,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SoilWaterMain%SoilLiqTmp,2)/=1 .or. ubound(workspaces(1)%SoilWaterMain%SoilLiqTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SurfaceEnergyFluxGlacier%SoilIceTmp,2)/=1 .or. ubound(workspaces(1)%SurfaceEnergyFluxGlacier%SoilIceTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SurfaceRadiation%RadSwAbsCanDir,2)/=1 .or. ubound(workspaces(1)%SurfaceRadiation%RadSwAbsCanDir,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%SurfaceRadiation%RadSwAbsCanDif,2)/=1 .or. ubound(workspaces(1)%SurfaceRadiation%RadSwAbsCanDif,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%TileDrainageHooghoudt%ThickSatZone,2)/=1 .or. ubound(workspaces(1)%TileDrainageHooghoudt%ThickSatZone,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%TileDrainageHooghoudt%LateralWatCondTmp,2)/=1 .or. ubound(workspaces(1)%TileDrainageHooghoudt%LateralWatCondTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%TileDrainageHooghoudt%WatExcFieldCapTmp,2)/=1 .or. ubound(workspaces(1)%TileDrainageHooghoudt%WatExcFieldCapTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%TileDrainageHooghoudt%SoilLiqWaterAftDrain,2)/=1 .or. ubound(workspaces(1)%TileDrainageHooghoudt%SoilLiqWaterAftDrain,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%TileDrainageSimple%DrainFracTmp,2)/=1 .or. ubound(workspaces(1)%TileDrainageSimple%DrainFracTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%TileDrainageSimple%SoilFieldCapLiq,2)/=1 .or. ubound(workspaces(1)%TileDrainageSimple%SoilFieldCapLiq,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%TileDrainageSimple%WatExcFieldCap,2)/=1 .or. ubound(workspaces(1)%TileDrainageSimple%WatExcFieldCap,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%WaterMainGlacier%SoilIceTmp,2)/=1 .or. ubound(workspaces(1)%WaterMainGlacier%SoilIceTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%WaterMainGlacier%SoilLiqWaterTmp,2)/=1 .or. ubound(workspaces(1)%WaterMainGlacier%SoilLiqWaterTmp,2)/=capacity) errors=errors+1
    if (lbound(workspaces(1)%WaterTableEquilibrium%DepthSoilFineLy,2)/=1 .or. ubound(workspaces(1)%WaterTableEquilibrium%DepthSoilFineLy,2)/=capacity) errors=errors+1
    workspaces(1)%WaterTableEquilibrium%DepthSoilFineLy(1,lane)=real(1000*lane+73,kind_noahmp)
    if (any(workspaces(1)%GlacierPhaseChange%EnergyRes(:,lane)/=real(1000*lane+1,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierPhaseChange%GlacierPhaseChg(:,lane)/=real(1000*lane+2,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierPhaseChange%MassWatTotInit(:,lane)/=real(1000*lane+3,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierPhaseChange%MassWatIceInit(:,lane)/=real(1000*lane+4,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierPhaseChange%MassWatLiqInit(:,lane)/=real(1000*lane+5,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierPhaseChange%MassWatIceTmp(:,lane)/=real(1000*lane+6,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierPhaseChange%MassWatLiqTmp(:,lane)/=real(1000*lane+7,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierPhaseChange%EnergyResLeft(:,lane)/=real(1000*lane+8,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierTemperatureMain%MatRight(:,lane)/=real(1000*lane+9,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierTemperatureMain%MatLeft1(:,lane)/=real(1000*lane+10,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierTemperatureMain%MatLeft2(:,lane)/=real(1000*lane+11,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierTemperatureMain%MatLeft3(:,lane)/=real(1000*lane+12,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierTemperatureSolver%MatRightTmp(:,lane)/=real(1000*lane+13,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierTemperatureSolver%MatLeft3Tmp(:,lane)/=real(1000*lane+14,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierThermalDiffusion%DepthSnowSoilInv(:,lane)/=real(1000*lane+15,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierThermalDiffusion%HeatCapacPerArea(:,lane)/=real(1000*lane+16,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierThermalDiffusion%TempGradDepth(:,lane)/=real(1000*lane+17,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GlacierThermalDiffusion%EnergyExcess(:,lane)/=real(1000*lane+18,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GroundWaterTopModel%DepthSoilMid(:,lane)/=real(1000*lane+19,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GroundWaterTopModel%ThicknessSoil(:,lane)/=real(1000*lane+20,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GroundWaterTopModel%SoilLiqTmp(:,lane)/=real(1000*lane+21,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GroundWaterTopModel%SoilEffPorosity(:,lane)/=real(1000*lane+22,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GroundWaterTopModel%SoilWatConductTmp(:,lane)/=real(1000*lane+23,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%GroundWaterTopModel%SoilMoisture(:,lane)/=real(1000*lane+24,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%RunoffSurfaceFreeDrain%SoilWatMaxHold(:,lane)/=real(1000*lane+25,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SnowLayerDivide%SnowIceTmp(:,lane)/=real(1000*lane+26,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SnowLayerDivide%SnowLiqTmp(:,lane)/=real(1000*lane+27,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SnowLayerDivide%TemperatureSnowTmp(:,lane)/=real(1000*lane+28,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SnowLayerDivide%SnowThickTmp(:,lane)/=real(1000*lane+29,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SnowThermalProperty%SnowDensBulk(:,lane)/=real(1000*lane+30,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilMoistureSolver%MatRightTmp(:,lane)/=real(1000*lane+31,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilMoistureSolver%MatLeft3Tmp(:,lane)/=real(1000*lane+32,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowTemperatureMain%MatRight(:,lane)/=real(1000*lane+33,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowTemperatureMain%MatLeft1(:,lane)/=real(1000*lane+34,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowTemperatureMain%MatLeft2(:,lane)/=real(1000*lane+35,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowTemperatureMain%MatLeft3(:,lane)/=real(1000*lane+36,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowTemperatureSolver%MatRightTmp(:,lane)/=real(1000*lane+37,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowTemperatureSolver%MatLeft3Tmp(:,lane)/=real(1000*lane+38,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowThermalDiffusion%DepthSnowSoilInv(:,lane)/=real(1000*lane+39,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowThermalDiffusion%HeatCapacPerArea(:,lane)/=real(1000*lane+40,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowThermalDiffusion%TempGradDepth(:,lane)/=real(1000*lane+41,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowThermalDiffusion%EnergyExcess(:,lane)/=real(1000*lane+42,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowWaterPhaseChange%EnergyRes(:,lane)/=real(1000*lane+43,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowWaterPhaseChange%WaterPhaseChg(:,lane)/=real(1000*lane+44,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowWaterPhaseChange%MassWatTotInit(:,lane)/=real(1000*lane+45,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowWaterPhaseChange%MassWatIceInit(:,lane)/=real(1000*lane+46,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowWaterPhaseChange%MassWatLiqInit(:,lane)/=real(1000*lane+47,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowWaterPhaseChange%MassWatIceTmp(:,lane)/=real(1000*lane+48,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilSnowWaterPhaseChange%MassWatLiqTmp(:,lane)/=real(1000*lane+49,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilThermalProperty%SoilIceTmp(:,lane)/=real(1000*lane+50,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilWaterDiffusionRichards%DepthSnowSoilInv(:,lane)/=real(1000*lane+51,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilWaterDiffusionRichards%SoilThickTmp(:,lane)/=real(1000*lane+52,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilWaterDiffusionRichards%SoilWaterGrad(:,lane)/=real(1000*lane+53,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilWaterDiffusionRichards%WaterExcess(:,lane)/=real(1000*lane+54,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilWaterDiffusionRichards%SoilMoistureTmp(:,lane)/=real(1000*lane+55,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilWaterMain%MatRight(:,lane)/=real(1000*lane+56,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilWaterMain%MatLeft1(:,lane)/=real(1000*lane+57,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilWaterMain%MatLeft2(:,lane)/=real(1000*lane+58,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilWaterMain%MatLeft3(:,lane)/=real(1000*lane+59,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SoilWaterMain%SoilLiqTmp(:,lane)/=real(1000*lane+60,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SurfaceEnergyFluxGlacier%SoilIceTmp(:,lane)/=real(1000*lane+61,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SurfaceRadiation%RadSwAbsCanDir(:,lane)/=real(1000*lane+62,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%SurfaceRadiation%RadSwAbsCanDif(:,lane)/=real(1000*lane+63,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%TileDrainageHooghoudt%ThickSatZone(:,lane)/=real(1000*lane+64,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%TileDrainageHooghoudt%LateralWatCondTmp(:,lane)/=real(1000*lane+65,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%TileDrainageHooghoudt%WatExcFieldCapTmp(:,lane)/=real(1000*lane+66,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%TileDrainageHooghoudt%SoilLiqWaterAftDrain(:,lane)/=real(1000*lane+67,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%TileDrainageSimple%DrainFracTmp(:,lane)/=real(1000*lane+68,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%TileDrainageSimple%SoilFieldCapLiq(:,lane)/=real(1000*lane+69,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%TileDrainageSimple%WatExcFieldCap(:,lane)/=real(1000*lane+70,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%WaterMainGlacier%SoilIceTmp(:,lane)/=real(1000*lane+71,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%WaterMainGlacier%SoilLiqWaterTmp(:,lane)/=real(1000*lane+72,kind_noahmp))) errors=errors+1
    if (any(workspaces(1)%WaterTableEquilibrium%DepthSoilFineLy(:,lane)/=real(1000*lane+73,kind_noahmp))) errors=errors+1
  enddo
  !$acc end parallel loop
  if (errors/=0) error stop 'Actual-type probe: scratch component mapping failed'
  print '(a)', 'PASS: all 73 block scratch fields and independent workspace lanes.'

  do p=1,ncolumns
    if (active(p)) columns(p)%forcing%TemperatureAirRefHeight=real(p+300,kind_noahmp)
  enddo
  call NoahmpGpuTraceStage('probe boundary upload: begin')
  call NoahmpGpuUploadRefresh(columns)
  call NoahmpGpuTraceStage('probe boundary upload: completed')
  ! Mutate host input storage so the download must restore the device values.
  do p=1,ncolumns
    if (active(p)) columns(p)%forcing%TemperatureAirRefHeight=-123.0_kind_noahmp
  enddo
  call NoahmpGpuTraceStage('probe boundary download: begin')
  call NoahmpGpuBeforeRefresh(columns)
  call NoahmpGpuTraceStage('probe boundary download: completed')
  do p=1,ncolumns
    if (.not. active(p)) cycle
    if (columns(p)%forcing%TemperatureAirRefHeight/=real(p+300,kind_noahmp)) &
      error stop 'Actual-type probe: boundary transfer failed'
  enddo
  device_flag=0
  !$acc serial copyout(device_flag)
  if (acc_on_device(acc_device_nvidia)) device_flag=1
  !$acc end serial
  if (device_flag/=1) error stop 'Actual-type probe did not run on NVIDIA'
  call NoahmpGpuTraceStage('probe independent column writes: begin')
  !$acc parallel loop present(columns,active)
  do p=1,ncolumns
    if (.not. active(p)) cycle
    columns(p)%water%state%SoilMoisture(:)=real(p+1000,kind_noahmp)
    columns(p)%energy%state%TemperatureSoilSnow(:)=real(p+2000,kind_noahmp)
  enddo
  !$acc end parallel loop
  call NoahmpGpuTraceStage('probe independent column writes: completed')
  errors=0
  !$acc serial present(columns) copy(errors)
  if (any(columns(1)%energy%flux%FracRadSwAbsSnowDir/=17.0_kind_noahmp)) errors=errors+1
  if (any(columns(1)%energy%flux%FracRadSwAbsSnowDif/=23.0_kind_noahmp)) errors=errors+1
  if (lbound(columns(1)%energy%flux%FracRadSwAbsSnowDir,1)/=-2) errors=errors+1
  if (lbound(columns(1)%energy%flux%FracRadSwAbsSnowDir,2)/=2) errors=errors+1
  if (size(columns(1)%water%flux%TranspWatLossSoilMean)/=0) errors=errors+1
  if (allocated(columns(1)%biochem%param%LeafDeathTempCoeffCrop)) errors=errors+1
  columns(1)%energy%flux%FracRadSwAbsSnowDir(:,:)=31.0_kind_noahmp
  columns(1)%energy%flux%FracRadSwAbsSnowDif(:,:)=37.0_kind_noahmp
  !$acc end serial
  !$acc serial present(workspaces) copy(errors)
  if (size(workspaces(1)%SoilMoistureSolver%MatRightTmp)/=0) errors=errors+1
  if (size(workspaces(1)%WaterMainGlacier%SoilLiqWaterTmp)/=0) errors=errors+1
  !$acc end serial
  if (errors/=0) error stop 'Actual-type probe: optional/empty/unallocated descriptors failed'

  call NoahmpGpuTraceStage('probe snapshot download: begin')
  call NoahmpGpuSyncHost(columns)
  call NoahmpGpuTraceStage('probe snapshot download: completed')
  do p=1,ncolumns
    if (.not. active(p)) then
      if (any(columns(p)%water%state%SoilMoisture/=real(p,kind_noahmp))) &
        error stop 'Actual-type probe changed a skipped column'
      cycle
    endif
    if (any(columns(p)%water%state%SoilMoisture/=real(p+1000,kind_noahmp))) &
      error stop 'Actual-type probe: independent column payload/snapshot failed'
    if (any(columns(p)%energy%state%TemperatureSoilSnow/=real(p+2000,kind_noahmp))) &
      error stop 'Actual-type probe: snow/soil payload/snapshot failed'
    if (lbound(columns(p)%energy%state%TemperatureSoilSnow,1)/=-nsnow+1) &
      error stop 'Actual-type probe: host snow bounds were corrupted'
  enddo
  if (any(columns(1)%energy%flux%FracRadSwAbsSnowDir/=31.0_kind_noahmp)) &
    error stop 'Actual-type probe: first matrix snapshot failed'
  if (any(columns(1)%energy%flux%FracRadSwAbsSnowDif/=37.0_kind_noahmp)) &
    error stop 'Actual-type probe: second matrix snapshot failed'
  if (lbound(columns(1)%energy%flux%FracRadSwAbsSnowDir,2)/=2) &
    error stop 'Actual-type probe: host matrix bounds lost'
  print '(a)', 'PASS: optional matrices, empty arrays, unallocated components and snapshots.'

  ! Exercise valid, invalid and endpoint albedo updates between output events.
  ! Start the next refresh cycle so the following device writes require a new
  ! full snapshot; the previous SyncHost marked the old host state current.
  call NoahmpGpuUploadRefresh(columns)
  errors=0
  !$acc parallel loop present(columns,active) reduction(+:errors)
  do p=1,ncolumns
    if (.not. active(p)) cycle
    if (columns(p)%energy%state%AlbedoSfcOutput/=0.125_kind_noahmp) errors=errors+1
    columns(p)%energy%state%AlbedoSfc=0.25_kind_noahmp
    call NoahmpUpdateOutputAlbedo(columns(p))
    call NoahmpResetStepTemporaries(columns(p))
    if (columns(p)%energy%state%AlbedoSfcOutput/=0.25_kind_noahmp) errors=errors+1
    columns(p)%energy%state%AlbedoSfc=-9999.0_kind_noahmp
    call NoahmpUpdateOutputAlbedo(columns(p))
    if (columns(p)%energy%state%AlbedoSfcOutput/=0.25_kind_noahmp) errors=errors+1
    columns(p)%energy%state%AlbedoSfc=1.1_kind_noahmp
    call NoahmpUpdateOutputAlbedo(columns(p))
    if (columns(p)%energy%state%AlbedoSfcOutput/=0.25_kind_noahmp) errors=errors+1
    columns(p)%energy%state%AlbedoSfc=0.0_kind_noahmp
    call NoahmpUpdateOutputAlbedo(columns(p))
    if (columns(p)%energy%state%AlbedoSfcOutput/=0.0_kind_noahmp) errors=errors+1
    columns(p)%energy%state%AlbedoSfc=1.0_kind_noahmp
    call NoahmpUpdateOutputAlbedo(columns(p))
    if (columns(p)%energy%state%AlbedoSfcOutput/=1.0_kind_noahmp) errors=errors+1
    columns(p)%energy%state%AlbedoSfc=real(p,kind_noahmp)/real(ncolumns,kind_noahmp)
    call NoahmpUpdateOutputAlbedo(columns(p))
    columns(p)%energy%state%AlbedoSfc=-1.0_kind_noahmp
    call NoahmpUpdateOutputAlbedo(columns(p))
  enddo
  !$acc end parallel loop
  if (errors/=0) error stop 'Community probe: albedo bounds or persistence failed'
  call NoahmpGpuSyncHost(columns)
  do p=1,ncolumns
    if (.not. active(p)) cycle
    if (columns(p)%energy%state%AlbedoSfcOutput/=real(p,kind_noahmp)/real(ncolumns,kind_noahmp)) &
      error stop 'Community probe: output albedo snapshot failed'
    if (columns(p)%energy%state%AlbedoSfc/=-1.0_kind_noahmp) &
      error stop 'Community probe: raw physics albedo changed'
  enddo
  print '(a)', 'PASS: community albedo safeguard, timestep persistence and GPU snapshot.'

  call NoahmpGpuTraceStage('probe teardown: begin')
  call NoahmpGpuFinish(columns,workspaces,active)
  call NoahmpGpuTraceStage('probe teardown: completed')
  if (size(workspaces(1)%SoilMoistureSolver%MatRightTmp)/=0) &
    error stop 'Actual-type probe: sparse workspace host bounds changed at teardown'
  if (size(workspaces(1)%WaterMainGlacier%SoilLiqWaterTmp)/=0) &
    error stop 'Actual-type probe: second sparse workspace host bounds changed at teardown'
  print '(a)', 'PASS: sparse workspace subgroups, complete sibling payloads and teardown.'

  do p=1,ncolumns
    if (.not. allocated(columns(p)%water%state%SoilMoisture)) error stop 'Host descriptor lost during GPU teardown'
  enddo
  deallocate(columns,workspaces,active)
  print '(a)', 'PASS: actual Noah-MP GPU mapping, bounds, transfers, masks, partial blocks, and teardown.'
end program noahmp_gpu_mapping_probe
