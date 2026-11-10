#ifdef WRF_HYDRO
#error "The offline CPU column port does not support WRF_HYDRO"
#endif
module module_NoahMP_hrldas_driver

  use NoahmpVarType, only: noahmp_type
  use Machine
  use NoahmpGpuDataMod, only: NoahmpGpuBeforeRefresh,NoahmpGpuSyncHost,NoahmpGpuFinish
  use NoahmpWorkspaceMod, only: NoahmpWorkspace_type, NoahmpAllocateWorkspacePool
  use HrldasContextMod
  use HrldasSetupMod
  use NoahmpReadNamelistMod
  use NoahmpReadTableMod
  use module_hrldas_netcdf_io
  use NoahmpInitMainMod
  use NoahmpDriverMainMod
  use HrldasColumnInitMod
  use HrldasOutputMod
  use NoahmpColumnResetMod
  use module_date_utilities

#ifdef MPP_LAND
  use module_mpp_land, only: MPP_LAND_PAR_INI, mpp_land_init, getLocalXY, mpp_land_bcast_char
  use module_mpp_land, only: check_land, my_id , node_info
#endif

  implicit none

#ifdef MPP_LAND
  include "mpif.h"
#endif


  character(len=9), parameter :: version = "v20220701"
  integer :: LDASIN_VERSION

#ifdef MPP_LAND
  integer ix_tmp, jx_tmp
#endif

  type(noahmp_type), allocatable :: noahmp(:)
  type(NoahmpWorkspace_type), allocatable :: workspaces(:)
  type(HrldasContext_type)                                    :: ctx

  integer                                                :: I
  integer                                                :: J
  integer                                                :: imode

contains

  subroutine land_driver_ini(NTIME_out,wrfits,wrfite,wrfjts,wrfjte)

    implicit none

    integer      :: NTIME_out
    real(kind=kind_noahmp) :: max_landtype2d
    integer, optional, intent(in) :: wrfits,wrfite,wrfjts,wrfjte

    !local
    integer ::   xstart, xend, ystart, yend, ixfull, jxfull

#ifdef MPP_LAND
    call  MPP_LAND_INIT()
#endif

!---------------------------------------------------------------------
!  read in input data from table and initial file
!---------------------------------------------------------------------

    allocate(ctx%setup)
    call NoahmpReadNamelist(ctx)


! Explicit urban models are excluded from the standalone column port.

!----------------------------------------------------------------------
! Initialize gridded domain
!----------------------------------------------------------------------

#ifdef MPP_LAND
       call read_dim(ctx%runtime%hrldas_setup_file,ix_tmp,jx_tmp)
       call MPP_LAND_PAR_INI(1,ix_tmp,jx_tmp,1)
       call getLocalXY(ix_tmp,jx_tmp,ctx%grid%xstart,ctx%grid%ystart,ctx%grid%xend,ctx%grid%yend)
#endif

    call read_hrldas_hdrinfo(ctx%runtime%hrldas_setup_file, ctx%grid%ix, ctx%grid%jx,           &
                             ctx%grid%xstart, ctx%grid%xend, ctx%grid%ystart, ctx%grid%yend, &
                             ctx%grid%iswater, ctx%grid%islake, ctx%grid%isurban,            &
                             ctx%grid%isice, ctx%grid%llanduse, ctx%grid%dx, ctx%grid%dy,    &
                             ctx%grid%truelat1, ctx%grid%truelat2, ctx%grid%cen_lon,         &
                             ctx%grid%lat1, ctx%grid%lon1, ctx%grid%igrid, ctx%grid%mapproj  &
                            )

    write (ctx%grid%hgrid,'(I1)') ctx%grid%igrid
    write (ctx%runtime%olddate,'(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2)') &
           ctx%runtime%start_year, ctx%runtime%start_month, ctx%runtime%start_day, ctx%runtime%start_hour, &
      & ctx%runtime%start_min, 0

    ctx%runtime%startdate = ctx%runtime%olddate

#ifdef MPP_LAND
    ctx%grid%ix = ix_tmp
    ctx%grid%jx = jx_tmp
#endif


    ctx%grid%ids = ctx%grid%xstart
    ctx%grid%ide = ctx%grid%xend
    ctx%grid%jds = ctx%grid%ystart
    ctx%grid%jde = ctx%grid%yend
    ctx%grid%kds = 1
    ctx%grid%kde = 2
    ctx%grid%its = ctx%grid%xstart
    ctx%grid%ite = ctx%grid%xend
    ctx%grid%jts = ctx%grid%ystart
    ctx%grid%jte = ctx%grid%yend
    ctx%grid%kts = 1
    ctx%grid%kte = 2
    ctx%grid%ims = ctx%grid%xstart
    ctx%grid%ime = ctx%grid%xend
    ctx%grid%jms = ctx%grid%ystart
    ctx%grid%jme = ctx%grid%yend
    ctx%grid%kms = 1
    ctx%grid%kme = 2


!---------------------------------------------------------------------
!  Allocate multi-dimension fields for subwindow calculation
!---------------------------------------------------------------------

    ctx%grid%ixfull    = ctx%grid%xend - ctx%grid%xstart+1
    ctx%grid%jxfull    = ctx%grid%yend - ctx%grid%ystart+1
    ctx%grid%ixpar     = ctx%grid%ixfull
    ctx%grid%jxpar     = ctx%grid%jxfull
    ctx%grid%xstartpar = 1
    ctx%grid%ystartpar = 1

!---------------------------------------------------------------------
!  initialize HrldasContextMod data type and variables
!---------------------------------------------------------------------

    call HrldasAllocateSetup(ctx)

!---------------------------------------------------------------------
!  read in Noah-MP Parameters from MPTABLE
!---------------------------------------------------------------------

    call NoahmpReadTable(ctx)


!----------------------------------------------------------------------
! Read Landuse Type and Soil Texture and Other Information
!----------------------------------------------------------------------

    call READLAND_HRLDAS(ctx%runtime%HRLDAS_SETUP_FILE, ctx%grid%XSTART, ctx%grid%XEND, ctx%grid%YSTART, ctx%grid%YEND, &
      & &
                         ctx%grid%ISWATER, ctx%grid%ISLAKE, ctx%grid%IVGTYP, ctx%grid%ISLTYP, ctx%setup%TERRAIN, ctx%grid%TMN, &
                         ctx%grid%XLAT, ctx%grid%XLONG, ctx%grid%XLAND, ctx%setup%SEAICE, ctx%setup%MSFTX, ctx%setup%MSFTY)

    where (ctx%setup%SEAICE > 0.0) ctx%grid%XICE = 1.0

!------------------------------------------------------------------------
! For spatially-varying soil parameters, read in necessary extra fields
!------------------------------------------------------------------------

    if ( ctx%runtime%IOPT_SOIL == 2 ) then
       call READ_SOIL_TEXTURE(ctx%runtime%HRLDAS_SETUP_FILE, ctx%grid%XSTART, ctx%grid%XEND, ctx%grid%YSTART, ctx%grid%YEND, &
                              ctx%runtime%NSOIL, ctx%grid%IVGTYP, ctx%grid%SOILCL1, ctx%grid%SOILCL2, ctx%grid%SOILCL3,      &
                              ctx%grid%SOILCL4, ctx%grid%ISICE, ctx%grid%ISWATER)
       ctx%grid%ISLTYP = nint(ctx%grid%SOILCL1(:,:))
    endif

    if ( ctx%runtime%IOPT_SOIL == 3 ) then
       call READ_SOIL_COMPOSITION(ctx%runtime%HRLDAS_SETUP_FILE, ctx%grid%XSTART, ctx%grid%XEND, ctx%grid%YSTART, &
      & ctx%grid%YEND, &
                                  ctx%runtime%NSOIL, ctx%grid%IVGTYP, ctx%grid%ISICE, ctx%grid%ISWATER, ctx%grid%SOILCOMP)
    endif

    if ( ctx%runtime%IOPT_SOIL == 4 ) then
       call READ_3D_SOIL(ctx%runtime%SPATIAL_FILENAME, ctx%grid%XSTART, ctx%grid%XEND, ctx%grid%YSTART, ctx%grid%YEND,       &
                        ctx%runtime%NSOIL, ctx%grid%BEXP_3D, ctx%grid%SMCDRY_3D, ctx%grid%SMCWLT_3D, ctx%grid%SMCREF_3D,     &
                        ctx%grid%SMCMAX_3D, ctx%grid%DKSAT_3D, ctx%grid%DWSAT_3D, ctx%grid%PSISAT_3D, ctx%grid%QUARTZ_3D, &
                        ctx%grid%REFDK_2D, ctx%grid%REFKDT_2D, ctx%grid%IRR_FRAC_2D, ctx%grid%IRR_HAR_2D,                 &
                        ctx%grid%IRR_LAI_2D, ctx%grid%IRR_MAD_2D, ctx%grid%FILOSS_2D, ctx%grid%SPRIR_RATE_2D,             &
                        ctx%grid%MICIR_RATE_2D, ctx%grid%FIRTFAC_2D, ctx%grid%IR_RAIN_2D, ctx%grid%BVIC_2D,               &
                        ctx%grid%AXAJ_2D, ctx%grid%BXAJ_2D, ctx%grid%XXAJ_2D, ctx%grid%BDVIC_2D, ctx%grid%GDVIC_2D,       &
                        ctx%grid%BBVIC_2D, ctx%grid%KLAT_FAC, ctx%grid%TDSMC_FAC, ctx%grid%TD_DC, ctx%grid%TD_DCOEF,      &
                        ctx%grid%TD_DDRAIN, ctx%grid%TD_RADI, ctx%grid%TD_SPAC)
    endif

!----------------------------------------------------------------------
! For spatially-varying irrigation parameters, read in necessary extra fields
! only if IOPT_IRR >= 1
!----------------------------------------------------------------------
    if ( ctx%runtime%IOPT_IRR >= 1 ) then
       call READ_AGRICULTURE_DATA(ctx%runtime%AGDATA_FLNM, ctx%grid%XSTART, ctx%grid%XEND, ctx%grid%YSTART, ctx%grid%YEND, &
                                  ctx%grid%IRFRACT, ctx%grid%SIFRACT, ctx%grid%MIFRACT, ctx%grid%FIFRACT)
    endif

!------------------------------------------------------------------------
! For IOPT_RUNSUB = 5 (MMF groundwater), read in necessary extra fields
! This option is not tested for parallel use in the offline driver
!------------------------------------------------------------------------


!------------------------------------------------------------------------
! For IOPT_ALB = 3 (SNICAR), read in necessary fileds and parameters
!------------------------------------------------------------------------

!------------------------------------------------------------------------
! For OPT_CROP=1 (crop model), read in necessary extra fields
!------------------------------------------------------------------------

    ctx%setup%CROPTYPE   = 0       ! make default 0% crops everywhere
    if (ctx%runtime%IOPT_CROP == 1) then ! use crop_option=1 for reading crop extra fields
       call READ_CROP_INPUT(ctx%runtime%HRLDAS_SETUP_FILE, ctx%grid%XSTART, ctx%grid%XEND, ctx%grid%YSTART, ctx%grid%YEND,&
                            ctx%setup%CROPTYPE, ctx%grid%PLANTING, ctx%grid%HARVEST, ctx%grid%SEASON_GDD)
    endif

!------------------------------------------------------------------------
! For IOPT_TDRN = 1 or 2 READ TILE DRAIN MAP
!------------------------------------------------------------------------

    ctx%grid%TD_FRACTION = 0.0
    if (ctx%runtime%IOPT_TDRN > 0) then
       call READ_TILE_DRAIN_MAP(ctx%runtime%TDINPUT_FLNM, ctx%grid%XSTART, ctx%grid%XEND, ctx%grid%YSTART, ctx%grid%YEND, &
      & ctx%grid%TD_FRACTION)
    endif

!------------------------------------------------------------------------
! For surface wetland scheme, IOPT_WETLAND = 1 or 2, read in necessary extra fields
!------------------------------------------------------------------------
    if (ctx%runtime%IOPT_WETLAND == 2) then ! use wetland_option=2 for reading wetland extra fields
       ctx%grid%FSATMX = 0.38
       ctx%grid%WCAP   = 0.1
       call READ_WETLAND_INPUT(ctx%runtime%HRLDAS_SETUP_FILE, ctx%grid%XSTART, ctx%grid%XEND, ctx%grid%YSTART, ctx%grid%YEND,&
                               ctx%grid%FSATMX,ctx%grid%WCAP)
    endif

!------------------------------------------------------------------------
! For SF_URBAN_PHYSICS > 0 read 2D map of urban parameters
!------------------------------------------------------------------------


!----------------------------------------------------------------------
! Initialize Model State
!----------------------------------------------------------------------

    ctx%runtime%SLOPETYP  =  1 ! it was 2 here and 1 in the noahmpdrv- pvk
    ctx%runtime%DZS       =  ctx%runtime%SOIL_THICK_INPUT(1:ctx%runtime%NSOIL)
    ctx%runtime%ITIMESTEP = 1
    xstart = ctx%grid%xstart
    xend   = ctx%grid%xend
    ystart = ctx%grid%ystart
    yend   = ctx%grid%yend
    ixfull = ctx%grid%ixfull
    jxfull = ctx%grid%jxfull

    if (ctx%runtime%restart_filename_requested /= " ") then
       ctx%runtime%restart_flag = .true.

       call find_restart_file(ctx%runtime%rank, trim(ctx%runtime%restart_filename_requested), ctx%runtime%startdate, &
                              ctx%runtime%khour, ctx%runtime%olddate, ctx%runtime%restart_flnm)

       call read_restart(trim(ctx%runtime%restart_flnm), xstart, xend, xstart, ixfull, jxfull, ctx%runtime%NSOIL, &
      & ctx%runtime%olddate)

#ifdef MPP_LAND
       call mpp_land_bcast_char(19,ctx%runtime%OLDDATE(1:19))
#endif

       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SOIL_T"  , ctx%setup%TSLB     )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNOW_T"  , ctx%setup%TSNOXY   )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SMC"     , ctx%setup%SMOIS    )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SH2O"    , ctx%setup%SH2O     )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ZSNSO"   , ctx%setup%ZSNSOXY  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNICE"   , ctx%setup%SNICEXY  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNLIQ"   , ctx%setup%SNLIQXY  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "FWET"    , ctx%setup%FWETXY   )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNEQVO"  , ctx%setup%SNEQVOXY )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "EAH"     , ctx%setup%EAHXY    )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "TAH"     , ctx%setup%TAHXY    )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ALBOLD"  , ctx%setup%ALBOLDXY )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "CM"      , ctx%setup%CMXY     )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "CH"      , ctx%setup%CHXY     )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ISNOW"   , ctx%setup%ISNOWXY  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "CANLIQ"  , ctx%setup%CANLIQXY )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "CANICE"  , ctx%setup%CANICEXY )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNEQV"   , ctx%setup%SNOW     )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNOWH"   , ctx%setup%SNOWH    )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "TV"      , ctx%setup%TVXY     )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "TG"      , ctx%setup%TGXY     )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ZWT"     , ctx%setup%ZWTXY    )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "WA"      , ctx%setup%WAXY     )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "WT"      , ctx%setup%WTXY     )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "WSLAKE"  , ctx%setup%WSLAKEXY )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "LFMASS"  , ctx%setup%LFMASSXY )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "RTMASS"  , ctx%setup%RTMASSXY )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "STMASS"  , ctx%setup%STMASSXY )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "CROPCAT" , ctx%grid%CROPCAT  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "WOOD"    , ctx%setup%WOODXY   )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "GRAIN"   , ctx%setup%GRAINXY  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "GDD"     , ctx%setup%GDDXY    )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "STBLCP"  , ctx%setup%STBLCPXY )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "FASTCP"  , ctx%setup%FASTCPXY )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "LAI"     , ctx%input%LAI      )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SAI"     , ctx%setup%XSAIXY   )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "VEGFRA"  , ctx%input%VEGFRA   )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "GVFMIN"  , ctx%grid%GVFMIN   )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "GVFMAX"  , ctx%grid%GVFMAX   )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACMELT"  , ctx%setup%ACSNOM   )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACSNOW"  , ctx%setup%ACSNOW   )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "TAUSS"   , ctx%setup%TAUSSXY  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QSFC"    , ctx%setup%QSFC     )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SFCRUNOFF",ctx%setup%SFCRUNOFF)
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "UDRUNOFF" ,ctx%setup%UDRUNOFF )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QTDRAIN"  ,ctx%setup%QTDRAIN  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_SSOIL" , ctx%setup%ACC_SSOILXY )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_QINSUR", ctx%setup%ACC_QINSURXY)
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_QSEVA" , ctx%setup%ACC_QSEVAXY )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_ETRANI", ctx%setup%ACC_ETRANIXY)
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_DWATER", ctx%setup%ACC_DWATERXY)
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_PRCP"  , ctx%setup%ACC_PRCPXY  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_ECAN"  , ctx%setup%ACC_ECANXY  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_ETRAN" , ctx%setup%ACC_ETRANXY )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_EDIR"  , ctx%setup%ACC_EDIRXY  )
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACC_GLAFLW", ctx%setup%ACC_GLAFLWXY)
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ALBSOILDIR", ctx%setup%ALBSOILDIRXY)
       call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ALBSOILDIF", ctx%setup%ALBSOILDIFXY)

       ! below for SNICAR snow albedo scheme

       ! below for irrigation scheme
       if ( ctx%runtime%IOPT_IRR > 0 ) then
          call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRNUMSI" , ctx%setup%IRNUMSI  )
          call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRNUMMI" , ctx%setup%IRNUMMI  )
          call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRNUMFI" , ctx%setup%IRNUMFI  )
          call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRWATSI" , ctx%setup%IRWATSI  )
          call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRWATMI" , ctx%setup%IRWATMI  )
          call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRWATFI" , ctx%setup%IRWATFI  )
          call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRSIVOL" , ctx%setup%IRSIVOL  )
          call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRMIVOL" , ctx%setup%IRMIVOL  )
          call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRFIVOL" , ctx%setup%IRFIVOL  )
          call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRELOSS" , ctx%setup%IRELOSS  )
          call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "IRRSPLH" , ctx%setup%IRRSPLH  )
       endif

       ! below for MMF groundwater scheme

       ! for wetland scheme
       if ( ctx%runtime%IOPT_WETLAND > 0 ) then
          call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "FSATXY"     , ctx%setup%FSATXY    )
          call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "WSURFXY"    , ctx%setup%WSURFXY   )
       endif

       ! below for urban model

       ctx%runtime%STEPWTD = nint(ctx%runtime%WTDDT*60.0/ctx%runtime%DTBL)
       ctx%runtime%STEPWTD = max(ctx%runtime%STEPWTD,1)

! Must still call NoahmpInitMain even in restart to set up parameter arrays (also done in WRF)

!------------------------------------------------------------------------------------------
!    Initialize Noah-MP states and variables
!------------------------------------------------------------------------------------------
       call NoahmpInitMain(ctx)



    else

       ctx%runtime%restart_flag = .false.

          ctx%runtime%inflnm = trim(ctx%runtime%indir)//"/"//&
          ctx%runtime%startdate(1:4)//ctx%runtime%startdate(6:7)//ctx%runtime%startdate(9:10)//ctx%runtime%startdate(12:13)//&
          ".LDASIN_DOMAIN"//ctx%grid%hgrid

          call READINIT_HRLDAS(ctx%runtime%HRLDAS_SETUP_FILE, ctx%grid%xstart, ctx%grid%xend,    &
                              ctx%grid%ystart, ctx%grid%yend, ctx%runtime%NSOIL, ctx%runtime%DZS,  &
                              ctx%runtime%OLDDATE, LDASIN_VERSION, ctx%setup%SMOIS,              &
                              ctx%setup%TSLB, ctx%setup%CANWAT, ctx%setup%TSK, ctx%setup%SNOW,   &
                              ctx%setup%SNOWH, ctx%runtime%FNDSNOWH)

          ctx%input%VEGFRA =  0.0
          call READVEG_HRLDAS(ctx%runtime%HRLDAS_SETUP_FILE, ctx%grid%xstart, ctx%grid%xend,        &
                              ctx%grid%ystart, ctx%grid%yend, ctx%runtime%OLDDATE, ctx%grid%IVGTYP, &
                              ctx%input%VEGFRA, ctx%input%LAI, ctx%grid%GVFMIN, ctx%grid%GVFMAX)


       !----------------------------------------------------------------------------------
       !   Initialize Noah-MP variables
       !----------------------------------------------------------------------------------
       call NoahmpInitMain(ctx)

       !----------------------------------------------------------------------------------
       !   Initialize MMF Groundwater variables
       !----------------------------------------------------------------------------------



    endif  ! restart flag

    ctx%runtime%NTIME       = (ctx%runtime%KHOUR)*3600.0/nint(ctx%runtime%dtbl)*(ctx%runtime%spinup_loops+1)
    ctx%runtime%spinup_loop = 0
    ctx%runtime%reset_spinup_date = .false.
    ctx%runtime%reset_spinup_datea = .false.

    print*, "NTIME = ", ctx%runtime%NTIME , "KHOUR=",ctx%runtime%KHOUR,"dtbl = ", ctx%runtime%dtbl

    call system_clock(count=ctx%runtime%clock_count_1)   ! Start a timer



    call HrldasPrepareBoundary(ctx,noahmp)
    call NoahmpAllocateWorkspacePool(workspaces,min(ctx%runtime%block_size,size(noahmp)), &
         ctx%runtime%nsoil,ctx%runtime%nsnow,ctx%runtime%numrad)
    NTIME_out = ctx%runtime%NTIME

  end subroutine land_driver_ini

!===============================================================================
  subroutine land_driver_exe(itime)
     implicit  none
     integer :: itime          ! timestep loop
     integer :: I, J, K, p
     logical :: first_columns,reset_soil,fatal_error
     real(kind_noahmp) :: console_soil_temperature,console_leaf_area

!---------------------------------------------------------------------------------
! Read the forcing data.
!---------------------------------------------------------------------------------
! For HRLDAS, we're assuming (for now) that each time period is in a
! separate file.  So we can open a new one right now.

     ctx%runtime%inflnm = trim(ctx%runtime%indir)//"/"//&
          ctx%runtime%olddate(1:4)//ctx%runtime%olddate(6:7)//ctx%runtime%olddate(9:10)//ctx%runtime%olddate(12:13)//&
          ".LDASIN_DOMAIN"//ctx%grid%hgrid

     ! Build a filename template
     ctx%runtime%inflnm_template = trim(ctx%runtime%indir)//"/<date>.LDASIN_DOMAIN"//ctx%grid%hgrid

#ifdef MPP_LAND
     call mpp_land_bcast_char(19,ctx%runtime%OLDDATE(1:19))
#endif

     call READFORC_HRLDAS(ctx%runtime%INFLNM_TEMPLATE, ctx%runtime%FORCING_TIMESTEP,                               &
                          ctx%runtime%OLDDATE, ctx%grid%XSTART, ctx%grid%XEND, ctx%grid%YSTART, ctx%grid%YEND,  &
                          ctx%runtime%forcing_name_T,ctx%runtime%forcing_name_Q,ctx%runtime%forcing_name_U,           &
                          ctx%runtime%forcing_name_V,ctx%runtime%forcing_name_P,ctx%runtime%forcing_name_LW,          &
                          ctx%runtime%forcing_name_SW,ctx%runtime%forcing_name_PR,ctx%runtime%forcing_name_SN,        &
                          ctx%runtime%forcing_name_DirFrac,ctx%runtime%forcing_name_VisFrac,                       &
                          ctx%input%T_PHY(:,:),ctx%input%QV_CURR(:,:),ctx%input%U_PHY(:,:),               &
                          ctx%input%V_PHY(:,:),ctx%input%P8W(:,:), ctx%input%GLW, ctx%input%SWDOWN,          &
                          ctx%input%RAINBL, ctx%input%SNOWBL, ctx%input%VEGFRA, ctx%runtime%update_veg,            &
                          ctx%input%LAI, ctx%runtime%update_lai, ctx%runtime%reset_spinup_date, ctx%runtime%startdate, &
                          ctx%input%RadSwDirFrac, ctx%input%RadSwVisFrac )

       if(maxval(ctx%input%VEGFRA) <= 1 ) then
           ctx%input%VEGFRA = ctx%input%VEGFRA * 100.   ! added for input vegfra as a fraction (0~1)
       endif



     ! Preserve the accepted input arithmetic and its rounding order.
     where(ctx%grid%XLAND>1.5)
       ctx%input%T_PHY(:,:)=0.0
       ctx%input%U_PHY(:,:)=0.0
       ctx%input%V_PHY(:,:)=0.0
       ctx%input%QV_CURR(:,:)=0.0
       ctx%input%P8W(:,:)=0.0
       ctx%input%GLW=0.0
       ctx%input%SWDOWN=0.0
       ctx%input%RAINBL=0.0
       ctx%input%SNOWBL=0.0
     endwhere
     ctx%input%QV_CURR(:,:)=ctx%input%QV_CURR(:,:)/(1.0-ctx%input%QV_CURR(:,:))
     ctx%input%RAINBL=ctx%input%RAINBL*ctx%runtime%DTBL
     ctx%input%SNOWBL=ctx%input%SNOWBL*ctx%runtime%DTBL
     do J=ctx%grid%YSTART,ctx%grid%YEND
       do I=ctx%grid%XSTART,ctx%grid%XEND
         call CALC_DECLIN(ctx%runtime%OLDDATE(1:19),ctx%grid%XLAT(I,J),ctx%grid%XLONG(I,J), &
                          ctx%input%COSZEN(I,J),ctx%runtime%JULIAN)
         if (ctx%runtime%IOPT_IRR>0) call LOCAL_TIME(ctx%runtime%OLDDATE(1:19),ctx%grid%XLONG(I,J),ctx%grid%LOCTIM(I,J))
       enddo
     enddo
     read(ctx%runtime%OLDDATE(1:4),*) ctx%runtime%YR
     ctx%runtime%YEARLEN=365
     if (mod(ctx%runtime%YR,4)==0) then
       ctx%runtime%YEARLEN=366
       if (mod(ctx%runtime%YR,100)==0) then
         ctx%runtime%YEARLEN=365
         if (mod(ctx%runtime%YR,400)==0) ctx%runtime%YEARLEN=366
       endif
     endif
     if (ITIME==1 .and. .not. ctx%runtime%RESTART_FLAG) then
       ctx%setup%EAHXY=(ctx%input%P8W(:,:)*ctx%input%QV_CURR(:,:))/(0.622+ctx%input%QV_CURR(:,:))
       ctx%setup%TAHXY=ctx%input%T_PHY(:,:)
       ctx%setup%CHXY=0.1
       ctx%setup%CMXY=0.1
     endif
     if (ITIME>0) then
       where(ctx%grid%XICE>=ctx%runtime%XICE_THRESHOLD) ctx%input%LAI=0.01
       ctx%runtime%ITIMESTEP=ITIME
       ctx%runtime%SOIL_UPDATE_STEPS=max(nint(ctx%runtime%SOILTSTEP/ctx%runtime%DTBL),1)
       ctx%runtime%CALCULATE_SOIL=mod(ITIME,ctx%runtime%SOIL_UPDATE_STEPS)==0
       ctx%grid%ZSOIL(1)=-ctx%runtime%DZS(1)
       do K=2,ctx%runtime%NSOIL
         ctx%grid%ZSOIL(K)=-ctx%runtime%DZS(K)+ctx%grid%ZSOIL(K-1)
       enddo
       first_columns=allocated(ctx%setup)
#ifdef NOAHMP_OPENACC
       if (.not. first_columns) call NoahmpGpuBeforeRefresh(noahmp)
#endif
       ! Retain initial inactive-cell values for the unchanged file schema.
       if (first_columns) then
         do J=ctx%grid%JTS,ctx%grid%JTE
           do I=ctx%grid%ITS,ctx%grid%ITE
             if ((ctx%grid%XLAND(I,J)-1.5)>=0.0) then
               ctx%setup%SMSTAV(I,J)=1.0
               ctx%setup%SMSTOT(I,J)=1.0
               ctx%setup%SMOIS(I,1:ctx%runtime%NSOIL,J)=1.0
               ctx%setup%TSLB(I,1:ctx%runtime%NSOIL,J)=273.16
             elseif (ctx%grid%XICE(I,J)==1.0) then
               ctx%setup%SMSTAV(I,J)=1.0
               ctx%setup%SMSTOT(I,J)=1.0
               ctx%setup%SMOIS(I,1:ctx%runtime%NSOIL,J)=1.0
             endif
             if (ctx%grid%XICE(I,J)>=ctx%runtime%XICE_THRESHOLD) then
               ctx%setup%SH2O(I,1:ctx%runtime%NSOIL,J)=1.0
               ctx%input%LAI(I,J)=0.01
             endif
           enddo
         enddo
       endif
       p=0
       do J=ctx%grid%JTS,ctx%grid%JTE
         do I=ctx%grid%ITS,ctx%grid%ITE
           p=p+1
           if (.not. ctx%active(p)) cycle
           if (first_columns) then
             call HrldasInitializeColumn(noahmp(p),ctx,I,J)
           else
             call NoahmpResetStepTemporaries(noahmp(p))
             call HrldasRefreshColumn(noahmp(p),ctx,I,J)
           endif
         enddo
       enddo
       if (first_columns) call HrldasReleaseSetup(ctx)
       reset_soil=ctx%runtime%SOIL_UPDATE_STEPS==1
       if (ctx%runtime%SOIL_UPDATE_STEPS>1) reset_soil=mod(ITIME,ctx%runtime%SOIL_UPDATE_STEPS)==1
       ctx%runtime%sflx_count_sum=0
       call system_clock(count=ctx%runtime%count_before_sflx,count_rate=ctx%runtime%clock_rate)
       call NoahmpDriverMain(noahmp,ctx%active,ctx%runtime%block_size,reset_soil,workspaces,fatal_error)
       if (fatal_error) then
#ifdef MPP_LAND
         call mpp_land_abort()
#endif
         error stop 'Noah-MP fatal column diagnostic; see report above'
       endif
       call system_clock(count=ctx%runtime%count_after_sflx,count_rate=ctx%runtime%clock_rate)
       ctx%runtime%sflx_count_sum=ctx%runtime%count_after_sflx-ctx%runtime%count_before_sflx
     endif
! Output for history
    OUTPUT_FOR_HISTORY: if (ctx%runtime%output_timestep > 0) then
        if (mod(ITIME*ctx%runtime%noah_timestep, ctx%runtime%output_timestep) == 0 &
            .and. .not. ctx%runtime%skip_first_output ) then

#ifdef NOAHMP_OPENACC
           call NoahmpGpuSyncHost(noahmp)
#endif
           call prepare_output_file (trim(ctx%runtime%outdir), version,                                   &
                                     ctx%grid%igrid, ctx%runtime%output_timestep, ctx%grid%llanduse,      &
                                     ctx%runtime%split_output_count, ctx%grid%hgrid,                      &
                                     ctx%grid%ixfull, ctx%grid%jxfull, ctx%grid%ixpar, ctx%grid%jxpar, &
                                     ctx%grid%xstartpar, ctx%grid%ystartpar,                           &
                                     ctx%grid%iswater, ctx%grid%mapproj, ctx%grid%lat1, ctx%grid%lon1, &
                                     ctx%grid%dx, ctx%grid%dy, ctx%grid%truelat1, ctx%grid%truelat2,   &
                                     ctx%grid%cen_lon, ctx%runtime%NSOIL, ctx%runtime%nsnow, ctx%runtime%NUMRAD,&
                                     ctx%runtime%dzs, ctx%runtime%startdate, ctx%runtime%olddate,               &
                                     ctx%runtime%spinup_loop,                                             &
                                     ctx%runtime%spinup_loops, ctx%grid%IVGTYP, ctx%grid%ISLTYP)

           DEFINE_MODE_LOOP : do imode = 1, 2

              call set_output_define_mode(imode)

              ! For 3D arrays, we need to know whether the Z dimension is snow layers, or soil layers.

            ! Properties - Assigned or predicted
              call HrldasGatherField(ctx,noahmp,'ivgtyp')
              call add_to_output(ctx%output_int, "IVGTYP"  , "Dominant vegetation category"         , "category"              )
              call HrldasGatherField(ctx,noahmp,'isltyp')
              call add_to_output(ctx%output_int, "ISLTYP"  , "Dominant soil category"               , "category"              )
              call HrldasGatherField(ctx,noahmp,'fvegxy')
              call add_to_output(ctx%output2, "FVEG"    , "Green Vegetation Fraction"            , "-"                     )
              call HrldasGatherField(ctx,noahmp,'lai')
              call add_to_output(ctx%output2, "LAI"     , "Leaf area index"                      , "m2/m2"                 )
              call HrldasGatherField(ctx,noahmp,'xsaixy')
              call add_to_output(ctx%output2, "SAI"     , "Stem area index"                      , "m2/m2"                 )
            ! Forcing
              call HrldasGatherField(ctx,noahmp,'swdown')
              call add_to_output(ctx%output2, "SWFORC"  , "Shortwave forcing"                    , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'coszen')
              call add_to_output(ctx%output2, "COSZ"    , "Cosine of zenith angle"               , "-"                     )
              call HrldasGatherField(ctx,noahmp,'glw')
              call add_to_output(ctx%output2, "LWFORC"  , "Longwave forcing"                     , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'rainbl')
              call add_to_output(ctx%output2, "RAINRATE", "Precipitation rate"                   , "mm/timestep"           )
            ! Grid energy budget terms
              call HrldasGatherField(ctx,noahmp,'emiss')
              call add_to_output(ctx%output2, "EMISS"   , "Grid emissivity"                      , "-"                     )
              call HrldasGatherField(ctx,noahmp,'fsaxy')
              call add_to_output(ctx%output2, "FSA"     , "Total absorbed SW radiation"          , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'firaxy')
              call add_to_output(ctx%output2, "FIRA"    , "Total net LW radiation to atmosphere" , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'grdflx')
              call add_to_output(ctx%output2, "GRDFLX"  , "Heat flux into the soil"              , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'hfx')
              call add_to_output(ctx%output2, "HFX"     , "Total sensible heat to atmosphere"    , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'lh')
              call add_to_output(ctx%output2, "LH"      , "Total latent heat to atmosphere"      , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'ecanxy')
              call add_to_output(ctx%output2, "ECAN"    , "Canopy water evaporation rate"        , "mm/s"                  )
              call HrldasGatherField(ctx,noahmp,'etranxy')
              call add_to_output(ctx%output2, "ETRAN"   , "Transpiration rate"                   , "mm/s"                  )
              call HrldasGatherField(ctx,noahmp,'edirxy')
              call add_to_output(ctx%output2, "EDIR"    , "Direct from soil evaporation rate"    , "mm/s"                  )
              call HrldasGatherField(ctx,noahmp,'albedo')
              call add_to_output(ctx%output2, "ALBEDO"  , "Surface albedo"                       , "-"                     )
            ! Grid water budget terms - in addition to above
              call HrldasGatherField(ctx,noahmp,'udrunoff')
              call add_to_output(ctx%output2, "UGDRNOFF", "Accumulated underground runoff"       , "mm"                    )
              call HrldasGatherField(ctx,noahmp,'sfcrunoff')
              call add_to_output(ctx%output2, "SFCRNOFF", "Accumulatetd surface runoff"          , "mm"                    )
              call HrldasGatherField(ctx,noahmp,'canliqxy')
              call add_to_output(ctx%output2, "CANLIQ"  , "Canopy liquid water content"          , "mm"                    )
              call HrldasGatherField(ctx,noahmp,'canicexy')
              call add_to_output(ctx%output2, "CANICE"  , "Canopy ice water content"             , "mm"                    )
              call HrldasGatherField(ctx,noahmp,'zwtxy')
              call add_to_output(ctx%output2, "ZWT"     , "Depth to water table"                 , "m"                     )
              call HrldasGatherField(ctx,noahmp,'waxy')
              call add_to_output(ctx%output2, "WA"      , "Water in aquifer"                     , "mm"                    )
              call HrldasGatherField(ctx,noahmp,'wtxy')
              call add_to_output(ctx%output2, "WT"      , "Water in aquifer and saturated soil"  , "mm"                    )
              call HrldasGatherField(ctx,noahmp,'qtdrain')
              call add_to_output(ctx%output2, "QTDRAIN" , "Accumulated tile drainage"            , "mm"                    )
              call HrldasGatherField(ctx,noahmp,'pondingxy')
              call add_to_output(ctx%output2, "PONDING" , "total surface ponding per time step"            , "mm/s"        )
            ! Additional needed to close the canopy energy budget
              call HrldasGatherField(ctx,noahmp,'savxy')
              call add_to_output(ctx%output2, "SAV"     , "Solar radiation absorbed by canopy"   , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'trxy')
              call add_to_output(ctx%output2, "TR"      , "Transpiration heat"                   , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'evcxy')
              call add_to_output(ctx%output2, "EVC"     , "Canopy evap heat"                     , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'ircxy')
              call add_to_output(ctx%output2, "IRC"     , "Canopy net LW rad"                    , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'shcxy')
              call add_to_output(ctx%output2, "SHC"     , "Canopy sensible heat"                 , "W/m2"                  )
            ! Additional needed to close the under canopy ground energy budget
              call HrldasGatherField(ctx,noahmp,'irgxy')
              call add_to_output(ctx%output2, "IRG"     , "below-canopy ground net LW rad"       , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'shgxy')
              call add_to_output(ctx%output2, "SHG"     , "below-canopy ground sensible heat"    , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'evgxy')
              call add_to_output(ctx%output2, "EVG"     , "below-canopy ground evap heat"        , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'ghvxy')
              call add_to_output(ctx%output2, "GHV"     , "below-canopy ground heat to soil"     , "W/m2"                  )
            ! Needed to close the bare ground energy budget
              call HrldasGatherField(ctx,noahmp,'sagxy')
              call add_to_output(ctx%output2, "SAG"     , "Solar radiation absorbed by ground"   , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'irbxy')
              call add_to_output(ctx%output2, "IRB"     , "Net LW rad to atm bare ground"        , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'shbxy')
              call add_to_output(ctx%output2, "SHB"     , "Sensible heat to atm bare ground"     , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'evbxy')
              call add_to_output(ctx%output2, "EVB"     , "Evaporation heat to atm bare ground"  , "W/m2"                  )
              call HrldasGatherField(ctx,noahmp,'ghbxy')
              call add_to_output(ctx%output2, "GHB"     , "Ground heat flux to soil bare ground" , "W/m2"                  )
            ! Above-soil temperatures
              call HrldasGatherField(ctx,noahmp,'tradxy')
              call add_to_output(ctx%output2, "TRAD"    , "Surface radiative temperature"        , "K"                     )
              call HrldasGatherField(ctx,noahmp,'tgxy')
              call add_to_output(ctx%output2, "TG"      , "Ground temperature"                   , "K"                     )
              call HrldasGatherField(ctx,noahmp,'tvxy')
              call add_to_output(ctx%output2, "TV"      , "Vegetation temperature"               , "K"                     )
              call HrldasGatherField(ctx,noahmp,'tahxy')
              call add_to_output(ctx%output2, "TAH"     , "Canopy air temperature"               , "K"                     )
              call HrldasGatherField(ctx,noahmp,'tgvxy')
              call add_to_output(ctx%output2, "TGV"     , "Ground surface Temp vegetated"        , "K"                     )
              call HrldasGatherField(ctx,noahmp,'tgbxy')
              call add_to_output(ctx%output2, "TGB"     , "Ground surface Temp bare"             , "K"                     )
              call HrldasGatherField(ctx,noahmp,'t2mvxy')
              call add_to_output(ctx%output2, "T2MV"    , "2m Air Temp vegetated"                , "K"                     )
              call HrldasGatherField(ctx,noahmp,'t2mbxy')
              call add_to_output(ctx%output2, "T2MB"    , "2m Air Temp bare"                     , "K"                     )
            ! Above-soil moisture
              call HrldasGatherField(ctx,noahmp,'q2mvxy')
              call add_to_output(ctx%output2, "Q2MV"    , "2m mixing ratio vegetated"            , "kg/kg"                 )
              call HrldasGatherField(ctx,noahmp,'q2mbxy')
              call add_to_output(ctx%output2, "Q2MB"    , "2m mixing ratio bare"                 , "kg/kg"                 )
              call HrldasGatherField(ctx,noahmp,'eahxy')
              call add_to_output(ctx%output2, "EAH"     , "Canopy air vapor pressure"            , "Pa"                    )
              call HrldasGatherField(ctx,noahmp,'fwetxy')
              call add_to_output(ctx%output2, "FWET"    , "Wetted fraction of canopy"            , "fraction"              )
            ! Snow and soil - 3D terms
              call HrldasGatherField(ctx,noahmp,'zsnsoxy')
              ctx%output_snow=ctx%output_snso(:,1:ctx%runtime%nsnow,:)
              call add_to_output(ctx%output_snow,"ZSNSO_SN","Snow layer depth from snow surface","m","SNOW")
              call HrldasGatherField(ctx,noahmp,'snicexy')
              call add_to_output(ctx%output_snow, "SNICE"   , "Snow layer ice"                       , "mm"             , &
      & "SNOW")
              call HrldasGatherField(ctx,noahmp,'snliqxy')
              call add_to_output(ctx%output_snow, "SNLIQ"   , "Snow layer liquid water"              , "mm"             , &
      & "SNOW")
              call HrldasGatherField(ctx,noahmp,'tslb')
              call add_to_output(ctx%output_soil, "SOIL_T"  , "soil temperature"                     , "K"              , &
      & "SOIL")
              call HrldasGatherField(ctx,noahmp,'smois')
              call add_to_output(ctx%output_soil, "SOIL_M"  , "volumetric soil moisture"             , "m3/m3"          , &
      & "SOIL")
              call HrldasGatherField(ctx,noahmp,'sh2o')
              call add_to_output(ctx%output_soil, "SOIL_W"  , "liquid volumetric soil moisture"      , "m3/m3"          , &
      & "SOIL")
              call HrldasGatherField(ctx,noahmp,'tsnoxy')
              call add_to_output(ctx%output_snow, "SNOW_T"  , "snow temperature"                     , "K"              , &
      & "SNOW")
              call HrldasGatherField(ctx,noahmp,'albsnowdirxy')
              call add_to_output(ctx%output_rad, "ALBSNOWDIR" , "Snow albedo (direct)"             , "-"              , "RADN")
              call HrldasGatherField(ctx,noahmp,'albsnowdifxy')
              call add_to_output(ctx%output_rad, "ALBSNOWDIF" , "Snow albedo (diffuse)"            , "-"              , "RADN")
              call HrldasGatherField(ctx,noahmp,'albsfcdirxy')
              call add_to_output(ctx%output_rad, "ALBSFCDIR"  , "Surface albedo (direct)"          , "-"              , "RADN")
              call HrldasGatherField(ctx,noahmp,'albsfcdifxy')
              call add_to_output(ctx%output_rad, "ALBSFCDIF"  , "Surface albedo (diffuse)"         , "-"              , "RADN")
              call HrldasGatherField(ctx,noahmp,'albsoildirxy')
              call add_to_output(ctx%output_rad, "ALBSOILDIR"  , "Soil albedo (direct)"            , "-"              , "RADN")
              call HrldasGatherField(ctx,noahmp,'albsoildifxy')
              call add_to_output(ctx%output_rad, "ALBSOILDIF"  , "Soil albedo (diffuse)"           , "-"              , "RADN")
            ! Snow - 2D terms
              call HrldasGatherField(ctx,noahmp,'snowh')
              call add_to_output(ctx%output2, "SNOWH"   , "Snow depth"                           , "m"                     )
              call HrldasGatherField(ctx,noahmp,'snow')
              call add_to_output(ctx%output2, "SNEQV"   , "Snow water equivalent"                , "mm"                    )
              call HrldasGatherField(ctx,noahmp,'qsnowxy')
              call add_to_output(ctx%output2, "QSNOW"   , "Snowfall rate on the ground"          , "mm/s"                  )
              call HrldasGatherField(ctx,noahmp,'qrainxy')
              call add_to_output(ctx%output2, "QRAIN"   , "Rainfall rate on the ground"          , "mm/s"                  )
              call HrldasGatherField(ctx,noahmp,'isnowxy')
              call add_to_output(ctx%output_int, "ISNOW"   , "Number of snow layers"                , "-"                     )
              call HrldasGatherField(ctx,noahmp,'snowc')
              call add_to_output(ctx%output2, "FSNO"    , "Snow-cover fraction on the ground"    , "-"                     )
              call HrldasGatherField(ctx,noahmp,'acsnow')
              call add_to_output(ctx%output2, "ACSNOW"  , "accumulated snow fall"                , "mm"                    )
              call HrldasGatherField(ctx,noahmp,'acsnom')
              call add_to_output(ctx%output2, "ACSNOM"  , "accumulated snow melt water"          , "mm"                    )
              call HrldasGatherField(ctx,noahmp,'qsnbotxy')
              call add_to_output(ctx%output2, "QSNBOT"  , "water (melt+rain through) out of snow bottom"   , "mm/s"        )
              call HrldasGatherField(ctx,noahmp,'qmeltxy')
              call add_to_output(ctx%output2, "QMELT"   , "snow melt due to phase change"                  , "mm/s"        )
            ! SNICAR snow albedo scheme
            ! Exchange coefficients
              call HrldasGatherField(ctx,noahmp,'cmxy')
              call add_to_output(ctx%output2, "CM"      , "Momentum drag coefficient"            , "m/s"                   )
              call HrldasGatherField(ctx,noahmp,'chxy')
              call add_to_output(ctx%output2, "CH"      , "Sensible heat exchange coefficient"   , "m/s"                   )
              call HrldasGatherField(ctx,noahmp,'chvxy')
              call add_to_output(ctx%output2, "CHV"     , "Exchange coefficient vegetated"       , "m/s"                   )
              call HrldasGatherField(ctx,noahmp,'chbxy')
              call add_to_output(ctx%output2, "CHB"     , "Exchange coefficient bare"            , "m/s"                   )
              call HrldasGatherField(ctx,noahmp,'chleafxy')
              call add_to_output(ctx%output2, "CHLEAF"  , "Exchange coefficient leaf"            , "m/s"                   )
              call HrldasGatherField(ctx,noahmp,'chucxy')
              call add_to_output(ctx%output2, "CHUC"    , "Exchange coefficient bare"            , "m/s"                   )
              call HrldasGatherField(ctx,noahmp,'chv2xy')
              call add_to_output(ctx%output2, "CHV2"    , "Exchange coefficient 2-m vegetated"   , "m/s"                   )
              call HrldasGatherField(ctx,noahmp,'chb2xy')
              call add_to_output(ctx%output2, "CHB2"    , "Exchange coefficient 2-m bare"        , "m/s"                   )
            ! Carbon allocation model
              call HrldasGatherField(ctx,noahmp,'lfmassxy')
              call add_to_output(ctx%output2, "LFMASS"  , "Leaf mass"                            , "g/m2"                  )
              call HrldasGatherField(ctx,noahmp,'rtmassxy')
              call add_to_output(ctx%output2, "RTMASS"  , "Mass of fine roots"                   , "g/m2"                  )
              call HrldasGatherField(ctx,noahmp,'stmassxy')
              call add_to_output(ctx%output2, "STMASS"  , "Stem mass"                            , "g/m2"                  )
              call HrldasGatherField(ctx,noahmp,'woodxy')
              call add_to_output(ctx%output2, "WOOD"    , "Mass of wood and woody roots"         , "g/m2"                  )
              call HrldasGatherField(ctx,noahmp,'grainxy')
              call add_to_output(ctx%output2, "GRAIN"   , "Mass of grain"                        , "g/m2"                  )
              call HrldasGatherField(ctx,noahmp,'gddxy')
              call add_to_output(ctx%output2, "GDD"     , "Growing degree days "                 , "-"                     )
              call HrldasGatherField(ctx,noahmp,'stblcpxy')
              call add_to_output(ctx%output2, "STBLCP"  , "Stable carbon in deep soil"           , "gC/m2"                 )
              call HrldasGatherField(ctx,noahmp,'fastcpxy')
              call add_to_output(ctx%output2, "FASTCP"  , "Short-lived carbon in shallow soil"   , "gC/m2"                 )
              call HrldasGatherField(ctx,noahmp,'neexy')
              call add_to_output(ctx%output2, "NEE"     , "Net ecosystem exchange"               , "gCO2/m2/s"             )
              call HrldasGatherField(ctx,noahmp,'gppxy')
              call add_to_output(ctx%output2, "GPP"     , "Net instantaneous assimilation"       , "gC/m2/s"               )
              call HrldasGatherField(ctx,noahmp,'nppxy')
              call add_to_output(ctx%output2, "NPP"     , "Net primary productivity"             , "gC/m2/s"               )
              call HrldasGatherField(ctx,noahmp,'psnxy')
              call add_to_output(ctx%output2, "PSN"     , "Total photosynthesis"                 , "umol CO2/m2/s"         )
              call HrldasGatherField(ctx,noahmp,'aparxy')
              call add_to_output(ctx%output2, "APAR"    , "Photosynthesis active energy by canopy", "W/m2"                 )
            ! additional NoahMP output
            if (ctx%runtime%noahmp_output > 0) then
              ! additional water budget terms
              call HrldasGatherField(ctx,noahmp,'qintsxy')
              call add_to_output(ctx%output2, "QINTS"   , "canopy interception (loading) rate for snowfall", "mm/s"        )
              call HrldasGatherField(ctx,noahmp,'qintrxy')
              call add_to_output(ctx%output2, "QINTR"   , "canopy interception rate for rain"              , "mm/s"        )
              call HrldasGatherField(ctx,noahmp,'qdripsxy')
              call add_to_output(ctx%output2, "QDRIPS"  , "drip (unloading) rate for intercepted snow"     , "mm/s"        )
              call HrldasGatherField(ctx,noahmp,'qdriprxy')
              call add_to_output(ctx%output2, "QDRIPR"  , "drip rate for canopy intercepted rain"          , "mm/s"        )
              call HrldasGatherField(ctx,noahmp,'qthrosxy')
              call add_to_output(ctx%output2, "QTHROS"  , "throughfall of snowfall"                        , "mm/s"        )
              call HrldasGatherField(ctx,noahmp,'qthrorxy')
              call add_to_output(ctx%output2, "QTHROR"  , "throughfall for rain"                           , "mm/s"        )
              call HrldasGatherField(ctx,noahmp,'qsnsubxy')
              call add_to_output(ctx%output2, "QSNSUB"  , "snow surface sublimation rate"                  , "mm/s"        )
              call HrldasGatherField(ctx,noahmp,'qsnfroxy')
              call add_to_output(ctx%output2, "QSNFRO"  , "snow surface frost rate"                        , "mm/s"        )
              call HrldasGatherField(ctx,noahmp,'qsubcxy')
              call add_to_output(ctx%output2, "QSUBC"   , "canopy snow sublimation rate"                   , "mm/s"        )
              call HrldasGatherField(ctx,noahmp,'qfrocxy')
              call add_to_output(ctx%output2, "QFROC"   , "canopy snow frost rate"                         , "mm/s"        )
              call HrldasGatherField(ctx,noahmp,'qevacxy')
              call add_to_output(ctx%output2, "QEVAC"   , "canopy snow evaporation rate"                   , "mm/s"        )
              call HrldasGatherField(ctx,noahmp,'qdewcxy')
              call add_to_output(ctx%output2, "QDEWC"   , "canopy snow dew rate"                           , "mm/s"        )
              call HrldasGatherField(ctx,noahmp,'qfrzcxy')
              call add_to_output(ctx%output2, "QFRZC"   , "refreezing rate of canopy liquid water"         , "mm/s"        )
              call HrldasGatherField(ctx,noahmp,'qmeltcxy')
              call add_to_output(ctx%output2, "QMELTC"  , "melting rate of canopy snow"                    , "mm/s"        )
              call HrldasGatherField(ctx,noahmp,'fpicexy')
              call add_to_output(ctx%output2, "FPICE"   , "snow fraction in precipitation"                 , "-"           )
              call HrldasGatherField(ctx,noahmp,'acc_qinsurxy')
              call add_to_output(ctx%output2,"ACC_QINSUR", "accumuated water flux to soil within soil timestep"     , &
      & "m/s*dt_soil/dt_main")
              call HrldasGatherField(ctx,noahmp,'acc_qsevaxy')
              call add_to_output(ctx%output2,"ACC_QSEVA" , "accumulated soil surface evap rate within soil timestep", &
      & "m/s*dt_soil/dt_main")
              call HrldasGatherField(ctx,noahmp,'acc_etranixy')
              call add_to_output(ctx%output_soil,"ACC_ETRANI", "accumualted transpiration rate within soil timestep"    , &
      & "m/s*dt_soil/dt_main","SOIL")
              call HrldasGatherField(ctx,noahmp,'acc_dwaterxy')
              call add_to_output(ctx%output2,"ACC_DWATER", "accumulated water storage change within soil timestep"  , "mm")
              call HrldasGatherField(ctx,noahmp,'acc_prcpxy')
              call add_to_output(ctx%output2,"ACC_PRCP"  , "accumulated precipitation within soil timestep"         , "mm")
              call HrldasGatherField(ctx,noahmp,'acc_ecanxy')
              call add_to_output(ctx%output2,"ACC_ECAN"  , "accumulated net canopy evaporation within soil timestep", "mm")
              call HrldasGatherField(ctx,noahmp,'acc_etranxy')
              call add_to_output(ctx%output2,"ACC_ETRAN" , "accumulated transpiration within soil timestep"         , "mm")
              call HrldasGatherField(ctx,noahmp,'acc_edirxy')
              call add_to_output(ctx%output2,"ACC_EDIR"  , "accumulated net ground evaporation within soil timestep", "mm")
              call HrldasGatherField(ctx,noahmp,'acc_glaflwxy')
              call add_to_output(ctx%output2,"ACC_GLAFLW", "accumuated glacier excessive flow per soil timestep"    , "mm")
              ! additional energy terms
              call HrldasGatherField(ctx,noahmp,'pahxy')
              call add_to_output(ctx%output2, "PAH"     , "Precipitation advected heat flux"                         , &
      & "W/m2"    )
              call HrldasGatherField(ctx,noahmp,'pahgxy')
              call add_to_output(ctx%output2, "PAHG"    , "Precipitation advected heat flux to below-canopy ground"  , &
      & "W/m2"    )
              call HrldasGatherField(ctx,noahmp,'pahbxy')
              call add_to_output(ctx%output2, "PAHB"    , "Precipitation advected heat flux to bare ground"          , &
      & "W/m2"    )
              call HrldasGatherField(ctx,noahmp,'pahvxy')
              call add_to_output(ctx%output2, "PAHV"    , "Precipitation advected heat flux to canopy"               , &
      & "W/m2"    )
              call HrldasGatherField(ctx,noahmp,'acc_ssoilxy')
              call add_to_output(ctx%output2, "ACC_SSOIL","accumulated heat flux into snow/soil within soil timestep", &
      & "W/m2"    )
              call HrldasGatherField(ctx,noahmp,'eflxbxy')
              call add_to_output(ctx%output2, "EFLXB"   , "accumulated heat flux through soil bottom"                , &
      & "J/m2"    )
              call HrldasGatherField(ctx,noahmp,'soilenergy')
              call add_to_output(ctx%output2, "SOILENERGY","energy content in soil relative to 273.16"               , &
      & "KJ/m2"   )
              call HrldasGatherField(ctx,noahmp,'snowenergy')
              call add_to_output(ctx%output2, "SNOWENERGY","energy content in snow relative to 273.16"               , &
      & "KJ/m2"   )
              call HrldasGatherField(ctx,noahmp,'canhsxy')
              call add_to_output(ctx%output2, "CANHS"   , "canopy heat storage change"                               , &
      & "W/m2"    )
              ! additional forcing terms
              call HrldasGatherField(ctx,noahmp,'rainlsm')
              call add_to_output(ctx%output2, "RAINLSM" , "lowest model liquid precipitation into LSM"               , &
      & "mm/s"    )
              call HrldasGatherField(ctx,noahmp,'snowlsm')
              call add_to_output(ctx%output2, "SNOWLSM" , "lowest model snowfall into LSM"                           , &
      & "mm/s"    )
              call HrldasGatherField(ctx,noahmp,'forctlsm')
              call add_to_output(ctx%output2, "FORCTLSM", "lowest model temperature into LSM"                        , &
      & "K"       )
              call HrldasGatherField(ctx,noahmp,'forcqlsm')
              call add_to_output(ctx%output2, "FORCQLSM", "lowest model specific humidty into LSM"                   , &
      & "kg/kg"   )
              call HrldasGatherField(ctx,noahmp,'forcplsm')
              call add_to_output(ctx%output2, "FORCPLSM", "lowest model pressure into LSM"                           , &
      & "Pa"      )
              call HrldasGatherField(ctx,noahmp,'forczlsm')
              call add_to_output(ctx%output2, "FORCZLSM", "lowest model forcing height into LSM"                     , &
      & "m"       )
              call HrldasGatherField(ctx,noahmp,'forcwlsm')
              call add_to_output(ctx%output2, "FORCWLSM", "lowest model wind speed into LSM"                         , &
      & "m/s"     )
              call HrldasGatherField(ctx,noahmp,'radswvisfrac')
              call add_to_output(ctx%output2, "SWVISFRAC", "Fraction of visible band downward solar radiation", "-"      )
              call HrldasGatherField(ctx,noahmp,'radswdirfrac')
              call add_to_output(ctx%output2, "SWDIRFRAC", "Fraction of downward solar direct radiation",   "-"          )
            endif

            ! Irrigation
            if ( ctx%runtime%IOPT_IRR > 0 ) then
              call HrldasGatherField(ctx,noahmp,'irnumsi')
              call add_to_output(ctx%output_int, "IRNUMSI" , "Sprinkler irrigation count"                               , &
      & "-"       )
              call HrldasGatherField(ctx,noahmp,'irnummi')
              call add_to_output(ctx%output_int, "IRNUMMI" , "Micro irrigation count"                                   , &
      & "-"       )
              call HrldasGatherField(ctx,noahmp,'irnumfi')
              call add_to_output(ctx%output_int, "IRNUMFI" , "Flood irrigation count"                                   , &
      & "-"       )
              call HrldasGatherField(ctx,noahmp,'ireloss')
              call add_to_output(ctx%output2, "IRELOSS" , "Accumulated sprinkler Evaporation"                        , &
      & "mm"      )
              call HrldasGatherField(ctx,noahmp,'irsivol')
              call add_to_output(ctx%output2, "IRSIVOL" , "Sprinkler irrigation amount"                              , &
      & "mm"      )
              call HrldasGatherField(ctx,noahmp,'irmivol')
              call add_to_output(ctx%output2, "IRMIVOL" , "Micro irrigation amount"                                  , &
      & "mm"      )
              call HrldasGatherField(ctx,noahmp,'irfivol')
              call add_to_output(ctx%output2, "IRFIVOL" , "Flood irrigation amount"                                  , &
      & "mm"      )
              call HrldasGatherField(ctx,noahmp,'irrsplh')
              call add_to_output(ctx%output2, "IRRSPLH" , "Accumulated latent heating due to sprinkler"              , &
      & "J/m2"    )
            endif
            ! MMF groundwater  model
            ! Wetland model
            if ( ctx%runtime%IOPT_WETLAND > 0 ) then
              call HrldasGatherField(ctx,noahmp,'fsatxy')
              call add_to_output(ctx%output2, "FSAT"    , "saturated fraction of the grid"                           , "-" )
              call HrldasGatherField(ctx,noahmp,'wsurfxy')
              call add_to_output(ctx%output2, "WSURF"   , "Wetland Water Storage"                                    , "mm")
            endif

           ! For now, no urban output variables included

           enddo DEFINE_MODE_LOOP

           call finalize_output_file(ctx%runtime%split_output_count,itime)

        endif

        if(ctx%runtime%skip_first_output) ctx%runtime%skip_first_output = .false.

     endif OUTPUT_FOR_HISTORY

     call HrldasConsoleValues(ctx,noahmp,console_soil_temperature,console_leaf_area)
     if (ctx%grid%IVGTYP(ctx%grid%xstart,ctx%grid%ystart)==ctx%grid%ISWATER) then
       write(*,'(" ***DATE=", A19)', advance="NO") ctx%runtime%olddate
     else
       write(*,'(" ***DATE=", A19, 6F10.5)', advance="NO")     &
       ctx%runtime%olddate, console_soil_temperature, &
       console_leaf_area
     endif

!------------------------------------------------------------------------
! Write Restart - timestamp equal to output will have same states
!------------------------------------------------------------------------

     if (ctx%runtime%restart_frequency_hours > 0) then
        if (mod(ITIME, int(ctx%runtime%restart_frequency_hours*3600./nint(ctx%runtime%dtbl))) == 0) then
           call lsm_restart()
        endif
     endif
     if (ctx%runtime%restart_frequency_hours <= 0) then
        if ((ctx%runtime%olddate( 9:10) == "01") .and. (ctx%runtime%olddate(12:13) == "00") .and. &
            (ctx%runtime%olddate(15:16) == "00") .and. (ctx%runtime%olddate(18:19) == "00")) then
            call lsm_restart()  ! jlm - i moved all the restart code to a subroutine.
        endif
     endif

!------------------------------------------------------------------------
! Advance the time
!------------------------------------------------------------------------

     call geth_newdate(ctx%runtime%newdate, ctx%runtime%olddate, nint(ctx%runtime%dtbl))
     ctx%runtime%olddate = ctx%runtime%newdate

     if(itime > 0 .and. ctx%runtime%spinup_loops > 0 .and. mod(itime,ctx%runtime%ntime/(ctx%runtime%spinup_loops+1)) == 0) then
       ctx%runtime%spinup_loop = ctx%runtime%spinup_loop + 1
       call geth_newdate(ctx%runtime%olddate, ctx%runtime%startdate, nint(ctx%runtime%dtbl))
       ctx%runtime%reset_spinup_date = .true.
       ctx%runtime%reset_spinup_datea = .true.
     endif

! update the timer
     call system_clock(count=ctx%runtime%clock_count_2, count_rate=ctx%runtime%clock_rate)
     ctx%runtime%timing_sum = ctx%runtime%timing_sum + &
                        float(ctx%runtime%clock_count_2-ctx%runtime%clock_count_1)/float(ctx%runtime%clock_rate)
     write(*,'("    Timing: ",f6.2," Cumulative:  ", f10.2, "  SFLX: ", f6.2 )') &
          float(ctx%runtime%clock_count_2-ctx%runtime%clock_count_1)/float(ctx%runtime%clock_rate), &
          ctx%runtime%timing_sum, real(ctx%runtime%sflx_count_sum) / real(ctx%runtime%clock_rate)

     ctx%runtime%clock_count_1 = ctx%runtime%clock_count_2


end subroutine land_driver_exe

!!===============================================================================
subroutine land_driver_finish()
#ifdef NOAHMP_OPENACC
  call NoahmpGpuFinish(noahmp,workspaces,ctx%active)
#endif
end subroutine land_driver_finish

subroutine lsm_restart()
  implicit none

  print*, 'Write restart at '//ctx%runtime%olddate(1:13)

#ifdef NOAHMP_OPENACC
  call NoahmpGpuSyncHost(noahmp)
#endif
  call prepare_restart_file (trim(ctx%runtime%outdir), version, ctx%grid%igrid, ctx%grid%llanduse, &
              ctx%runtime%olddate, ctx%runtime%startdate, ctx%grid%ixfull, ctx%grid%jxfull,  &
              ctx%grid%ixpar, ctx%grid%jxpar, ctx%grid%xstartpar, ctx%grid%ystartpar,  &
              ctx%runtime%NSOIL, ctx%runtime%nsnow, ctx%runtime%NUMRAD,                         &
              ctx%grid%dx, ctx%grid%dy, ctx%grid%truelat1, ctx%grid%truelat2,          &
              ctx%grid%mapproj, ctx%grid%lat1, ctx%grid%lon1, ctx%grid%cen_lon,        &
              ctx%grid%iswater, ctx%grid%ivgtyp,               &
              ctx%runtime%num_urban_hi,  ctx%runtime%urban_map_zrd,  &
              ctx%runtime%urban_map_zwd, ctx%runtime%urban_map_gd,   &
              ctx%runtime%urban_map_zd,  ctx%runtime%urban_map_zdf,  &
              ctx%runtime%urban_map_bd,  ctx%runtime%urban_map_wd,   &
              ctx%runtime%urban_map_gbd, ctx%runtime%urban_map_fbd,  &
              ctx%runtime%urban_map_zgrd,ctx%runtime%num_urban_ndm,  &
              ctx%runtime%SF_URBAN_PHYSICS)

  call HrldasGatherField(ctx,noahmp,'tslb')

  call add_to_restart(ctx%output_soil, "SOIL_T", layers="SOIL")
  call HrldasGatherField(ctx,noahmp,'tsnoxy')
  call add_to_restart(ctx%output_snow, "SNOW_T", layers="SNOW")
  call HrldasGatherField(ctx,noahmp,'smois')
  call add_to_restart(ctx%output_soil, "SMC"   , layers="SOIL")
  call HrldasGatherField(ctx,noahmp,'sh2o')
  call add_to_restart(ctx%output_soil, "SH2O"  , layers="SOIL")
  call HrldasGatherField(ctx,noahmp,'zsnsoxy')
  call add_to_restart(ctx%output_snso, "ZSNSO" , layers="SOSN")
  call HrldasGatherField(ctx,noahmp,'snicexy')
  call add_to_restart(ctx%output_snow, "SNICE" , layers="SNOW")
  call HrldasGatherField(ctx,noahmp,'snliqxy')
  call add_to_restart(ctx%output_snow, "SNLIQ" , layers="SNOW")
  call HrldasGatherField(ctx,noahmp,'fwetxy')
  call add_to_restart(ctx%output2, "FWET"  )
  call HrldasGatherField(ctx,noahmp,'sneqvoxy')
  call add_to_restart(ctx%output2, "SNEQVO")
  call HrldasGatherField(ctx,noahmp,'eahxy')
  call add_to_restart(ctx%output2, "EAH"   )
  call HrldasGatherField(ctx,noahmp,'tahxy')
  call add_to_restart(ctx%output2, "TAH"   )
  call HrldasGatherField(ctx,noahmp,'alboldxy')
  call add_to_restart(ctx%output2, "ALBOLD")
  call HrldasGatherField(ctx,noahmp,'cmxy')
  call add_to_restart(ctx%output2, "CM"    )
  call HrldasGatherField(ctx,noahmp,'chxy')
  call add_to_restart(ctx%output2, "CH"    )
  call HrldasGatherField(ctx,noahmp,'isnowxy')
  call add_to_restart(ctx%output_int, "ISNOW" )
  call HrldasGatherField(ctx,noahmp,'canliqxy')
  call add_to_restart(ctx%output2, "CANLIQ")
  call HrldasGatherField(ctx,noahmp,'canicexy')
  call add_to_restart(ctx%output2, "CANICE")
  call HrldasGatherField(ctx,noahmp,'snow')
  call add_to_restart(ctx%output2, "SNEQV" )
  call HrldasGatherField(ctx,noahmp,'snowh')
  call add_to_restart(ctx%output2, "SNOWH" )
  call HrldasGatherField(ctx,noahmp,'tvxy')
  call add_to_restart(ctx%output2, "TV"    )
  call HrldasGatherField(ctx,noahmp,'tgxy')
  call add_to_restart(ctx%output2, "TG"    )
  call HrldasGatherField(ctx,noahmp,'zwtxy')
  call add_to_restart(ctx%output2, "ZWT"   )
  call HrldasGatherField(ctx,noahmp,'waxy')
  call add_to_restart(ctx%output2, "WA"    )
  call HrldasGatherField(ctx,noahmp,'wtxy')
  call add_to_restart(ctx%output2, "WT"    )
  call HrldasGatherField(ctx,noahmp,'wslakexy')
  call add_to_restart(ctx%output2, "WSLAKE")
  call HrldasGatherField(ctx,noahmp,'lfmassxy')
  call add_to_restart(ctx%output2, "LFMASS")
  call HrldasGatherField(ctx,noahmp,'rtmassxy')
  call add_to_restart(ctx%output2, "RTMASS")
  call HrldasGatherField(ctx,noahmp,'stmassxy')
  call add_to_restart(ctx%output2, "STMASS")
  call HrldasGatherField(ctx,noahmp,'cropcat')
  call add_to_restart(ctx%output_int, "CROPCAT")
  call HrldasGatherField(ctx,noahmp,'woodxy')
  call add_to_restart(ctx%output2, "WOOD"  )
  call HrldasGatherField(ctx,noahmp,'grainxy')
  call add_to_restart(ctx%output2, "GRAIN" )
  call HrldasGatherField(ctx,noahmp,'gddxy')
  call add_to_restart(ctx%output2, "GDD"   )
  call HrldasGatherField(ctx,noahmp,'stblcpxy')
  call add_to_restart(ctx%output2, "STBLCP")
  call HrldasGatherField(ctx,noahmp,'fastcpxy')
  call add_to_restart(ctx%output2, "FASTCP")
  call HrldasGatherField(ctx,noahmp,'lai')
  call add_to_restart(ctx%output2, "LAI"   )
  call HrldasGatherField(ctx,noahmp,'xsaixy')
  call add_to_restart(ctx%output2, "SAI"   )
  call HrldasGatherField(ctx,noahmp,'vegfra')
  call add_to_restart(ctx%output2, "VEGFRA")
  call HrldasGatherField(ctx,noahmp,'gvfmin')
  call add_to_restart(ctx%output2, "GVFMIN")
  call HrldasGatherField(ctx,noahmp,'gvfmax')
  call add_to_restart(ctx%output2, "GVFMAX")
  call HrldasGatherField(ctx,noahmp,'acsnom')
  call add_to_restart(ctx%output2, "ACMELT")
  call HrldasGatherField(ctx,noahmp,'acsnow')
  call add_to_restart(ctx%output2, "ACSNOW")
  call HrldasGatherField(ctx,noahmp,'taussxy')
  call add_to_restart(ctx%output2, "TAUSS" )
  call HrldasGatherField(ctx,noahmp,'qsfc')
  call add_to_restart(ctx%output2, "QSFC"  )
  call HrldasGatherField(ctx,noahmp,'sfcrunoff')
  call add_to_restart(ctx%output2, "SFCRUNOFF")
  call HrldasGatherField(ctx,noahmp,'udrunoff')
  call add_to_restart(ctx%output2, "UDRUNOFF" )
  call HrldasGatherField(ctx,noahmp,'qtdrain')
  call add_to_restart(ctx%output2, "QTDRAIN"  )
  call HrldasGatherField(ctx,noahmp,'acc_ssoilxy')
  call add_to_restart(ctx%output2,"ACC_SSOIL" )
  call HrldasGatherField(ctx,noahmp,'acc_qinsurxy')
  call add_to_restart(ctx%output2,"ACC_QINSUR")
  call HrldasGatherField(ctx,noahmp,'acc_qsevaxy')
  call add_to_restart(ctx%output2,"ACC_QSEVA" )
  call HrldasGatherField(ctx,noahmp,'acc_etranixy')
  call add_to_restart(ctx%output_soil,"ACC_ETRANI", layers="SOIL")
  call HrldasGatherField(ctx,noahmp,'acc_dwaterxy')
  call add_to_restart(ctx%output2,"ACC_DWATER")
  call HrldasGatherField(ctx,noahmp,'acc_prcpxy')
  call add_to_restart(ctx%output2,"ACC_PRCP"  )
  call HrldasGatherField(ctx,noahmp,'acc_ecanxy')
  call add_to_restart(ctx%output2,"ACC_ECAN"  )
  call HrldasGatherField(ctx,noahmp,'acc_etranxy')
  call add_to_restart(ctx%output2,"ACC_ETRAN" )
  call HrldasGatherField(ctx,noahmp,'acc_edirxy')
  call add_to_restart(ctx%output2,"ACC_EDIR"  )
  call HrldasGatherField(ctx,noahmp,'acc_glaflwxy')
  call add_to_restart(ctx%output2,"ACC_GLAFLW")
  call HrldasGatherField(ctx,noahmp,'albsoildirxy')
  call add_to_restart(ctx%output_rad, "ALBSOILDIR", layers="RADN")
  call HrldasGatherField(ctx,noahmp,'albsoildifxy')
  call add_to_restart(ctx%output_rad, "ALBSOILDIF", layers="RADN")

  ! SNICAR snow albedo scheme

 ! irrigation scheme
  if ( ctx%runtime%IOPT_IRR > 0 ) then
     call HrldasGatherField(ctx,noahmp,'irnumsi')
     call add_to_restart(ctx%output_int, "IRNUMSI")
     call HrldasGatherField(ctx,noahmp,'irnummi')
     call add_to_restart(ctx%output_int, "IRNUMMI")
     call HrldasGatherField(ctx,noahmp,'irnumfi')
     call add_to_restart(ctx%output_int, "IRNUMFI")
     call HrldasGatherField(ctx,noahmp,'irwatsi')
     call add_to_restart(ctx%output2, "IRWATSI")
     call HrldasGatherField(ctx,noahmp,'irwatmi')
     call add_to_restart(ctx%output2, "IRWATMI")
     call HrldasGatherField(ctx,noahmp,'irwatfi')
     call add_to_restart(ctx%output2, "IRWATFI")
     call HrldasGatherField(ctx,noahmp,'irsivol')
     call add_to_restart(ctx%output2, "IRSIVOL")
     call HrldasGatherField(ctx,noahmp,'irmivol')
     call add_to_restart(ctx%output2, "IRMIVOL")
     call HrldasGatherField(ctx,noahmp,'irfivol')
     call add_to_restart(ctx%output2, "IRFIVOL")
     call HrldasGatherField(ctx,noahmp,'ireloss')
     call add_to_restart(ctx%output2, "IRELOSS")
     call HrldasGatherField(ctx,noahmp,'irrsplh')
     call add_to_restart(ctx%output2, "IRRSPLH")
  endif

  ! below for MMF groundwater scheme

  ! for wetland scheme
  if ( ctx%runtime%IOPT_WETLAND > 0 ) then
     call HrldasGatherField(ctx,noahmp,'fsatxy')
     call add_to_restart(ctx%output2, "FSATXY"     )
     call HrldasGatherField(ctx,noahmp,'wsurfxy')
     call add_to_restart(ctx%output2, "WSURFXY"    )
  endif

  ! below for urban model

  call finalize_restart_file()

end subroutine lsm_restart

!------------------------------------------------------------------------------------

SUBROUTINE CALC_DECLIN ( NOWDATE, LATITUDE, LONGITUDE, COSZ, JULIAN, &
                         HRANG_OUT, DECLIN_OUT, GMT_OUT, JULDAY_OUT )

  USE MODULE_DATE_UTILITIES
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------

! !ARGUMENTS:
   CHARACTER(LEN=19),                   INTENT(IN)  :: NOWDATE    ! YYYY-MM-DD_HH:MM:SS
   REAL(kind=kind_noahmp),              INTENT(IN)  :: LATITUDE
   REAL(kind=kind_noahmp),              INTENT(IN)  :: LONGITUDE
   REAL(kind=kind_noahmp),              INTENT(OUT) :: COSZ
   REAL(kind=kind_noahmp),              INTENT(OUT) :: JULIAN
   REAL(kind=kind_noahmp),    OPTIONAL, INTENT(OUT) :: HRANG_OUT
   REAL(kind=kind_noahmp),    OPTIONAL, INTENT(OUT) :: DECLIN_OUT
   REAL(kind=kind_noahmp),    OPTIONAL, INTENT(OUT) :: GMT_OUT
   INTEGER,                   OPTIONAL, INTENT(OUT) :: JULDAY_OUT
   REAL(kind=kind_noahmp)                           :: OBECL
   REAL(kind=kind_noahmp)                           :: SINOB
   REAL(kind=kind_noahmp)                           :: SXLONG
   REAL(kind=kind_noahmp)                           :: ARG
   REAL(kind=kind_noahmp)                           :: TLOCTIM
   REAL(kind=kind_noahmp)                           :: HRANG
   REAL(kind=kind_noahmp)                           :: DECLIN
   REAL(kind=kind_noahmp)                           :: GMT
   INTEGER                                          :: JULDAY
   INTEGER                                          :: IDAY
   INTEGER                                          :: IHOUR
   INTEGER                                          :: IMINUTE
   INTEGER                                          :: ISECOND

   REAL(kind=kind_noahmp), PARAMETER                :: DEGRAD = 3.14159265/180.
   REAL(kind=kind_noahmp), PARAMETER                :: DPD    = 360./365.

   call GETH_IDTS(NOWDATE(1:10), NOWDATE(1:4)//"-01-01", IDAY)
   READ(NOWDATE(12:13), *) IHOUR
   READ(NOWDATE(15:16), *) IMINUTE
   READ(NOWDATE(18:19), *) ISECOND
   GMT = REAL(IHOUR) + IMINUTE/60.0 + ISECOND/3600.0
   JULIAN = REAL(IDAY) + GMT/24.
   JULDAY = IDAY

!
! FOR SHORT WAVE RADIATION

   DECLIN=0.

!-----OBECL : OBLIQUITY = 23.5 DEGREE.

   OBECL=23.5*DEGRAD
   SINOB=SIN(OBECL)

!-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:

   IF(JULIAN >= 80.0) SXLONG = DPD*(JULIAN-80.0)*DEGRAD
   IF(JULIAN <  80.0) SXLONG = DPD*(JULIAN+285.0)*DEGRAD
   ARG    = SINOB * SIN(SXLONG)
   DECLIN = ASIN(ARG)

   TLOCTIM = REAL(IHOUR) + REAL(IMINUTE)/60.0 + REAL(ISECOND)/3600.0 + LONGITUDE/15.0 ! LOCAL TIME IN HOURS
   TLOCTIM = AMOD(TLOCTIM+24.0, 24.0)
   HRANG=15.*(TLOCTIM-12.)*DEGRAD
   COSZ=SIN(LATITUDE*DEGRAD)*SIN(DECLIN)+COS(LATITUDE*DEGRAD)*COS(DECLIN)*COS(HRANG)

   IF (PRESENT(HRANG_OUT ))  HRANG_OUT = HRANG
   IF (PRESENT(DECLIN_OUT)) DECLIN_OUT = DECLIN
   IF (PRESENT(GMT_OUT   ))    GMT_OUT = GMT
   IF (PRESENT(JULDAY_OUT)) JULDAY_OUT = JULDAY

 END SUBROUTINE CALC_DECLIN

!---------------------------------------------------------------------

 SUBROUTINE LOCAL_TIME(NOWDATE, LONGITUDE, TLOCTIM)

   IMPLICIT NONE
!---------------------------------------------------------------------
   CHARACTER(LEN=19),      INTENT(IN)  :: NOWDATE    ! YYYY-MM-DD_HH:MM:SS
   REAL(kind=kind_noahmp), INTENT(IN)  :: LONGITUDE
   REAL(kind=kind_noahmp), INTENT(OUT) :: TLOCTIM
   INTEGER                             :: IHOUR
   INTEGER                             :: IMINUTE
   INTEGER                             :: ISECOND

   READ(NOWDATE(12:13), *) IHOUR
   READ(NOWDATE(15:16), *) IMINUTE
   READ(NOWDATE(18:19), *) ISECOND

   TLOCTIM = REAL(IHOUR) + REAL(IMINUTE)/60.0 + REAL(ISECOND)/3600.0 + LONGITUDE/15.0 ! LOCAL TIME IN HOURS
   TLOCTIM = AMOD(TLOCTIM+24.0, 24.0)

 END SUBROUTINE LOCAL_TIME


end module module_NoahMP_hrldas_driver

