! ==============================================================================

!                               Model Overview
! Model Name: TECO-CNP
! Purpose: Develop a coupled Carbon-Nitrogen-Phosphorus (C-N-P) biogeochemical model.
! Programmer: Fangxiu Wan
! Submitted Version: TECO-CNP V1.0
! Version Date: February 2025
! Reference: 
!   TECO-CNP v1.0 - A coupled carbon-nitrogen-phosphorus model with data assimilation 
!   capabilities for subtropical forests, Wan et al., 2025
! ==============================================================================

!                              Important Statements
! ------------------------------------------------------------------------------
!
! This code is constructed based on the work of Weng & Luo (2008).
! 
! The core processes of the carbon cycle, particularly the photosynthesis module,
! primarily inherit from the 2008 version. The nitrogen and phosphorus cycles,
! along with their interactions with carbon, were developed by WFX.
!
! The model currently operates in single precision and may be updated
! according to future research requirements.
!
! For detailed interpretation of the code, please refer to Wan et al. (2025) 
!
! Copyright notice is required for any research utilizing this code.

! Module statements:
! The main module - TECO-CNP_main.f90 - Implements the coupled C-N-P cycle for terrestrial ecosystems.
! Included Main Modules:
! - FileSize.f90        : Defines array sizes and file sizes
! - IntersVariables.f90 : Defines global variables
! - outputs_mod.f90     : Manages output file generation
! - FOR_MCMC.f90        : Contains the data assimilation module
! - vars_site.f90       : Holds site-specific parameters

! ==============================================================================


PROGRAM MAIN

! Main program for TECO-CNP
! Overview:
! This program executes various simulations based on the command-line arguments provided.
! The configuration options include:
! 1. Normal Simulation [CYCLE_CNP]:
!    - Options: (a) C-only, (b) CN, (c) CNP
!    - Related subroutine: \CNP_simu\
!    - Note: A spin-up procedure is included within \CNP_simu\
!
! 2. Sensitivity Analysis [SensTest]:
!    - Customize variables of interest according to specific research objectives.
!    - Related subroutine: \SensitivityTest\
!
! 3. Data Assimilation [MCMC]:
!    - Related subroutine: \MCMC_simu\

! Programmer: Fangxiu Wan
! Version Date: February 12 2025

    USE FileSize
    USE IntersVariables
    USE outputs_mod
    USE FOR_MCMC
    USE vars_site
    IMPLICIT NONE

    INTEGER istat1,m,n,nn,i,j,upgraded
    INTEGER :: nseed
    INTEGER, ALLOCATABLE :: seed(:)
    CHARACTER(len=99) commts,line
    CHARACTER(len=50) arg
    NPaddition = 0   ! 默认：不进入施肥实验情景（非常关键，避免随机值）
    IF (command_argument_count() == 7 .OR. command_argument_count() == 8 .OR. command_argument_count() == 9) THEN

    CALL get_command_argument(1, arg)
    READ(arg,*) start_year

    CALL get_command_argument(2, arg)
    READ(arg,*) end_year

    CALL get_command_argument(3, arg)
    READ(arg,*) CYCLE_CNP

    CALL get_command_argument(4, arg)
    READ(arg,*) MCMC

    CALL get_command_argument(5, arg)
    READ(arg,*) NDSPINUP

    CALL get_command_argument(6, arg)
    READ(arg,*) nspinup

    CALL get_command_argument(7, arg)
    READ(arg,*) SensTest

    IF (command_argument_count() >= 8) THEN
        CALL get_command_argument(8, arg)
        READ(arg,*) NPaddition
    ELSE
        NPaddition = 0
    ENDIF

    IF (command_argument_count() >= 9) THEN
        CALL get_command_argument(9, arg)
        READ(arg,*) P_add_rate_kg_ha_yr
    ELSE
        P_add_rate_kg_ha_yr = 0.0
    ENDIF

ELSE
    WRITE(*,*) "Error: Expected 7, 8 or 9 arguments."
    WRITE(*,*) "Usage: ./teco_cnp.exe start_year end_year CYCLE_CNP MCMC NDSPINUP nspinup SensTest"
    WRITE(*,*) "Optional: [NPaddition] [P_add_rate_kg_ha_yr]"
    STOP
ENDIF

WRITE(*,*) "NPaddition = ", NPaddition
WRITE(*,*) "P_add_rate_kg_ha_yr = ", P_add_rate_kg_ha_yr

    CALL Filespath() ! Define the path for output based on simulation type
    CALL ReadInitialFile() ! Read the initial data from input files

    !START SIMULATION

    IF (MCMC.eq.0) THEN ! START BIOGEOCHEMICAL CYCLE SIMULATION
        WRITE(*,*) 'starting simulation,choosed coupled scheme is',CYCLE_CNP
        ispinup = ispinup + 1
        !1- REGULAR SIMULATION
        CALL CNP_simu(upgraded)          
        !2- SENSITIVITY ANALYSIS - alternative action
        IF(SensTest .eq. 1) THEN
            CALL SensitivityTest()
            PRINT*,'NOW RUN SENSITIVITY TEST'
        ENDIF       
        WRITE(*,*)'Here, finish the simulation, iloops, CYCLE_CNP:',ispinup, CYCLE_CNP
        ! OUTPUTS
        CALL FilesForOutput ()   !Configuring the file path for model output
        CALL WriteFiles_noMCMC() !Customize the outputs based on research requirements
    ELSE!START MCMC
        WRITE(*,*) 'STARTING MCMC'
        ! SET SEED
        CALL RANDOM_SEED(size = nseed)
        ALLOCATE(seed(nseed))
        seed = 123
        CALL RANDOM_SEED(put = seed)
        ! OPEN FILE
        OPEN(63,file=TRIM(ADJUSTL(fileplace))//'Paraest_1000'//filesignal)
        OPEN(72,file=TRIM(ADJUSTL(fileplace))//'covvariance_temp'//filesignal)
        ! RUN MCMC - GENERATE ONE CHAIN        
        CALL MCMC_simu(upgraded)
        DEALLOCATE(seed)
    ENDIF!END of the simulation or MCMC 
    stop
END PROGRAM MAIN

SUBROUTINE CNP_simu(upgraded)

! Biogeochemical cycle of TECO-CNP
!
! Overview:
! This subroutine simulates the biogeochemical cycles for various coupling schemes.
! The carbon cycle is based on the work of Weng & Luo (2008) and Wan et al. (2025).
! The nitrogen and phosphorus cycles are modeled according to Wan et al. (2025).
!
! Programmer: Fangxiu Wan
! Version Date: February 12 2025

    USE FileSize
    USE IntersVariables
    USE DaysHours
    USE LIMITATION
    USE SPINUP_mod
    USE vars_site
    USE vars_consts
    USE outputs_mod
    USE update_traits
    IMPLICIT NONE 

    INTEGER iiterms
    REAL co2ca,co2
    REAL,DIMENSION(10):: wupl,evapl
    REAL runoff,rain
    ! variables for canopy model
    REAL ET,G,wind,eairp,esat,rnet,Esoil
    REAL,DIMENSION(3):: reffbm,reffdf,extkbm,extkdm
    REAL,DIMENSION(2):: radabv
    REAL Qcan(3,2)
    ! parameters for photosynthesis model
    ! additional arrays to allow output of info for each layer
    REAL,DIMENSION(5):: RnStL,QcanL,RcanL,AcanL,EcanL,HcanL
    REAL,DIMENSION(5):: GbwcL,GswcL,hG,hIL
    REAL,DIMENSION(5):: Gaussx,Gaussw,Gaussw_cum 
    ! for phenology
    REAL totlivbiom
    REAL L_fall,L_add,litter,seeds
    REAL stor_use
    ! respiratiom
    REAL RmLeaf,RmStem,RmRoot,RmRe    ! maintanence respiration
    REAL RgLeaf,RgStem,RgRoot,RgRe    ! growth respiration
    REAL RaLeaf,RaStem,RaRoot,RaRe   
    REAL Ressoil,Rtotal !Ressoil - hetero respiration
    ! climate variables for every day
    REAL Ta,Tair,Ts,Tsoil,TaDaily
    REAL doy,hour,Dair,Rh,radsol,PAR
    ! For loops
    INTEGER year,i,j,m,n,upgraded   
    INTEGER idoy,ileaf
    !for nitrogen sub-model
    REAL CNmin,CNmax,NSNmax,NSNmin
    REAL QNplant
    REAL SNvcmax,SNgrowth,SNRauto,SNrs
    REAL SNvcmax_use,SNgrowth_use,SPvcmax_use,SPgrowth_use
    !for phosphorus sub-model
    REAL P_leach  !p_LEACH = P_LOSS
    REAL QPplant
    REAL ksye,Vcmax0_update,SLAx_update
    !INTEGER Days,Hours
    INTEGER simu_year,simu_days,EndDays
    INTEGER simu_hours,EndHours,start_day,end_day,count
    REAL (SGL) :: b1,b2,b3,b4
    REAL (SGL),DIMENSION(10,1) :: Xct,dXdt,BaseTau,EcoTau,Xpt_yr,delta_QC_yr
    REAL (SGL) k1_lastloop,ksi_lastloop
    TYPE(ResultType) :: result
    !!!

    INTEGER :: ihour_in_day
    REAL    :: P_annual, P_pulse

    CALL consts()
    CALL site_value()
    count = 0
    ! DEBUG
    Flag_Growth = 0
    Flag_C = 0
    Flag_N = 0 
    Flag_P = 0
    Flag_GPP = 0
    Flag_QP = 0
    Flag_NSC = 0
    rejectnow = 0
    !TEST FOR NUTRIENT ADDITION EXP.
    SLAx_update = 0.0
    Vcmax0_update= 0.0
    ! TEST FOR LEAF PHENOLOGY
    Tairmax = maxval(forcing_data(4,:))
    Tairmax_Loc = maxloc(forcing_data(4,:))
    Tairmax_Loc_value = Tairmax_Loc(1)

DO  ! CYCLE FOR SPIN-UP
    CALL loop_initial_value(iiterms,simu_year,EndDays,EndHours,stor_use)
    WRITE(*,*) 'P parameters inside model:', &
        s_sorb_to_ssorb, s_ssorb_to_occ
    DO itime=1,iiterms
!       Reset the variables each year 每年重置变量
        IF (EndHours .EQ. 0)THEN
            storage=accumulation
            stor_use=storage/times_storage_use
            accumulation=0.0
            GDD5 = 0.
            onset = 0
            iyear = iyear+1
            simu_year = start_year+iyear
            result= Days_cal(simu_year)
            EndDays=result%Days
            EndHours=result%Hours
            WRITE(*,*)'EndDays,EndHours,simu_year,idays',EndDays,EndHours,simu_year,idays
        ENDIF  ! End reset  

        N_fert = 0.0
        P_fert = 0.0

        ! ----------------------------------------------------------
        ! P addition experiment
        ! P_annual unit: kg P ha-1 -> g P m-2
        ! ----------------------------------------------------------
        IF (NPaddition .EQ. 1) THEN

            P_annual = P_add_rate_kg_ha_yr * 0.1
            P_pulse  = P_annual / 4.0
            P_fert   = 0.0

            ! Directly obtain date from forcing data, so the code does
            ! not depend on whether the simulation begins in 2021 or 2022.
            year = INT(forcing_data(1, itime))
            idoy = INT(forcing_data(2, itime))
            ihour_in_day = INT(forcing_data(3, itime))

            ! Apply P only at the first hour of the selected day.
            IF (ihour_in_day .EQ. 1) THEN

                ! 2022: one full annual dose on 1 July.
                IF (year .EQ. 2022 .AND. idoy .EQ. 182) THEN
                    P_fert = P_annual
                ENDIF

                ! 2023: four equal additions.
                IF (year .EQ. 2023) THEN
                    IF (idoy .EQ. 1   .OR. idoy .EQ. 91  .OR. &
                        idoy .EQ. 182 .OR. idoy .EQ. 274) THEN
                        P_fert = P_pulse
                    ENDIF
                ENDIF

                ! 2024 is a leap year.
                IF (year .EQ. 2024) THEN
                    IF (idoy .EQ. 1   .OR. idoy .EQ. 92  .OR. &
                        idoy .EQ. 183 .OR. idoy .EQ. 275) THEN
                        P_fert = P_pulse
                    ENDIF
                ENDIF

                IF (P_fert .GT. 0.0) THEN
                    WRITE(*,*) 'P addition: year=', year,       &
                               ' doy=', idoy,                   &
                               ' rate=', P_add_rate_kg_ha_yr,   &
                               ' P_fert=', P_fert
                ENDIF

            ENDIF
        ENDIF

        
        StemSap=AMIN1(Stemmax,SapS*bmStem)
        RootSap=AMIN1(Rootmax,SapR*bmRoot)
        NSCmax=0.01*(StemSap+RootSap+QC(1))

        IF(NSC.LE.NSCmin)fnsc=0.0
        IF(NSC.ge.NSCmax)fnsc=1.0
        IF((NSC.LT.NSCmax).and.(NSC.GT.NSCmin))THEN 
            fnsc=(NSC-NSCmin)/(NSCmax-NSCmin)
        ENDIF     
!------------------------------------------------
        ! Nutrient limitation factors
        CALL xnp(alphaN,alphaP)
        x_leaf_limit = min(xnp_leaf_limit,1.)
        x_stem_limit = min(xnp_stem_limit,1.)
        x_root_limit = min(xnp_root_limit,1.)
        xuptake      = min(xNPuptake,1.)
        !forcing data
        CALL step_forcing (doy,hour,Tair,Tsoil,RH,Dair,Rain,wind,&
        & PAR,radsol,co2,co2ca)
        ! Ajust some unreasonable values
        RH=AMAX1(0.01,AMIN1(99.99,RH))                      
        eairP = esat(Tair)*RH/100.              ! Added for SPRUCE, due to lack of VPD data
        !Dair=esat(Tair)-eairP                  ! wan- commented out, since we have VPD data in TianTong
        radsol=AMAX1(radsol,0.01)

        IF(radsol.GT.10.0) THEN
                G=-25.0  !wan-see Garratt, pp116, G - the heat flux into soil            
        ELSE
                G=20.5
        ENDIF
        Esoil=0.05*radsol
        IF(radsol.LE.10.0) Esoil=0.5*G		    ! If solar radiation less than 10...
  
        ! PHOTOSYNTHESIS FROM LEAF TO CANOPY
        CALL canopy(doy,hour,radsol,Tair,Dair,eairP,wind,rain,Rnet,G,Esoil,&
            &              Tsoil,co2ca,pi,tauL,rhoL,rhoS,emleaf,emsoil, &
            &              Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,&
            &              H2OMw,chi,Dheat,wleaf,gsw0,theta,&
            &              conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,&
            &              Edjm,Entrpy,gam0,gam1,gam2,AcanL)
        rain = rain
        ! SOIL WATER
        CALL soilwater(Tair,RH,rain,runoff,upgraded)                
            ET=evap+transp
        ! PLANT RESPIRATION - auto-respiration
        CALL respiration(Tair,Tsoil,SNRauto,RmLeaf,RmStem,RmRoot,RmRe)
        ! PLANT GROWTH - allocation, phenology, growth
        CALL plantgrowth(Tair,stor_use,SNgrowth,L_fall,& 
            &                RgLeaf,RgStem,RgRoot,RgRe,count,upgraded)
        ! CARBON TRANSFER MODULE
        CALL TCS_CNP(Tair,Tsoil,runoff,L_fall,&
        &               CNmin,CNmax,NSNmax,NSNmin,&       ! nitrogen
        &               SNvcmax,SNgrowth,SNRauto,SNrs,&
        &               P_leach)
    !This loop for: if a new set of paras get wrong simulation,
    !CNP_simu will stop run and goto MCMC module to regenerate new paras 

        IF(CYCLE_CNP .eq. 1)THEN
            Rnitrogen = 0.0 
            Rphosphorus=0.0
        ENDIF
        IF(CYCLE_CNP .eq. 2) Rphosphorus = 0.0 

        EC_out = EC_store * exitEC
        EC_store = EC_store + ExcessC - EC_out
        Rauto = Rmain+Rgrowth!+Rnitrogen+Rphosphorus+ExcessC
        Rauto_new = Rauto + Ec_out

        IF (ISNAN(NSC))THEN
            write(*,*)'NSC is NA'
            IF(MCMC .EQ. 0)stop
            Flag_NSC = 1
            !stop
        ENDIF

        NSC    =NSC+GPP-Rauto-(NPP-add)-store-Rnitrogen-Rphosphorus-ExcessC 

        IF(NSC<0)THEN            !wan- NSC can be smaller then 0 20230516
                bmstem=bmstem+NSC/0.48
                NPP=NPP+NSC
                NSN=NSN-NSC/CN(2)
                NSP=NSP-NSC/CP(2) 
                NSC=0.
        ENDIF
        !update
        RaLeaf = RgLeaf+ RmLeaf  
        RaStem = RgStem + RmStem
        RaRoot = RgRoot + RmRoot !+ Rnitrogen
        RaRe = RgRe+ RmRe

        Rhetero= Rh_pools(1)+Rh_pools(2)+Rh_pools(3)&
        &         +Rh_pools(4)+Rh_pools(5)
        Ressoil  =Rhetero+RmRoot+RgRoot+Rnitrogen
        
        NEE= Rauto+Rhetero-GPP
        NEE_new= Rauto+Rhetero+Rnitrogen+Rphosphorus+store+ExcessC-GPP-add
        NEE_new= Rauto+Rhetero+Rnitrogen+Rphosphorus+ExcessC-GPP

        !NEE = NEE_new

        Reco=Rauto+Rhetero
        bmleaf=QC(1)/fbmC
        bmstem=QC(2)/fbmC
        bmroot=QC(3)/fbmC
        bmRe = QC(4)/fbmC
        bmplant=bmleaf+bmroot+bmstem+bmRe
        LAI=bmleaf*SLA     

        IF(isnan(gpp))THEN
            WRITE(*,*)'gpp is nan,itime',itime
            Flag_GPP = 1
        ENDIF

!-----------Reset after 1 year
!           Update the phenological variables after each year simulation
            !update the GDD5 every day (Calculate the cumulative temperature)
            simu_hours = EndHours -1
            EndHours = simu_hours
            IF(mod(itime,24).eq. 0)THEN
                Ta = sum(forcing_data(4,(itime-23):itime))/24
                IF(Ta.GT.5.0)GDD5 = GDD5+Ta-5.0   !wan 20230707
                simu_days = EndDays -1
                EndDays = simu_days
            ENDIF  
            IF (MCMC .eq. 1)THEN
                IF (Flag_Growth .eq. 1) rejectnow = 1
                IF (Flag_C .eq. 1) rejectnow = 2  !QC
                IF (Flag_N .eq. 1) rejectnow = 3  !QNMINER
                IF (Flag_P .eq. 1) rejectnow = 4  !QPLAB
                IF (Flag_GPP .eq. 1) rejectnow = 5
                IF (Flag_QP .eq. 1) rejectnow = 6
                IF (Flag_NSC .eq. 1) rejectnow = 7
            ENDIF      

            IF (rejectnow .ne. 0) goto 111 

    ! The results were recorded at hourly and daily time scalar
    ! after one steps (/hour) is finishe
    CALL ForRecords(AcanL)
        IF (EndDays .EQ. 0) THEN
            write(*,*)'simu_year,nyear',simu_year,nyear
            if (check_leap(simu_year) .EQ. 1) then
                ndays = 366
            else
                ndays = 365
            end if
            start_day = idays-ndays+1
            end_day   = idays
            write(*,*)'idays',idays,idays-ndays
            !CALL ForYearRecords(start_day,end_day,simu_year)
        ENDIF
    ENDDO !END one [CNP_simu], ndays = 365 or 365*nspinup (if there is no leap year)

    
    !---------------------------------------------------------------------------
    ! After one loop, i.e. the simulation diven by one period of climate forcing
    ! The nspinup is not used, the whole simulation stop when meet the conditions

    IF (NDSPINUP .EQ. 1) THEN
        CALL Spinup_output(ispinup)  !Sunbroutine for recording the output of NDspinup
        IF (ispinup .EQ. nspinup) THEN
            WRITE (*,*)'QC',QC
            exit
        ENDIF
        ispinup =ispinup+1
    ENDIF

    IF (NDSPINUP .EQ. 0) exit

ENDDO  !END spinup

WRITE(*,*)'END SIMU'

111 continue ! END simulation for MCMC
! Flag for MCMC
IF(rejectnow .ne. 0)THEN
    write(*,*)'QC',QC
    write(*,*)'QPlab,QNminer,GrowthP,GPP',&
    & QPlab,QNminer,GrowthP,GPP,QC(1),QN(1),QP(1),LfactorP,x_leaf_limit,Vcmax0,SLA
    write(*,*)'QP',QP
    write(*,*)'itime,rejectnow',itime,rejectnow
    !stop
ENDIF
RETURN
END SUBROUTINE CNP_simu


SUBROUTINE canopy(doy,hour,radsol,tair,Dair,eairP,&
    &               wind,rain,&
    &               Rnet,G,Esoil,&
    &               Tsoil,& !constants specific to soil and plant
    &               co2ca,&
    &               pi,tauL,rhoL,rhoS,emleaf,emsoil,&
    &               Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,&
    &               H2OMw,chi,Dheat,wleaf,gsw0,theta,&
    &               conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,&
    &               Edjm,Entrpy,gam0,gam1,gam2,AcanL)

    ! Overview: simulate canopy photosynthesis based on two-leaf model

     USE IntersVariables
     USE vars_site

      REAL doy
      REAL tauL(3),rhoL(3),rhoS(3),reffbm(3),reffdf(3)
      REAL extkbm(3),extkdm(3)
      REAL radabv(2),Qcan(3,2)
!     extra variables used to run the model for the wagga data
      INTEGER idoy,ileaf
      INTEGER i,j,k
!     additional arrays to allow output of info for each layer
      REAL RnStL(5),QcanL(5),RcanL(5),AcanL(5),EcanL(5),HcanL(5)
      REAL GbwcL(5),GswcL(5),hG(5),hIL(5)
      REAL Gaussx(5),Gaussw(5),Gaussw_cum(5)
      CHARACTER*80 commts
!     Normalised Gaussian points and weights (Goudriaan & van Laar, 1993, P98)
!     5-point
      data Gaussx/0.0469101,0.2307534,0.5,0.7692465,0.9530899/
      data Gaussw/0.1184635,0.2393144,0.2844444,0.2393144,0.1184635/
      data Gaussw_cum/0.11846,0.35777,0.64222,0.88153,1.0/

!     calculate beam fraction in incoming solar radiation
      CALL  yrday(doy,hour,lat,radsol,fbeam)
      idoy=int(doy)
      hours=idoy*1.0+hour/24.0
      coszen=sinbet(doy,lat,pi,hour)             !cos zenith angle of sun

!     set windspeed to the minimum speed to avoid zero Gb
      IF(wind.LT.0.01) wind=0.01
!     calculate soil albedo for NIR as a function of soil water (Garratt pp292)
      IF(topfws.GT.0.5) THEN
            rhoS(2)=0.18            !wan- original, 0.18
      ELSE
            rhoS(2)=0.52-0.68*topfws  !wan- original, 0.52-0.68*topfws
      ENDIF
!        assign plant biomass and leaf area index at time t
!        assume leaf biomass = root biomass
      FLAIT =LAI 
      !wan - eairP=esat(Tair)-Dair                !air water vapour pressure
      !wan - eairP=esat(Tair)*(RH/100)
      radabv(1)=0.5*radsol                 !(1) - solar radn
      radabv(2)=0.5*radsol                 !(2) - NIR
!     CALL multilayer model of Leuning - uses Gaussian integration but radiation scheme
!     is that of Goudriaan
      CALL xlayers(Tair,Dair,radabv,G,Esoil,fbeam,eairP,          &
     &           wind,co2ca,fwsoil,wcl,LAI,coszen,idoy,hours,     &
     &           tauL,rhoL,rhoS,xfang,extkd,extkU,wleaf,          &
     &           Rconst,sigma,emleaf,emsoil,theta,a1,Ds0,         &
     &           cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,       &
     &           gsw0,alpha,stom_n,wsmax,wsmin,                   &
     &           Vcmx0,eJmx0,conKc0,conKo0,Ekc,Eko,o2ci,          &
     &           Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,       &
     &           extKb,Rsoilabs,Hsoil,Acan1,Acan2,Ecan1,Ecan2,    &
     &           RnStL,QcanL,RcanL,AcanL,EcanL,HcanL,GbwcL,GswcL,gddonset,itime,&
     &           Aleaf,Tlk1,Tlk2,fslt,fshd,                       &
     &           XX_rec,XX_med_rec,gscx_rec,gscx_med_rec,EqStomata)
         Acanop=Acan1+Acan2
         Ecanop=Ecan1+Ecan2
         gpp=Acanop*3600.0*12.0                           ! every hour, conv mol m-2 s-1 to g C m-2 h-1
	     transp=AMAX1(Ecanop*3600.0/(1.0e6*(2.501-0.00236*Tair)),0.) ! mm H2O /hour
         evap=AMAX1(Esoil*3600.0/(1.0e6*(2.501-0.00236*Tair)),0.)    ! wan- 1.0e6*(2.501-0.00236*Tair ) lamda
!        H2OLv0=2.501e6               !latent heat H2O (J/kg)

      RETURN
 END
!****************************************************************************
!   autotrophic respiration
SUBROUTINE respiration(Tair,Tsoil,SNRauto,RmLeaf,RmStem,RmRoot,RmRe)

    ! Overview: calculate plant and soil respiration by the following equation:
    !           RD=BM*Rd*Q10**((T-25)/10) (Sun et al. 2005. Acta Ecologica Sinica)
!     
    USE vars_site
    USE IntersVariables
    IMPLICIT NONE
    REAL Tair,Tsoil
    REAL RmLeaf,RmStem,RmRoot,RmRe
    REAL SNRauto,SNRauto_use
    REAL SPRauto,SPRauto_use
    REAL conv                  ! converter from "umol C /m2/s" to "gC/m2/hour"

    conv=3600.*12./1000000.    ! umol C /m2/s--> gC/m2/hour

!-----------------------------
!    IF(LAI.GT.LAIMIN) THEN
    IF(GPP .GT. 0) THEN !wan 2024 0306
        !RmLeaf=Rl0*bmleaf*0.48*SLA*0.1     &               
        !   &   *Q10**((Tair-10.)/10.)*0.02*conv
        !RmStem=Rs0*StemSap*0.001*Q10**((Tair-10.)/10.)*0.02*conv   !wan- change (Tair-25) to (Tair -10)
        !RmRoot=Rr0*RootSap*0.001*Q10**((Tair-10.)/10.)*0.02*conv   !wan- change (Tair-25) to (Tair -10)
        !RmRe=  Rre0*bmRe*0.8*0.001*Q10**((Tair-10.)/10.)*0.02*conv
        !wan- replace above equations, from TECO.for. similarity results between two methods. 2023.
        Rl0=bmleaf*SLA*0.1            !umol/m2/s
	    Rs0=bmstem*0.006*0.1          !umol/m2/s
	    Rr0=bmroot*0.006*0.1          !umol/m2/s
	    RmLeaf=Rl0*exp(0.084*Tair)*fnsc*conv 
	    RmStem=Rs0*exp(0.084*Tair)*fnsc*conv 
	    RmRoot=Rr0*exp(0.084*Tair)*fnsc*conv 
    ELSEIF (GPP .EQ. 0.) THEN !wan, 2024 0306
        ! Plant will use the carbon in NSC pool to maintain the life activity
        RmLeaf=0.3*0.001*NSC    
        RmStem=0.3*0.001*NSC
        RmRoot=0.3*0.001*NSC
        RmRe  =0.1*0.001*NSC
    ENDIF
    Rmain=Rmleaf+Rmstem+RmRoot+RmRe
    !rescale
    IF(Rmain > 0.0015*NSC)THEN            !If Total autotropic respiration greater than 0.15% of Nonstructure Carbon
        Rmleaf=Rmleaf/Rmain*0.0015*NSC    !in TECO.for RaL=RaL*NSC*0.6/GPP_Ra
        Rmstem=Rmstem/Rmain*0.0015*NSC
        RmRoot=RmRoot/Rmain*0.0015*NSC
        RmRe=RmRoot/Rmain*0.0015*NSC
        Rmain=Rmleaf+Rmstem+RmRoot+RmRe
    ENDIF
    !For debug
    IF(Rmain .LT. 0.)THEN
        WRITE(*,*)'Q10**((Tair-15.)/10.)',Q10**((Tair-15.)/10.)
        WRITE(*,*)'NSC',NSC
    ENDIF

    RETURN
END ! sub respiration

SUBROUTINE soilwater(Tair,RH,rain,runoff,upgraded)  !outputs

    ! Note: the unit of water is 'mm', soil moisture or soil water content is a ratio

    USE vars_site
    USE IntersVariables
    IMPLICIT NONE
!   soil traits
    REAL wsmaxL(10),wsminL(10) !from input percent x%
!   plant traits
    INTEGER nfr,upgraded
!   climate conditions
    REAL rain ! mm/hour
!   output from canopy model
    
!   output variables
    REAL fw(10),ome(10)
    REAL depth(10),WUPL(10),EVAPL(10),SRDT(10)
    REAL plantup(10)
    REAL Tsrdt
    
!  REAL fwcln(10) !  fraction of water in layers, like field capacity
    REAL DWCL(10),Tr_ratio(10)
    REAL wtadd,twtadd,runoff,tr_allo

    REAL exchangeL,supply,demand,omegaL(10)
    INTEGER i,j,k
    REAL infilt_max
!   CHARACTERs annotation for water table module  -MS
!   wan- water table module is neglectable for TECO_CNP
    REAL vtot,zmax,thetasmin,zthetasmin,az
    REAL zwt1,zwt2,zwt3
!   water table CHARACTERs annotation END here  -MS
!   wan- calculate soil evaporation, [Weng & Luo,2008]    
    REAL wtdeficit(10),roff_layer
    REAL Rsoil,Rd,P,density,la,sp_heat,psychro,Tair,RH,esat,evap_soil

!   infilt_max=15.           !the value should be larger in peatland with higher porosity
!   infilt_max=50.           !changed by Mary on 11/14/2016   optimization to the warming scenarios
    infilt_max=30.           !wan- for Tiantong , ! no change in results
    WILTPT =wsmin/100.0      !ratio
    FILDCP =wsmax/100.0      !ratio
    
    DO i=1,10
        wtdeficit(i)=0.0
        dwcl(i)=0.0   
        evapl(i)=0.0
        WUPL(i)=0.0
        SRDT(i)=0.0
        DEPTH(i)=0.0
    ENDDO

!   Determine which layers are reached by the root system. 
!   Layer volume (cm3)
    DEPTH(1)=10.0
    DO i=2,10
        DEPTH(i)=DEPTH(i-1)+thksl(i)
    ENDDO
    DO i=1,10
        IF(rdepth.GT.DEPTH(i)) nfr=i+1
    ENDDO
    IF (nfr.GT.10) nfr=10
!   wan- water deficite (in layers)
    do i=1,10
	    if(FILDCP.GT.wcl(i)) wtdeficit(i)=FILDCP-wcl(i)
	enddo

! water infiltration through layers
!   infilt=infilt+rain  !mm/hour   !wan- changed equation, infilt = precp, 2023.
    infilt = rain
    TWTADD=0
    roff_layer=0.0
!   Loop over all soil layers.
DO i =1,10   
    IF(infilt.GT.0.0)THEN           ! corrected by Mary on 11/14/2016 to fix the bug in 'WTADD 'NAN''
!       Add water to this layer, pass extra water to the next.
        WTADD=AMIN1(INFILT,wtdeficit(i)*thksl(i)*10.0)
        !wan- commented out,  WTADD=AMIN1(INFILT,infilt_max,AMAX1((FILDCP-wcl(1))*thksl(1)*10.0,0.0)) !wan- thksl*10?
        !changed by Mary ON 11/14/2016 AMAX doesn't seem make sense and the result does not change with this correction
        !wan- commented out,  WTADD=AMIN1(INFILT,infilt_max,(FILDCP-wcl(1))*thksl(1)*10.0)   
!       change water content of this layer
        wcl(i)=(wcl(i)*(thksl(i)*10.0)+WTADD)/(thksl(i)*10.0)    !wan- wcl unit % relative soil content
!       FWCLN(I)=wcl(I)           !  /VOLUM(I)! update fwcln of this layer
        TWTADD=TWTADD+WTADD       !calculating total added water to soil layers (mm)
        INFILT=INFILT-WTADD       !update infilt
    ENDIF

    if(infilt.GT.0.0)THEN  !runoff [TECO.for]
		roff_layer=roff_layer + INFILT*0.05*(i-1) 
        INFILT=INFILT -  INFILT*0.05*(i-1)
	endif
ENDDO
    runoff=INFILT + roff_layer
!wan- commented out
!   runoff
!   runoff=INFILT*0.0019
!   runoff=0.0
    !runoff=INFILT*0.0019    ! change by zhoujian
    !infilt = infilt-runoff
!wan- end 

!wan- unresonable mechanism, add by wan, from TECO.for
    !supply=0
    !if(rain.GT.0.0.and.wcl(1).GT.wcl(2))then
	!	supply=(wcl(1)-wcl(2))/3.0
	!	wcl(1)=wcl(1)-2.0*supply
	!	wcl(2)=wcl(2)+supply
	!endif
!********************************************************
!wan- commented out, intricate equations, can't find the source
    !IF (transp .GT. 0.2 .and. transp .LE. 0.22) THEN
    !    infilt = infilt+transp*0.4
    !ELSE IF (transp .GT. 0.22) THEN
!        infilt = infilt+infilt*0.0165
   !     infilt = infilt+transp*0.9
!        infilt = infilt+0.22*0.4+(transp-0.22)*0.9
    !ELSE
    !    infilt = infilt+transp*0.001
    !ENDIF
!    
    !IF (evap .ge. 0.1 .and. evap .LE. 0.15) THEN
    !    infilt = infilt+evap*0.4
    !ELSE IF (evap .GT. 0.15) THEN
    !    infilt = infilt+evap*0.9
    !ELSE
    !    infilt = infilt+evap*0.001
    !ENDIF  
!wan- end

!   water redistribution among soil layers
    DO i=1,10
        wsc(i)=Amax1(0.00,(wcl(i)-WILTPT)*thksl(i)*10.0)      !wan, 2023/3/21 *10, convert thksl [cm] to mm
        omegaL(i)=Amax1(0.001,(wcl(i)-WILTPT)/(FILDCP-WILTPT))
    ENDDO
    
    supply=0.0
    demand=0.0
    DO i=1,9
        IF(omegaL(i).GT.0.3)THEN
            supply=wsc(i)*omegaL(i)
            demand=(FILDCP-wcl(i+1))*thksl(i+1)*10.0      &
                &               *(1.0-omegaL(i+1))
            exchangeL=AMIN1(supply,demand)
            wsc(i)=wsc(i)- exchangeL
            wsc(i+1)=wsc(i+1)+ exchangeL
            wcl(i)=wsc(i)/(thksl(i)*10.0)+WILTPT
            wcl(i+1)=wsc(i+1)/(thksl(i+1)*10.0)+WILTPT
        ENDIF
    ENDDO
!wan- commented out, no need for another runoff
    wsc(10)=wsc(10)-wsc(10)*0.00005
    runoff = runoff+wsc(10)*0.00005
    wcl(10)=wsc(10)/(thksl(10)*10.0)+WILTPT
!    END of water redistribution among soil layers
!wan- add from TECO.for [Weng & Luo, 2008], eq.A9
!calculate evap demand by eq.24 of Seller et al. 1996 (demand)
	if(wcl(1).LT.WILTPT)then
	    evap_soil=0.0
	ELSE
		Rsoil=10.1*exp(1.0/wcl(1))
		Rd=   20.5 !*exp(LAI/1.5)!LAI is added by Weng
		P=101325.0  !Pa, atmospheric pressure
		density=1.204 !kg/m3
		la=(2.501-0.00236*Tair)*1000000.0 !J/kg  !wan-lamuda
		sp_heat=1012.0  !J/kg/K
		psychro=1628.6*P/la                      

		evap_soil=1.0*esat(Tair)*(1.0-RH/100.0)/ &
     &         (Rsoil+Rd)*density*sp_heat/psychro/la*3600.0
	endif
!wan- end add

!   Redistribute evaporation among soil layers
    Tsrdt=0.0
    DO i=1,10
!   Fraction of SEVAP supplied by each soil layer
    SRDT(I)=EXP(-6.73*(DEPTH(I)-thksl(I)/2.0)/100.0) !/1.987       
    Tsrdt=Tsrdt+SRDT(i)  ! to normalize SRDT(i)               
    ENDDO

!wan- merged in the eq.EVAPL
    !do i=1,10
	!	SRDT(i)=SRDT(i)/Tsrdt                             
	!enddo
    
    DO i=1,10
        EVAPL(i)=Amax1(AMIN1(evap*SRDT(i)/Tsrdt,wsc(i)),0.0)    !mm  wan- change evap_soil to evap
        DWCL(i)=EVAPL(i)/(thksl(i)*10.0) !ratio
        wcl(i)=wcl(i)-DWCL(i)
    ENDDO

    evap=0.0       
    DO i=1,10
        evap=evap+EVAPL(I)
    ENDDO

!   Redistribute transpiration according to root biomass
!   and available water in each layer
    tr_allo=0.0                     !transpiration total
    DO i=1,nfr
        tr_ratio(i)=FRLEN(i)*wsc(i) !no big change when commented*(wcl(i)-WILTPT)) !*thksl(I))
        tr_allo=tr_allo+tr_ratio(i)
    ENDDO

    DO i=1,nfr
        plantup(i)=AMIN1(transp*tr_ratio(i)/tr_allo, wsc(i)) !mm              
        wupl(i)=plantup(i)/(thksl(i)*10.0)
        wcl(i)=wcl(i)-wupl(i)
    ENDDO

    transp=0.0
    DO i=1,nfr
        transp=transp+plantup(i)
    ENDDO

!******************************************************    
!   water table module starts here
    vtot = MAX(145.,wsc(1)+wsc(2)+wsc(3)+infilt)!+wsc(4)+wsc(5)   ! total amount of water in top 500mm of soil  mm3/mm2 infilt here is standing water   infilt has distributed to wsc?
!   infilt means standing water according to jiangjiang
!   vtot = MAX(145.,vtot+145.+rain-evap-transp-runoff)         ! vtot should not be smaller than 145, which is the water content when wt is at -300mm
!   phi = 0.95   !wan-soil porosity   mm3/mm3   the same unit with theta, Granberg et al. (1999)
    zmax = 300   !maximum water table depth   mm
    thetasmin = 0.25    !minimum volumetric water content at the soil surface   cm3/cm3
    zthetasmin = 100    !maximum depth where evaporation influences soil moisture   mm
    az = (phi-thetasmin)/zthetasmin     ! gradient in soil moisture resulting from evaporation at the soil surface    mm-1
    
    zwt1 = -sqrt(3.0*(phi*zmax-vtot)/(2.0*az))
    zwt2 = -(3.0*(phi*zmax-vtot)/(2.0*(phi-thetasmin)))
    zwt3 = vtot-phi*zmax                                   
    IF ((zwt1 .ge. -100) .and. (zwt1 .LE. 0))   zwt = zwt1  !the non-linear part of the water table changing line
    IF (zwt2 .LT. -100)                         zwt = zwt2  !the linear part of the water table changing line

!   IF ((zwt2 .LT. -100) .and. (zwt2 .ge. -300))zwt = zwt2 !the linear part of the water table changing line valid when Vtot>145mm
!   IF (zwt2 .LE. -300)                         zwt = -300
    IF (phi*zmax .LT. vtot)                     zwt = zwt3  !the linear part when the water table is above the soil surface 
        
!   water table module ends here
!******************************************************    
    
!   Output fwsoil, omega, and topfws
    DO i=1,nfr       
        ome(i)=(wcl(i)-WILTPT)/(FILDCP-WILTPT)  !wan- (W_soil-W_min)/(W_max-W_min), wcl: relative soil water content
        ! wan 2023/3/21: wcl, unit % or mm3/mm3 or 1
        ome(i)=AMIN1(1.0,AMAX1(0.0,ome(i)))
        fw(i)=amin1(1.0,3.333*ome(i))           !wan- Weng and Luo 2008 (A12),omega -- normalized soil moisture
        fw(i)=amin1(1.0,2*ome(i))               !wan 0628 change this equation due to the water limitation
        !assume water limite plant growth when ome lower than 0.5
    ENDDO

    topfws=amin1(1.0,(wcl(1)-WILTPT)/((FILDCP-WILTPT))) 
	if(topfws.LT.0.0) topfws=0.000001

    fwsoil=0.0
    omega=0.0
    DO i=1,nfr
        fwsoil= fwsoil+fw(i)*frlen(i)    
        omega = omega+ome(i)*frlen(i)
    ENDDO
    fwsoil_rec(itime) = fwsoil
    omega_rec(itime) = omega
    RETURN
    END ! END soilwater subroutine
    
!**********************************************************************
!     plant growth model
SUBROUTINE plantgrowth(Tair,stor_use,SNgrowth,L_fall,&
    &                RgLeaf,RgStem,RgRoot,RgRe,count,upgraded)
     USE LIMITATION 
     USE IntersVariables
     USE vars_site
     IMPLICIT NONE
     REAL nsCN
     REAL SnscnL,SnscnS,SnscnR         !w: not use in original code
     REAL stor_use
     REAL tauLeaf
     REAL GrowthL,GrowthR,GrowthS
     REAL Tair
!     biomass
     REAL CNP0           !w: CNP0 - CNp0 - plant CN ratio
     REAL la0,GPmax,acP,c1,c2,b1,b2,b3,b4
     REAL QCbmL,QCbmR,QCbmP,QCbmS,QCbmRe
     REAL Rgroot,Rgleaf,Rgstem,RgRe
!     scalars
     REAL Ss,Sn,SL_rs,SR_rs,Slai,SNgrowth,phiN,SNgrowth_use
     REAL RSw
     REAL gamma_Wmax,gamma_Tmax,gamma_N
     REAL beta_T,Tcold,Twarm,Topt
     REAL bW,bT
     REAL L_fall,L_add,NL_fall,NL_add
     REAL alpha_St
     INTEGER i
     REAL SPgrowth,SPgrowth_use,CPp0  
     REAL LFactor
     REAL L,W2,k_n,eS,eR,eL,ww
     REAL xnpGrowthL,xnpGrowthR,xnpGrowthS,xnpGrowthRe
     REAL xnpGrowthL2,xnpGrowthR2,xnpGrowthS2,xnpGrowthRe2
     INTEGER count,upgraded

    Twarm=35.0       ! wan- Not use
    !Tcold=5.0
    !Tcold=0.0       ! For SPRUCE
    Tcold =5.0 !5.0  ! wan- change to 5 for TianTong， Aug 13
    Topt=30.         ! wan- Not use
    phiN=0.33

    QCbmL=bmleaf*fbmC ! Carbon
    QCbmR=bmRoot*fbmC
    QCbmS=bmStem*fbmC
    QCbmRe=bmRe*fbmC
    
    IF(QCbmL.LT.NSC/0.3)THEN
        QCbmL=NSC/0.3
    ENDIF
    IF(QCbmR.LT.NSC/0.3)QCbmR=NSC/0.3
    IF(QCbmS.LT.NSC/0.3)QCbmS=NSC/0.3
    IF(QCbmRe.LT.NSC/0.1)QCbmRe=NSC/0.1
    
    StemSap=SapS*QCbmS  ! Weng 12/05/2008
    RootSap=SapR*QCbmR
    IF(StemSap.LT.0.001)StemSap=0.001
    IF(RootSap.LT.0.001)RootSap=0.00111

    QCbmP=QCbmL+QCbmR+QCbmS+QCbmRe					! Plant C biomass 
    acP=QCbmL+StemSap+RootSap+0.5*QCbmRe 			! Plant available sapwood C  
    la0=0.2 
    ht=hmax*(1.-exp(-hl0*QCbmP))				  ! Scaling plant C biomass to height
    LAIMAX=AMAX1(LAIMAX0*(1.-exp(-la0*ht)),LAIMIN+0.1)  ! Scaling plant height to maximum LAI

!   Phenology
!   Only update the storage pool, and the phenology state
    gddonset = 488.
    IF((GDD5.GT.gddonset).and.onset.eq.0.and.storage.GT.stor_use) THEN  !wan- add the 'onset eq 0' 
        onset=1 
    ENDIF  
    
    IF((onset.eq.1).and.(storage.GT.stor_use))THEN 
        IF(LAI.LT.LAIMAX) add=stor_use      ! wan-2022.10: use the C from storage pool, if LAI .LT. than LAIMAX
        storage=storage-add
        !stop
    ELSE
        add=0.0
        onset=0
        ! add by zhoujian
    ENDIF
    !IF(accumulation.LT.(NSCmax*0.1+0.005*RootSap))THEN   !accumulation, initial value = 0
    IF(accumulation.LT.(NSCmax*0.4))THEN    !accumulation, 2024 01 23
        store=AMAX1(0.,0.005*NSC)			! 0.5% of nonstructure carbon is stored
    ELSE
        store=0.0  
    ENDIF
    accumulation = accumulation+store !wan- 'accumulation' help judge when to add carbon to store

    St = AMAX1(0.0, 1.0-exp(-(Tair-Tcold)/10))  ! wan 0320 change /15 to 10
!   Sw=AMIN1(0.5,AMAX1(0.333, 0.333+omega))      !Jiang modified June 2016
    Sw=omega  !wan- 20230205
    W = AMIN1(1.0,3.333*omega)                   !same equation in Weng & Luo 2008

    GPmax=(GLmax*QCbmL+GSmax*StemSap+GRmax*QCbmR+GRemax*QCbmRe)  !/acP 
    !GrowthP = AMIN1(GPmax*St*Sw*fnsc,0.04*NSC)  !wan- to exclude the N, same with results of above equation almostly
    GrowthP = GPmax*St*Sw*fnsc


    IF (CYCLE_CNP .GT. 1) THEN
        LFactor = AMIN1 (LFactorP,LFactorN)
        xnpGrowthP = GrowthP*LFactor*x_leaf_limit
        ExcessC = GrowthP - xnpGrowthP
    ELSE
        LFactor = 1.
        xnpGrowthP = GrowthP
        ExcessC =0.
    ENDIF

    xnpGrowthP_temp =xnpGrowthP
    ! For debug
    IF(GrowthP .GT. 200 .or. isnan(GrowthP)) THEN
        WRITE(*,*)'GrowthP,fnsc,nsc',GrowthP,fnsc,nsc
        Flag_Growth = 1
    ENDIF

!  Original code from zhoujian
!  Fixed allocation coefficience, wan- what's the rationale?????
!      GrowthL=MAX(0.0,GrowthP*0.4)  ! add by zhoujian 
      !GrowthR=MIN(GrowthP*0.3,MAX(0.0,0.75/Sw*QCbmL-QCbmR))  ! *c1/(1.+c1+c2)
!      GrowthR=MAX(0.0,GrowthP*0.3)
!      GrowthS=MAX(0.0,GrowthP - (GrowthL+GrowthR) )         ! *c2/(1.+c1+c2)
!      NPP = GrowthL + GrowthR + GrowthS + add       ! Modified by Jiang Jiang 2015/10/13

!  Luo et al., 1995, Denisom and Loomis 1989, Shevliakova et al., 2009, Weng and Luo, 2008
!  wan-2022.10 : for the flexible allocation coefficience, REF -- Du doctoral dissertation, P55
      
    !c1=QCbmL/QCbmR*CN(1)/CN0(1) !bm -- biomass of leaf or fine root
!          !c2 = 0  --- not use
      !c2=0.5*250e3*SLA*0.00021*ht*2
      !b1 = 1./(1.+c1+c2) !allocate to leaf
      !b2 = c2/(1.+c1+c2) !allocate to stem/wood
      !b3 = c1/(1.+c1+c2) !allocate to root

    !   GrowthL=GrowthP*b1
    !   GrowthS=GrowthP*b2
    !   GrowthR=GrowthP*b3
    !   NPP = GrowthL + GrowthR + GrowthS + add  ! +add / +store?  +add is more reasonalble
    
!--------Arora & Boer
    k_n = 0.5
    ww = 0.36  !origin: 0.8
    eL  = 0.35
    eS  = 0.05 
    eR  = 0.60
  
    IF(LAI .LE. LAIMIN)THEN
        a_Re=0.0
    ELSE
        a_Re=0.12 ! wan assume a fixed allocation fraction for fecundity 12% of new carbon
    ENDIF

    L = exp(-k_n*LAI)
    W2 = AMAX1(0.0, Amin1(1.,omega))
    a_S = (eS*x_stem_limit+ww*(1.0-L))/(1.+ww*(2.0-L-W2))
    a_L = (eL*x_leaf_limit)/(1.+ww*(2.-L-W2))
    a_R=1-a_L-a_S   

    ! 10% C allocate to reproduction  
    xnpGrowthRe=xnpGrowthP*a_Re
    xnpGrowthP = xnpGrowthP-xnpGrowthRe

    xnpGrowthL=xnpGrowthP*a_L 
    xnpGrowthS=xnpGrowthP*a_S 
    xnpGrowthR=xnpGrowthP*a_R
    
    !NPP_FORNSC = GrowthL + GrowthR + GrowthS + add   !GrowthP+add
    !xnpNPP = xnpGrowthL + xnpGrowthR + xnpGrowthS + add + xnpGrowthRe
    xnpNPP = (xnpGrowthL+add)+xnpGrowthR+xnpGrowthS+xnpGrowthRe
    NPP = xnpNPP
    NPP_noadd = NPP - add


    !wan- I just set x_leaf_limit,xuptake to Vcmax! try 2022/12/4
    !xnpNPP = x_leaf_limit*xuptake*NPP  !Note, see the code in casa_cnp, confused euqation!
    !NPP=xnpNPP    !2023/3/22

    ! This part is different from Wang 2010, AppendixB, wan- 2022/12/4, move after xnpNPP
    IF(NPP.eq.0.0)THEN
        !count = 1+count
        alpha_L=0.3
        alpha_W=0.3
        alpha_R=0.3
        alpha_Re=0.1
    ELSE
        alpha_L=(xnpGrowthL+add)/NPP       !wan 20230517 change GrowthL to xnpgrowthL
        alpha_W=xnpGrowthS/NPP
        alpha_R=xnpGrowthR/NPP
        !alpha_Re=xnpGrowthRe/NPP
        alpha_Re = 1-(alpha_L+alpha_W+alpha_R)
    ENDIF

!   Carbon cost for growth
!   Rgrowth,Rgroot,Rgleaf,Rgstem, 0.5 is from IBIS and Amthor, 1984
    !Rgleaf=0.5*GrowthL    ! 2024 02 26 xnpGrowth
    !Rgstem=0.5*GrowthS
    !Rgroot=0.5*GrowthR

    ! growth respiration cost per unit of new plant tissue carbon (GR ) 
    ! set as a constant for all tissues (GRc = 0.3; Larcher, 1995)
    ! Thornton & Rosenbloom 2005
    Rgleaf=0.2*xnpGrowthL       !change 0.3 to 0.2
    Rgstem=0.2*xnpGrowthS
    Rgroot=0.2*xnpGrowthR
    RgRe=0.2*xnpGrowthRe
    Rgrowth=Rgleaf+Rgstem+Rgroot+RgRe

!   Leaf litter 
!   Maximum rates of leaf fall induced by low T and drought respectively
    !gamma_Wmax=0.12/24. ! maxmum leaf fall rate per hour
    !gamma_Tmax=0.12/24.

    gamma_Wmax=0.005/24. ! maxmum leaf fall rate per hour
    gamma_Tmax=0.005/24.


    bW=4.0   
    bT=2.0
    beta_T = 1

    IF(Tair .LE. Tcold)THEN
        beta_T= 0.  
!Second stage, rapid growth, suitable temperature
    ELSEIF(Tair .GT. Tcold .and. Tair .LT. (Tairmax-5.) .and. itime .LT. Tairmax_Loc_value)THEN
        beta_T = 1.
!Third stage, keep，but a litter decrease because Too hot, 
    ELSEIF(Tair .GE. (Tairmax-5.).and. Tair .LE. Tairmax)THEN
        beta_T=(Tair-Tcold)/10.  !100
        beta_T = AMAX1(0.,AMIN1(1., beta_T))
!Fourth stage, prepare for winter 
    ELSEIF(Tair .GT. Tcold .and. Tair .LT. Tairmax .and. itime .GT. Tairmax_Loc_value )THEN
        beta_T=(Tair-Tcold)/20.
        beta_T = AMAX1(0.,AMIN1(1., beta_T))
    ENDIF



    !IF (tauC(1) < 8760.)THEN                     !wan- why add this condition?
        gamma_W = (1. - W)     **bW * gamma_Wmax  !gamma_W always be 1
        !gamma_W = (1. -  S_omega)**bW * gamma_Wmax  ! 2024 01 22 wan
        gamma_T = (1. - beta_T)  **bT * gamma_Tmax
    !ELSE
    !    gamma_W = 0.
    !    gamma_T = 0.
    !ENDIF                                         !wan- delete the condition according to  Arora & Boer 2005  
    !gamma_N = 1.0/(tauC(1)*Sw)   !gamma_N=1.0/Tau_L*Sw      ! Modify by Jiang Jiang 2015/10/20
                                 !tauC modified by Zhoujian? tauC = Tau_L
    gamma_N = 1.0/tauC(1) !wan- 2022/12/14 Ref Arora & Boer 2005.
    gamma_N = exitK(1)

    !since the unit of leaf loss in Arora & Boer is ... d-1, they multiple 365 there, we dont need this. 
    !keep them houly, QC(1)*gamma_N = QC(1)/tau(1)

    IF(LAI < LAIMIN) THEN
        gamma_W = 0.
        gamma_T = 0.
        !gamma_N = 0.   
    ENDIF

    ! L_fall = bmleaf*0.45*AMIN1((gamma_W+gamma_T+gamma_N),0.99) 
    ! L_fall=bmleaf*0.48*gamma_N
    L_tau_new = gamma_W+gamma_T+gamma_N !wan 20230708
    L_fall = QC(1)*L_tau_new!*(gamma_W+gamma_T+gamma_N) 
    RETURN
    END
!************************************************************************
! carbon transfer module
! wan- Wang & Luo. 2008
SUBROUTINE TCS_CNP(Tair,Tsoil,runoff,L_fall,&
    &               CNmin,CNmax,NSNmax,NSNmin,&        ! nitrogen
    &               SNvcmax,SNgrowth,SNRauto,SNrs,&
    &               P_leach) 
         
    !USE TransFraction !partitioning coefficients
    USE NPmin_imm
    USE NutrientsUptake
    USE NPleaching
    USE P_specific_process
    USE IntersVariables
    USE vars_site

    IMPLICIT NONE
 
    REAL L_fall,L_add
    REAL Tair,Tsoil,runoff
    REAL Q_plant
    REAL Q10h(5) ! Q10 of the litter and soil C pools
!     the fraction of C-flux which enters the atmosphere from the kth pool
    REAL f_CO2_fine,f_CO2_coarse,f_CO2_Micr,f_CO2_Slow,f_CO2_Pass
!---------------------------------------------------------------------------
!     for nitrogen sub-model
    REAL CNmin,CNmax,NSNmax,NSNmin,fnsn
    REAL CN_foliage
    REAL N_imm(5),N_imm2(3)
    REAL N_fixation_x,Nfix0
    REAL N_loss,ScNloss
    REAL Qroot0,Cfix0
    REAL Scalar_N_flow,Scalar_N_T
    REAL ksye,kappaVcmax
    REAL SNvcmax,SNgrowth,SNRauto,SNrs
    REAL SNleaf,SNwood,SNroot
    REAL SNfine,SNcoarse,SNmicr,SNslow,SNpass
    REAL costCuptake,costCfix,costCreuse
    REAL Creuse0,Nup0,N_deN0,LDON0
    REAL Nreqmax(3),Nreqmin(3)

!---------------------------------------------------------------------------
!     the variables relative to soil moisture calcualtion
    !REAL S_omega    !  average values of the moisture scaling functions
    !REAL S_t(5)     !  average values of temperature scaling functions
    !REAL S_w_min    !  minimum decomposition rate at zero available water
!     For test
    REAL totalC1,totalN1,C_in,C_out,N_in,N_out,totalC2,totalN2
!---------------------------------------------------------------------
    INTEGER i,j,k,n,m
    INTEGER day,week,month,year

!---------------------------------------------------------------------------

    REAL CP_foliage,QPplant,QNplant
    !REAL P_imm(5) ! immobilization P, for - QP(4,5,6,7,8)
    REAL P_imm2(3),P_imm(5)
    !REAL Fptase ! biochemical immobilization, for - QP(7,8)
    REAL SNleaf_use, SNwood_use, SNroot_use, SNfine_use
    REAL SNcoarse_use, SNmicr_use, SNslow_use, SNpass_use
    REAL xDsoillab,Dsoilsorb
       !P biochemical mineralization
    REAL prodptase !Phosphate production, biome-specific in CABLE
      ! same as u_pmax, the maximun specific biochemical P mineralization rate(d-1)
    REAL costNpup,Pup0 !N cost of P uptake 40gN/P for tropical biomes
      ! and 25gN/g P for other biomes
    REAL  P_biominer 
    REAL  Psorbmax,kmlabP,xkplab,xkpsorb,xkpocc !biome-specific,Parameters for Langmuir equation
    REAL  kplab,kpsorb,kpocc  !rate
    REAL  Pdep,Pwea,P_fer!P input
    REAL  Scalar_P_flow,LDOP0
    REAL  P_leach
    REAL  xkoptsoil,xksoil
    REAL  Preqmax(3),Preqmin(3)                     !CABLE
    REAL  xde_n,xde_p !N,P limitation on decomposition
    REAL  xde_n2,xde_p2,xde_n22,xde_p22
    REAL  MaxUptakeCapN,MaxUptakeCapP
    REAL  fnsp,f_PNplant,N_cost
    REAL  QC_o(8),QC_update(8),delta_QC_inProcess(8)
    REAL (SGL) :: delta_QC_TSC(9),delta_QN_TSC(9),delta_QP_TSC(9)
    REAL  P_in_out
    REAL  NfluxP2StrucL,NfluxP2MetL,CfluxP2L


    !Initial value for C cycle
    Q10h=(/Q10,Q10,Q10,Q10,Q10/)

    !Shared parameters of NP
    Qroot0=500.
    

    !Initial value/parameters for N cycle
    Nfix0=1./60.   ! maximum N fix ratio, N/C
    Nfix0=1./600.  ! J.Whiltshire 0.0016gN / gC
    Nup0 =0.02     ! nitrogen uptake rate  wan-the maximum rate of N absorption per step when Root_total
    ! approaches infinity
    Cfix0=12.      ! C cost per N for fixation
    ksye=0.05      ! C cost per N for uptake
    Creuse0=2.     ! C cost per N for resorption
    ScNloss=1    ! 1
    N_deN0=1.E-4*ScNloss   ! 1.E-3, 5.E-3, 10.E-3, 20.E-3
    LDON0=1.E-4*ScNloss
    
    Rnitrogen=0.
    !for N scalars
    CNmin=40.0
    CNmax=200.0
!   Max and min NSN pool
    !NSNmax = QN(1) + 0.2*QN(2) + QN(3)  ! 15.0
    NSNmax = 0.5*(QN(1) + QN(2)+ QN(3))  !wan 2023/3/22, value = 17.33
    NSNmin = 0.05
    !Initial value/parameters for P cycle
    !alphaP = 0.4    !YanEnRong, DoctoralDissertation, 2006
    ! P Input --- biome specific
    ! P weathering rate --- TianTong, soil orders: Ultisol, value = 0.005 gPm-2year-1
    ! Ref: [Wang et al., 2010, Table2]
    !Pwea = 0.005/8760 !covert gPm-2year-1 to gP m-2 hour-1， Ultisol,calibrated
    Pwea = 0.01/8760 ! wan20241224
    ! [Zhu Jianxing 2016] mean value of wet dep for forest ecosystem 
    ! *2 dry + wet 
    ! dry P deposition contributed 40–75% of the total atmospheric P deposition
    Pdep = 0.1/8760!6.5E-6!7.5E-6!3.28E-6*2!6.5E-6*1.4
    ! Initialize...
    P_miner  = 0.0
    P_imm2    = 0.0
    P_immob  = 0.0
    !P_net    = 0.0
    P_leach  = 0.0
    P_loss   = 0.0
    P_biominer= 0.0

    xde = 1.0
    xde_n = 1.0
    xde_p = 1.0
    xde2 = 1.0
    xde_n2 = 1.0
    xde_p2 = 1.0

    xde22 = 1.0
    xde_n22 = 1.0
    xde_p22 = 1.0
      
    ! Calculating soil scaling factors, S_omega and S_tmperature
    S_w_min=0.08 !0.08 !minimum decomposition rate at zero soil moisture  !0.03, 0522 wan 
    !S_omega=S_w_min + (1.-S_w_min) * Amin1(1.0,0.3*omega)   ! Jiang modified on June 2016
    S_omega=S_w_min + (1.-S_w_min) * Amin1(1.0,0.3*omega) !wan 2022/12/14s

    DO i=1,5
        S_t(i)=Q10h(i)**((Tsoil-10.)/10.)  !wan: use this equation according to GaoJie, 20230516
    ENDDO

!   Calculating NPP allocation and changes of each C pools
    NPP_L = alpha_L * NPP           ! NPP allocation, !wan 2023/3/17 change NPP to xnpNPP
    NPP_W = alpha_W * NPP
    NPP_R = alpha_R * NPP
    NPP_Re = alpha_Re * NPP

!Scalar for limitating on decomposition wan- 2022/12/5 Ref Wang et al., 2010

IF(CYCLE_CNP .GT. 1) THEN
    xde_n22 = AMIN1(AMAX1(0.001,1.0+(N_net22/QNminer)),1.0)
    IF(CYCLE_CNP .GT. 2) THEN
        xde_p22 = AMIN1(AMAX1(0.001,1.0+(P_net22/QPlab)),1.0)
    ENDIF
ENDIF

xde22= MIN(xde_n22,xde_p22) 
 
!     the carbon leaving the pools

    OutC(1)=L_fall!*SNleaf_use*SPleaf_use*xde_n*xde_p
    OutC(2)=QC(2)*exitK(2)
    OutC(3)=QC(3)*exitK(3)
    OutC(4)=QC(4)*exitK(4)
    OutC(5)=QC(5)*exitK(5)*S_omega* S_t(1)*xde22
    OutC(6)=QC(6)*exitK(6)*S_omega* S_t(2)*xde22
    OutC(7)=QC(7)*exitK(7)*S_omega* S_t(3)*xde22
    OutC(8)=QC(8)*exitK(8)*S_omega* S_t(4)*xde22
    OutC(9)=QC(9)*exitK(9)*S_omega* S_t(5)*xde22

!   This part for nutrient cycle
    OutC2(1:4) = OutC(1:4)
    OutC2(5)=QC(5)*exitK(5)*S_omega* S_t(1)
    OutC2(6)=QC(6)*exitK(6)*S_omega* S_t(2)
    OutC2(7)=QC(7)*exitK(7)*S_omega* S_t(3)
    OutC2(8)=QC(8)*exitK(8)*S_omega* S_t(4)
    OutC2(9)=QC(9)*exitK(9)*S_omega* S_t(5)
    
!-----------------------------------------


!     heterotrophic respiration from each pool
    Rh_pools(1)=OutC(5)* (1. - f_F2M) ! f_F2M -- partition parameter
    Rh_pools(2)=OutC(6)* (1. - f_C2M - f_C2S) 
    Rh_pools(3)=OutC(7)* (1. - f_M2S - f_M2P) 
    Rh_pools(4)=OutC(8)* (1. - f_S2P - f_S2M)
    Rh_pools(5)=OutC(9)* (1. - f_P2M)
    ! For debug
    IF(Rh_pools(1) .EQ. 0) THEN 
        WRITE(*,*)'Rh_pools,QNminer',sum(Rh_pools),QNminer
        IF(MCMC .EQ.0) THEN
            stop
        ELSE
            Flag_N = 1
        ENDIF
    ENDIF

    etaL2 = 1-etaL
    etaW2 = 1-etaW
    etaR2 = 1-etaR
    etaRe2= 1-etaRe

IF(CYCLE_CNP .GT. 1) THEN
    DO i=1,9
        OutN(i) = OutC(i)/CN(i)
    ENDDO
    DO i=1,9
        OutN2(i) = OutC2(i)/CN(i)
    ENDDO

    IF(CYCLE_CNP .GT. 2) THEN
    !OutP
    !Ref: TECO,nitrogen leaving the pools
    !OUTPUT: OutP(8), gPm-2hour-1
        DO i=1,9
            OutP(i) = OutC(i)/CP(i)            
        ENDDO
        DO i=1,9
            OutP2(i) = OutC2(i)/CP(i)            
        ENDDO
    ENDIF

!NP, MINERALIZATION
    CALL N_mineralization()

    IF(CYCLE_CNP .GT. 2.) THEN
    !P mineralization
    !Ref   : TECO,nitrogen mineralization
    !INPUT : OutP, transfer fraction
    !OUTPUT: P_miner
    !CABLE : Decomposition rate of litter/soil * P pools of litter/soil
        CALL P_mineralization()
    ENDIF

!NP, IMMOBILIZATION
    CALL N_immobilization(N_imm2)

    IF(CYCLE_CNP .GT. 2.)THEN
        CALL P_immobilization(P_imm2)
    ENDIF

    ! Before plant uptake, we close an count  of net_min
    ! 2024 01 25, wan

    N_immob2 = N_imm2(1)+N_imm2(2)+N_imm2(3)
    N_net2 = N_miner2 - N_immob2

!NP, UPTAKE 

    !IF (NutrientsUp .eq. 1.0) THEN
    CALL N_UptakeDemand(NSNmax,NSNmin,Nfix0,&
                        N_fixation_x,Tsoil)

    IF (NSN .LE. (0.001 * NSNmin))THEN
    ! if there is no nutrient support growth, there is no carbon allocate to tissue
        NPP_L = 0.0
        NPP_W = 0.0
        NPP_R = 0.0
        NPP_Re = 0.0
        NPP = 0.0
        CALL N_UptakeDemand(NSNmax,NSNmin,Nfix0,&
                            N_fixation_x,Tsoil)
    ENDIF

    !ELSEIF (NutrientsUp .eq. 2.0) THEN
    !CALL CABLE_Nuptake(NPP,NPP_L,NPP_W,NPP_R,alphaN,OutN,&
    !            CN,CN0,QC,&
    !            QNminer,Nreqmax,Nreqmin,NSN,N_uptake,&
    !            N_leaf,N_wood,N_root)
    !ENDIF

    IF(CYCLE_CNP .GT. 2.)THEN
        !IF (NutrientsUp .eq. 1.0) THEN
        CALL P_UptakeDemand() 
        IF(NSP .LE. (0.001 * NSPmin))THEN 
            NPP_L = 0.0
            NPP_W = 0.0
            NPP_R = 0.0
            NPP_Re = 0.0
            NPP = 0.0
            CALL N_UptakeDemand(NSNmax,NSNmin,Nfix0,&
                                N_fixation_x,Tsoil)
            CALL P_UptakeDemand() 
        ENDIF
    ENDIF

    QNminer=QNminer+N_deposit+N_fert+N_fixation_x+N_net2-N_uptake
    N_fixation = N_fixation_x

    CALL Nleaching(runoff,rdepth,N_deN0,Tsoil,LDON0,N_loss,&
                    & Scalar_N_flow,Scalar_N_T)
    !     update QNminer
    QNminer=QNminer-N_loss
!labile P dynamic
    IF(CYCLE_CNP .GT. 2) THEN
        CALL SorbedPDynamic(runoff,rdepth,QPplant,Pwea,Pdep,N_cost,P_fert)
        QNminer = QNminer - N_cost
        IF (QNminer .LT. 0) THEN
            Flag_N = 1
            write(*,*)'QNminer .LT. 0',QNminer
            !stop
        ENDIF
    ENDIF   
ENDIF !CYCLE_Nitrogen & Phosphorus

P_in_out = (Pwea+Fptase_rec(itime)+Pdep+P_net2_rec(itime))-&
    & (P_uptake_rec(itime)+asorb_rec(itime)+P_loss_rec(itime))

P_in_out_rec(itime) =  P_in_out

! update carbon pools, ! steply change of each pool size

    delta_QC_TSC(1) = - OutC(1) + NPP_L
    delta_QC_TSC(2) = - OutC(2) + NPP_W
    delta_QC_TSC(3) = - OutC(3) + NPP_R 
    delta_QC_TSC(4) = - OutC(4) + NPP_Re
    delta_QC_TSC(5) = - OutC(5) + etaL*OutC(1)+etaW*OutC(2)+etaR*OutC(3)+etaRe*OutC(4)
    delta_QC_TSC(6) = - OutC(6) + etaW2*OutC(2) &
                        + etaL2*OutC(1) + etaR2*OutC(3) & 
                        + etaRe2*OutC(4)
    delta_QC_TSC(7) = - OutC(7) + f_F2M*OutC(5)+f_C2M*OutC(6) &         
                        + f_S2M*OutC(8)+f_P2M*OutC(9)
    delta_QC_TSC(8) = - OutC(8)+f_C2S*OutC(6)+f_M2S*OutC(7)
    delta_QC_TSC(9) = - OutC(9)+f_M2P*OutC(7)+f_S2P*OutC(8)

    QC(1) = QC(1) + delta_QC_TSC(1)
    QC(2) = QC(2) + delta_QC_TSC(2)
    QC(3) = QC(3) + delta_QC_TSC(3)
    QC(4) = QC(4) + delta_QC_TSC(4)
    QC(5) = QC(5) + delta_QC_TSC(5)
    QC(6) = QC(6) + delta_QC_TSC(6)
    QC(7) = QC(7) + delta_QC_TSC(7)
    QC(8) = QC(8) + delta_QC_TSC(8)
    QC(9) = QC(9) + delta_QC_TSC(9)

    Q_plant =QC(1) + QC(2) + QC(3) + QC(4)      !wan, 2023/3/22 for NSN,NSP
    DO i=1,9
        IF(QC(i) .LT. 0.0)Flag_C = 1
    ENDDO

    IF(CYCLE_CNP > 1) THEN
    !   update nitrogen pools

        delta_QN_TSC(1) = - OutN(1) + N_leaf
        delta_QN_TSC(2) = - OutN(2) + N_wood
        delta_QN_TSC(3) = - OutN(3) + N_root
        delta_QN_TSC(4) = - OutN(4) + N_re
        delta_QN_TSC(5) = - OutN(5)  &
        &            + (etaL*OutN(1) + etaW*OutN(2) + etaR*OutN(3))*(1.0-alphaN) &   !*（1-alphaN）here
        &            + etaRe*OutN(4)
        delta_QN_TSC(6) =  - OutN(6) &
        &       + ((1.0-etaL)*OutN(1)+(1.0-etaW)*OutN(2)+(1.0-etaR)*OutN(3)) &
        &       * (1.0-alphaN) + (1.0-etaRe)*OutN(4)

        IF (CN(7) .GT. (1+0.05)*CN0(7))THEN
            N_loss_mic=0.0
        ELSE
            N_loss_mic=Scalar_N_flow*QN(7)*LDON0
        ENDIF
        delta_QN_TSC(7) = - OutN(7) + N_imm2(1)-N_loss_mic
        delta_QN_TSC(8) = - OutN(8) + N_imm2(2)
        delta_QN_TSC(9) = - OutN(9) + N_imm2(3)

        QN(1)=QN(1)+delta_QN_TSC(1)
        QN(2)=QN(2)+delta_QN_TSC(2) 
        QN(3)=QN(3)+delta_QN_TSC(3) 
        QN(4)=QN(4)+delta_QN_TSC(4) 
        QN(5)=QN(5)+delta_QN_TSC(5) 
        QN(6)=QN(6)+delta_QN_TSC(6) 
        !Divide the flow from plant to structural and metabolic litter
        !CfluxP2L = (OutC(1)+OutC(2)+OutC(3)+OutC(4))*(1.0-alphaN)
        !NfluxP2StrucL = CfluxP2L*(1.0/CN(6))
        !NfluxP2MetL = CfluxP2L - NfluxP2StrucL
        !QN(5)=QN(5) - OutN(5) + NfluxP2MetL
        !QN(6)=QN(6) - OutN(6) + NfluxP2StrucL 
        QN(7)=QN(7)+delta_QN_TSC(7) 
        QN(8)= QN(8)+delta_QN_TSC(8) 
        QN(9)= QN(9)+delta_QN_TSC(9) 
        
        CN=QC/QN
        QNplant =QN(1) + QN(2) + QN(3) + QN(4) 
        CN_foliage=(QC(1)+QC(3))/(QN(1)+QN(3))
        
    ENDIF  
    !Update organic P pools
    IF(CYCLE_CNP >2) THEN
        delta_QP_TSC(1) = - OutP(1) + P_leaf
        delta_QP_TSC(2) = - OutP(2) + P_wood
        delta_QP_TSC(3) = - OutP(3) + P_root
        delta_QP_TSC(4) = - OutP(4) + P_re
        delta_QP_TSC(5) = - OutP(5) & !alphaP, the transfer of N before littering
        &       + (etaL*OutP(1) + etaW*OutP(2) + etaR*OutP(3))*(1.0-alphaP) & 
        &       + etaRe*OutP(4)
        delta_QP_TSC(6) =  - OutP(6) &
        &       + ((1.0-etaL)*OutP(1)+(1.0-etaW)*OutP(2)+(1.0-etaR)*OutP(3)) &
        &       * (1.0-alphaP) + (1.0-etaRe)*OutP(4)
        delta_QP_TSC(7) = - OutP(7) + P_imm2(1)   

        delta_QP_TSC(8)=-OutP(8)+P_imm2(2)-Fptase_slow
        delta_QP_TSC(9)=-OutP(9)+P_imm2(3)-Fptase_pass
        

        QP(1) = QP(1)+delta_QP_TSC(1)
        QP(2) = QP(2)+delta_QP_TSC(2)
        QP(3) = QP(3)+delta_QP_TSC(3) 
        QP(4) = QP(4)+delta_QP_TSC(4) 
        QP(5) = QP(5)+delta_QP_TSC(5) 
        QP(6) = QP(6)+delta_QP_TSC(6) 
        QP(7) = QP(7)+delta_QP_TSC(7) 
        QP(8) = QP(8)+delta_QP_TSC(8)  !producted the phosphate
        QP(9) = QP(9)+delta_QP_TSC(9)  !producted the phosphate
        
    !---------------------------------------------------------
    !update C/P ratio
        CP=QC/QP
        QPplant = QP(1) + QP(2)+ QP(3)+QP(4)
        CP_foliage=(QP(1)+QP(3))/(QP(1)+QP(3))
        DO i=1,9
            IF(QP(i) .LT. 0.0)Flag_QP =1
        ENDDO 
    ENDIF 

END  SUBROUTINE
!------------------------------------------------------------------------


SUBROUTINE xlayers(Tair,Dair,radabv,G,Esoil,fbeam,eairP,  &
     &             wind,co2ca,fwsoil,wcl,FLAIT,coszen,idoy,hours,   &
     &           tauL,rhoL,rhoS,xfang,extkd,extkU,wleaf,            &
     &           Rconst,sigma,emleaf,emsoil,theta,a1,Ds0,           &
     &           cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,         &
     &           gsw0,alpha,stom_n,wsmax,wsmin,                     &
     &           Vcmx0,eJmx0,conKc0,conKo0,Ekc,Eko,o2ci,            &
     &           Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,         &
     &           extKb,Rsoilabs,Hsoil,Acan1,Acan2,Ecan1,Ecan2,      &
     &           RnStL,QcanL,RcanL,AcanL,EcanL,HcanL,GbwcL,GswcL,gddonset,itime,&
     &           Aleaf,Tlk1,Tlk2,fslt,fshd,                         &
     &           XX_rec,XX_med_rec,gscx_rec,gscx_med_rec,EqStomata)


!    the multi-layered canopy model developed by 
!    Ray Leuning with the new radiative transfer scheme   
!    implemented by Y.P. Wang (from Sellers 1986)
!    12/Sept/96 (YPW) correction for mean surface temperature of sunlit
!    and shaded leaves
!    Tleaf,i=sum{Tleaf,i(n)*fslt*Gaussw(n)}/sum{fslt*Gaussw(n)} 
!    
      REAL Gaussx(5),Gaussw(5)
      REAL layer1(5),layer2(5)
      REAL tauL(3),rhoL(3),rhoS(3),Qabs(3,2),radabv(2),Rnstar(2)
      REAL Aleaf(2),Eleaf(2),Hleaf(2),Tleaf(2),co2ci(2)
      REAL gbleaf(2),gsleaf(2),QSabs(3,2),Qasoil(2)
      INTEGER ng,nw
      REAL rhoc(3,2),reff(3,2),kpr(3,2),scatt(2)       !Goudriaan

      REAL rsoil,rlai,raero
      REAL wsmax,wsmin,WILTPT,FILDCP,wcl(10)
      REAL gddonset
!    additional arrays to allow output of info for each Layer
      REAL RnStL(5),QcanL(5),RcanL(5),AcanL(5),EcanL(5),HcanL(5)
      REAL GbwcL(5),GswcL(5)
      REAL XX_rec,XX_med_rec,gscx_rec,gscx_med_rec
      INTEGER EqStomata

      !INTEGER itime
      
! Normalised Gaussian points and weights (Goudriaan & van Laar, 1993, P98)
!* 5-point
      data Gaussx/0.0469101,0.2307534,0.5,0.7692465,0.9530899/
      data Gaussw/0.1184635,0.2393144,0.2844444,0.2393144,0.1184635/

!     soil water conditions
      WILTPT=wsmin/100.
      FILDCP=wsmax/100.
!     reset the vairables
      Rnst1=0.0        !net rad, sunlit
      Rnst2=0.0        !net rad, shaded
      Qcan1=0.0        !vis rad
      Qcan2=0.0
      Rcan1=0.0        !NIR rad
      Rcan2=0.0
      Acan1=0.0        !CO2
      Acan2=0.0
      Ecan1=0.0        !Evap
      Ecan2=0.0
      Hcan1=0.0        !Sens heat
      Hcan2=0.0
      Gbwc1=0.0        !Boundary layer conductance
      Gbwc2=0.0
      Gswc1=0.0        !Canopy conductance
      Gswc2=0.0
      Tleaf1=0.0       !Leaf Temp
      Tleaf2=0.0  
	  
!     aerodynamic resistance                                                
      raero=50./wind    !wan- aerodynamic resistance between the ground and the canopy air space                   
      !raero=5./wind    !wan- changed raero, to reduce the GPP when high weed speed, 2023
!     Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13; Wang & Leuning 1998 B7)
      xphi1 = 0.5 - 0.633*xfang -0.33*xfang*xfang
      xphi2 = 0.877 * (1.0 - 2.0*xphi1)
      funG=xphi1 + xphi2*coszen                             !G-function: Projection of unit leaf area in direction of beam
      
      IF(coszen.GT.0) THEN                                  !check IF day or night
        extKb=funG/coszen                                   !beam extinction coeff - black leaves
      ELSE
        extKb=100.
      END IF

!     Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
!     Effective extinction coefficient for diffuse radiation Goudriaan & van Laar Eq 6.6)
      pi180=3.1416/180.
      cozen15=cos(pi180*15)
      cozen45=cos(pi180*45)
      cozen75=cos(pi180*75)
      xK15=xphi1/cozen15+xphi2
      xK45=xphi1/cozen45+xphi2
      xK75=xphi1/cozen75+xphi2
      transd=0.308*exp(-xK15*FLAIT)+0.514*exp(-xK45*FLAIT)+     &
     &       0.178*exp(-xK75*FLAIT)
      extkd=(-1./FLAIT)*alog(transd)
      extkn=extkd                        !N distribution coeff 

!canopy reflection coefficients (Array indices: first;  1=VIS,  2=NIR
!                                               second; 1=beam, 2=diffuse
      DO nw=1,2                                                      !nw:1=VIS, 2=NIR
       scatt(nw)=tauL(nw)+rhoL(nw)                      !scattering coeff
       IF((1.-scatt(nw))<0.0)scatt(nw)=0.9999           ! Weng 10/31/2008
       kpr(nw,1)=extKb*sqrt(1.-scatt(nw))               !modified k beam scattered (6.20)
       kpr(nw,2)=extkd*sqrt(1.-scatt(nw))               !modified k diffuse (6.20)
       rhoch=(1.-sqrt(1.-scatt(nw)))/(1.+sqrt(1.-scatt(nw)))            !canopy reflection black horizontal leaves (6.19)
       rhoc15=2.*xK15*rhoch/(xK15+extkd)                                !canopy reflection (6.21) diffuse
       rhoc45=2.*xK45*rhoch/(xK45+extkd)
       rhoc75=2.*xK75*rhoch/(xK75+extkd)
       rhoc(nw,2)=0.308*rhoc15+0.514*rhoc45+0.178*rhoc75
       rhoc(nw,1)=2.*extKb/(extKb+extkd)*rhoch                          !canopy reflection (6.21) beam 
       reff(nw,1)=rhoc(nw,1)+(rhoS(nw)-rhoc(nw,1))   &                   !effective canopy-soil reflection coeff - beam (6.27)
     &            *exp(-2.*kpr(nw,1)*FLAIT) 
       reff(nw,2)=rhoc(nw,2)+(rhoS(nw)-rhoc(nw,2))   &                   !effective canopy-soil reflection coeff - diffuse (6.27)
     &            *exp(-2.*kpr(nw,2)*FLAIT)  
      ENDDO


!     isothermal net radiation & radiation conductance at canopy top - needed to calc emair
      CALL Radiso(flait,flait,Qabs,extkd,Tair,eairP,cpair,Patm, &
     &            fbeam,airMa,Rconst,sigma,emleaf,emsoil,       &
     &            emair,Rnstar,grdn)
      TairK=Tair+273.2

!     below      
      DO ng=1,5
         flai=gaussx(ng)*FLAIT
!        radiation absorption for visible and near infra-red
         CALL goudriaan(FLAI,coszen,radabv,fbeam,reff,kpr,      &   
     &                  scatt,xfang,Qabs) 
!        isothermal net radiation & radiation conductance at canopy top
         CALL Radiso(flai,flait,Qabs,extkd,Tair,eairP,cpair,Patm,   &
     &               fbeam,airMa,Rconst,sigma,emleaf,emsoil,        &
     &               emair,Rnstar,grdn)
         windUx=wind*exp(-extkU*flai)               !windspeed at depth xi !
         scalex=exp(-extkn*flai)                    !scale Vcmx0 & Jmax0   !wan- scaling up for single leaf to big leaves
         Vcmxx=Vcmx0*scalex
         eJmxx=eJmx0*scalex

         !IF(radabv(1).ge.10.0) THEN                          !check solar Radiation > 10 W/m2
         IF(radabv(1).ge.20.0) THEN      !wan 20240616 change according to the situation in TT
!           leaf stomata-photosynthesis-transpiration model - daytime
            CALL agsean_day(Qabs,Rnstar,grdn,windUx,Tair,Dair,      &
     &               co2ca,wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,  &
     &               Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,  &
     &               gsw0,alpha,stom_n,                                 &
     &               Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,            &
     &               Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,         &
     &               Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci,gddonset,itime,&
     &               XX_rec,XX_med_rec,gscx_rec,gscx_med_rec,Dleaf,EqStomata)    

         ELSE

            CALL agsean_ngt(Qabs,Rnstar,grdn,windUx,Tair,Dair,      &
     &               co2ca,wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,  &
     &               Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,  &
     &               gsw0,alpha,stom_n,                                 &
     &               Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,            &
     &               Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,         &
     &               Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci,Dleaf,&
     &               XX_rec,XX_med_rec,gscx_rec,gscx_med_rec)
         ENDIF  

         fslt=exp(-extKb*flai)                        !fraction of sunlit leaves
         fshd=1.0-fslt                                !fraction of shaded leaves
         IF (fslt .LE. 0.1) fshd=(0.1)-fslt     ! wan add 20240627
         Rnst1=Rnst1+fslt*Rnstar(1)*Gaussw(ng)*FLAIT  !Isothermal net rad`
         Rnst2=Rnst2+fshd*Rnstar(2)*Gaussw(ng)*FLAIT
         RnstL(ng)=Rnst1+Rnst2
!
         Qcan1=Qcan1+fslt*Qabs(1,1)*Gaussw(ng)*FLAIT  !visible
         Qcan2=Qcan2+fshd*Qabs(1,2)*Gaussw(ng)*FLAIT
         QcanL(ng)=Qcan1+Qcan2
!
         Rcan1=Rcan1+fslt*Qabs(2,1)*Gaussw(ng)*FLAIT  !NIR
         Rcan2=Rcan2+fshd*Qabs(2,2)*Gaussw(ng)*FLAIT
         RcanL(ng)=Rcan1+Rcan2
!
         IF(Aleaf(1).LT.0.0)Aleaf(1)=0.0      !Weng 2/16/2006
         IF(Aleaf(2).LT.0.0)Aleaf(2)=0.0      !Weng 2/16/2006
         Acan1=Acan1+fslt*Aleaf(1)*Gaussw(ng)*FLAIT*stom_n    !amphi/hypostomatous
         Acan2=Acan2+fshd*Aleaf(2)*Gaussw(ng)*FLAIT*stom_n
         AcanL(ng)=Acan1+Acan2

         layer1(ng)=Aleaf(1)
         layer2(ng)=Aleaf(2)

         Ecan1=Ecan1+fslt*Eleaf(1)*Gaussw(ng)*FLAIT
         Ecan2=Ecan2+fshd*Eleaf(2)*Gaussw(ng)*FLAIT
         EcanL(ng)=Ecan1+Ecan2
!
         Hcan1=Hcan1+fslt*Hleaf(1)*Gaussw(ng)*FLAIT
         Hcan2=Hcan2+fshd*Hleaf(2)*Gaussw(ng)*FLAIT
         HcanL(ng)=Hcan1+Hcan2
!
         Gbwc1=Gbwc1+fslt*gbleaf(1)*Gaussw(ng)*FLAIT*stom_n
         Gbwc2=Gbwc2+fshd*gbleaf(2)*Gaussw(ng)*FLAIT*stom_n
!
         Gswc1=Gswc1+fslt*gsleaf(1)*Gaussw(ng)*FLAIT*stom_n
         Gswc2=Gswc2+fshd*gsleaf(2)*Gaussw(ng)*FLAIT*stom_n
!
         Tleaf1=Tleaf1+fslt*Tleaf(1)*Gaussw(ng)*FLAIT
         Tleaf2=Tleaf2+fshd*Tleaf(2)*Gaussw(ng)*FLAIT
      ENDDO  ! 5 layers

      FLAIT1=(1.0-exp(-extKb*FLAIT))/extkb
      Tleaf1=Tleaf1/FLAIT1
      Tleaf2=Tleaf2/(FLAIT-FLAIT1)

!     Soil surface energy and water fluxes
!    Radiation absorbed by soil
      Rsoilab1=fbeam*(1.-reff(1,1))*exp(-kpr(1,1)*FLAIT)        &
     &         +(1.-fbeam)*(1.-reff(1,2))*exp(-kpr(1,2)*FLAIT)          !visible
      Rsoilab2=fbeam*(1.-reff(2,1))*exp(-kpr(2,1)*FLAIT)        &
     &         +(1.-fbeam)*(1.-reff(2,2))*exp(-kpr(2,2)*FLAIT)          !NIR
      Rsoilab1=Rsoilab1*radabv(1)
      Rsoilab2=Rsoilab2*radabv(2)
!  
      Tlk1=Tleaf1+273.2
      Tlk2=Tleaf2+273.2

!      temp1=-extkd*FLAIT
      QLair=emair*sigma*(TairK**4)*exp(-extkd*FLAIT)
      QLleaf=emleaf*sigma*(Tlk1**4)*exp(-extkb*FLAIT)           &
     &      +emleaf*sigma*(Tlk2**4)*(1.0-exp(-extkb*FLAIT))
      QLleaf=QLleaf*(1.0-exp(-extkd*FLAIT)) 
      QLsoil=emsoil*sigma*(TairK**4)
      Rsoilab3=(QLair+QLleaf)*(1.0-rhoS(3))-QLsoil

!    Net radiation absorbed by soil
!    the old version of net long-wave radiation absorbed by soils 
!    (with isothermal assumption)
!     Rsoil3=(sigma*TairK**4)*(emair-emleaf)*exp(-extkd*FLAIT)         !Longwave
!     Rsoilab3=(1-rhoS(3))*Rsoil3

!    Total radiation absorbed by soil    
      Rsoilabs=Rsoilab1+Rsoilab2+Rsoilab3 

!    thermodynamic parameters for air
      TairK=Tair+273.2
      rhocp=cpair*Patm*AirMa/(Rconst*TairK)     !wan- rhocp=c_p? the specific heat of the air?
      H2OLv=H2oLv0-2.365e3*Tair
      slope=(esat(Tair+0.1)-esat(Tair))/0.1
      psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar=Patm/(Rconst*TairK)
      fw1=AMIN1(AMAX1((FILDCP-wcl(1))/(FILDCP-WILTPT),0.05),1.0)  !wan- why use FILDCP-wcl(1)?
      Rsoil=30.*exp(0.2/fw1)
      !Rsoil=exp(8.206-4.255*wcl(1))  !similar results with above equation, !wan- add, Sellers,1996
      rLAI=exp(FLAIT)
!     latent heat flux into air from soil
!           Eleaf(ileaf)=1.0*
!     &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    !2* Weng 0215
!     &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
      Esoil=(slope*(Rsoilabs-G)+rhocp*Dair/(raero+rLAI))/       & !wan- soil Latent heat, why use this equation???
     &      (slope+psyc*(rsoil/(raero+rLAI)+1.))
!     sensible heat flux into air from soil
      Hsoil=Rsoilabs-Esoil-G

      RETURN
END 
SUBROUTINE goudriaan(FLAI,coszen,radabv,fbeam,reff,kpr,   &
     &                  scatt,xfang,Qabs)
     
!    for spheric leaf angle distribution only
!    compute within canopy radiation (PAR and near infra-red bands)
!    using two-stream approximation (Goudriaan & vanLaar 1994)
!    tauL: leaf transmittance
!    rhoL: leaf reflectance
!    rhoS: soil reflectance
!    sfang XiL function of Ross (1975) - allows for departure from spherical LAD
!         (-1 vertical, +1 horizontal leaves, 0 spherical)
!    FLAI: canopy leaf area index
!    funG: Ross' G function
!    scatB: upscatter parameter for direct beam
!    scatD: upscatter parameter for diffuse
!    albedo: single scattering albedo
!    output:
!    Qabs(nwave,type), nwave=1 for visible; =2 for NIR,
!                       type=1 for sunlit;   =2 for shaded (W/m2)

      REAL radabv(2)
      REAL Qabs(3,2),reff(3,2),kpr(3,2),scatt(2)
      xu=coszen                                         !cos zenith angle
      
!     Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
      xphi1 = 0.5 - 0.633*xfang -0.33*xfang*xfang      !xfang = x^2 Weng&Leuning 1998, B7
      xphi2 = 0.877 * (1.0 - 2.0*xphi1)
      funG=xphi1 + xphi2*xu                             !G-function: Projection of unit leaf area in direction of beam
      
      IF(coszen.GT.0) THEN                                  !check IF day or night
        extKb=funG/coszen                                   !beam extinction coeff - black leaves
      ELSE
        extKb=100.
      END IF
                       
! Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
      DO nw=1,2
       Qd0=(1.-fbeam)*radabv(nw)                                          !diffuse incident radiation
       Qb0=fbeam*radabv(nw)                                               !beam incident radiation
       Qabs(nw,2)=Qd0*(kpr(nw,2)*(1.-reff(nw,2))*exp(-kpr(nw,2)*FLAI))+  & !absorbed radiation - shaded leaves, diffuse
     &            Qb0*(kpr(nw,1)*(1.-reff(nw,1))*exp(-kpr(nw,1)*FLAI)-   & !beam scattered
     &            extKb*(1.-scatt(nw))*exp(-extKb*FLAI))
       Qabs(nw,1)=Qabs(nw,2)+extKb*Qb0*(1.-scatt(nw))                     !absorbed radiation - sunlit leaves 
      END DO
      RETURN
END

!****************************************************************************
SUBROUTINE Radiso(flai,flait,Qabs,extkd,Tair,eairP,cpair,Patm,    &
     &                  fbeam,airMa,Rconst,sigma,emleaf,emsoil,         &
     &                  emair,Rnstar,grdn)
!     output
!     Rnstar(type): type=1 for sunlit; =2 for shaded leaves (W/m2)
!     23 Dec 1994
!     calculates isothermal net radiation for sunlit and shaded leaves under clear skies
!     implicit REAL (a-z)
      REAL Rnstar(2)
      REAL Qabs(3,2)
      TairK=Tair+273.2

! thermodynamic properties of air
      rhocp=cpair*Patm*airMa/(Rconst*TairK)   !volumetric heat capacity (J/m3/K) wan- cpair

! apparent atmospheric emissivity for clear skies (Brutsaert, 1975)
      emsky=0.642*(eairP/Tairk)**(1./7)       !note eair in Pa
      emsky=0.533*(eairP/100)**(1./7)          !wan- e0 is in Pa, should change to millibars, (Brutsaert, 1975), equation 11

! apparent emissivity from clouds (Kimball et al 1982)
      ep8z=0.24+2.98e-12*eairP*eairP*exp(3000/TairK)     !wan. eq. 5, 6a, 6b
      tau8=amin1(1.0,1.0-ep8z*(1.4-0.4*ep8z))            !ensure tau8<1
      emcloud=0.36*tau8*(1.-fbeam)*(1-10./TairK)**4      !10 from Tcloud = Tair-10   !wan- cant fine the reference equation.

! apparent emissivity from sky plus clouds      
!      emair=emsky+emcloud
! 20/06/96
      emair=emsky

      IF(emair.GT.1.0) emair=1.0
      
! net isothermal outgoing longwave radiation per unit leaf area at canopy
! top & thin layer at flai (Note Rn* = Sn + Bn is used rather than Rn* = Sn - Bn in Leuning et al 1985)
      Bn0=sigma*(TairK**4.)
      Bnxi=Bn0*extkd*(exp(-extkd*flai)*(emair-emleaf)       &
     &    + exp(-extkd*(flait-flai))*(emsoil-emleaf))
!     isothermal net radiation per unit leaf area for thin layer of sunlit and shaded leaves
      Rnstar(1)=Qabs(1,1)+Qabs(2,1)+Bnxi
      Rnstar(2)=Qabs(1,2)+Qabs(2,2)+Bnxi
      !wan, change '+Bnxi' to '-Bnxi' 20240625
      !Rnstar(1)=Qabs(1,1)+Qabs(2,1)-Bnxi 
      !Rnstar(2)=Qabs(1,2)+Qabs(2,2)-Bnxi
      !wan test 2
      IF (Rnstar(1) .LT. 0) Rnstar(1) = 0.0
      IF (Rnstar(2) .LT. 0) Rnstar(2) = 0.0
!     radiation conductance (m/s) @ flai
      grdn=4.*sigma*(TairK**3.)*extkd*emleaf*               &       ! corrected by Jiang Jiang 2015/9/29
     &    (exp(-extkd*flai)+exp(-extkd*(flait-flai)))       &       ! wan- coorect what? same as TECO.for
     &    /rhocp
      RETURN
      END

!***********************************************************************
SUBROUTINE agsean_day(Qabs,Rnstar,grdn,windUx,Tair,Dair,      &
     &               co2ca,wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,  &
     &               Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,  &
     &               gsw0,alpha,stom_n,                                 &
     &               Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,            &
     &               Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,         &
     &               Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci,gddonset,itime, &
     &               XX_rec,XX_med_rec,gscx_rec,gscx_med_rec,Dleaf,EqStomata)

!    implicit REAL (a-z)
      INTEGER kr1,ileaf
      REAL Aleaf(2),Eleaf(2),Hleaf(2),Tleaf(2),co2ci(2)
      REAL gbleaf(2), gsleaf(2)
      REAL Qabs(3,2),Rnstar(2)
      REAL XX_rec,XX_med_rec,gscx_rec,gscx_med_rec
      INTEGER EqStomata
!    thermodynamic parameters for air
      TairK=Tair+273.2
      rhocp=cpair*Patm*AirMa/(Rconst*TairK)
      H2OLv=H2oLv0-2.365e3*Tair
      slope=(esat(Tair+0.1)-esat(Tair))/0.1
      psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar=Patm/(Rconst*TairK)
      weighJ=1.0
!    boundary layer conductance for heat - single sided, forced convection
!    (Monteith 1973, P106 & notes dated 23/12/94)
      IF(windUx/wleaf>=0.0)THEN
          gbHu=0.003*sqrt(windUx/wleaf)    !m/s
      ELSE
          gbHu=0.003 !*sqrt(-windUx/wleaf)
      ENDIF         ! Weng 10/31/2008
!     raero=0.0                        !aerodynamic resistance s/m
      DO ileaf=1,2              ! loop over sunlit and shaded leaves
!        first estimate of leaf temperature - assume air temp
         Tleaf(ileaf)=Tair
         Tlk=Tleaf(ileaf)+273.2    !Tleaf to deg K
!        first estimate of deficit at leaf surface - assume Da
         IF (EqStomata .EQ. 2) Dair = max(Dair,50*0.001)
         Dleaf=Dair                !Pa
!        first estimate for co2cs
         co2cs=co2ca               !mol/mol
         Qapar = (4.6e-6)*Qabs(1,ileaf)
!    ********************************************************************
         kr1=0                     !iteration counter for LE
!        RETURN point for evaporation iteration
         DO               !iteration for leaf temperature
!          single-sided boundary layer conductance - free convection (see notes 23/12/94)
           Gras=1.595e8*ABS(Tleaf(ileaf)-Tair)*(wleaf**3.)     !Grashof
           gbHf=0.5*Dheat*(Gras**0.25)/wleaf
           gbH=gbHu+gbHf                         !m/s
           rbH=1./gbH                            !b/l resistance to heat transfer
           rbw=0.93*rbH                          !b/l resistance to water vapour
!          Y factor for leaf: stom_n = 1.0 for hypostomatous leaf;  stom_n = 2.0 for amphistomatous leaf
           rbH_L=rbH*stom_n/2.                   !final b/l resistance for heat  
           rrdn=1./grdn
           Y=1./(1.+ (rbH_L+raero)/rrdn)
!          boundary layer conductance for CO2 - single side only (mol/m2/s)
           gbc=Cmolar*gbH/1.32            !mol/m2/s
           gsc0=gsw0/1.57                 !convert conductance for H2O to that for CO2
           varQc=0.0
           weighR=1.0
           CALL photosyn(co2ca,CO2Cs,Dleaf,Tlk,Qapar,Gbc,   &   !Qaparx<-Qapar,Gbcx<-Gsc0
                &         a1,Ds0,fwsoil,varQc,weighR,                  &
                &         gsc0,alpha,Vcmxx,eJmxx,weighJ,               &   !wan- gsc0-g0
                &         Aleafx,gscx,gddonset,itime,                  &
                &         XX_rec,XX_med_rec,gscx_rec,gscx_med_rec,Dair,EqStomata)  !outputs

!          choose smaller of Ac, Aq
           Aleaf(ileaf) = Aleafx      !0.7 Weng 3/22/2006          !mol CO2/m2/s
		   !          calculate new values for gsc, cs (Lohammer model)
           co2cs = co2ca-Aleaf(ileaf)/gbc
           co2Ci(ileaf) = co2cs-Aleaf(ileaf)/gscx
!          scale variables
!           gsw=gscx*1.56      !gsw in mol/m2/s, oreginal:gsw=gsc0*1.56,Weng20060215
           gsw=gscx*1.56       !gsw in mol/m2/s, oreginal:gsw=gscx*1.56,Weng20090226
           gswv=gsw/Cmolar                           !gsw in m/s
           rswv=1./gswv
!          calculate evap'n using combination equation with current estimate of gsw
           Eleaf(ileaf)=1.0*                                    &
     &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    &   !2* Weng 0215
     &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero)) 
!          calculate sensible heat flux
           Hleaf(ileaf)=Y*(Rnstar(ileaf)-Eleaf(ileaf))
!          calculate new leaf temperature (K)
           Tlk1=273.2+Tair+Hleaf(ileaf)*(rbH/2.+raero)/rhocp

!          calculate Dleaf use LE=(rhocp/psyc)*gsw*Ds
           Dleaf=psyc*Eleaf(ileaf)/(rhocp*gswv)    !wan- as input for phosyn module
           ! Dleafx
           gbleaf(ileaf)=gbc*1.32*1.075
           gsleaf(ileaf)=gsw
!          compare current and previous leaf temperatures

           IF(abs(Tlk1-Tlk) .GT. 3)Then            !Wan 2022/12/2
                Tlk=TairK
                exit
            ENDIF

           IF(abs(Tlk1-Tlk).LE.0.1) exit ! original is 0.05 C Weng 10/31/2008
!          update leaf temperature  ! leaf temperature calculation has many problems! Weng 10/31/2008
           Tlk=Tlk1
           Tleaf(ileaf)=Tlk1-273.2
           kr1=kr1+1
           IF(kr1 > 500)THEN
               Tlk=TairK
               exit
           ENDIF
           !IF(Tlk < 200.)THEN   !wan- wrong? kr1<200 Tlk<200  Matter!!!!! change TLK to kr1! nonono!
           !     Tlk=TairK
           !     exit 
           !ENDIF                     ! Weng 10/31/2008
!        goto 100                    ! solution not found yet
         ENDDO

      ENDDO
      RETURN
      END
!     ****************************************************************************
SUBROUTINE agsean_ngt(Qabs,Rnstar,grdn,windUx,Tair,Dair,co2ca,    &
     &               wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,            &
     &               Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,      &
     &               gsw0,alpha,stom_n,                                     &
     &               Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,                &
     &               Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,             &
     &               Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci,Dleaf,     &
     &               XX_rec,XX_med_rec,gscx_rec,gscx_med_rec)
!    implicit REAL (a-z)
      INTEGER kr1,ileaf
      REAL Aleaf(2),Eleaf(2),Hleaf(2),Tleaf(2),co2ci(2)
      REAL gbleaf(2), gsleaf(2)
      REAL Qabs(3,2),Rnstar(2)
!    thermodynamic parameters for air
      TairK=Tair+273.2
      rhocp=cpair*Patm*AirMa/(Rconst*TairK)
      H2OLv=H2oLv0-2.365e3*Tair
      slope=(esat(Tair+0.1)-esat(Tair))/0.1
      psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar=Patm/(Rconst*TairK)
      weighJ=1.0

!     boundary layer conductance for heat - single sided, forced convection
!     (Monteith 1973, P106 & notes dated 23/12/94)
      gbHu=0.003*sqrt(windUx/wleaf)    !m/s
!     raero=0.0                        !aerodynamic resistance s/m

      DO ileaf=1,2                  ! loop over sunlit and shaded leaves
!        first estimate of leaf temperature - assume air temp
         Tleaf(ileaf)=Tair
         Tlk=Tleaf(ileaf)+273.2    !Tleaf to deg K
!        first estimate of deficit at leaf surface - assume Da
         Dleaf=Dair                !Pa
!        first estimate for co2cs
         co2cs=co2ca               !mol/mol
         Qapar = (4.6e-6)*Qabs(1,ileaf)
!        ********************************************************************
         kr1=0                     !iteration counter for LE
         DO
!100        continue !    RETURN point for evaporation iteration
!           single-sided boundary layer conductance - free convection (see notes 23/12/94)
            Gras=1.595e8*abs(Tleaf(ileaf)-Tair)*(wleaf**3)     !Grashof
            gbHf=0.5*Dheat*(Gras**0.25)/wleaf
            gbH=gbHu+gbHf                         !m/s
            rbH=1./gbH                            !b/l resistance to heat transfer
            rbw=0.93*rbH                          !b/l resistance to water vapour
!           Y factor for leaf: stom_n = 1.0 for hypostomatous leaf;  stom_n = 2.0 for amphistomatous leaf
            rbH_L=rbH*stom_n/2.                   !final b/l resistance for heat  
            rrdn=1./grdn
            Y=1./(1.+ (rbH_L+raero)/rrdn)
!           boundary layer conductance for CO2 - single side only (mol/m2/s)
            gbc=Cmolar*gbH/1.32            !mol/m2/s
            gsc0=gsw0/1.57                        !convert conductance for H2O to that for CO2
            varQc=0.0                  
            weighR=1.0
!           respiration      
            Aleafx=-0.0089*Vcmxx*exp(0.069*(Tlk-293.2))
            gsc=gsc0
!           choose smaller of Ac, Aq
            Aleaf(ileaf) = Aleafx                     !mol CO2/m2/s
!           calculate new values for gsc, cs (Lohammer model)
            co2cs = co2ca-Aleaf(ileaf)/gbc
            co2Ci(ileaf) = co2cs-Aleaf(ileaf)/gsc
!           scale variables
            gsw=gsc*1.56                              !gsw in mol/m2/s
            gswv=gsw/Cmolar                           !gsw in m/s
            rswv=1./gswv
!           calculate evap'n using combination equation with current estimate of gsw
            Eleaf(ileaf)=                                       &
     &      (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/   &
     &      (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
!           calculate sensible heat flux
            Hleaf(ileaf)=Y*(Rnstar(ileaf)-Eleaf(ileaf))
!           calculate new leaf temperature (K)
            Tlk1=273.2+Tair+Hleaf(ileaf)*(rbH/2.+raero)/rhocp
!           calculate Dleaf use LE=(rhocp/psyc)*gsw*Ds
            Dleaf=psyc*Eleaf(ileaf)/(rhocp*gswv)
            gbleaf(ileaf)=gbc*1.32*1.075
            gsleaf(ileaf)=gsw

!          compare current and previous leaf temperatures
            IF(abs(Tlk1-Tlk) .GT. 3)Then            !Wan 2022/12/2
                Tlk=TairK
                exit
            ENDIF

            IF(abs(Tlk1-Tlk).LE.0.1)exit
            IF(kr1.GT.500)exit
!           update leaf temperature
            Tlk=Tlk1 
            Tleaf(ileaf)=Tlk1-273.2
            kr1=kr1+1
         ENDDO                          !solution not found yet
10    continue
      ENDDO
      RETURN
      END
!     ****************************************************************************
SUBROUTINE ciandA(Gma,Bta,g0,XX,Rd,co2Cs,gammas,ciquad,Aquad,itime)
!     calculate coefficients for quadratic equation for ci
      b2 = g0+XX*(Gma-Rd)
      b1 = (1.-co2cs*XX)*(Gma-Rd)+g0*(Bta-co2cs)-XX*(Gma*gammas+Bta*Rd)
      b0 = -(1.-co2cs*XX)*(Gma*gammas+Bta*Rd)-g0*Bta*co2cs
      bx=b1*b1-4.*b2*b0
      IF(bx.GT.0.0) THEN 
!       calculate larger root of quadratic
!       Weng & Lenuning, 1998, E7
        ciquad = (-b1+sqrt(bx))/(2.*b2)
      ENDIF
      IF(ciquad.LT.0.or.bx.LT.0.) THEN
        Aquad = 0.0
        ciquad = 0.7 * co2Cs
      ELSE
        Aquad = Gma*(ciquad-gammas)/(ciquad+Bta)
      ENDIF
      RETURN
END

!****************************************************************************
SUBROUTINE goud1(FLAIT,coszen,radabv,fbeam,               &
     &                  Tair,eairP,emair,emsoil,emleaf,sigma,   &
     &                  tauL,rhoL,rhoS,xfang,extkb,extkd,       &
     &                  reffbm,reffdf,extkbm,extkdm,Qcan)
!    wan- un-call subroutine
!    use the radiation scheme developed by
!    Goudriaan (1977, Goudriaan and van Larr 1995)
!=================================================================
!    Variable      unit      defintion
!    FLAIT         m2/m2     canopy leaf area index       
!    coszen                  cosine of the zenith angle of the sun
!    radabv(nW)    W/m2      incoming radiation above the canopy
!    fbeam                   beam fraction
!    fdiff                   diffuse fraction
!    funG(=0.5)              Ross's G function
!    extkb                   extinction coefficient for beam PAR
!    extkd                   extinction coefficient for diffuse PAR
!    albedo                  single scattering albedo
!    scatB                   upscattering parameter for beam
!    scatD                   upscattering parameter for diffuse
! ==================================================================
!    all intermediate variables in the calculation correspond
!    to the variables in the Appendix of of Seller (1985) with
!    a prefix of "x".
      INTEGER nW
      REAL radabv(3)
      REAL rhocbm(3),rhocdf(3)
      REAL reffbm(3),reffdf(3),extkbm(3),extkdm(3)
      REAL tauL(3),rhoL(3),rhoS(3),scatL(3)
      REAL Qcan(3,2)
!
!     for PAR: using Goudriann approximation to account for scattering
      fdiff=1.0-fbeam
      xu=coszen
      xphi1 = 0.5 -0.633*xfang - 0.33*xfang*xfang
      xphi2 = 0.877 * (1.0 - 2.0*xphi1)
      funG = xphi1 + xphi2*xu
      extkb=funG/xu
                       
!     Effective extinction coefficient for diffuse radiation Goudriaan & van Laar Eq 6.6)
      pi180=3.1416/180.
      cozen15=cos(pi180*15)
      cozen45=cos(pi180*45)
      cozen75=cos(pi180*75)
      xK15=xphi1/cozen15+xphi2
      xK45=xphi1/cozen45+xphi2
      xK75=xphi1/cozen75+xphi2
      transd=0.308*exp(-xK15*FLAIT)+0.514*exp(-xK45*FLAIT)+     &
     &       0.178*exp(-xK75*FLAIT)
      extkd=(-1./FLAIT)*alog(transd)

!     canopy reflection coefficients (Array indices: 1=VIS,  2=NIR
      DO nw=1,2                                                         !nw:1=VIS, 2=NIR
         scatL(nw)=tauL(nw)+rhoL(nw)                                    !scattering coeff
         IF((1.-scatL(nw))<0.0) scatL(nw)=0.9999                        !Weng 10/31/2008
         extkbm(nw)=extkb*sqrt(1.-scatL(nw))                            !modified k beam scattered (6.20)
         extkdm(nw)=extkd*sqrt(1.-scatL(nw))                            !modified k diffuse (6.20)
         rhoch=(1.-sqrt(1.-scatL(nw)))/(1.+sqrt(1.-scatL(nw)))          !canopy reflection black horizontal leaves (6.19)
         rhoc15=2.*xK15*rhoch/(xK15+extkd)                              !canopy reflection (6.21) diffuse
         rhoc45=2.*xK45*rhoch/(xK45+extkd)
         rhoc75=2.*xK75*rhoch/(xK75+extkd)   
       
         rhocbm(nw)=2.*extkb/(extkb+extkd)*rhoch                        !canopy reflection (6.21) beam 
         rhocdf(nw)=0.308*rhoc15+0.514*rhoc45+0.178*rhoc75
         reffbm(nw)=rhocbm(nw)+(rhoS(nw)-rhocbm(nw))        &               !effective canopy-soil reflection coeff - beam (6.27)
     &             *exp(-2.*extkbm(nw)*FLAIT)                              
         reffdf(nw)=rhocdf(nw)+(rhoS(nw)-rhocdf(nw))        &            !effective canopy-soil reflection coeff - diffuse (6.27)
     &             *exp(-2.*extkdm(nw)*FLAIT)  

!        by the shaded leaves
         abshdn=fdiff*(1.0-reffdf(nw))*extkdm(nw)                       &           !absorbed NIR by shaded
     &      *(funE(extkdm(nw),FLAIT)-funE((extkb+extkdm(nw)),FLAIT))    &
     &      +fbeam*(1.0-reffbm(nw))*extkbm(nw)                          &
!    &      *(funE(extkbm(nw),FLAIT)-funE((extkb+extkdm(nw)),FLAIT))    ! error found by De Pury
     &      *(funE(extkbm(nw),FLAIT)-funE((extkb+extkbm(nw)),FLAIT))    &
     &      -fbeam*(1.0-scatL(nw))*extkb                                &
     &      *(funE(extkb,FLAIT)-funE(2.0*extkb,FLAIT))
!        by the sunlit leaves
         absltn=fdiff*(1.0-reffdf(nw))*extkdm(nw)                       &  !absorbed NIR by sunlit
     &      *funE((extkb+extkdm(nw)),FLAIT)                             &
     &      +fbeam*(1.0-reffbm(nw))*extkbm(nw)                          &
!    &      *funE((extkb+extkdm(nw)),FLAIT)                         ! error found by De Pury
     &      *funE((extkb+extkbm(nw)),FLAIT)                             &
     &      +fbeam*(1.0-scatL(nw))*extkb                                &
     &      *(funE(extkb,FLAIT)-funE(2.0*extkb,FLAIT))

!        scale to REAL flux 
!        sunlit    
          Qcan(nw,1)=absltn*radabv(nw)
!        shaded
          Qcan(nw,2)=abshdn*radabv(nw)
      ENDDO
!     
!    calculate the absorbed (iso)thermal radiation
      TairK=Tair+273.2
      
!     apparent atmospheric emissivity for clear skies (Brutsaert, 1975)
      emsky=0.642*(eairP/Tairk)**(1./7)      !note eair in Pa

!     apparent emissivity from clouds (Kimball et al 1982)
      ep8z=0.24+2.98e-12*eairP*eairP*exp(3000.0/TairK)
      tau8=amin1(1.0,1-ep8z*(1.4-0.4*ep8z))                !ensure tau8<1
      emcloud=0.36*tau8*(1.-fbeam)*(1-10./TairK)**4        !10 from Tcloud = Tair-10 

!     apparent emissivity from sky plus clouds      
!     emair=emsky+emcloud
! 20/06/96
      emair=emsky
      IF(emair.GT.1.0) emair=1.0                             

      Bn0=sigma*(TairK**4)
      QLW1=-extkd*emleaf*(1.0-emair)*funE((extkd+extkb),FLAIT)      &
     &     -extkd*(1.0-emsoil)*(emleaf-emair)*exp(-2.0*extkd*FLAIT) &
     &     *funE((extkb-extkd),FLAIT)
      QLW2=-extkd*emleaf*(1.0-emair)*funE(extkd,FLAIT)              &
     &     -extkd*(1.0-emsoil)*(emleaf-emair)                       &
     &     *(exp(-extkd*FLAIT)-exp(-2.0*extkd*FLAIT))/extkd         &
     &     -QLW1
      Qcan(3,1)=QLW1*Bn0
      Qcan(3,2)=QLW2*Bn0
      RETURN
END

!****************************************************************************
SUBROUTINE photosyn(co2ca,CO2Csx,Dleafx,Tlkx,Qaparx,Gbcx, &
     &         a1,Ds0,fwsoil,varQc,weighR,                &
     &         g0,alpha,Vcmx1,eJmx1,weighJ,               &
     &         Aleafx,gscx,gddonset,itime,        &
     &         XX_rec,XX_med_rec,gscx_rec,gscx_med_rec,Dair,EqStomata)

    USE vars_consts
    IMPLICIT NONE

    REAL co2ca,CO2Csx,Dleafx,Tlkx,Qaparx,Gbcx
    REAL a1,Ds0,fwsoil,varQc,weighR
    REAL g0,alpha,Vcmx1,eJmx1,weighJ
    REAL Aleafx,gscx,gddonset
    REAL Acx,Aqx,Gma,Bta,XX,Rd,gammas,co2ci2,co2ci4
    REAL VcmxT,conKcT,conKoT,eJ,eJmxT,gamma,sps
    REAL Tdiff,Tlf,TmaxV,TminJ,TminV,TmaxJ,ToptV,ToptJ
    INTEGER itime,EqStomata
    REAL VJtemp,EnzK,fJQres ! Declare function
    REAL XX_rec,XX_med_rec,gscx_rec,gscx_med_rec,Dair,g1_med
    REAL XX_med,gscx_med
     
!     calculate Vcmax, Jmax at leaf temp (Eq 9, Harley et al 1992)
!     turned on by Weng, 2012-03-13
!     VcmxT = Vjmax(Tlkx,Trefk,Vcmx1,Eavm,Edvm,Rconst,Entrpy)
!     eJmxT = Vjmax(Tlkx,Trefk,eJmx1,Eajm,Edjm,Rconst,Entrpy)
      CO2Csx=AMAX1(CO2Csx,0.6*co2ca)


      IF(Qaparx.LE.0.) THEN
        !WRITE(*,*)'This is night'                    !night, umol quanta/m2/s
        !stop!改过
        Aleafx=-0.0089*Vcmx1*exp(0.069*(Tlkx-293.2))   ! original: 0.0089 Weng 3/22/2006
        gscx=g0
      ENDIF
!     calculate  Vcmax, Jmax at leaf temp using Reed et al (1976) function J appl Ecol 13:925
      TminV=-5.
      !TminV=gddonset/10.  ! original -5.        !-Jiang Jiang 2015/10/13

      TmaxV=50.
      ToptV=35.
      
      TminJ=TminV
      TmaxJ=TmaxV
      ToptJ=ToptV 
      
      Tlf=Tlkx-273.2

      VcmxT=VJtemp(Tlf,TminV,TmaxV,ToptV,Vcmx1)

      eJmxT=VJtemp(Tlf,TminJ,TmaxJ,ToptJ,eJmx1)      
!     calculate J, the asymptote for RuBP regeneration rate at given Q
      eJ = weighJ*fJQres(eJmxT,alpha,Qaparx,theta)
!     calculate Kc, Ko, Rd gamma*  & gamma at leaf temp
      conKcT = EnzK(Tlkx,Trefk,conKc0,Rconst,Ekc)
      conKoT = EnzK(Tlkx,Trefk,conKo0,Rconst,Eko)
!     following de Pury 1994, eq 7, make light respiration a fixed proportion of
!     Vcmax
      Rd = 0.0089*VcmxT*weighR                              !de Pury 1994, Eq7 !wan- de Pury, 1997
      Tdiff=Tlkx-Trefk
      gammas = gam0*(1.+gam1*Tdiff+gam2*Tdiff*Tdiff)       !gamma*
!     gamma = (gammas+conKcT*(1.+O2ci/conKoT)*Rd/VcmxT)/(1.-Rd/VcmxT)
      gamma = 25*1.0E-6!0.0   !wan- Walter Hill, in Algal Ecology, 1996 & 
      ! & Park S. Nobel, in Physicochemical and Environmental Plant Physiology (Fifth Edition), 2020
!     ***********************************************************************
!     Analytical solution for ci. This is the ci which satisfies supply and demand
!     functions simultaneously
!     calculate XX using Lohammer model, and scale for soil moisture
      a1= 11!1./(1.-0.7) 
      g1_med = 3 !3 !unit kpa^0.5 = 17.3 kpa
      Ds0 = 1000
      XX = a1*fwsoil/((co2csx - gamma)*(1.0 + Dleafx/Ds0))   !for stomatal conductance

      !MED model
      XX_med = 1.6*(1+g1_med/(SQRT(Dleafx/1000)))/co2csx
      IF(Dair .EQ. (50*0.001)) THEN
        XX_med = 0.
      ENDIF

      XX_rec = XX
      XX_med_rec = XX_med
      IF (EqStomata .EQ. 2) XX = XX_med

      !10 to 40 μmol m−2 s−1
      ! XX -- bulk stomatal and residual conductance for water vapour
      ! gamma -- CO2 compensation point 
      ! Dleafx -- water vapour mol fraction deficits (VPD) at the laef surface

!     calculate solution for ci when Rubisco activity limits A
      Gma = VcmxT  
      Bta = conKcT*(1.0+ o2ci/conKoT)
      CALL ciandA(Gma,Bta,g0,XX,Rd,co2Csx,gammas,co2ci2,Acx,itime)
!     calculate +ve root for ci when RuBP regeneration limits A
      Gma = eJ/4.
      Bta = 2.*gammas
!    calculate coefficients for quadratic equation for ci
      CALL ciandA(Gma,Bta,g0,XX,Rd,co2Csx,gammas,co2ci4,Aqx,itime)

!     choose smaller of Ac, Aq
      sps=AMAX1(0.001,sps)                  !Weng, 3/30/2006
      Aleafx = (amin1(Acx,Aqx) - Rd) !*sps     ! Weng 4/4/2006
!      IF(Aleafx.LT.0.0) Aleafx=0.0    ! by Weng 3/21/2006
!    calculate new values for gsc, cs (Lohammer model)
      CO2csx = co2ca-Aleafx/Gbcx
      gscx=g0 + XX*Aleafx  ! revised by Weng, wan- Weng & Luo 2008, eq.A5
      gscx_med = g0 + XX_med*Aleafx
      gscx_rec = gscx
      gscx_med_rec = gscx_med
      !IF (EqStomata .EQ. 2) gscx = gscx_med
      !gscx = gscx_med when Eqstomata eq 2 since I have replace the XX by XX_med
      RETURN
      END
!***********************************************************************
REAL function funeJ(alpha,eJmxT,Qaparx)
      funeJ=alpha*Qaparx*eJmxT/(alpha*Qaparx+2.1*eJmxT)
      RETURN
END
!****************************************************************************
REAL function esat(T)
!     returns saturation vapour pressure in Pa
      esat=610.78*exp(17.27*T/(T+237.3))     
      RETURN
END

!****************************************************************************
REAL function evapor(Td,Tw,Patm)
!* returns vapour pressure in Pa from wet & dry bulb temperatures
      gamma = (64.6 + 0.0625*Td)/1.e5
      evapor = esat(Tw)- gamma*(Td-Tw)*Patm
      RETURN
END

!****************************************************************************
REAL function Vjmax(Tk,Trefk,Vjmax0,Eactiv,Edeact,Rconst,Entrop)
      anum = Vjmax0*EXP((Eactiv/(Rconst*Trefk))*(1.-Trefk/Tk))
      aden = 1. + EXP((Entrop*Tk-Edeact)/(Rconst*Tk))
      Vjmax = anum/aden
      RETURN
      END
!****************************************************************************
REAL function funE(extkbd,FLAIT)
      funE=(1.0-exp(-extkbd*FLAIT))/extkbd
      RETURN
END

!     ****************************************************************************
!     Reed et al (1976, J appl Ecol 13:925) equation for temperature response
!     used for Vcmax and Jmax
REAL function VJtemp(Tlf,TminVJ,TmaxVJ,ToptVJ,VJmax0)
!constrain leaf temperatures between min and max
      !IF(Tlf.LT.TminVJ) Tlf=TminVJ
      IF(Tlf.GT.TmaxVJ) Tlf=TmaxVJ
      pwr=(TmaxVJ-ToptVJ)/(ToptVj-TminVj)
      VJtemp=VJmax0*((Tlf-TminVJ)/(ToptVJ-TminVJ))*     &
     &       ((TmaxVJ-Tlf)/(TmaxVJ-ToptVJ))**pwr 
      RETURN
END

!     ****************************************************************************
REAL function fJQres(eJmx,alpha,Q,theta)
    AX = theta                                 !a term in J fn
    BX = alpha*Q+eJmx                          !b term in J fn
    CX = alpha*Q*eJmx                          !c term in J fn
    IF((BX*BX-4.*AX*CX)>=0.0)THEN
        fJQres = (BX-SQRT(BX*BX-4.*AX*CX))/(2*AX)
    ELSE
        fJQres = (BX)/(2*AX)                   !Weng 10/31/2008
    ENDIF

    RETURN
END

!     *************************************************************************
REAL function EnzK(Tk,Trefk,EnzK0,Rconst,Eactiv)

    temp1=(Eactiv/(Rconst* Trefk))*(1.-Trefk/Tk)
!   IF (temp1<50.)THEN
    EnzK = EnzK0*EXP((Eactiv/(Rconst* Trefk))*(1.-Trefk/Tk))
!   ELSE
!   EnzK = EnzK0*EXP(50.)                                          ! Weng 10/31/2008
!   ENDIF

    RETURN
END

!     *************************************************************************
REAL function sinbet(doy,lat,pi,timeh)
    REAL lat
!   sin(bet), bet = elevation angle of sun
!   calculations according to Goudriaan & van Laar 1994 P30
    rad = pi/180.
!   sine and cosine of latitude
    sinlat = sin(rad*lat)
      coslat = cos(rad*lat)
!     sine of maximum declination
      sindec=-sin(23.45*rad)*cos(2.0*pi*(doy+10.0)/365.0)
      cosdec=sqrt(1.-sindec*sindec)
!     terms A & B in Eq 3.3
      A = sinlat*sindec
      B = coslat*cosdec
      sinbet = A+B*cos(pi*(timeh-12.)/12.)
      RETURN
END

!     *************************************************************************
SUBROUTINE yrday(doy,hour,lat,radsol,fbeam)
      REAL lat
      pi=3.14159256
      pidiv=pi/180.0
      slatx=lat*pidiv
      sindec=-sin(23.4*pidiv)*cos(2.0*pi*(doy+10.0)/365.0)             !wan- sin_delta Garratt, pp30, eqn. 3.3
      cosdec=sqrt(1.-sindec*sindec)
      a=sin(slatx)*sindec
      b=cos(slatx)*cosdec
      sinbet=a+b*cos(2*pi*(hour-12.)/24.)                              !wan- sine of the height of the sun
      solext=1370.0*(1.0+0.033*cos(2.0*pi*(doy-10.)/365.0))*sinbet     !wan- incident global radiation,1375 solar constant
      
      tmprat=radsol/solext                                             !wan- tau_a???

      tmpR=0.847-1.61*sinbet+1.04*sinbet*sinbet
      tmpK=(1.47-tmpR)/1.66

      IF(tmprat.LE.0.22) fdiff=1.0                                     
      !IF(tmprat .LE. 0.3) fdiff=1.0                                   !wan- Garratt, pp30, Figure3.2
      IF(tmprat.GT.0.22.and.tmprat.LE.0.35) THEN
        fdiff=1.0-6.4*(tmprat-0.22)*(tmprat-0.22)
      ENDIF
      !IF(tmprat.GT.0.3.and.tmprat.LE.0.7) THEN                        !wan- Garratt, pp30, Figure3.2     
      !    fdiff=1.0-6.4*(tmprat-0.3)*(tmprat-0.3)
      !ENDIF
      IF(tmprat.GT.0.35.and.tmprat.LE.tmpK) THEN
        fdiff=1.47-1.66*tmprat
      ENDIF
      IF(tmprat.ge.tmpK) THEN
        fdiff=tmpR
      ENDIF
      !IF(tmprat.ge.0.7) THEN                                          !wan- Garratt, pp30, Figure3.2
      !  fdiff=0.2
      !ENDIF

      fbeam=1.0-fdiff
      IF(fbeam.LT.0.0) fbeam=0.0
      RETURN
END

SUBROUTINE getCov(gamma,covfile,npara)    
    IMPLICIT NONE
    INTEGER npara,i,k
    REAL gamma(npara,npara)    
    CHARACTER(len=80) covfile
    
    OPEN(14,file=covfile,status='old')

    DO i=1,npara
        READ (14,*)(gamma(i,k),k=1,npara)
    ENDDO  
    RETURN
    END        
!MCMC
    SUBROUTINE costFObsNee(J_last, upgraded, isimu)

        USE FileSize
        IMPLICIT NONE

        !==========================================================
        ! 输入和输出参数
        !==========================================================

        REAL,    INTENT(INOUT) :: J_last
        INTEGER, INTENT(INOUT) :: upgraded
        INTEGER, INTENT(IN)    :: isimu

        !==========================================================
        ! 局部变量
        !==========================================================

        INTEGER :: i
        INTEGER :: j
        INTEGER :: nvalid
        INTEGER :: day_idx
        INTEGER :: vcode

        INTEGER :: litter_base_day

        REAL :: J_new
        REAL :: delta_J
        REAL :: r_num

        REAL :: J_leaf
        REAL :: J_stem
        REAL :: J_root
        REAL :: J_litter_cum
        REAL :: J_mbc
        REAL :: J_plab
        REAL :: J_psorb
        REAL :: J_pss
        REAL :: J_pocc
        REAL :: J_other

        REAL :: obs
        REAL :: obs_used
        REAL :: sd
        REAL :: wt
        REAL :: sim
        REAL :: d

        REAL :: litter_base_obs
        REAL :: cost_i

        LOGICAL :: litter_base_found

        !==========================================================
        ! 初始化
        !==========================================================

        J_new  = 0.0
        nvalid = 0
        J_leaf         = 0.0
        J_stem         = 0.0
        J_root         = 0.0
        J_litter_cum   = 0.0
        J_mbc          = 0.0
        J_plab         = 0.0
        J_psorb        = 0.0
        J_pss          = 0.0
        J_pocc         = 0.0
        J_other        = 0.0

        litter_base_day   = -999
        litter_base_obs   = 0.0
        litter_base_found = .FALSE.

        !==========================================================
        ! 找到累计凋落量 vcode=101 的第一条观测
        !
        ! 后续累计凋落量将转换为相对于该观测日的增量：
        !
        ! obs_used = 当前累计观测 - 第一条累计观测
        ! sim      = 第一条观测日之后的 OutC(1) 累计量
        !==========================================================

        DO j = 1, nobs_assim

            IF (INT(obs_TianTong(2,j)) .EQ. 101) THEN

                IF (.NOT. litter_base_found) THEN

                    litter_base_day = INT(obs_TianTong(1,j))
                    litter_base_obs = obs_TianTong(3,j)
                    litter_base_found = .TRUE.

                ELSE IF (INT(obs_TianTong(1,j)) .LT. litter_base_day) THEN

                    litter_base_day = INT(obs_TianTong(1,j))
                    litter_base_obs = obs_TianTong(3,j)

                ENDIF

            ENDIF

        ENDDO

        !==========================================================
        ! 遍历所有观测值并计算目标函数
        !==========================================================

        DO i = 1, nobs_assim

            !------------------------------------------------------
            ! 读取第 i 条观测
            !------------------------------------------------------

            day_idx = INT(obs_TianTong(1,i))
            vcode   = INT(obs_TianTong(2,i))

            obs = obs_TianTong(3,i)
            sd  = obs_TianTong(4,i)
            wt  = obs_TianTong(5,i)

            ! 默认情况下，实际进入代价函数的观测值就是原观测
            obs_used = obs

            !------------------------------------------------------
            ! 检查观测数据是否有效
            !------------------------------------------------------

            ! 标准差小于或等于0，无法标准化残差
            IF (sd .LE. 0.0) CYCLE

            ! 权重没有设置或小于等于0时，默认使用1
            IF (wt .LE. 0.0) wt = 1.0

            ! 观测日超出模型输出范围则跳过
            IF (day_idx .LT. 1 .OR. day_idx .GT. output_ndays) CYCLE

            !------------------------------------------------------
            ! 根据vcode提取相应的模型模拟值
            !------------------------------------------------------

            SELECT CASE (vcode)

            !======================================================
            ! 植物碳库
            !======================================================

            CASE (1)

                ! 叶片碳库
                sim = outputd_ccycle_Cpools(1, day_idx)

            CASE (2)

                ! 木质部碳库
                sim = outputd_ccycle_Cpools(2, day_idx)

            CASE (3)

                ! 细根碳库
                sim = outputd_ccycle_Cpools(3, day_idx)

            !======================================================
            ! 凋落物现存碳库
            !======================================================

            CASE (101)

                !------------------------------------------------------
                ! 地上总累计凋落物：
                ! OutC1 = 叶片凋落
                ! OutC2 = 木质部/枝条凋落
                ! OutC4 = 繁殖体凋落
                ! 不包括 OutC3，因为它是细根凋落
                !------------------------------------------------------

                IF (.NOT. litter_base_found) CYCLE

                IF (day_idx .LT. litter_base_day) CYCLE

                ! 实测值转换为相对于第一条观测的累计增量
                obs_used = obs - litter_base_obs

                ! 第一条累计观测只作为基线，不产生残差
                IF (day_idx .EQ. litter_base_day) CYCLE

                ! 模型地上总凋落物累计增量
                sim = SUM( &
                    outputd_outc(1, litter_base_day + 1 : day_idx) &
                    + outputd_outc(2, litter_base_day + 1 : day_idx) &
                    + outputd_outc(4, litter_base_day + 1 : day_idx) &
                )

            !======================================================
            ! 微生物生物量碳
            !======================================================

            CASE (107)

                ! 第7个碳库是微生物生物量碳
                sim = outputd_ccycle_Cpools(7, day_idx)

            !======================================================
            ! 无机磷库
            !======================================================

            CASE (201)

                ! 活性无机磷库 QPlab
                sim = outputd_Pdynamic(6, day_idx)

            CASE (202)

                ! 吸附态无机磷库 QPsorb
                sim = outputd_Pdynamic(7, day_idx)

            CASE (203)

                ! 强吸附态无机磷库 QPss
                sim = outputd_Pdynamic(8, day_idx)

            CASE (204)

                ! 闭蓄态无机磷库 QPocc
                sim = outputd_Pdynamic(9, day_idx)

            !======================================================
            ! 未定义的变量编码
            !======================================================

            CASE DEFAULT

                CYCLE

            END SELECT

            !------------------------------------------------------
            ! 检查模型值是否有效
            !------------------------------------------------------

            ! NaN判断：NaN是唯一不等于自身的浮点数
            IF (sim .NE. sim) CYCLE

            !------------------------------------------------------
            ! 输出前20条观测的检查信息
            !------------------------------------------------------

            IF (i .LE. 20) THEN

                WRITE(*,'(A,I5,A,I6,A,I5,A,F14.4,A,F14.4,A,F12.4,A,F8.3)') &
                    'DEBUG i=', i, &
                    ' day=', day_idx, &
                    ' vcode=', vcode, &
                    ' obs=', obs_used, &
                    ' sim=', sim, &
                    ' sd=', sd, &
                    ' wt=', wt

            ENDIF
            ! ----------------------------------------------
            ! 计算相对误差
            ! sim 和 obs_used 必须是相同定义下的量
            ! ----------------------------------------------

            ! IF (ABS(obs_used) .LT. 1.0E-8) CYCLE

            d = (sim - obs_used) / ABS(obs_used)

            ! cost_i = wt * 0.5 * d * d

            ! J_new = J_new + cost_i

            ! ----------------------------------------------
            
            !d = (sim - obs_used) / sd

            
            ! 当前这一条观测对总代价的贡献
            cost_i = wt * 0.5 * d * d

            ! 累加总代价
            J_new = J_new + cost_i
            
            IF (vcode .EQ. 201 .AND. &
    (isimu .LE. 3 .OR. MOD(isimu,10) .EQ. 0)) THEN

    WRITE(*,'(A,I7,A,I6,A,F14.6,A,F14.6,A,F12.6,A,F8.3,A,F14.5)') &
        'PLAB DETAIL isimu=', isimu, &
        ' day=', day_idx, &
        ' obs=', obs, &
        ' sim=', sim, &
        ' sd=', sd, &
        ' wt=', wt, &
        ' cost_i=', cost_i

ENDIF

            ! 根据vcode分别累计不同变量的代价
            SELECT CASE (vcode)

            CASE (1)
                J_leaf = J_leaf + cost_i

            CASE (2)
                J_stem = J_stem + cost_i

            CASE (3)
                J_root = J_root + cost_i

            CASE (101)
                J_litter_cum = J_litter_cum + cost_i

            CASE (107)
                J_mbc = J_mbc + cost_i

            CASE (201)
                J_plab = J_plab + cost_i

            CASE (202)
                J_psorb = J_psorb + cost_i

            CASE (203)
                J_pss = J_pss + cost_i

            CASE (204)
                J_pocc = J_pocc + cost_i

            CASE DEFAULT
                J_other = J_other + cost_i

            END SELECT

            nvalid = nvalid + 1
        ENDDO

        ! 前5次模拟每次都打印，之后每100次打印一次
        IF (isimu .LE. 5 .OR. MOD(isimu,10) .EQ. 0) THEN

            WRITE(*,'(A,I8,A,I8,A,ES16.7,A,ES16.7,A,I8)') &
                'COST CHECK isimu=', isimu, &
                ' nvalid=', nvalid, &
                ' J_new=', J_new, &
                ' J_last=', J_last, &
                ' upgraded=', upgraded

            WRITE(*,'(A,ES16.7)') '  J_leaf         = ', J_leaf
            WRITE(*,'(A,ES16.7)') '  J_stem         = ', J_stem
            WRITE(*,'(A,ES16.7)') '  J_root         = ', J_root
            WRITE(*,'(A,ES16.7)') '  J_litter_cum   = ', J_litter_cum
            WRITE(*,'(A,ES16.7)') '  J_mbc          = ', J_mbc
            WRITE(*,'(A,ES16.7)') '  J_plab         = ', J_plab
            WRITE(*,'(A,ES16.7)') '  J_psorb        = ', J_psorb
            WRITE(*,'(A,ES16.7)') '  J_pss          = ', J_pss
            WRITE(*,'(A,ES16.7)') '  J_pocc         = ', J_pocc
            WRITE(*,'(A,ES16.7)') '  J_other        = ', J_other

            WRITE(*,'(A,ES16.7)') &
                '  J_component_sum = ', &
                J_leaf + J_stem + J_root + &
                 J_litter_cum + J_mbc + &
                J_plab + J_psorb + J_pss + J_pocc + J_other

        ENDIF

        !==========================================================
        ! 此处继续保留你原有的MCMC接受/拒绝代码
        ! 例如使用J_new、J_last、delta_J和随机数r_num
        !==========================================================
                !==========================================================
        ! 检查目标函数是否有效
        !==========================================================

        IF (nvalid .LE. 0) THEN
            WRITE(*,*) 'REJECT: no valid observations'
            RETURN
        ENDIF

        ! NaN检查
        IF (J_new .NE. J_new) THEN
            WRITE(*,*) 'REJECT: J_new is NaN, isimu=', isimu
            RETURN
        ENDIF

        ! 如果模型运行过程中出现异常，则拒绝该参数组
        IF (rejectnow .NE. 0) THEN
            WRITE(*,*) 'REJECT: model failure, rejectnow=', rejectnow
            RETURN
        ENDIF

        !==========================================================
        ! Metropolis-Hastings接受或拒绝
        !==========================================================

        delta_J = J_new - J_last

        CALL RANDOM_NUMBER(r_num)

        IF (delta_J .LE. 0.0) THEN

            ! 新参数使目标函数变小：直接接受
            J_last = J_new
            upgraded = upgraded + 1

            WRITE(*,'(A,I8,A,I8,A,ES16.7,A,ES16.7)') &
                'ACCEPT BETTER: isimu=', isimu, &
                ' upgraded=', upgraded, &
                ' J_last=', J_last, &
                ' delta_J=', delta_J

        ELSE

            ! 新参数目标函数较大：按概率接受
            !
            ! 当delta_J很大时，EXP(-delta_J)接近0，
            ! 为避免数值下溢，直接拒绝。

            IF (delta_J .LT. 80.0) THEN

                IF (r_num .LT. EXP(-delta_J)) THEN

                    J_last = J_new
                    upgraded = upgraded + 1

                    WRITE(*,'(A,I8,A,I8,A,ES16.7,A,F10.6,A,ES16.7)') &
                        'ACCEPT WORSE: isimu=', isimu, &
                        ' upgraded=', upgraded, &
                        ' delta_J=', delta_J, &
                        ' random=', r_num, &
                        ' probability=', EXP(-delta_J)

                ELSE

                    WRITE(*,'(A,I8,A,ES16.7,A,F10.6,A,ES16.7)') &
                        'REJECT: isimu=', isimu, &
                        ' delta_J=', delta_J, &
                        ' random=', r_num, &
                        ' probability=', EXP(-delta_J)

                ENDIF

            ELSE

                WRITE(*,'(A,I8,A,ES16.7)') &
                    'REJECT LARGE DELTA: isimu=', isimu, &
                    ' delta_J=', delta_J

            ENDIF

        ENDIF

        RETURN

    END SUBROUTINE costFObsNee


!     SUBROUTINE costFObsNee (J_last,upgraded,isimu)
!     !(output_daily_mcmc,obs_TianTong,len1,J_last,upgraded) len1=nobs
!     USE FileSize
!     IMPLICIT NONE
!     INTEGER day,hour
!     REAL J_new,J_last,delta_J
!     REAL J_gpp,J_nee,J_er,J_foliage,J_fnpp
!     REAL J_sw
!     REAL J_wood,J_wnpp,J_root,J_rnpp,J_soilc,J_pheno
!     REAL tmp1, dObsSim, random_harvest
!     INTEGER i,len1,upgraded,m,isimu,count_999
!     INTEGER j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11
!     REAL r_num
! 	REAL J_lai,J_biomass
! 	REAL std_obs_gpp,std_obs_er,std_obs_nee,std_obs_lai,std_obs_biomass  !add by zhoujian
!     REAL, dimension(8760) :: GPP_mod,GPP_obs,NEE_mod,NEE_obs
!     REAL slope_GPP,slope_NEE

! !    compute J_obs
!     J_gpp=0.0
!     J_nee=0.0
!     J_er=0.0
!     j1=0
!     j2=0
!     j3=0

!     ! GPP _ interp
!     std_obs_gpp=0.348!0.292
! 	std_obs_er=0.147!0.149
! 	std_obs_nee=0.355!0.269

!     !len1 = nobs
!     count_999 = 0
!     len1 = 8760 !nobs
!     DO i=1,len1		!len1=nobs=365
!         !day=int(obs_TianTong(1,i))        !index for simudailyflux data 
!         hour = int(obs_TianTong(1,i))  
!         IF(obs_TianTong(2,i).GT.-999)THEN !2 for gpp  
!             j1=j1+1
!             !dObsSim=output_daily_mcmc(1,i)-obs_TianTong(2,i)  !output 1 --gpp, obs 2--gpp
!             !IF(obs_TianTong(2,i) .lt. 0) obs_TianTong(2,i) = 0 
!             dObsSim=output_mcmc(1,i)-obs_TianTong(2,i)  !output 1 --gpp, obs 2--gpp
!             J_gpp=J_gpp+(dObsSim*dObsSim)  
!             !IN SPRUCE: J_gpp=J_gpp+(dObsSim*dObsSim)/(2*std(2,i)*std(2,i))
!         ENDIF
! !        			
!         IF(obs_TianTong(3,i).GT.-999)THEN
!             j3=j3+1
!             !dObsSim=output_daily_mcmc(2,i)-obs_TianTong(3,i) !output 2 --RECO, obs 3--RECO
!             dObsSim=output_mcmc(2,i)-obs_TianTong(3,i) !output 2 --RECO, obs 3--RECO
!             J_er=J_er+(dObsSim*dObsSim)
!         ENDIF

!         IF(obs_TianTong(4,i).GT.-999)THEN
!             j2=j2+1
!             !dObsSim=output_daily_mcmc(3,i)-obs_TianTong(4,i) !output 3 --NEE, obs 4--NEE
!             dObsSim=output_mcmc(3,i)-obs_TianTong(4,i) !output 2 --RECO, obs 3--RECO
!             J_nee=J_nee+(dObsSim*dObsSim)
!         ENDIF
        
!     ENDDO

!     GPP_mod =  output_mcmc(1,:)
!     GPP_obs =  obs_TianTong(2,:)
!     !CALL linear_regression(GPP_obs,GPP_mod,slope_GPP)
!     !write(*,*)'sum_GPP', sum(output_mcmc(1,:))!,sum(obs_TianTong(2,:))

!     NEE_mod =  output_mcmc(3,:)
!     NEE_obs =  obs_TianTong(4,:)
!     !CALL linear_regression(NEE_obs,NEE_mod,slope_NEE)
!     !write(*,*)'sum_NEE', sum(output_mcmc(3,:)),sum(obs_TianTong(4,:))
    
! 	J_gpp=J_gpp/(2*std_obs_gpp*std_obs_gpp)
!     J_nee=J_nee/(2*std_obs_nee*std_obs_nee)
!     J_er=J_er/(2*std_obs_er*std_obs_er)
	
!     J_gpp=J_gpp/REAL(j1)
!     J_nee=J_nee/REAL(j2)
!     J_er=J_er/REAL(j3)

!     !J_gpp=ABS(1-slope_GPP)*J_gpp/REAL(j1)
!     !J_nee=ABS(1-slope_NEE)*J_nee/REAL(j2)
!     !J_er=J_er/REAL(j3)

!     J_new=J_er+J_gpp+J_nee
!     J_new=J_new*600 

!     IF(isimu .eq. 1) J_last = J_new  ! wan, avoid infinity at the first simulation
!     delta_J=J_new-J_last           !    delta_J=(J_new-J_last)/J_last
!     CALL random_number(r_num)

! !     delta_J=-1     !accept all samples, No data
!     IF(ISNAN(J_new))THEN
!         WRITE(*,*)'In cosNEE: NaN RETURN, upgraded', upgraded 
!         RETURN
!     ENDIF

!     IF (rejectnow .eq. 0) THEN
!         IF(AMIN1(1.0,exp(-delta_J)).gt.r_num)THEN 
!         ! IF(AMIN1(1.0,exp(-delta_J)).gt.0.95)THEN  ! wan change this accepet probability
!             upgraded=upgraded+1
!             J_last=J_new    
!         ENDIF
!     ENDIF

!     RETURN
!     END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Square root of a matrix							  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE racine_mat(M, Mrac,npara)

    INTEGER npara,i
    REAL M(npara,npara),Mrac(npara,npara)
    REAL valpr(npara),vectpr(npara,npara)
    Mrac=0.
    CALL jacobi(M,npara,npara,valpr,vectpr,nrot)
    DO i=1,npara
	IF(valpr(i).ge.0.) THEN
            Mrac(i,i)=sqrt(valpr(i))
	ELSE
            print*, 'WARNING!!! Square root of the matrix is undefined.'
            print*, ' A negative eigenvalue has been set to zero - results may be wrong'
            Mrac=M
            RETURN
	ENDIF
    ENDDO
    Mrac=matmul(matmul(vectpr, Mrac),transpose(vectpr))

END SUBROUTINE racine_mat      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Extraction of the eigenvalues and the eigenvectors !!
!! of a matrix (Numerical Recipes)					  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE jacobi(a,n,np,d,v,nrot)
INTEGER :: n,np,nrot
REAL :: a(np,np),d(np),v(np,np)
INTEGER, PARAMETER :: NMAX=500
INTEGER :: i,ip,iq,j
REAL :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      
DO ip=1,n
	DO iq=1,n
		v(ip,iq)=0.
	END DO
	v(ip,ip)=1.
END DO

DO ip=1,n
	b(ip)=a(ip,ip)
	d(ip)=b(ip)
	z(ip)=0.
END DO

nrot=0
DO i=1,50
	sm=0.
	DO ip=1,n-1
		DO iq=ip+1,n
			sm=sm+abs(a(ip,iq))
		END DO
	END DO
	IF(sm.eq.0.)RETURN
	IF(i.LT.4)THEN
		tresh=0.2*sm/n**2
	ELSE
		tresh=0.
	ENDIF
	DO ip=1,n-1
		DO iq=ip+1,n
			g=100.*abs(a(ip,iq))
			IF((i.GT.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))THEN
				a(ip,iq)=0.
			ELSE IF(abs(a(ip,iq)).GT.tresh)THEN
				h=d(iq)-d(ip)
				IF(abs(h)+g.eq.abs(h))THEN
					t=a(ip,iq)/h
				ELSE
					theta=0.5*h/a(ip,iq)
					t=1./(abs(theta)+sqrt(1.+theta**2))
					IF(theta.LT.0.) THEN
						t=-t
					ENDIF
				ENDIF
				c=1./sqrt(1+t**2)
				s=t*c
				tau=s/(1.+c)
				h=t*a(ip,iq)
				z(ip)=z(ip)-h
				z(iq)=z(iq)+h
				d(ip)=d(ip)-h
				d(iq)=d(iq)+h
				a(ip,iq)=0.
				DO j=1,ip-1
					g=a(j,ip)
					h=a(j,iq)
					a(j,ip)=g-s*(h+g*tau)
					a(j,iq)=h+s*(g-h*tau)
				END DO
				DO j=ip+1,iq-1
					g=a(ip,j)
					h=a(j,iq)
					a(ip,j)=g-s*(h+g*tau)
					a(j,iq)=h+s*(g-h*tau)
				END DO
				DO j=iq+1,n
					g=a(ip,j)
					h=a(iq,j)
					a(ip,j)=g-s*(h+g*tau)
					a(iq,j)=h+s*(g-h*tau)
				END DO
				DO j=1,n
					g=v(j,ip)
					h=v(j,iq)
					v(j,ip)=g-s*(h+g*tau)
					v(j,iq)=h+s*(g-h*tau)
				END DO
				nrot=nrot+1
			ENDIF
		END DO
	END DO
	DO ip=1,n
		b(ip)=b(ip)+z(ip)
		d(ip)=b(ip)
		z(ip)=0.
	END DO
END DO
print*, 'too many iterations in jacobi' 
RETURN
END SUBROUTINE jacobi


!===================================================
!       generate new coefficents
        SUBROUTINE coefgenerate(isimu,coefac,coefmax,coefmin,coef,search_length,npara)
        
        INTEGER npara
        REAL coefac(npara),coefmax(npara),coefmin(npara),coef(npara)
        REAL r,coefmid,random_harvest
        INTEGER i
        REAL search_length
        DO i=1,npara
999         continue
!w!rite(*,*)'continue',i
            CALL random_number(random_harvest)
            r=random_harvest-0.5
            coef(i)=coefac(i)+r*(coefmax(i)-coefmin(i))*search_length
            IF(coef(i).GT.coefmax(i).or.coef(i).LT.coefmin(i))goto 999
        ENDDO
        !WRITE(*,*)isimu,(coef(i),i=1,npara)
!! ********** added to debug memo error in MCMC *************  11/13/2017           
!            WRITE(91,901)isimu,(coef(i),i=1,npara)
!901         format(I12,",",8(F15.4,","))        
!! ********** added to debug memo error in MCMC *************  11/13/2017          
        
        RETURN
        END
!============
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Generation of a random vector from a multivariate  !!
!! normal distribution with mean zero and covariance  !!
!! matrix gamma.									  !!
!! Beware!!! In order to improve the speed of the	  !!
!! algorithms, the SUBROUTINE use the Square root	  !!
!! matrix of gamma									  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE gengaussvect(gamma_racine,xold,xnew,npara)

INTEGER npara
REAL gamma_racine(npara,npara)
REAL x(npara),xold(npara),xnew(npara)

DO i=1,npara
    x(i)=rangauss(25)
ENDDO

x = matmul(gamma_racine, x)
xnew = xold + x
END SUBROUTINE gengaussvect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Generation of a random number from a standard	  !!
!! normal distribution. (Numerical Recipes)           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function rangauss(idum)


INTEGER idum
REAL v1, v2, r, fac, gset
REAL r_num

data iset/0/
IF(iset==0) THEN
1	CALL random_number(r_num)
        v1=2.*r_num-1
        CALL random_number(r_num)
	v2=2.*r_num-1
	r=(v1)**2+(v2)**2
	IF(r>=1) go to 1
	fac=sqrt(-2.*log(r)/r)
	gset=v1*fac
	rangauss=v2*fac
	iset=1
ELSE
	rangauss=gset
	iset=0
END IF

RETURN
END function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! variance matrix of a matrix of data				  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE varcov(tab,varcovar,npara,ncov)

INTEGER npara,ncov
REAL tab(ncov,npara),tab2(ncov,npara)
REAL varcovar(npara,npara)

CALL centre(tab,tab2,npara,ncov)

varcovar = matmul(transpose(tab2), tab2)*(1./REAL(ncov))

END SUBROUTINE varcov

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute the centered matrix, ie. the matrix minus  !!
!! the column means									  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE centre(mat,mat_out,npara,ncov)

    INTEGER npara,ncov
    REAL mat(ncov,npara),mat_out(ncov,npara)
    REAL mean

DO i=1,npara
    mat_out(:,i) = mat(:,i) - mean(mat(:,i),ncov)
ENDDO

END SUBROUTINE centre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mean of a vector									  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Function mean(tab,ncov)
    INTEGER ncov
REAL tab(ncov)
REAL mean,mean_tt
mean_tt=0.
DO i=1,ncov	
mean_tt=mean_tt+tab(i)/REAL(ncov)
ENDDO
mean=mean_tt
END Function

SUBROUTINE linear_regression (x1,y1,slope)
    implicit none
    real :: x1(8760), y1(8760)
    real :: slope, intercept, sum_x, sum_y, sum_xy, sum_x_squared
    real :: n, denominator
    integer :: i, valid_count

    ! calculated the sum, skip the NA value (999)
    sum_x = 0.0
    sum_y = 0.0
    sum_xy = 0.0
    sum_x_squared = 0.0
    valid_count = 0
    do i = 1, 8760
        if (x1(i) /= -999) then ! x is observation, have the NA value
            sum_x = sum_x + x1(i)
            sum_y = sum_y + y1(i)
            sum_xy = sum_xy + x1(i) * y1(i)
            sum_x_squared = sum_x_squared + x1(i) * x1(i)
            valid_count = valid_count + 1
        end if
    end do

    ! calculate the slope and intercept
    n = real(valid_count) ! valid data point
    denominator = n * sum_x_squared - sum_x * sum_x
    slope = (n * sum_xy - sum_x * sum_y) / denominator
    intercept = (sum_y * sum_x_squared - sum_x * sum_xy) / denominator
    print *, "slop: ", slope,valid_count
    print *, "interp: ", intercept
END  SUBROUTINE linear_regression

real function sum_without_na(arr)
    real, intent(in) :: arr(:)
    integer :: i
    sum_without_na = 0.0
    do i = 1, size(arr)
        if (arr(i) /= NA) then
            sum_without_na = sum_without_na + arr(i)
        end if
    end do
end function sum_without_na
