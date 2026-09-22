MODULE IntersVariables
IMPLICIT NONE
SAVE
TYPE :: ResultType2
INTEGER :: nleap_rere
INTEGER :: force_nhours_rere
INTEGER :: output_ndays_rere
INTEGER :: nyear_rere
END TYPE ResultType2

INTEGER,PARAMETER :: SGL = SELECTED_REAL_KIND(p=6)  !SINGLE
INTEGER,PARAMETER :: DBL = SELECTED_REAL_KIND(p=13) !DOUBLE
INTEGER,PARAMETER :: intn = SELECTED_INT_KIND(9) 
INTEGER,PARAMETER :: pools = 9
REAL,PARAMETER :: fbmC = 0.55
INTEGER start_year,end_year
INTEGER :: nobs_assim

REAL (SGL) :: Fptase_slow,Fptase_pass,N_loss_mic
REAL (SGL):: GPP,NEE,Reco,NPP,xnpNPP,NPP_noadd,NPP_FORNSC,NEE_new
REAL (SGL):: Acan1,Acan2,fslt,fshd
REAL (SGL):: Rgrowth,Rnitrogen,Rphosphorus,Rmain,Rauto,Rhetero
REAL (SGL):: GrowthP,xnpGrowthP,bmleaf,bmstem,bmroot,bmRe,bmplant,xnpGrowthP_temp
REAL (SGL):: NSC,fnsc,storage,add,accumulation,store
REAL (SGL):: LAI,ht,L_tau_new
! allocation ratio to Leaf, stem, and Root
REAL (SGL):: alpha_L,alpha_W,alpha_R,alpha_Re,a_L,a_S,a_R,a_Re,NPP_L,NPP_W,NPP_R,NPP_Re
REAL (SGL):: OutC(pools),OutN(pools),OutP(pools)
REAL (SGL):: QC(pools),QN(pools),QP(pools),CN(pools),CP(pools),CN0(pools),CP0(pools)
REAL (SGL):: N_miner,N_immob,N_leach,N_vol,N_net,N_uptake,N_fixation,N_transfer,N_deficit
REAL (SGL):: P_transfer,P_deficit,P_deficit_test,P_demand_test
REAL (SGL):: xde,xnp_leaf_limit,xnp_stem_limit,xnp_root_limit,xNPuptake,LfactorN,LfactorP
! This for new xde
REAL (SGL):: xde2,OutC2(pools),OutN2(pools),OutP2(pools),N_miner2,N_immob2,N_net2
REAL (SGL):: P_miner2,P_immob2,P_net2
REAL (SGL):: P_miner22,P_immob22,P_net22,N_miner22,N_immob22,N_net22,xde22
!-
REAL (SGL):: P_leaf,P_wood,P_root,P_re,N_leaf,N_wood,N_root,N_re,NSP,NSN
REAL (SGL):: N_demand,P_demand
REAL (SGL):: Dsoillab,P_net,Fptase,P_loss,P_uptake,P_miner,P_immob,asorb,bss
REAL (SGL):: QPlab,QPsorb,QPss,QNminer,QPocc
REAL (SGL):: S_omega,W,Sw,wsc(10),wcl(10),omega,zwt,S_t(5),st !wcl, volum ratio
REAL (SGL):: gamma_T,gamma_W
REAL (SGL):: GDD5
REAL (SGL):: tauC(pools),exitK(pools)
REAL (SGL):: x_stem_limit,x_root_limit,x_leaf_limit,xuptake
REAL (SGL) :: etaL,etaW,etaR,etaRe,f_F2M,f_C2M,f_C2S,f_M2S,f_M2P,f_S2M,f_S2P,f_P2M
REAL (SGL) :: etaL2,etaW2,etaR2,etaRe2
REAL (SGL) :: s_w_min,wsmin,wsmax
REAL tau_L,tau_W,tau_R,tau_Re
REAL tau_F,tau_C,tau_Micr,tau_Slow,tau_Pass
! canopy
REAL Acanop,Hcanop
REAL evap,transp        !wan- recalculater soil evaporation accord. to Weng & Luo 2008, eq.A9.
REAL (SGL):: Aleaf(2),Tlk1,Tlk2
REAL Rsoilabs!Total radiation absorbed by soil 
INTEGER (SGL) itime,idays,Flag_N,Flag_P,Flag_C,Flag_Growth,Flag_GPP,Flag_QP,Flag_NSC,rejectnow
REAL Tairmax
INTEGER, Dimension(1) :: Tairmax_Loc
INTEGER Tairmax_Loc_value
INTEGER onset !flag of phenological stage
REAL (SGL) :: aN,aP,bN,bP
INTEGER, DIMENSION(:), allocatable :: year_all 
INTEGER, DIMENSION(:), allocatable :: uni_Year
INTEGER num_cols
INTEGER force_nhours,output_ndays,nobs,ndays,nyear,lastyear_days,lastyear_hours !joutput = the number of output varible
INTEGER, DIMENSION(:), allocatable :: leap_year
TYPE (ResultType2) :: rere
REAL (SGL) :: P_reserve, N_reserve,NSPmax,NSPmin
REAL (SGL) :: rootSurfaceP,rate_maxN,rate_maxP
INTEGER, PARAMETER :: npara_all = 61
REAL (SGL) :: site_paras(npara_all),InitialCNP(pools,3) !change for paras
INTEGER iyear,resoforce_hr !INTEGER resoforce_hr =temporal resolution
INTEGER n_leap_year,nspinup,ispinup
INTEGER CYCLE_CNP,NutrientsUp,MCMC,NDSPINUP,SensTest,EQstomata,NPaddition
CHARACTER(len=99) fileplace,filesignal
REAL infilt,fwsoil,topfws
REAL f_rootp
REAl Traceresults(10,10)
REAL Rh_pools(5)            !fine lit.,coarse lit.,Micr,Slow,Pass
REAL ExcessC, EC_store, EC_out, exitEC, Rauto_new
REAL XX_rec,XX_med_rec,gscx_rec,gscx_med_rec,XX_med,g1_med,Gscx_med
REAL,DIMENSION(3)::newsoilNCmax,newsoilNCmin,newsoilNC
!REAL, DIMENSION(8760) :: 

REAL, DIMENSION(:), ALLOCATABLE :: MaxUpCapP_rec,bP_rec,aP_rec,fwsoil_rec,omega_rec
REAL, DIMENSION(:), ALLOCATABLE :: P_loss_rec,Fptase_rec,asorb_rec,Dsoillab_rec
REAL, DIMENSION(:), ALLOCATABLE :: P_miner2_rec,P_immob2_rec,P_net2_rec
REAL, DIMENSION(:), ALLOCATABLE :: Dsoillab_rec2, P_uptake_rec,P_in_out_rec,P_leach_rec
REAL, DIMENSION(:), ALLOCATABLE :: FdiffP_rec,delta_QLabP_root_rec,D_rec
REAL, DIMENSION(:), ALLOCATABLE :: swc_p_rec,tf_rec,rdiff_rec,tf1_rec,tf2_rec
REAL, DIMENSION(:), ALLOCATABLE :: f_rootp_rec,QC3_rec,LfactorP_rec
REAL, DIMENSION(:), ALLOCATABLE :: P_deficit_rec,P_demand_rec,P_deficit_rec_test,P_demand_rec_test
REAL, DIMENSION(:), ALLOCATABLE :: N_demand_rec,N_uptake_rec,N_deficit_rec
REAL, DIMENSION(:), ALLOCATABLE :: N_fixation_x_rec
REAL, DIMENSION(:,:), ALLOCATABLE :: Xss_new_step_rec
!!!!!!
! REAL(SGL) :: s_Psorbmax      = 2.94
! REAL(SGL) :: s_kmlabP        = 0.21
! REAL(SGL) :: s_sorb_to_ssorb = 0.84
! REAL(SGL) :: s_ssorb_to_occ  = 1.15
! REAL(SGL) :: s_P_loss        = 0.21

REAL(SGL) :: s_Psorbmax      = 1.00
REAL(SGL) :: s_kmlabP        = 1.00
REAL(SGL) :: s_sorb_to_ssorb = 1.00
REAL(SGL) :: s_ssorb_to_occ  = 1.00
<<<<<<< /home/user/fangwei/FangxiuWan-TECO-CNP-Sv1.0-edff198/TECO-CNP Sv1.0 20250316/main/FileSize.f90
REAL(SGL) :: s_P_loss        = 1.00
REAL :: P_add_rate_kg_ha_yr
||||||| /home/user/fangwei/TECO_sensitivity_server/baseline_main/FileSize.f90
REAL(SGL) :: s_P_loss        = 1.00
REAL :: P_add_rate_kg_ha_yr
=======
REAL(SGL) :: s_P_loss        = 1.00
REAL :: P_add_rate_kg_ha_yr

! ------------------------------------------------------------------
! External sensitivity-analysis controls.  Every multiplier defaults
! to one, so a normal TECO run is bit-for-bit unchanged unless the
! corresponding TECO_* environment variable is explicitly supplied.
! ------------------------------------------------------------------
INTEGER :: teco_sensitivity_mode = 0
REAL(SGL) :: sens_alloc_kn       = 1.0
REAL(SGL) :: sens_alloc_ww       = 1.0
REAL(SGL) :: sens_alloc_el       = 1.0
REAL(SGL) :: sens_alloc_es       = 1.0
REAL(SGL) :: sens_alloc_repro    = 1.0
REAL(SGL) :: sens_kpmin          = 1.0
REAL(SGL) :: sens_kpmin_low      = 1.0
REAL(SGL) :: sens_fptase_capacity= 1.0
>>>>>>> /home/user/fangwei/TECO_sensitivity_server/patched_main/FileSize.f90


END MODULE IntersVariables

!-----------------------------

MODULE DaysHours
    !Module to decalre Translate fraction data to share between routines
    USE IntersVariables
    IMPLICIT NONE

    TYPE :: ResultType
    INTEGER :: Days
    INTEGER :: Hours
    END TYPE ResultType

    CONTAINS

    FUNCTION check_leap(year)
        IMPLICIT NONE
        INTEGER year,check_leap,n

        if (mod(year,4) == 0 .and. mod(year,100)/=0) then ! 能被4整除且不能被100整除为闰年
                write(*,*) " leap year",year 
                check_leap = 1
        elseif (mod(year,400) == 0) then ! 能被400整除也为闰年
                write(*,*) " leap year",year 
                check_leap = 1
        else
                check_leap =0
        end if
    END FUNCTION check_leap

        FUNCTION YearDayHour() RESULT(rere1)

        IMPLICIT NONE

        INTEGER :: n, nleap1, year, check_leap1
        INTEGER :: force_nhours1, output_ndays1, nyear1
        TYPE(ResultType2) :: rere1

        output_ndays1 = 0
        force_nhours1 = 0
        nleap1 = 0

        nyear1 = end_year - start_year + 1

        IF (nyear1 .LT. 1) THEN
            WRITE(*,*) 'ERROR: invalid simulation period'
            WRITE(*,*) 'start_year=', start_year
            WRITE(*,*) 'end_year=', end_year
            STOP 1
        ENDIF

        IF (ALLOCATED(uni_Year)) DEALLOCATE(uni_Year)
        ALLOCATE(uni_Year(nyear1))

        DO n = 1, nyear1
            uni_Year(n) = start_year + n - 1
        ENDDO

        WRITE(*,*) 'uni_Year', uni_Year, nyear1

        DO n = 1, nyear1

            year = uni_Year(n)
            check_leap1 = check_leap(year)

            IF (check_leap1 .EQ. 1) THEN
                nleap1 = nleap1 + 1
                ndays = 366
            ELSE
                ndays = 365
            ENDIF

            force_nhours1 = force_nhours1 + ndays * 24
            output_ndays1 = output_ndays1 + ndays

            IF (n .EQ. nyear1) THEN
                lastyear_days = ndays
                lastyear_hours = ndays * 24
            ENDIF

        ENDDO

        rere1%nleap_rere = nleap1
        rere1%force_nhours_rere = force_nhours1
        rere1%output_ndays_rere = output_ndays1
        rere1%nyear_rere = nyear1

    END FUNCTION YearDayHour


    FUNCTION Days_cal(simu_year) RESULT(result)
        IMPLICIT NONE
        INTEGER simu_year,m,Days1,Hours1,check_leap1
        TYPE (ResultType) :: result

        Days1 = 365
        
        check_leap1 = check_leap(simu_year)
        IF (check_leap1 .EQ. 1) THEN
            Days1 = 366
        ENDIF

        Hours1 = Days1*24
        result%Days = Days1
        result%Hours = Hours1

    END FUNCTION Days_cal

END MODULE DaysHours

MODULE FileSize
    !------------------------------------------------------------------
    ! Module for managing global variables and subroutines:
    ! 1. `Allosize` - Allocate memory
    ! 2. `Deallo_mod` - Deallocate memory
    ! 3. `FilesForOutput` - Configure output files
    ! (Fangxiu, 20230626)
    !-------------------------------------------------------------------
    USE IntersVariables
    USE DaysHours
    IMPLICIT NONE
    SAVE
    
    REAL (SGL), dimension(:,:), allocatable ::  forcing_data !change for the forcing iiterms
    !REAL, dimension(:),   allocatable ::  output_data  !hourly output record, same raws of force data
    REAL (SGL), dimension(:,:), allocatable ::  output_daily
    REAL (SGL), dimension(:,:), allocatable ::  output_yr
    REAL (SGL), dimension(:,:), allocatable ::  outputd_ccycle_Cpools
    REAL (SGL), dimension(:,:), allocatable ::  outputd_outc
    REAL (SGL), dimension(:,:), allocatable ::  outputd_outn
    REAL (SGL), dimension(:,:), allocatable ::  outputd_outp
    REAL (SGL), dimension(:,:), allocatable ::  outputd_Pcycle_Ppools
    REAL (SGL), dimension(:,:), allocatable ::  outputd_Pcycle_CPratios
    REAL (SGL), dimension(:,:), allocatable ::  outputd_Ncycle_Npools
    REAL (SGL), dimension(:,:), allocatable ::  outputd_Ncycle_CNratios
    REAL (SGL), dimension(:,:), allocatable ::  outputd_Pdynamic
    REAL (SGL), dimension(:,:), allocatable ::  outputd_Ndynamic
    REAL (SGL), dimension(:,:), allocatable ::  output_P_record 
    REAL (SGL), dimension(:,:), allocatable ::  output_N_record
    REAL (SGL), dimension(:,:), allocatable ::  output_OutC_record
    REAL (SGL), dimension(:,:), allocatable ::  output_OutN_record
    REAL (SGL), dimension(:,:), allocatable ::  output_OutP_record

    !check GPP
    REAL (SGL), dimension(:,:), allocatable ::  output_Aleaf_record   !hourly
    REAL (SGL), dimension(:,:), allocatable ::  outputd_Aleaf

    !for SpinUp/SAS xia et al., 2012
    REAL (SGL), dimension(:), allocatable :: LoopNPP_d,LoopNPP_yr,LoopNPP_yr_noadd
    REAL (SGL), dimension(:), allocatable :: LoopNEE_d,LoopNEE_yr
    REAL (SGL) :: DeltNPP,LastLoopNPP,EndLoopNPP
    REAL (SGL), dimension(:,:), allocatable :: LoopQC
    !REAL (SGL), dimension(:,:), allocatable :: LoopQC_yrmean
    REAL (SGL), dimension(:,:), allocatable :: Dynamic_trace_Xct
    REAL (SGL), dimension(:,:), allocatable :: Dynamic_trace_Xpt
    REAL (SGL), dimension(:,:), allocatable :: Dynamic_trace_EcoTau
    REAL (SGL), dimension(:,:), allocatable :: Dynamic_trace_BaseTau
    REAL (SGL), dimension(:,:), allocatable :: Dynamic_trace_deltaQC
    REAL (SGL), dimension(:,:), allocatable :: Rh_pools_record

    REAL (SGL), dimension(:,:), allocatable :: Dynamic_trace_Xct_day
    REAL (SGL), dimension(:,:), allocatable :: Dynamic_trace_Xpt_day
    REAL (SGL), dimension(:,:), allocatable :: Dynamic_trace_EcoTau_day
    REAL (SGL), dimension(:,:), allocatable :: Dynamic_trace_BaseTau_day
    REAL (SGL), dimension(:,:), allocatable :: Dynamic_trace_deltaQC_day


    REAL (SGL), dimension(:), allocatable :: LoopQC9
    REAL (SGL), dimension(:,:), allocatable :: Loop_variables

    REAL (SGL), dimension(:,:), allocatable :: LoopQN
    REAL (SGL), dimension(:,:), allocatable :: LoopQP
    REAL (SGL), dimension(:,:), allocatable :: Loop_variablesNNN
    REAL (SGL), dimension(:,:), allocatable :: Loop_variablesPPP
    
    !    !environmental scalar
    REAL (SGL), dimension(:,:), allocatable :: Somega_output
    REAL (SGL), dimension(:,:), allocatable :: ST_NP_output
    REAL (SGL),dimension(:,:), allocatable :: ksi_record
    REAL (SGL),DIMENSION (:,:),allocatable :: ksi_np_output_daily
    REAL (SGL), dimension(:,:), allocatable ::NPPallo
    REAL (SGL), dimension(:,:), allocatable ::output_nppallo_record

    REAL (SGL), dimension(:,:), allocatable ::  QC_record

    
    ! for MCMC
    REAL (SGL), dimension(:,:), allocatable ::  output_daily_mcmc
    REAL (SGL), dimension(:,:), allocatable ::  obs_TianTong
    !for hourly record
    REAL (SGL), dimension(:,:), allocatable ::  output_mcmc         !hourly for mcmc
    REAL (SGL), dimension(:,:), allocatable ::  output_record        !hourly for no mcmc
    
    INTEGER, PARAMETER :: joutput = 12 !the colnums of output varsiables, daily
    INTEGER, PARAMETER :: jrecord = 18 !the colnums of hourly simitations
    INTEGER, PARAMETER :: jforce=11    !force_nhours=8760, iforce = force_nhours
    INTEGER, PARAMETER :: jobs=5       !the colnums of the observation for MCMC
    INTEGER, PARAMETER :: jmcmc=3      !the colnums of the constraind results after MCMC
 
    CONTAINS
    
    SUBROUTINE AlloSize()

    ! Subroutine for allocate size for variables

    IMPLICIT NONE
    INTEGER m

    !Consdering the leap year
        ispinup = 0
        rere = YearDayHour()
        force_nhours = rere%force_nhours_rere
        output_ndays = rere%output_ndays_rere
        nyear        = rere%nyear_rere

        ! ALLO _REC 
        allocate(MaxUpCapP_rec(force_nhours))
        allocate(bP_rec(force_nhours))
        allocate(aP_rec(force_nhours))
        allocate(fwsoil_rec(force_nhours))
        allocate(omega_rec(force_nhours))
        allocate(P_loss_rec(force_nhours))
        allocate(Fptase_rec(force_nhours))
        allocate(asorb_rec(force_nhours))
        allocate(P_miner2_rec(force_nhours))
        allocate(P_immob2_rec(force_nhours))
        allocate(P_net2_rec(force_nhours))
        allocate(Dsoillab_rec(force_nhours))
        allocate(Dsoillab_rec2(force_nhours))
        allocate(P_uptake_rec(force_nhours))
        allocate(P_in_out_rec(force_nhours))
        allocate(FdiffP_rec(force_nhours))
        allocate(delta_QLabP_root_rec(force_nhours))
        allocate(D_rec(force_nhours))
        allocate(swc_p_rec(force_nhours))
        allocate(tf_rec(force_nhours))
        allocate(rdiff_rec(force_nhours))
        allocate(tf1_rec(force_nhours))
        allocate(tf2_rec(force_nhours))
        allocate(f_rootp_rec(force_nhours))
        allocate(QC3_rec(force_nhours))
        allocate(LfactorP_rec(force_nhours))
        allocate(P_deficit_rec(force_nhours))
        allocate(P_deficit_rec_test(force_nhours))
        allocate(P_demand_rec(force_nhours))
        allocate(P_demand_rec_test(force_nhours))
        allocate(N_demand_rec(force_nhours))
        allocate(N_uptake_rec(force_nhours))
        allocate(N_deficit_rec(force_nhours))
        allocate(N_fixation_x_rec(force_nhours))
        allocate(P_leach_rec(force_nhours))
        allocate(Xss_new_step_rec(force_nhours,9))


        allocate(forcing_data(11,force_nhours))
        write(*,*)'forcing_data size:',shape(forcing_data),nyear

        nobs = output_ndays  !nobs = the raws of observation

        !output
        !allocate(output_data(force_nhours))
        allocate(output_daily(5,output_ndays))  !note the dimension1 of output_daily isnt the joutput!
        allocate(output_yr(14,nyear))
        allocate(outputd_outc(pools,output_ndays))
        allocate(outputd_outn(pools,output_ndays))
        allocate(outputd_outp(pools,output_ndays))
        allocate(outputd_ccycle_Cpools(pools,output_ndays)) !pools -- 8 carbon pools
        allocate(outputd_Pcycle_Ppools(pools,output_ndays))
        allocate(outputd_Ncycle_Npools(pools,output_ndays))
        allocate(outputd_Pcycle_CPratios(pools,output_ndays))
        allocate(outputd_Ncycle_CNratios(pools,output_ndays))
        allocate(outputd_Pdynamic(9,output_ndays))
        allocate(outputd_Ndynamic(7,output_ndays))
      
        allocate(QC_record(pools,force_nhours))
        allocate(Dynamic_trace_Xct(pools,force_nhours))
        allocate(Dynamic_trace_Xpt(pools,force_nhours))
        allocate(Dynamic_trace_EcoTau(pools,force_nhours))
        allocate(Dynamic_trace_BaseTau(pools,force_nhours))
        allocate(Dynamic_trace_deltaQC(pools,force_nhours))
        allocate(Rh_pools_record(5,force_nhours))

        allocate(Dynamic_trace_Xct_day(pools,output_ndays))
        allocate(Dynamic_trace_Xpt_day(pools,output_ndays))
        allocate(Dynamic_trace_EcoTau_day(pools,output_ndays))
        allocate(Dynamic_trace_BaseTau_day(pools,output_ndays))
        allocate(Dynamic_trace_deltaQC_day(pools,output_ndays))

        !for spinup/SAS
        IF (NDSPINUP .EQ. 1)  THEN
            allocate(LoopNPP_d(nspinup))
            allocate(LoopNPP_yr(nspinup))
            allocate(LoopNPP_yr_noadd(nspinup))
            allocate(LoopNEE_d(nspinup))
            allocate(LoopNEE_yr(nspinup))
            allocate(LoopQC(pools,nspinup))
            !allocate(LoopQC_yrmean(pools,nspinup))
            allocate(LoopQC9(nspinup))
            allocate(Loop_variables(14,nspinup))
            allocate(LoopQN(14,nspinup))
            allocate(Loop_variablesNNN(14,nspinup))
            allocate(LoopQP(17,nspinup))
            allocate(Loop_variablesPPP(22,nspinup))
            !allocate(Limite_NP_loops(7,nspinup))
        ENDIF
        
        allocate(Somega_output(24,output_ndays))
        allocate(ST_NP_output(15,output_ndays))  !wan 0519 change to force_nhours
        allocate(NPPallo(4,output_ndays))
        allocate(ksi_np_output_daily(1,output_ndays))
        
        

        !for GPP modify
        allocate(output_Aleaf_record(12,force_nhours))
        allocate(outputd_Aleaf(5,output_ndays))
    
        !for mcmc
        allocate(output_daily_mcmc(jmcmc,nobs))
        !allocate(obs_TianTong(jobs,nobs))
        !allocate(obs_TianTong(4,8760))
        allocate(obs_TianTong(5,nobs))
        obs_TianTong = -999

        !for hourly record
        allocate(ksi_record(6,force_nhours))

        allocate(output_record(12,force_nhours))
        allocate(output_mcmc(jmcmc,force_nhours))
        allocate(output_P_record(5,force_nhours))
        allocate(output_N_record(7,force_nhours))
        allocate(output_nppallo_record(4,force_nhours))
        allocate(output_OutC_record(pools,force_nhours))
        allocate(output_OutN_record(pools,force_nhours))
        allocate(output_OutP_record(pools,force_nhours))

        write(*,*)'allsize called sucessfully'
        print *, shape(outputd_Ndynamic),shape(forcing_data),shape(output_daily)
        write(*,*) 'force_nhours,output_ndays',force_nhours,output_ndays

    END SUBROUTINE AlloSize


    SUBROUTINE Deallo_mod
    !-----------------------------------
    !Subroutine for deallocate the space
    !-----------------------------------
    IMPLICIT NONE
        deallocate(forcing_data)
        !deallocate(output_data)
        deallocate(output_daily)  !note the dimension1 of output_daily isnt the joutput!
        deallocate(outputd_outc)
        deallocate(outputd_ccycle_Cpools) !8 -- 8 carbon pools
        deallocate(outputd_Pcycle_Ppools)
        deallocate(outputd_Ncycle_Npools)
        deallocate(outputd_Pcycle_CPratios)
        deallocate(outputd_Ncycle_CNratios)
        deallocate(outputd_Pdynamic)
        deallocate(outputd_Ndynamic)
        
        !for spinup/SAS
        deallocate(LoopNPP_d)
        deallocate(LoopNPP_yr)
        deallocate(LoopNPP_yr_noadd)
        deallocate(LoopNEE_d)
        deallocate(LoopNEE_yr)
        deallocate(LoopQC)
        deallocate(LoopQC9)
        deallocate(Loop_variables)


        deallocate(LoopQN)
        deallocate(Loop_variablesNNN)
        deallocate(LoopQP)
        deallocate(Loop_variablesPPP)
        

        deallocate(Somega_output)
        deallocate(ST_NP_output)  !wan 0519 change to force_nhours
        deallocate(NPPallo) 

        !for GPP modify
        deallocate(output_Aleaf_record)
        deallocate(outputd_Aleaf)
    
        !for mcmc
        deallocate(output_daily_mcmc)
        deallocate(obs_TianTong)

        !for hourly record
        deallocate(output_record)
        deallocate(output_mcmc)
        deallocate(output_P_record)
        deallocate(output_N_record)
        deallocate(output_nppallo_record)
        deallocate(output_OutC_record)

        write(*,*)'deallsize called sucessfully'

    END SUBROUTINE Deallo_mod

END MODULE FileSize

MODULE outputs_mod
    USE FileSize
    CONTAINS
    SUBROUTINE ForRecords(AcanL)

        IMPLICIT NONE

        INTEGER :: m
        REAL :: Ta
        REAL :: AcanL(5)
        INTEGER, PARAMETER :: cc=1

        !===========================================================
        ! MCMC mode
        !===========================================================
        IF(MCMC .EQ. 1) THEN

            !-------------------------------------------------------
            ! Hourly carbon fluxes
            !-------------------------------------------------------
            output_mcmc(1,itime)=GPP
            output_mcmc(2,itime)=Reco
            output_mcmc(3,itime)=NEE

            ! Carbon pools and litter inputs
            QC_record(:,itime)=QC

            output_OutC_record(1:pools,itime)=OutC(1:pools)
            output_OutN_record(1:pools,itime)=OutN(1:pools)
            output_OutP_record(1:pools,itime)=OutP(1:pools)

            !-------------------------------------------------------
            ! Hourly P-cycle variables
            !-------------------------------------------------------
            IF(CYCLE_CNP .GT. 2) THEN

                output_P_record(1,itime)=Dsoillab
                output_P_record(2,itime)=P_net2
                output_P_record(3,itime)=Fptase
                output_P_record(4,itime)=P_loss
                output_P_record(5,itime)=P_uptake

            ENDIF

            !-------------------------------------------------------
            ! Daily outputs
            !-------------------------------------------------------
            IF(MOD(itime,24) .EQ. 0) THEN

                ! Global day index increases exactly once per day
                idays=idays+cc

                IF(idays .GT. output_ndays) THEN
                    WRITE(*,*) 'ERROR: idays exceeds output_ndays'
                    WRITE(*,*) 'itime=',itime
                    WRITE(*,*) 'idays=',idays
                    WRITE(*,*) 'output_ndays=',output_ndays
                    STOP
                ENDIF

                ! Daily GPP, Reco and NEE
                output_daily_mcmc(1,idays)=                 &
                    SUM(output_mcmc(1,(itime-23):itime))

                output_daily_mcmc(2,idays)=                 &
                    SUM(output_mcmc(2,(itime-23):itime))

                output_daily_mcmc(3,idays)=                 &
                    SUM(output_mcmc(3,(itime-23):itime))

                ! Daily instantaneous carbon pools
                outputd_ccycle_Cpools(1:pools,idays)=QC

                ! Daily litter input
                DO m=1,pools

                    outputd_outc(m,idays)=                  &
                        SUM(output_OutC_record(             &
                        m,(itime-23):itime))

                    outputd_outn(m,idays)=                  &
                        SUM(output_OutN_record(             &
                        m,(itime-23):itime))

                    outputd_outp(m,idays)=                  &
                        SUM(output_OutP_record(             &
                        m,(itime-23):itime))

                ENDDO

                ! Daily P variables and pools
                IF(CYCLE_CNP .GT. 2) THEN

                    outputd_Pdynamic(1,idays)=              &
                        SUM(output_P_record(                &
                        1,(itime-23):itime))

                    outputd_Pdynamic(2,idays)=              &
                        SUM(output_P_record(                &
                        2,(itime-23):itime))

                    outputd_Pdynamic(3,idays)=              &
                        SUM(output_P_record(                &
                        3,(itime-23):itime))

                    outputd_Pdynamic(4,idays)=              &
                        SUM(output_P_record(                &
                        4,(itime-23):itime))

                    outputd_Pdynamic(5,idays)=              &
                        SUM(output_P_record(                &
                        5,(itime-23):itime))

                    ! Soil inorganic P pools
                    outputd_Pdynamic(6,idays)=QPlab
                    outputd_Pdynamic(7,idays)=QPsorb
                    outputd_Pdynamic(8,idays)=QPss
                    outputd_Pdynamic(9,idays)=QPocc

                    ! Ecosystem P pools and C:P ratios
                    outputd_Pcycle_Ppools(1:pools,idays)=QP
                    outputd_Pcycle_CPratios(1:pools,idays)=CP

                ENDIF

            ENDIF

        !===========================================================
        ! Normal forward-simulation mode
        !===========================================================
        ELSE

            !-------------------------------------------------------
            ! Hourly ecosystem outputs
            !-------------------------------------------------------
            output_record(1,itime)=GPP
            output_record(2,itime)=Reco
            output_record(3,itime)=NEE
            output_record(4,itime)=LAI
            output_record(5,itime)=bmplant
            output_record(6,itime)=NPP
            output_record(7,itime)=Rhetero
            output_record(8,itime)=Rauto
            output_record(9,itime)=bmleaf
            output_record(10,itime)=bmstem
            output_record(11,itime)=bmroot
            output_record(12,itime)=bmRe

            Rh_pools_record(1:5,itime)=Rh_pools(1:5)

            ! NPP allocation
            output_nppallo_record(1,itime)=alpha_L
            output_nppallo_record(2,itime)=alpha_W
            output_nppallo_record(3,itime)=alpha_R
            output_nppallo_record(4,itime)=alpha_Re

            ! Leaf photosynthesis
            output_Aleaf_record(1,itime)=Aleaf(1)
            output_Aleaf_record(2,itime)=Aleaf(2)

            !-------------------------------------------------------
            ! Hourly N-cycle variables
            !-------------------------------------------------------
            IF(CYCLE_CNP .GT. 1) THEN

                output_N_record(1,itime)=N_leach
                output_N_record(2,itime)=N_vol
                output_N_record(3,itime)=N_net2
                output_N_record(4,itime)=N_uptake
                output_N_record(5,itime)=N_fixation
                output_N_record(6,itime)=N_transfer
                output_N_record(7,itime)=QNminer

            ENDIF

            !-------------------------------------------------------
            ! Hourly P-cycle variables
            !-------------------------------------------------------
            IF(CYCLE_CNP .GT. 2) THEN

                output_P_record(1,itime)=Dsoillab
                output_P_record(2,itime)=P_net2
                output_P_record(3,itime)=Fptase
                output_P_record(4,itime)=P_loss
                output_P_record(5,itime)=P_uptake

            ENDIF

            !-------------------------------------------------------
            ! Hourly pools and litter input
            !-------------------------------------------------------
            QC_record(:,itime)=QC

            output_OutC_record(1:pools,itime)=OutC(1:pools)
            output_OutN_record(1:pools,itime)=OutN(1:pools)
            output_OutP_record(1:pools,itime)=OutP(1:pools)

            !-------------------------------------------------------
            ! Daily output
            !-------------------------------------------------------
            IF(MOD(itime,24) .EQ. 0) THEN

                ! Global day index increases exactly once per day
                idays=idays+cc

                IF(idays .GT. output_ndays) THEN
                    WRITE(*,*) 'ERROR: idays exceeds output_ndays'
                    WRITE(*,*) 'itime=',itime
                    WRITE(*,*) 'idays=',idays
                    WRITE(*,*) 'output_ndays=',output_ndays
                    STOP
                ENDIF

                output_daily(1,idays)=                     &
                    SUM(output_record(1,(itime-23):itime))

                output_daily(2,idays)=                     &
                    SUM(output_record(2,(itime-23):itime))

                output_daily(3,idays)=                     &
                    SUM(output_record(3,(itime-23):itime))

                output_daily(4,idays)=                     &
                    SUM(output_record(6,(itime-23):itime))

                output_daily(5,idays)=                     &
                    SUM(output_record(4,(itime-23):itime))/24.0

                ! Carbon pools
                outputd_ccycle_Cpools(1:pools,idays)=QC

                ! Daily litter input
                DO m=1,pools

                    outputd_outc(m,idays)=                 &
                        SUM(output_OutC_record(             &
                        m,(itime-23):itime))

                    outputd_outn(m,idays)=                 &
                        SUM(output_OutN_record(             &
                        m,(itime-23):itime))

                    outputd_outp(m,idays)=                 &
                        SUM(output_OutP_record(             &
                        m,(itime-23):itime))

                ENDDO

                !---------------------------------------------------
                ! Environmental scalars
                !---------------------------------------------------
                Somega_output(1,idays)=S_omega
                Somega_output(2,idays)=W
                Somega_output(3,idays)=Sw
                Somega_output(4:13,idays)=wsc(:)
                Somega_output(14:23,idays)=wcl(:)
                Somega_output(24,idays)=omega

                ST_NP_output(1:5,idays)=S_t
                ST_NP_output(6,idays)=St
                ST_NP_output(7,idays)=xde
                ST_NP_output(8,idays)=xnp_leaf_limit
                ST_NP_output(9,idays)=xnp_stem_limit
                ST_NP_output(10,idays)=xnp_root_limit
                ST_NP_output(11,idays)=xNPuptake
                ST_NP_output(12,idays)=LfactorN
                ST_NP_output(13,idays)=LfactorP

                ! ST_NP_output(14,idays)=                   &
                !     SUM(output_N_record(                  &
                !     9,(itime-23):itime))

                ! ST_NP_output(15,idays)=                   &
                !     SUM(output_P_record(                  &
                !     10,(itime-23):itime))

                !---------------------------------------------------
                ! Daily NPP allocation
                !---------------------------------------------------
                NPPallo(1,idays)=                         &
                    SUM(output_nppallo_record(            &
                    1,(itime-23):itime))/24.0

                NPPallo(2,idays)=                         &
                    SUM(output_nppallo_record(            &
                    2,(itime-23):itime))/24.0

                NPPallo(3,idays)=                         &
                    SUM(output_nppallo_record(            &
                    3,(itime-23):itime))/24.0

                NPPallo(4,idays)=                         &
                    SUM(output_nppallo_record(            &
                    4,(itime-23):itime))/24.0

                !---------------------------------------------------
                ! Daily N cycle
                !---------------------------------------------------
                IF(CYCLE_CNP .GT. 1) THEN

                    outputd_Ndynamic(1,idays)=            &
                        SUM(output_N_record(              &
                        1,(itime-23):itime))

                    outputd_Ndynamic(2,idays)=            &
                        SUM(output_N_record(              &
                        2,(itime-23):itime))

                    outputd_Ndynamic(3,idays)=            &
                        SUM(output_N_record(              &
                        3,(itime-23):itime))

                    outputd_Ndynamic(4,idays)=            &
                        SUM(output_N_record(              &
                        4,(itime-23):itime))

                    outputd_Ndynamic(5,idays)=            &
                        SUM(output_N_record(              &
                        5,(itime-23):itime))

                    outputd_Ndynamic(6,idays)=            &
                        SUM(output_N_record(              &
                        6,(itime-23):itime))

                    outputd_Ndynamic(7,idays)=            &
                        SUM(output_N_record(              &
                        7,(itime-23):itime))

                    outputd_ncycle_Npools(1:pools,idays)=QN
                    outputd_Ncycle_CNratios(1:pools,idays)=CN

                ENDIF

                !---------------------------------------------------
                ! Daily P cycle
                !---------------------------------------------------
                IF(CYCLE_CNP .GT. 2) THEN

                    outputd_Pdynamic(1,idays)=            &
                        SUM(output_P_record(              &
                        1,(itime-23):itime))

                    outputd_Pdynamic(2,idays)=            &
                        SUM(output_P_record(              &
                        2,(itime-23):itime))

                    outputd_Pdynamic(3,idays)=            &
                        SUM(output_P_record(              &
                        3,(itime-23):itime))

                    outputd_Pdynamic(4,idays)=            &
                        SUM(output_P_record(              &
                        4,(itime-23):itime))

                    outputd_Pdynamic(5,idays)=            &
                        SUM(output_P_record(              &
                        5,(itime-23):itime))

                    outputd_Pdynamic(6,idays)=QPlab
                    outputd_Pdynamic(7,idays)=QPsorb
                    outputd_Pdynamic(8,idays)=QPss
                    outputd_Pdynamic(9,idays)=QPocc

                    outputd_Pcycle_Ppools(1:pools,idays)=QP
                    outputd_Pcycle_CPratios(1:pools,idays)=CP

                ENDIF

            ENDIF

        ENDIF

    END SUBROUTINE ForRecords
    
    SUBROUTINE ForYearRecords(start_day,end_day,simu_year)
        IMPLICIT NONE

        REAL GPP_yr,Reco_yr,NEE_yr,NPP_yr,NEEnew_yr
        REAL Rhetero_yr,Rauto_yr,Rmain_yr,Rnitrogen_yr,Rphosphorus_yr,Rauto_new_yr
        INTEGER start_day,end_day,simu_year
        REAL QC_yr(9)

        GPP_yr = sum(output_daily(1,start_day:end_day))
        Reco_yr= sum(output_daily(2,start_day:end_day)) 
        NEE_yr = sum(output_daily(3,start_day:end_day))
        NPP_yr = sum(output_daily(7,start_day:end_day))
        NEEnew_yr = sum(output_daily(28,start_day:end_day))

        QC_yr(:)= outputd_ccycle_Cpools(:,end_day) 

        Rhetero_yr = sum(output_daily(14,start_day:end_day))
        Rauto_yr  = sum(output_daily(15,start_day:end_day))
        Rauto_new_yr  = sum(output_daily(27,start_day:end_day))
        Rmain_yr  = sum(output_daily(16,start_day:end_day))
        Rnitrogen_yr = sum(output_daily(18,start_day:end_day))
        Rphosphorus_yr = sum(output_daily(19,start_day:end_day))

        output_yr(1,iyear) = simu_year
        output_yr(2,iyear) = GPP_yr
        output_yr(3,iyear) = Reco_yr
        output_yr(4,iyear) = NEE_yr
        output_yr(5,iyear) = NPP_yr
        output_yr(6:14,iyear) = QC_yr(:)

    END SUBROUTINE ForYearRecords


    SUBROUTINE Filespath()

     ! Define the path for output based on simulation type

        IMPLICIT NONE
        CHARACTER(len=99) fold

        ! Determine the folder path based on the simulation conditions
        IF (MCMC .EQ. 0) THEN
            IF (NDSPINUP .EQ. 1) THEN
                fold = 'NDSPINUP/' ! Path for the spin-up procedure
            ELSEIF (SensTest .EQ. 0 .and. MCMC .EQ. 0 .and. &
                &   NDSPINUP .EQ. 0 )THEN
                fold = 'sim/'      ! Path for normal simulation
            ELSEIF (SensTest .EQ. 1 .and. MCMC .EQ. 0 .and. &
                &   NDSPINUP .EQ. 0 )THEN
                fold = 'SensTest/' ! Path for sensitivity testing
            ENDIF
        ELSE
            fold = 'MCMC/'         ! Path for MCMC simulations
        ENDIF

        IF(CYCLE_CNP .eq. 1) THEN
        fileplace = '../output/'//TRIM(ADJUSTL(fold))//'teco_c/'
        filesignal = '_c.csv'
        ELSEIF(CYCLE_CNP .eq. 2)THEN
        fileplace = '../output/'//TRIM(ADJUSTL(fold))//'teco_cn/'
        filesignal = '_cn.csv'
        ELSEIF(CYCLE_CNP .eq. 3)THEN
        fileplace = '../output/'//TRIM(ADJUSTL(fold))//'teco_cnp/'
        filesignal='_cnp.csv'
        ENDIF


    END SUBROUTINE Filespath
    
    
    SUBROUTINE FilesForOutput()

    ! Subroutine for configuring the file path for model output.

        IMPLICIT NONE

        OPEN(211, file=TRIM(ADJUSTL(fileplace))//'teco_simu_cflux'//filesignal) 
        OPEN(212, file=TRIM(ADJUSTL(fileplace))//'teco_8cpools'//filesignal)
        OPEN(213, file=TRIM(ADJUSTL(fileplace))//'teco_outc'//filesignal)
        OPEN(21313, file=TRIM(ADJUSTL(fileplace))//'teco_outn'//filesignal)
        OPEN(2131313, file=TRIM(ADJUSTL(fileplace))//'teco_outp'//filesignal)
        OPEN(214, file=TRIM(ADJUSTL(fileplace))//'NPP_allocation_fraction'//filesignal)
        OPEN(998, file=TRIM(ADJUSTL(fileplace))//'S_omega_record'//filesignal)
        OPEN(997, file=TRIM(ADJUSTL(fileplace))//'ST_NP'//filesignal)
        OPEN(996, file=TRIM(ADJUSTL(fileplace))//'Rh_pools'//filesignal) 
        open(215, file=TRIM(ADJUSTL(fileplace))//'Aleaf'//filesignal)
        OPEN (191919,file=TRIM(ADJUSTL(fileplace))//'YearOutPut'//filesignal)
        open(8881,file=TRIM(ADJUSTL(fileplace))//'QC_record'//filesignal)

        !CHANGE the output filepath for different configuration

        IF(CYCLE_CNP .gt. 1)THEN
            OPEN(771,file=TRIM(ADJUSTL(fileplace))//'CNP_Ncycle_simu_pools'//filesignal)
            OPEN(772,file=TRIM(ADJUSTL(fileplace))//'CNP_Ncycle_simu_ratios'//filesignal)
            OPEN(773,file=TRIM(ADJUSTL(fileplace))//'CNP_NNNcycle_dynamics'//filesignal)
            OPEN(774,file=TRIM(ADJUSTL(fileplace))//'Limit_NP_record'//filesignal)
            OPEN(775,file=TRIM(ADJUSTL(fileplace))//'Limit_NP_loops'//filesignal)
            open(665,file=TRIM(ADJUSTL(fileplace))//'plant_NP_allo_record'//filesignal)
        ENDIF

        IF(CYCLE_CNP .gt. 2)THEN

            OPEN(661,file=TRIM(ADJUSTL(fileplace))//'CNP_Pcycle_simu_pools'//filesignal)
            OPEN(662,file=TRIM(ADJUSTL(fileplace))//'CNP_Pcycle_simu_ratios'//filesignal)
            OPEN(663,file=TRIM(ADJUSTL(fileplace))//'CNP_PPPcycle_dynamics'//filesignal)
        ENDIF
        
        IF(NDSPINUP .EQ. 1) THEN
            OPEN(999,file=TRIM(ADJUSTL(fileplace))//'spinup_LoopVariables'//filesignal)
            IF (CYCLE_CNP .GT. 1) OPEN(9991,file=TRIM(ADJUSTL(fileplace))//'spinup_LoopVariablesNNN'//filesignal)
            IF (CYCLE_CNP .GT. 2) OPEN(9992,file=TRIM(ADJUSTL(fileplace))//'spinup_LoopVariablesPPP'//filesignal)
        ENDIF
        
        open(1212,file=TRIM(ADJUSTL(fileplace))//'simu_cflux_record'//filesignal)

        IF(nyear .gt. 1 .OR. (nyear .eq. 1 .AND. NDSPINUP .eq. 1 ))THEN
            open(216,file=TRIM(ADJUSTL(fileplace))//'simu_cflux_lastyear'//filesignal)
            open(12,file=TRIM(ADJUSTL(fileplace))//'simu_cflux_lastyear_record'//filesignal)
            open(217,file=TRIM(ADJUSTL(fileplace))//'simu_8cpools_lastyear'//filesignal)
            IF(CYCLE_CNP .gt. 1)THEN
                OPEN(7711,file=TRIM(ADJUSTL(fileplace))//'CNP_Ncycle_simu_pools_lastyear'//filesignal)
                OPEN(7722,file=TRIM(ADJUSTL(fileplace))//'CNP_Ncycle_simu_ratios_lastyear'//filesignal)
                OPEN(7733,file=TRIM(ADJUSTL(fileplace))//'CNP_NNNcycle_dynamics_lastyear'//filesignal)
            ENDIF
            IF(CYCLE_CNP .gt. 2)THEN
                OPEN(6611,file=TRIM(ADJUSTL(fileplace))//'CNP_Pcycle_simu_pools_lastyear'//filesignal)
                OPEN(6622,file=TRIM(ADJUSTL(fileplace))//'CNP_Pcycle_simu_ratios_lastyear'//filesignal)
                OPEN(6633,file=TRIM(ADJUSTL(fileplace))//'CNP_PPPcycle_dynamics_lastyear'//filesignal)
            ENDIF
        ENDIF

    END SUBROUTINE FilesForOutput

    SUBROUTINE WriteFiles_noMCMC()
        
    ! Customize the outputs

        IMPLICIT NONE
        INTEGER m,n
        DO m=1,nyear
            WRITE(191919,'(*( G0.8, :, ",", X))')output_yr(:,m)
        ENDDO

        DO m=1,output_ndays
            WRITE(211,'(*( G0.8, :, ",", X))')output_daily(:,m)
            WRITE(212,'(*( G0.8, :, ",", X))')outputd_ccycle_Cpools(:,m)                 

            WRITE(213,'(*( G0.8, :, ",", X))')outputd_outc(:,m)
            WRITE(21313,'(*( G0.8, :, ",", X))')outputd_outn(:,m)
            WRITE(2131313,'(*( G0.8, :, ",", X))')outputd_outp(:,m)
            WRITE(214,'(*( G0.8, :, ",", X))')NPPallo(:,m)
            WRITE(998,'(*( G0.8, :, ",", X))')Somega_output(:,m)
            WRITE(997,'(*( G0.8, :, ",", X))')ST_NP_output(:,m)

            IF(CYCLE_CNP .gt. 1)THEN
                WRITE(771,'(*( G0.8, :, ",", X))')outputd_ncycle_Npools(:,m)
                WRITE(772,'(*( G0.8, :, ",", X))')outputd_Ncycle_CNratios(:,m)    
                WRITE(773,'(*( G0.8, :, ",", X))')outputd_Ndynamic(:,m)
            ENDIF
            IF(CYCLE_CNP .gt. 2)THEN
                WRITE(661,'(*( G0.8, :, ",", X))')outputd_Pcycle_Ppools(:,m)
                WRITE(662,'(*( G0.8, :, ",", X))')outputd_Pcycle_CPratios(:,m)
                WRITE(663,'(*( G0.8, :, ",", X))')outputd_Pdynamic(:,m)
            ENDIF
        ENDDO

        DO n=1,force_nhours
            WRITE(8881,'(*( G0.8, :, ",", X))')QC_record(:,n)
            WRITE(1212,'(*( G0.8, :, ",", X))')output_record(:,n)
            WRITE(996,'(*( G0.8, :, ",", X))')Rh_pools_record(:,n)
            WRITE(215,'(*( G0.8, :, ",", X))')output_Aleaf_record(:,n)
        ENDDO

    IF(nyear .gt. 1 .OR. (nyear .eq. 1 .AND. NDSPINUP .eq. 1 ))THEN
    ! if we have long-term forcing, we want to keep the results of the last year for our convenient 
        DO n =output_ndays-lastyear_days+1,output_ndays
            WRITE(216,'(*( G0.8, :, ",", X))')output_daily(:,n)
            WRITE(217,'(*( G0.8, :, ",", X))')outputd_ccycle_Cpools(:,n)
            IF(CYCLE_CNP .gt. 1)THEN
                    WRITE(7711,'(*( G0.8, :, ",", X))')outputd_ncycle_Npools(:,n)
                    WRITE(7722,'(*( G0.8, :, ",", X))')outputd_Ncycle_CNratios(:,n)    
                    WRITE(7733,'(*( G0.8, :, ",", X))')outputd_Ndynamic(:,n)
            ENDIF
            IF(CYCLE_CNP .gt. 2)THEN
                    WRITE(6611,'(*( G0.8, :, ",", X))')outputd_Pcycle_Ppools(:,n)
                    WRITE(6622,'(*( G0.8, :, ",", X))')outputd_Pcycle_CPratios(:,n)
                    WRITE(6633,'(*( G0.8, :, ",", X))')outputd_Pdynamic(:,n)
            ENDIF 

        ENDDO

        DO n= force_nhours-lastyear_hours+1,force_nhours
            WRITE(12,'(*( G0.8, :, ",", X))')output_record(:,n)
        ENDDO
    ENDIF

    IF (NDSPINUP .EQ. 1 ) THEN
        DO m=1,ispinup   !or nspinup-1?? or  = ispinup ??
            WRITE(999,'(*( G0.8, :, ",", X))') Loop_variables(:,m)
            IF (CYCLE_CNP .GT. 1) WRITE(9991,'(*( G0.8, :, ",", X))')Loop_variablesNNN(:,m)
            IF (CYCLE_CNP .GT. 2) WRITE(9992,'(*( G0.8, :, ",", X))')Loop_variablesPPP(:,m)
        ENDDO
    ENDIF  ! END the output when SpinUp is activated

<<<<<<< /home/user/fangwei/FangxiuWan-TECO-CNP-Sv1.0-edff198/TECO-CNP Sv1.0 20250316/main/FileSize.f90
    END SUBROUTINE WriteFiles_noMCMC
END MODULE outputs_mod
||||||| /home/user/fangwei/TECO_sensitivity_server/baseline_main/FileSize.f90
    END SUBROUTINE WriteFiles_noMCMC

END MODULE outputs_mod
=======
    END SUBROUTINE WriteFiles_noMCMC

    SUBROUTINE WriteSensitivitySummary()
        IMPLICIT NONE
        INTEGER :: u, ios, nday
        CHARACTER(len=1024) :: outfile
        REAL :: npp_leaf_sum, npp_wood_sum, npp_root_sum, npp_repro_sum

        CALL get_environment_variable('TECO_SUMMARY_FILE', outfile, status=ios)
        IF (ios .NE. 0 .OR. LEN_TRIM(outfile) .EQ. 0) THEN
            outfile = 'teco_sensitivity_summary.csv'
        ENDIF

        nday = output_ndays
        IF (nday .LT. 1) THEN
            WRITE(*,*) 'ERROR: no daily output available for sensitivity summary'
            STOP 91
        ENDIF

        npp_leaf_sum  = SUM(output_daily(4,:)*NPPallo(1,:))
        npp_wood_sum  = SUM(output_daily(4,:)*NPPallo(2,:))
        npp_root_sum  = SUM(output_daily(4,:)*NPPallo(3,:))
        npp_repro_sum = SUM(output_daily(4,:)*NPPallo(4,:))

        OPEN(newunit=u, file=TRIM(outfile), status='replace', action='write', iostat=ios)
        IF (ios .NE. 0) THEN
            WRITE(*,*) 'ERROR: cannot create sensitivity summary: ', TRIM(outfile)
            STOP 92
        ENDIF

        WRITE(u,'(A)') 'gpp_sum,npp_sum,lai_mean,reco_sum,rauto_sum,rhetero_sum,'// &
            'npp_leaf_sum,npp_wood_sum,npp_root_sum,npp_repro_sum,litterfall_c_sum,'// &
            'leaf_c_end,wood_c_end,root_c_end,plant_c_end,soil_c_end,nsc_end,'// &
            'p_demand_sum,p_uptake_sum,p_retrans_sum,leaf_p_end,root_p_end,plant_p_end,'// &
            'leaf_cp_end,root_cp_end,p_net_sum,fptase_sum,p_loss_sum,'// &
            'labile_p_end,sorbed_p_end,strongly_sorbed_p_end,occluded_p_end,'// &
            'microbial_p_end,organic_p_end,total_inorganic_p_end'

        WRITE(u,'(*(G0.10,:,","))') &
            SUM(output_daily(1,:)), SUM(output_daily(4,:)), &
            SUM(output_daily(5,:))/REAL(nday), SUM(output_daily(2,:)), &
            SUM(output_record(8,:)), SUM(output_record(7,:)), &
            npp_leaf_sum, npp_wood_sum, npp_root_sum, npp_repro_sum, &
            SUM(outputd_outc(1:4,:)), &
            outputd_ccycle_Cpools(1,nday), outputd_ccycle_Cpools(2,nday), &
            outputd_ccycle_Cpools(3,nday), SUM(outputd_ccycle_Cpools(1:4,nday)), &
            SUM(outputd_ccycle_Cpools(5:9,nday)), NSC, &
            SUM(P_demand_rec), SUM(outputd_Pdynamic(5,:)), &
            SUM((outputd_outp(1,:)+outputd_outp(2,:)+outputd_outp(3,:))*alphaP), &
            outputd_Pcycle_Ppools(1,nday), outputd_Pcycle_Ppools(3,nday), &
            SUM(outputd_Pcycle_Ppools(1:4,nday)), &
            outputd_Pcycle_CPratios(1,nday), outputd_Pcycle_CPratios(3,nday), &
            SUM(outputd_Pdynamic(2,:)), SUM(outputd_Pdynamic(3,:)), &
            SUM(outputd_Pdynamic(4,:)), outputd_Pdynamic(6,nday), &
            outputd_Pdynamic(7,nday), outputd_Pdynamic(8,nday), &
            outputd_Pdynamic(9,nday), outputd_Pcycle_Ppools(7,nday), &
            SUM(outputd_Pcycle_Ppools(5:9,nday)), SUM(outputd_Pdynamic(6:9,nday))

        CLOSE(u)
        WRITE(*,*) 'Sensitivity summary written: ', TRIM(outfile)
    END SUBROUTINE WriteSensitivitySummary
END MODULE outputs_mod
>>>>>>> /home/user/fangwei/TECO_sensitivity_server/patched_main/FileSize.f90






