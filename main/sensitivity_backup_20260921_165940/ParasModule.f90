!-----------------------------

MODULE update_traits
    IMPLICIT NONE
    CONTAINS
    REAL FUNCTION calc_vcmax(t, Vcmax0)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: t !iitime, update every hour
    REAL, INTENT(IN) :: Vcmax0
    REAL :: k, vcmax_2011, vcmax_2020
    INTEGER t1

    ! parameters
    vcmax_2011 = Vcmax0  ! initial
    vcmax_2020 = Vcmax0 * (1+0.28)     ! end（+response ratio）
    k = abs(vcmax_2020-vcmax_2011)/10         ! rate
    ! cal
    t1 = t - 2010
    calc_vcmax = vcmax_2011 + (vcmax_2020 - vcmax_2011) * (1.0 - EXP(-k * t1))

    END FUNCTION calc_vcmax

    REAL FUNCTION calc_sla(t, SLAx)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: t !simu_year, update annually
        REAL, INTENT(IN) :: SLAx    ! 
        REAL :: k, sla_2011,sla_2020
        INTEGER t1
        ! parameters
        sla_2011 = SLAx             ! initial
        sla_2020 = SLAx * (1-0.12)       ! end（+response ratio）
        k = abs(sla_2020-sla_2011)/10 ! update every year
        ! cal
        t1 = t - 2010
        calc_sla = sla_2011 + (sla_2020 - sla_2011) * (1.0 - EXP(-k * t1))
    END FUNCTION calc_sla

END MODULE update_traits

MODULE vars_site

    USE IntersVariables
    USE FileSize

    IMPLICIT NONE
    SAVE

    REAL, PARAMETER :: lat=22.68 !heshan !29.48 TianTong
    REAL longi
    REAL N_deposit,N_fert,P_fert
    REAL NSCmin,NSCmax   ! range of non-structural carbon pool
    REAL,DIMENSION(10):: thksl,FRLEN !fraction of root length in every layer
    REAL rdepth
    REAL LAIMAX,LAIMIN,SLA,SLAx,hmax 
    REAL, PARAMETER :: LAIMAX0 = 5.8!8.
    REAL, PARAMETER :: hl0  = 0.0008
    REAL alphaN,alphaP
    REAL LNC,LPC
    REAL gddonset
    REAL GRmax,GLmax,GSmax,GRemax
    REAL Q10,Rl0,Rs0,Rr0,Rre0         ! parameters for auto respiration
    REAL Rootmax,Stemmax,SapS,SapR,StemSap,RootSap
        ! parameters for photosynthesis model
    REAL stom_n,a1,Ds0,Vcmax0,extkU,xfang,alpha
    REAL phi
    REAL Vcmx0, eJmx0,eJmax0
    REAL WILTPT,FILDCP
    INTEGER ii
    REAL, parameter:: times_storage_use=720. ! for 30 days - rapid growth

    CONTAINS

    SUBROUTINE ReadInitialFile()

    !   Read the initial file for simulations

        IMPLICIT NONE
        INTEGER :: m, n, k, nvalid, day_idx, vcode
        REAL    :: sd, wt
        !INTEGER m,n
        CHARACTER(len=99) commts
        
        ! ! Site-level specific parameters
        ! OPEN(112,file='../input/pre-HeShanSpecificParameters.csv',status='old',ACTION='READ')
        ! ! Forcing data
        ! OPEN(111,file='../input/',status='old',ACTION='READ')
        ! ! Initial C pools, CN, CP ratios
        ! OPEN(113,file='../input/Initialstate_heshan.csv',status='old',ACTION='READ')
        
        ! IF (MCMC .eq. 1)THEN
        !     OPEN(114,file='../input/obs_for_MCMC_P0_2022_2024.txt',status='old',ACTION='READ')    !Observed data for data assimilation
        !     OPEN(115,file='../input/mcmc_range.csv',status='old',ACTION='READ')
        ! ENDIF
                ! Site-level parameters
        OPEN(112, file='../input/pre-HeShanSpecificParameters.csv', &
             status='old', action='read')

        ! Active forcing file, replaced by the shell script for each stage
        OPEN(111, file='../input/TECO forcing 2021-2024.txt', &
             status='old', action='read')

        ! Active initial state, replaced by the shell script for each stage
        OPEN(113, file='../input/Initialstate_heshan zero.csv', &
             status='old', action='read')

        IF (MCMC .EQ. 1) THEN
            OPEN(114, file='../input/obs_for_MCMC_P0_2022_2024.txt', &
                 status='old', action='read')
            OPEN(115, file='../input/mcmc_range.csv', &
                 status='old', action='read')
        ENDIF

        CALL AlloSize()  !Also need deallocate, TO BE IMPROVED 

    !   Read data from readed files

        WRITE(*,*) 'reading parameters ...'
        READ(112,'(a150)') commts
        write(*,*)commts
        READ(112,*) site_paras
        
        READ(111,'(a150)') commts
        write(*,*) commts

        m=0
        DO 
            m=m+1
            READ(111,*,IOSTAT=n) forcing_data(:,m)
            IF(n<0)exit 
        ENDDO !END READ the forcing data

        ! Check the forcing data
        write(*,*)'force_nhour,forcing_data',force_nhours,forcing_data(:,8760)
        WRITE(*,*)'END read forcing data, number of cols', num_cols

        m = 0
        Do 
            m = m+1
            READ(113,*,IOSTAT=n) InitialCNP(:,m)
            IF(n<0)exit 
        ENDDO !END READ the initial C,CN,CP
        
        WRITE(*,*)'QC',InitialCNP(:,1)
        close(111)
        close(112)
        close(113)

    !   Read observed data for MCMC
        IF (MCMC .eq. 1) THEN
            WRITE(*,*) 'reading obs for MCMC'
            READ(114,'(a160)') commts

             m = 0
            DO
                m = m + 1
                READ(114,*,IOSTAT=n) obs_TianTong(:,m)
                IF (n < 0) EXIT
            ENDDO
            CLOSE(114)

            nobs_assim = m - 1
            WRITE(*,*) 'obs_for_MCMC header: ', TRIM(commts)
            WRITE(*,*) 'nobs_assim(read) = ', nobs_assim

            ! ---------- 过滤 valid，并把 valid 压缩到前面 ----------
            nvalid = 0
            DO k = 1, nobs_assim
                day_idx = INT(obs_TianTong(1,k))
                vcode   = INT(obs_TianTong(2,k))
                sd      = obs_TianTong(4,k)

                ! weight<=0 自动设为 1
                IF (obs_TianTong(5,k) .LE. 0.0) obs_TianTong(5,k) = 1.0
                wt = obs_TianTong(5,k)

                IF (day_idx < 1 .OR. day_idx > output_ndays) CYCLE
                IF (sd <= 0.0) CYCLE
                IF (wt <= 0.0) CYCLE
                IF (.NOT. (                        &
                    vcode == 1   .OR.              &
                    vcode == 2   .OR.              &
                    vcode == 3   .OR.              &
                    vcode == 100 .OR.              &
                    vcode == 101 .OR.              &
                    vcode == 107 .OR.              &
                    vcode == 201 .OR.              &
                    vcode == 202 .OR.              &
                    vcode == 203 .OR.              &
                    vcode == 204)) CYCLE
                

                nvalid = nvalid + 1
                ! 压缩：把第 k 条 valid 放到第 nvalid 行
                obs_TianTong(:, nvalid) = obs_TianTong(:, k)
            ENDDO

            WRITE(*,*) 'nobs_assim(valid) = ', nvalid

            ! ---------- 打印前 5 条 valid ----------
            WRITE(*,*) 'First valid obs rows: day_idx, vcode, value, sd, weight'
            DO k = 1, MIN(5, nvalid)
                WRITE(*,'(I6,1X,I6,1X,F12.4,1X,F12.4,1X,F8.3)') &
                    INT(obs_TianTong(1,k)), INT(obs_TianTong(2,k)), &
                    obs_TianTong(3,k), obs_TianTong(4,k), obs_TianTong(5,k)
            ENDDO
            CALL FLUSH(6)

            ! 最关键：后续 MCMC 用 valid 数量
            nobs_assim = nvalid
        END IF


        ! IF (MCMC .eq. 1) THEN
        !     WRITE(*,*) 'reading obs for MCMC'
        !     READ(114,'(a160)') commts
        !     m=0
        !     DO 
        !         m=m+1
        !         READ(114,*,IOSTAT=n)obs_TianTong(:,m)
        !         IF(n<0)exit  
        !     ENDDO
        !     !!!!!!!!!!!!!
        !     INTEGER :: nread, k
        !     nread = m-1
        !     WRITE(*,*) 'obs_for_MCMC header: ', TRIM(commts)
        !     WRITE(*,*) 'nobs_assim(read) = ', nread

        !     ! 打印前 5 行（或不足 5 行就打印全部）
        !     WRITE(*,*) 'First obs rows: day_idx, vcode, value, sd, weight'
        !     DO k = 1, MIN(5, nread)
        !         WRITE(*,'(I6,1X,I6,1X,F12.4,1X,F12.4,1X,F8.3)') &
        !         INT(obs_TianTong(1,k)), INT(obs_TianTong(2,k)), &
        !         obs_TianTong(3,k), obs_TianTong(4,k), obs_TianTong(5,k)
        !     ENDDO
        !     CALL FLUSH(6)
        !     !!!!!!!!!!!!!
        !     close(114)
        ! ENDIF

    END SUBROUTINE ReadInitialFile

    SUBROUTINE site_value()

!Define variables before simulation, these variable DONT update every loop
!Fixed initial values
        !thickness of every soil layer
        thksl = [10.,10.,10.,10.,10.,10.,20.,20.,20.,20.]            !TianTong data_0515
        !ratio of roots in every layer
        FRLEN = [0.5,0.3,0.15,0.03,0.015,0.005,0.0,0.0,0.0,0.0]      !TianTong

        N_deposit=3.602/8760.  !(gN/m2/h)    
        !N_fixation = 0. !2.25/8760. !(gN/h/m2)  from value of Goll,2017
        N_fert=0.0 !5.6 ! (11.2 gN m-2 yr-1, in spring, Duke Forest FACE) 
        P_fert=0.0

!Paras read
        wsmax	   = 	35!site_paras(3)
        !#20% - 40% wan20241224
        wsmin	   = 	site_paras(4)!
        LAIMAX     = 	site_paras(5)!LAI最大值
        LAIMIN     = 	site_paras(6)!LAI最小值
        rdepth     = 	site_paras(7)!根系最大深度
        Rootmax    = 	site_paras(8)!最大根系边材量
        Stemmax    =  	site_paras(9)!最大茎边材量
        SapR	   = 	site_paras(10)!根边材比例
        SapS	   = 	site_paras(11)!茎边材比例
        SLAx	   = 	site_paras(12)!SLA cm2/g
        GLmax	   = 	site_paras(13)!叶最大生长潜力
        GRmax	   = 	site_paras(14)
        Gsmax	   = 	site_paras(15)
        GRemax     =    site_paras(58)
        stom_n     = 	site_paras(16)!气孔特性，单面/双面Vcmx0 = 18.28185
        a1	       = 	11!site_paras(17)
        Ds0	       = 	1000!site_paras(18)
        Vcmax0     = 	site_paras(19)
        extkU	   = 	site_paras(20)
        xfang	   = 	site_paras(21)!叶片角度分布参数
        alpha	   = 	site_paras(22)
        tau_L	   = 	site_paras(23)!叶周转时间/年
        tau_W	   =    site_paras(24)!茎周转时间/年
        tau_R	   = 	site_paras(25)!根周转时间/年
        tau_Re     =    site_paras(56)!繁殖周转时间/年
        tau_F	   = 	site_paras(26)!代谢周转时间/年
        tau_C	   = 	site_paras(27)!结构周转时间/年
        tau_Micr   = 	site_paras(28)!快周转时间/年
        tau_Slow   = 	site_paras(29)!慢周转时间/年
        tau_Pass   =    site_paras(30)!被动周转时间/年
        gddonset   = 	site_paras(31)
        Q10	       = 	site_paras(32) 
        Rl0	       = 	site_paras(33)
        Rs0	       = 	site_paras(34)
        Rr0	       = 	site_paras(35)
        Rre0       =    site_paras(59)
        phi        =    site_paras(36)
        hmax       =    site_paras(37)
        LPC        =    site_paras(38) 
        LNC        =    site_paras(39)
        alphaP     =    site_paras(40)
        alphaN     =    site_paras(41)
        etaL       =    site_paras(42)!叶凋落物进入代谢库的比例
        etaW       =    site_paras(43)
        etaR       =    site_paras(44)
        etaRe      =    site_paras(57)
        f_F2M      =    site_paras(45)!5 to 7
        f_C2M      =    site_paras(46)!6 to 7
        f_C2S      =    site_paras(47)!6 to 8
        f_M2S      =    site_paras(48)!7 to 8
        f_M2P      =    site_paras(49)!7 to 9
        f_S2M      =    site_paras(50)!8 to 7
        f_S2P      =    site_paras(51)!8 to 9
        f_P2M      =    site_paras(52)!9 to 7
        QPlab      =    0.021!site_paras(53) 
        QNminer    =    site_paras(54)!土壤无机氮库
        s_w_min    =    site_paras(55)
        rate_maxN  =    site_paras(60)
        rate_maxP  =    site_paras(61)

! UNIT TRANSFER
        !   the unit of residence time is transformed from yearly to hourly
        tauC=(/tau_L,tau_W,tau_R,tau_Re,tau_F,tau_C,tau_Micr,tau_Slow,tau_Pass/)*8760.  !hourly
        exitK = 1./tauC   !wan 20230708

        exitEC = 1.0/(3*8760)
        EC_store = 0.0

        !   growth rates of plant
        GLmax=GLmax/8760.
        GRmax=GRmax/8760.
        GSmax=GSmax/8760.
        GRemax=GRemax/8760.
        !   Initialize parameters and initial state:
        WILTPT=wsmin/100.0
        FILDCP=wsmax/100.0
        
        !   gddonset=320.0
        eJmax0 = 2.04*Vcmax0-10.49 ! Equation derive from measurements in TT. wan20241222
        Vcmx0 = Vcmax0*1.0e-6  !μmol to mol
        eJmx0 = eJmax0*1.0e-6 
        ! eJmx0 = 1.67*Vcmx0 ! Weng 02/21/2011 Medlyn et al. 2002
        !eJmx0 = 1.58 *Vcmx0    !wan Aug 28
        SLA=SLAx/10000.         ! Convert unit from cm2/g to m2/g

! INITIAL VALUES
! initialize variables
        !   pools and ratios
        QC  = InitialCNP(:,1)
        CN0 = InitialCNP(:,2)   
        CP0 = InitialCNP(:,3)

        newsoilNCmax(1) = 1.0/6.0!1.0/3.0 !5.95
        newsoilNCmax(2) = 1.0/10.0
        newsoilNCmax(3) = 1.0/10.0
            
        newsoilNCmin(1) = 1.0/10.0 !29.96
        newsoilNCmin(2) = 1.0/30.0
        newsoilNCmin(3) = 1.0/30.0

! test for new dataset 0902 

        P_deficit=0.0

        CN=CN0
        QN=QC/CN0
        CP=CP0
        QP=QC/CP0
        
        DO ii=1,10
            wcl(ii)=wsmax/100. !wsmin/100.  ! wsmax - site parameters,initial value  wsmax/100
        ENDDO 

        NSPmax = 0.2*(QP(1) + QP(2)+ QP(3) + QP(4))   !wan 2023/3/22 value = 2.16
        NSPmin = 0.001

        QPsorb = 1.268!133.0*QPlab/(64.0+QPlab)
        QPss= 18.09!54
        QPocc= 131.93!396.5
        P_reserve = 0.5*(QP(1) + QP(2)+ QP(3) + QP(4))

        ispinup = 1   
        onset = 0

        !   20230710
        !   let storage = accumulation, means initial the [storage] by [accumulation] 
        accumulation=0
        storage= 63.44!63.38!58.36!74.09           !g C/m2
        NSC= 9.64!10.22!18.88!25  !85.3          !wan- g C/m2  from 10 change to 23, 2023/3/22,initial value
        LAI=LAIMIN

        infilt=0.

        NSN = NSC/200. 
        NSP=NSC/2300. !change to 0.009 2023/3/21 !  change 1 to 0.5 2023/3/12

    END SUBROUTINE site_value

    SUBROUTINE loop_initial_value(iiterms,simu_year,EndDays,EndHours,stor_use)
        USE DaysHours
        IMPLICIT NONE
        INTEGER iiterms,i
        INTEGER simu_year,EndDays
        INTEGER simu_hours,EndHours
        TYPE(ResultType) :: result
        REAL stor_use

        idays=0               ! the whole days in simulated year (1 yr or more than 1 yr) in one spin-up loop
        iiterms= force_nhours ! the whole hours in simulated year (1 yr or more than 1 yr) in one spin-up loop
        iyear = 0
        simu_year = start_year+iyear
        result= Days_cal(simu_year)
        EndDays=result%Days
        EndHours=result%Hours

!   VARIABLES
        N_demand = 0.
        P_demand = 0.
        LFactorP = 1.
        LFactorN = 1.
        P_net = 0.
        NPP = 0.
        GPP = 0.
        Rnitrogen = 0.
        Rphosphorus = 0.
        Rauto = 0.
        ExcessC=0.

        !QPsorb=0.8
        !QPss= 9.4 ! reset the Qpss since we dont consider the occluded pools
        !QPocc=18.4

        !QPsorb=2.01!0.03!4!0.3
        !QPss=9.4!5!37.25!0.8!37.25!1.2
    !   define soil for export variables for satisfying usage of canopy submodel first time
        DO i=1,10
            wcl(i)=wsmax/100!wsmin/100.  ! wsmax - site parameters,initial value
        ENDDO 

        !fwsoil=0.5!1.0  !w! Parameter fw represents the relative availability of soil water for plants
        !topfws=0.25!1.0
        !omega=0.3!1.0
        !infilt=0.

        bmleaf=QC(1)/0.55
        bmstem=QC(2)/0.55
        bmroot=QC(3)/0.55
        bmRe  =QC(4)/0.55
        bmplant=bmstem+bmroot+bmleaf

        IF(ispinup .GT. 1)THEN
            storage = accumulation
            accumulation = 0.0
        ENDIF
        stor_use = storage/times_storage_use

        GDD5=0.0   !wan- initial value, debug
        onset = 0

        !NSN=6.0       !labile N pools, both for reserve and transfer station, wan, 2023/3/16
        
        N_deficit=0.
        
        OutN = 0.0
        N_immob=0.
        OutP =0.0
        P_immob=0.
        f_rootp = 0.02 !(this initial value shoule be smaller!)

        StemSap=AMIN1(Stemmax,SapS*bmStem)   
    ! Stemmax and SapS were input from parameter file, what are they? Unit? Maximum stem biomass? -JJJJJJJJJJJJJJJJJJJJJJ 
        RootSap=AMIN1(Rootmax,SapR*bmRoot)
        NSCmin=2. 
        NSCmax=0.05*(StemSap+RootSap+QC(1))  !wan- change 0.05 to 0.005 for Tiantong 20230322,change to 0.08

        !   calculating scaling factor of NSC
        IF(NSC.le.NSCmin)fnsc=0.0
        IF(NSC.ge.NSCmax)fnsc=1.0
        IF((NSC.lt.NSCmax).and.(NSC.gt.NSCmin))THEN 
            fnsc=(NSC-NSCmin)/(NSCmax-NSCmin)
        ENDIF

        
    END SUBROUTINE loop_initial_value


    SUBROUTINE step_forcing (doy,hour,Tair,Tsoil,RH,Dair,Rain,wind,&
                & PAR,radsol,co2,co2ca)

        IMPLICIT NONE

        REAL doy,hour
        REAL Tair,Tsoil,RH,Dair,Rain,wind
        REAL PAR,radsol,co2,co2ca

                !!! input the force_data	
        doy     =forcing_data(2,itime)
        hour    =forcing_data(3,itime)	  
        Tair    =forcing_data(4,itime)          ! Tair
        Tsoil   =forcing_data(5,itime)          ! SLT
        RH      =forcing_data(6,itime)
        Dair    =forcing_data(7,itime)          !air water vapour defficit? Unit Pa  !wan- VPD
        rain    =forcing_data(8,itime)          ! rain fal per hour          
        wind    =ABS(forcing_data(9,itime))     ! wind speed m s-1
        PAR     =forcing_data(10,itime)         ! Unit ? umol/s/m-2
        radsol  =forcing_data(10,itime)         ! unit ? PAR actually
        co2     =forcing_data(11,itime)
        co2ca=co2*1.0E-6 


    END SUBROUTINE

END MODULE vars_site


MODULE vars_consts

    IMPLICIT NONE
    SAVE
    REAL, PARAMETER :: pi = 3.1415926
!     physical constants
    REAL tauL(3), rhoL(3), rhoS(3)
    REAL, PARAMETER :: emleaf=0.96
    REAL, PARAMETER :: emsoil=0.94
    REAL, PARAMETER :: Rconst=8.314                 ! universal gas constant (J/mol)
    REAL, PARAMETER :: sigma=5.67e-8                ! Steffan Boltzman constant (W/m2/K4)
    REAL, PARAMETER :: cpair=1010.                  ! heat capacity of air (J/kg/K)   wan- need change to 1006?
    REAL, PARAMETER :: Patm=101325. !1.e5           ! atmospheric pressure  (Pa)
    REAL, PARAMETER :: Trefk=293.2                  !reference temp K for Kc, Ko, Rd
    REAL, PARAMETER :: H2OLv0=2.501e6               !latent heat H2O (J/kg)
    REAL, PARAMETER :: AirMa=29.e-3                 !mol mass air (kg/mol)
    REAL, PARAMETER :: H2OMw=18.e-3                 !mol mass H2O (kg/mol)
    REAL, PARAMETER :: chi=0.93                     !gbH/gbw
    REAL, PARAMETER :: Dheat=21.5e-6                !molecular diffusivity for heat
!     plant parameters
    REAL, PARAMETER :: gsw0 = 1.0e-2                !g0 for H2O in BWB model
    !REAL, PARAMETER :: eJmx0 = Vcmx0*1.67 ! 2.7            !@20C Leuning 1996 from Wullschleger (1993)    !wan 0523 change 2,7 to 1.67
    REAL, PARAMETER :: theta = 0.9
    REAL, PARAMETER :: wleaf=0.01                   !leaf width (m)

!     thermodynamic parameters for Kc and Ko (Leuning 1990)
    REAL, PARAMETER :: conKc0 = 302.e-6                !mol mol^-1
    REAL, PARAMETER :: conKo0 = 256.e-3                !mol mol^-1
    REAL, PARAMETER :: Ekc = 59430.                    !J mol^-1
    REAL, PARAMETER :: Eko = 36000.                    !J mol^-1
!     Erd = 53000.                    !J mol^-1
    REAL, PARAMETER :: o2ci= 210.e-3                   !mol mol^-1

!     thermodynamic parameters for Vcmax & Jmax (Eq 9, Harley et al, 1992; #1392)
    REAL, PARAMETER :: Eavm = 116300.               !J/mol  (activation energy)
    REAL, PARAMETER :: Edvm = 202900.               !J/mol  (deactivation energy)
    REAL, PARAMETER :: Eajm = 79500.                !J/mol  (activation energy) 
    REAL, PARAMETER :: Edjm = 201000.               !J/mol  (deactivation energy)
    REAL, PARAMETER :: Entrpy = 650.                !J/mol/K (entropy term, for Jmax & Vcmax)

!     parameters for temperature dependence of gamma* (revised from von Caemmerer et al 1993)
    REAL, PARAMETER :: gam0 = 28.0e-6               !mol mol^-1 @ 20C = 36.9 @ 25C
    REAL, PARAMETER :: gam1 = .0509
    REAL, PARAMETER :: gam2 = .0010
    REAL, PARAMETER :: Pa_air=101325.0         !Pa

    CONTAINS
    
        SUBROUTINE consts() 

            tauL(1)=0.1                  ! leaf transmittance for vis
            rhoL(1)=0.1                  ! leaf reflectance for vis
            rhoS(1)=0.1                  ! soil reflectance for vis
            tauL(2)=0.425                ! for NIR
            rhoL(2)=0.425                ! for NIR
            rhoS(2)=0.3                  ! for NIR - later function of soil water content
            tauL(3)=0.00                 ! for thermal
            rhoL(3)=0.00                 ! for thermal
            rhoS(3)=0.00                 ! for thermal
            !eJmx0 = Vcmx0*1.67           ! !@20C Leuning 1996 from Wullschleger (1993)    !wan 0523 change 2,7 to 1.67

        END SUBROUTINE consts

END MODULE vars_consts


!MODULE TransFraction
    !Module to decalre Translate fraction data to share between routines
!        USE IntersVariables
!        IMPLICIT NONE
!        SAVE
    !     partitioning coefficients
        !INTEGER,PARAMETER :: SGL= 4
!        REAL (SGL) :: etaL=0.6          ! 60% of foliage litter is fine, didn't use
!        REAL (SGL) :: etaW=0.15         ! 15% of woody litter is fine
!        REAL (SGL) :: etaR=0.85         ! 85% of root litter is fine  , didn't use    
!        REAL (SGL) :: f_F2M=0.55!0.55        ! *exp((CN0(4)-CN_fine)*0.1)
!        REAL (SGL) :: f_C2M=0.275!0.275       ! *exp((CN0(5)-CN_coarse)*0.1)
!        REAL (SGL) :: f_C2S=0.275!0.275       ! *exp((CN0(5)-CN_coarse)*0.1)
!        REAL (SGL) :: f_M2S=0.3!0.3
!        REAL (SGL) :: f_M2P=0.03!0.03!0.1  wan 0523 change to 0.01 
!        REAL (SGL) :: f_S2P=0.08!0.05!0.2        !0.03 Change by Jiang Jiang 10/10/2015   wan change to 0.03
!        REAL (SGL) :: f_S2M=0.5
!        REAL (SGL) :: f_P2M=0.45

!END MODULE TransFraction 