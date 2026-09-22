
!----------------------------------------------------
Module UptakeScalars

  USE IntersVariables
  IMPLICIT NONE
  CONTAINS
! Step 1, cal
  SUBROUTINE FDiffP_cal ()

  ! Apply the Fiki'law to calculate the diffusion flux between root surface and soil
  ! around root zone. F_diff = -D*delta_P_sol
  ! Variables description:
  ! FDiffP - the diffusion P flux of following Fick law [g m-2 h-1]

      IMPLICIT NONE
      REAL soilbulk,wsctot,swc_p,wcltot,waterD,swc_p_ref
      REAL sawp,eq1P,eq1
      REAL D_ref,rd,rr,f1,f2,pi,conv_D_ref
      REAL RLD,rdiff,tf,FdiffP,delta_QlabP_root,D
      REAL conv_fac_vmaxP!,!wcltot_rec(8760),wsctot_rec(8760)
      INTEGER i

    !molar mass for P P的相对原子量
      sawp = 31 ![(gP mol-1)]
    !relative water content 相对含水量参考值 扩散开始明显受限的含水阈值
      swc_p_ref = 0.12 !m3/m3 same as mm3/mm3
  ! Diffusion coefficient of phosphate in free water at 25degreeC  (Mollier et al 2008))
  !自由水中磷酸根扩散系数基准
      !D_ref = 1.581E-2  ![m2/30min]  (Table 2) 
      !D_ref = D_ref * 2 ![m2/h]
      D_ref = 0.759/24.*1E-4  !wan 2024 0417 m2/h 
  !parameters for diffusion path
  !Bonan 2014
      ! fine root radius, 0.1-0.5 mm 细根半径
      rr = 0.29E-3  !m 
      ! Specific root density 200-500 kg/m3, r^−1 = 12.2 m g−1 × πr2 根组织密度
      rd = 0.31E6   !g biomass /m3 root Bonan 2014, different with Goll 2017
  !Barraclough & Tinker 1981
      !The best fit line is given in Fig. 1 (fl = 1.58theta - 0.17).
      f1 = 1.58
      f2 = -0.17
  !constant pi
      pi=3.14159256

  !1 -- calculate the soil water 
      wcltot = 0.0
      DO i = 1,6  !I do not need ten layers
          !  relative soilwater content [m3/m3] for calculate root surface concentration
          wcltot = wcltot+wcl(i)   !! total volumetric soilwater content 6个土层的含水量求和的总含水量
      ENDDO
      swc_p = wcltot / 6  ! cm3/cm3 or m3/m3
      waterD = 1. !g/cm3 water density
      !thikness = 10 cm
      !*10 g/cm2 to kg/m2
      wsctot = swc_p*waterD*10*10 !unit kg/m2
      !wcltot_rec(itime) = wcltot
      !wsctot_rec(itime) = wsctot

  !2 -- calculate the root length density 单位体积土壤中根系总长度
      ! = root biomass density / (rd*pi*rr**2)
      ! root biomass density (g/m3) = root biomass * fraction / thickness (m)
      ! fraction 1: the soil contains f of the total root biomass
      ! root biomass: Bonan 2014 fine root biomass - 500 g m-2 -- here is QC(3)*coef
      ! thickness 0.1m
      ! since we have no thickness in soil cabron, we set the thickness is 1.5m
      RLD = ((QC(3)*2)/1.5)/(rd*pi*rr**2)
  !2 -- one-half the distance between roots (m)
      ! calculated with the assumption of uniform root spacing and 
      ! assuming the soil is divided into cylinders with the root along the middle axis
      !rdiff = MIN(0.1,(pi*RLD)**(-0.5))
      ! 0.1: restrict the diffusion path to 10 cm when there is dedcious tree
      rdiff = (pi*RLD)**(-0.5)
    
  !3 -- calculate the tortuosity factor

      IF (swc_p .ge. swc_p_ref) THEN
            tf = f1*swc_p+f2
      ELSE
            tf = (swc_p*(f1*swc_p+f2)/swc_p_ref)   
            !in orchidee code, only use this equation
            !some times, swc_p<swc_p_ref, the f_rootp equal 0! 0520
      ENDIF

  !4 -- calculate the diffusion coefficient
      !1.E-3 from Goll et al., 2017. 基准系数*平均相对含水量*曲折度因子*
      D = D_ref*swc_p*tf*(1/rdiff)!*1.E-3
      tf1_rec(itime) = f1*swc_p+f2
      tf2_rec(itime) = (swc_p*(f1*swc_p+f2)/swc_p_ref)
  !5 -- calculate the changes in the phosphorus concentration in the root zone
      ! Goll et al., 2017 Eq.23
      delta_QLabP_root = QPlab*(f_rootp-1)/swc_p
      ! in ORCHIDEE code:
      ! delta_QLabP_root = ((f_rootp-1)*QPlab)/(wsctot*1.E-3)
      ! 1.E-3 convert kg/m2 to g/m2
      ! weird: in orchidee tmc_pft(:,j)*1E-3, the unit of tmc_pft is kg/m2
  !5 -- calculate the diffusion flux
      FdiffP = -D * delta_QLabP_root
      FdiffP_rec(itime) = FdiffP
      delta_QLabP_root_rec(itime) = delta_QlabP_root
      D_rec(itime) = D
      swc_p_rec(itime) = swc_p
      tf_rec(itime) = tf
      rdiff_rec(itime) = rdiff
  !6 -- calculate the fraction of P concentration at root surface with F_diffP
      f_rootp = AMAX1(AMIN1((f_rootp*QPlab+FdiffP)/QPlab,1.0),0.0)
      !stop

  END SUBROUTINE FDiffP_cal

! Step 2 update the soil nutrient concentration of root surface using f_rootp
  !tmpvar = f_Pdissolved*soil_p_min(:,:,ipdissolved)

! Step 3 calulate the root maximal uptake capacity using the updated root surface 
  ! nutrient concentration
SUBROUTINE MaxUptakeCap_P (MaxUpCapP)

      IMPLICIT NONE
      REAL kpmin,kpmin_low,wcmax,sawp
      REAl conv_fac_concentP,MaxUpCapP

      !molar mass for P 
      sawp = 31 ![(gP mol-1)]
      !linear factor, chosen to match the observed rate of increase in overall phosphorus uptake
          ! at high dissolved labile phosphorus concentration μmolP L-1
      kPmin = 3 * sens_kpmin
      kPmin_low = 0.01 * sens_kpmin_low
      !maximal uptake capacity of roots  (Goll et al., 2017 Table 2),t = 30min
      !rate_maxP = 4.31E-6  ![g(P)/g(C)/t]
      !rate_maxP = rate_maxP * 2 ![g(P)/g(C)/h]
      !water holding capacity as approximation of pore space
      wcmax = 150 ![kg/m2]
      !Conversion factor from [umol(nutrient) per L] to [g(nutrient) m-2] 
      conv_fac_concentp=  1.e-6 * sawp  * 1.E3 * wcmax*1.E-3
  ! calculate the maximal uptake capacity
      MaxUpCapP = rate_maxP *rootSurfaceP  &
                  * (kPmin_low/conv_fac_concentp &
                      + (1/(rootSurfaceP+conv_fac_concentp*kPmin)))

  END SUBROUTINE MaxUptakeCap_P

  SUBROUTINE MaxUptakeCap_N (MaxUpCapN)
      IMPLICIT NONE
      REAL sawn,knmin,knmin_low
      REAL conv_fac_vmaxn,conv_fac_concentn,eq1N
      REAL wcmax,MaxUpCapN

      !rate_maxN = 5.4   !umol (g DryWeight_root)-1 h-1
      knmin_low = 0.0002
      knmin = 98
      wcmax = 150 ![kg/m2]
      
      sawn = 14 !g/mol
      conv_fac_vmaxn= 1.e-6 * sawn  
      conv_fac_concentn = sawn * 1.e3 * 1.e-6 * wcmax/1.e3   
      eq1N = knmin_low/conv_fac_concentn + (1/(QNminer+conv_fac_concentn*knmin))
      MaxUpCapN = rate_maxN*conv_fac_vmaxn*QNminer*eq1N

  END SUBROUTINE MaxUptakeCap_N

SUBROUTINE SCALAR_PN(f_PNplant)
      IMPLICIT NONE
      REAL pn_leaf_min,pn_leaf_max
      REAL n_plant,p_plant,pn_plant,f_PNplant
  
      !pn_leaf_max = 1/10.84 ! Goll 2017, Table A1
      !pn_leaf_min = 1/14.67

      !pn_leaf_max = 1/16.68 ! Goll 2017, Table A1
      !pn_leaf_min = 1/22.57

      pn_leaf_max = 1/20.87 ! wan
      pn_leaf_min = 1/42.99
  
      n_plant = QC(1)/CN(1)+QC(2)/CN(2)+QC(3)/CN(3)+QC(4)/CN(4)
      p_plant = QC(1)/CP(1)+QC(2)/CP(2)+QC(3)/CP(3)+QC(4)/CP(4)

      !n_plant = QC(1)/CN(1)+QC(3)/CN(3)+NSN
      !p_plant = QC(1)/CP(1)+QC(3)/CP(3)+NSP

      pn_plant = p_plant/n_plant
  
      f_PNplant = max((pn_plant - pn_leaf_max)/(pn_leaf_min - pn_leaf_max),0.01)
      f_PNplant = min(f_PNplant,1.)
  

END SUBROUTINE !SCALAR_PN


SUBROUTINE SCALAR_NC(f_NCplant)

      IMPLICIT NONE
      REAL nc_leaf_min,nc_leaf_max,c_plant,n_plant,nc_plant
      REAL f_NCplant

      !nc_leaf_max = 1/16.0
      !nc_leaf_min = 1/45.0

      nc_leaf_max = 1/26.19  !wan
      nc_leaf_min = 1/46.69
      
      c_plant = QC(1)+QC(2)+QC(3)+QC(4)
      n_plant = QC(1)/CN(1)+QC(2)/CN(2)+QC(3)/CN(3)+QC(4)/CN(4)

      !c_plant = QC(1)+QC(3)+NSC
      !n_plant = QC(1)/CN(1)+QC(3)/CN(3)+NSN

      nc_plant = n_plant/c_plant
      f_NCplant = max((nc_plant - nc_leaf_max)/(nc_leaf_min - nc_leaf_max),0.01)
      f_NCplant = min(f_NCplant,1.)

  
END SUBROUTINE SCALAR_NC

END MODULE UptakeScalars


MODULE NutrientsUptake
 
      !USE NutrientsDemand
  USE IntersVariables
  USE LIMITATION
  USE vars_site
  USE UptakeScalars

  IMPLICIT NONE
  CONTAINS

  SUBROUTINE N_UptakeDemand(NSNmax,NSNmin,Nfix0,&
                            N_fixation_x,Tsoil)

    !The plant nutrient uptake (g N m−2 s−1) from the soil mineral N pool 
    !is a function of the root biomass density
    !Ref: N uptake  Du2018, (10), (11); CABLE: Wang et al., 2010, (D12)
    
    IMPLICIT NONE
    REAL N_fixation_x
    REAL costCuptakeN,costCfix,costCreuse,NSNmax,NSNmin
    REAL Cfix,Tsoil,Cuptake,MaxUpCapN,f_NCplant
    REAL Qroot0,Nup0
    REAL ksye,Nfix0,Creuse0N1,Creuse0N2,Creuse0N3
    REAL a,b,FN_reserve,Nmax

    Nmax =12.0

    N_transfer=0.
    N_uptake=0.
    N_fixation_x =0.
    costCuptakeN=0.
    costCfix=0.
    costCreuse=0.

    ! Carbon cost for N resorption
    Creuse0N1 = 1/(QC(1)/CN(1))      !wan - [Fisher, 2010, FUN-model], 1 gC m-2, leaf. original Cresuse0N =2
    Creuse0N2 = 1/(QC(2)/CN(2))      !wood
    Creuse0N3 = 1/(QC(3)/CN(3))      !root
    ! Carbon cost for fixation and uptake
    Cfix = -6.25*(exp(-3.62+0.27*Tsoil*(1-0.5*(Tsoil/25.15)))-2)   !Cfix bwtween 7.5 ~ 12.5
    Cuptake = (1/QNminer)*(1/QC(3))

    ! Step 1 -- calculate N demand according to the new growth of tissues 
    ! and N_deficit from last loop
    N_demand = NPP_L/CN(1)+NPP_W/CN(2)+NPP_R/CN(3)+NPP_Re/CN(4) + N_deficit
    
    ! Step 2 -- calculate the resorted N with fixed ratio of leaf, wood and root
    !IF (NPP .gt. 0) THEN  ! NPP, available C for nutrient activity
    N_transfer =(OutN(1) + OutN(2) +OutN(3))*alphaN 
    ! Step 3 -- change Creuse from fixed value to tissues [N] depended vlaue.
    costCreuse =(Creuse0N1*OutN(1)+Creuse0N2*OutN(2)+Creuse0N3*OutN(3))*alphaN     
    !ENDIF
    N_demand   = AMAX1(N_demand-N_transfer,0.0)
    
    ! step 4 -- root uptake N from soil if there exist demand
    ! select for uptake or fixation according to the C cost
    If(N_demand > 0.0)THEN   !add by wan, for debug
      !     N uptake  Du2018, (10), (11)
      !     wan: ksye/QNminer is a standardized processing,
      !     When QNminer is big enough, the uptake occur, otherwise, plant can't uptake nutrients from soil.
            !         N_uptake=AMIN1(N_demand+N_deficit,      &
            !&                       QNminer*QC(3)/(QC(3)+Qroot0), &
            !&                       Nup0*NSC/(ksye/QNminer))  ！ Du et al., 2018
      !     ksye - C cost per N for uptake, Cfix0 - C cost per N for fixation
      !IF(Cuptake .lt. Cfix)THEN  现在默认开启N的固定
          CALL SCALAR_NC(f_NCplant)
          CALL MaxUptakeCap_N(MaxUpCapN)
          N_uptake = AMIN1 (MAX(0.0, N_demand), &
          &                 MaxUpCapN*(QC(3)/fbmC)*S_t(1)*f_NCplant)

          aN = MAX(0.0, N_demand+N_deficit)
          bN = MaxUpCapN*(QC(3)/fbmC)*S_t(1)*f_NCplant

          ! calculate the limitation factor - LfactorN
          CALL xnp_soil()
        ! ensure soil P min doesn't get depleted
          !IF(N_uptake .gt. N_miner2)THEN
          !      N_uptake=N_miner2  ! From CABLE!!!! WAN   no change the results
          !      LFactorN = N_miner2/aN
          !ENDIF
          costCuptakeN=N_uptake*Cuptake
          N_deficit=N_demand - N_uptake
      !ELSE
          ! Nitrogen fixation
          !N_fixation_x=Amin1(AMAX1(0.0,N_demand),fnsc*Nfix0*NSC)  !wan: fnsc - nsc-limiting factor
          !costCfix=Cfix*N_fixation_x
          !N_demand=N_demand-N_fixation_x
          !aN = AMAX1(0.0, N_demand+N_deficit)
          !bN = fnsc*Nfix0*NSC
          !CALL xnp_soil()
          !N_deficit=N_demand - N_fixation_x
      !ENDIF

    ELSE  !when NPP is zero at night, the P_demand is 0, so no need for uptaking P
      N_demand = 0.0
      N_uptake = 0.0
      costCuptakeN = 0.0    
      !costCfix = 0.0
      N_deficit = N_demand
      LfactorN =1.0
    ENDIF ! the uptake of fixation

    ! wan add this part to avoid no fixation 
    ! the C-fix always bigger than C-uptake according to the equation
    !IF ((N_deficit .gt. 0).and. (Cuptake .lt. Cfix))THEN
        N_demand = N_deficit
        !N_fixation_x=Amin1(N_demand,fnsc*Nfix0*NSC)  !wan: fnsc - nsc-limiting factor
        N_fixation_x=fnsc*Nfix0*NSC*((Nmax-QNminer)/Nmax)
        costCfix=Cfix*N_fixation_x
        !N_deficit=N_demand-N_fixation_x
        !aN=N_demand
        !bN=fnsc*Nfix0*NSC
        !CALL xnp_soil
    !ENDIF

    ! step 5 -- the N allocation according to new growth and the initial CN ratio

    N_leaf = NPP_L/CN(1)+QC(1)/CN0(1)-QC(1)/CN(1)
    N_wood = NPP_W /CN(2)+QC(2)/CN0(2)-QC(2)/CN(2)
    N_root = NPP_R/CN(3)+QC(3)/CN0(3)-QC(3)/CN(3)
    N_re =   NPP_Re/CN(4)+QC(4)/CN0(4)-QC(4)/CN(4)
    
    ! step 6 -- update the NSN
    NSN=NSN+N_transfer+N_uptake-(N_leaf+N_wood+N_root+N_re) !+N_fixation_x

    IF (NSN .LE. (0.6 * NSNmax))THEN
      FN_reserve = AMIN1((NSNmax-NSN),0.1*N_reserve)
      NSN = NSN + FN_reserve
      N_reserve = N_reserve - FN_reserve
    ELSE
      FN_reserve = NSN - NSNmax
      NSN = NSN - FN_reserve
      N_reserve = N_reserve + FN_reserve
    ENDIF

  !     Total C cost for nitrogen   
    Rnitrogen=costCuptakeN+costCfix+costCreuse
    N_deficit_rec(itime) = N_deficit
    N_demand_rec(itime) = N_demand
    N_uptake_rec(itime) = N_uptake
    N_fixation_x_rec(itime) = N_fixation_x


  END SUBROUTINE N_UptakeDemand
      !---------------------------------------------------
      ! Plant demand-uptake P
      
  SUBROUTINE P_UptakeDemand() 
  
  ! Calculte the Plant uptake P, and update the NSP
  ! same flow with N cycle

    IMPLICIT NONE
    REAL Creuse0P !alphaP - trasfer fraction befor littering
    REAL LDOP0
    REAL costCuptakeP,costCreuseP
    REAL Creuse0P1,Creuse0P2,Creuse0P3,CPuptake
    REAL MaxUpCapP,f_PNplant,a,b
    REAL FP_reserve
    
    P_transfer=0.
    P_uptake=0.
    costCuptakeP=0. 
    costCreuseP=0.
    Creuse0P=1 !gC/gP
    !P_fixation=0.      ! no P fixation,

    Creuse0P1 = 1/(QC(1)/CP(1))      !wan - [Fisher, 2010, FUN-model], 1 gC m-2, leaf. original Cresuse0N =2
    Creuse0P2 = 1/(QC(2)/CP(2))      !wood
    Creuse0P3 = 1/(QC(3)/CP(3))      !root

    CPuptake = (1/QPlab)*(1/QC(3))
    
    ! Step 1 -- P demand
    P_demand=NPP_L/CP(1)+NPP_W/CP(2)+NPP_R/CP(3)+NPP_Re/CP(4)+P_deficit
    P_demand_test=NPP_L/CP(1)+NPP_W/CP(2)+NPP_R/CP(3)+NPP_Re/CP(4)

    ! Step 2 -- P resorption
    !IF (NPP .gt. 0) THEN
    P_transfer=(OutP(1) + OutP(2) +OutP(3))*alphaP
    ! Step 3 C cost for resorption
    costCreuseP = (Creuse0P1*OutP(1)+Creuse0P2*OutP(2)+Creuse0P3*OutP(3))*alphaP
    !ENDIF
    P_demand=AMAX1 (P_demand-P_transfer, 0.0)
    !P_demand_test = AMAX1 (P_demand_test-P_transfer, 0.0)
    
    ! Step 4 -- Uptake
    ! When the resorption < demand, the plant need uptake P from soil additionally
    ! Goll et al., 2017: the P concertation on root surface
    IF(P_demand > 0.0)THEN   
        ! (1) calculate the fraction of P concentration at root surface
        ! using Fiki's law
        CALL FDiffP_cal()
        ! f_rootp -- fraction of P concentration at root surface
        ! (2) update the P concentration on the root surface
        f_rootp_rec(itime) = f_rootp
        rootSurfaceP = f_rootp * QPlab
        ! (3) caculate the maximal root uptake capacity
        CALL MaxUptakeCap_P(MaxUpCapP)
        ! (4) calculate the N scalar on the P uptake
        CALL SCALAR_PN(f_PNplant)
        
        ! (5) calculate the actual uptake with the maximal uptake rate and root biomass
        ! The dynamic of P uptake has big impact on NPP, 
        ! thus a feasible MaxUpCapP will lead a increased loops of spin-up to reach the steady state
        ! MaxUpCapP = 1.E-5 ～ 8.E-6 is a safe rate. too low will lead a too low NPP
        ! MaxUpCapP = 8.E-6
        ! QPlab = 0.3
        !f_PNplant = 0.02
        P_uptake = AMIN1 (MAX(0.0, P_demand),MaxUpCapP*(QC(3)/fbmC)*S_t(1)*f_PNplant)
        aP = MAX(0.0, P_demand+P_deficit)
        bP = MaxUpCapP*(QC(3)/fbmC)*S_t(1)*f_PNplant
        MaxUpCapP_rec(itime) = MaxUpCapP
        QC3_rec(itime) = QC(3)
        aP_rec(itime) = aP
        bP_rec(itime) = bP
        ! calculate the limitation factor - LfactorP
        CALL xnp_soil()

        ! ensure soil P min doesn't get depleted

        P_demand_rec(itime) = P_demand
        P_demand_rec_test(itime) = P_demand_test
        P_deficit_rec(itime) = P_deficit
        P_uptake_rec(itime) = P_uptake
        P_miner2_rec(itime) = P_miner2
        LfactorP_rec(itime) = LfactorP
        
        ! (7) update the reduction fraction of [P] again
        f_rootp = (f_rootp*QPlab-P_uptake)/QPlab

        ! (8) C cost for uptake
        costCuptakeP=P_uptake*CPuptake !ksye/QPlab
          
    ELSE  !when NPP is zero at night, the P_demand is 0, so no need for uptaking P
        P_demand = 0.0
        P_uptake = 0.0
        costCuptakeP =0.0
        P_deficit = P_demand
        LfactorP = 1.0
    ENDIF

    P_deficit=Amax1(0.0, P_demand-P_uptake)
    Rphosphorus=costCuptakeP+costCreuseP
    
    !Nutrients allocation
    ! step 5 -- the N allocation according to new growth and the initial CN ratio

    P_leaf = NPP_L/CP(1)+QC(1)/CP0(1)-QC(1)/CP(1)
    P_wood = NPP_W/CP(2)+QC(2)/CP0(2)-QC(2)/CP(2)
    P_root = NPP_R/CP(3)+QC(3)/CP0(3)-QC(3)/CP(3)
    P_Re =   NPP_Re/CP(4)+QC(4)/CP0(4)-QC(4)/CP(4)

    ! step 6 -- update the NSP
    NSP=NSP+(P_transfer+P_uptake)-(P_leaf+P_wood+P_root+P_re)

    ! A reserve pool server the long-term variation
    IF (NSP .LE. (0.6 * NSPmax))THEN
          FP_reserve = AMIN1((NSPmax-NSP),0.1*P_reserve)
          NSP = NSP + FP_reserve
          P_reserve = P_reserve - FP_reserve
    ELSE
          FP_reserve = NSP - NSPmax
          NSP = NSP - FP_reserve
          P_reserve = P_reserve + FP_reserve
    ENDIF
          

  END SUBROUTINE P_UptakeDemand

END MODULE NutrientsUptake
      
      
      !------------------------------------------------
      
