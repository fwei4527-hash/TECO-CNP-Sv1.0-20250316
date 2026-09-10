
Module LIMITATION

    IMPLICIT NONE
    CONTAINS
    !CASA_CNP
    SUBROUTINE xnp(alphaN,alphaP)
    
    USE IntersVariables
    IMPLICIT NONE
    REAL xncleaf,xncstem,xncroot,xnlimit
    REAL xpcleaf,xpcstem,xpcroot,xplimit
    REAL alphaN,alphaP
    REAL xNuptake,xPuptake
    REAL Nreqmax(3), Nreqmin(3),Ntrans(3),totNreqmax,totNreqmin
    REAL xnCnpp,xpCnpp
    REAL Preqmax(3), Preqmin(3),Ptrans(3),totPreqmax,totPreqmin
    REAL xn_stem_limit,xp_stem_limit
    REAL xn_root_limit,xp_root_limit
    REAL xnpmax
    REAL tot_N_demand,tot_P_demand

    xnpmax = 1.28 !TECO !1.28 !CABLE

    NPP_L=NPP*alpha_L
    NPP_W=NPP*alpha_W
    NPP_R=NPP*alpha_R
    !write(*,*)'NPP',NPP

    ! nutrient concentration limiting factor
    xnlimit  = 1.0 ! CN limitation
    xplimit  = 1.0 ! CP limitation
    xnp_leaf_limit=1.0 !This factor calculated from CABLE
    xnp_stem_limit=1.0
    xnp_root_limit=1.0

    IF(CYCLE_CNP == 1.0) THEN
        xnp_leaf_limit=1.0
    ENDIF

    IF(CYCLE_CNP == 2.0) THEN

        xncleaf = QN(1)/(QC(1)+1.0e-10)  
        xnlimit = xncleaf/(xncleaf+0.01)
        !0.01,kn, an empirical parameter for nitrogen limitation on NPP
        ! = 0.01gN/gC
        xplimit=1.0
        !write(*,*)'xnlimit_leaf,QN(1),QC(1)',xnlimit*xnpmax,QN(1),QC(1)
        xnp_leaf_limit=AMIN1(xnlimit,xplimit)*xnpmax 

        !---for stem
        xncstem = QN(2)/(QC(2)+1.0e-10)
        xn_stem_limit = xncstem/(xncstem+0.01)
        xp_stem_limit = 1.0
        xnp_stem_limit=AMIN1(xn_stem_limit,xp_stem_limit)*xnpmax

        !---for root
        xncroot = QN(3)/(QC(3)+1.0e-10)
        xn_root_limit = xncroot/(xncroot+0.01)
        xp_root_limit = 1.0
        xnp_root_limit=AMIN1(xn_root_limit,xp_root_limit)*xnpmax
    ENDIF

    IF(CYCLE_CNP == 3.0) THEN

        xncleaf = QN(1)/(QC(1)+1.0e-10)  
        xnlimit=xncleaf/(xncleaf+0.01)
        xpcleaf=QP(1)/(QC(1)+1.0e-10)
        xplimit=xpcleaf/(xpcleaf+0.0006)
        !0.0006 kp
        xnp_leaf_limit=AMIN1(xnlimit,xplimit)*xnpmax ! why *xnpmax in casa

        !write(*,*)'QN(1),QC(1),xncleaf,xnlimit',QN(1),QC(1),xncleaf,xnlimit
        !write(*,*)'QP(1),QC(1),xpcleaf,xplimit',QP(1),QC(1),xpcleaf,xplimit

        !---for stem
        xncstem = QN(2)/(QC(2)+1.0e-10)
        xn_stem_limit = xncstem/(xncstem+0.01)  !0.01
        !xn_stem_limit = xncstem/(xncstem+0.001)
        xpcstem = QP(2)/(QC(2)+1.0e-10)
        xp_stem_limit = xpcstem/(xpcstem+0.0006) !0.0006
        !xp_stem_limit = xpcstem/(xpcstem+0.0002) !0.0006

        !write(*,*)xpcstem,xncstem,xn_stem_limit*xnpmax,xp_stem_limit*xnpmax

        xnp_stem_limit=AMIN1(xn_stem_limit,xp_stem_limit)*xnpmax
        !write(*,*)'QP(2),QC(2),xp_stem_limit',QP(2),QC(2),xp_stem_limit
        

        !---for root
        xncroot = QN(3)/(QC(3)+1.0e-10)
        xn_root_limit = xncroot/(xncroot+0.01)
        xpcroot = QP(3)/(QC(3)+1.0e-10)
        xp_root_limit = xpcroot/(xpcroot+0.0006)

        xnp_root_limit=AMIN1(xn_root_limit,xp_root_limit)*xnpmax
        !write(*,*)'xnlimit,xn_stem_limit,xn_root_limit',xnlimit,xn_stem_limit,xn_root_limit
        !write(*,*)'xplimit,xp_stem_limit,xp_root_limit',xplimit,xp_stem_limit,xp_root_limit

    ENDIF

    END SUBROUTINE xnp
    
    SUBROUTINE xnp_soil()
    
        USE IntersVariables
        IMPLICIT NONE
        !REAL aN, aP !actual_demand
        !REAL bN, bP !pre_uptake
        IF(CYCLE_CNP .gt. 1)THEN
            IF(aN>bN)THEN
                !LfactorN = bN/aN
                LfactorN = 1/(1 +exp(-12*(bN/aN)+6))
            ELSE
                LFactorN = 1.
            ENDIF
            LfactorP = 1.
        ENDIF

        IF(CYCLE_CNP .gt.2)THEN
            IF (aP>bP)THEN
                !LfactorP = bP/aP
                LfactorP = 1/(1 +exp(-12*(bP/aP)+6)) ! wan20241223
            ELSE
                LfactorP = 1.
            ENDIF
        ENDIF
 
    END SUBROUTINE xnp_soil

END MODULE LIMITATION



