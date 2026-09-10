MODULE SPINUP_mod
    IMPLICIT NONE
    CONTAINS

    ! This module will be call if ispinup > 2

    SUBROUTINE Spinup_output(m)
    !This subtoutine for the output of traditional spinup
    USE FileSize
    IMPLICIT NONE
    INTEGER m,i

    !variables: NPP and passive SOM for criteria of SASU
    !output_ndays is cumulative variable

    LoopNPP_d(m)=sum(output_daily(7,:))/output_ndays    ! [/output_ndays], ref: see xia et al., 2012, Fig.3 the unit is gC m-2 d-1
    LoopNPP_yr(m)=sum(output_daily(7,:))                ! annual NPP, add for test wrong: output_daily not for one year
    LoopNPP_yr_noadd(m)=sum(output_daily(25,:))
    LoopNEE_d(m)=sum(output_daily(3,:))/output_ndays
    LoopNEE_yr(m)=sum(output_daily(3,:))
    LoopQC9(m) = outputd_ccycle_Cpools(9,output_ndays)   !the slowest pool
    LoopQC(:,m)= outputd_ccycle_Cpools(:,output_ndays)   !the m row = the [output_ndays] row C pools
    
    ! allocate the value for temporary space [Loop_variables]
    Loop_variables(1,m)  = LoopNPP_d(m)
    Loop_variables(2,m)  = LoopNPP_yr(m)
    Loop_variables(3,m)  = LoopNEE_d(m)
    Loop_variables(4,m)  = LoopNEE_yr(m)
    Loop_variables(5:13,m)= LoopQC(:,m)
    Loop_variables(14,m)  = LoopNPP_yr_noadd(m)

    IF (CYCLE_CNP .gt. 1)THEN
        LoopQN(1:9,m)= outputd_ncycle_Npools(:,output_ndays)
        Loop_variablesNNN(:,m)=LoopQN(:,m)
    ENDIF
    IF (CYCLE_CNP .gt. 2)THEN
        LoopQP(1:9,m) = outputd_Pcycle_Ppools(:,output_ndays) 
        LoopQP(10,m)  = sum(outputd_Pdynamic(1,:))  !Dsoillab
        LoopQP(11,m)  = sum(outputd_Pdynamic(2,:))
        LoopQP(12,m)  = sum(outputd_Pdynamic(3,:))
        LoopQP(13,m)  = sum(outputd_Pdynamic(4,:))
        LoopQP(14,m)  = sum(outputd_Pdynamic(5,:))
        LoopQP(15,m)  = sum(outputd_Pdynamic(6,:))
        LoopQP(16,m)  = sum(outputd_Pdynamic(7,:))
        LoopQP(17,m)  = sum(outputd_Pdynamic(8,:))

        Loop_variablesPPP(:,m)=LoopQP(:,m)
    ENDIF

END SUBROUTINE Spinup_output

FUNCTION AssignMatrixD1(sourceD1,nloops) 
    IMPLICIT NONE

    REAL, DIMENSION(:) :: sourceD1     ! 源矩阵
    REAL, DIMENSION(nloops) :: AssignMatrixD1
    INTEGER :: i,n,ndim,nloops

    ndim = SIZE (SHAPE(sourceD1))
    n = SIZE (sourceD1)
    DO i = 1,n
        AssignMatrixD1 (i) = sourceD1 (i)
        IF((nloops - n) .EQ. 1) AssignMatrixD1 (nloops) = -999.
    ENDDO
END FUNCTION AssignMatrixD1

FUNCTION AssignMatrixD2(sourceD2,nloops) 
    IMPLICIT NONE
    REAL, DIMENSION(:,:) :: sourceD2   
    REAL, DIMENSION(SIZE(sourceD2,1),nloops) :: AssignMatrixD2  ! 结果矩阵
    INTEGER :: i, j,m,n,ndim,nloops

    ndim = SIZE (SHAPE(sourceD2))
    m = SIZE(sourceD2,1)
    n = SIZE(sourceD2,2)
    !write(*,*)'m,n,nloops',m,n,nloops

    DO i = 1,n     !row
        DO j = 1,m !col
            AssignMatrixD2 (j,i) = sourceD2 (j,i)
            IF((nloops - n) .EQ. 1) AssignMatrixD2 (j,nloops) = -999.
        ENDDO
    ENDDO
END FUNCTION AssignMatrixD2

END MODULE








