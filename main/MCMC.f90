MODULE FOR_MCMC

IMPLICIT NONE
CONTAINS

SUBROUTINE MCMC_simu (upgraded)

      USE FileSize
      USE IntersVariables
      IMPLICIT NONE
      INTEGER IDUM,upgraded,isimu
      REAL search_length
      REAL J_last
      INTEGER,PARAMETER :: npara=16  !para for mcmc
      REAL coefminmax(npara,2)
      REAL coef(npara),coefmin(npara),coefmax(npara) !set paras for MCMC
      REAL coefac(npara),coefnorm(npara),a(npara)
      REAL gamma(npara,npara),gamnew(npara,npara)     ! covariance matrix
      INTEGER k1,k2,rejet,paraflag,k3,m
      INTEGER, PARAMETER :: nc=100
      INTEGER, PARAMETER :: ncov=1000
	
      REAL coefhistory(ncov,npara)
      REAL covexist   
      REAL r,fact_rejet
      INTEGER i,j,new,reject,tmp_up,n,flag_new
      character(len=512) covfile,outfile,outfile1,outfile2
 
	!coef=site_paras !now, we choose all paras to MCMC
      coef(1)=site_paras(32) !Q10=2.105
      !coef(2)=site_paras(12)!SLA=150!!!!
      coef(2)=site_paras(19)!Vcmax0=33.721!!!!
      coef(3)=site_paras(23)!T1                              
      coef(4)=site_paras(24)!T2
      coef(5)=site_paras(25)!T3  
      coef(6)=site_paras(56)!T4     
      coef(7)=site_paras(26)!T5    
      coef(8)=site_paras(27)!T6   
      coef(9)=site_paras(28)!T7    
      coef(10)=site_paras(29)!T8  
      coef(11)=site_paras(30)!T9 
      coef(12)=2
      coef(13)=0.5
      coef(14)=2.09
      coef(15)=0.9
      coef(16)=0.78
      !coef(13)=site_paras(42)!etaL
      !coef(14)=site_paras(43)!etaW
      !coef(15)=site_paras(44)!etaR
      !coef(16)=site_paras(57)!etaRe
      !coef(13)=site_paras(45)!f_F2M
      !coef(14)=site_paras(46)!f_C2M
      !coef(15)=site_paras(47)!f_C2S
      !!coef(16)=site_paras(48)!f_M2S
      !coef(17)=site_paras(49)!f_M2P
      !coef(18)=site_paras(50)!f_S2M
      !coef(19)=site_paras(51)!f_S2P
      !coef(20)=site_paras(52)!f_P2M
      READ(115,*,IOSTAT=n) coefminmax(:,1)
        IF (n /= 0) STOP "ERROR: cannot read coefmin line from mcmc_range.csv"

      READ(115,*,IOSTAT=n) coefminmax(:,2)
        IF (n /= 0) STOP "ERROR: cannot read coefmax line from mcmc_range.csv"

    !   m = 0
    ! DO 
    !     m=m+1
    !     READ(115,*,IOSTAT=n) coefminmax(:,m)
    !     IF(n<0)exit 
    ! ENDDO
    write(*,*)'coefminmax',coefminmax
    coefmin = coefminmax(:,1)
    coefmax = coefminmax(:,2)
    WRITE(*,*)'coefmin',coefmin
    WRITE(*,*)'coefmax',coefmax
    CLOSE(115)      
    !stop
	!initialize covariance matrix
	covexist=0
	IF (covexist.eq.1)THEN ! If prior covariance exists, read from file
	    covfile='covariance.txt'
	    CALL getCov(gamma,covfile,npara)
        CALL racine_mat(gamma,gamnew,npara)      ! square root of covariance matrix
        gamma=gamnew
        DO k1=1,npara
!           coefnorm(k1)=(coef(k1)-coefmin(k1))/(coefmax(k1)-coefmin(k1))
            coefnorm(k1)=0.5
            coefac(k1)=coefnorm(k1)
        ENDDO
    ELSE
        coefac=coef
    ENDIF
	
    fact_rejet=0.2 !2.4/sqrt(REAL(npara))
    search_length=0.02 !0.05  !!0.05 ! wan change this value from 0.05 to 0.2 (08 Sep.)
    rejet = 0
	
	J_last=9000000.0
	IDUM = 542    ! IDUM is a random number to generate different random number sequence
	upgraded=0
	new=0
    flag_new = 0
    k3=0
	!OPEN(71,file='output/Paraest.csv')
    OPEN(71,file=TRIM(ADJUSTL(fileplace))//'Paraest'//filesignal)
    OPEN(73,file=TRIM(ADJUSTL(fileplace))//'Paraest_cov'//filesignal)

      !MCCM NOW
    isimu = 0
	DO !isimu=1,50000
        isimu = isimu+1
        IF(covexist.eq.1)THEN   
            paraflag=1
            DO while(paraflag.gt.0)
                CALL gengaussvect(fact_rejet*gamma,coefac,coefnorm,npara)
                paraflag=0
                DO k1=1,npara
                    IF(coefnorm(k1).lt.0. .or. coefnorm(k1).gt.1.)THEN
                    paraflag=paraflag+1
                    WRITE(*,*)'out of range,k1',paraflag,k1
                    ENDIF
                ENDDO
            ENDDO
            DO k1=1,npara
                coef(k1)=coefmin(k1)+coefnorm(k1)*(coefmax(k1)-coefmin(k1))
            ENDDO
            !stop
        ELSE
            !write(*,*)'call coefgenerate'
            CALL coefgenerate(isimu,coefac,coefmax,coefmin,coef,search_length,npara)       
        ENDIF
!!!END of regenerating paras
	
!updated paras
        site_paras(32)=coef(1) !Q10=2.105
        !site_paras(12)=coef(2)!SLA=150!!!!
        site_paras(19)=coef(2)!Vcmax0=33.721!!!!
        site_paras(23)=coef(3)!T1                              
        site_paras(24)=coef(4)!T2
        site_paras(25)=coef(5)!T3    
        site_paras(56)=coef(6)!T4   
        site_paras(26)=coef(7)!T5    
        site_paras(27)=coef(8)!T6   
        site_paras(28)=coef(9)!T7    
        site_paras(29)=coef(10)!T8  
        site_paras(30)=coef(11)!T9
        s_Psorbmax      = coef(12)
        s_kmlabP        = coef(13)
        s_sorb_to_ssorb = coef(14)   
        s_ssorb_to_occ  = coef(15)
        s_P_loss        = coef(16)   
        !site_paras(42)=coef(13)!etaL
        !site_paras(43)=coef(14)!etaW
        !site_paras(44)=coef(15)!etaR
        !site_paras(57)=coef(16)!etaRe
        !site_paras(45)=coef(13)!f_F2M
        !site_paras(46)=coef(14)!f_C2M
        !site_paras(47)=coef(15)!f_C2S
        !site_paras(48)=coef(16)!f_M2S
        !site_paras(49)=coef(17)!f_M2P
        !site_paras(50)=coef(18)!f_S2M
        !site_paras(51)=coef(19)!f_S2P
        !site_paras(52)=coef(20)!f_P2M

! start running the CNP_simu
        WRITE(*,*) 'P parameters before CNP_simu:', &
    s_sorb_to_ssorb, s_ssorb_to_occ
   CALL CNP_simu()
! M-H algorithm
	tmp_up=upgraded  !wan: temporary upgraded
    WRITE(*,*) 'BEFORE COST day366 pools:', &
    outputd_ccycle_Cpools(1,366),       &
    outputd_ccycle_Cpools(2,366),       &
    outputd_ccycle_Cpools(3,366)
	!calculate the cost
	CALL costFObsNee (J_last,upgraded,isimu)
    !(output_daily_mcmc,obs_TianTong,nobs,J_last,upgraded)
!+++********insert parts that was not included
    IF(upgraded.gt.tmp_up)THEN !IF upgraged, THEN...
        new=new+1 ! add a accepted value
        IF(covexist.eq.1)THEN  !covexist 
            coefac=coefnorm
            coefhistory(new,:)=coefnorm
        ELSE ! no covexist THEN ...
            coefac=coef
            DO k1=1,npara
                coefnorm(k1)=(coef(k1)-coefmin(k1))/(coefmax(k1)-coefmin(k1)) !norm the accepted values
            ENDDO
        ENDIF
        coefhistory(new,:)=coefnorm !add the normed accepted values to the history file

        IF(new.ge.ncov)THEN
            new=0 !! ncov is the covariance matrix gamma will be updated  !!
                    !! every ncov iterations 1000
            flag_new = 1 !wan add this to generate new para file
        ENDIF
        WRITE(71,701) upgraded,(coef(i),i=1,npara)
        IF(flag_new .eq. 1 .and. new  .gt. 0) WRITE(73,701) new,(coef(i),i=1,npara)

    701 format(I12,",",35(F15.4,",")) 
            IF(upgraded.gt.5000 .and. k3.lt.1000)THEN        !Mary changed the value from 1000 to 2500 concerning burn-in period 
            !2500  --- 250 for test
                CALL random_number(r) 
                IF(r.gt.0.75)THEN    !this part for 500 times simulation result to save      
                    ! change the 0.9 to 0.8      
                    ! 20% acceptance probability 
                    k3=k3+1
                    !WRITE(outfile1,'(A30,I4.4,A10)')  "simu_results/MCMC_dailyflux",k3,filesignal ! 根据K3定义文件名
                    !outfile1 = TRIM(ADJUSTL(outfile1))
                    !outfile1 = TRIM(ADJUSTL(fileplace))//outfile1
                    !write(*,*)outfile
                    !OPEN(62,file=outfile1,status='replace')
                    !DO i=1,365
                        !WRITE(62,602)i,(output_daily_mcmc(j,i),j=1,3) ! no need
!602                    !format((i7),",",14(f15.4,","))
                    !    WRITE(62,'(*( G0.8, :, ",", X))')(output_daily_mcmc(j,i),j=1,3)
                    !ENDDO
                    ! 确保这里的目录最后有 /
                    outfile2 = TRIM(ADJUSTL(fileplace)) // 'simu_results_hourly/hourlyflux'

                    WRITE(outfile2, '(A, I4.4, A)') TRIM(outfile2), k3, TRIM(filesignal)
                    ! 现在 outfile2 就是：../output/MCMC/teco_cnp/simu_results_hourly/hourlyflux0001_cnp.csv

                    OPEN(6622, file=TRIM(outfile2), status='replace', action='write', iostat=n)
                    IF (n /= 0) THEN
                    WRITE(*,*) 'ERROR: cannot open outfile2 = ', TRIM(outfile2), ' iostat=', n
                        STOP
                    ENDIF


                !     WRITE(outfile2,'(A30,I4.4,A10)')"simu_results_hourly/hourlyflux",k3,filesignal ! 根据K3定义文件名
                !     outfile2 = TRIM(ADJUSTL(outfile2))
                !     outfile2 = TRIM(ADJUSTL(fileplace))//outfile2
                !     write(*,*)outfile
                !     OPEN(6622,file=outfile2,status='replace')
                !     DO i=1,force_nhours
                !         WRITE(6622,'(*( G0.8, :, ",", X))')(output_mcmc(j,i),j=1,3)
                !     ENDDO
                                   
                !     CLOSE(62)  
                !     WRITE(63,'(*( G0.8, :, ",", X))')coef
                ENDIF
            ENDIF
    ELSE
        reject=reject+1
    ENDIF
   
    WRITE(*,*)'isimu,upgraded',isimu,upgraded

    IF(covexist.eq.0 .and. upgraded.gt.0 .and. mod(upgraded,ncov).eq.0)THEN
        covexist=1
        coefac=coefnorm
        CALL varcov(coefhistory,gamnew,npara,ncov)
        IF (.not.(all(gamnew==0.))) THEN
            gamma=gamnew
            CALL racine_mat(gamma,gamnew,npara)
            gamma=gamnew
        ENDIF
    ENDIF
   	
	IF (covexist.eq.1 .and. mod(upgraded,ncov).eq.0) THEN
        CALL varcov(coefhistory,gamnew,npara,ncov)
        IF (.not.(all(gamnew==0.))) THEN
            gamma=gamnew
            CALL racine_mat(gamma,gamnew,npara)
            gamma=gamnew
        ENDIF
	ENDIF               
    
    IF (upgraded .eq. 10000) exit
                     
    ENDDO ! END of isimu
    DO i=1,npara
        WRITE(72,*) (gamma(j,i),j=1,npara)
    ENDDO
    CLOSE(72)    
    CLOSE(71)
    WRITE(*,*)'run MCMC,CYCLE_CNP:',CYCLE_CNP        
END SUBROUTINE MCMC_simu


SUBROUTINE SensitivityTest()
    USE IntersVariables
    USE FileSize
    USE outputs_mod
    IMPLICIT NONE

    INTEGER,PARAMETER :: npara=24
    REAL coef(npara),coef_lower,coef_higher
    CHARACTER(len=80) outfile,outfile_state
    INTEGER n,i,k

        coef(1)=site_paras(32) !Q10=2.105
        coef(2)=site_paras(12)!SLA=150!!!!
        coef(3)=site_paras(19)!Vcmax0=33.721!!!!
        coef(4)=site_paras(23) !tau_Leaf=1.868             0.5,0.6,0.7                                 !para in 
        coef(5)=site_paras(24) !tau_Wood=40                                          	!para in 
        coef(6)=site_paras(25) !tau_Root=1.342      
        coef(7)=site_paras(56)                                  	
        coef(8)=site_paras(26) !tau_F=0.3                                          
        coef(9)=site_paras(27) !tau_C=5.86                                      	
        coef(10)=site_paras(28) !tau_Micro=0.356
        coef(11)=site_paras(29) !tau_SlowSOM=567.756!!!
        coef(12)=site_paras(30) !tau_Passive=2050!!!
        coef(13)=site_paras(42)!etaL
        coef(14)=site_paras(43)!etaW
        coef(15)=site_paras(44)!etaR
        coef(16)=site_paras(57)
        coef(17)=site_paras(31)!F2M
        coef(18)=site_paras(32)!C2M
        coef(19)=site_paras(33)!C2s
        coef(20)=site_paras(34)!M2S
        coef(21)=site_paras(35)!M2P
        coef(22)=site_paras(36)!S2M
        coef(23)=site_paras(37)!S2P
        coef(24)=site_paras(38)!P2M
        !coef(25)=InitialCNP(1,1)
        !coef(26)=InitialCNP(2,1)
        !coef(27)=InitialCNP(3,1)
        !coef(28)=InitialCNP(4,1)
        !coef(29)=InitialCNP(5,1)
        !coef(30)=InitialCNP(6,1)
        !coef(31)=InitialCNP(7,1)
        !coef(32)=InitialCNP(8,1)
        !coef(33)= InitialCNP(9,1)
        !coef(34)= site_paras(13)
        !coef(35)= site_paras(14)
        !coef(36)= site_paras(15)
        !coef(37)= site_paras(58)


        DO i=1,npara
            coef_lower =  coef(i)-coef(i)*0.25
            coef_higher = coef(i)+coef(i)*0.25

            coef(i) = coef_lower

            site_paras(32)=coef(1) !Q10=2.105
            site_paras(12)=coef(2)!SLA=150!!!!
            site_paras(19)=coef(3)!Vcmax0=33.721!!!!
            !site_paras(37)=coef(4)!hmax
            site_paras(23)=coef(4) !tau_Leaf=1.868             0.5,0.6,0.7                                 !para in 
            site_paras(24)=coef(5) !tau_Wood=40                                          	!para in 
            site_paras(25)=coef(6) !tau_Root=1.342     
            site_paras(56)=coef(7) !tau_Re                             	
            site_paras(26)=coef(8) !tau_F=0.3                                          
            site_paras(27)=coef(9) !tau_C=5.86                                      	
            site_paras(28)=coef(10) !tau_Micro=0.356
            site_paras(29)=coef(11) !tau_SlowSOM=567.756!!!
            site_paras(30)=coef(12) !tau_Passive=2050!!!
            site_paras(42)=coef(13)!etaL
            site_paras(43)=coef(14)!etaW
            site_paras(44)=coef(15)!etaR
            site_paras(57)=coef(16)!etaRe
            site_paras(45)=coef(17)!F2M
            site_paras(46)=coef(18)!C2M
            site_paras(47)=coef(19)!C2s
            site_paras(48)=coef(20)!M2S
            site_paras(49)=coef(21)!M2P
            site_paras(50)=coef(22)!S2M
            site_paras(51)=coef(23)!S2P
            site_paras(52)=coef(24)!P2M
            !InitialCNP(1,1)=coef(25)
            !InitialCNP(2,1)=coef(26)
            !InitialCNP(3,1)=coef(27)
            !InitialCNP(4,1)=coef(28)
            !InitialCNP(5,1)=coef(29)
            !InitialCNP(6,1)=coef(30)
            !InitialCNP(7,1)=coef(31)
            !InitialCNP(8,1)=coef(32)
            !InitialCNP(9,1)=coef(33)

            CALL CNP_simu()
            WRITE(*,'(A30,I2.2,A10)')"L/L_dailyflux",i,filesignal
            WRITE(outfile,'(A30,I2.2,A10)') "L/L_dailyflux",i,filesignal
            outfile = TRIM(ADJUSTL(outfile))
            outfile = TRIM(ADJUSTL(fileplace))//outfile

            write(*,*)'L - outfile',outfile
            
            OPEN(1100,file=outfile,status='replace')
            DO n=1,365
                WRITE(1100,'(*( G0.8, :, ",", X))')(output_daily(k,n),k=1,3)
            ENDDO
            close(1100)
            ! =========================================================
            ! L: write daily STATE (QC(1:9) + QPlab QPsorb QPss QPocc)
            !     -> output path is the SAME folder "L/"
            ! =========================================================
            WRITE(outfile_state,'(A30,I2.2,A10)') "L/L_state", i, filesignal
            outfile_state = TRIM(ADJUSTL(outfile_state))
            outfile_state = TRIM(ADJUSTL(fileplace))//outfile_state

            write(*,*) 'L - outstate', outfile_state

            OPEN(1200, file=outfile_state, status='replace')
            DO n = 1, 365
            WRITE(1200,'(*( G0.8, :, ",", X))')  &
                    (outputd_ccycle_Cpools(k,n), k=1,9), &   ! QC(1:9)
                    outputd_Pdynamic(6,n), &                 ! QPlab
                    outputd_Pdynamic(7,n), &                 ! QPsorb
                    outputd_Pdynamic(8,n), &                 ! QPss
                    outputd_Pdynamic(9,n)                    ! QPocc
            ENDDO
            CLOSE(1200)


            coef(i) = coef_higher

            site_paras(32)=coef(1) !Q10=2.105
            site_paras(12)=coef(2)!SLA=150!!!!
            site_paras(19)=coef(3)!Vcmax0=33.721!!!!
            !site_paras(37)=coef(4)!hmax
            site_paras(23)=coef(4) !tau_Leaf=1.868             0.5,0.6,0.7                                 !para in 
            site_paras(24)=coef(5) !tau_Wood=40                                          	!para in 
            site_paras(25)=coef(6) !tau_Root=1.342     
            site_paras(56)=coef(7) !tau_Re                             	
            site_paras(26)=coef(8) !tau_F=0.3                                          
            site_paras(27)=coef(9) !tau_C=5.86                                      	
            site_paras(28)=coef(10) !tau_Micro=0.356
            site_paras(29)=coef(11) !tau_SlowSOM=567.756!!!
            site_paras(30)=coef(12) !tau_Passive=2050!!!
            site_paras(42)=coef(13)!etaL
            site_paras(43)=coef(14)!etaW
            site_paras(44)=coef(15)!etaR
            site_paras(57)=coef(16)!etaRe
            site_paras(45)=coef(17)!F2M
            site_paras(46)=coef(18)!C2M
            site_paras(47)=coef(19)!C2s
            site_paras(48)=coef(20)!M2S
            site_paras(49)=coef(21)!M2P
            site_paras(50)=coef(22)!S2M
            site_paras(51)=coef(23)!S2P
            site_paras(52)=coef(24)!P2M
            !InitialCNP(1,1)=coef(25)
            !InitialCNP(2,1)=coef(26)
            !InitialCNP(3,1)=coef(27)
            !InitialCNP(4,1)=coef(28)
            !InitialCNP(5,1)=coef(29)
            !InitialCNP(6,1)=coef(30)
            !InitialCNP(7,1)=coef(31)
            !InitialCNP(8,1)=coef(32)
            !InitialCNP(9,1)=coef(33)
            

            CALL CNP_simu()
            WRITE(*,'(A30,I2.2,A10)')"H/H_dailyflux",i,filesignal
            WRITE(outfile,'(A30,I2.2,A10)')  "H/H_dailyflux",i,filesignal
            outfile = TRIM(ADJUSTL(outfile))
            outfile = TRIM(ADJUSTL(fileplace))//outfile
            write(*,*)'H - outfile',outfile
            OPEN(110,file=outfile,status='replace')
            DO n=1,365
                WRITE(110,'(*( G0.8, :, ",", X))')(output_daily(k,n),k=1,3)
            ENDDO
            close(110)

            WRITE(outfile_state,'(A30,I2.2,A10)') "H/H_state", i, filesignal
            outfile_state = TRIM(ADJUSTL(outfile_state))
            outfile_state = TRIM(ADJUSTL(fileplace))//outfile_state
            write(*,*) 'H - outstate', outfile_state

            OPEN(1300, file=outfile_state, status='replace')
            DO n = 1, 365
            WRITE(1300,'(*( G0.8, :, ",", X))')  &
                (outputd_ccycle_Cpools(k,n), k=1,9), &
                outputd_Pdynamic(6,n), &
                outputd_Pdynamic(7,n), &
                outputd_Pdynamic(8,n), &
                outputd_Pdynamic(9,n)
            ENDDO
            CLOSE(1300)

            PRINT*,i

        ENDDO

END SUBROUTINE SensitivityTest

END MODULE FOR_MCMC