

*     max. number of function evaluations
*     restart      


      SUBROUTINE LVMINI(NPAR,MVEC,NFCN, AUX) ! Initialization
*     ... of Limited-memory Variable-metrik Minimization
      DOUBLE PRECISION X(*),F,AUX(*),ALPHA,DGOUT,DGINIT
      DOUBLE PRECISION DBDOT,FCLAST,FPITER,FCDACT,AGRADL,PGRADL
      DOUBLE PRECISION FOPT,FPLUS,FMINUS,FDER,ALPHA0
      LOGICAL PRINT
      COMMON/LMCOMM/N,M,NB,IPAR,IPGR,IAPP,ISCA,ISEA,IFUN,IRHO,ISPT,IYPT,
     +              INGR,INDI,INHS,LAST,NCPT,ICPT,ICP,NPT,NCAL

      PARAMETER (NDIMLS=12)                         ! line search data
      COMMON/COLISE/NLS,XLS(NDIMLS),YLS(NDIMLS),TLS(NDIMLS),KLS(NDIMLS)
      PARAMETER (DFEPS=1.0E-4,DFWLF1=1.0E-4,DFWLF2=0.9) ! default const
      DOUBLE PRECISION TINY
      PARAMETER (TINY=1.0D-10)
      DATA CONVEG/DFEPS/,WOLFC1/DFWLF1/,WOLFC2/DFWLF2/
      SAVE
*     ...
      N=ABS(NPAR)             ! number of parameters 
      PRINT=.FALSE.
      IF(NPAR.LT.0) PRINT=.TRUE.

      LAST=LVMDIM(NPAR,MVEC)  ! necessary dimension

      M =MAX(1,ABS(MVEC))     ! number of history pairs, at least =1
      NF=NFCN                 ! maximum number of function evaluations
      IF(NF.LE.0) THEN        ! default value
         NF=91+(909*N*N)/(100+N*N)  ! 100 ... 1000 function evaluations 
         IF(MVEC.LT.0) NF=NF+N+N    ! + 2 * n for cov. matrix
      END IF 
      NCAL=0                  ! count number of function evaluations      

      ITER=0                  ! iteration counter 
      ICPT=0                  ! current pair-index (and index k)
      NCPT=0                  ! number of stored sy pairs

      IDER=0
 
      FPITER=0.0D0
      FCDACT=0.0D0
      ANGLGS=0.0
      ALPHA =1.0
      AGRADL=0.0D0
      PGRADL=0.0
      COSA  =1.0  
      LITER=-1

      DO I=1,LAST  
       AUX(I)=0.0D0 ! reset aux complete array
      END DO 
      RETURN
      

*     MTLINE search initialization -------------------------------------

*     (re)initialize     "Linesearch Algorithm with Sufficient Decrease"
 10   CALL MTLDEF(WOLFC1,WOLFC2,1.0E-17,1.0E-20,1.0E+20, 1,10)
      DGINIT=DBDOT(N,AUX,AUX(ISEA+1))  ! initial slope
      GOTO 40


*     ------------------------------------------------------------------
      ENTRY LVMEPS(EPS,WLF1,WLF2) ! change fit constants     
      IF(EPS.LE.0.0) THEN
         CONVEG=DFEPS   ! default value
      ELSE
         CONVEG=EPS     ! user value
      END IF
      IF(WLF1.LE.0.0) THEN
         WOLFC1=DFWLF1  ! default value
      ELSE
         WOLFC1=WLF1    ! user value
      END IF
      IF(WLF2.LE.0.0) THEN
         WOLFC2=DFWLF2  ! default value
      ELSE
         WOLFC2=WLF2    ! user value
      END IF
      RETURN
     
*     ------------------------------------------------------------------
      ENTRY LVMFUN(X,F,JRET,AUX)        ! called with new F and gradient

      IF(NCAL.EQ.0) THEN        ! function etc calculated the first time
 
         IF((.NOT.F+TINY.GT.0.0D0.AND..NOT.F-EPS.LT.0.0D0).OR.
     +       F.GT.1.0D+30.OR.F.LT.(-1.0D-30)) THEN  ! nan or infty
            IF(PRINT) WRITE(*,*) 'Bad first function value: ',F
            JRET=7
         END IF
       
         NB=0
         DO I=1,N
          IF(AUX(N+I).NE.0) NB=1 ! diagonal element not zero
          AUX(ISCA+I)=1.0D0      ! default scale factor 
          AUX(ISPT+I)=AUX(I)     ! for min derivative         
          AUX(IYPT+I)=AUX(I)     ! for max derivative
          AUX(IPAR+I)  =X(I)     ! initla parameter values
         END DO

         IF(PRINT) THEN
            WRITE(*,*) ' '
            WRITE(*,*) 'Limited-memory BFGS routine LVMINI initialized',
     +                 '========================'     
            WRITE(*,*)   '   with N =',N, '  dimension parameter'
            WRITE(*,*)   '        M =',M, '  number of history pairs'
            WRITE(*,191) '       NF =',NF,
     +                   '   max number of function calls'   
            WRITE(*,101) '      eps =',CONVEG,'   convergence parameter'
            WRITE(*,101) '      c_1 =',WOLFC1,'   first Wolfe parameter'
            WRITE(*,101)'      c_2 =',WOLFC2,'   second Wolfe parameter'
            IF(NB.EQ.1) THEN
               WRITE(*,*)   'Diagonal of Hessian provided' 
            END IF
            WRITE(*,*) 'Requires AUX array with at least',LAST,' words'
            CALL LVMPRT(6,AUX,-1) ! print start values
            WRITE(*,*) ' '
            WRITE(*,112) 
         END IF

      END IF ! function etc. calculated the first time 

      IF(NCPT.EQ.0) THEN     ! (s,y) pairs not yet stored 
         IF(NB.NE.0) THEN    ! diagonal elements present
            DO I=1,N
             IF(AUX(N+I).GT.0.0) THEN
                AUX(ISCA+I)=SQRT((1.0D0+ABS(F))/AUX(N+I)) 
             END IF
            END DO
         ELSE                ! no diagonal elements   
            DO I=1,N
             AUX(ISPT+I)=MIN(AUX(ISPT+I),AUX(I)) ! min
             AUX(IYPT+I)=MAX(AUX(IYPT+I),AUX(I)) ! max
             FDER=0.5D0*(ABS(AUX(ISPT+I))+ABS(AUX(IYPT+I)))
             IF(FDER.NE.0.0) AUX(ISCA+I)=1.0D0/FDER
            END DO
         END IF 
      END IF ! (s,y) pairs not yet stored

      NCAL=NCAL+1 ! count nr of function evaluations

      IF(IDER.EQ.0) GOTO 30 ! derivative caculation

*     derivatives ------------------------------------------------------
*
*     using: IPAR ... optimal parameter values, with FOPT
*            IPGR ... previous gradient (optimal gradient)
*            ISEA ... search vector (gradient at +)  
*            INGR ... gradient, numerical
*            INHS ... Hessian, numerical

 20   IF(IDER.EQ.0) THEN              ! start numerical derivatives

*        prepare error calculation

         KRET=JRET ! save return code
         JRET=-1   ! set return code to function evaluation

         ETA =1.0E-5         ! relative function error
         NDE =1              ! counter for derivatives - not used

         FOPT=F              ! save optimal function value
         DO I=1,N
          AUX(IPAR+I)=X(I)   ! save optimal values
          AUX(IGPR+I)=AUX(I) ! save gradient at optimal values
          AUX(INGR+I)=0.0D0  ! reset gradient, numerical  
         END DO
         DO I=1,(N*N+N)/2
          AUX(INHS+I)=0.0D0  ! reset Hessian, numerical
         END DO
         IF(PRINT) THEN
            WRITE(*,*) ' '
            WRITE(*,*) 'derivative calculation:'
            WRITE(*,*) '  par_i   step(+-)     secder   dF_minus',
     +      '    dF_plus'
         END IF
      ELSE IF(IDER.GT.0) THEN
         FPLUS=F                           ! F_+
         DO I=1,N
          AUX(ISEA+I)=AUX(I)               ! g_+                   
         END DO
         X(IDER)=AUX(IPAR+IDER)-HSTEP      ! X - h_j
         IDER=-IDER
         GOTO 90
      ELSE IF(IDER.LT.0) THEN
         FMINUS=F                          ! F_-
         IDER=-IDER
         X(IDER)=AUX(IPAR+IDER)            ! restore X 

         DO I=1,N
          HIJ=0.5D0*(AUX(ISEA+I)-AUX(I))/HSTEP
          IF(I.LT.IDER) THEN
             IJ=INHS+I+(IDER*IDER-IDER)/2
             AUX(IJ)=0.5D0*(AUX(IJ)+HIJ)     ! average
          ELSE IF(I.EQ.IDER) THEN
             AUX(INHS+(I*I+I)/2)=HIJ         ! diagonal element
          ELSE IF(I.GT.IDER) THEN
             AUX(INHS+IDER+(I*I-I)/2)=HIJ    ! store
          END IF
         END DO
         HJJ=(FPLUS-2.0D0*FOPT+FMINUS)/HSTEP**2
         AUX(INDI+IDER)=HJJ 

         IF(PRINT) THEN 
            WRITE(*,121) IDER,HSTEP,HJJ,FMINUS-FOPT,FPLUS-FOPT
         END IF

         IF(IDER.EQ.N) THEN                  ! numerical work done
*           restore
            F=FOPT               ! optimal function value
            DO I=1,N
             AUX(I)=AUX(IPGR+I)  ! optimal gradient
            END DO
*                                           diagonal    pointer
            CALL LVMINV(AUX(INHS+1),N,NRANK,AUX(ISEA+1),AUX(IPGR+1))
            DO I=1,N
             DEN=AUX(INHS+(I*I+I)/2)             ! diagonal element
             AUX(INGR+I)=0.0D0
             IF(DEN.GT.0.0D0) AUX(INGR+I)=SQRT(DEN)  ! error
             DEN=DEN*AUX(ISEA+I)          
             AUX(INDI+I)=0.0D0
             IF(DEN.GT.0.0D0) AUX(INDI+I)=SQRT(1.0D0-1.0D0/DEN) !glcorr
            END DO

            IOPT=7
            IF(N.GT.12) IOPT=3 

            IF(PRINT) THEN
               CALL LVMPRT(6,AUX,6)
               WRITE(*,*) 'Limited-memory BFGS routine ending ',
     +                 '==================================='
               WRITE(*,*) ' '
            END IF     
            JRET=KRET
            GOTO 90
         END IF 

      END IF 

      IDER=IDER+1                     ! next derivative
      HSTEP=AUX(IAPP+IDER) 
      X(IDER)=AUX(IPAR+IDER)+HSTEP    ! X + h_j
      GOTO 90

*     minimization -----------------------------------------------------

 30   FCDACT=FCLAST-F                       ! actual function difference
      AGRADL=SQRT(DBDOT(N,AUX,AUX))         ! actual |gradient| length
      IF(ITER.EQ.0) FCDACT=0.0D0
      IF(PRINT) THEN
         IF(ITER.NE.LITER) THEN
            WRITE(*,111) ITER,NCAL,NCPT,F,0.5*FPITER,FCDACT,
     +                   ALPHA,ANGLGS
         ELSE
            WRITE(*,113)      NCAL,NCPT,F,   FCDACT,ALPHA
         END IF
      END IF  

      IF(AGRADL.EQ.0.0D0) THEN 
         IF(PRINT) WRITE(*,*) 'gradient zero'
         CALL LVMOUT(N,M,AUX)
         JRET=+2   ! gradient zero 
         GOTO 90    ! error !
      END IF

      LITER=ITER

      JRET=-1

      IF(NCAL.EQ.1) GOTO 50 ! first function value: calculate searchdir

*     continue line search

 40   IF(DGINIT.GE.0.0D0) THEN         ! ... is not negative 
         WRITE(*,*) 'Huch ',DGINIT
         WRITE(*,*) 'AUX ',AUX(1),AUX(2)
         WRITE(*,*) 'D   ',AUX(ISEA+1),AUX(ISEA+2)
         CALL LVMOUT(N,M,AUX)
         JRET=+3
         GOTO 90 ! error !
      END IF
      CALL MTLINE(N,X,F,AUX(1),AUX(ISEA+1),ALPHA, AUX(IPAR+1),
     +            DGOUT,IRET)
*     value of IRET=-1   return for another function evaluation
*                  = 0   N or STP bad
*                  = 1   conditions hold
*                  = 2   interval at limit
*                  = 3   max number of calls
*                  = 4   lower bound
*                  = 5   upper bound  
*                  = 6   rounding error limitation
*                  = 7   all below 10^-5  

      IF(IRET.LT.0) GOTO 90 ! ... for another function evaluation    

#ifdef graline 
*     print line search   
c      WRITE(*,103)
c      DO I=1,NLS
c       WRITE(*,102) I,XLS(I),YLS(I),TLS(I),KLS(I)
c      END DO
c      WRITE(*,*) 'Parameter values',X(1),X(2) 
 103  FORMAT('     I     X    Y    T    N ')
 102  FORMAT(I6,G15.5,G15.5,G15.5,I8)

      IF(NLS.GE.2) THEN  ! plot line search 
         CALL QFINIT(0.0,0.0)
         CALL QUAT('red circle marker')
         CALL QUAT('red full line')
         DO I=1,NLS
          CALL QDNXY(1,XLS(I),YLS(I),'ml')
         END DO
         CALL QUAT('blue bars line')
         DX=0.05*XLS(NLS)
         DO I=1,NLS
          CALL QD2XY(XLS(I)-DX,YLS(I)-DX*TLS(I),
     +               XLS(I)+DX,YLS(I)+DX*TLS(I),'l')
         END DO
         CALL QUPDAT
      END IF
#endif 

      IF(PRINT.AND.IRET.NE.1) CALL MTINFO(IRET) ! print error info 
      IF(IRET.EQ.0) GOTO 90                     ! bad N or STP 
      IF(IRET.NE.1) GOTO 60 ! conditions do NOT hold

*     ******************************************************************
*     calculation of new search vector by the L-BFGS algorithm ---------        
*     ****************************************************************** 


*     store s- and y-vector, calculate rho from the two vectors -------- 

      IF(NCAL.GT.1) THEN
         CALL LVMUPD(ALPHA,AUX(ISEA+1),AUX,AUX(IPGR+1),AUX) ! store next pair
      END IF  

*     L-BFGS two-loop recursion ========================================

 50   DO I=1,N
       AUX(ISEA+I)=AUX(I)             ! copy gradient  
      END DO
      IF(ITER.EQ.0) THEN
         DO I=1,N
          AUX(ISEA+I)=AUX(ISEA+I)/AGRADL  ! adjust gradient length 
         END DO
      ELSE
         CALL LVMSOL(AUX(ISEA+1),AUX)     ! get next search vector
      END IF

*     ==================================================================

      DO I=1,N
       AUX(ISEA+I)=-AUX(ISEA+I)   ! reverse sign of search vector
       AUX(IPGR+I)=AUX(I)         ! store last gradient
      END DO
      PGRADL=AGRADL              ! store last gradient length
      ANGLGS=57.2957795*
     +       ACOS(-DBDOT(N,AUX(ISEA+1),AUX)/
     +      (AGRADL*SQRT(DBDOT(N,AUX(ISEA+1),AUX(ISEA+1))))) 

*     start new line search iteration ----------------------------------

      FPLAST=FPITER              ! previous gradient * step vector 
      FPITER=-DBDOT(N,AUX,AUX(ISEA+1))     ! gradient * step vector

      ALPHA=0.1D0                     ! default start value
      IF(ITER.GT.0) THEN 
         ALPHA0=2.0D0*(FCLAST-F)/FPITER
         IF(ALPHA0.GT.0.0D0) ALPHA=MIN(1.0D0,1.01D0*ALPHA0)        
      END IF
      ITER=ITER+1
      FCLAST=F        ! previous function value
      IF(ITER.LE.M.OR.MAX(FPITER,FCDACT).GT.CONVEG) GOTO 10

 60   JRET=+1                  ! iteration ending

      AUX(IFUN)=F              ! save optimal function value

*     approximate errors from pair info -------------------------------- 

      DO J=1,N ! use pair-info to estimate the inverse Hessian
       DO I=1,N
        AUX(ISEA+I)=0.0
       END DO                       ! unit vector e_j
       AUX(ISEA+J)=1.0
       CALL LVMSOL(AUX(ISEA+1),AUX) ! solve -> matrix column 
       AUX(IAPP+J)=SQRT(ABS(AUX(ISEA+J))) ! store sqrt(diagonal element)
      END DO  
      IF(INGR.NE.0) GOTO 20 ! derivative calculation
      IF(PRINT) THEN
         CALL LVMPRT(6,AUX,0) ! final function value ...  
         CALL LVMPRT(6,AUX,2) ! and parameters with approx. errors  
      END IF 
*     convergence reached ---------------------------------------------- 

 90   CONTINUE
      IF(JRET.LT.0) THEN
         DO I=1,N+N
          AUX(I)=0.0D0   ! reset AUX
         END DO 
      END IF 

*     Format statements ------------------------------------------------ 

 101  FORMAT(1X,A,G10.3,A)
 191  FORMAT(1X,A,I10,A)
 112  FORMAT('Iter Nfcn Npair',6X,
     +       'Fcn_value   dFcn__ex   dFcn__ac   fstep',3X,' degr')
 121  FORMAT(I8,G11.2,E11.2,2E11.2)
 111  FORMAT(I4,I5,I6,F15.4,    2E11.2,F8.3,F8.2) 
 113  FORMAT(4X,I5,I6,F15.4,11X,1E11.2,F8.3)
      END 

      SUBROUTINE LVMOUT(N,M,AUX) ! test printout
      DOUBLE PRECISION AUX(1000)
*
      WRITE(*,*) ' '
      WRITE(*,*) '-----------------------------'
      IND=0
      WRITE(*,*) '   Gradient:'
      WRITE(*,*) (AUX(IND+I),I=1,N)
      IND=IND+N
      WRITE(*,*) '   Diagonal elements:'
      WRITE(*,*) (AUX(IND+I),I=1,N)
      IND=IND+N
      WRITE(*,*) '   Parameter:'
      WRITE(*,*) (AUX(IND+I),I=1,N)
      IND=IND+N

      WRITE(*,*) '   Pgradient:'
      WRITE(*,*) (AUX(IND+I),I=1,N)
      IND=IND+N

      WRITE(*,*) '   Scale factors:'
      WRITE(*,*) (AUX(IND+I),I=1,N)
      IND=IND+N

      WRITE(*,*) '   Search:'
      WRITE(*,*) (AUX(IND+I),I=1,N)
      IND=IND+N

      WRITE(*,*) '   Rho(k):'
      WRITE(*,*) (AUX(IND+I),I=1,M)
      IND=IND+M
      WRITE(*,*) '   Alpha(k):'
      WRITE(*,*) (AUX(IND+I),I=1,M)
      IND=IND+M

      WRITE(*,*) '   s vectors:'
      DO J=1,M
       WRITE(*,*) J,(AUX(IND+I),I=1,N)
       IND=IND+N
      END DO

      WRITE(*,*) '   y vectors:'
      DO J=1,M
       WRITE(*,*) J,(AUX(IND+I),I=1,N)
       IND=IND+N
      END DO
      WRITE(*,*) '-----------------------------'
      WRITE(*,*) ' '
      END

      FUNCTION LVMDIM(NPAR,MVEC) ! pointer and total mem calculation
      COMMON/LMCOMM/N,M,NB,IPAR,IPGR,IAPP,ISCA,ISEA,IFUN,IRHO,ISPT,IYPT,
     +              INGR,INDI,INHS,LAST,NCPT,ICPT,ICP,NPT,NCAL
      SAVE
*     ... 
      NP=ABS(NPAR)      
      MV=ABS(MVEC)

*     offset                  nr_of_words   description
*     ------------------      -----------   -------------------------- 
*          0                !     N         gradient 
*          NP               !     N         diagonal elements
      IPAR=NP+NP            !     N         parameter vector
      IPGR=IPAR+NP          !     N         previous gradient
      IAPP=IPGR+NP          !     N         approximate errors
      ISCA=IAPP+NP          !     N         scale factors
      ISEA=ISCA+NP          !     N         search vector
      IFUN=ISEA+NP+1        !     1         function value  
      IRHO=IFUN             !     M         rho(k)
!          IRHO+MV          !     M         alpha(k)
      ISPT=IRHO+MV+MV       !     M*N       s vectors
      IYPT=ISPT+MV*NP       !     M*N       y vectors

      INGR=0                !     unused
      INDI=0                !     unused
      INHS=0                !     unused
      LAST=IYPT+MV*NP       !     total space

      IF(MVEC.LT.0) THEN
         INGR=IYPT+MV*NP    !     N         gradient, numerical
         INDI=INGR+NP       !     N         diagonal, numerical 
         INHS=INDI+NP       !  (N*N+N)/2    Hessian, numerical              
         LAST=INHS+(NP*NP+NP)/2
      END IF
      LVMDIM=LAST
      END 


      FUNCTION LVMIND(IARG) ! return index for ...
*     IARG = 0   function value
*          = 1   parameter vector
*          = 2   approximate errors
*          = 3   accurate errors
*          = 4   global correlations
*          = 5   covariance matrix
      COMMON/LMCOMM/N,M,NB,IPAR,IPGR,IAPP,ISCA,ISEA,IFUN,IRHO,ISPT,IYPT,
     +              INGR,INDI,INHS,LAST,NCPT,ICPT,ICP,NPT,NCAL
*     ...
      LVMIND=0         
      IF(IARG.EQ.(-1)) THEN
         LVMIND=N
      ELSE IF(IARG.EQ.(-2)) THEN
         LVMIND=NCAL
      ELSE IF(IARG.EQ.0) THEN
         LVMIND=IFUN
      ELSE IF(IARG.EQ.1) THEN 
         LVMIND=IPAR
      ELSE IF(IARG.EQ.2) THEN 
         LVMIND=IAPP
      ELSE IF(IARG.EQ.3) THEN 
         LVMIND=INGR
      ELSE IF(IARG.EQ.4) THEN 
         LVMIND=INDI
      ELSE IF(IARG.EQ.5) THEN 
         LVMIND=INHS
      END IF
      END

      SUBROUTINE LVMPRT(LUN,AUX,JARG) ! print on unit LUN ... 
      DOUBLE PRECISION AUX(*)
      INTEGER MC(10)
      SAVE
*          =-1  initial parameter values
*     IARG = 0  function value
*          = 1  parameter values
*          = 2  parameter values with approximate errors
*          = 3  parameter values with errors
*          = 4  parameter values with global correlations and errors
*          = 5  correlation matrix
*          = 6  means 0, 2, 4, 5  
*     ...
      LUP=LUN
      IF(LUP.EQ.0) LUP=6 ! standard output
      NFUN=LVMIND(-2)   
      NPAR=LVMIND(-1)
      IFUN=LVMIND(0)   ! get indices
      IPAR=LVMIND(1) 
      IAPP=LVMIND(2)
      INGR=LVMIND(3) 
      INDI=LVMIND(4)
      INHS=LVMIND(5) 

      IF(JARG.NE.6) THEN      
         IARG1=JARG
         IARG2=JARG
      ELSE
         IARG1=0
         IARG2=6  
      END IF

      DO IARG=IARG1,IARG2,2
      WRITE(LUP,*) ' '
      IF(IARG.EQ.0) THEN
         WRITE(LUP,*) NPAR,'-parameter function value =',AUX(IFUN),
     +                     ' after',NFUN,' function evaluations' 
      ELSE IF(IARG.EQ.(-1)) THEN  ! parameter values
         WRITE(LUP,*) '  Initial parameter values:'
         WRITE(LUP,101)    
         WRITE(LUP,102) (I,AUX(IPAR+I),I=1,NPAR)
      ELSE IF(IARG.EQ.1) THEN  ! parameter values 
         WRITE(LUP,101) 
         WRITE(LUP,102) (I,AUX(IPAR+I),I=1,NPAR)
      ELSE IF(IARG.EQ.2) THEN  ! parameter values with appr. errors
         WRITE(LUP,*) '  Optimization results ',
     +                 'with approximate errors:'
         WRITE(LUP,103)
         WRITE(LUP,104) (I,AUX(IPAR+I),AUX(IAPP+I),I=1,NPAR)
      ELSE IF(IARG.EQ.3.OR.IARG.EQ.4) THEN
         WRITE(LUP,*) '    Optimization results:'  
         IF(IARG.EQ.3) THEN
            WRITE(LUP,*)'    par_i         value             error'
         ELSE IF(IARG.EQ.4) THEN
            WRITE(LUP,*)'    par_i         value             error',
     +      '        glcorr'
         END IF
         DO I=1,NPAR 
          IF(IARG.EQ.3) THEN
             WRITE(LUP,105) I,AUX(IPAR+I),AUX(INGR+I)
          ELSE IF(IARG.EQ.4) THEN
             WRITE(LUP,105) I,AUX(IPAR+I),AUX(INGR+I),AUX(INDI+I) 
          END IF
         END DO 
      ELSE IF(IARG.EQ.5.OR.IARG.EQ.6) THEN         
         WRITE(LUP,*) '    Correlation coefficients [per mille]',
     +              ' between parameters i and j:'
         WRITE(LUP,106) (I,I=1,MIN(10,NPAR))
         DO I=1,NPAR
          DEI=AUX(INGR+I)
          IF(DEI.GT.0.0D0) THEN
             M=0
             DO J=1,I
              RHO=0.0D0
              DEJ=AUX(INGR+J)
              IF(DEI*DEJ.GT.0.0D0) THEN
                 RHO=AUX(INHS+J+(I*I-I)/2)/(DEI*DEJ)
              END IF
              M=M+1
              IF(RHO.GT.0.0D0) THEN
                 MC(M)=1000.0D0*RHO+0.5D0
              ELSE
                 MC(M)=-1000.0D0*RHO+0.5D0
                 MC(M)=-MC(M)
              END IF
              IF(M.EQ.10) THEN
                 IF(J.LE.10) THEN
                    WRITE(LUP,107) I,MC
                 ELSE
                    WRITE(LUP,108)   MC
                 END IF
                 M=0
              END IF
             END DO
             IF(M.NE.0) THEN
                IF(J.LE.10) THEN
                   WRITE(LUP,107) I,(MC(L),L=1,M)
                ELSE
                   WRITE(LUP,108)   (MC(L),L=1,M)
                END IF
                M=0
             END IF
          ELSE
             WRITE(LUP,107) I
          END IF
         END DO
      END IF

      END DO 
      WRITE(LUP,*) ' '

*     Formats 
 101  FORMAT(3('   par_i       value    '))
 102  FORMAT(3(I8,G16.6))
 103  FORMAT(2('   par_i      value         error    '))
 104  FORMAT(2(I8,G15.6,' +-',G11.3))
 105  FORMAT(I10,G18.6,'+-',G16.6,F10.3)
 106  FORMAT('     par_i  \ j=',10I6)
 107  FORMAT(I10,6X,10I6)
 108  FORMAT(16X,   10I6)
      END 


      SUBROUTINE LVMUPD(FAC,SVEC,YVECP,YVECM,AUX)          ! pair update
      DOUBLE PRECISION FAC,SVEC(*),YVECP(*),YVECM(*),AUX(*),DBDOT,D(*)
      DOUBLE PRECISION DOTSY,DOTYY,DOTSBY
      COMMON/LMCOMM/N,M,NB,IPAR,IPGR,IAPP,ISCA,ISEA,IFUN,IRHO,ISPT,IYPT,
     +              INGR,INDI,INHS,LAST,NCPT,ICPT,ICP,NPT,NCAL
      SAVE
*     ... 

*     update set of vector pairs ---------------------------------------

      NPT =ICPT*N           ! offset before increasing ICPT
      NCPT=MIN(NCPT+1,M)    ! number of stored vectors
      ICPT=ICPT+1           ! increase pointer   ICPT=POINT
      IF(ICPT.EQ.M) ICPT=0       
      ICP =ICPT             ! pointer
      IF(ICP.EQ.0) ICP=M

      DO I=1,N    ! new vectors with scaling factors
       AUX(ISPT+NPT+I)=FAC*SVEC(I)      /AUX(ISCA+I)   ! param diff.
       AUX(IYPT+NPT+I)=(YVECP(I)-YVECM(I))*AUX(ISCA+I) ! gradient diff. 
      END DO

      DOTYS=DBDOT(N,AUX(ISPT+NPT+1),AUX(IYPT+NPT+1))
      DOTYY=DBDOT(N,AUX(IYPT+NPT+1),AUX(IYPT+NPT+1))
      DOTSBY=0.0D0 
      IF(DOTYY.NE.0.0D0) DOTSBY=DOTYS/DOTYY

      AUX(IRHO+ICP)=0.0D0 ! rho(k)=1/(s(k)*y(k))
      IF(DOTYS.NE.0.0D0) AUX(IRHO+ICP)=1.0D0/DOTYS
      RETURN

*     ------------------------------------------------------------------

      ENTRY LVMSOL(D,AUX)       ! calculate B (vector d)

      DO I=1,N
       D(I)=D(I)*AUX(ISCA+I)    ! scaling
      END DO 

*     backward: alpha(i) := rho(i) * s(i)T*y(i); q := q-alpha(i)*y(i)

      ICP=ICPT 
      DO J=1,NCPT               ! backward loop
       ICP=ICP-1                ! ... with step -1 
       IF(ICP.LT.0) ICP=M-1     ! calculate alpha(j)    
       AUX(IRHO+M+ICP+1)=AUX(IRHO+ICP+1)*DBDOT(N,AUX(ISPT+ICP*N+1),D)
       CALL DBAXPY(N,-AUX(IRHO+M+ICP+1),AUX(IYPT+ICP*N+1),D)
      END DO

*     solve with approx Hessian ---------------------------------------- 

      IF(NB.EQ.0) THEN          ! no second derivative    
         DO I=1,N
          D(I)=DOTSBY*D(I)      ! scale gradient
         END DO 
      ELSE IF(NB.EQ.1) THEN     ! diagonal matrix 
         DO I=1,N               ! solve for approx. search direction
          IF(AUX(N+I).EQ.0.0D0) THEN
             D(I)=0.0D0         ! diagonal element zero 
          ELSE
             D(I)=D(I)/ABS(AUX(N+I))  ! using aprox. Hessian ...
             D(I)=D(I)/AUX(ISCA+I)**2 ! with scaling 
          END IF
         END DO 
      ELSE IF(NB.GT.1) THEN     ! band matrix
         STOP 'not implemented'
      END IF

*     forward loop -----------------------------------------------------

      DO J=1,NCPT               ! forward loop
       BETA=AUX(IRHO+ICP+1)*DBDOT(N,AUX(IYPT+ICP*N+1),D) 
       CALL DBAXPY(N,AUX(IRHO+M+ICP+1)-BETA,AUX(ISPT+ICP*N+1),D)
       ICP=ICP+1                ! ... with step +1 
       IF(ICP.EQ.M) ICP=0    
      END DO

      DO I=1,N
       D(I)=D(I)*AUX(ISCA+I)    ! scaling back
      END DO 

      END




*     matrix/vector products -------------------------------------------

      DOUBLE PRECISION FUNCTION DBDOT(N,DX,DY)   ! dot product
*     Dot product of two vectors: DX * DY
      DOUBLE PRECISION DX(*),DY(*),DTEMP
*     ...
      DTEMP=0.0D0
      DO I = 1,MOD(N,5)
       DTEMP=DTEMP+DX(I)*DY(I)
      END DO
      DO I =MOD(N,5)+1,N,5
       DTEMP=DTEMP+DX(I)*DY(I)+DX(I+1)*DY(I+1)+DX(I+2)*DY(I+2)
     +                        +DX(I+3)*DY(I+3)+DX(I+4)*DY(I+4)
      END DO
      DBDOT=DTEMP
      END


      SUBROUTINE DBAXPY(N,DA,DX,DY)              ! multiply, addition
*     Constant times vector added to a vector: DY:=DY+DA*DX
      DOUBLE PRECISION DX(*),DY(*),DA
*     ...
      DO I=1,MOD(N,4)
        DY(I)=DY(I)+DA*DX(I)
      END DO
      DO I=MOD(N,4)+1,N,4
        DY(I  )=DY(I  )+DA*DX(I  )
        DY(I+1)=DY(I+1)+DA*DX(I+1)
        DY(I+2)=DY(I+2)+DA*DX(I+2)
        DY(I+3)=DY(I+3)+DA*DX(I+3)
      END DO
      END

*     matrix inversion ------------------------------------------------- 

      SUBROUTINE LVMINV(V,N,NRANK,DIAG,NEXT)   ! matrix inversion  
*           V = symmetric N-by-N matrix in symmetric storage mode
*               V(1) = V11, V(2) = V12, V(3) = V22, V(4) = V13, . . .
*               replaced by inverse matrix
*
*     DIAG(N) =  double precision scratch array
*     NEXT(N) =  aux array (double precision instead of integer)
*
      DOUBLE PRECISION V(*),DIAG(N),VKK,VJK,EPS,NEXT(N)
      PARAMETER (EPS=1.0D-10)
*     ...
      NEXT0=1                     ! start index
      DO I=1,N
       NEXT(I)=I+1                ! set "next" pointer
       DIAG(I)=ABS(V((I*I+I)/2))  ! save abs of diagonal elements
      END DO
      NEXT(N)=-1                  ! end flag

      NRANK=0
      DO I=1,N                    ! start of loop
       K  =0
       VKK=0.0D0

       J=NEXT0
       LAST=0 
 05    IF(J.GT.0) THEN
          JJ=(J*J+J)/2
          IF(ABS(V(JJ)).GT.MAX(ABS(VKK),EPS*DIAG(J))) THEN
             VKK=V(JJ)              
             K=J
             L=LAST
          END IF
          LAST=J
          J=NEXT(LAST)
          GOTO 05   
       END IF 

       IF(K.NE.0) THEN            ! pivot found
          KK=(K*K+K)/2
          IF(L.EQ.0) THEN
             NEXT0=NEXT(K)
          ELSE 
             NEXT(L)=NEXT(K)
          END IF
          NEXT(K)=0.0D0           ! index is used, reset
          NRANK=NRANK+1           ! increase rank and ...
          VKK    =1.0/VKK
          V(KK)  =-VKK
          JK     =KK-K
          JL     =0
          DO J=1,N                ! elimination
           IF(J.EQ.K) THEN
              JK=KK
              JL=JL+J
           ELSE
              IF(J.LT.K) THEN
                 JK=JK+1
              ELSE
                 JK=JK+J-1
              END IF
              VJK  =V(JK)
              V(JK)=VKK*VJK
              LK   =KK-K
              DO L=1,J
               JL=JL+1
               IF(L.EQ.K) THEN
                  LK=KK
               ELSE
                  IF(L.LT.K) THEN
                     LK=LK+1
                  ELSE
                     LK=LK+L-1
                  END IF
                  V(JL)=V(JL)-V(LK)*VJK
               END IF
              END DO
           END IF
          END DO
       ELSE
          DO K=1,N
           IF(NEXT(K).NE.0.0D0) THEN
              DO J=1,K
               IF(NEXT(J).NE.0.0D0) V((K*K-K)/2+J)=0.0D0  ! clear row/col
              END DO
           END IF
          END DO
          GOTO 10
       END IF
      END DO             ! end of loop
 10   DO IJ=1,(N*N+N)/2
       V(IJ)=-V(IJ)      ! finally reverse sign of all matrix elements
      END DO
      END


      SUBROUTINE LVMIDI(NPAR,PAR,FUNCT, AUX) ! minimization of [Minuit] function
      DOUBLE PRECISION PAR(*),FSUM,AUX(*)
      SAVE
      EXTERNAL FUNCT
*     ...      
c      CALL LVMEPS(1.0E-4,0.0,0.9)

      MVEC=MIN(ABS(NPAR),6)                ! number of delta pairs
      CALL LVMINI(NPAR,-MVEC,0, AUX) ! initialization

      NFCNS=0
*     function and gradient evaluation
 10   CALL FUNCT(ABS(NPAR),AUX,FSUM,PAR,2)            ! calculate F, gradient
      NFCNS=NFCNS+1
      CALL LVMFUN(PAR,FSUM, IRET,AUX)
      IF(IRET.LT.0)          GOTO 10
      END


*     line search subprograms ------------------------------------------

      SUBROUTINE MTLINE(N,X,F,G,S,STP, WA, DGOUT,INFO)
*     Linesearch Algorithms wih Guaranteed Sufficient Decrease"
*     by JJ More and D Thuente, ACM Transactions on Mathenatical
*     Software 20 (1994), pp 286-307 - modified by V.B.
*
*        ===============================================================
*        The algorithm is designed to find a step which satisfies 
*        the sufficient decrease condition
*              F(X+STEP*S) <= F(X) + FTOL*STP*[gradF(X).S]
*        and the curvature condition
*              ABS[gradF(X+STP*S).S] <= GTOL*ABS[gradF(X).S]
*
*        Typical values for Constants: 
*                               FTOL = 1.0E-4
*                               GTOL = 0.9  for Newton, Quasi-Newton
*                                    = 0.1  for CG
*
*        Conditions: 0 < FTOL < GTOL <  1   for Newton, Quasi-Newton
*                    0 < FTOL < GTOL < 1/2  for CG
*        ===============================================================
*      
*     CALL MTLDEF(...)      ! initialization, see end of routine
*     NC=0                  ! reset
*     define STP =  (in units of S(.), usually 1.0) 
* 01  calculate     (note: first value before first MTLINE call) 
*             F    = function value,
*             G(.) = gradient,
*             S(.) = search vector, at X(.) 
*     NC=NC+1               ! count function calls  
*     CALL MTLINE(N,X,F,G,S,STP, WA, DGOUT,INFO)
*        WA(.) =work array, of dimension N
*        DGOUT =value of directional derivative when the Wolfe 
*               conditions hold; may be used to check descent     
*        INFO  =returned flag (not to be modified) 
*     IF(INFO.LT.0) GOTO 01 ! recalculate function ...  
*     IF(INFO.NE.1)         ! error condition
*     IF(INFO.EQ.1)         ! sufficient decrease and deriv. condition
*                           ! i.e. Wolfe conditions hold 
*     CALL MTINFO(INFO)     ! to print end condition
*     ------------------------------------------------------------------ 

      INTEGER N,MAXF,INFO,NFEV
      DOUBLE PRECISION F,STP,FTOL,GTOL,XTOL,STPMIN,STPMAX
      DOUBLE PRECISION X(N),G(N),S(N),WA(N),DBDOT

      DOUBLE PRECISION SFD(3,5) ! step|function_value|slope 
      DOUBLE PRECISION A1,A2,A3,A4,FOPT,SOPT

      SAVE

      INTEGER INFOC,J
      LOGICAL BRACKT,STAGE1
      DOUBLE PRECISION DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM,
     *       FINIT,FTEST1,FM,FX,FXM,FY,FYM,P5,P66,STX,STY,
     *       STMIN,STMAX,WIDTH,WIDTH1,XTRAPF,ZERO,DGOUT,STPAC

      PARAMETER (NDIMLS=12)
      COMMON/COLISE/NLS,XLS(NDIMLS),YLS(NDIMLS),TLS(NDIMLS),KLS(NDIMLS)
 
      DATA SFD/15*0.0D0/

      DATA P5,P66,XTRAPF,ZERO /0.5D0,0.66D0,4.0D0,0.0D0/
c      DATA P5,P66,XTRAPF,ZERO /0.5D0,0.66D0,10.0D0,0.0D0/
c      DATA P5,P66,XTRAPF,ZERO /0.5D0,0.66D0,2.618D0,0.0D0/
  
      SPAN(A1,A2,A3,A4)=MAX(A1,A2,A3,A4)-MIN(A1,A2,A3,A4)

*     ...

*     -------------------------  initialization  ----------------------
      IF(NFEV.GE.0) GOTO 10
      NFEV =0     ! reset  
      INFOC=1
      INFO =0 
      IF(N.LE.0.OR.STP.LE.ZERO) THEN 
         WRITE(*,*) 'MTLINE call with N=',N,', STP=',STP
         RETURN ! check N and step value
      END IF 
      DGINIT=DBDOT(N,G,S)
 
      NLS=1         ! store initial value
      XLS(1)=0.0
      YLS(1)=0.0
      TLS(1)=DGINIT
      KLS(1)=0 

      IF(DGINIT.GE.ZERO) THEN
         WRITE(*,*) 'Initial DG is ge zero:',DGINIT,N
         RETURN ! check initial search direction (INFO=0)
      END IF 

      BRACKT= .FALSE.         ! local variables
      STAGE1= .TRUE.

      FINIT = F               ! initial function value
      FOPT  = F
      SOPT   =0.0D0 
      DGTEST= FTOL*DGINIT     ! initial direct. derivative
      WIDTH = STPMAX - STPMIN ! width 
      WIDTH1= WIDTH/P5        
      DO J = 1, N
       WA(J) = X(J)           ! save initial parameter values
      END DO

      STX=ZERO    ! step is zero             at best step
      FX =FINIT   ! function value
      DGX=DGINIT  ! directional derivative

      STY=ZERO    ! step is zero             at other endpoint
      FY =FINIT   ! function value
      DGY=DGINIT  ! directional derivative

      SFD(1,1)=0.0D0  ! step 
      SFD(2,1)=0.0D0  ! function value
      SFD(3,1)=DGINIT ! slope

      GOTO 30

*     ------------------ new function and gradient ---------------------  
 10   IF(INFO.GT.0) GOTO 20
      INFO=0
      NFEV=NFEV+1
       
      IF(F.LT.FOPT) THEN
         FOPT=F
         SOPT=STP 
      END IF



c      WRITE(*,*) NFEV,'.th function value',F,' with step STP',STP
      STPAC=STP  ! actual step
      DG=DBDOT(N,G,S)
 
*     move up

      DO I=3,4      
       DO J=1,3
        SFD(J,I)=SFD(J,I-1)
       END DO
      END DO
      SFD(1,4)=STP
      SFD(2,4)=F-FINIT
      SFD(3,4)=DG
      DO J=1,3  
       SFD(J,5)=SPAN(SFD(J,1),SFD(J,2),SFD(J,3),SFD(J,4))
      END DO 
       
      IF(NFEV.GE.4) THEN
c         WRITE(*,*) NFEV,SFD(1,5),SFD(2,5),SFD(3,5)
         IF(SFD(1,5).LT.1.0E-5.AND.
     +      SFD(2,5).LT.1.0E-5.AND.
     +      SFD(3,5).LT.1.0E-3.AND.
     +      SOPT.EQ.0.0D0) INFO=7
         IF(SFD(1,5).LT.1.0E-5.AND.
     +      SFD(2,5).LT.1.0E-5.AND.
     +      SFD(3,5).LT.1.0E-5) INFO=7
         SFD1=SFD(1,5)  
         SFD2=SFD(2,5)  
         SFD3=SFD(3,5)  
c         IF(INFO.EQ.7) CALL RUCORR(SFD1,SFD2)        
*        INFO=7 
      END IF 


      IF(NFEV.LT.NDIMLS) THEN  ! store line search data
         DO I=1,NLS
          IF(STP.LT.XLS(I)) THEN
             DO J=NLS,I,-1
              XLS(J+1)=XLS(J)
              YLS(J+1)=YLS(J)
              TLS(J+1)=TLS(J)
              KLS(J+1)=KLS(J) 
             END DO
             ILS=I
             GOTO 15
          END IF
         END DO
         ILS=NLS+1
 15      NLS=NLS+1   
         XLS(ILS)=STP
         YLS(ILS)=F-FINIT  
         TLS(ILS)=DG
         KLS(ILS)=NFEV
      END IF



*     STP=        ! current step             at current step
*     F=          ! function value
*     DG=         ! directional derivative       
      FTEST1 = FINIT + STP*DGTEST

      IF((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     +      .OR. INFOC .EQ. 0)
     +                            INFO = 6 ! rounding error limitation 
      IF(STP .EQ. STPMAX .AND.
     +       F .LE. FTEST1 .AND. DG .LE. DGTEST)
     +                            INFO = 5 ! step at the upper bound
      IF(STP .EQ. STPMIN .AND.
     +       (F .GT. FTEST1 .OR. DG .GE. DGTEST))
     +                            INFO = 4 ! step at the lower bound
      IF(NFEV .GE. MAXF)          INFO = 3 ! max number of calls reached

      IF(INFO.EQ.3) THEN
c         WRITE(*,*) '[INFO=3] ',' DGINIT,DG ',DGINIT,DG
c         INFO=1 
      END IF

      IF(BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX)
     +                            INFO = 2 ! relative interval width limit

      IF(F.LE. FTEST1 .AND. ABS(DG) .LE. GTOL*(-DGINIT).
     *       AND.NFEV.GE.MINF)    INFO = 1 ! min conditions hold ! 

      IF(INFO.NE.0) THEN ! termination 
         DGOUT=DG
         IF(INFO.EQ.1) RETURN
c            WRITE(*,*) 'Error ',INFO,' FOPT,SOPT=',FOPT,SOPT
            DO J=1,N
             X(J)=WA(J)+SOPT*S(J)
            END DO 
         RETURN
      END IF

*     ------------------------ define STP ------------------------------

*     first stage" seek a step for which the modified function has
*     a non-positive value and non-negative derivative
 20   IF(STAGE1.AND.F.LE.FTEST1.AND.
     +   DG.GE.MIN(FTOL,GTOL)*DGINIT) STAGE1=.FALSE.

c      WRITE(*,*) 'TEST STAGE1 F,FTEST1,DG,FTOL,GTOL,DGINIT,M*DGINIT'
c      WRITE(*,*) STAGE1,F,FTEST1,DG,FTOL,GTOL,DGINIT,
c     +           MIN(FTOL,GTOL)*DGINIT
c      WRITE(*,*) 'STAGE1 FX FTEST1',STAGE1,FX,FTEST1

      IF(STAGE1.AND.F.LE.FX.AND.F.GT.FTEST1) THEN

*        use modified function to predict a step, if no step was
*        obtained for which the modified function has a non-positive
*        function value and non-negative derivative, and if   
*        a lower function value has been obtained, but with insufficient
*        decrease

         FM = F - STP*DGTEST    ! modify function value ...
         FXM = FX - STX*DGTEST
         FYM = FY - STY*DGTEST
         DGM = DG - DGTEST      ! ... and gradient values
         DGXM = DGX - DGTEST
         DGYM = DGY - DGTEST

*        update the interval and compute new step

         CALL CSTEPM(STX,FXM,DGXM,
     +               STY,FYM,DGYM,
     +               STP,FM ,DGM,
     +               BRACKT,STMIN,STMAX,INFOC)

         FX = FXM + STX*DGTEST  ! reset function value ...
         FY = FYM + STY*DGTEST
         DGX = DGXM + DGTEST    ! ... and gradient values
         DGY = DGYM + DGTEST

      ELSE

*        update the interval and compute new step

         CALL CSTEPM(STX,FX,DGX,
     +               STY,FY,DGY,
     +               STP,F ,DG,                       ! update interval
     +               BRACKT,STMIN,STMAX,INFOC)        ! and get new step
      END IF

c         WRITE(*,*) 'Update by CSTEPM: STX,FXM,DGXM,STY,FYM,DGYM,
c     +               STP,FM ,DGM,BRACKT,STMIN,STMAX,INFOC'
c         WRITE(*,*) STX,FXM,DGXM,STY,FYM,DGYM,
c     +              STP,FM ,DGM,BRACKT,STMIN,STMAX,INFOC

      IF(BRACKT) THEN ! force sufficient decrease of size of interval
         IF(ABS(STY-STX).GE.P66*WIDTH1) STP=STX+P5*(STY-STX)
         WIDTH1 = WIDTH
         WIDTH = ABS(STY-STX)
      END IF

*     ---------- prepare calculation of next function value ------------

   30 IF(BRACKT) THEN               ! interval of uncertainty defined      
         STMIN=MIN(STX,STY)             ! minimum step
         STMAX=MAX(STX,STY)             ! maximum step
      ELSE                          ! interval of uncertainty undefined 
         STMIN=STX                      ! minimum step
c         WRITE(*,*) 'vor XTRAPF STP STX XTRAPF STMAX',
c     +              STP,STX,XTRAPF,STMAX     
         STMAX=STP+XTRAPF*(STP-STX)     ! maximum step
c         WRITE(*,*) 'nach XTRAPF STMAX=',STMAX         
      END IF

      STP=MAX(STPMIN,MIN(STP,STPMAX)) ! STPMIN < STP < STPMAX

      IF((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     +   .OR. NFEV .GE. MAXF-1 .OR. INFOC .EQ. 0
     +   .OR. (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX)) THEN
*        unusual termination
c         WRITE(*,*) 'NFEV MAXF INFOC',NFEV,MAXF,INFOC
c         WRITE(*,*) 'Unusual termination, new,best STP  is',STP,STX
         STP = STX
c         WRITE(*,*) 'STX STPAC',STX,STPAC
         IF(NFEV.GE. MAXF-1.AND.STX.EQ.STPAC) THEN
            INFO=3
            STP=STPAC
            DGOUT=DG
c            WRITE(*,*) 'Special return ! DGINIT DGOUT ',
c     +                 DGINIT,DGOUT,STP
            RETURN
         END IF
      END IF

 90   DO J = 1, N          ! calculate the parameter values
       X(J)=WA(J)+STP*S(J) ! for step STP
      END DO
      INFO=-1 ! return to compute function value
      RETURN

      ENTRY MTINFO(INFO)  ! print information about errors
      IF(INFO.GE.0) WRITE(*,*) ' '
      IF(INFO.EQ.0) THEN
         WRITE(*,101) INFO,'N or STP lt 0 or step not descending'
      ELSE IF(INFO.EQ.1) THEN
         WRITE(*,101) INFO,'Linesearch convergence conditions reached'
      ELSE IF(INFO.EQ.2) THEN
         WRITE(*,101) INFO,'interval of uncertainty at lower limit'
      ELSE IF(INFO.EQ.3) THEN
         WRITE(*,101) INFO,'max nr of line search calls reached'
      ELSE IF(INFO.EQ.4) THEN
         WRITE(*,101) INFO,'step at the lower bound' 
      ELSE IF(INFO.EQ.5) THEN
         WRITE(*,101) INFO,'step at the upper bound'
      ELSE IF(INFO.EQ.6) THEN
         WRITE(*,101) INFO,'rounding error limitation'  
      ELSE IF(INFO.EQ.7) THEN
         WRITE(*,101) INFO,'function noise limit'
      END IF 
      IF(INFO.GE.0) WRITE(*,*) ' '
 101  FORMAT('     Info =',I3,3X,A)
      RETURN

*     MTLINE search initialization -------------------------------------

      ENTRY MTLDEF(FTOLE, GTOLE, XTOLE,  STMINE, STMAXE, MINFE,MAXFE)
*        example:  1.0E-4 1.0E-1 1.0E-17 1.0E-20 1.0E+20 1     40  
*     initialize "Linesearch Algorithm with Sufficient Decrease" 

      FTOL=FTOLE    ! tolerance C1  for sufficient decrease of function value
      GTOL=GTOLE    ! tolerance C2  for directional derivative condition
      XTOL=XTOLE    ! tolerance for width of intercal
      STPMIN=STMINE ! lower bound for the step
      STPMAX=STMAXE ! upper bound for the step
      MINF=MINFE    ! min number of linesearch function calls
      MAXF=MAXFE    ! max number of linesearch function calls (>1)
      IF(FTOL.LT.ZERO.OR.GTOL.LT.ZERO.OR.XTOL.LT.ZERO.OR.
     +   STMIN.LT.ZERO.OR.STPMAX.LT.STPMIN.OR.
     +   MINF.GT.MAXF.OR.MAXF.LE.1) THEN
         WRITE(*,*) 'Improper MTLDEF arguments - stop'
         STOP
      END IF 
      NFEV=-1               ! reset counter
      END


      SUBROUTINE CSTEPM(STX,FX,DX,  ! step, fvalue, der at best step
     +                  STY,FY,DY,  !                   at other endpt
     +                  STP,FP,DP,  !                   at current step 
     +                  BRACKT,        ! bracketed true/false
     +                  STPMIN,STPMAX, ! bounds for ste step
     +                  INFO)          ! returned: INFO=0 means error 
*                       INFO=1 ... 4   ! used method of interpolation
*     compute a safeguarded step for a linesearch
*     update an interval of uncertainty for the function minimum
      INTEGER INFO
      DOUBLE PRECISION STX,FX,DX,STY,FY,DY,STP,FP,DP,STPMIN,STPMAX
      LOGICAL BRACKT,BOUND
      DOUBLE PRECISION GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA
      DOUBLE PRECISION TINY,REDUCE
      PARAMETER (TINY=1.0D-10)
      PARAMETER (REDUCE=0.1D0) 
*     ...
      INFO=0
      IF((BRACKT.AND.(STP.LE.MIN(STX,STY).OR.      ! check errors
     *    STP.GE.MAX(STX,STY))).OR.
     *    DX*(STP-STX).GE.0.0.OR.STPMAX.LT.STPMIN) RETURN

      SGND=DP*SIGN(1.0D0,DX)

      IF((.NOT.FP+TINY.GT.0.0D0.AND..NOT.FP-EPS.LT.0.0D0).OR.
     +   FP.GT.1.0D+30.OR.FP.LT.(-1.0D-30)) THEN  ! nan or infty
         INFO=5 
         BOUND=.TRUE.
         BRACKT=.TRUE. 
         STPF=REDUCE*STP
         FP=1.0D10
      ELSE IF(FP.GT.FX) THEN
         INFO = 1
*        higher function value: minimum is bracketed.
*        If the cubic step is closer to STX than the quadratic step,
*        the cubic step is taken, else the average of the two steps. 
         BOUND=.TRUE.
         THETA= 3.0*(FX - FP)/(STP - STX) + DX + DP
         S=MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA=S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF(STP.LT.STX) GAMMA = -GAMMA
         P=(GAMMA-DX)+THETA
         Q=((GAMMA-DX)+GAMMA)+DP
         R=P/Q
         STPC=STX+R*(STP-STX)
         STPQ=STX+((DX/((FX-FP)/(STP-STX)+DX))/2)*(STP-STX)
         IF(ABS(STPC-STX).LT.ABS(STPQ-STX)) THEN
            STPF=STPC
         ELSE
            STPF=STPC+(STPQ-STPC)/2
         END IF
         BRACKT = .TRUE.
      ELSE IF(SGND.LT.0.0) THEN
         INFO = 2
*        lower function value and opp. derivative signs: minimum bracketed.
*        If the cubic step is closer to STX than the quadratic step,
*        the cubic step is taken, else the quadratic step.          
         BOUND = .FALSE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = ((GAMMA - DP) + GAMMA) + DX
         R = P/Q
         STPC = STP + R*(STX - STP)
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF(ABS(STPC-STP) .GT. ABS(STPQ-STP)) THEN
            STPF = STPC
         ELSE
            STPF = STPQ
         END IF
         BRACKT = .TRUE.
      ELSEIF(ABS(DP).LT.ABS(DX)) THEN
         INFO = 3
*        lower function value and same derivative signs, der. decreases:
*        The cubic step is only used, if the cubic tends to infinity in
*        the driection of the step, or if the minimum of the cuic is beyond
*        STP. Otherwise the cubic step is defined to be either STPMIN or
*        STPMAX. The quadratic step also computed and if the minimum is
*        bracketed, then the step closest to STX is taken, else the step
*        farthest away is taken.
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT(MAX(0.0D0,(THETA/S)**2 - (DX/S)*(DP/S)))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = (GAMMA + (DX - DP)) + GAMMA
         R = P/Q
         IF(R .LT. 0.0 .AND. GAMMA .NE. 0.0) THEN
            STPC = STP + R*(STX - STP)
         ELSE IF (STP .GT. STX) THEN
            STPC = STPMAX
         ELSE
            STPC = STPMIN
         END IF
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF(BRACKT) THEN
            IF(ABS(STP-STPC) .LT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
            END IF
         ELSE
            IF(ABS(STP-STPC) .GT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
            END IF
         END IF
      ELSE
         INFO = 4
*        lower function value and same derivative signs, no der. decrease:
*        If the minimum is not bracketed, the step is either STPMIN or
*        STPMAX, else the cubic step is taken. 
         BOUND = .FALSE.
         IF (BRACKT) THEN
            THETA = 3*(FP - FY)/(STY - STP) + DY + DP
            S = MAX(ABS(THETA),ABS(DY),ABS(DP))
            GAMMA = S*SQRT((THETA/S)**2 - (DY/S)*(DP/S))
            IF (STP .GT. STY) GAMMA = -GAMMA
            P = (GAMMA - DP) + THETA
            Q = ((GAMMA - DP) + GAMMA) + DY
            R = P/Q
            STPC = STP + R*(STY - STP)
            STPF = STPC
         ELSE IF (STP .GT. STX) THEN
            STPF = STPMAX
         ELSE
            STPF = STPMIN
         END IF
      END IF

*     update the interval of uncertainty (does not depend on new step) 
      IF(FP.GT.FX) THEN
         STY=STP
         FY =FP            !!!!!!!!!!!!!
         DY =DP
      ELSE
         IF(SGND.LT.0.0) THEN
            STY=STX
            FY =FX
            DY =DX
         END IF
         STX=STP
         FX =FP
         DX =DP
      END IF

      STPF=MAX(STPMIN,MIN(STPF,STPMAX))      
      STP =STPF                       ! new step
      IF(BRACKT.AND.BOUND) THEN       ! safeguard the step
         IF(STY.GT.STX) THEN
            STP= MIN(STX+0.66*(STY-STX),STP)
         ELSE
            STP= MAX(STX+0.66*(STY-STX),STP)
         END IF
      END IF
      END

