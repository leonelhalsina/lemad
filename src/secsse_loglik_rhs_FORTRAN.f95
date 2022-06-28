! Helper function: 
! fill vec with N elements from parms, starting at position ii
!==========================================================================

      SUBROUTINE secsse_fill1d (vec, DIMP, parms, II)
      IMPLICIT NONE
      INTEGER DIMP, II, I
      DOUBLE PRECISION vec(DIMP), parms(*)
     
        DO I = 1, DIMP
          II = II + 1
          vec(I) = parms(II)
        ENDDO
        
      END SUBROUTINE secsse_fill1d

!==========================================================================
! module with declarations
!==========================================================================

      MODULE secsse_dimmod

      ! length of the vector -  decided in R-code
      INTEGER  :: N
      
      ! 1 parameter vectors with unknown length
      DOUBLE PRECISION, ALLOCATABLE  :: P(:)                                
      
      ! Boolean: will become TRUE if the parameters have a value
      LOGICAL :: initialised = .FALSE.

      END MODULE secsse_dimmod

!==========================================================================
!==========================================================================
! Initialisation: name of this function as passed by "initfunc" argument
! Sets the fixed parameter vector, and allocates memory
!==========================================================================
!==========================================================================

      SUBROUTINE secsse_initmod (steadyparms)
      USE secsse_dimmod 

      IMPLICIT NONE
      EXTERNAL steadyparms

      INTEGER, PARAMETER :: nparsmall = 1  ! constant-length parameters
      
      DOUBLE PRECISION parms(nparsmall)
      COMMON /XCBPar/parms                 ! common block 

! Set the fixed parameters obtained from R
      CALL steadyparms(nparsmall, parms)

! first parameter has the length of the vector       
      N = INT(parms(1) + 1e-6)  

! Allocate variable size arrays (state variables, derivatives and parameters)

      IF (ALLOCATED(P)) DEALLOCATE(P)  
      ALLOCATE(P((N**3)/8 + (N**2)/4 + N/2))

      initialised = .FALSE.
       
      END SUBROUTINE secsse_initmod
      
!==========================================================================
!==========================================================================
! Dynamic routine: name of this function as passed by "func" argument
! variable parameter values are passed via yout
!==========================================================================
!==========================================================================
       
      SUBROUTINE secsse_runmod (neq, t, Conc, dConc, yout, ip)
      USE secsse_dimmod
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i, ii
      DOUBLE PRECISION  :: t, Conc(N), dConc(N), yout(*), FF1, FF2


! parameters - named here
      DOUBLE PRECISION rn
      COMMON /XCBPar/rn


! local variables
      CHARACTER(len=80) msg

!............................ statements ..................................

      IF (.NOT. Initialised) THEN
        ! check memory allocated to output variables
        IF (ip(1) < 1) CALL rexit("nout not large enough") 

        ! save parameter values in yout
        ii = ip(1)   ! Start of parameter values
        CALL secsse_fill1d(P, (N**2) / 4 + N, yout, ii) ! ii is updated in Fill1D
        Initialised = .TRUE.      ! to prevent from initialising more than once
      ENDIF

! dynamics

!  R code
!  dE<- mus-(lambdas + mus + Q %*% (rep(1,d))) * Es + lambdas * Es * Es + ( Q %*% Es )

      DO I = 1, N/2
        FF1 = P(N/2 + I) - (P(I) + P(N/2 + I)) * Conc(I)
        FF2 = P(I) * Conc(I) * Conc(I)
        dConc(I) = FF1 + FF2
        DO II = 1, N/2
           FF1 = Conc(II) - Conc(I)
           FF2 = P(N + (I - 1) * N/2 + II) * FF1
           dConc(I) = dConc(I) + FF2
        ENDDO
      ENDDO

!  R code
!  dD<- -(lambdas + mus + Q %*% (rep(1,d))) * Ds + 2 * lambdas * Es * Ds + ( Q %*% Ds )

      DO I = 1, N/2
        FF1 = (-P(I) - P(N/2 + I) + 2 * P(I) * Conc(I)) * Conc(N/2 + I)
        dConc(N/2 + I) = FF1
        DO II = 1, N/2
           FF1 = Conc(N/2 + II) - Conc(N/2 + I)
           FF2 = P(N + (I - 1) * N/2 + II) * FF1
           dConc(N/2 + I) = dConc(N/2 + I) + FF2
        ENDDO
      ENDDO

      END SUBROUTINE secsse_runmod
      
!==========================================================================
       
      SUBROUTINE cla_secsse_runmod (neq, t, Conc, dConc, yout, ip)
      USE secsse_dimmod
      IMPLICIT NONE

!......................... declaration section.............................

      INTEGER           :: neq, ip(*), i, ii, iii, arraydim, matdim
      DOUBLE PRECISION  :: t, Conc(N), dConc(N), yout(*), FF1, FF2
      REAL(16)          :: lambdas(N/2,N/2,N/2), mus(N/2), Qs(N/2,N/2)
      REAL(16)          :: lamEE(N/2,N/2,N/2), lamDE(N/2,N/2,N/2)

! parameters - named here
      DOUBLE PRECISION rn
      COMMON /XCBPar/rn


! local variables
      CHARACTER(len=80) msg

!............................ statements ..................................

      IF (.NOT. Initialised) THEN
        ! check memory allocated to output variables
        IF (ip(1) < 1) CALL rexit("nout not large enough") 

        ! save parameter values in yout
        ii = ip(1)   ! Start of parameter values
        CALL secsse_fill1d(P, (N**3)/8 + N/2 + (N**2)/4, yout, ii)   ! ii is updated in Fill1D
        Initialised = .TRUE.          ! to prevent from initialising more than once
      ENDIF

!lambdas = P(1 ... (N**3)/8)
!mus = P((N**3)/8 + 1 ... (N**3)/8 + N/2)
!Qs = P((N**3)/8 + N/2 + 1 ... (N**3)/8 + N/2 + (N**2)/4)

      arraydim = (N**3)/8
      matdim = (N**2)/4
      DO I = 1, N/2
        mus(I) = P(arraydim + I)
        DO II = 1, N/2
           Qs(I,II) = P(arraydim + N/2 + (II - 1) * N/2 + I)
           DO III = 1,N/2
              lambdas(I,II,III) = P((I - 1) * matdim + (III - 1) * N/2 + II)
              FF1 = Conc(II) * Conc(III)
              lamEE(I,II,III) = lambdas(I,II,III) * FF1
              FF1 = Conc(N/2 + II) * Conc(III)
              FF2 = Conc(N/2 + III) * Conc(II)
              lamDE(I,II,III) = lambdas(I,II,III) * (FF1 + FF2)
           ENDDO
        ENDDO
        Qs(I,I) = 0
      ENDDO


! dynamics

!  R code
!  all_states<-cbind(Ds,Es)
!  a<-cbind(all_states[,2],all_states[,1])
!  b<-t(all_states)
!  cross_D_E <- a%*%b
!  dE <-  -((unlist(lapply(lambdas,sum))) + mus + Q %*% (rep(1,d))) * Es + ( Q %*% Es ) + mus + unlist(lapply(lapply(lambdas,"*",Es%*%t(Es)),sum))

      DO I = 1, N/2
        FF1 = mus(I) - (SUM(lambdas(I,:,:)) + mus(I)) * Conc(I)
        dConc(I) = FF1 + SUM(lamEE(I,:,:))
        DO II = 1, N/2
           FF1 = Conc(II) - Conc(I)
           dConc(I) = dConc(I) + Qs(I, II) * FF1
        ENDDO
      ENDDO

!  R code
!  dD <-  -((unlist(lapply(lambdas,sum))) + mus + Q %*% (rep(1,d))) * Ds +( Q %*% Ds ) + unlist(lapply(lapply(lambdas,"*",cross_D_E),sum))

      DO I = 1, N/2
       FF1 = (-SUM(lambdas(I,:,:)) - mus(I)) * Conc(N/2 + I) 
       dConc(N/2 + I) = FF1 + SUM(lamDE(I,:,:))
        DO II = 1, N/2
           FF1 = Conc(N/2 + II) - Conc(N/2 + I)
           dConc(N/2 + I) = dConc(N/2 + I) + Qs(I, II) * FF1
        ENDDO
      ENDDO

      END SUBROUTINE cla_secsse_runmod
      
