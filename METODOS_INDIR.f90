MODULE METODOS_INDIR
    USE VYM_MANIP
    USE VYM_IO
    IMPLICIT NONE
CONTAINS    
    !----------Métodos Indirectos----------!
    !Usa A y B para aproximar X dentro de una TOL.
    !No revisa si el SEL converge o diverge, si se quiere se puede revisar antes de llamar.
    !Por parametro puede venir un vector XI (X inicial) que sirve de "semilla".
    SUBROUTINE MET_JACOBI(AIN, BIN, TOL, X, XI)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: AIN
        REAL(8), DIMENSION(:), INTENT(IN) :: BIN
        REAL(8), INTENT(IN) :: TOL
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: X
        REAL(8), DIMENSION(:), INTENT(IN), OPTIONAL :: XI
        !
        REAL(8) :: ERROR
        REAL(8), ALLOCATABLE :: XANT(:), B(:), A(:,:)
        
        INTEGER :: I, J, N, ITER
        INTEGER, PARAMETER :: MAXITER = 100 !Como para que no se quede en bucle infinito
        
        N = SIZE(AIN,1)
        IF (.NOT.PRESENT(XI)) THEN
            ALLOCATE(X(N))
            X = 0.
        ELSE
            X = XI
        END IF
        
        ALLOCATE(A(N,N), B(N))
        A = AIN
        B = BIN
        ALLOCATE(XANT(N))
        ITER = 0; ERROR = 2.*TOL !Valor imposible
        DO WHILE(ITER < MAXITER .AND. ERROR >= TOL)
            XANT = X
            DO I = 1, N
                X(I) = B(I)
                DO J = 1, I-1
                    X(I) = X(I) - A(I,J)*XANT(J)
                END DO
                
                DO J = I+1, N
                    X(I) = X(I) - A(I,J)*XANT(J)
                END DO
                X(I) = X(I)/A(I,I)
            END DO
            ERROR = VEC_NORMAM(X-XANT)
!            ERROR = VEC_NORMAM(VEC_RESIDUO(A, B, X))
            WRITE(*,*) 'Error', ERROR
            ITER = ITER + 1
        END DO
        IF (ITER >= MAXITER) PRINT *, 'Jacobi diverge'
        DEALLOCATE(XANT)
    END SUBROUTINE
    
    SUBROUTINE MET_GS(AIN, BIN, TOL, X, XI)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: AIN
        REAL(8), DIMENSION(:), INTENT(IN) :: BIN
        REAL(8), INTENT(IN) :: TOL
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: X
        REAL(8), DIMENSION(:), INTENT(IN), OPTIONAL :: XI
        !
        REAL(8) :: ERROR
        REAL(8), ALLOCATABLE :: XANT(:), B(:), A(:,:)
        INTEGER :: I, J, N, ITER
        INTEGER, PARAMETER :: MAXITER = 100 !Como para que no se quede en bucle infinito
        
        
        N = SIZE(AIN,1)
        IF (.NOT.PRESENT(XI)) THEN
            ALLOCATE(X(N))
            X = 0.
        ELSE
            X = XI
        END IF
        
        ALLOCATE(A(N,N), B(N))
        A = AIN
        B = BIN
        PRINT *,'A'
        CALL MAT_MOSTRAR(A)
        PRINT *,'B'
        CALL VEC_MOSTRAR(B)
!        ALLOCATE(XANT(N))
        ITER = 0; ERROR = 2.*TOL !Valor imposible
        DO WHILE(ITER < MAXITER .AND. ERROR >= TOL)
!            XANT = X
            DO I = 1, N
                X(I) = B(I)
        
                DO J = 1, I-1
                    X(I) = X(I) - A(I,J)*X(J)
                END DO
                
                DO J = I+1, N
                    X(I) = X(I) - A(I,J)*X(J)
                END DO
        
                X(I) = X(I)/A(I,I)
            END DO
!            ERROR = VEC_NORMAM(X-XANT)
            ERROR = VEC_NORMAM(VEC_RESIDUO(A, B, X))
            WRITE(*,*) 'Error', ERROR
            ITER = ITER + 1
        END DO
        IF (ITER >= MAXITER) PRINT *, 'Gauss-Seidel diverge'
!        DEALLOCATE(XANT)
    END SUBROUTINE
    
    FUNCTION VEC_RESIDUO(A, B, X)
        REAL(8), DIMENSION(:), ALLOCATABLE :: VEC_RESIDUO
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A
        REAL(8), DIMENSION(:), INTENT(IN) :: B, X
        !
        ALLOCATE(VEC_RESIDUO(SIZE(X)))
        VEC_RESIDUO = MATMUL(A,X) - B
    END FUNCTION
END MODULE
