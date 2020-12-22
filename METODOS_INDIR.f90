MODULE METODOS_INDIR
    USE VYM_MANIP
    USE VYM_IO
!    USE VYM_
    IMPLICIT NONE
CONTAINS    
    !----------Métodos Indirectos----------!
    !Usa A y B para aproximar X dentro de una TOL.
    !No revisa si el SEL converge o diverge, si se quiere se puede revisar antes de llamar.
    !Se por parametro viene una X inicial que sirve de ""semilla"".
    SUBROUTINE MET_JACOBI(AIN, BIN, TOL, X, XI)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: AIN
        REAL(8), DIMENSION(:), INTENT(IN) :: BIN
        REAL(8), INTENT(IN) :: TOL
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: X
        REAL(8), DIMENSION(:), INTENT(IN), OPTIONAL :: XI
        !
        REAL(8) :: ERROR
        REAL(8), ALLOCATABLE :: XANT(:), BMAT(:,:), B(:), A(:,:), MATAMP(:,:)
        
        INTEGER :: I, J, N, ITER
        INTEGER, PARAMETER :: MAXITER = 100 !Como para que no se quede en bucle infinito
        
        N = SIZE(AIN,1)
        IF (.NOT.PRESENT(XI)) THEN
            ALLOCATE(X(N))
            X = 0.
        ELSE
            X = XI
        END IF
        
!        ALLOCATE(BMAT(N,1))
!        BMAT(:,1) = BIN
        
!        CALL MATRIZAMPLIADA(AIN, BMAT, MATAMP)
        
!        PRINT *, 'MATAMP'; CALL MAT_MOSTRAR(MATAMP)
!        CALL SUBRUT_PIVOTEOMATAMP(MATAMP)
!        PRINT *, 'MATAMP'; CALL MAT_MOSTRAR(MATAMP)
        ALLOCATE(A(N,N), B(N))
!        A = MATAMP(:N,:N)
!        B = MATAMP(:,N+1)
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
            !En la teoria habla de la variacion proporcional.
            !Pero en el codigo solo hace la diferencia, asi que lo dejo asi.
            ERROR = VEC_NORMAM(X-XANT)
!            ERROR = VEC_NORMAM(VEC_RESIDUO(A, B, X))
            WRITE(*,*) 'Error', ERROR
            ITER = ITER + 1
        END DO
        IF (ITER >= MAXITER) PRINT *, 'Jacobi diverge'
        
        PRINT *, 'Iteraciones: ', ITER
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
        REAL(8), ALLOCATABLE :: XANT(:), BMAT(:,:), B(:), A(:,:), MATAMP(:,:)
        INTEGER :: I, J, N, ITER
        INTEGER, PARAMETER :: MAXITER = 100 !Como para que no se quede en bucle infinito
        
        
        N = SIZE(AIN,1)
        IF (.NOT.PRESENT(XI)) THEN
            ALLOCATE(X(N))
            X = 0.
        ELSE
            X = XI
        END IF
        
!        ALLOCATE(BMAT(N,1))
!        BMAT(:,1) = BIN
        
!        CALL MATRIZAMPLIADA(AIN, BMAT, MATAMP)
        
!        PRINT *, 'MATAMP'; CALL MAT_MOSTRAR(MATAMP)
!        CALL SUBRUT_PIVOTEOMATAMP(MATAMP)
!        CALL MAT_INTERFILAS(MATAMP, 3, 4)
!        PRINT *, 'MATAMP'; CALL MAT_MOSTRAR(MATAMP)
        ALLOCATE(A(N,N), B(N))
!        A = MATAMP(:N,:N)
!        B = MATAMP(:,N+1)
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
!            !En la teoria habla de la variacion proporcional.
!            !Pero en el codigo solo hace la diferencia, asi que lo dejo asi.
!            ERROR = VEC_NORMAM(X-XANT)
            ERROR = VEC_NORMAM(VEC_RESIDUO(A, B, X))
            WRITE(*,*) 'Error', ERROR
            ITER = ITER + 1
        END DO
        IF (ITER >= MAXITER) PRINT *, 'Gauss-Seidel diverge'
        WRITE(*,*) 'Iteraciones: ', ITER
        WRITE(*,*) 'Error: ', ERROR
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
    
    !---no entiendo bien el metodo y no se por que no anda---!
    !Es un método que hace sobre Gauss-Seidel.
    !Si W (omega) > 2 diverge.
    !Si W < 1 es sub-relajacion.
    !Si 1 < w < 2 es sobre-relajacion (el OR de SOR)
    !La idea sería entonces que W sobre-relaje.
    SUBROUTINE MET_SOR(A, B, TOL, W, X, XI)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A
        REAL(8), DIMENSION(:), INTENT(IN) :: B
        REAL(8), INTENT(IN) :: TOL, W
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: X
        REAL(8), DIMENSION(:), INTENT(IN), OPTIONAL :: XI
        !
        REAL(8), DIMENSION(:), ALLOCATABLE :: R, XANT
        REAL(8) :: ERROR
        INTEGER :: I, N, ITER
        INTEGER, PARAMETER :: MAXITER = 1000 !Como para que no se quede en bucle infinito
        
        N = SIZE(A,1)
        IF (.NOT.PRESENT(XI)) THEN
            ALLOCATE(X(N))
            X = 0.
        ELSE
            X = XI
        END IF
        ALLOCATE(R(N), XANT(N))
        ERROR = 2.*TOL !Valor imposible
        ITER = 0; ERROR = 2.*TOL !Valor imposible
        DO WHILE(ITER < MAXITER .AND. ERROR >= TOL)
            I = 0
            XANT = X
            DO I = 1, N
                !Por alguna razon, en la clase 5 el residuo es ax-b, pero en la clase 6 es b-ax.
                !Voy a seguir la clase 6 asi que lo voy a poner negativo.
                R = VEC_RESIDUO(A, X, B)
                X(I) = XANT(I) + W*R(I)/A(I,I)
                CALL VEC_MOSTRAR(X)
                READ(*,*)
            END DO
            ITER = ITER + 1
        END DO
        DEALLOCATE(R)
    END SUBROUTINE
END MODULE
