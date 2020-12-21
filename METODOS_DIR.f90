MODULE METODOS_DIR
    IMPLICIT NONE
CONTAINS
    FUNCTION MATRIZAMPLIADA(A, B)
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: MATRIZAMPLIADA
        REAL(8), DIMENSION(:,:) :: A, B
        
        INTEGER :: N, M, Q
        
        N = SIZE(A,1); M = SIZE(A,2);
        Q = SIZE(B,2);
        
        ALLOCATE(MATRIZAMPLIADA(N,M+Q))
        
        MATRIZAMPLIADA(:,:M) = A
        MATRIZAMPLIADA(:,M+1:) = B
    END FUNCTION
    
    SUBROUTINE PIVOTEOMATAMP(MAT, J, B)
        REAL(8), DIMENSION(:,:), INTENT(INOUT) :: MAT
        REAL(8), DIMENSION(:), OPTIONAL :: B
        INTEGER, INTENT(IN) :: J
        !
        REAL(8) :: MAXIMO
        INTEGER :: I, INUEVO, N
        N = SIZE(MAT,1)
        MAXIMO = 0.
        INUEVO = J
        DO I = J+1, N
            IF (ABS(MAT(I,J))>MAXIMO) THEN
                MAXIMO = ABS(MAT(I,J))
                INUEVO = I
            END IF
        END DO
        
        IF (INUEVO /= J) THEN
!            PRINT *, 'INTERCAMBIO'
!            PRINT *, 'ANTES:'
!            CALL MAT_MOSTRAR(MAT, '(F8.4)')
            CALL MAT_INTERFILAS(MAT, INUEVO, J)
            IF (PRESENT(B)) CALL VEC_INTERELEM(B, INUEVO, J)
!            WRITE(*,*) 'Se intercambio la fila ', J, ' por la fila ', INUEVO
!            PRINT *, 'DESPUES:'
!            CALL MAT_MOSTRAR(MAT, '(F8.4)')
        END IF
    END SUBROUTINE
    
    SUBROUTINE SUST_REGRESIVA(MATAMP, X)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MATAMP
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: X
        !
        INTEGER :: I, K, N, M
        REAL(8), DIMENSION(:), ALLOCATABLE :: SUMA
        
        N = SIZE(MATAMP,1); M = SIZE(MATAMP,2)
        ALLOCATE(X(N,M-N), SUMA(N))
        X = MATAMP(:,N+1:)
        !Paso 1 dividir la ultima fila para tener el resultado del final.
        X(N,:) = X(N,:)/MATAMP(N,N)
        !Paso 2 sustitucion hacia arriba.
        DO I = N-1, 1, -1
            SUMA(:) = 0.
            DO K = I+1, N
                SUMA = SUMA + MATAMP(I,K)*X(K,:)
            END DO
            X(I,:) = (MATAMP(I,N+1:) - SUMA(:))/MATAMP(I,I)
        END DO
        
        DEALLOCATE(SUMA)
    END SUBROUTINE
    
    !Devuelve una matriz ampliada en GAUSS, con los resultados correspondientes a B.
    !PUEDE REALIZAR PIVOTEO.
    SUBROUTINE MET_GAUSS(A, B, GAUSS)
        REAL(8), ALLOCATABLE, DIMENSION(:,:) :: GAUSS
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A, B
        !
        INTEGER :: J, FILA, ORDEN
        
        ORDEN = SIZE(A,1)
        
        GAUSS = MATRIZAMPLIADA(A, B)
        DO J = 1, ORDEN
            CALL PIVOTEOMATAMP(GAUSS, J)
            DO FILA = J+1, ORDEN
                GAUSS(FILA,J+1:) = GAUSS(FILA,J+1:) - GAUSS(J,J+1:)*GAUSS(FILA,J)/GAUSS(J,J)
                GAUSS(FILA,J) = 0.
            END DO
        END DO
    END SUBROUTINE
    
    !Devuelve una matriz ampliada en GJ, con los resultados correspondientes a B.
    !PUEDE REALIZAR PIVOTEO.
    SUBROUTINE MET_GAUSSJORDAN(A, B, GJ)
        REAL(8), ALLOCATABLE, DIMENSION(:,:) :: GJ
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A, B
        !
        INTEGER :: J, FILA, ORDEN
        
        ORDEN = SIZE(A,1)
        
        GJ = MATRIZAMPLIADA(A, B)
        
        DO J = 1, ORDEN
            CALL PIVOTEOMATAMP(GJ, J)
            DO FILA = 1, J-1
                GJ(FILA,J+1:) = GJ(FILA,J+1:) - GJ(J,J+1:)*GJ(FILA,J)/GJ(J,J)
                GJ(FILA,J) = 0.
            END DO
            
            DO FILA = J+1, ORDEN
                GJ(FILA,J+1:) = GJ(FILA,J+1:) - GJ(J,J+1:)*GJ(FILA,J)/GJ(J,J)
                GJ(FILA,J) = 0.
            END DO
            !poner 1 en la diagonal, dividiendo toda la fila por el elemento de la diagonal.
            !equivalente a hacer sustitucion regresiva al final.
            GJ(J,J+1:) = GJ(J,J+1:)/GJ(J,J)
            GJ(J,J) = 1.
        END DO
    END SUBROUTINE
    
    SUBROUTINE MET_LU_CROUT(A, B, CROUT)
        REAL(8), INTENT(IN) :: A(:,:), B(:)
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: CROUT
        !
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: RED
        REAL(8), DIMENSION(:), ALLOCATABLE :: C
        INTEGER :: I, K, N
        
        CALL LU_REDUCCION(A,RED)
        N = SIZE(A,1)
        ALLOCATE(C(N), CROUT(N))
        !Copiado
        !Sustitucion progresiva para calcular C
        C(1) = B(1)/RED(1,1)
        DO I = 2, N
            C(I) = B(I)
            DO K = 1, I-1
                C(I) = C(I) - RED(I,K)*C(K)
            END DO
            C(I) = C(I)/RED(I,I)
        END DO
        !Sustitucion regresiva para calcular X
        CROUT(N) = C(N)
        DO I = N-1, 1, -1
            CROUT(I) = C(I)
            DO K = I+1, N
                CROUT(I) = CROUT(I) - RED(I,K)*CROUT(K)
            END DO
        END DO
        !Fin copiado
        
        DEALLOCATE(C)
    END SUBROUTINE
    
    !A es la combinacion de L y U1 de MAT.
    SUBROUTINE LU_REDUCCION(MAT, RED)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: RED
        !
        INTEGER :: I, FILA, COL, K, N
        N = SIZE(MAT,1)
        !copio
        RED = MAT
        RED(1,2:) = RED(1,2:)/RED(1,1)
        
        DO I = 2, N
            !"Calcula la columna i"
            DO FILA = I, N
                DO K = 1, I-1
                    RED(FILA,I) = RED(FILA,I) - RED(FILA,K)*RED(K,I)
                END DO
            END DO
            !"Calcula fila i"
            DO COL = I+1, N
                DO K = 1, I-1
                    RED(I,COL) = RED(I,COL) - RED(I,K)*RED(K,COL)
                END DO
                RED(I,COL) = RED(I,COL)/RED(I,I)
            END DO
        END DO
    END SUBROUTINE
END MODULE
