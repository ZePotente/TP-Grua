MODULE METODOS_DIR
    USE VYM_MANIP
    IMPLICIT NONE
CONTAINS    
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
        !Paso 2 sustitucion hacia arriba ("regresiva").
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
    !Realida pivoteo parcial por filas
    SUBROUTINE MET_GAUSS(A, B, GAUSS)
        REAL(8), ALLOCATABLE, DIMENSION(:,:) :: GAUSS
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A, B
        !
        INTEGER :: J, FILA, ORDEN
        
        ORDEN = SIZE(A,1)
        
        CALL MATRIZAMPLIADA(A, B, GAUSS)
        DO J = 1, ORDEN
            CALL PIVOTEOMATAMP(GAUSS, J)
            DO FILA = J+1, ORDEN
                GAUSS(FILA,J+1:) = GAUSS(FILA,J+1:) - GAUSS(J,J+1:)*GAUSS(FILA,J)/GAUSS(J,J)
                GAUSS(FILA,J) = 0.
            END DO
        END DO
    END SUBROUTINE
    
    !Devuelve una matriz ampliada en GJ, con los resultados correspondientes a B.
    !Realiza pivoteo parcial por filas
    SUBROUTINE MET_GAUSSJORDAN(A, B, GJ)
        REAL(8), ALLOCATABLE, DIMENSION(:,:) :: GJ
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A, B
        !
        INTEGER :: J, FILA, ORDEN
        
        ORDEN = SIZE(A,1)
        
        CALL MATRIZAMPLIADA(A, B, GJ)
        
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
        
        DEALLOCATE(C)
    END SUBROUTINE
    
    !A es la combinacion de L y U1 de MAT.
    SUBROUTINE LU_REDUCCION(MAT, RED)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: RED
        !
        INTEGER :: I, FILA, COL, K, N
        N = SIZE(MAT,1)
        
        RED = MAT
        RED(1,2:) = RED(1,2:)/RED(1,1)
        
        DO I = 2, N
            !Calcula la columna i
            DO FILA = I, N
                DO K = 1, I-1
                    RED(FILA,I) = RED(FILA,I) - RED(FILA,K)*RED(K,I)
                END DO
            END DO
            !Calcula fila i
            DO COL = I+1, N
                DO K = 1, I-1
                    RED(I,COL) = RED(I,COL) - RED(I,K)*RED(K,COL)
                END DO
                RED(I,COL) = RED(I,COL)/RED(I,I)
            END DO
        END DO
    END SUBROUTINE
END MODULE
