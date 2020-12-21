PROGRAM GRUA
    !Modulo
    USE SEL_MET
    USE VYM_IO
    
    IMPLICIT NONE
    
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: A
    REAL(8), DIMENSION(:), ALLOCATABLE :: B, XGAUSS, XGJ, XCROUT
    INTEGER :: BANDERAA, BANDERAB, BANDERA
    CHARACTER(*), PARAMETER :: FORMATO = '(F25.15)' !Es el mismo formato del módulo, pero por si se quiere cambiar.
    
    CALL MAT_LEER(A, BANDERAA, 'APV.txt')
    CALL VEC_LEER(B, BANDERAB, 'BPV.txt')
    !prueba
    CALL PRUEBA_PIVOTEO(A)
    
    IF (.NOT. VERIF_LECT(BANDERAA, BANDERAB)) GOTO 20
    CALL MOSTRAR_DATOS(A, B)
    
    CALL PRUEBA_METODOS_DIR(A, B, XGAUSS, XGJ, XCROUT)
    !Se guardan los resultados de la prueba
    CALL VEC_GUARDAR(XGAUSS, BANDERA, 'Resultados por Gauss.txt')
    CALL VEC_GUARDAR(XGJ, BANDERA, 'Resultados por Gauss-Jordan.txt')
    CALL VEC_GUARDAR(XCROUT, BANDERA, 'Resultados por LUCrout.txt')
    
    GOTO 10
20  PRINT *, 'Error al leer de archivo.'

10  PRINT *, 'Fin del programa.'
CONTAINS

    SUBROUTINE MOSTRAR_DATOS(A, B)
        REAL(8), INTENT(IN) :: A(:,:), B(:)
        !
        CHARACTER(*), PARAMETER :: FORMATO = '(F10.4)'
        
        PRINT *, 'Matriz inicial:'
        CALL MAT_MOSTRAR(A, FORMATO)
        
        PRINT *, 'Vector de términos independientes:'
        CALL VEC_MOSTRAR(B, FORMATO)
    
    END SUBROUTINE
    
    FUNCTION VERIF_LECT(BA, BB)
        LOGICAL :: VERIF_LECT
        INTEGER, INTENT(IN) :: BA, BB
        
        VERIF_LECT = (BA /= 1) .AND. (BB /= 1)
    END FUNCTION
    
    !Se devuelven los resultados obtenidos por los diferentes métodos directos.
    !B es un vector, asi que se devuelven resultados vectoriales.
    SUBROUTINE PRUEBA_METODOS_DIR(A, B, XGAUSS, XGJ, XCROUT)
        REAL(8), INTENT(IN) :: A(:,:), B(:)
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: XGAUSS, XGJ, XCROUT
        !
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: MATAMPGAUSS, MATGSR, MATAMPGJ
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: BMAT
        INTEGER :: N
        
        N = SIZE(B) !Cantidad de filas del sistema (cantidad de fuerzas).
        ALLOCATE(BMAT(N,1))
        BMAT(:,1) = B
        
        !Resolución por los métodos directos
        CALL MET_GAUSS(A, BMAT, MATAMPGAUSS)
        CALL SUST_REGRESIVA(MATAMPGAUSS, MATGSR) !Matriz de Gauss con Sustitucion Regresiva
        CALL MET_GAUSSJORDAN(A, BMAT, MATAMPGJ)
        CALL MET_LU_CROUT(A, B, XCROUT)
        
        !Se pasan los resultados por los métodos de Gauss y Gauss-Jordan a vectores.
        ALLOCATE(XGAUSS(N), XGJ(N))
        XGAUSS = MATGSR(:,1) !Tiene sólo los resultados, no es la matriz ampliada
        XGJ = MATAMPGJ(:,N+1)
        
        DEALLOCATE(MATAMPGAUSS, MATGSR, MATAMPGJ, BMAT)
    END SUBROUTINE
    
    SUBROUTINE PRUEBA_PIVOTEO(MAT)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT
        !
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: MATPIV
        INTEGER :: N, M
        N = SIZE(MAT,1); M = SIZE(MAT,2)
        ALLOCATE(MATPIV(N,M))
        MATPIV = MAT
        CALL SUBRUT_PIVOTEOMATAMP(MATPIV)
        PRINT *, 'Matriz original leida:'
        CALL MAT_MOSTRAR(MAT)
        PRINT *, 'Matriz pivoteada al inicio:'
        CALL MAT_MOSTRAR(MATPIV)
        DEALLOCATE(MATPIV)
    END SUBROUTINE
END PROGRAM
