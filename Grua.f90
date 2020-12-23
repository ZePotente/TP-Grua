PROGRAM GRUA
    !Modulo
    USE METODOS_DIR
    USE METODOS_INDIR
    USE VYM_IO
    
    IMPLICIT NONE
    
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: A, ACROUT
    REAL(8), DIMENSION(:), ALLOCATABLE :: B, BCROUT,           BMEDIO, BCERCA, BMAX
    REAL(8), DIMENSION(:), ALLOCATABLE :: XGAUSS, XGJ, XCROUT, XMEDIO, XCERCA, XMAX, XGS
    INTEGER :: BANDERAA, BANDERAB, BANDERAC, BANDERAD !Banderas reutilizables.
    CHARACTER(*), PARAMETER :: FORMATO = '(F10.4)' !Formato global.
    
    !Prueba inicial
    PRINT *, 'Lectura de los datos iniciales.'
    CALL MAT_LEER(A, BANDERAA, 'A.txt')
    CALL VEC_LEER(B, BANDERAB, 'B.txt')
    !Inicialmente era el sistema pivoteado manualmente para Crout, pero luego se utilizó también para probar los métodos indirectos 
    !(por eso el nombre de las variables).
    CALL MAT_LEER(ACROUT, BANDERAC, 'ACrout.txt')
    CALL VEC_LEER(BCROUT, BANDERAD, 'BCrout.txt')
    
    IF (.NOT. VERIF_ARCH(BANDERAA, BANDERAB, BANDERAC) .OR. BANDERAD == 1) GOTO 20
    PRINT *, 'Datos leidos correctamente.'
    CALL MOSTRAR_DATOS(A, B, ACROUT, BCROUT)
    
    CALL MET_JACOBI(ACROUT, BCROUT, 1D-10, XGS)
    PRINT *, 'Resultado por Jacobi con matriz pivoteada manualmente:'
    CALL VEC_MOSTRAR(XGS)
    CALL MET_GS(ACROUT, BCROUT, 1D-10, XGS)
    PRINT *, 'Resultado por Gauss-Seidel'
    CALL VEC_MOSTRAR(XGS)
    
    CALL CONFIRMAR()
    CALL SYSTEM("Clear")
    
    PRINT *, 'Prueba de métodos en la posición inicial.'
    CALL PRUEBA_METODOS_DIR(A, B, ACROUT, BCROUT, XGAUSS, XGJ, XCROUT)
    PRINT *, 'Métodos finalizados.'
    
    PRINT *, 'Guardando resultados en archivo...'
    CALL VEC_GUARDAR(XGAUSS, BANDERAA, 'Resultados por Gauss.txt')
    CALL VEC_GUARDAR(XGJ, BANDERAB, 'Resultados por Gauss-Jordan.txt')
    CALL VEC_GUARDAR(XCROUT, BANDERAC, 'Resultados por LUCrout.txt')
    
    IF (.NOT. VERIF_ARCH(BANDERAA, BANDERAB, BANDERAC)) GOTO 20 
    PRINT *, 'Resultados guardados correctamente.'
    PRINT *, 'Prueba inicial terminada.'
    
    CALL CONFIRMAR()
    CALL SYSTEM("Clear")
    
    !Prueba en distintas posiciones y carga con el Método de Gauss-Jordan.
    PRINT *, 'Resolución en diferentes condiciones mediante el método de Gauss-Jordan.'
    PRINT *, 'Leyendo términos independientes...'
    CALL VEC_LEER(BMEDIO, BANDERAA, 'B - Posicion intermedia.txt')
    CALL VEC_LEER(BCERCA, BANDERAB, 'B - Posicion cercana.txt')
    CALL VEC_LEER(BMAX,   BANDERAC, 'B - Maximo.txt')
    
    IF (.NOT. VERIF_ARCH(BANDERAA, BANDERAB, BANDERAC)) GOTO 20 
    PRINT *, 'Datos leidos correctamente.'
    CALL MOSTRAR_TERMIND(BMEDIO, BCERCA, BMAX)
    
    PRINT *, 'Resolviendo sistemas mediante Gauss-Jordan.'
    CALL RESOLUCION_GJ(A, BMEDIO, BCERCA, BMAX, XMEDIO, XCERCA, XMAX)
    PRINT *, 'Métodos finalizados.'
    
    PRINT *, 'Guardando resultados en archivo...'
    CALL VEC_GUARDAR(XMEDIO, BANDERAA, 'Resultados en el centro.txt')
    CALL VEC_GUARDAR(XCERCA, BANDERAB, 'Resultados en posicion cercana.txt')
    CALL VEC_GUARDAR(XMAX, BANDERAC, 'Resultados del peso maximo.txt')
    
    IF (.NOT. VERIF_ARCH(BANDERAA, BANDERAB, BANDERAC)) GOTO 20 
    PRINT *, 'Resultados guardados correctamente.'
    
    GOTO 10
20  PRINT *, 'Error al tratar con algún archivo.'

10  PRINT *, 'Fin del programa.'
CONTAINS
    SUBROUTINE CONFIRMAR()
        PRINT *, 'Presione enter para continuar.'
        READ(*,*)
    END SUBROUTINE
    
    SUBROUTINE MOSTRAR_DATOS(A, B, ACROUT, BCROUT)
        REAL(8), INTENT(IN) :: A(:,:), ACROUT(:,:), B(:), BCROUT(:)
        !
        
        PRINT *, 'Matriz inicial:'
        CALL MAT_MOSTRAR(A, FORMATO)
        
        PRINT *, 'Vector de términos independientes:'
        CALL VEC_MOSTRAR(B, FORMATO)
        
        PRINT *, 'Matriz inicial para método de Crout:'
        CALL MAT_MOSTRAR(ACROUT, FORMATO)
        
        PRINT *, 'Vector de términos independientes para método de Crout:'
        CALL VEC_MOSTRAR(BCROUT, FORMATO)
    
    END SUBROUTINE
    
    SUBROUTINE MOSTRAR_TERMIND(BMEDIO, BCERCA, BMAX)
        REAL(8), DIMENSION(:) :: BMEDIO, BCERCA, BMAX
        !
        
        PRINT *, 'Terminos independientes en el medio:'
        CALL VEC_MOSTRAR(BMEDIO, FORMATO)
        
        PRINT *, 'Terminos independientes en el punto más cercano:'
        CALL VEC_MOSTRAR(BCERCA, FORMATO)
        
        PRINT *, 'Terminos independientes con carga máxima:'
        CALL VEC_MOSTRAR(BMAX, FORMATO)
    END SUBROUTINE
    
    FUNCTION VERIF_ARCH(BA, BB, BC)
        LOGICAL :: VERIF_ARCH
        INTEGER, INTENT(IN) :: BA, BB, BC
        
        VERIF_ARCH = (BA /= 1) .AND. (BB /= 1) .AND. (BC /= 1)
    END FUNCTION
    
    !Se devuelven los resultados obtenidos por los diferentes métodos directos.
    !B es un vector, asi que se devuelven resultados vectoriales.
    SUBROUTINE PRUEBA_METODOS_DIR(A, B, ACROUT, BCROUT, XGAUSS, XGJ, XCROUT)
        REAL(8), INTENT(IN) :: A(:,:), B(:), ACROUT(:,:), BCROUT(:)
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
        CALL MET_LU_CROUT(ACROUT, BCROUT, XCROUT)
        
        !Se pasan los resultados por los métodos de Gauss y Gauss-Jordan a vectores.
        ALLOCATE(XGAUSS(N), XGJ(N))
        XGAUSS = MATGSR(:,1) !Tiene sólo los resultados, no es la matriz ampliada
        XGJ = MATAMPGJ(:,N+1)
        
        DEALLOCATE(MATAMPGAUSS, MATGSR, MATAMPGJ, BMAT)
    END SUBROUTINE
    
    SUBROUTINE RESOLUCION_GJ(A, BMEDIO, BCERCA, BMAX, XMEDIO, XCERCA, XMAX)
        REAL(8), DIMENSION(:,:) :: A
        REAL(8), DIMENSION(:) :: BMEDIO, BCERCA, BMAX
        REAL(8), DIMENSION(:), ALLOCATABLE :: XMEDIO, XCERCA, XMAX
        !
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: MATAMPMEDIO, MATAMPCERCA, MATAMPMAX, BMATMEDIO, BMATCERCA, BMATMAX
        INTEGER :: N
        
        N = SIZE(B) !Cantidad de filas del sistema (cantidad de fuerzas).
        ALLOCATE(BMATMEDIO(N,1), BMATCERCA(N,1), BMATMAX(N,1))
        BMATMEDIO(:,1) = BMEDIO
        BMATCERCA(:,1) = BCERCA
        BMATMAX(:,1) = BMAX
        
        CALL MET_GAUSSJORDAN(A, BMATMEDIO, MATAMPMEDIO)
        
        CALL MET_GAUSSJORDAN(A, BMATCERCA, MATAMPCERCA)
        
        CALL MET_GAUSSJORDAN(A, BMATMAX, MATAMPMAX)
        
        !Se pasan los resultados a vectores.
        ALLOCATE(XMEDIO(N), XCERCA(N), XMAX(N))
        XMEDIO = MATAMPMEDIO(:,N+1)
        XCERCA = MATAMPCERCA(:,N+1)
        XMAX = MATAMPMAX(:,N+1)
        
        DEALLOCATE(MATAMPMEDIO, MATAMPCERCA, MATAMPMAX, BMATMEDIO, BMATCERCA, BMATMAX)
    END SUBROUTINE
END PROGRAM
