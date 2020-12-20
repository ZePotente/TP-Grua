PROGRAM SEL
    !Modulos
    USE VYM_IO
    USE VYM_CALCULOS
    USE SEL_MET
!    USE Test_SEL
    
    IMPLICIT NONE
    
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: A, BMAT, XMAT, XMAT2, MATAMP, AINV, XG
    REAL(8), DIMENSION(:), ALLOCATABLE :: B, x
    REAL(8) :: TOLERANCIA
    INTEGER :: BANDERA
    
    CALL MAT_LEER(A, BANDERA, 'Aacomodado.txt')
    IF (BANDERA /= 0) PRINT *, 'Error al leer A'
!    CALL MAT_LEER(BMAT, BANDERA, '1)b)B.txt')
!    IF (BANDERA /= 0) PRINT *, 'Error al leer BMAT'
    CALL VEC_LEER(B, BANDERA, 'B.txt')
    IF (BANDERA /= 0) PRINT *, 'Error al leer B'
    
    ALLOCATE(BMAT(SIZE(B),1))
    BMAT(:,1) = B
    
    PRINT *, 'A'
    CALL MAT_MOSTRAR(A)
    PRINT *, 'B'
    CALL VEC_MOSTRAR(B)
    PRINT *, 'BMAT'
    CALL MAT_MOSTRAR(BMAT)
    CALL MET_GAUSSJORDAN(A, BMAT, XMAT)
    CALL MAT_MOSTRAR(XMAT)
    CALL MAT_GUARDAR(XMAT, BANDERA, 'GJ.txt')
    DEALLOCATE(XMAT)
    
    CALL MET_GAUSS(A, BMAT, MATAMP)
    CALL MAT_MOSTRAR(MATAMP, '(F10.5)')
    CALL SUST_REGRESIVA(MATAMP, XG)
    CALL MAT_MOSTRAR(XG)
    CALL MAT_GUARDAR(XG, BANDERA, 'G.txt')
    CALL MAT_GUARDAR(MATAMP, BANDERA, 'GDIAG.txt')
    
!    PRINT *, 'RESIDUO'
!    CALL VEC_MOSTRAR(VEC_RESIDUO(A, B, X(:,1)))
    
!    CALL MAT_INVERSA(A, AINV)
!    CALL MAT_MOSTRAR(AINV, '(F7.2)')
    
!    CALL SUB_GAUSS(A, BMAT)
    
!    CALL MET_LU_CROUT(A, B, X)
!    CALL VEC_MOSTRAR(X)
!    CALL VEC_GUARDAR(X, BANDERA, 'LUCrout.txt')
!    DEALLOCATE(X)

!    PRINT *, SENS_NROCOND(A)
!    PRINT *, MAT_NORMAM(A)
!    PRINT *, (0.02*MAT_NORMAM(A))/SENS_NROCOND(A)
!    CALL PRUEBA(A, BMAT)
!    TOLERANCIA = 0.00025
!    CALL SUB_JACOBI(A, B, TOLERANCIA)
!    CALL SUB_GS(A, B, TOLERANCIA)
CONTAINS
    SUBROUTINE SUB_GAUSS(A, B)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A, B
        !
        REAL(8), ALLOCATABLE, DIMENSION(:,:) :: MATAMP, X
        INTEGER :: BANDERA

        CALL MET_GAUSS(A, B, MATAMP)

        !CALL MAT_GUARDAR(MATAMP,BANDERA, ARCHIVO = 'Gauss.txt'); IF (BANDERA /= 0) PRINT *, 'Error al guardar la matriz ampliada de gauss resuelta.'
        
        CALL SUST_REGRESIVA(MATAMP, X)
        
        CALL MAT_GUARDAR(X,BANDERA, ARCHIVO = 'Gauss resuelto.txt')
        IF (BANDERA /= 0) PRINT *, 'Error al guardar la matriz ampliada de gauss resuelta.'
    END SUBROUTINE
    
    SUBROUTINE SUB_GAUSSJORDAN(A, B)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A, B
        !
        REAL(8), ALLOCATABLE, DIMENSION(:,:) :: MATAMP
        INTEGER :: BANDERA
        
        CALL MET_GAUSSJORDAN(A, B, MATAMP)
        
        CALL MAT_GUARDAR(MATAMP,BANDERA, ARCHIVO = 'Gauss-Jordan.txt')
        IF (BANDERA /= 0) PRINT *, 'Error al guardar la matriz ampliada de gauss-jordan.'
    END SUBROUTINE
    
!----------Métodos Indirectos----------!
    SUBROUTINE SUB_JACOBI(A, B, TOL)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A
        REAL(8), DIMENSION(:), INTENT(IN) :: B
        REAL(8), INTENT(IN) :: TOL
        !
        REAL(8), DIMENSION(:), ALLOCATABLE :: X
        INTEGER :: BANDERA
        
        !Pivoteo
        !CALL PIVOTEO
        !Chequear si es edd
!        IF (.NOT. MAT_ISEDD(A)) PRINT *, 'La matriz de entrada en Jacobi no es EDD.'
        !Si no es EDD habria que verificar las combinatorias de filas hasta que alguna dé.
        CALL MET_JACOBI(A, B, TOL, X)
!        CALL VEC_GUARDAR(X, BANDERA, ARCHIVO = 'Jacobi.txt')
!        IF (BANDERA /= 0) PRINT *, 'Error al guardar en SUB_JACOBI.'
        
        DEALLOCATE(X)
    END SUBROUTINE
        
    SUBROUTINE SUB_GS(A, B, TOL)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A
        REAL(8), DIMENSION(:), INTENT(IN) :: B
        REAL(8), INTENT(IN) :: TOL
        !
        REAL(8), DIMENSION(:), ALLOCATABLE :: X
        INTEGER :: BANDERA
        
        !Pivoteo
        !CALL PIVOTEO
        !Chequear si es edd
        IF (.NOT. MAT_ISEDD(A)) PRINT *, 'La matriz de entrada en Gauss Seidel no es EDD.'
        !Si no es EDD habria que verificar las combinatorias de filas hasta que alguna dé.
        CALL MET_GS(A, B, TOL, X)
        CALL VEC_GUARDAR(X, BANDERA, ARCHIVO = 'Gauss-Seidel.txt')
        IF (BANDERA /= 0) PRINT *, 'Error al guardar en SUB_GS.'
        
        DEALLOCATE(X)
    END SUBROUTINE
END PROGRAM
