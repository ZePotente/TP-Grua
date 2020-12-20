PROGRAM GRUA
    !Modulo
    USE SEL_MET
    USE VYM_IO
    
    IMPLICIT NONE
    
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: A
    REAL(8), DIMENSION(:), ALLOCATABLE :: B
    INTEGER :: BANDERAA, BANDERAB
    
    CALL MAT_LEER(A, BANDERAA, 'A.txt')
    CALL VEC_LEER(B, BANDERAB, 'B.txt')
    
    IF (.NOT. VERIF_LECT(BANDERAA, BANDERAB)) GOTO 20
    CALL MOSTRAR_DATOS(A, B)
    
    CALL PRUEBA_METODOS(A, B)
    
20  PRINT *, 'Error al leer de archivo.'
CONTAINS

    SUBROUTINE MOSTRAR_DATOS(A, B)
        REAL(8), INTENT(IN) :: A(:,:), B(:)
        !
        CHARACTER(*), PARAMETER :: FORMATO = '(F10.5)'
        
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
    
    SUBROUTINE PRUEBA_METODOS(A, B)
        REAL(8), INTENT(IN) :: A(:,:), B(:)
        !
        REAL(8), DIMENSION(:), ALLOCATABLE :: XGAUSS, XGJ!, XCROUT
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: BGJ
        
        
    END SUBROUTINE
END PROGRAM
