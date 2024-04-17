PROGRAM ParticulasConResortes
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: N = 10 ! Número de partículas
    REAL :: dt = 0.01            ! Paso de tiempo
    REAL :: t_max = 10.0         ! Tiempo máximo de simulación
    INTEGER :: i, j
    REAL, PARAMETER :: k = 1.0   ! Constante del resorte
    REAL :: g = 9.81             ! Aceleración gravitatoria
    
    REAL, DIMENSION(N, 2) :: x, v, k1, k2, k3, k4
    REAL :: t
    
    ! Condiciones iniciales
    x(:, 1) = 1.0
    x(:, 2) = 1.0
    v(:, 1) = 1.0
    v(:, 2) = 1.0
    
    OPEN(UNIT=1, FILE='posiciones_x.txt', STATUS='UNKNOWN')
    OPEN(UNIT=2, FILE='posiciones_y.txt', STATUS='UNKNOWN')
    OPEN(UNIT=3, FILE='posiciones_xy.txt', STATUS='UNKNOWN')
    
    t = 0.0
    
    ! Loop temporal
    DO WHILE (t < t_max)
        ! Escribir posiciones en los archivos
        WRITE(1, *) t, (x(i, 1), i = 1, N)
        WRITE(2, *) t, (x(i, 2), i = 1, N)
        WRITE(3, *) t, ((x(i, j), j=1,2), i = 1, N)
        
        ! Calcular k1
        CALL F(N, x, v, k1)
        ! Calcular k2
        CALL F(N, x + 0.5*dt*k1, v, k2)
        ! Calcular k3
        CALL F(N, x + 0.5*dt*k2, v, k3)
        ! Calcular k4
        CALL F(N, x + dt*k3, v, k4)
        
        ! Actualizar posición y velocidad
        x = x + dt/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)
        
        t = t + dt
    END DO
    
    CLOSE(UNIT=1)
    CLOSE(UNIT=2)
    CLOSE(UNIT=3)
    
    END PROGRAM ParticulasConResortes

SUBROUTINE F(N, x, v, F_out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    REAL, DIMENSION(N, 2), INTENT(IN) :: x, v
    REAL, DIMENSION(N, 2), INTENT(OUT) :: F_out
    INTEGER :: i
    
    REAL :: g,k
    
    k = 1
    g = 9.81
    
    F_out(1, :) = 0.0 ! Fuerza en la primera partícula
    F_out(N, :) = 0.0 ! Fuerza en la última partícula
    
    DO i = 2, N-1
        ! Fuerza en dirección x
        F_out(i, 1) = -k*(x(i, 1) - x(i-1, 1)) + k*(x(i+1, 1) - x(i, 1))
        ! Fuerza en dirección y (incluyendo la gravedad)
        F_out(i, 2) = -k*(x(i, 2) - x(i-1, 2)) + k*(x(i+1, 2) - x(i, 2)) - g
    END DO
    
END SUBROUTINE F

