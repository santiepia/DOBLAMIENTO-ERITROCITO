module funcion_volumen
contains

    real function f(x, A)
        implicit none
        real, parameter :: pi = 3.141592653589793
        real :: B, C   
        real :: x, A

        B = 18.73
        C = 1.05
      
        f = 2 * pi * (B * ((C**4.0 + 4.0 * A**2.0 * x**2.0)**0.5 - A**2.0 - x**2.0)**0.5)**2.0
                
        return
    end function f
    
end module funcion_volumen
    

module funcion_area_superficial
contains

    real function g(x, A)
        implicit none
        real, parameter :: pi = 3.141592653589793
        real :: B, C   
        real :: x, A
        real :: R, S, Y, dY_dx
        real :: termino1, termino2
        
        B = 18.73
        C = 1.05
        
        R = sqrt(C**4 + 4.0 * A**2 * x**2)
        S = R - A**2 - x**2
        
        if (S <= 0.0) then
            g = 0.0
            return
        end if
        
        Y = B * sqrt(S)
        termino1 = (4.0 * A**2 * x) / R - 2.0 * x
        termino2 = 2.0 * sqrt(S)
        dY_dx = B * termino1 / termino2
        
        g = 4.0 * pi * Y * sqrt(1.0 + dY_dx**2) !!recordar que la formula del area superficial es 4*pi*Y*sqrt(1+(dY/dx)^2) para integrales pares, por lo que no se hace de -l a l sino de 0 a l
        
        return
    end function g
    
end module funcion_area_superficial



program Simpson

    use funcion_volumen
    use funcion_area_superficial

    integer :: N, i, j, n_A
    real :: a_limite, b_limite, h
    real, allocatable :: x(:)
    real :: suma1_volumen, suma2_volumen, suma1_area, suma2_area
    real :: AreaSimpson1, AreaSimpson2
    real :: A_valor, A_inicio, A_fin, paso_A
    real :: C
    
    ! Parámetros fijos
    C = 1.05           ! Parámetro C constante
    h = 0.0001         ! Paso de integración
    a_limite = 0.0     ! Límite inferior (siempre 0 por simetría)
    
    ! Rango de A (de 1 a 10, paso 0.1)
    A_inicio = 1.0
    A_fin = 1.5
    paso_A = 0.001
    n_A = nint((A_fin - A_inicio) / paso_A) + 1
    
    ! Encabezado de la tabla
    write(*,*) '============================================================'
    write(*,*) '   A       b_limite      Volumen           Area Superficial'
    write(*,*) '============================================================'
    
    ! Bucle sobre A
    do j = 0, n_A - 1
        A_valor = A_inicio + j * paso_A
        
        ! Calcular el límite superior dependiendo de A
        b_limite = sqrt(A_valor**2 + C**2)
        
        ! Calcular N para la integración
        N = (b_limite - a_limite) / h
        if (mod(N,2) /= 0) N = N - 1
        b_limite = a_limite + N * h
        
        ! Reservar memoria para x
        allocate(x(0:N))
        
        ! Generar puntos x
        x(0) = a_limite
        do i = 1, N
            x(i) = x(i-1) + h
        end do
        
        ! Inicializar sumas
        suma1_volumen = 0.0
        suma2_volumen = 0.0
        suma1_area = 0.0
        suma2_area = 0.0
        
        ! Sumas para Simpson (puntos impares -> peso 4)
        do i = 1, N-1, 2
            suma1_volumen = suma1_volumen + f(x(i), A_valor)
            suma1_area = suma1_area + g(x(i), A_valor)
        end do
        
        ! Sumas para Simpson (puntos pares -> peso 2)
        do i = 2, N-2, 2
            suma2_volumen = suma2_volumen + f(x(i), A_valor)
            suma2_area = suma2_area + g(x(i), A_valor)
        end do
        
        ! Regla de Simpson
        AreaSimpson1 = (h/3.0) * (f(a_limite, A_valor) + f(b_limite, A_valor) + 4.0 * suma1_volumen + 2.0 * suma2_volumen)
        AreaSimpson2 = (h/3.0) * (g(a_limite, A_valor) + g(b_limite, A_valor) + 4.0 * suma1_area + 2.0 * suma2_area)
        
        ! Imprimir fila de la tabla
        write(*, '(F6.2, F10.4, 2ES18.6)') A_valor, b_limite, AreaSimpson1, AreaSimpson2
        
        deallocate(x)
        
    end do
    
    write(*,*) '============================================================'

end program Simpson