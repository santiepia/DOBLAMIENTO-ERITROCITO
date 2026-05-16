module funcion_volumen
contains

    real function f(x, A, B, C)
        implicit none
        real, parameter :: pi = 3.141592653589793   
        real :: x, A, B, C

  

      
        f = 2 * pi * (B * ((C**4.0 + 4.0 * A**2.0 * x**2.0)**0.5 - A**2.0 - x**2.0)**0.5)**2.0
                
        return
    end function f
    
end module funcion_volumen
    

module funcion_area_superficial
contains

    real function g(x, A, B,C)
        implicit none
        real, parameter :: pi = 3.141592653589793

        real :: x, A, B, C
        real :: R, S, Y, dY_dx
        real :: termino1, termino2
        

        
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

    integer :: N, i, i_A, i_B, n_A, n_B, i_C, n_C
    real :: a_limite, b_limite, h
    real, allocatable :: x(:)
    real :: suma1_volumen, suma2_volumen, suma1_area, suma2_area
    real :: AreaSimpson1, AreaSimpson2
    real :: A_valor, A_inicio, A_fin, paso_A
    real :: B_valor, B_inicio, B_fin, paso_B
    real :: C_valor, C_inicio, C_fin, paso_C

    
    ! Parámetros fijos

    h = 0.1              ! Paso de integración
    a_limite = 0.0     ! Límite inferior (siempre 0 por simetría)
    
    ! Rango de A (de 1 a 10, paso 0.1)
    A_inicio = 1.0
    A_fin = 1.1
    paso_A = 0.1
    n_A = nint((A_fin - A_inicio) / paso_A) + 1

    ! Rango de B
    B_inicio = 1.0
    B_fin = 1.1
    paso_B = 0.1
    n_B = nint((B_fin - B_inicio) / paso_B) + 1


    C_inicio = 1.0
    C_fin = 1.1
    paso_C = 0.1
    n_C = nint((C_fin - C_inicio) / paso_C) + 1


    open(10, file='resultados.txt', status='replace')  ! Abrir archivo para escribir resultados
    ! Encabezado de la tabla
    write(10,*) '============================================================'
    write(10,*) '   A      B      C      b_limite      Volumen           Area Superficial'
    write(10,*) '============================================================'
    
    ! Bucle sobre A
    do i_A = 0, n_A - 1

        A_valor = A_inicio + i_A * paso_A
        
        ! Calcular el límite superior dependiendo de A
        b_limite = sqrt(A_valor**2 + C_valor**2)
        
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
        
!INICIA BUCLE PARA B

         do i_B = 0, n_B - 1
            B_valor = B_inicio + i_B * paso_B
            




        !INICIANDO BUCLE PARA C
            do i_C = 0, n_C - 1

                C_valor = C_inicio + i_C * paso_C

                ! Inicializar sumas
                suma1_volumen = 0.0
                suma2_volumen = 0.0
                suma1_area = 0.0
                suma2_area = 0.0

                ! Puntos impares
                do i = 1, N-1, 2
                    suma1_volumen = suma1_volumen + f(x(i), A_valor, B_valor,C_valor)
                    suma1_area = suma1_area + g(x(i), A_valor, B_valor,C_valor)
                end do
                
                ! Puntos pares
                do i = 2, N-2, 2
                    suma2_volumen = suma2_volumen + f(x(i), A_valor, B_valor,C_valor)
                    suma2_area = suma2_area + g(x(i), A_valor, B_valor,C_valor)
                end do
                
                ! Regla de Simpson (líneas corregidas)
                AreaSimpson1 = (h/3.0) * (f(a_limite, A_valor, B_valor,C_valor) + &
                            f(b_limite, A_valor, B_valor,C_valor) + 4.0 * suma1_volumen + &
                            2.0 * suma2_volumen)
                
                AreaSimpson2 = (h/3.0) * (g(a_limite, A_valor, B_valor,C_valor) + &
                            g(b_limite, A_valor, B_valor,C_valor) + 4.0 * suma1_area + &
                            2.0 * suma2_area)

                ! Imprimir
                write(10, '(F6.2, F8.2, F10.4, 2ES18.6)') A_valor, B_valor, C_valor, &
                    b_limite, AreaSimpson1, AreaSimpson2
                
                end do  ! Fin del bucle para C
            
            
        end do
        
        deallocate(x)
        write(10,*) '----------------------------------------------------------------------'
        
    end do
    
    write(10,*) '======================================================================'
    close(10)

end program Simpson


        


