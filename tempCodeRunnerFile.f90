            write(*,*) 'ERROR: N <= 0 para A =', A_valor
            stop
        end if

        if (mod(N,2) /= 0) N = N - 1
        if (N <= 0) then
            write(*,*) 'ERROR: N <= 0 después de hacerlo par'
            stop
        end if