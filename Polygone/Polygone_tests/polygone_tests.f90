module polygone_tests    
implicit none
    
    contains
    
    
    
    
    subroutine test_star5
    use polygone, only: STAR5_new
    implicit none

    integer(8), parameter :: ND=5, N=3, M=2
    real(8) :: s1, s2
    complex(8) :: A(ND,N), B1(ND,M), B2(ND,M)
    complex(8) :: C1(ND,N), C2(ND,N)
    complex(8) :: R1(ND,M+3), R2(ND,M+3)

    integer(8) :: i, j
    real(8) :: tol
    logical :: ok

    tol = 1.0d-10
    ok = .true.

    ! =========================
    ! 1. »Õ»÷»¿À»«¿÷»ﬂ ƒ¿ÕÕ€’
    ! =========================

    call random_seed()

    do i = 1, N
        do j = 1, N
            call random_number(s1)
            call random_number(s2)
            A(i,j) = cmplx(s1, s2, kind=8)
        enddo
    enddo

    do i = 1, N
        do j = 1, M
            call random_number(s1)
            call random_number(s2)
            B1(i,j) = cmplx(s1, s2, kind=8)
        enddo
    enddo

    B2 = B1   ! ÍÓÔËˇ ‰Îˇ ‚ÚÓÓ„Ó ‡Î„ÓËÚÏ‡

    ! =========================
    ! 2. ¬€«Œ¬ Œ¡Œ»’ ¿À√Œ–»“ÃŒ¬
    ! =========================

    call STAR5(A, B1, C1, R1, ND, N, M, 1)
    call STAR5_new(A, B2, C2, R2, ND, N, M, 1)

    ! =========================
    ! 3. —–¿¬Õ≈Õ»≈ –≈«”À‹“¿“Œ¬
    ! =========================

    do i = 1, N
        do j = 1, M
            if (abs(B1(i,j) - B2(i,j)) > tol) then
                print *, "Mismatch at (", i, ",", j, ")"
                print *, "Old:", B1(i,j)
                print *, "New:", B2(i,j)
                ok = .false.
            endif
        enddo
    enddo

    ! =========================
    ! 4. »“Œ√
    ! =========================

    if (ok) then
        print *, "TEST PASSED: solutions match"
    else
        print *, "TEST FAILED"
    endif

    end subroutine test_star5
    
    
    
    
    
    
    
    
    subroutine test_speed_star5
    use polygone, only: STAR5_new
    implicit none

    integer(8), parameter :: ND=100, N=80, M=10
    integer(8), parameter :: n_iter = 200

    complex(8) :: A(ND,N), B1(ND,M), B2(ND,M)
    complex(8) :: C(ND,N), R(ND,M+3)

    real(8) :: t1, t2, time_old, time_new
    real(8) :: s1, s2
    integer :: i, j, iter

    ! =========================
    ! »Õ»÷»¿À»«¿÷»ﬂ
    ! =========================

    call random_seed()

    do i = 1, N
        do j = 1, N
            call random_number(s1)
            call random_number(s2)
            A(i,j) = cmplx(s1, s2)
        enddo
    enddo

    do i = 1, N
        do j = 1, M
            call random_number(s1)
            call random_number(s2)
            B1(i,j) = cmplx(s1, s2)
        enddo
    enddo

    ! =========================
    ! “≈—“ —“¿–Œ… ¬≈–—»»
    ! =========================

    call cpu_time(t1)

    do iter = 1, n_iter
        B2 = B1
        call STAR5(A, B2, C, R, ND, N, M, 1)
    enddo

    call cpu_time(t2)
    time_old = t2 - t1

    ! =========================
    ! “≈—“ ÕŒ¬Œ… ¬≈–—»»
    ! =========================

    call cpu_time(t1)

    do iter = 1, n_iter
        B2 = B1
        call STAR5_new(A, B2, C, R, ND, N, M, 1)
    enddo

    call cpu_time(t2)
    time_new = t2 - t1

    ! =========================
    ! –≈«”À‹“¿“
    ! =========================

    print *, "OLD TIME:", time_old
    print *, "NEW TIME:", time_new
    print *, "SPEEDUP (old/new):", time_old / time_new

    end subroutine test_speed_star5
    
    
    
    
    
    
    
    
    
    
    
    
end module polygone_tests