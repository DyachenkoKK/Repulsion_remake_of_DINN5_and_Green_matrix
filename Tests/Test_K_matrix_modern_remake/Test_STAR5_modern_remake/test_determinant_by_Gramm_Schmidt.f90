module test_determinant_by_Gramm_Schmidt    
implicit none
    
contains
   subroutine test_get_determinant_by_Gramm_Schmidt
    use determinant_by_Gramm_Schmidt, only: get_determinant_by_Gramm_Schmidt
    implicit none
    
    ! Параметры тестов
    integer, parameter :: num_tests = 6
    integer(8) :: ND, N, KORT
    integer :: test_num
    complex(8), allocatable :: A(:, :), A_copy(:, :)
    complex(8), allocatable :: C1(:, :), C2(:, :), S1(:), S2(:), SM1(:), SM2(:)
    complex(8) :: det_old, det_new
    real(8) :: relative_error
    logical :: test_passed
    
    ! Константы
    real(8), parameter :: EPS = 1.0D-10
    integer(8), parameter :: KORT_VAL = 3
    
    print *, "============================================="
    print *, "Testing DSTAR vs DSTAR_new"
    print *, "============================================="
    print *, ""
    
    ! ========================================================================
    ! ТЕСТ 1: Матрица 2x2 с известным определителем
    ! ========================================================================
    test_num = 1
    N = 2
    ND = 2
    
    allocate(A(ND, N), A_copy(ND, N))
    allocate(C1(ND, N), C2(ND, N))
    allocate(S1(N), S2(N), SM1(N), SM2(N))
    
    ! Матрица: [1+2i, 3+4i; 5+6i, 7+8i]
    ! Определитель = -16i
    A(1,1) = (1.0d0, 2.0d0); A(1,2) = (3.0d0, 4.0d0)
    A(2,1) = (5.0d0, 6.0d0); A(2,2) = (7.0d0, 8.0d0)
    A_copy = A
    
    call test_dstar_comparison(A, A_copy, C1, C2, S1, S2, SM1, SM2, &
                               ND, N, KORT_VAL, det_old, det_new, &
                               test_num, EPS)
    
    deallocate(A, A_copy, C1, C2, S1, S2, SM1, SM2)
    
    ! ========================================================================
    ! ТЕСТ 2: Единичная матрица 3x3 (определитель = 1)
    ! ========================================================================
    test_num = 2
    N = 3
    ND = 3
    
    allocate(A(ND, N), A_copy(ND, N))
    allocate(C1(ND, N), C2(ND, N))
    allocate(S1(N), S2(N), SM1(N), SM2(N))
    
    call identity_matrix_complex(A, ND, N)
    A_copy = A
    
    call test_dstar_comparison(A, A_copy, C1, C2, S1, S2, SM1, SM2, &
                               ND, N, KORT_VAL, det_old, det_new, &
                               test_num, EPS)
    
    deallocate(A, A_copy, C1, C2, S1, S2, SM1, SM2)
    
    ! ========================================================================
    ! ТЕСТ 3: Диагональная матрица 4x4
    ! ========================================================================
    test_num = 3
    N = 4
    ND = 4
    
    allocate(A(ND, N), A_copy(ND, N))
    allocate(C1(ND, N), C2(ND, N))
    allocate(S1(N), S2(N), SM1(N), SM2(N))
    
    call diagonal_matrix_complex(A, ND, N)
    A_copy = A
    
    call test_dstar_comparison(A, A_copy, C1, C2, S1, S2, SM1, SM2, &
                               ND, N, KORT_VAL, det_old, det_new, &
                               test_num, EPS)
    
    deallocate(A, A_copy, C1, C2, S1, S2, SM1, SM2)
    
    ! ========================================================================
    ! ТЕСТ 4: Матрица 5x5 с целочисленным определителем
    ! ========================================================================
    test_num = 4
    N = 5
    ND = 5
    
    allocate(A(ND, N), A_copy(ND, N))
    allocate(C1(ND, N), C2(ND, N))
    allocate(S1(N), S2(N), SM1(N), SM2(N))
    
    call integer_matrix_test(A, ND, N)
    A_copy = A
    
    call test_dstar_comparison(A, A_copy, C1, C2, S1, S2, SM1, SM2, &
                               ND, N, KORT_VAL, det_old, det_new, &
                               test_num, EPS)
    
    deallocate(A, A_copy, C1, C2, S1, S2, SM1, SM2)
    
    ! ========================================================================
    ! ТЕСТ 5: Случайная матрица 10x10
    ! ========================================================================
    test_num = 5
    N = 10
    ND = 10
    
    allocate(A(ND, N), A_copy(ND, N))
    allocate(C1(ND, N), C2(ND, N))
    allocate(S1(N), S2(N), SM1(N), SM2(N))
    
    call random_matrix_complex(A, ND, N)
    A_copy = A
    
    call test_dstar_comparison(A, A_copy, C1, C2, S1, S2, SM1, SM2, &
                               ND, N, KORT_VAL, det_old, det_new, &
                               test_num, EPS)
    
    deallocate(A, A_copy, C1, C2, S1, S2, SM1, SM2)
    
    ! ========================================================================
    ! ТЕСТ 6: Вырожденная матрица (определитель = 0)
    ! ========================================================================
    test_num = 6
    N = 3
    ND = 3
    
    allocate(A(ND, N), A_copy(ND, N))
    allocate(C1(ND, N), C2(ND, N))
    allocate(S1(N), S2(N), SM1(N), SM2(N))
    
    call singular_matrix(A, ND, N)
    A_copy = A
    
    call test_dstar_comparison(A, A_copy, C1, C2, S1, S2, SM1, SM2, &
                               ND, N, KORT_VAL, det_old, det_new, &
                               test_num, EPS)
    
    deallocate(A, A_copy, C1, C2, S1, S2, SM1, SM2)
    
    print *, "============================================="
    print *, "All tests completed!"
    print *, "============================================="

contains

    !==========================================================================
    ! Вспомогательные подпрограммы
    !==========================================================================
    
    subroutine identity_matrix_complex(A, ND, N)
        integer(8), intent(in) :: ND, N
        complex(8), intent(out) :: A(ND, N)
        integer(8) :: i, j
        
        do i = 1, N
            do j = 1, N
                if (i == j) then
                    A(i, j) = (1.0d0, 0.0d0)
                else
                    A(i, j) = (0.0d0, 0.0d0)
                end if
            end do
        end do
    end subroutine identity_matrix_complex
    
    subroutine diagonal_matrix_complex(A, ND, N)
        integer(8), intent(in) :: ND, N
        complex(8), intent(out) :: A(ND, N)
        integer(8) :: i, j
        
        do i = 1, N
            do j = 1, N
                if (i == j) then
                    A(i, j) = cmplx(dble(i) * 2.0d0, 0.0d0, kind=8)
                else
                    A(i, j) = (0.0d0, 0.0d0)
                end if
            end do
        end do
    end subroutine diagonal_matrix_complex
    
    subroutine integer_matrix_test(A, ND, N)
        integer(8), intent(in) :: ND, N
        complex(8), intent(out) :: A(ND, N)
        integer(8) :: i, j
        
        ! Матрица Гильберта (известный определитель)
        do i = 1, N
            do j = 1, N
                A(i, j) = cmplx(1.0d0 / dble(i + j - 1), 0.0d0, kind=8)
            end do
        end do
    end subroutine integer_matrix_test
    
    subroutine random_matrix_complex(A, ND, N)
        integer(8), intent(in) :: ND, N
        complex(8), intent(out) :: A(ND, N)
        integer(8) :: i, j
        real(8) :: r1, r2
        
        call init_random_seed_test()
        
        do i = 1, N
            do j = 1, N
                call random_number(r1)
                call random_number(r2)
                A(i, j) = cmplx((r1 - 0.5d0) * 20.0d0, (r2 - 0.5d0) * 20.0d0, kind=8)
            end do
        end do
    end subroutine random_matrix_complex
    
    subroutine singular_matrix(A, ND, N)
        integer(8), intent(in) :: ND, N
        complex(8), intent(out) :: A(ND, N)
        integer(8) :: i, j
        
        ! Создаем матрицу с двумя линейно зависимыми строками
        do i = 1, N
            do j = 1, N
                A(i, j) = cmplx(dble(i + j), 0.0d0, kind=8)
            end do
        end do
        ! Делаем последнюю строку равной сумме первых двух
        if (N >= 3) then
            A(N, :) = A(1, :) + A(2, :)
        end if
    end subroutine singular_matrix
    
    subroutine init_random_seed_test()
        integer :: i, n, clock
        integer, allocatable :: seed(:)
        
        call random_seed(size=n)
        allocate(seed(n))
        
        call system_clock(count=clock)
        seed = clock + 37 * (/(i-1, i=1, n)/)
        call random_seed(put=seed)
        
        deallocate(seed)
    end subroutine init_random_seed_test
    
    subroutine test_dstar_comparison(A, A_copy, C1, C2, S1, S2, SM1, SM2, &
                                     ND, N, KORT, det_old, det_new, &
                                     test_num, EPS)
        integer(8), intent(in) :: ND, N, KORT, test_num
        complex(8), intent(in) :: A(ND, N), A_copy(ND, N)
        complex(8), intent(out) :: C1(ND, N), C2(ND, N)
        complex(8), intent(out) :: S1(N), S2(N), SM1(N), SM2(N)
        complex(8), intent(out) :: det_old, det_new
        real(8), intent(in) :: EPS
        
        ! Вызываем старую версию DSTAR
        call DSTAR(A, C1, S1, SM1, det_old, ND, N, KORT)
        
        ! Вызываем новую версию DSTAR_new
        call get_determinant_by_Gramm_Schmidt(A_copy, C2, S2, SM2, det_new, ND, N, KORT)
        
        ! Вычисляем относительную ошибку
        if (abs(det_old) > EPS) then
            relative_error = abs(det_new - det_old) / abs(det_old)
        else
            relative_error = abs(det_new - det_old)
        end if
        
        test_passed = relative_error < EPS
        
        print '(A, I2, A)', " Test ", test_num, ":"
        print '(A, 2F12.6)', "   Old determinant: ", det_old
        print '(A, 2F12.6)', "   New determinant: ", det_new
        print '(A, ES12.4)', "   Relative error:  ", relative_error
        print '(A, L1)',      "   Test passed:     ", test_passed
        print *, ""
        
    end subroutine test_dstar_comparison

end subroutine test_get_determinant_by_Gramm_Schmidt

    
    
end module test_determinant_by_Gramm_Schmidt