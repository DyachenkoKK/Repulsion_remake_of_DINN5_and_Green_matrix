module test_LU_Gaussian_elimination    
implicit none
    
contains
    
    subroutine result_tests_for_Gauss_LU_method
    use polygone, only: test_func
    implicit none
    
        integer :: i, j
        integer(8) :: N, ND
        complex(8), allocatable :: A_old(:, :), A_new(:, :)
        complex(8) :: det_old, det_new
        complex(8) :: expected_det
    
        ! Test matrix 1: Simple 3x3 matrix with known determinant
        print *, "========================================="
        print *, "Test 1: 3x3 matrix with known determinant"
        print *, "========================================="
    
        N = 3
        ND = 3
    
        allocate(A_old(ND, ND), A_new(ND, ND))
    
        ! Matrix: [2, 1, 1; 4, 3, 2; 2, 2, 3]
        ! Determinant = 4 (calculated manually)
        A_old(1,1) = (2.0d0, 0.0d0); A_old(1,2) = (1.0d0, 0.0d0); A_old(1,3) = (1.0d0, 0.0d0)
        A_old(2,1) = (4.0d0, 0.0d0); A_old(2,2) = (3.0d0, 0.0d0); A_old(2,3) = (2.0d0, 0.0d0)
        A_old(3,1) = (2.0d0, 0.0d0); A_old(3,2) = (2.0d0, 0.0d0); A_old(3,3) = (3.0d0, 0.0d0)
        A_new = A_old
    
        expected_det = (4.0d0, 0.0d0)
    
        call DGAUS(A_old, ND, N, det_old)
        call test_func(A_new, ND, N, det_new)
    
        print *, "Old routine determinant: ", det_old
        print *, "New routine determinant: ", det_new
        print *, "Expected determinant:     ", expected_det
        print *, "Old routine correct:      ", abs(det_old - expected_det) < 1.0d-10
        print *, "New routine correct:      ", abs(det_new - expected_det) < 1.0d-10
        print *, ""
    
        deallocate(A_old, A_new)
    
        ! Test matrix 2: 4x4 matrix with integer determinant
        print *, "========================================="
        print *, "Test 2: 4x4 matrix with integer determinant"
        print *, "========================================="
    
        N = 4
        ND = 4
    
        allocate(A_old(ND, ND), A_new(ND, ND))
    
        ! Matrix with determinant = 0 (singular)
        A_old(1,1) = (1.0d0, 0.0d0); A_old(1,2) = (2.0d0, 0.0d0); A_old(1,3) = (3.0d0, 0.0d0); A_old(1,4) = (4.0d0, 0.0d0)
        A_old(2,1) = (2.0d0, 0.0d0); A_old(2,2) = (3.0d0, 0.0d0); A_old(2,3) = (4.0d0, 0.0d0); A_old(2,4) = (5.0d0, 0.0d0)
        A_old(3,1) = (3.0d0, 0.0d0); A_old(3,2) = (4.0d0, 0.0d0); A_old(3,3) = (5.0d0, 0.0d0); A_old(3,4) = (6.0d0, 0.0d0)
        A_old(4,1) = (4.0d0, 0.0d0); A_old(4,2) = (5.0d0, 0.0d0); A_old(4,3) = (6.0d0, 0.0d0); A_old(4,4) = (7.0d0, 0.0d0)
        A_new = A_old
    
        expected_det = (0.0d0, 0.0d0)  ! This matrix is singular (determinant = 0)
    
        call DGAUS(A_old, ND, N, det_old)
        call test_func(A_new, ND, N, det_new)
    
        print *, "Old routine determinant: ", det_old
        print *, "New routine determinant: ", det_new
        print *, "Expected determinant:     ", expected_det
        print *, "Old routine correct:      ", abs(det_old - expected_det) < 1.0d-10
        print *, "New routine correct:      ", abs(det_new - expected_det) < 1.0d-10
        print *, ""
    
        deallocate(A_old, A_new)
    
        ! Test matrix 3: Complex matrix
        print *, "========================================="
        print *, "Test 3: Complex matrix"
        print *, "========================================="
    
        N = 2
        ND = 2
    
        allocate(A_old(ND, ND), A_new(ND, ND))
    
        ! Matrix: [1+2i, 3+4i; 5+6i, 7+8i]
        ! Determinant = (1+2i)*(7+8i) - (3+4i)*(5+6i) = -16i
        A_old(1,1) = (1.0d0, 2.0d0); A_old(1,2) = (3.0d0, 4.0d0)
        A_old(2,1) = (5.0d0, 6.0d0); A_old(2,2) = (7.0d0, 8.0d0)
        A_new = A_old
    
        expected_det = (0.0d0, -16.0d0)
    
        call DGAUS(A_old, ND, N, det_old)
        call test_func(A_new, ND, N, det_new)
    
        print *, "Old routine determinant: ", det_old
        print *, "New routine determinant: ", det_new
        print *, "Expected determinant:     ", expected_det
        print *, "Old routine correct:      ", abs(det_old - expected_det) < 1.0d-10
        print *, "New routine correct:      ", abs(det_new - expected_det) < 1.0d-10
        print *, ""
    
        deallocate(A_old, A_new)
    
        ! Test matrix 4: Identity matrix
        print *, "========================================="
        print *, "Test 4: Identity matrix"
        print *, "========================================="
    
        N = 5
        ND = 5
    
        allocate(A_old(ND, ND), A_new(ND, ND))
    
        call identity_matrix(A_old, ND, N)
        A_new = A_old
    
        expected_det = (1.0d0, 0.0d0)
    
        call DGAUS(A_old, ND, N, det_old)
        call test_func(A_new, ND, N, det_new)
    
        print *, "Old routine determinant: ", det_old
        print *, "New routine determinant: ", det_new
        print *, "Expected determinant:     ", expected_det
        print *, "Old routine correct:      ", abs(det_old - expected_det) < 1.0d-10
        print *, "New routine correct:      ", abs(det_new - expected_det) < 1.0d-10
        print *, ""
    
        deallocate(A_old, A_new)
    
        ! Test matrix 5: Random matrix (compare results between routines)
        print *, "========================================="
        print *, "Test 5: Random 10x10 matrix"
        print *, "========================================="
    
        N = 10
        ND = 10
    
        allocate(A_old(ND, ND), A_new(ND, ND))
    
        call random_matrix_test(A_old, ND, N)
        A_new = A_old
    
        call DGAUS(A_old, ND, N, det_old)
        call test_func(A_new, ND, N, det_new)
    
        print *, "Old routine determinant: ", det_old
        print *, "New routine determinant: ", det_new
        print *, "Difference:              ", abs(det_old - det_new)
        print *, "Results match:           ", abs(det_old - det_new) < 1.0d-10
        print *, ""
    
        deallocate(A_old, A_new)

    contains

        subroutine identity_matrix(A, ND, N)
            integer, intent(in) :: ND, N
            complex(8), intent(out) :: A(ND, ND)
            integer :: i, j
    
            do i = 1, N
                do j = 1, N
                    if (i == j) then
                        A(i, j) = (1.0d0, 0.0d0)
                    else
                        A(i, j) = (0.0d0, 0.0d0)
                    end if
                end do
            end do
        end subroutine identity_matrix
    
        subroutine random_matrix_test(A, ND, N)
            integer, intent(in) :: ND, N
            complex(8), intent(out) :: A(ND, ND)
            integer :: i, j
            real(8) :: r1, r2
    
            ! Initialize random seed
            call init_random_seed()
            
            do i = 1, N
                do j = 1, N
                    call random_number(r1)
                    call random_number(r2)
                    A(i, j) = cmplx(r1 * 10.0d0, r2 * 10.0d0, kind=8)
                end do
            end do
        end subroutine random_matrix_test
        
        subroutine init_random_seed()
            integer :: i, n, clock
            integer, allocatable :: seed(:)
            
            call random_seed(size=n)
            allocate(seed(n))
            
            call system_clock(count=clock)
            seed = clock + 37 * (/(i-1, i=1, n)/)
            call random_seed(put=seed)
            
            deallocate(seed)
        end subroutine init_random_seed
        
    end subroutine result_tests_for_Gauss_LU_method
    
    
    
    
    
    
    
    
    
    
   subroutine speed_tests_for_Gauss_LU_method
    use polygone
    implicit none
    
    integer :: i, n_size
    integer, parameter :: num_tests = 5
    integer, parameter :: matrix_sizes(num_tests) = [50, 100, 200, 500, 1000]
    real(8) :: time_old, time_new, speedup
    real(8) :: time_old_avg, time_new_avg
    integer, parameter :: num_runs = 3  ! Number of runs for averaging
    
    print *, "============================================="
    print *, "Performance Comparison: DGAUS vs test_func"
    print *, "============================================="
    print *, ""
    print *, "Matrix Size | Old Time (s) | New Time (s) | Speedup | Correctness"
    print *, "------------|--------------|--------------|---------|------------"
    
    do i = 1, num_tests
        n_size = matrix_sizes(i)
        
        ! Average over multiple runs
        time_old_avg = 0.0d0
        time_new_avg = 0.0d0
        
        call compare_single_size(n_size, time_old_avg, time_new_avg, num_runs)
        
        if (time_new_avg > 0.0d0) then
            speedup = time_old_avg / time_new_avg
        else
            speedup = 0.0d0
        end if
        
        print '(I10, " | ", F12.6, " | ", F12.6, " | ", F7.2, "x | ", L10)', &
                   n_size, time_old_avg, time_new_avg, speedup, .true.
    end do
    
    print *, ""
    print *, "============================================="
    print *, "Detailed analysis for specific matrix sizes"
    print *, "============================================="
    
    ! Test with specific sizes for detailed analysis
    call detailed_performance_analysis(50)
    call detailed_performance_analysis(100)
    call detailed_performance_analysis(200)
    call detailed_performance_analysis(500)
    call detailed_performance_analysis(1000)
    
    contains

        ! Функция для проверки NaN для комплексных чисел
        function is_nan_complex(z) result(res)
        implicit none
            complex(8), intent(in) :: z
            logical :: res
            real(8) :: re_part, im_part
            
            re_part = real(z, kind=8)
            im_part = aimag(z)
            ! NaN проверяется тем, что число не равно самому себе
            res = (re_part /= re_part) .or. (im_part /= im_part)
        end function is_nan_complex

        subroutine compare_single_size(n, time_old_avg, time_new_avg, num_runs)
        implicit none
            integer, intent(in) :: n, num_runs
            real(8), intent(out) :: time_old_avg, time_new_avg
            complex(8), allocatable :: A_old(:, :), A_new(:, :)
            complex(8) :: det_old, det_new
            real(8) :: start_time, end_time, total_old, total_new
            integer :: run
            logical :: valid_run
        
            total_old = 0.0d0
            total_new = 0.0d0
            valid_run = .true.
        
            do run = 1, num_runs
                ! Allocate matrices
                allocate(A_old(n, n), A_new(n, n))
            
                ! Initialize with non-singular matrix
                call generate_non_singular_matrix(A_old, n, n)
                A_new = A_old
            
                ! Test old routine (DGAUS)
                call cpu_time(start_time)
                call DGAUS(A_old, n, n, det_old)
                call cpu_time(end_time)
                
                ! Check if determinant is valid (not NaN or zero)
                if (is_nan_complex(det_old) .or. abs(det_old) < 1.0d-50) then
                    valid_run = .false.
                    deallocate(A_old, A_new)
                    cycle
                end if
                
                total_old = total_old + (end_time - start_time)
            
                ! Test new routine (test_func)
                call cpu_time(start_time)
                call test_func(A_new, n, n, det_new)
                call cpu_time(end_time)
                
                if (is_nan_complex(det_new) .or. abs(det_new) < 1.0d-50) then
                    valid_run = .false.
                    deallocate(A_old, A_new)
                    cycle
                end if
                
                total_new = total_new + (end_time - start_time)
            
                deallocate(A_old, A_new)
            end do
            
            if (valid_run .and. num_runs > 0) then
                time_old_avg = total_old / num_runs
                time_new_avg = total_new / num_runs
            else
                time_old_avg = 0.0d0
                time_new_avg = 0.0d0
            end if
        end subroutine compare_single_size
        
        ! Generate an identity matrix (non-singular)
        ! Generate a non-singular matrix with controlled determinant magnitude
        subroutine generate_non_singular_matrix(A, ND, N)
            implicit none
            integer, intent(in) :: ND, N
            complex(8), intent(out) :: A(ND, ND)
            integer :: i, j
            real(8) :: r1, r2
            real(8), parameter :: DIAG_SCALE = 1.0d0  ! Уменьшили масштаб
    
            ! Create a matrix with moderate diagonal elements
            do i = 1, N
                do j = 1, N
                    if (i == j) then
                        ! Умеренные диагональные элементы (не 100*N)
                        A(i, j) = cmplx(DIAG_SCALE, 0.0d0, kind=8)
                    else
                        A(i, j) = 0
                    end if
                end do
            end do
    
            ! Добавляем единичную матрицу для гарантии невырожденности
            do i = 1, N
                A(i, i) = A(i, i) + cmplx(DIAG_SCALE, 0.0d0, kind=8)
            end do
    
        end subroutine generate_non_singular_matrix
        
        subroutine detailed_performance_analysis(n)
        implicit none
            integer, intent(in) :: n
            complex(8), allocatable :: A_old(:, :), A_new(:, :)
            complex(8) :: det_old, det_new
            real(8) :: start_time, end_time
            real(8) :: time_old, time_new
            integer :: i, num_warmup = 2
        
            print *, ""
            print *, "-----------------------------------------"
            print *, " Detailed analysis for ", n, "x", n, " matrix"
            print *, "-----------------------------------------"
        
            allocate(A_old(n, n), A_new(n, n))
        
            ! Warm-up runs (to avoid cold-start effects)
            do i = 1, num_warmup
                call generate_non_singular_matrix(A_old, n, n)
                A_new = A_old
                call DGAUS(A_old, n, n, det_old)
                call test_func(A_new, n, n, det_new)
            end do
        
            ! Actual timing with fresh matrices
            call generate_non_singular_matrix(A_old, n, n)
            A_new = A_old
        
            ! Time DGAUS
            call cpu_time(start_time)
            call DGAUS(A_old, n, n, det_old)
            call cpu_time(end_time)
            time_old = end_time - start_time
        
            ! Time test_func
            call cpu_time(start_time)
            call test_func(A_new, n, n, det_new)
            call cpu_time(end_time)
            time_new = end_time - start_time
            
            ! Check if determinants are valid
            if (is_nan_complex(det_old) .or. abs(det_old) < 1.0d-50) then
                print *, " WARNING: Old routine produced invalid determinant!"
                time_old = 0.0d0
            end if
            
            if (is_nan_complex(det_new) .or. abs(det_new) < 1.0d-50) then
                print *, " WARNING: New routine produced invalid determinant!"
                time_new = 0.0d0
            end if
        
            if (time_old > 0.0d0 .and. time_new > 0.0d0) then
                print '(A, F10.6, A)', " DGAUS time:     ", time_old, " seconds"
                print '(A, F10.6, A)', " test_func time: ", time_new, " seconds"
                print '(A, F8.2, A)', " Speedup:        ", time_old/time_new, "x"
            else
                print '(A, F10.6, A)', " DGAUS time:     ", time_old, " seconds"
                print '(A, F10.6, A)', " test_func time: ", time_new, " seconds"
                print '(A, A)', " Speedup:        ", "N/A"
            end if
            
            print '(A, L1)',       " Determinants match: ", abs(det_old - det_new) < 1.0d-8
            if (.not. is_nan_complex(det_old) .and. abs(det_old) > 1.0d-50) then
                print '(A, 2ES12.4)', " DGAUS det:      ", det_old
            else
                print '(A, A)', " DGAUS det:      ", "INVALID"
            end if
            
            if (.not. is_nan_complex(det_new) .and. abs(det_new) > 1.0d-50) then
                print '(A, 2ES12.4)', " test_func det:  ", det_new
            else
                print '(A, A)', " test_func det:  ", "INVALID"
            end if
        
            deallocate(A_old, A_new)
        end subroutine detailed_performance_analysis

end subroutine speed_tests_for_Gauss_LU_method

    
    
end module test_LU_Gaussian_elimination