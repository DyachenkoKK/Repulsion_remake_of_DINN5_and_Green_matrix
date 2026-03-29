module polygone
implicit none

    contains
    
    subroutine test_func(A, A_size, n, result_)
    implicit none
    
    integer(8), intent(in) :: n, A_size
    complex(8), intent(inout) :: A(A_size, A_size)
    complex(8), intent(out) :: result_
    
    integer(8) :: i, i1, j, j_max, k, l
    real(8) :: A_current, A_temp
    complex(8) :: tmp
    
    result_ = (1.0d0, 0.0d0)
    
    do i = 1, n
        j = i
        j_max = j
        A_current = abs(A(i, j))
        
        if (i /= n) then
            i1 = i + 1
            
            ! Search for main matrix element
            do l = i1, n
                A_temp = abs(A(i, l))
                if (A_temp > A_current) then
                    A_current = A_temp
                    j_max = l
                end if
            end do
            
            ! Checking matrix for degeneracy
            if (A_current <= 1.0d-20) then
                result_ = (0.0d0, 0.0d0)
                return
            end if
            
            ! Swap rows if the current and maximum elements are different
            if (j_max /= j) then
                do k = 1, n
                    tmp = A(k, j)
                    A(k, j) = A(k, j_max)
                    A(k, j_max) = tmp
                end do
                result_ = -result_
            end if
            
            ! Gaussian elimination
            do k = i1, n
                tmp = A(k, j) / A(i, j)      
                do l = i1, n
                    A(k, l) = A(k, l) - tmp * A(i, l)
                end do
            end do
        end if
        
        result_ = result_ * A(i, j)
        
    end do
    
    return
end subroutine test_func
    
    

end module polygone
