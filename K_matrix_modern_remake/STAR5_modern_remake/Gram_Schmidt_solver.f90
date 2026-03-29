module Gram_Schmidt_solver
implicit none

    contains
    
    subroutine get_solution_by_Gram_Schmidt_solver(A, B, C, R, ND, N, M, kort)
    implicit none
    
    ! Входные/выходные параметры
    integer(8), intent(in) :: ND, N, M, kort
    complex(8), intent(inout) :: A(ND, N)     ! Входная матрица
    complex(8), intent(inout) :: B(ND, M)     ! Правая часть/решение
    complex(8), intent(out) :: C(ND, N)       ! Ортонормированный базис
    complex(8), intent(out) :: R(ND, M+3)     ! Рабочая память
    
    ! Локальные переменные
    integer(8) :: i, k, k1, kw, jr, im, col
    integer(8) :: ms, msm, mt
    real(8) :: norm, residual_norm
    complex(8) :: proj, diag_elem
    
    ! Константы
    real(8), parameter :: EPS = 1.0D-13
    real(8), parameter :: TINY = 1.0D-12
    
    ! Индексы для рабочего массива R
    ms = M + 1
    msm = M + 2
    mt = M + 3
    
    !=========================================================================
    ! 1. ОРТОГОНАЛИЗАЦИЯ ГРАМА-ШМИДТА
    !    Создание ортонормированного базиса из столбцов матрицы A
    !=========================================================================
    
    do k = 1, N
        k1 = k - 1
        kw = 0
        
        ! Копируем k-й столбец A в рабочий массив R (в столбцы ms и msm)
        do i = 1, N
            R(i, ms) = A(i, k)
            R(i, msm) = A(i, k)
        enddo
        
        ! Вычисляем Евклидову норму столбца (корень из суммы попарно перемноженных элементов 'R(1:N,ms)' и 'conjg(R(1:N, ms))')
        norm = sqrt(sum(R(1:N, ms) * conjg(R(1:N, ms))))
        
        ! Первый столбец - просто нормируем
        if (k == 1) then
            R(1, mt) = norm
            C(1:N, 1) = R(1:N, ms) / norm
            cycle
        endif
        
        ! Для остальных столбцов: ортогонализация относительно предыдущих
        do jr = 1, k1
            ! Вычисляем проекцию на jr-й ортонормированный вектор
            proj = (0.0d0, 0.0d0)
            do i = 1, N
                proj = proj + R(i, ms) * conjg(C(i, jr))
            enddo
            
            ! Вычитаем проекцию
            do i = 1, N
                R(i, msm) = R(i, msm) - proj * C(i, jr)
            enddo
        enddo
        
        ! Вычисляем норму после ортогонализации
        residual_norm = 0.0d0
        do i = 1, N
            residual_norm = residual_norm + abs(R(i, msm))**2
        enddo
        residual_norm = sqrt(residual_norm)
        
        ! Проверка на вырожденность
        if (residual_norm > EPS) then
            kw = kw + 1
            ! Нормируем ортогональный вектор
            do i = 1, N
                R(i, ms) = R(i, msm) / residual_norm
                R(i, msm) = R(i, ms)  ! Копируем для следующей итерации
            enddo
            
            ! Повторяем ортогонализацию при необходимости
            if (kw < kort) cycle
            
            ! Вычисляем элемент для обратной подстановки
            proj = (0.0d0, 0.0d0)
            do i = 1, N
                proj = proj + R(i, ms) * conjg(A(i, k))
            enddo
            R(k, mt) = proj
            
            ! Сохраняем ортонормированный вектор
            C(1:N, k) = R(1:N, ms)
            
        else
            ! Вектор линейно зависим - обнуляем
            do i = 1, N
                R(i, ms) = (0.0d0, 0.0d0)
            enddo
            proj = (0.0d0, 0.0d0)
            R(k, mt) = proj
            C(1:N, k) = R(1:N, ms)
        endif
        
   enddo
    
    !=========================================================================
    ! 2. ПРЯМОЙ ХОД
    !    Копирование правых частей в рабочий массив
    !=========================================================================
    
    do im = 1, M
        R(1:N, im) = B(1:N, im)
    enddo
    
    !=========================================================================
    ! 3. ОБРАТНЫЙ ХОД
    !    Решение системы с использованием ортонормированного базиса
    !=========================================================================
    
    do col = N, 1, -1
        ! Загружаем col-й ортонормированный вектор
        R(1:N, ms) = C(1:N, col)
        
        diag_elem = R(col, mt)
        
        ! Если диагональный элемент слишком мал, решение = 0
        if (abs(diag_elem) <= TINY) then
            B(col, 1:M) = (0.0d0, 0.0d0)
            cycle
        endif
        
        ! Вычисляем решение для текущей строки
        do im = 1, M
            proj = (0.0d0, 0.0d0)
            do i = 1, N
                proj = proj + R(i, ms) * conjg(R(i, im))
            enddo
            B(col, im) = conjg(proj / diag_elem)
        enddo
        
        if (col == 1) then
            exit 
        endif
        
        ! Обновляем правые части для следующих итераций
        do i = 1, N
            do im = 1, M
                R(i, im) = R(i, im) - B(col, im) * A(i, col)
            enddo
        enddo
        
    enddo 
    
    return
    end subroutine get_solution_by_Gram_Schmidt_solver
    

end module Gram_Schmidt_solver
