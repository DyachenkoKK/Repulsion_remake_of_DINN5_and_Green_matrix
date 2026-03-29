        !COMPILER-GENERATED INTERFACE MODULE: Tue Feb 10 19:15:36 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EV_EVALUATE__genmod
          INTERFACE 
            SUBROUTINE EV_EVALUATE(A,N,E_VALUES,E_VECTORS)
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: A(N,N)
              COMPLEX(KIND=8) :: E_VALUES(N)
              COMPLEX(KIND=8) :: E_VECTORS(N,N)
            END SUBROUTINE EV_EVALUATE
          END INTERFACE 
        END MODULE EV_EVALUATE__genmod
