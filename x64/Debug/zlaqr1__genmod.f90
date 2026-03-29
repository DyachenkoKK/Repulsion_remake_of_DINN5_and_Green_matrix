        !COMPILER-GENERATED INTERFACE MODULE: Tue Feb 10 19:15:39 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZLAQR1__genmod
          INTERFACE 
            SUBROUTINE ZLAQR1(N,H,LDH,S1,S2,V)
              INTEGER(KIND=4) :: LDH
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: H(LDH,*)
              COMPLEX(KIND=8) :: S1
              COMPLEX(KIND=8) :: S2
              COMPLEX(KIND=8) :: V(*)
            END SUBROUTINE ZLAQR1
          END INTERFACE 
        END MODULE ZLAQR1__genmod
