subroutine K_amp
use Mult_Glob
implicit none
integer i, j, alfa1_num, alfa2_num
real*8 alfa1_min, alfa1_max, alfa1_step, alfa2_min, alfa2_max, alfa2_step, alfa1, alfa2, gamma, alfa
complex*16 alfa_c
    print*, "K_amp has been started "
    f = 0.01d0; w=f*2d0*pi;
    
    alfa1_num = 1000; alfa2_num = 1000;
    alfa1_min = -0.06d0; alfa2_min = -0.06d0;
    alfa1_max = 0.06d0; alfa2_max = 0.06d0;
    key=1
    
    alfa1_step = (alfa1_max - alfa1_min)/alfa1_num; alfa2_step = (alfa2_max - alfa2_min)/alfa2_num;
    open(unit=1033,file='K_amp.dat',form ='formatted')
    do i = 1, alfa1_num
        do j = 1, alfa2_num
            alfa1 = alfa1_min + alfa1_step*i; alfa2 = alfa2_min + alfa2_step*j;
            alfa = sqrt(alfa1**2+alfa2**2); alfa_c = alfa*(1d0, 0d0); gamma = atan(alfa2, alfa1);
            
            call MultiK_An(alfa_c, gamma, 0d0)
            write(1033, *) alfa1, alfa2, abs(Kaz(1,1))
        enddo
    enddo
    close(1033)      
end
    
subroutine K_pol
use Mult_Glob
implicit none

integer j, alfa1_num, alfa2_num
real*8 alfa1_min, alfa1_max, alfa1_step, alfa2_min, alfa2_max, alfa2_step, alfa1, alfa2, alfa

integer i, gm_num, Ndz
real*8 gm_step, gamma, alfa_max, pols(10), RDabs_alf, gm
complex*16 alfa_c
common/gm/gm
external RDabs,RDabs_alf


    print*, "K_pol has been started "
    f = 0.01d0; w=f*2d0*pi;
    
    gm_num = 500; alfa_max = 10d0;
    gm_step = 1.999d0*pi/gm_num;
    open(unit=1033,file='K_pol.dat',form ='formatted')
    do i = 1, gm_num
        
           
        !alfa1_num = 1000; alfa2_num = 1000;
        !alfa1_min = -0.06d0; alfa2_min = -0.06d0;
        !alfa1_max = 0.06d0; alfa2_max = 0.06d0;
        !alfa1 = alfa1_min + alfa1_step*i; alfa2 = alfa2_min + alfa2_step*j;
        !alfa = sqrt(alfa1**2+alfa2**2); alfa_c = alfa*(1d0, 0d0); gamma = atan(alfa2, alfa1);
        !    
        !call MultiK_An(alfa_c, gamma, 0d0)
        !
        !write(1033, *) abs(Kaz(1,1))
        !
        
        
        gm = gm_step*i
        call Hamin(RDabs_alf, 1d-3, alfa_max, 1d-2, 1d-6, 10, pols, Ndz)
        if (Ndz>0) then
            write(1033, *) gm, pols(1)
        else 
            print*, 'No poles!'
        endif
        
    enddo
    
    close(1033)
        
end