!  Repulsion -- calculation of GW dispersion curves and   
!    frequency responses of a multilayered elastic waveguide 
!    for investigating repulsion and ZGV-resonance phenomena
 
    program Repulsion
    use Mult_Glob
    !use globals
    use RealPoles
    use integrate
    use from_zero_to_infinity_integrating
    use compare_integration_value_and_real_model
    
    use test_materials
    use test_gauss_integrating
    
    implicit none
  
    character(len=200) fortran_file, comsol_file, output_file
    complex(8) :: tst
    
    
    
	call Init_RP
    
    
    !tst = integrate_simpson(my_func, 2000, 0d0, 2d0)
    !print *, 'Integral = ', tst
    !pause
    
    !call plot_Ktrace_in_area()
    !call PolesAtAngle(1d0,1d0,1d0)
    !call plot_poles_at_angles()
             
    
    fortran_file = 'My_files\My_Resource_Files\integral_by_area.dat'
    
    comsol_file = 'My_files\My_Resource_Files\LineOnAnisotropicRectangle.txt'
    
    output_file = 'My_files\My_Resource_Files\merged_file.dat'
    
    !call plot_wavefield(fortran_file)
    !call test_material
    call run_tests_for_gauss()
    
    !call merge_data_files(fortran_file, comsol_file, output_file)
    
    pause
    !contains
    !    complex(8) function my_func(x) result(fx)
    !       complex(8), intent(in) :: x
    !       fx = exp(-x) / (x - 1d0)
    !   end function
    
	end 
