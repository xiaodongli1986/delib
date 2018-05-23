
!This is a test for your code.
!
!Just type these commands,
!
! ifort test5.f90 -mkl -$delm 
! ./a.out

program main

use de_model_init
use de_chisqs

implicit none

	double precision :: y, y2
 	double precision :: omegam, z, a, ez, qz, q_ez_fun, H_residual, amin
 	integer :: i,j,k, ia
 
 	pr_chisq_info = .false.
 
!#############################################
! Compute rszd
if(.true.) then
    	do i = 1, 999
        de_model_lab = de_lcdm_lab
        de_CP%Ob0hsq    =  0.02204530992
        omegam = (i+0.5)*0.001
        de_CP%h         =  0.6777
        de_CP%alpha     =  0.142358E+01
        de_CP%beta      =  0.325629E+01
        de_CP%Odm0      =  omegam - de_CP%Ob0hsq/de_CP%h**2.0
	call de_init()
	y = de_chisq_sdssdr12_APinlcdm(dble(omegam))
	print *, omegam, y
	enddo
endif

if(.false.) then
        de_model_lab = de_lcdm_lab
        de_CP%Ob0hsq    =  0.022d0
        omegam = 0.31d0
        de_CP%h         =  0.676d0
        de_CP%alpha     =  0.142358E+01
        de_CP%beta      =  0.325629E+01
        de_CP%Odm0      =  omegam - de_CP%Ob0hsq/de_CP%h**2.0
	call de_init()
	print *, 'rszd = ', de_CP%rszd*de_const_c
	print *, 'rszd / rszd(from camb) = ', de_CP%rszd*de_const_c / 147.78d0
endif

if(.true.) then
        de_model_lab = de_lcdm_lab
        de_CP%Ob0hsq    =  0.022d0
        omegam = 0.31d0
        de_CP%h         =  0.676d0
        de_CP%alpha     =  0.142358E+01
        de_CP%beta      =  0.325629E+01
        de_CP%Odm0      =  omegam - de_CP%Ob0hsq/de_CP%h**2.0
	call de_init()
	print *, 'rszd = ', de_CP%rszd*de_const_c
	print *, 'rszd / rszd(from camb) = ', de_CP%rszd*de_const_c / 147.78d0
endif

end program
