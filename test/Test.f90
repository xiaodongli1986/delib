
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
 	double precision :: omegam, z, a, ez, qz, q_ez_fun, H_residual
 	integer :: i, ia
 
 	pr_chisq_info = .false.
 
!==========================================
! LCDM
!==========================================
if(.true.) then
	print *
	print *
	de_model_lab = de_lcdm_lab
	
	de_CP%Ob0hsq    =  0.02253
	omegam 		= 0.264936E+00
	de_CP%h		=  0.711833E+00
	de_CP%alpha	=  0.142125E+01
	de_CP%beta	=  0.325121E+01  
	de_CP%Odm0 	=  omegam - de_CP%Ob0hsq/de_CP%h**2.0
	de_CP%Ok0	= 0
	
	call de_init()

!	y = de_chisq_g06(.false., .true.)
!	y = de_chisq_jla()
	y = de_chisq_snls3() + de_chisq_sdssdr7_old() + de_chisq_wmap7() + de_chisq_h_Riess()
	
	write(*,*) ""
	write(*,*) "====================================="
	write(*,*) "  Resut of base flat LCDM:"
	write(*,*) "     Total chisq (lnlike) = ", y, y/2.0
	write(*,*) "     Expected chisq = ", 424.911450626633d0
	write(*,*) "====================================="	

!	stop

!==========================================
! XCDM
!==========================================

!                0.264213E+00
!               -0.108350E+01
!                0.723622E+00
!                0.142359E+01
!                0.325647E+01

	print *
	print *
	de_model_lab = de_wcdm_lab

	de_CP%Ob0hsq    =  0.02253
	omegam 		=  0.264213E+00
	de_CP%wcdm%w	=  -0.108350E+01
	de_CP%h		=  0.723622E+00
	de_CP%alpha	=  0.142359E+01
	de_CP%beta	=  0.325647E+01
	de_CP%Odm0 	=  omegam - de_CP%Ob0hsq/de_CP%h**2.0
	de_CP%Ok0	= 0
	
	call de_init()

	y = de_chisq_snls3() + de_chisq_sdssdr7_old() + de_chisq_wmap7() + de_chisq_h_Riess()
	
	write(*,*) ""
	write(*,*) "====================================="
	write(*,*) "  Resut of base flat XCDM:"
	write(*,*) "     Total chisq (lnlike) = ", y, y/2.0
	write(*,*) "     Expected chisq = ", 423.444613315747 
	write(*,*) "====================================="


!==========================================
! CPL
!==========================================

!                0.264624E+00
!               -0.106659E+01
!               -0.911151E-01
!                0.724149E+00
!                0.142358E+01
!                0.325629E+01
	
	print *
	print *
	de_model_lab = de_CPL_lab
	de_CP%Ob0hsq    =  0.02253
	omegam = 0.264624E+00
	de_CP%CPL%w0 =  -0.106659E+01
	de_CP%CPL%wa = 	-0.911151E-01
	de_CP%h		=  0.724149E+00
	de_CP%alpha	=  0.142358E+01
	de_CP%beta	=  0.325629E+01  
	de_CP%Odm0 	=  omegam - de_CP%Ob0hsq/de_CP%h**2.0
	de_CP%Ok0	= 0

	call de_init()
	y = de_chisq_snls3() + de_chisq_impwig()+de_chisq_dr11()+de_chisq_6dfGS() + de_chisq_wmap9() + de_chisq_planck() + de_chisq_h_Riess()
	
	write(*,*) ""
	write(*,*) "====================================="
	write(*,*) "  Resut of flat CPL:"
	write(*,*) "     Total chisq (lnlike) = ", y, y/2.0
	write(*,*) "     Expected chisq = ", 448.620325726934
	write(*,*) "====================================="


!==========================================
! srom
!==========================================

	print *
        print *
        pr_chisq_info = .false.
	de_model_lab = de_srom_lab
	de_CP%Ob0hsq    =  0.02253d0
	de_CP%Odm0 	=  0.26d0
	de_CP%srom%w0 =  -1.0d0
	de_CP%srom%w1 =   0.0d0
	de_CP%srom%xi =	 3.0d0
	de_CP%h		=  0.7d0
	de_CP%alpha	=  0.142358E+01
	de_CP%beta	=  0.325629E+01  
	de_CP%Ok0	= 0.0d0

	de_model_lab = de_srom_lab
	call de_init()
	y = de_chisq_snls3() + de_chisq_impwig()+de_chisq_dr11()+de_chisq_6dfGS() + de_chisq_wmap9() + de_chisq_planck() + de_chisq_h_Riess()

	write(*,*) ""
	write(*,*) "====================================="
	write(*,*) "  Resut of srom:"
	write(*,*) "     Total chisq (lnlike) = ", y, y/2.0
	write(*,*) "     Expected chisq = ", 1002.72980461404
	write(*,*) "====================================="

	do i = 1, 5
		z = i*10000		
		print *, 'z, H(z) (in SROM) = ', z, de_CP%H0/de_inv_e(z)
		de_model_lab = de_lcdm_lab
		call de_init()
		print *, 'z, H(z) (in LCDM) = ', z, de_CP%H0/de_inv_e(z)
		de_model_lab = de_srom_lab
		call de_init()
	enddo
	
!==========================================
! ICG	
!==========================================
	
	print *
        print *
        pr_chisq_info = .false.
	de_model_lab = de_ICG_lab
	de_CP%Ob0hsq    =  0.02253d0
	de_CP%Odm0 	=  0.26d0
	de_CP%h		=  0.7d0
	de_CP%alpha	=  0.142358E+01
	de_CP%beta	=  0.325629E+01  
	de_CP%Ok0	= 0.0d0

	de_CP%ICG%xi =  3.0d0
	de_CP%ICG%A =   (1.0-de_CP%Odm0-de_CP%Ob0hsq/de_CP%h**2.0 - de_CP%Or0)**2.0
	call de_init()
	de_CP%ICG%A =   (1.0-de_CP%Odm0-de_CP%Ob0hsq/de_CP%h**2.0 - de_CP%Or0)**2.0	
	call de_init()
	
	y = de_chisq_snls3() + de_chisq_impwig()+de_chisq_dr11()+de_chisq_6dfGS() + de_chisq_wmap9() + de_chisq_planck() + de_chisq_h_Riess()

	write(*,*) "====================================="
	write(*,*) "  Resut of ICG:"
	write(*,*) "     Total chisq (lnlike) = ", y, y/2.0
	write(*,*) "     Expected chisq = ", 1002.72980461404
	write(*,*) "====================================="

	do i = 1, 5
		z = i*10000		
		print *, 'z, H(z) (in ICG ) = ', z, de_CP%H0/de_inv_e(z)
		de_model_lab = de_lcdm_lab
		call de_init()
		print *, 'z, H(z) (in LCDM) = ', z, de_CP%H0/de_inv_e(z)
		de_model_lab = de_ICG_lab
		call de_init()
	enddo

!==========================================
! qz
!==========================================
	
	print *
        print *
        pr_chisq_info = .false.
	de_model_lab = de_qz_lab
	de_CP%Ob0hsq    =  0.02253d0
	de_CP%Odm0 	=  0.26d0
	de_CP%h		=  0.7d0
	de_CP%alpha	=  0.142358E+01
	de_CP%beta	=  0.325629E+01  
	de_CP%Ok0	= 0.0d0

	de_CP%qz%q0 =  -0.8d0
	de_CP%qz%q1 =  2.0d0
	call de_init()
	
	y = de_chisq_snls3() !+ de_chisq_impwig()+de_chisq_dr11()!+de_chisq_6dfGS() !+ de_chisq_wmap9() + de_chisq_planck() + de_chisq_h_Riess()

	write(*,*) "====================================="
	write(*,*) "  Resut of qz:"
	write(*,*) "     Total chisq (lnlike) = ", y, y/2.0
	write(*,*) "     Expected chisq = ", 1002.72980461404
	write(*,*) "====================================="

	! Check of union2p1 zc (z cutted)
	if (.false.) then
		do i = 15, 1, -3
			print *
			z = i*0.1
	!		y = de_chisq_snls3()/2.0
	!		print *, 'z<', z, '; snls3 like = ', y

			y = de_chisq_union2p1()/2.0
			print *, 'z<', z, '; union2p1 like = ', y

			y = de_chisq_union2p1zcs(z, 2.0d0)/2.0
			print *, 'z<', z, '; union2p1 zcs like = ', y
			y = de_chisq_union2p1zc(z)/2.0
			print *, 'z<', z, '; union2p1 zc like = ', y
		enddo
	endif
endif

!==========================================
! de_mauricehde
!==========================================

	print *
	print *
	
        DO i = 1, 30
	de_model_lab = de_mauricehde_lab
	de_CP%Ob0hsq    =  0.02253
	omegam 		=  0.26+i*0.002!0.284936E+00
	!de_CP%h		=  0.711833E+00
	de_CP%h		=  0.71
	de_CP%alpha	=  0.142125E+01
	de_CP%beta	=  0.325121E+01  
	de_CP%Odm0 	=  omegam - de_CP%Ob0hsq/de_CP%h**2.0
	de_CP%Ok0	= 0

	call de_init()

!	y = de_chisq_g06(.false., .true.)
!	y = de_chisq_jla()
!	y = de_chisq_snls3() !+ de_chisq_sdssdr7_old() + de_chisq_wmap7() + de_chisq_h_Riess()
	y = de_chisq_union2p1() !+ de_chisq_wmap7()
	y = de_chisq_wmap7()
	y = de_chisq_planck()

	open(unit=987,file='qz.txt')
	do ia = 1, 150
		z=de_zi(ia*100+1)
		!z = ia * 0.1d0
		a=1.0/(1.0+z)
		ez=1.0/de_inv_e(z)
		qz=de_mauricehde_q(a,ez)
		q_ez_fun=(1-qz)*ez*ez
		H_residual=de_CP%Om0/a**3 + de_CP%Or0/a**4 + (1-qz)*ez*ez / 3.0 - ez*ez
		!print *, 'z / h / q = ', real(z), real(ez), real(de_mauricehde_q(a,ez))
		write(987,*)  real(z), real(ez), real(qz), real(q_ez_fun), real(H_residual)
	enddo
	close(987)
	
	de_model_lab = de_lcdm_lab
	call de_init()
!	y2 = de_chisq_snls3() !+ de_chisq_sdssdr7_old() + de_chisq_wmap7() + de_chisq_h_Riess()
	y2 = de_chisq_union2p1() !+ de_chisq_wmap7()
	y2 = de_chisq_wmap7()
	y2 = de_chisq_planck()

	open(unit=987,file='qz_lcdm.txt')
	do ia = 1, 150
		z=de_zi(ia*100+1)
		!z = ia * 0.1d0
		a=1.0/(1.0+z)
		ez=1.0/de_inv_e(z)
		qz= -1.0 - (a/ez) * (-3*de_CP%Om0*a**(-4) -4*de_CP%Or0*a**(-5))/(2*ez)
		q_ez_fun=(1-qz)*ez*ez
		!print *, 'z / h / q = ', real(z), real(ez), real(qz)
		write(987,*)  real(z), real(ez), real(qz), real(q_ez_fun)
	enddo
	close(987)

	print *, 'omegam /chisqs of HDE/LambdaCDM =', real(omegam), real(y), real(y2)
	ENDDO
	write(*,*) ""
	write(*,*) "====================================="
	write(*,*) "  Resut of Maurice's HDE:"
	write(*,*) "     Total chisq (lnlike) = ", y, y/2.0
	write(*,*) "     Expected chisq = ", 424.911450626633d0
	write(*,*) "====================================="	

end program main
