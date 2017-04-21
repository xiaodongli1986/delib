!######################################################################################
!This module choose the dark energy model for calculation.
!
!For your own dark energy model, you shall give:
!
! A. An init subroutine
! B. Its expansion history inv_e(z).
! C. Its normalized dark energy density rho_de(z)
! D. Its equation-of-state w(z)
!
!To caluate the chisqs defined in the chisqs.f90, only A and B are necessary.
!To play with CAMB, C and D are required.
!
!Please define these functions/subroutins in your MODELNAME.f90,
!and add thie file into the Makefile.
!
!Please give a name of the model, and define 
!
!######################################################################################

MODULE de_model_init

USE de_settings
USE de_tools
USE de_types
USE de_hde
USE de_wcdm3
USE de_srom
USE de_ICG
USE de_qz
USE de_mauricehde

IMPLICIT NONE

	!
	!Please tell the program which model shall be initied
	!
	INTEGER, PARAMETER :: de_lcdm_lab = 1
	INTEGER, PARAMETER :: de_wcdm_lab = 2
	INTEGER, PARAMETER :: de_CPL_lab =  3
	INTEGER, PARAMETER :: de_hde_lab =  4 !Holographic dark energy
	INTEGER, PARAMETER :: de_wcdm3_lab = 5 !binned w(z), 3 bins
	INTEGER, PARAMETER :: de_srom_lab = 6
	INTEGER, PARAMETER :: de_ICG_lab = 7
	INTEGER, PARAMETER :: de_qz_lab = 8 ! deceleration parameter as a function of redshift
	INTEGER, PARAMETER :: de_Rhct_lab = 9 ! H(z) = (1+z)^alpha H0
	INTEGER, PARAMETER :: de_mauricehde_lab = 10 ! H(z) = (1+z)^alpha H0
	
	
	!model to be used
!	INTEGER :: de_model_lab = de_hde_lab
	INTEGER :: de_model_lab = de_wcdm3_lab
	
	LOGICAL :: de_inited = .FALSE.

	CONTAINS


  !------------------------------------------
  ! Initialize the Calculation
  !------------------------------------------
	SUBROUTINE de_init()
		INTEGER :: i
		DOUBLE PRECISION :: b1, b2, g1, g2
		DOUBLE PRECISION :: z, x, lnH, F ! for test

		IF(.not. de_tools_inited) THEN
			CALL de_tools_init()
			de_tools_inited = .TRUE.
		ENDIF

		!Common initialization for all models.
		!Assume that Ob0hsq, h, Odm0, alpha, beta in de_CP are known.
		de_CP%H0	= de_CP%h*100.0d0
	
		de_CP%Ob0	= de_CP%Ob0hsq / de_CP%h**2
		de_CP%Odm0hsq 	= de_CP%Odm0 * de_CP%h**2
		de_CP%Om0    	= de_CP%Odm0+de_CP%Ob0
		de_CP%Om0hsq 	= de_CP%Om0*de_CP%h**2.0d0
		de_CP%Og0	= 2.469d-5 / de_CP%h**2
		de_CP%Or0 	= de_CP%Og0 * (1.0d0 + 0.2271d0 * 3.04d0)

		de_CP%Ode0   	= 1.0d0 - de_CP%Odm0 - de_CP%Ok0 - de_CP%Ob0 - de_CP%Or0

		IF(de_model_lab .EQ. de_lcdm_lab .or. &
			de_model_lab .EQ. de_wcdm_lab .or.&
			de_model_lab .EQ. de_cpl_lab .or. &
			de_model_lab .EQ. de_Rhct_lab) THEN
			CONTINUE
			!No intialization need to be done
		ELSEIF(de_model_lab .EQ. de_hde_lab) THEN
			CALL de_hde_init()
		ELSEIF(de_model_lab .EQ. de_wcdm3_lab) THEN
			CALL de_wcdm3_init()
		ELSEIF(de_model_lab .EQ. de_srom_lab) then
			CALL de_srom_init()
		ELSEIF(de_model_lab .EQ. de_ICG_lab) then
			CALL de_ICG_init()
		ELSEIF(de_model_lab .EQ. de_qz_lab) then
			CALL de_qz_init()
		ELSEIF(de_model_lab .EQ. de_mauricehde_lab) then
			CALL de_mauricehde_init()
		ELSE
			WRITE(*,*) "Error! Initializaton of model ", de_model_lab, "not found!"
			STOP
		ENDIF
		
		!calculating values of zd, rs(zd)		
		b1 = 0.313d0*de_CP%Om0hsq**(-0.419d0)*(1.0d0+0.607d0*de_CP%Om0hsq**0.674d0)
		b2 = 0.238d0*de_CP%Om0hsq**(0.223d0)

		de_CP%zd = 1291.0d0*de_CP%Om0hsq**(0.251d0)* (1.0d0+b1*de_CP%Ob0hsq**b2) / (1.0d0+0.659d0 *de_CP%Om0hsq**(0.828d0)) 
		de_CP%rszd = de_rs(de_CP%zd)

		!calculating values of zstar, DAzstar, rszstar, R, lA		
		g1 = 0.0783d0*de_CP%Ob0hsq**(-0.238d0) / (1.0d0+39.5d0*de_CP%Ob0hsq**(0.763d0));
		g2 = 0.56d0 / (1.0d0+21.1d0*de_CP%Ob0hsq**(1.81d0));   
		de_CP%zstar = 1048d0*(1.0d0+0.00124d0*de_CP%Ob0hsq**(-0.738d0))*(1.0d0+g1*((de_CP%Odm0+de_CP%Ob0)*de_CP%h**2.0d0)**g2);
		de_CP%DAzstar = de_DA(de_CP%zstar)
		de_CP%R       = sqrt(de_CP%Odm0+de_CP%Ob0)*de_CP%H0*(1.0d0+de_CP%zstar)*de_CP%DAzstar
		de_CP%rszstar = de_rs(de_CP%zstar)
		de_CP%lA      = (1.0d0+de_CP%zstar)*de_Pi*de_CP%DAzstar/de_CP%rszstar
		
		de_inited = .TRUE. 

!		WRITE(*,*) de_CP%Ob0hsq, de_CP%Ob0, de_CP%Om0, de_CP%Ode0, de_CP%h
!		WRITE(*,*) Inv_e(1.0d0), Inv_e(2.0d0)		
!		WRITE(*,*) g(1.0d-3), g(2.0d-3)
!		WRITE(*,*) rs(1.0d3), rs(2.0d3)
	END SUBROUTINE de_init
	

  !------------------------------------------
  ! inv_e(z) = 1 / e(z)
  !------------------------------------------
  	DOUBLE PRECISION FUNCTION de_intpl_inv_e(z)    
		INTEGER :: i, i1, i2, i3
		DOUBLE PRECISION :: z, z1, z2, z3, f1, f2, f3, Ez

		IF(z > 5.0d8) THEN
			Ez = sqrt(de_CP%Or0 * (1.0d0+z)**4.0d0)
			de_intpl_inv_e = 1.0d0 / Ez
			RETURN
		ENDIF

		i = de_iz(z)

		IF(i < 2) i=2
		IF(i > de_num_intpl - 1) i = de_num_intpl-1

		i1 = i-1; i2 = i; i3 = i+1
		z1 = de_zdata(i1); z2 = de_zdata(i2);	z3 = de_zdata(i3)
		f1 = de_ezdata(i1); f2 = de_ezdata(i2);	f3 = de_ezdata(i3)

		Ez = de_intpl_vl(z, z1, f1, z2, f2, z3, f3)
		de_intpl_inv_e = 1.0d0 / Ez
	END FUNCTION de_intpl_inv_e

	
  !------------------------------------------
  ! inv_e(z) = 1 / e(z)
  !------------------------------------------
	DOUBLE PRECISION FUNCTION de_inv_e(z)    
		DOUBLE PRECISION z

		IF(de_model_lab .EQ. de_lcdm_lab) THEN
			de_inv_e = 1.0d0 / sqrt(de_CP%Om0*(1.0+z)**3.0  &
				+ de_CP%Ok0*(1.0+z)**2.0 &
				+ de_CP%Or0*(1.0+z)**4.0 &
				+ de_CP%Ode0)
		ELSEIF(de_model_lab .EQ. de_wcdm_lab) THEN
			de_inv_e = 1.0d0 / sqrt(de_CP%Om0*(1.0+z)**3.0  &
				+ de_CP%Ok0*(1.0+z)**2.0 &
				+ de_CP%Or0*(1.0+z)**4.0 &
				+ de_CP%Ode0*(1.0+z)**(3.0*(1.0+de_CP%wcdm%w)))
		ELSEIF(de_model_lab .EQ. de_CPL_lab) THEN
			de_inv_e = 1.0d0 / sqrt(de_CP%Om0*(1.0+z)**3.0  &
				+ de_CP%Ok0*(1.0+z)**2.0 &
				+ de_CP%Or0*(1.0+z)**4.0 &
				+ de_CP%Ode0*(1.0+z)**(3.0*(1.0+de_CP%CPL%w0+de_CP%CPL%wa)) * exp(-3.0*z*de_CP%CPL%wa/(1.0d0+z)))
		ELSEIF(de_model_lab .EQ. de_Rhct_lab) then
			de_inv_e = 1.0d0 / sqrt((1.0+z)**de_CP%Rhct%alpha)
		ELSEIF(de_model_lab .EQ. de_hde_lab) THEN
			de_inv_e = de_hde_inv_e(z)	
		ELSEIF(de_model_lab .EQ. de_wcdm3_lab) THEN
			de_inv_e = de_wcdm3_inv_e(z)
		ELSEIF(de_model_lab .EQ. de_srom_lab .OR. de_model_lab .EQ. de_ICG_lab .OR. de_model_lab .EQ. de_qz_lab &
		       .OR. de_model_lab .EQ. de_mauricehde_lab) THEN
			de_inv_e = de_intpl_inv_e(z)
		ELSE
			WRITE(*,*) "Error! inv_e(z) of model ", de_model_lab, "not found!"
			STOP
		ENDIF
	END FUNCTION de_inv_e
	
	
  !------------------------------------------
  ! rho_de(z),
  !  propotional to de energy density, 
  !  normalized to 1 at z = 0
  !------------------------------------------	
	DOUBLE PRECISION FUNCTION de_rhode(z)
		DOUBLE PRECISION :: z
		
		IF(de_model_lab .EQ. de_lcdm_lab) THEN
			de_rhode = 1.0d0
		ELSEIF(de_model_lab .EQ. de_wcdm_lab) THEN
			de_rhode = (1.0+z)**(3.0*(1.0+de_CP%wcdm%w))
		ELSEIF(de_model_lab .EQ. de_CPL_lab) THEN
			de_rhode = (1.0+z)**(3.0*(1.0+de_CP%CPL%w0+de_CP%CPL%wa)) * exp(-3.0*z*de_CP%CPL%wa/(1.0d0+z))
		ELSEIF(de_model_lab .EQ. de_hde_lab) THEN
			de_rhode = de_hde_rhode(z)
		ELSEIF(de_model_lab .EQ. de_wcdm3_lab) THEN
			de_rhode = de_wcdm3_rhode(z)
		ELSEIF(de_model_lab .EQ. de_qz_lab .or. de_model_lab .EQ. de_Rhct_lab &
			.or. de_model_lab .EQ. de_mauricehde_lab) then
			print *, 'ERROR! No rhode(z) available for q(z), Rhct model!'
			stop
		ELSE	
			WRITE(*,*) "Error! rho_de(z) of model ", de_model_lab, "not found!"
			STOP
		ENDIF
	END FUNCTION de_rhode

  !------------------------------------------
  ! w_hde(z)
  !------------------------------------------	
	DOUBLE PRECISION FUNCTION de_w(z)
		DOUBLE PRECISION :: z

		IF(de_model_lab .EQ. de_lcdm_lab) THEN
			de_w = 1.0d0
		ELSEIF(de_model_lab .EQ. de_wcdm_lab) THEN
			de_w = de_CP%wcdm%w
		ELSEIF(de_model_lab .EQ. de_CPL_lab) THEN
			de_w = de_CP%CPL%w0 + de_CP%CPL%wa * z / (1.0d0+z)
		ELSEIF(de_model_lab .EQ. de_hde_lab) THEN
			de_w = de_hde_w(z)
		ELSEIF(de_model_lab .EQ. de_wcdm3_lab) THEN
			de_w = de_wcdm3_w(z)
		ELSEIF(de_model_lab .EQ. de_qz_lab .or. de_model_lab .EQ. de_Rhct_lab .or.  &
			de_model_lab .EQ. de_mauricehde_lab) then
			print *, 'ERROR! No w(z) available for q(z), Rhct model!'
			stop
		ELSE
			WRITE(*,*) "Error! w(z) of model ", de_model_lab, "not found!"
			STOP		
		ENDIF
	END FUNCTION de_w



!!!!!!!!!! Useful Functions (Basically independent with DE model) !!!!!!!!!!
		
  !----------------------------------------------------------------
  ! \int ^ zright _ 0 f(z) dz
  !----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION de_Inte(zright)
		DOUBLE PRECISION :: zright
		DOUBLE PRECISION :: z1,z2,f1,f2
		DOUBLE PRECISION :: aright
		INTEGER :: i, N
		z1 = 0.0d0
		z2 = 0.5d0
		de_Inte = 0.0d0
		N  = 32
		DO WHILE(z2<zright)
			de_Inte = de_Inte + de_Simpson(de_inv_e,z1,z2,N)
			z1 = z2
			z2 = z2*2.d0
		ENDDO
		de_Inte = de_Inte + de_Simpson(de_inv_e,z1,zright,N)
	END FUNCTION de_Inte
   

  !----------------------------------------------------------------
  ! The FUNCTION to be integrated by rs(z)
  !----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION de_g(a)
		DOUBLE PRECISION :: a, z
		z = 1.0d0/a - 1.0d0
		de_g = de_inv_e(z)/( a**2 * sqrt( 1.0d0+3.0d0*de_CP%Ob0*a/(4.0d0*de_CP%Og0) ) )
	END FUNCTION de_g

   
  !----------------------------------------------------------------
  ! rs(z, vars)
  !----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION de_rs(z)
		DOUBLE PRECISION :: z
		de_rs = de_Simpson(de_g,1.0d-10,1.0d0/(1.0d0+z), 128) / ( sqrt(3.0d0) * de_CP%H0 )
	END FUNCTION de_rs 
	
	
  !----------------------------------------------------------------
  ! rstoDv(z, vars)
  !----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION de_rstoDv(z)
		DOUBLE PRECISION :: z
		de_rstoDV = de_CP%rszd / de_DV(z)
	END FUNCTION de_rstoDV

  !----------------------------------------------------------------
  ! DA/rd
  !----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION de_DAtord(z)
		DOUBLE PRECISION :: z
		de_DAtord = de_DA(z)/de_CP%rszd
	END FUNCTION de_DAtord


  !----------------------------------------------------------------
  ! H*rd
  !----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION de_Hrd(z)
		DOUBLE PRECISION :: z
		de_Hrd = (de_CP%H0/de_inv_e(z))*(de_CP%rszd) 
	END FUNCTION de_Hrd


  !----------------------------------------------------------------
  !renamed as de_DVtors
  !----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION de_DVtors(z)
		DOUBLE PRECISION :: z
		de_DVtors = de_DV(z)/de_CP%rszd
	END FUNCTION de_DVtors

  !------------------------------------------
  ! fk(vars, distance)
  !------------------------------------------
	DOUBLE PRECISION FUNCTION de_fk(distance)
		DOUBLE PRECISION :: distance, sqrtOk0
		sqrtOk0		= sqrt(ABS(de_CP%Ok0))
		IF(de_CP%Ok0 < 0.0d0) THEN
			de_fk = sin(sqrtOk0*distance)/sqrtOk0
			ELSEIF(de_CP%Ok0 > 0.0d0) THEN
			de_fk = sinh(sqrtOk0*distance)/sqrtOk0
			ELSE
			de_fk = distance
		ENDIF
	END FUNCTION de_fk
	
  !----------------------------------------------------------------
  ! Angular Diameter Distance
  !----------------------------------------------------------------	
	DOUBLE PRECISION FUNCTION de_DA(z)
		DOUBLE PRECISION :: z
		de_DA = de_fk(de_Inte(z)) / (de_CP%H0*(1.0d0+z))
	END FUNCTION de_DA


  !----------------------------------------------------------------
  ! Function DV
  !----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION de_DV(z)
		DOUBLE PRECISION :: z
		de_DV = ((1.0d0+z)**2 * de_DA(z)**2 * z * de_inv_e(z) / de_CP%H0)**(1.0d0/3.0d0)
	END FUNCTION de_DV

  !----------------------------------------------------------------
  ! A parameter in the bao data
  !----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION de_A_bao(z)
		DOUBLE PRECISION :: z
		de_A_bao = de_CP%H0 * de_DV(z) * sqrt(de_CP%Odm0+de_CP%Ob0) / z
	END FUNCTION de_A_bao
END MODULE de_model_init
