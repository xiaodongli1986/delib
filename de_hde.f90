

!-----------------------------------------------------------------------
!------------------ BEGINNING OF THE MODULE Union2 ---------------------
!-----------------------------------------------------------------------
MODULE de_hde

USE de_settings
USE de_tools
USE de_types

IMPLICIT NONE
	
	PRIVATE

	PUBLIC :: de_hde_init, de_hde_inv_e, de_hde_rhode, de_hde_w
!	DOUBLE PRECISION, PRIVATE 	:: Omz, Orz, Obz, Odez, Okz, Ode_iz, Ezsq

  CONTAINS

!########################################################################
!########################################################################
!		PUBLIC FUNCTIONS 					#
!########################################################################
!########################################################################

  !------------------------------------------
  ! inv_e(z) = 1 / e(z)
  !------------------------------------------
  	DOUBLE PRECISION FUNCTION de_hde_inv_e(z)    
		INTEGER :: i, i1, i2, i3
		DOUBLE PRECISION :: z, z1, z2, z3, f1, f2, f3, Ez

		IF(z > 5.0d5) THEN
			Ez = sqrt(de_CP%Or0 * (1.0d0+z)**4.0d0)
			de_hde_inv_e = 1.0d0 / Ez
			RETURN
		ENDIF

		i = de_iz(z)

		IF(i < 2) i=2
		IF(i > de_num_intpl - 1) i = de_num_intpl-1

		i1 = i-1; i2 = i; i3 = i+1
		z1 = de_zdata(i1); z2 = de_zdata(i2);	z3 = de_zdata(i3)
		f1 = de_ezdata(i1); f2 = de_ezdata(i2);	f3 = de_ezdata(i3)

		Ez = de_intpl_vl(z, z1, f1, z2, f2, z3, f3)
		de_hde_inv_e = 1.0d0 / Ez
	END FUNCTION de_hde_inv_e


  !------------------------------------------
  ! rho_hde(z),
  !  propotional to hde energy density, normalized to 1 at z = 0
  !------------------------------------------	
	DOUBLE PRECISION FUNCTION de_hde_rhode(z)
		DOUBLE PRECISION :: z
		IF(z < 2000.0) THEN
			de_hde_rhode = de_get_odez(z) * de_get_ez(z)**2.0d0 / de_CP%Ode0
			RETURN
		ENDIF
		
		de_hde_rhode = 0.0d0
	END FUNCTION de_hde_rhode

  !------------------------------------------
  ! w_hde(z)
  !------------------------------------------	
	DOUBLE PRECISION FUNCTION de_hde_w(z)
		DOUBLE PRECISION :: z
		IF(z < 2000.0) THEN
			de_hde_w = - 1.0/3.0 - 2.0 / 3.0 * &
				sqrt(de_get_odez(z)/de_CP%hde%c**2.0 + de_CP%Ok0*(1.0d0+z)**2.0/(de_get_Ez(z)**2.0))
			RETURN
		ENDIF
		
		de_hde_w = - 1.0d0/3.0d0		
	END FUNCTION de_hde_w

  !------------------------------------------
  ! Initialize the HDE model
  !------------------------------------------
	SUBROUTINE de_hde_Init()
		INTEGER :: i
		DOUBLE PRECISION :: z1, Ez1, Odez1, z2, Ez2, Odez2

		i	= 1
		z1 	= 0.0d0
		Ez1 	= 1.0d0
		Odez1 	= de_CP%Ode0

		DO i = 1, de_num_intpl
			z2 = de_zi(i)

			CALL de_RK2(DEzDz, DOdezDz, z1, Ez1, Odez1, z2, Ez2, Odez2, 1)

			de_ezdata(i) = Ez2
			de_odezdata(i) = Odez2

			z1 = z2; Ez1 = Ez2; Odez1 = Odez2
		ENDDO		
	END SUBROUTINE de_hde_Init
	

!########################################################################
!########################################################################
!		PRIVATE FUNCTIONS 					#
!########################################################################
!########################################################################

  !----------------------------------------------------------------
  ! d E(z) / dz
  !----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION DEzDz(z, Ez, Odez)
		DOUBLE PRECISION :: z, Ez, Odez
		DOUBLE PRECISION :: Okz, Orz, Obz, Omz, Ode_iz
		DOUBLE PRECISION :: Ezsq
		Ezsq = Ez**2.0d0
		Okz  = de_CP%Ok0*(1.0d0+z)**2.0d0/Ezsq
		Orz  = de_CP%Or0*(1.0d0+z)**4.0d0/Ezsq
		Obz  = de_CP%Ob0*(1.0d0+z)**3.0d0/Ezsq
		Omz  = 1.0d0 - Okz - Orz - Obz - Odez
		Ode_iz = de_CP%hde%Gamma * Ez**(2.0d0*(de_CP%hde%a+de_CP%hde%b)-3.0d0)  * Omz**de_CP%hde%a * Odez**de_CP%hde%b
 
		DEzDz = (3.0d0 * Odez + Okz - Orz - 3.0d0 + Ode_iz) / (2.0d0*Odez)
		DEzDz = DEzDz - 1.0d0 + sqrt(Odez/de_CP%hde%c**2.0d0+Okz)
		DEzDz = - Odez/(1.0d0+z) * DEzdz
		DEzDz = DEzDz * Ez
	END FUNCTION DEzDz


  !----------------------------------------------------------------
  ! d \Omega_{de}(z) / dz
  !----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION DOdezDz(z, Ez, Odez)
		DOUBLE PRECISION :: z, Ez, Odez
		DOUBLE PRECISION :: Okz, Orz, Obz, Omz, Ode_iz
		DOUBLE PRECISION :: Ezsq
		Ezsq = Ez**2.0d0
		Okz  = de_CP%Ok0*(1.0d0+z)**2.0d0/Ezsq
		Orz  = de_CP%Or0*(1.0d0+z)**4.0d0/Ezsq
		Obz  = de_CP%Ob0*(1.0d0+z)**3.0d0/Ezsq
		Omz  = 1.0d0 - Okz - Orz - Obz - Odez
		Ode_iz = de_CP%hde%Gamma * Ez**(2.0d0*(de_CP%hde%a+de_CP%hde%b)-3.0d0)  * Omz**de_CP%hde%a * Odez**de_CP%hde%b
 
		DOdezDz = sqrt(Odez/de_CP%hde%c**2.0d0+Okz) - 1.0d0
		DOdezDz = DOdezDz - (3.0d0*Odez-Orz+Okz-3.0d0+Ode_iz) / (2.0d0 - 2.0d0*Odez)
		DOdezDz = - (2.0d0*Odez*(1.0d0-Odez)/(1.0d0+z)) * DOdezDz 
	END FUNCTION DOdezDz  
END MODULE de_hde
