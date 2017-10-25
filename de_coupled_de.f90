

!-----------------------------------------------------------------------
!------------------ BEGINNING OF THE MODULE Union2 ---------------------
!-----------------------------------------------------------------------
MODULE de_mauricehde

USE de_settings
USE de_tools
USE de_types

IMPLICIT NONE
	
	PRIVATE

	PUBLIC :: de_mauricehde_init, de_mauricehde_dezda, de_mauricehde_q !, de_srom_init, de_srom_dFdx, de_srom_dlnHdx, de_srom_rhobtorhox
!	DOUBLE PRECISION, PRIVATE 	:: Omz, Orz, Obz, Odez, Okz, Ode_iz, Ezsq

!	DOUBLE PRECISION :: de_sromlnHdata(de_num_intpl), de_sromFdata(de_num_intpl)
  CONTAINS

!########################################################################
!########################################################################
!		PUBLIC FUNCTIONS 					#
!########################################################################
!########################################################################


  !------------------------------------------
  ! Initialize the Maurice's hde model
  !------------------------------------------
	SUBROUTINE de_mauricehde_init()
		INTEGER :: i
		DOUBLE PRECISION :: z1,z2, a1,a2, ez1,ez2
			
		i	= 1
		z1 	= 0.0d0
		a1 	= 1.0d0/(1.0d0+z1)
		ez1 	= 1.0d0

		DO i = 1, de_num_intpl
			a2 = 1.0d0/(1.0d0+de_zi(i))
			ez2 = de_RK(de_mauricehde_dezda, a1, ez1, a2, 1)

			de_ezdata(i) = ez2
!			if(mod(i,100).eq.1) then
!				print *, i, x2, exp(lnH2)
!			endif
			a1 = a2; ez1 = ez2;
		ENDDO
	
	END SUBROUTINE de_mauricehde_init
	

!########################################################################
!########################################################################
!		PRIVATE FUNCTIONS 					#
!########################################################################
!########################################################################

	DOUBLE PRECISION FUNCTION de_mauricehde_dezda(a,ez)
		DOUBLE PRECISION :: a,ez

		!de_mauricehde_dezda = ez/a - 3.0d0/(a*ez) * (de_CP%Om0 * a**(-3.0d0) )
		!de_mauricehde_dezda = ez/a - 3.0d0/(a*ez) * (de_CP%Om0 * a**(-3.0d0) + de_CP%Or0 * a**(-4.0d0) * 2.0d0)
		de_mauricehde_dezda = ez/a - 3.0d0/(a*ez) * (de_CP%Om0 * a**(-3.0d0) + de_CP%Or0 * a**(-4.0d0) )
		!de_mauricehde_dezda =  - 1.5d0*de_CP%Om0 / a**(4.0d0) / ez - 2.0d0*de_CP%Or0 / a**(5.0d0) / ez 
	END FUNCTION de_mauricehde_dezda

	DOUBLE PRECISION FUNCTION de_mauricehde_q(a,ez)
		DOUBLE PRECISION :: a,ez,dezda
		dezda = de_mauricehde_dezda(a,ez)
		de_mauricehde_q = -1.0 - a/(ez) * dezda
	END FUNCTION de_mauricehde_q
		

END MODULE de_mauricehde
