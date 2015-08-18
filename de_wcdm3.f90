

!-----------------------------------------------------------------------
!------------------ BEGINNING OF THE MODULE Union2 ---------------------
!-----------------------------------------------------------------------
MODULE de_wcdm3

USE de_settings
USE de_tools
USE de_types

IMPLICIT NONE
	
	PRIVATE

	PUBLIC :: de_wcdm3_init, de_wcdm3_inv_e, de_wcdm3_rhode, de_wcdm3_w
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
  	DOUBLE PRECISION FUNCTION de_wcdm3_inv_e(z)    
		DOUBLE PRECISION :: z
		
		de_wcdm3_inv_e = 1.0d0 / sqrt( &
			de_CP%Om0*(1.0+z)**3.0 &
			+ de_CP%Or0*(1.0+z)**4.0 &
			+ de_CP%Ok0*(1.0+z)**2.0 &
			+ de_wcdm3_rhode(z) * de_CP%Ode0 )
			
	END FUNCTION de_wcdm3_inv_e


  !------------------------------------------
  ! rho_wcdm3(z),
  !  propotional to wcdm3 energy density, normalized to 1 at z = 0
  !------------------------------------------	
	DOUBLE PRECISION FUNCTION de_wcdm3_rhode(z)
		INTEGER :: i, n
		DOUBLE PRECISION :: z, zi(3), wi(3)
		
		zi(1) = de_CP%wcdm3%z1
		zi(2) = de_CP%wcdm3%z2
		wi(1) = de_CP%wcdm3%w1
		wi(2) = de_CP%wcdm3%w2
		wi(3) = de_CP%wcdm3%w3
		
		IF(z .le. zi(1)) THEN
			n = 1
			ELSEIF(z .le. zi(2)) THEN
			n = 2
			ELSE
			n = 3
		ENDIF
		
		de_wcdm3_rhode = (1.0+z) ** (3.0 * (1.0+wi(n)))
		
		DO i = 1, n-1
			de_wcdm3_rhode = de_wcdm3_rhode * (1+zi(i))**(3.0*(wi(i)-wi(i+1)))
		ENDDO

	END FUNCTION de_wcdm3_rhode

  !------------------------------------------
  ! w_wcdm3(z)
  !------------------------------------------	
	DOUBLE PRECISION FUNCTION de_wcdm3_w(z)
		DOUBLE PRECISION :: z

		IF(z .le. de_CP%wcdm3%z1) THEN
			de_wcdm3_w = de_CP%wcdm3%w1
			return
		ENDIF
		
		IF(z .le. de_CP%wcdm3%z2) THEN
			de_wcdm3_w = de_CP%wcdm3%w2
			return
		ENDIF
		
		de_wcdm3_w = de_CP%wcdm3%w3
		
	END FUNCTION de_wcdm3_w

  !------------------------------------------
  ! Initialize the wcdm3 model
  !------------------------------------------
	SUBROUTINE de_wcdm3_Init()
	
	END SUBROUTINE de_wcdm3_Init
	

!########################################################################
!########################################################################
!		PRIVATE FUNCTIONS 					#
!########################################################################
!########################################################################
END MODULE de_wcdm3
