

!-----------------------------------------------------------------------
!------------------ BEGINNING OF THE MODULE Union2 ---------------------
!-----------------------------------------------------------------------
MODULE de_qz

USE de_settings
USE de_tools
USE de_types

IMPLICIT NONE
	
	PRIVATE

	DOUBLE PRECISION, PARAMETER :: de_qz_sublab_linear=0, de_qz_sublab_2bin = 1, de_qz_sublab_deltazA=2, de_qz_sublab_poly2=3, &
		de_qz_sublab_linearzt=4, de_qz_sublab_CPL=5
	DOUBLE PRECISION :: de_qz_sublab = de_qz_sublab_linear

	PUBLIC :: de_qz_init, de_qz_sublab, de_qz_sublab_linear, de_qz_sublab_2bin, de_qz_sublab_deltazA, de_qz_sublab_poly2, &
		de_qz_sublab_linearzt, de_qz_sublab_CPL
!	DOUBLE PRECISION, PRIVATE 	:: Omz, Orz, Obz, Odez, Okz, Ode_iz, Ezsq

  CONTAINS

!########################################################################
!########################################################################
!		PUBLIC FUNCTIONS 					#
!########################################################################
!########################################################################
	

  !------------------------------------------
  ! Initialize the qz model
  !------------------------------------------
	SUBROUTINE de_qz_Init()
		INTEGER :: i
		DOUBLE PRECISION :: z1, Ez1, z2, Ez2

		i	= 1
		z1 	= 0.0d0
		Ez1 	= 1.0d0

		DO i = 1, de_num_intpl
			z2 = de_zi(i)

			! de_RK(dfundz, zleft, funleft, zright, N)
			Ez2 = de_RK(DEzDz, z1, Ez1, z2, 1)
			de_ezdata(i) = Ez2

			z1 = z2; Ez1 = Ez2
		ENDDO		
	END SUBROUTINE de_qz_Init
	

!########################################################################
!########################################################################
!		PRIVATE FUNCTIONS 					#
!########################################################################
!########################################################################

  !----------------------------------------------------------------
  ! qzs
  !----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION qz_linear(z)
		DOUBLE PRECISION :: z
		qz_linear = de_CP%qz%q0 + de_CP%qz%q1 * z
	END FUNCTION qz_linear
	DOUBLE PRECISION FUNCTION qz_poly2(z)
		DOUBLE PRECISION :: z
		qz_poly2 = de_CP%qz%q0 + de_CP%qz%q1 * z + de_CP%qz%q2 * z**2.0
	END FUNCTION qz_poly2
	DOUBLE PRECISION FUNCTION qz_linearzt(z)
		DOUBLE PRECISION :: z
		qz_linearzt = de_CP%qz%q0 + de_CP%qz%q1 * (z - de_CP%qz%zt)
	END FUNCTION qz_linearzt
	DOUBLE PRECISION FUNCTION qz_CPL(z)
		DOUBLE PRECISION :: z
		qz_CPL = de_CP%qz%q0 + de_CP%qz%q1 * z / (1.0+z)
	END FUNCTION qz_CPL
	DOUBLE PRECISION FUNCTION qz_2bin(z)
		DOUBLE PRECISION :: z, q0, q1, zt
		zt = de_CP%qz%zt
		if(z<=zt) then
			qz_2bin = de_CP%qz%q0 
		else
			qz_2bin = de_CP%qz%q1
		endif
	END FUNCTION qz_2bin
	DOUBLE PRECISION FUNCTION qz_deltazA(z)
		DOUBLE PRECISION :: z, q0, q1, zt, deltaz
		zt = de_CP%qz%zt
		deltaz = de_CP%qz%deltaz
		if(z<=zt-deltaz) then
			qz_deltazA=de_CP%qz%q0
		elseif(z<=zt+deltaz .and. z>zt-deltaz) then
			qz_deltazA=(de_CP%qz%q1-de_CP%qz%q0) / (2*de_CP%qz%deltaz) * (z - de_CP%qz%zt - de_CP%qz%deltaz) + de_CP%qz%q1
		elseif(z<zt-deltaz) then
			qz_deltazA=de_CP%qz%q1
		endif
	END FUNCTION qz_deltazA

  !----------------------------------------------------------------
  ! d E(z) / dz
  !----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION DEzDz(z, Ez)
		DOUBLE PRECISION :: z, Ez, qz
		if(de_qz_sublab.eq.de_qz_sublab_linear) then
			qz = qz_linear(z)
		elseif(de_qz_sublab.eq.de_qz_sublab_deltazA) then
			qz = qz_deltazA(z)
		elseif(de_qz_sublab.eq.de_qz_sublab_2bin) then
			qz = qz_2bin(z)
		elseif(de_qz_sublab.eq.de_qz_sublab_poly2) then
			qz = qz_poly2(z)
		elseif(de_qz_sublab.eq.de_qz_sublab_linearzt) then
			qz = qz_linearzt(z)
		elseif(de_qz_sublab .eq. de_qz_sublab_CPL) then
			qz = qz_CPL(z)
		endif
		DEzDz = Ez * (1.0+qz) / (1.0+z)
	END FUNCTION DEzDz

END MODULE de_qz
