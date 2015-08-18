

!-----------------------------------------------------------------------
!------------------ BEGINNING OF THE MODULE Union2 ---------------------
!-----------------------------------------------------------------------
MODULE de_ICG

USE de_settings
USE de_tools
USE de_types

IMPLICIT NONE
	
	PRIVATE

	PUBLIC :: de_ICG_init, de_ICG_dudx, de_ICG_dvdx
!	DOUBLE PRECISION, PRIVATE 	:: Omz, Orz, Obz, Odez, Okz, Ode_iz, Ezsq

!	DOUBLE PRECISION :: de_sromlnHdata(de_num_intpl), de_sromFdata(de_num_intpl)
  CONTAINS

!########################################################################
!########################################################################
!		PUBLIC FUNCTIONS 					#
!########################################################################
!########################################################################
	
	DOUBLE PRECISION FUNCTION de_ICG_dudx(x,u,v)
		DOUBLE PRECISION :: x, u, v, w, a, r0
		w = -v**(-2.0) * de_CP%ICG%A
		a = exp(x) ! scale fator as function of x=ln a
		r0 = de_CP%Odm0/(1.0-de_CP%Odm0-de_CP%Ob0-de_CP%Or0)
		de_ICG_dudx = (de_CP%ICG%xi/3.0 + w) / (1.0d0 + r0*a**(-de_CP%ICG%xi))
		de_ICG_dudx = -3.0 * u * (1.0 + de_ICG_dudx)
	END FUNCTION de_ICG_dudx

	DOUBLE PRECISION FUNCTION de_ICG_dvdx(x,u,v)
		DOUBLE PRECISION :: x, u, v, w, a, r0
		w = -v**(-2.0) * de_CP%ICG%A
		a = exp(x) ! scale fator as function of x=ln a
		r0 = de_CP%Odm0/(1.0-de_CP%Odm0-de_CP%Ob0-de_CP%Or0)
		de_ICG_dvdx = (de_CP%ICG%xi/3.0 + w) / (1.0d0 + r0*a**(-de_CP%ICG%xi))
		de_ICG_dvdx = -3.0*v*(1.0+w) + 3.0*u*de_ICG_dvdx 
	END FUNCTION de_ICG_dvdx

  !------------------------------------------
  ! Initialize the srom model
  !------------------------------------------
	SUBROUTINE de_ICG_init()
		INTEGER :: i
		DOUBLE PRECISION :: z1,z2, x1,x2, u1,v1,u2,v2
			
		i	= 1
		z1 	= 0.0d0
		x1 	= log(1.0d0/(1.0d0+z1))
		u1 	= de_CP%Odm0
		v1 	= 1.0-de_CP%Odm0-de_CP%Ob0-de_CP%Or0

		DO i = 1, de_num_intpl
			z2 = de_zi(i)
			x2 = log(1.0d0/(1.0d0+z2))
			CALL de_RK2(de_ICG_dudx, de_ICG_dvdx, x1,u1,v1, x2,u2,v2,1)

			de_ezdata(i) = (u2+v2+de_CP%Ob0*(1.0+z2)**3.0 + de_CP%Or0*(1.0+z2)**4.0)**0.5
!			if(mod(i,100).eq.1) then
!				print *, i, x2, exp(lnH2)
!			endif
			z1 = z2; x1=x2; u1 = u2; v1 = v2
		ENDDO
	
	END SUBROUTINE de_ICG_Init
	

!########################################################################
!########################################################################
!		PRIVATE FUNCTIONS 					#
!########################################################################
!########################################################################
END MODULE de_ICG
