

!-----------------------------------------------------------------------
!------------------ BEGINNING OF THE MODULE Union2 ---------------------
!-----------------------------------------------------------------------
MODULE de_srom

USE de_settings
USE de_tools
USE de_types

IMPLICIT NONE
	
	PRIVATE

	PUBLIC :: de_srom_init, de_srom_dFdx, de_srom_dlnHdx, de_srom_rhobtorhox
!	DOUBLE PRECISION, PRIVATE 	:: Omz, Orz, Obz, Odez, Okz, Ode_iz, Ezsq

!	DOUBLE PRECISION :: de_sromlnHdata(de_num_intpl), de_sromFdata(de_num_intpl)
  CONTAINS

!########################################################################
!########################################################################
!		PUBLIC FUNCTIONS 					#
!########################################################################
!########################################################################

	!!! TBC
	DOUBLE PRECISION FUNCTION de_srom_rhobtorhox(x,F)
		DOUBLE PRECISION :: x, F, a, r0
		a = exp(x) ! scale fator as function of x=ln a
		r0 = de_CP%Odm0/(1.0-de_CP%Odm0-de_CP%Ob0-de_CP%Or0)
		de_srom_rhobtorhox = de_CP%Ob0/de_CP%Odm0*(a**-3.0)*r0*exp(F)
	END FUNCTION de_srom_rhobtorhox
	
	!!! TBC
	DOUBLE PRECISION FUNCTION de_srom_dlnHdx(x,lnH,F)
		DOUBLE PRECISION :: x, lnH, F, a, r0, y
		a = exp(x) ! scale fator as function of x=ln a
		r0 = de_CP%Odm0/(1.0-de_CP%Odm0-de_CP%Ob0-de_CP%Or0)
		de_srom_dlnHdx = (de_CP%srom%w0+de_CP%srom%w1*(1.0-a)) / (1.0+r0*a**(-de_CP%srom%xi)+ (1.0d0+de_CP%Or0/de_CP%Ob0/a) * de_srom_rhobtorhox(x,F))
		
		y =  (1.0+r0*a**(-de_CP%srom%xi)) * (de_CP%Odm0/de_CP%Or0) * (a**4.0d0) / r0 * Exp(-F)
		y = (1.0/3.0d0) / (1.0+de_CP%Ob0/de_CP%Or0*a+y)
		
		de_srom_dlnHdx = -1.5d0 * (1.0 + de_srom_dlnHdx + y)
	END FUNCTION de_srom_dlnHdx

	DOUBLE PRECISION FUNCTION de_srom_dFdx(x,lnH,F)
		DOUBLE PRECISION :: x, lnH, F, a, r0, r0axi
		a = exp(x) ! scale fator as function of x=ln a
		r0 = de_CP%Odm0/(1.0-de_CP%Odm0-de_CP%Ob0-de_CP%Or0)
		r0axi = r0*(a**(-de_CP%srom%xi))
		de_srom_dFdx = (de_CP%srom%xi/3.0d0 + de_CP%srom%w0 + de_CP%srom%w1*(1.0d0-a)) * (r0axi/(1.0d0+r0axi))
		de_srom_dFdx = 3.0d0 * (1.0d0 + de_CP%srom%w0 + de_CP%srom%w1*(1.0d0-a) - de_srom_dFdx) 
	END FUNCTIOn de_srom_dFdx

  !------------------------------------------
  ! Initialize the srom model
  !------------------------------------------
	SUBROUTINE de_srom_init()
		INTEGER :: i
		DOUBLE PRECISION :: z1,z2, x1,x2, lnH1,F1,lnH2,F2
			
		i	= 1
		z1 	= 0.0d0
		x1 	= log(1.0d0/(1.0d0+z1))
		lnH1 	= log(de_CP%H0)
		F1 	= 0.0d0

		DO i = 1, de_num_intpl
			z2 = de_zi(i)
			x2 = log(1.0d0/(1.0d0+z2))
			CALL de_RK2(de_srom_dlnHdx, de_srom_dFdx, x1, lnH1, F1, x2, lnH2, F2, 1)

			de_ezdata(i) = exp(lnH2) / de_CP%H0
!			if(mod(i,100).eq.1) then
!				print *, i, x2, exp(lnH2)
!			endif
			z1 = z2; x1=x2; lnH1 = lnH2; F1 = F2
		ENDDO
	
	END SUBROUTINE de_srom_Init
	

!########################################################################
!########################################################################
!		PRIVATE FUNCTIONS 					#
!########################################################################
!########################################################################
END MODULE de_srom
