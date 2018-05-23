
MODULE de_types

USE de_settings
USE de_tools

IMPLICIT NONE

	!Free parameters for different models; hde, wcdm3, wcdm, cpl, ...
	TYPE :: hde_para
		DOUBLE PRECISION :: c
		!Interaction parameters. See arXiv:1204.6135 for their meanings.
		DOUBLE PRECISION :: Gamma=0
		DOUBLE PRECISION :: a, b  
	END TYPE

	TYPE :: wcdm3_para
		DOUBLE PRECISION :: z1 = 0.5d0, z2 = 1.0d0
		DOUBLE PRECISION :: w1, w2, w3
	END TYPE
	
	INTEGER, PARAMETER :: de_zbinned_type=1, de_abinned_type=2
	TYPE :: w_binned_para
		DOUBLE PRECISION :: zmax = 1.5d0
		INTEGER :: nbins = 30 ! should be nbins < 1000
		DOUBLE PRECISION :: ws(1000)
		DOUBLE PRECISION :: whighz = -1.0d0
		INTEGER :: binnedmethod = de_zbinned_type
	END TYPE	
	
	
	TYPE :: wcdm_para
		DOUBLE PRECISION :: w
	END TYPE
	
	TYPE :: CPL_para
		DOUBLE PRECISION :: w0, wa
	END TYPE
	
	TYPE :: srom_para
		DOUBLE PRECISION :: xi, w0, w1
	END TYPE
	
	TYPE :: ICG_para
		DOUBLE PRECISION :: xi, A
	END TYPE

	TYPE :: qz_para
		DOUBLE PRECISION :: q0, q1, q2, q3, zt, deltaz
	END TYPE

	TYPE :: Rhct_para
		DOUBLE PRECISION :: alpha
	END TYPE
	
	TYPE :: coupled_de_para
		LOGICAL :: use_xi1 
		DOUBLE PRECISION :: wde, xi1,xi2
	END TYPE


	TYPE :: de_para

		!common parameter parameters for most models
		DOUBLE PRECISION :: Ob0hsq, Odm0, h
		DOUBLE PRECISION :: Ok0 = 0.0 ! Default flat
		DOUBLE PRECISION :: alpha, beta !nuisance paramters for snls

		!Parameters for specific de models.
		!Once you add a new de model, add its parameters in this structure.
		
		!constant w
		TYPE(wcdm_para)  :: wcdm
		!CPL
		TYPE(CPL_para)   :: CPL
		!holographic dark energy
		TYPE(hde_para) 	 :: hde
		!Binned w(z) (3 bins)
		TYPE(wcdm3_para) :: wcdm3
		! srom model
		TYPE(srom_para) :: srom
		! ICG model
		TYPE(ICG_para) :: ICG
		! qz model
		TYPE(qz_para) :: qz
		! Rhct model
		TYPE(Rhct_para) :: Rhct
		
		TYPE(w_binned_para) :: w_binned
		TYPE(coupled_de_para) :: coupled_de


		!derived parameters
		DOUBLE PRECISION :: Ob0
		DOUBLE PRECISION :: Odm0hsq
		DOUBLE PRECISION :: Om0hsq, Om0
		DOUBLE PRECISION :: Ode0, Or0, Og0 !dark energy, radiation, and photon
		DOUBLE PRECISION :: H0
		
		DOUBLE PRECISION :: zd, rszd
		DOUBLE PRECISION :: zstar, DAzstar, rszstar
		DOUBLE PRECISION :: R, lA
		DOUBLE PRECISION :: ns
	END TYPE de_para
	
	!Global colleection of parameters
	TYPE(de_para) :: de_CP 

CONTAINS



END MODULE de_types	
