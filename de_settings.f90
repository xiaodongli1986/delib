
MODULE de_settings

IMPLICIT NONE

	INTEGER, PARAMETER :: de_char_len = 100
	INTEGER, PARAMETER :: de_dl = KIND(1.0d0)
	INTEGER, PARAMETER :: de_sp = KIND(1.0)

	REAL(de_dl),PARAMETER :: de_epsion = 0.0000001d0
	
	DOUBLE PRECISION,   PARAMETER  ::  de_const_c =  2.99792458d5
	DOUBLE PRECISION,   PARAMETER  ::  de_Pi = 3.1415926535897932384626433832795d0
	DOUBLE PRECISION,   PARAMETER  ::  de_inv_twoPI = 1.0d0 / ( 2.0d0 * 3.1415926535897932384626433832795d0 ) 

!PATH 
	CHARACTER(len=de_char_len), PARAMETER :: de_data_path = '~/software/delib/data/'

! print chisq
	LOGICAL, PUBLIC	:: pr_chisq_info = .TRUE.
	
	CONTAINS
	
END MODULE de_settings


