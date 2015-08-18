
This will build a library for dark energy.

-----------------
Installation
-----------------

To install it, firstly modify the path for data 

	CHARACTER(len=de_char_len), PARAMETER :: de_data_path = ...
		
defined in the file 'de_chisqs.f90'. Then, compile the codes using

	make
	
After make, please modify the first line of the file 'de.sh'

	DElibPATH=...
	
and add the following command to ~/.bashrc:

	source YOUR_SHFILE_PATH/de.sh
	
Then you can use a command like

	ifort Test.f90 $delm -mkl
	
or equivalently

	ifort Test.f90 -lde -I$demods -mkl
	
to compile f90 file with functions defined in delib.

#######################
As an example, test/Test.f90 provides a test:
	ifort Test.f90 $delm -mkl
	./a.out


#################################################################
List of modules
#################################################################

de_chisqs
de_model_init
de_settings
de_tools
de_types

de_hde
de_wcm3

#################################################################
List of Tools functions
#################################################################

  !------------------------------------------
  ! Get z from index
  !------------------------------------------
	DOUBLE PRECISION FUNCTION de_zi(i)
		INTEGER :: i

  !------------------------------------------
  ! Get index from z
  !------------------------------------------
	INTEGER FUNCTION de_iz(z)
		DOUBLE PRECISION :: z

  !------------------------------------------
  ! 2th order interpolating function.
  !------------------------------------------
	DOUBLE PRECISION FUNCTION de_intpl_vl(x, x1, f1, x2, f2, x3, f3)
		DOUBLE PRECISION :: x, x1, f1, x2, f2, x3, f3
		
  !---------------------------------------------------------------
  ! Simpson Integration
  !---------------------------------------------------------------
	DOUBLE PRECISION FUNCTION de_Simpson(fun,xleft,xright,N)   
		DOUBLE PRECISION, external :: fun
		DOUBLE PRECISION :: xleft,xright
		INTEGER :: N

  !---------------------------------------------------------------
  ! 4rd order Runge-Kutta 
  !---------------------------------------------------------------
	DOUBLE PRECISION FUNCTION de_RK(dfundz, zleft, funleft, zright, N)
		DOUBLE PRECISION, EXTERNAL :: dfundz
		DOUBLE PRECISION :: zleft, funleft, zright, nowz, nowfun

  !---------------------------------------------------------------
  ! 4rd order Runge-Kutta for two functions.
  !---------------------------------------------------------------
	SUBROUTINE de_RK2(dfunAdz, dfunBdz, zleft, funAleft, funBleft, zright, funAright, funBright, N)
		DOUBLE PRECISION, EXTERNAL :: dfunAdz, dfunBdz
		DOUBLE PRECISION, INTENT(IN) :: zleft, funAleft, funBleft, zright
		DOUBLE PRECISION, INTENT(OUT):: funAright, funBright
		
  !-----------------------------------------------------------
  ! Get the number of lines of a file.
  !-----------------------------------------------------------
	SUBROUTINE de_count_line_number (file_name, line_number)
		INTEGER :: line_number
		CHARACTER(LEN=de_char_len) :: file_name
		
  !-----------------------------------------------------------
  ! Output a two dimensional table to a file.
  !-----------------------------------------------------------
	SUBROUTINE de_output_2d (file_name, output_data)
		CHARACTER(LEN=de_char_len) :: file_name
		DOUBLE PRECISION :: output_data(:,:)
			
  !-----------------------------------------------------------
  ! Read in a file with given name and number of columns.
  !-----------------------------------------------------------
	SUBROUTINE de_read_in (file_name, n_column, n_lines, total_data)
		CHARACTER(LEN=de_char_len), INTENT(IN) :: file_name
		INTEGER, INTENT(IN)  :: n_column
		INTEGER, INTENT(OUT) :: n_lines
		DOUBLE PRECISION, ALLOCATABLE :: total_data(:,:)

	
		
#################################################################
List of chisq functions
#################################################################		

	!supernovae chisq functions
	PUBLIC :: de_chisq_snls3, de_chisq_union2p1, de_chisq_bao_ver3
	
	!bao chisq functions
	PUBLIC :: de_chisq_6dFGS!, de_chisq_2dFGS
	PUBLIC :: de_chisq_SDSSDR7_old, de_chisq_SDSSDR7_new, de_chisq_BOSSDR9, de_chisq_dr11
	PUBLIC :: de_chisq_WigZ_rsDV, de_chisq_WigZ_A, de_chisq_impwig
	PUBLIC :: de_chisq_bao_ver1, de_chisq_bao_ver2
	
	!CMB chisq functions
	PUBLIC :: de_chisq_wmap7, de_chisq_wmap9, de_chisq_planck
	
	!H0 chisq functions
	PUBLIC :: de_chisq_h_Riess, de_chisq_h_Carnegie
