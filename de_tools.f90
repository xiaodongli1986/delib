!-----------------------------------------------------------------------
!------------------ Useful Tools. Including Simpson Integral, ----------
!------------------      RK, interpolation, read_in/output, and so on --
!-----------------------------------------------------------------------
MODULE de_tools

USE de_settings

IMPLICIT NONE

	!number of points used for interpolations in redshift.
	INTEGER,  PARAMETER		:: de_num_intpl = 10650
	
	!If we take the bound in integral of rs to be 1.0d-9,
	! then 2700  for basenumber 128;
	!      5350  for basenumber 256; 
	!      10650 for basenumber 512.

	!If we take cut (5.0d5),
	! then	850   for basenumber 64  (typical speed: 0.47s / 100p) 
	!	1700  for basenumber 128 (typical speed: 0.60s / 100p)
	!	3400  for basenumber 256 (typical speed: 0.83s / 100p) 
	!	6750  for basenumber 512 (typical speed: 1.26s / 100p)

	DOUBLE PRECISION,  PARAMETER	:: de_basenumber = 1.0d0 + 1.0d0/512.0d0
	!DOUBLE PRECISION,  PARAMETER	:: de_basenumber = 1.0d0 + 1.0d0/128.0d0
	DOUBLE PRECISION,  PARAMETER	:: de_logbasenumber = log(de_basenumber)
	DOUBLE PRECISION :: de_maxintplz, de_minintpla
	

	
	!Common used array, saving the interpolating data.
	! Frequently used in many dark energy models, so put them here as common variables.
	! Array for redshift.
	DOUBLE PRECISION :: de_zdata(de_num_intpl)
	! Array for ez, defined as H(z)/H(0)
	DOUBLE PRECISION :: de_ezdata(de_num_intpl)
	! Array for \Omega_de(z), defined as rho_de(z)/rho_crit(z)
	DOUBLE PRECISION :: de_odezdata(de_num_intpl)
	! Array for dark energy equation of state
	DOUBLE PRECISION :: de_wzdata(de_num_intpl)
	! Array for rhode (z) / rho(z=0)
	DOUBLE PRECISION :: de_rhodezdata(de_num_intpl)
	
	LOGICAL :: de_tools_inited = .FALSE.
	
	
CONTAINS

  !------------------------------------------
  ! Initialization for de_tools
  !------------------------------------------
	SUBROUTINE de_tools_init()
		INTEGER :: i
		
		WRITE(*,*) "Initializing de_tools..."
		WRITE(*,*) " Redshift interpolating numbers = ", de_num_intpl

		!Initialize the redshift data
		DO i = 1, de_num_intpl
			de_zdata(i) = de_zi(i)
		ENDDO
		
		!Range of the interpolation (maximal redshift)
		de_maxintplz = de_zdata(de_num_intpl)
		de_minintpla = 1.0/(1.0+de_maxintplz)
		WRITE(*,*) " Maximal redshift in interpolating = ", de_maxintplz
				
		de_tools_inited = .TRUE.
	END SUBROUTINE de_tools_init


  !------------------------------------------
  ! Get z from index
  !------------------------------------------
	DOUBLE PRECISION FUNCTION de_zi(i)
		INTEGER :: i
		de_zi = de_basenumber**DBLE(i-1) - 1.0d0
	END FUNCTION de_zi

  !------------------------------------------
  ! Get index from z
  !------------------------------------------
	INTEGER FUNCTION de_iz(z)
		DOUBLE PRECISION :: z, temp
		de_iz = Ceiling(log(1.0d0+z)/de_logbasenumber + 1.0d0)
		IF(de_iz < 2) de_iz = 2
	END FUNCTION de_iz

  !------------------------------------------
  ! 2th order interpolating function.
  !------------------------------------------
	DOUBLE PRECISION FUNCTION de_intpl_vl(x, x1, f1, x2, f2, x3, f3)
		DOUBLE PRECISION :: x, x1, f1, x2, f2, x3, f3
		DOUBLE PRECISION :: d1, d2, d3, d12, d13, d23
		d1  = x - x1
		d2  = x - x2
		d3  = x - x3
		d12 = x1 - x2
		d13 = x1 - x3
		d23 = x2 - x3
		de_intpl_vl = f1*d2*d3/(d12*d13) - f2*d1*d3/(d12*d23) + f3*d1*d2/(d13*d23)
	END FUNCTION de_intpl_vl


  !------------------------------------------
  ! get the value of e(z) from the ezdata
  !------------------------------------------
	DOUBLE PRECISION FUNCTION de_get_ez(z)    
		INTEGER :: i, i1, i2, i3
		DOUBLE PRECISION :: z, z1, z2, z3, f1, f2, f3

		i = de_iz(z)

		IF(i < 2) i=2
		IF(i > de_num_intpl - 1) i = de_num_intpl -1

		i1 = i-1 
		i2 = i 
		i3 = i+1
		
		z1 = de_zdata(i1)
		z2 = de_zdata(i2)
		z3 = de_zdata(i3)

		f1 = de_ezdata(i1); 
		f2 = de_ezdata(i2); 
		f3 = de_ezdata(i3)

		de_get_ez = de_intpl_vl(z, z1, f1, z2, f2, z3, f3)
	END FUNCTION de_get_ez


  !------------------------------------------
  ! get the value of omega_de(z) from odezdata
  !------------------------------------------
	DOUBLE PRECISION FUNCTION de_get_odez(z)    
		INTEGER :: i, i1, i2, i3
		DOUBLE PRECISION :: z, z1, z2, z3, f1, f2, f3

		i = de_iz(z)

		IF(i < 2) i=2
		IF(i > de_num_intpl - 1) i = de_num_intpl -1

		i1 = i-1 
		i2 = i 
		i3 = i+1
		
		z1 = de_zdata(i1)
		z2 = de_zdata(i2)
		z3 = de_zdata(i3)
		
		f1 = de_odezdata(i1) 
		f2 = de_odezdata(i2)
		f3 = de_odezdata(i3)

		de_get_odez = de_intpl_vl(z, z1, f1, z2, f2, z3, f3)
	END FUNCTION de_get_odez


  !------------------------------------------
  ! get the value of w(z) from wzdata
  !------------------------------------------
	DOUBLE PRECISION FUNCTION de_get_wz(z)    
		INTEGER :: i, i1, i2, i3
		DOUBLE PRECISION :: z, z1, z2, z3, f1, f2, f3

		i = de_iz(z)

		IF(i < 2) i=2
		IF(i > de_num_intpl - 1) i = de_num_intpl -1

		i1 = i-1 
		i2 = i 
		i3 = i+1
		
		z1 = de_zdata(i1)
		z2 = de_zdata(i2)
		z3 = de_zdata(i3)
		
		f1 = de_wzdata(i1) 
		f2 = de_wzdata(i2)
		f3 = de_wzdata(i3)

		de_get_wz = de_intpl_vl(z, z1, f1, z2, f2, z3, f3)
	END FUNCTION de_get_wz
	
  !------------------------------------------
  ! get the value of w(z) from wzdata
  !------------------------------------------
	DOUBLE PRECISION FUNCTION de_get_rhodez(z)    
		INTEGER :: i, i1, i2, i3
		DOUBLE PRECISION :: z, z1, z2, z3, f1, f2, f3

		i = de_iz(z)

		IF(i < 2) i=2
		IF(i > de_num_intpl - 1) i = de_num_intpl -1

		i1 = i-1 
		i2 = i 
		i3 = i+1
		
		z1 = de_zdata(i1)
		z2 = de_zdata(i2)
		z3 = de_zdata(i3)
		
		f1 = de_rhodezdata(i1) 
		f2 = de_rhodezdata(i2)
		f3 = de_rhodezdata(i3)

		de_get_rhodez = de_intpl_vl(z, z1, f1, z2, f2, z3, f3)
	END FUNCTION de_get_rhodez


  !---------------------------------------------------------------
  ! Simpson Integration
  !---------------------------------------------------------------
	DOUBLE PRECISION FUNCTION de_Simpson(fun,xleft,xright,N)   
		DOUBLE PRECISION, external :: fun
		DOUBLE PRECISION :: x1,x2,xleft,xright,BC
		DOUBLE PRECISION :: f1,f2
		INTEGER :: i, N
		BC=(xright-xleft)/DBLE(N)
		x1=xleft;x2=x1+BC;
		f1=fun(x1);f2=fun(x2);de_Simpson=(f1+fun((x1+x2)*0.5d0)*4.0d0+f2)*BC/6.0d0;
	
		DO i = 2,N
			x1=x2;f1=f2;x2=x2+BC;
			f2=fun(x2);de_Simpson=de_Simpson+(f1+fun((x1+x2)*0.5d0)*4.0d0+f2)*BC/6.0d0; 
		ENDDO
	END FUNCTION de_Simpson

  !---------------------------------------------------------------
  ! 4rd order Runge-Kutta 
  !---------------------------------------------------------------
	DOUBLE PRECISION FUNCTION de_RK(dfundz, zleft, funleft, zright, N)
		DOUBLE PRECISION, EXTERNAL :: dfundz
		DOUBLE PRECISION :: zleft, funleft, zright, nowz, nowfun
		DOUBLE PRECISION :: BC, K1, K2, K3, K4
		INTEGER :: i, N
		BC = (zright-zleft) / N
		nowz = zleft
		nowfun = funleft
		DO i = 1, N
			K1 = dfundz(nowz, nowfun)
			K2 = dfundz(nowz + 0.5d0*BC, nowfun + 0.5d0*BC*K1)
			K3 = dfundz(nowz + 0.5d0*BC, nowfun + 0.5d0*BC*K2)
			K4 = dfundz(nowz + BC, nowfun + BC*K3)
			nowfun = nowfun + BC * (K1+2.0d0*K2+2.0d0*K3+K4) / 6.0d0
			nowz = nowz + BC
		ENDDO
		de_RK = nowfun
	END FUNCTION de_RK


  !---------------------------------------------------------------
  ! 4rd order Runge-Kutta for two functions.
  !---------------------------------------------------------------
	SUBROUTINE de_RK2(dfunAdz, dfunBdz, zleft, funAleft, funBleft, zright, funAright, funBright, N)
		DOUBLE PRECISION, EXTERNAL :: dfunAdz, dfunBdz
		DOUBLE PRECISION, INTENT(IN) :: zleft, funAleft, funBleft, zright
		DOUBLE PRECISION, INTENT(OUT):: funAright, funBright
		DOUBLE PRECISION :: nowz, nowfunA, nowfunB
		DOUBLE PRECISION :: BC, K1A, K1B, K2A, K2B, K3A, K3B, K4A, K4B
		INTEGER :: i, N
		BC = (zright-zleft) / N
		nowz = zleft
		nowfunA = funAleft
		nowfunB = funBleft
		DO i = 1, N
			K1A = dfunAdz(nowz, nowfunA, nowfunB)
			K1B = dfunBdz(nowz, nowfunA, nowfunB)
			K2A = dfunAdz(nowz + 0.5d0*BC, nowfunA + 0.5d0*BC*K1A, nowfunB + 0.5d0*BC*K1B)
			K2B = dfunBdz(nowz + 0.5d0*BC, nowfunA + 0.5d0*BC*K1A, nowfunB + 0.5d0*BC*K1B)
			K3A = dfunAdz(nowz + 0.5d0*BC, nowfunA + 0.5d0*BC*K2A, nowfunB + 0.5d0*BC*K2B)
			K3B = dfunBdz(nowz + 0.5d0*BC, nowfunA + 0.5d0*BC*K2A, nowfunB + 0.5d0*BC*K2B)
			K4A = dfunAdz(nowz + BC, nowfunA + BC*K3A, nowfunB + 0.5d0*BC*K3B)
			K4B = dfunBdz(nowz + BC, nowfunA + BC*K3A, nowfunB + 0.5d0*BC*K3B)
			nowfunA = nowfunA + BC * (K1A+2.0d0*K2A+2.0d0*K3A+K4A) / 6.0d0
			nowfunB = nowfunB + BC * (K1B+2.0d0*K2B+2.0d0*K3B+K4B) / 6.0d0
			nowz = nowz + BC
		ENDDO
		funAright = nowfunA
		funBright = nowfunB
	END SUBROUTINE de_RK2
	
	
  !-----------------------------------------------------------
  ! Get the number of lines of a file.
  !-----------------------------------------------------------
	SUBROUTINE de_count_line_number (file_name, line_number)
		INTEGER :: line_number
		CHARACTER(len=*) :: file_name
		CHARACTER(len=de_char_len) :: inline
    
		OPEN(UNIT=1456,FILE=file_name,ERR=2)
    
		line_number = 0
		DO While(1 .EQ. 1)
			READ(1456,*,ERR=3,END=100) inline
			line_number = line_number + 1
		ENDDO

2		WRITE(*,*) "Error occurs when opening the file ", file_name
		CLOSE(1456)
		STOP

3		WRITE(*,*) "Error occurs when couting the lines of the file ", file_name
		CLOSE(1456)
		STOP

100		CLOSE(1456)
	END SUBROUTINE de_count_line_number 


  !-----------------------------------------------------------
  ! Output a two dimensional table to a file.
  !-----------------------------------------------------------
	SUBROUTINE de_output_2d (file_name, output_data)
		CHARACTER(LEN=de_char_len) :: file_name
		DOUBLE PRECISION :: output_data(:,:)
		INTEGER :: d1, d2 
		INTEGER  :: i,j
		CHARACTER(LEN=de_char_len) :: tmpstr1, tmpstr2

		d1 = SIZE(output_data, 1)
		d2 = SIZE(output_data, 2)

		OPEN(UNIT=9873,FILE=file_name,ERR=22)

!		DO I=1, d1
!			WRITE(9873,'(<d2>(e14.7,2x))',ERR=33) output_data(i,1:d2)
!		ENDDO
		DO I=1, d1
			tmpstr1 = ''
			DO j=1, d2
				WRITE(tmpstr2,'(e14.7)',ERR=33) output_data(i,j)
				tmpstr1 = trim(adjustl(tmpstr1))//' '//trim(adjustl(tmpstr2))
			ENDDO
			WRITE(9873,'(A)',ERR=33) trim(adjustl(tmpstr2))
		ENDDO

		CLOSE(9873)
		RETURN

22		WRITE(*,*) "Error occurs when opening the file ", file_name
		CLOSE(9873)
		STOP

33		WRITE(*,*) "Error occurs when writing into the file ", file_name
		CLOSE(9873)
		STOP
	END SUBROUTINE de_output_2d

  !-----------------------------------------------------------
  ! Read in a file with given name and number of columns.
  !-----------------------------------------------------------
	SUBROUTINE de_read_in (file_name, n_column, n_lines, total_data)
		CHARACTER(LEN=de_char_len), INTENT(IN) :: file_name
		INTEGER, INTENT(IN)  :: n_column
		INTEGER, INTENT(OUT) :: n_lines
		DOUBLE PRECISION, ALLOCATABLE :: total_data(:,:)
		INTEGER :: i
		CALL de_count_line_number(file_name, n_lines)
		ALLOCATE(total_data(n_lines, n_column))
		OPEN(UNIT=4321,FILE=file_name)
		DO i = 1, n_lines
			READ(4321, *) total_data(i,1:n_column)
		ENDDO
		CLOSE(4321)
	END SUBROUTINE de_read_in

  !---------------------------------------------------------------
  ! inversion of a matrix
  !---------------------------------------------------------------      
        subroutine de_nizhen(aa,b,n)
                ! Dummy
                double precision, intent(in) :: aa(n,n)
                integer, intent(in) :: n
                double precision, intent(out) :: b(n,n)
                ! Local
                integer :: i,j,k
                double precision :: a(n,n)
                a=aa
                b=0.0d0
                do i=1,n
                        b(i,i)=1
                enddo
                do i=1,n
                        b(i,:)=b(i,:)/a(i,i)
                        a(i,i:n)=a(i,i:n)/a(i,i)
                        do j=i+1,n
                                do k=1,n
                                        b(j,k)=b(j,k)-b(i,k)*a(j,i)
                                enddo
                                a(j,i:n)=a(j,i:n)-a(i,i:n)*a(j,i)
                        enddo
                enddo
                do i=n,1,-1
                do j=i-1,1,-1
                do k=1,n
                        b(j,k)=b(j,k)-b(i,k)*a(j,i)
                enddo
                enddo
                enddo
        end subroutine de_nizhen

END MODULE de_tools
