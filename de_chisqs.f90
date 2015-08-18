

!-----------------------------------------------------------------------
!------------------ BEGINNING OF THE MODULE Union2 ---------------------
!-----------------------------------------------------------------------
MODULE de_chisqs

USE de_model_init
USE de_chisqs_JLA

IMPLICIT NONE

	PRIVATE
	PUBLIC :: de_chisq_all
	
	!supernovae chisq functions
	PUBLIC :: de_chisq_snls3, de_chisq_union2p1, de_chisq_union2p1zc, de_chisq_union2p1zcs, &
		de_chisq_g06, de_chisq_jla, de_chisq_bao_ver3
	
	!bao chisq functions
	PUBLIC :: de_chisq_6dFGS!, de_chisq_2dFGS
	PUBLIC :: de_chisq_SDSSDR7_old, de_chisq_SDSSDR7_new, de_chisq_BOSSDR9, de_chisq_dr11
	PUBLIC :: de_chisq_WigZ_rsDV, de_chisq_WigZ_A, de_chisq_impwig
	PUBLIC :: de_chisq_bao_ver1, de_chisq_bao_ver2
	
	!CMB chisq functions
	PUBLIC :: de_chisq_wmap7, de_chisq_wmap9, de_chisq_planck
	
	!H0 chisq functions
	PUBLIC :: de_chisq_h_Riess, de_chisq_h_Carnegie


!Some useful constants
	DOUBLE PRECISION,   PARAMETER  ::  logzero = 1.0d10
	DOUBLE PRECISION,   PARAMETER  ::  zfacsq = 25.0d0/(LOG(10.0d0))**2

!Some snls parameters
	DOUBLE PRECISION,   PARAMETER  ::  de_alphatol = 1.0d-10, betatol = 1.0d-10
	DOUBLE PRECISION,   PARAMETER  ::  scriptmcut = 10.0d0
	DOUBLE PRECISION,   PARAMETER  ::  intrinsicsq(4) = (/ 0.0675d0**2, 0.1133d0**2, 0.0815d0**2, 0.0989d0**2 /)
	DOUBLE PRECISION,   PARAMETER  ::  pecz = 0.0005d0
	CHARACTER,          PARAMETER  ::  uplo = 'U' !For LAPACK

!Supernova data TYPE
	TYPE, PRIVATE :: supernova
		CHARACTER(LEN=20) :: name  !The name of the SN
		DOUBLE PRECISION  :: zhel, zcmb    !The heliocentric and CMB frame redshIFts
		DOUBLE PRECISION  :: z_var         !The variance of the redshIFt
		DOUBLE PRECISION  :: mag           !The K-corrected peak magnitude
		DOUBLE PRECISION  :: mag_var       !The variance of mag
		DOUBLE PRECISION  :: stretch       !The light-curve fit stretch PARAMETER
		DOUBLE PRECISION  :: stretch_var   !The variance in the stretch
		DOUBLE PRECISION  :: colour        !The colour of the SN
		DOUBLE PRECISION  :: colour_var    !The variance of colour
		DOUBLE PRECISION  :: thirdvar      !Third variable for scripm split
		DOUBLE PRECISION  :: thirdvar_var  !Variance in thirdvar
		DOUBLE PRECISION  :: cov_mag_stretch !Covariance between mag and stretch
		DOUBLE PRECISION  :: cov_mag_colour  !Covariance between mag and colour
		DOUBLE PRECISION  :: cov_stretch_colour !Covariance between stretch and colour
		INTEGER           :: dataset       !Subset identIFier IF subset depENDent intrinsic disp is used
	END TYPE supernova

!---------------------------------------------------
!!! SNLS3
!Number of sn
	INTEGER, PARAMETER   ::   nsnls3  =  472

!Supernova data
	TYPE( supernova )    ::   sndata(nsnls3)

!The covarariance matrix 
	DOUBLE PRECISION, PRIVATE     ::   m_covmat(nsnls3,nsnls3), s_covmat(nsnls3,nsnls3)
	DOUBLE PRECISION, PRIVATE     ::   c_covmat(nsnls3,nsnls3), m_s_covmat(nsnls3,nsnls3)
	DOUBLE PRECISION, PRIVATE     ::   m_c_covmat(nsnls3,nsnls3), s_c_covmat(nsnls3,nsnls3)

!The inverse covariance matrix
	DOUBLE PRECISION, PRIVATE     ::   invcovmat(nsnls3,nsnls3)

!Some useful arrays
	DOUBLE PRECISION, PRIVATE     ::   pre_vars(nsnls3), A1(nsnls3), A2(nsnls3)
	DOUBLE PRECISION, PRIVATE     ::   snls3_lumdists(nsnls3), diffmag(nsnls3), invvars(nsnls3)
  
!Some useful variables
	LOGICAL, PUBLIC     :: snls_read       =  .FALSE.
	LOGICAL, PUBLIC     :: snls_prepped    =  .FALSE.
	LOGICAL		    :: snls_cov_check  =  .FALSE.
	LOGICAL, PRIVATE    :: first_inversion =  .TRUE.
	LOGICAL, PRIVATE    :: use_ab_prev     =  .TRUE.
	DOUBLE PRECISION    :: alpha_prev, beta_prev

!---------------------------------------------------
!!! Gold06
	INTEGER, PRIVATE, PARAMETER :: g06num = 292
	DOUBLE PRECISION :: g06z(g06num), g06mu(g06num), g06dmu(g06num)
	INTEGER :: g06mark(g06num) ! 1 means gold; 0 means silver
	INTEGER :: g06goodbad(g06num) ! 1 means good; 0 means bad
	LOGICAL, PUBLIC :: g06_read = .false.		
				


!---------------------------------------------------
!!! Union2.1
!Union2.1 data
	INTEGER, PRIVATE, PARAMETER   	:: union2p1_num  =  580
	DOUBLE PRECISION, PRIVATE 	:: union2p1_z(union2p1_num), union2p1_moduli(union2p1_num), union2p1_modulierr(union2p1_num), union2p1_plow(union2p1_num), union2p1_sumninv
	DOUBLE PRECISION, PRIVATE 	:: union2p1_Ninv(union2p1_num,union2p1_Num)

!Union2 data
	INTEGER, PRIVATE, PARAMETER   	:: union2_num  =  558
	DOUBLE PRECISION, PRIVATE 	:: union2_z(union2_num), union2_moduli(union2_num), union2_modulierr(union2_num), union2_plow(union2_num), union2_sumninv
	DOUBLE PRECISION, PRIVATE 	:: union2_Ninv(union2_num,union2_Num)



!Some useful variables
	LOGICAL, PUBLIC	    	:: union2p1_syscovmat  =  .TRUE.
	LOGICAL, PUBLIC     	:: union2p1_inited     =  .FALSE.
	DOUBLE PRECISION    	:: union2p1_lumdists(union2p1_num), union2p1_diffs(union2p1_num)
	INTEGER, PRIVATE    	:: count_chisq     =   0

	LOGICAL, PUBLIC	    	:: union2_syscovmat  =  .TRUE.
	LOGICAL, PUBLIC     	:: union2_inited     =  .FALSE.
	DOUBLE PRECISION    	:: union2_lumdists(union2p1_num), union2_diffs(union2_num)



  CONTAINS


  !------------------------------------------
  ! The total chisq. 
  !------------------------------------------
	DOUBLE PRECISION FUNCTION de_chisq_all()
	
		IF(de_CP%Odm0<0.0d0 .or. de_CP%Odm0>1.0d0 &
			.or. de_CP%h<0.5d0 .or. de_CP%h>1.0d0) THEN
			de_chisq_all = logzero
			RETURN
		ENDIF
		
		CALL de_Init()

		de_chisq_all =  de_chisq_union2p1() + de_chisq_snls3()  + de_chisq_wmap7() + de_chisq_h_Carnegie() + de_chisq_h_Riess() + de_chisq_bao_ver1() + de_chisq_bao_ver2()
		!+ de_chisq_union2p1() + chisq_bao_old()
		If(de_chisq_all .ge. logzero) de_chisq_all = logzero
		count_chisq = count_chisq + 1
	END FUNCTION de_chisq_all

!#####################################################
!#####################################################
  ! Supernovae chisq
!#####################################################
!#####################################################

  ! Gold06
  ! Gold06 read_in
	SUBROUTINE g06_init()
		INTEGER :: i, i_gold, i_silver
		CHARACTER(LEN=100) :: tmpstr, tmpstr1
		OPEN(UNIT=870,FILE=trim(adjustl(de_data_path))//'sn_Gold06.txt')
		i_gold = 0; i_silver = 0
		DO i = 1, g06num
			read(870, *) tmpstr, g06z(i), g06mu(i), g06dmu(i), tmpstr, tmpstr1
			if(trim(adjustl(tmpstr)).eq.'Gold') then
				g06mark(i) = 1
				i_gold = i_gold + 1 
			elseif(trim(adjustl(tmpstr)).eq.'Silver') then
				g06mark(i) = 0
				i_silver = i_silver + 1
			else
				print *, 'ERROR!!! ', i, trim(adjustl(tmpstr))
			endif
			if(trim(adjustl(tmpstr1)).eq.'Bad') then
				g06goodbad(i) = 0
				print *, 'bad found:', i, g06z(i)
			elseif(trim(adjustl(tmpstr1)).eq.'\') then
				g06goodbad(i) = 1
			else
				print *, 'ERROR!!! ', i, trim(adjustl(tmpstr1))
			endif
		enddo
		write(*,'(A,i4,A,i4,i4)'), 'Finishing in reading in ', g06num, 'G06 SNe; # of gold/silver = ', i_gold, i_silver
		g06_read = .true.
	END SUBROUTINE g06_init
  ! Gold06 Chisq
	DOUBLE PRECISION FUNCTION de_chisq_g06(ignoresilver, nobadsilver)
		! Dummy 
		LOGICAL, INTENT(IN) :: ignoresilver 
		LOGICAL, INTENT(IN), optional :: nobadsilver
		! Local
		INTEGER :: i, num_bin
		DOUBLE PRECISION :: z,mu,dmu,lumdist,muth,diffmu,dmufac, A, B, C

		IF(g06_read.eq..FALSE.) then
			CALL g06_init()
		ENDIF

		A = 0
		B = 0
		C = 0
		DO i = 1, g06num
			if(ignoresilver .and. g06mark(i).eq.0) then
				cycle ! ignore silver snIa
			endif
			if(present(nobadsilver)) then
				if(nobadsilver.and.g06goodbad(i).eq.0) then
					cycle
				endif
			endif
			z = g06z(i); mu = g06mu(i); dmu = g06dmu(i)
			num_bin = Ceiling(z*32.0d0)
			lumdist =  de_Simpson(de_inv_e, 0.0d0, z, num_bin)
			muth = 5.0d0 * LOG10(  (1.0d0+z) * lumdist )
			diffmu = mu - muth 
			dmufac = 1.0d0 / dmu**2.0
			A = A + diffmu**2.0 * dmufac
			B = B + diffmu * dmufac
			C = C + dmufac
		ENDDO

		de_chisq_g06 = A - B*B / C
	END FUNCTION de_chisq_g06


  ! SNLS3
  ! SNLS3 read_in
	  !---------------------------------------------------------------
	  ! READ in the SNLS data. Including 1 data FILE, 6 covmats.     
	  ! The global variable snls_READ will be setted to .TRUE.
	  !---------------------------------------------------------------
	SUBROUTINE read_snls_dataset   
		INTEGER  :: I
		DOUBLE PRECISION :: dz, dm, ds, dc, dt

		OPEN(UNIT=70,FILE=trim(adjustl(de_data_path))//'snls_3rdyear_lcparams.txt')
		OPEN(UNIT=71,FILE=trim(adjustl(de_data_path))//'snls3_v0_covmatrix.dat')
		OPEN(UNIT=72,FILE=trim(adjustl(de_data_path))//'snls3_va_covmatrix.dat')
		OPEN(UNIT=73,FILE=trim(adjustl(de_data_path))//'snls3_vb_covmatrix.dat')
		OPEN(UNIT=74,FILE=trim(adjustl(de_data_path))//'snls3_v0a_covmatrix.dat')
		OPEN(UNIT=75,FILE=trim(adjustl(de_data_path))//'snls3_v0b_covmatrix.dat')
		OPEN(UNIT=76,FILE=trim(adjustl(de_data_path))//'snls3_vab_covmatrix.dat')

		DO I=1, nsnls3
			READ(70,*)   &
			sndata(I)%name, sndata(I)%zcmb, sndata(I)%zhel,&
			dz, sndata(I)%mag, dm, sndata(I)%stretch, ds, &
			sndata(I)%colour,dc,sndata(I)%thirdvar, dt,&
			sndata(I)%cov_mag_stretch,&
			sndata(I)%cov_mag_colour,sndata(I)%cov_stretch_colour,&
			sndata(I)%dataset

			sndata(I)%z_var = dz**2
			sndata(I)%mag_var = dm**2
			sndata(I)%stretch_var = ds**2
			sndata(I)%colour_var = dc**2
			sndata(I)%thirdvar_var = dt**2

			READ(71,*) m_covmat(I,1:nsnls3)
			READ(72,*) s_covmat(I,1:nsnls3)
			READ(73,*) c_covmat(I,1:nsnls3)
			READ(74,*) m_s_covmat(I,1:nsnls3)
			READ(75,*) m_c_covmat(I,1:nsnls3)
			READ(76,*) s_c_covmat(I,1:nsnls3)
		ENDDO
		CLOSE(70); CLOSE(71); CLOSE(72); CLOSE(73); CLOSE(74); CLOSE(75); CLOSE(76);

		snls_read = .TRUE.
		first_inversion = .TRUE.
		snls_prepped = .FALSE.
		WRITE(*,*) "snls read_in complete"   
	END SUBROUTINE read_snls_dataset
  ! SNLS3 preparation
	  !----------------------------------------------------------------
	  ! Prepares the data for fitting by pre-calculating the parts of  
	  !  the errors that can be done ahead of time.                   
	  ! READ_snls_dataset() must have been CALLed before CALLing this       
	  ! The global variable snls_prepped will be setted to .TRUE.
	  !----------------------------------------------------------------
	SUBROUTINE snls_prep()
		INTEGER :: i

		pre_vars = sndata%mag_var + intrinsicsq(sndata%dataset+1)
		DO I=i,nsnls3
		pre_vars(i) = pre_vars(i) + &
			zfacsq * ( sndata(i)%z_var + pecz**2 ) * &
          		( (1.0d0 + sndata(i)%zcmb)/&
		        (sndata(i)%zcmb*(1.0d0 + 0.5d0*sndata(i)%zcmb)) )**2
		ENDDO

		DO i=1, nsnls3
			IF (sndata(i)%thirdvar .LE. scriptmcut ) THEN
			A1(i) = 1.0d0;  A2(i) = 0.0d0
			ELSE
			A1(i) = 0.0d0;  A2(i) = 1.0d0
			ENDIF
		ENDDO
!		WRITE(*,*) "test pre_vars 2: ", pre_vars(1), pre_vars(100), pre_vars(300)
		WRITE(*,*) "snls preparation complete"   
		snls_prepped = .TRUE.
		first_inversion = .TRUE.
	END SUBROUTINE snls_prep 
  ! SNLS3 inv covmat
	SUBROUTINE inv_cov_mat( alpha, beta, STATUS )     
		DOUBLE PRECISION  ::  alpha, beta
		INTEGER :: STATUS, i, j
		DOUBLE PRECISION  ::  alphasq, betasq, alphabeta 

!		IF (use_ab_prev .EQ. .TRUE. .AND. .NOT. first_inversion .AND. (ABS(alpha-alpha_prev) .LT. alphatol) .AND. &
!			( ABS(beta-beta_prev) .LT. betatol )) THEN           !Previous invcovmatrix is close enough
!			STATUS = 0
!			RETURN
!		ENDIF

		alphasq = alpha * alpha
		betasq = beta * beta
		alphabeta = alpha * beta

!		Build the covariance matrix    
 		invcovmat = m_covmat &
			+ alphasq * s_covmat &
			+ betasq * c_covmat &
			+ 2.0d0 * alpha * m_s_covmat &  
			- 2.0d0 * beta * m_c_covmat &
			- 2.0d0 * alphabeta * s_c_covmat

		if(snls_cov_check) then
			WRITE(*,*) "SNLS covmat check:"
			WRITE(*,*) m_covmat(1,1), m_covmat(235,256), m_covmat(nsnls3, nsnls3)			
			WRITE(*,*) s_covmat(1,1), s_covmat(235,256), s_covmat(nsnls3, nsnls3)
			WRITE(*,*) c_covmat(1,1), c_covmat(235,256), c_covmat(nsnls3, nsnls3)
			WRITE(*,*) m_s_covmat(1,1), m_s_covmat(235,256), m_s_covmat(nsnls3, nsnls3)	
			WRITE(*,*) m_c_covmat(1,1), m_c_covmat(235,256), m_c_covmat(nsnls3, nsnls3)	
			WRITE(*,*) s_c_covmat(1,1), s_c_covmat(235,256), s_c_covmat(nsnls3, nsnls3)	
			WRITE(*,*) invcovmat(1,1), invcovmat(235,256), invcovmat(nsnls3, nsnls3)
		endif
		
!		Update the diagonal terms
		DO I=1, nsnls3
			invcovmat(I,I) = invcovmat(I,I) + pre_vars(I) &
				+ alphasq * sndata(I)%stretch_var &
				+ betasq  * sndata(I)%colour_var &
				+ 2.0d0 * alpha * sndata(I)%cov_mag_stretch &
				- 2.0d0 * beta * sndata(I)%cov_mag_colour &
				- 2.0d0 * alphabeta * sndata(I)%cov_stretch_colour
		ENDDO
	
		CALL DPOTRF(uplo,nsnls3,invcovmat,nsnls3,STATUS)
		CALL DPOTRI(uplo,nsnls3,invcovmat,nsnls3,STATUS)
	
		first_inversion = .FALSE.
		alpha_prev = alpha
		beta_prev  = beta
	END SUBROUTINE inv_cov_mat
  ! SNLS3 chisq
	DOUBLE PRECISION FUNCTION de_chisq_snls3()
		DOUBLE PRECISION :: alpha, alphasq, beta, betasq, alphabeta
		DOUBLE PRECISION :: ogamma, Omegar, Omegade1, Omegade2
		DOUBLE PRECISION :: now_broader
		DOUBLE PRECISION :: estimated_scriptm, wtval
		DOUBLE PRECISION :: zhel, zcmb, chisq
		DOUBLE PRECISION :: amarg_A, amarg_B, amarg_C, amarg_D, amarg_E, amarg_F, tempG
		INTEGER :: i, num_bin, STATUS

		IF(snls_read .EQ. .FALSE.)	CALL read_snls_dataset()
		IF(snls_prepped .EQ. .FALSE. )  CALL snls_prep()

		alpha = de_CP%alpha
		beta  = de_CP%beta

		alphasq   = alpha*alpha
		betasq    = beta*beta
		alphabeta = alpha*beta

		DO I = 1, nsnls3
			zcmb = sndata(I)%zcmb
			num_bin = Ceiling(zcmb*32.0d0)
			snls3_lumdists(i) = de_Simpson(de_inv_e, 0.0d0, zcmb, num_bin)
			snls3_lumdists(i) = de_fk(snls3_lumdists(i))
			zhel = sndata(I)%zhel
			snls3_lumdists(I) = 5.0d0 * LOG10( (de_const_c/(100.0d0)) * (1.0d0+zhel) * snls3_lumdists(i) )
		ENDDO


!		Calculate estimated_scriptm and dIFfmag
		invvars = 1.0d0 / ( pre_vars + alphasq * sndata%stretch_var &
			+ betasq * sndata%colour_var &
			+ 2.0d0 * alpha * sndata%cov_mag_stretch &
			- 2.0d0 * beta * sndata%cov_mag_colour &
			- 2.0d0 * alphabeta * sndata%cov_stretch_colour )
!		WRITE(*,*) alpha, beta, alphasq, betasq, alphabeta

		wtval = SUM( invvars )
		estimated_scriptm= SUM( (sndata%mag - snls3_lumdists)*invvars ) / wtval
		diffmag = sndata%mag - snls3_lumdists + alpha*( sndata%stretch - 1.0d0 ) &
			- beta * sndata%colour - estimated_scriptm 

		CALL inv_cov_mat(alpha,beta,STATUS)

!		Now find the amarg_ PARAMETERs
!		We re-use the invvars variable to hold the intermediate product
!		which is sort of naughty
!		invvars = V^-1 * dIFfmag (invvars = 1.0*invcovmat*dIFfmag+0*invvars)
		CALL DSYMV(uplo,nsnls3,1.0d0,invcovmat,nsnls3,diffmag,1,0.0d0,invvars,1)

		amarg_A = DOT_PRODUCT( diffmag, invvars ) ! dIFfmag*V^-1*dIFfmag
		amarg_B = DOT_PRODUCT( invvars, A1 ) !dIFfmag*V^-1*A1
		amarg_C = DOT_PRODUCT( invvars, A2 ) !dIFfmag*V^-1*A2

!		Be naughty again and stick V^-1 * A1 in invvars
		CALL DSYMV(uplo,nsnls3,1.0d0,invcovmat,nsnls3,A1,1,0.0d0,invvars,1)
		amarg_D = DOT_PRODUCT( invvars, A2 ) !A2*V^-1*A1
		amarg_E = DOT_PRODUCT( invvars, A1 ) !A1*V^-1*A1

!		now V^-1 * A2
		CALL DSYMV(uplo,nsnls3,1.0d0,invcovmat,nsnls3,A2,1,0.0d0,invvars,1)
		amarg_F = DOT_PRODUCT( invvars, A2 ) !A2*V^-1*A2
		tempG = amarg_F - amarg_D*amarg_D/amarg_E; !g/e

!		Marginalized chisq
		chisq = amarg_A + LOG( amarg_E*de_inv_twopi ) + &
			LOG( tempG * de_inv_twopi ) - amarg_C*amarg_C/tempG - &
			amarg_B*amarg_B*amarg_F / ( amarg_E*tempG ) + &
			2.0d0*amarg_B*amarg_C*amarg_D/(amarg_E*tempG )

		de_chisq_snls3 = chisq 
		IF(pr_chisq_info) THEN
			WRITE(*,*) "     chisq_snls3 (lnlike) = ", de_chisq_snls3, de_chisq_snls3/2.0
			WRITE(*,*)
		ENDIF
	END FUNCTION de_chisq_snls3
  ! Union2.1
  ! Union2.1 Read in
	SUBROUTINE union2p1_init()
		CHARACTER(LEN=20) :: union2p1_name
		INTEGER :: i

		OPEN(UNIT=101, FILE=trim(adjustl(de_data_path))//'sn_z_mu_dmu_plow_union2.1.txt')
		DO i = 1, union2p1_num
			READ(101,*) union2p1_name, union2p1_z(i), union2p1_moduli(i), union2p1_modulierr(i), union2p1_plow(i)
		ENDDO
		CLOSE(101)

		IF(union2p1_syscovmat) THEN
			OPEN(UNIT=102, FILE=trim(adjustl(de_data_path))//'sn_wmat_sys_union2.1.txt')
			WRITE(*,*) "We are using sys in SN data."
			ELSE
			OPEN(UNIT=102, FILE=trim(adjustl(de_data_path))//'sn_wmat_nosys_union2.1.txt')
		ENDIF

		DO i = 1, union2p1_num
			READ(102,*) union2p1_ninv(i,1:union2p1_num)
		ENDDO
		CLOSE(102)
!		union2p1_sumninv = sum(union2p1_ninv)
	END SUBROUTINE union2p1_init
  ! Union2.1 chisq
	DOUBLE PRECISION FUNCTION de_chisq_union2p1()
		! Local
		INTEGER :: i, num_bin
		
		IF(union2p1_inited .EQ. .FALSE.) THEN
			CALL union2p1_init()
			WRITE(*,*) "Union2.1 read in completed"
			union2p1_inited = .TRUE.
		ENDIF

		DO i = 1, union2p1_num
			num_bin = Ceiling(union2p1_z(i)*16.0d0)
			union2p1_lumdists(i) = de_Simpson(de_inv_e, 0.0d0, union2p1_z(i), num_bin) 
			union2p1_lumdists(i) = de_fk(union2p1_lumdists(i))/ de_CP%H0
			union2p1_lumdists(i) = 5.0 * LOG10((1.0d0+union2p1_z(i))*union2p1_lumdists(i)) + 25	
			union2p1_diffs(i) = union2p1_lumdists(i) - union2p1_moduli(i)
		ENDDO

!		write(*,'(<560>(f10.5,1x))') union2p1_diffs
		de_chisq_union2p1 = dot_product(union2p1_diffs,matmul(union2p1_ninv,union2p1_diffs)) 

		IF(pr_chisq_info) THEN
			WRITE(*,*) "     chisq_Union2.1 = ", de_chisq_union2p1
			WRITE(*,*)
		ENDIF
	END FUNCTION de_chisq_union2p1


	DOUBLE PRECISION FUNCTION de_chisq_union2p1zcs(zcut_max, SNIaM)
		! Dummy
		DOUBLE PRECISION, INTENT(IN) :: zcut_max, SNIaM
		! Local
		INTEGER :: i, num_bin
		
		IF(union2p1_inited .EQ. .FALSE.) THEN
			CALL union2p1_init()
			WRITE(*,*) "Union2.1 read in completed"
			union2p1_inited = .TRUE.
		ENDIF

		DO i = 1, union2p1_num
			num_bin = Ceiling(union2p1_z(i)*16.0d0)
			union2p1_lumdists(i) = de_Simpson(de_inv_e, 0.0d0, union2p1_z(i), num_bin) 
			union2p1_lumdists(i) = de_fk(union2p1_lumdists(i))
			union2p1_lumdists(i) = 5.0 * LOG10((1.0d0+union2p1_z(i))*union2p1_lumdists(i)) + 25 + SNIaM
			union2p1_diffs(i) = union2p1_lumdists(i) - union2p1_moduli(i)
			if(union2p1_z(i) > zcut_max) then
				union2p1_diffs(i) = 0.0
			endif
		ENDDO
		de_chisq_union2p1zcs = dot_product(union2p1_diffs,matmul(union2p1_ninv,union2p1_diffs)) 

	END FUNCTION de_chisq_union2p1zcs
	! chisq of Union2.1, with maximal cut of redshift
	! Actually it is very demanding; I do a search in H direction
	DOUBLE PRECISION FUNCTION de_chisq_union2p1zc(zcut_max)
		! Dummy
		DOUBLE PRECISION, INTENT(IN) :: zcut_max
		! Local
		DOUBLE PRECISION :: H1, H2, Hmid, chisqA, chisqB, chisqmid, deltaH = 0.0001

		H1 = -100; H2 = 100.0;
		do while(abs(H2-H1).ge.3.0*deltaH)
			Hmid = (H1+H2)/2.0
			chisqA = de_chisq_union2p1zcs(zcut_max, Hmid)
			chisqB = de_chisq_union2p1zcs(zcut_max, Hmid+max((H2-H1)*0.01,deltaH))
!			print *, real(H1), real(H2), real(chisqA), real(chisqB)
			if(chisqB>chisqA) then
				H2 = Hmid
			else
				H1 = Hmid
			endif
		enddo
		de_chisq_union2p1zc = min(chisqA, chisqB, de_chisq_union2p1zcs(zcut_max, Hmid+deltaH/2.0))
	end function de_chisq_union2p1zc



!#####################################################
!#####################################################
  ! CMB chisq
!#####################################################
!#####################################################
  ! WMAP7
	DOUBLE PRECISION FUNCTION de_chisq_wmap7()
		DOUBLE PRECISION :: zstar, R, lA
		DOUBLE PRECISION :: g1, g2, DAzstar, rszstar
		DOUBLE PRECISION :: zstar_ML, R_ML, lA_ML, DVecCMB(3), CovCMB(3,3)

		CovCMB(1,1)=2.305d0;   CovCMB(1,2)=29.698d0;     CovCMB(1,3)=-1.333d0;
		CovCMB(2,1)=29.698d0;  CovCMB(2,2)=6825.27d0;    CovCMB(2,3)=-113.180d0;
		CovCMB(3,1)=-1.333d0;  CovCMB(3,2)=-113.180d0;   CovCMB(3,3)=3.414d0;
		lA_ML=302.09d0; R_ML=1.725d0; zstar_ML=1091.3d0;
   
		DVecCMB   = (/de_CP%lA-lA_ML, de_CP%R-R_ML, de_CP%zstar-zstar_ML/)
		de_chisq_wmap7 = DOT_PRODUCT(DVecCMB,matmul(CovCMB,DVecCMB))
		
		IF(pr_chisq_info) THEN
			WRITE(*,*) "     chisq_wmap7 (lnlike) = ", de_chisq_wmap7, de_chisq_wmap7/2.0
			WRITE(*,*)
		ENDIF
	END FUNCTION de_chisq_wmap7
  ! planck 1yr
	DOUBLE PRECISION FUNCTION de_chisq_planck()   !de_chisq_planck   !de_chisq_wmap7()
		DOUBLE PRECISION :: lA, R, Ob0hsq    !zstar
		DOUBLE PRECISION :: g1, g2, DAzstar, rszstar
		DOUBLE PRECISION :: lA_ML, R_ML, Ob0hsq_ML, DVecCMB(3), CovCMB(3,3)   !zstar_ML

		lA     = de_CP%lA
		R      = de_CP%R 
		Ob0hsq = de_CP%Ob0hsq

		CovCMB(1,1)=43.0179647667078;   CovCMB(1,2)=-366.771828690038;     CovCMB(1,3)=2972.52785899757;   !-----------------
		CovCMB(2,1)=-366.771828690038;  CovCMB(2,2)=24872.6578215152;      CovCMB(2,3)=446498.465930772;   ! planck's covmat
		CovCMB(3,1)=2972.52785899758;   CovCMB(3,2)=446498.465930772;      CovCMB(3,3)=21554701.3253821;   ! and mean value
		lA_ML=301.57d0; R_ML=1.7407d0; Ob0hsq_ML=0.02228d0                                                 !----------------- 

		DVecCMB         = (/lA-lA_ML, R-R_ML, Ob0hsq-Ob0hsq_ML/)
		de_chisq_planck = DOT_PRODUCT(DVecCMB,matmul(CovCMB,DVecCMB))

		IF(pr_chisq_info) THEN
			WRITE(*,*) "     chisq_planck (lnlike) = ", de_chisq_planck, de_chisq_planck/2.0
			WRITE(*,*)
		ENDIF

		!write(*,*) "la,r,ob0hsq", lA, R, Ob0hsq    !bossli
	END FUNCTION de_chisq_planck   !de_chisq_planck   !de_chisq_wmap7
  ! WMAP9
	DOUBLE PRECISION FUNCTION de_chisq_wmap9()   !de_chisq_planck   !de_chisq_wmap7()
		DOUBLE PRECISION :: lA, R, Ob0hsq    !zstar
		DOUBLE PRECISION :: g1, g2, DAzstar, rszstar
		DOUBLE PRECISION :: lA_ML, R_ML, Ob0hsq_ML, DVecCMB(3), CovCMB(3,3)   !zstar_ML

		lA     = de_CP%lA
		R      = de_CP%R 
		Ob0hsq = de_CP%Ob0hsq

		CovCMB(1,1)=3.68712321402858;   CovCMB(1,2)=-14.1725878878498;     CovCMB(1,3)=2566.01640514695;   !-----------------
		CovCMB(2,1)=-14.1725878878498;  CovCMB(2,2)=5179.04970797377;      CovCMB(2,3)=73212.4449017310;   ! wmap9's covmat
		CovCMB(3,1)=2566.01640514695;   CovCMB(3,2)=73212.4449017310;      CovCMB(3,3)=6692540.05113778;   ! and mean value
		lA_ML=302.02d0; R_ML=1.7327d0; Ob0hsq_ML=0.02260d0                                                 !----------------- 

		DVecCMB         = (/lA-lA_ML, R-R_ML, Ob0hsq-Ob0hsq_ML/)
		de_chisq_wmap9 = DOT_PRODUCT(DVecCMB,matmul(CovCMB,DVecCMB))

		IF(pr_chisq_info) THEN
			WRITE(*,*) "     chisq_wmap9 (lnlike) = ", de_chisq_wmap9, de_chisq_wmap9/2.0
			WRITE(*,*)
		ENDIF

		!write(*,*) "la,r,ob0hsq", lA, R, Ob0hsq    !bossli
	END FUNCTION de_chisq_wmap9   !de_chisq_planck   !de_chisq_wmap7

!#####################################################
!#####################################################
  ! BAO chisq
!#####################################################
!#####################################################
  ! 6dFGS rstodv(0.106) = 0.336 +/- 0.015 arXiv:1106.3366
	DOUBLE PRECISION FUNCTION de_chisq_6dFGS()
		de_chisq_6dFGS = ((de_rstodv(0.106d0) - 0.336d0)/0.015d0)**2.0 
		IF(pr_chisq_info) THEN
			WRITE(*,*) "        de_rstodv(0.106) = ", &
				de_rstodv(0.106d0)
			WRITE(*,*) "     chisq_6dFGS (lnlike) = ", &
				de_chisq_6dFGS, de_chisq_6dFGS/2.0
			WRITE(*,*)
		ENDIF
	END FUNCTION de_chisq_6dFGS
  ! BOSS DR9 dvtors(0.57) = 13.67 +/- 0.22 arXiv:1203.6594
  	DOUBLE PRECISION FUNCTION de_chisq_BOSSDR9()
  		DOUBLE PRECISION :: rstoDV0p57
  		
  		rstoDV0p57 = de_rstodv(0.57d0)
  		de_chisq_BOSSDR9 = ((1.0/rstoDV0p57 - 13.67d0)/0.22d0)**2.0
  		
  			
  		IF(pr_chisq_info) THEN
  			WRITE(*,*) "        de_rstodv(0.57) = ", &
  				rstoDV0p57
  			WRITE(*,*) "     chisq_BOSSDR9 (lnlike) = ", &
  				de_chisq_BOSSDR9, de_chisq_BOSSDR9/2.0
  			WRITE(*,*)
  		ENDIF
  	END FUNCTION de_chisq_BOSSDR9
  ! SDSS DR7 reconstructed dvtors(0.35) = 8.88 +/- 0.17 arXiv:1202.0090
	DOUBLE PRECISION FUNCTION de_chisq_SDSSDR7_new()
		de_chisq_SDSSDR7_new = ((1.0d0/de_rstodv(0.35d0)-8.88d0)/0.17d0)**2.0 
		IF(pr_chisq_info) THEN
			WRITE(*,*) "        de_rstodv(0.35) = ", &
				de_rstodv(0.35d0)
			WRITE(*,*) "     chisq_SDSSDR7_new (lnlike) = ", &
				de_chisq_SDSSDR7_new, de_chisq_SDSSDR7_new/2.0
			WRITE(*,*)
		ENDIF
	END FUNCTION de_chisq_SDSSDR7_new  
  ! SDSS DR7 arXiv:0907.1660 see wmap9 paper arXiv:1212.5226
  	DOUBLE PRECISION FUNCTION de_chisq_SDSSDR7_old()
  		DOUBLE PRECISION :: baodiff(2), covinv(2,2)
  		
  		baodiff(1) = de_rstodv(0.2d0)-0.1905d0
  		baodiff(2) = de_rstodv(0.35d0)-0.1097d0
  		
  		covinv(1,1) =  30124.0d0
                covinv(2,2) =  86977.0d0 
                
		covinv(1,2) = -17227.0d0
  		covinv(2,1) = covinv(1,2)
  		
  		de_chisq_SDSSDR7_old = dot_product(baodiff,matmul(covinv,baodiff))
  		
  		IF(pr_chisq_info) THEN
			WRITE(*,*) "        de_rstodv(0.2) = ", &
				de_rstodv(0.2d0)			
			WRITE(*,*) "        de_rstodv(0.35) = ", &
				de_rstodv(0.35d0)
			WRITE(*,*) "     chisq_SDSSDR7_old (lnlike) = ", &
				de_chisq_SDSSDR7_old, de_chisq_SDSSDR7_old/2.0
			WRITE(*,*)			
		ENDIF
	END FUNCTION de_chisq_SDSSDR7_old
  ! WiggleZ rstoDV see wmap9 paper arXiv:1212.5226
  	DOUBLE PRECISION FUNCTION de_chisq_WigZ_rsDV()
  		DOUBLE PRECISION :: baodiff(3), covinv(3,3)

		baodiff(1) = de_rstodv(0.44d0) - 0.0916d0
		baodiff(2) = de_rstodv(0.60d0) - 0.0726d0
		baodiff(3) = de_rstodv(0.73d0) - 0.0592d0

		covinv(1,1) = 24532.1d0
		covinv(2,2) = 134598.4d0
		covinv(3,3) = 128837.6d0
  			
  		covinv(1,2) = -25137.7d0
		covinv(2,1) = covinv(1,2)

		covinv(1,3) = 12099.1d0
		covinv(3,1) = covinv(1,3)

		covinv(2,3) = -64783.9d0
		covinv(3,2) = covinv(2,3)
  		
  		de_chisq_WigZ_rsDV = dot_product(baodiff,matmul(covinv,baodiff))
  		
		IF(pr_chisq_info) THEN
			WRITE(*,*) "        de_rstodv(0.44) = ", &
				de_rstodv(0.44d0)			
			WRITE(*,*) "        de_rstodv(0.60) = ", &
				de_rstodv(0.6d0)
			WRITE(*,*) "        de_rstodv(0.73) = ", &
				de_rstodv(0.73d0)			
			WRITE(*,*) "     chisq_WiggleZ_rstoDV (lnlike) = ", &
				de_chisq_WigZ_rsDV, de_chisq_WigZ_rsDV/2.0
			WRITE(*,*)			
		ENDIF
	END FUNCTION de_chisq_WigZ_rsDV
  ! WiggleZ A  see wmap9 paper arXiv:1212.5226
  	DOUBLE PRECISION FUNCTION de_chisq_WigZ_A()
  		DOUBLE PRECISION :: baodiff(3), covinv(3,3)

		baodiff(1) = de_A_bao(0.44d0) - 0.474d0
		baodiff(2) = de_A_bao(0.60d0) - 0.442d0
		baodiff(3) = de_A_bao(0.73d0) - 0.424d0
  		
  		covinv(1,1) = 1040.3d0
		covinv(2,2) = 3720.3d0
		covinv(3,3) = 2914.9d0
  		
  		covinv(1,2) = -807.5d0
		covinv(2,1) = covinv(1,2)

		covinv(1,3) = 336.8d0
		covinv(3,1) = covinv(1,3)

		covinv(2,3) = -1551.9d0
		covinv(3,2) = covinv(2,3)
  		
  		de_chisq_WigZ_A = dot_product(baodiff,matmul(covinv,baodiff))
  		
		IF(pr_chisq_info) THEN
			WRITE(*,*) "        de_A_bao(0.44) = ", &
				de_rstodv(0.44d0)			
			WRITE(*,*) "        de_A_bao(0.60) = ", &
				de_rstodv(0.6d0)
			WRITE(*,*) "        de_A_bao(0.73) = ", &
				de_rstodv(0.73d0)			
			WRITE(*,*) "     chisq_WiggleZ_A (lnlike) = ", &
				de_chisq_WigZ_A, de_chisq_WigZ_A/2.0
			WRITE(*,*)			
		ENDIF
	END FUNCTION de_chisq_WigZ_A
	
  !----------------------------------------------------------------
  ! used in wmap9 paper arXiv:1212.5226
  !  (6dFGS + WiggleZ rsDV + BOSS)
	DOUBLE PRECISION FUNCTION de_chisq_bao_ver2()
		DOUBLE PRECISION :: baodiff(6), covinv(6,6)

		baodiff(1) = de_rstodv(0.106d0) - 0.336d0
		baodiff(2) = 1.0d0 / de_rstodv(0.35d0) - 8.88d0
		baodiff(3) = 1.0d0 / de_rstodv(0.57d0) - 13.67d0
		baodiff(4) = de_rstodv(0.44d0) - 0.0916d0
		baodiff(5) = de_rstodv(0.60d0) - 0.0726d0
		baodiff(6) = de_rstodv(0.73d0) - 0.0592d0

		covinv = 0.0d0
		covinv(1,1) = 4444.4d0
		covinv(2,2) = 34.602d0
		covinv(3,3) = 20.661157d0

		covinv(4,4) = 24532.1d0
		covinv(5,5) = 134598.4d0
		covinv(6,6) = 128837.6d0

		covinv(4,5) = -25137.7d0
		covinv(5,4) = covinv(4,5)

		covinv(4,6) = 12099.1d0
		covinv(6,4) = covinv(4,6)

		covinv(5,6) = -64783.9d0
		covinv(6,5) = covinv(5,6)
		
		de_chisq_bao_ver2 = dot_product(baodiff,matmul(covinv,baodiff))
		IF(pr_chisq_info) THEN
			WRITE(*,*) "        de_rstodv(0.10) = ", de_rstodv(0.1d0)		
			WRITE(*,*) "        de_rstodv(0.35) = ", de_rstodv(0.35d0)
			WRITE(*,*) "        de_rstodv(0.57) = ", de_rstodv(0.57d0)			
			WRITE(*,*) "        de_rstodv(0.44) = ", de_rstodv(0.44d0)			
			WRITE(*,*) "        de_rstodv(0.60) = ", de_rstodv(0.6d0)
			WRITE(*,*) "        de_rstodv(0.73) = ", de_rstodv(0.73d0)			
			WRITE(*,*) "     chisq_bao_ver2 (lnlike) = ", de_chisq_bao_ver2, de_chisq_bao_ver2/2.0
			WRITE(*,*)
		ENDIF
	END FUNCTION de_chisq_bao_ver2
	
  !----------------------------------------------------------------
  ! SDSS DR7 + WiggleZ A + 6dFGS
	DOUBLE PRECISION FUNCTION de_chisq_bao_ver1()
		DOUBLE PRECISION :: z1, z2, z3, z4, z5, z6
		DOUBLE PRECISION :: p1(2), p2(3)
 		DOUBLE PRECISION :: cov1(2,2)=(/ 30124.0d0,-17227.0d0, &
                                   -17227.0d0,86977.0d0 /)
		DOUBLE PRECISION :: cov2(3,3)=(/ 1040.3d0,-807.5d0,336.8d0, &
                                      -807.5d0,3720.3d0,-1551.9d0, &
                                       336.8d0,-1551.9d0,2914.9d0 /)
		DOUBLE PRECISION :: chisqSDSS, chisqWiggleZ, chisq6dFGS
		INTEGER :: i
		
		z1 = 0.2d0
		z2 = 0.35d0
		z3 = 0.106d0
		z4 = 0.44d0
		z5 = 0.6d0
		z6 = 0.73d0

		p1 = (/ de_rstodv(z1)-0.1905d0, de_rstodv(z2)-0.1097d0 /)
		p2 = (/ de_A_bao(z4)-0.474d0, de_A_bao(z5)-0.442d0, de_A_bao(z6)-0.424d0 /)

		chisqSDSS = DOT_PRODUCT(p1, matmul(Cov1,p1))
		chisqWiggleZ = DOT_PRODUCT(p2, matmul(Cov2,p2))
		chisq6dFGS = (de_rstodv(z3)-0.336d0)**2.0d0/0.015d0**2.0d0
		
		de_chisq_bao_ver1  = chisqSDSS + chisqWiggleZ + chisq6dFGS
		
		IF(pr_chisq_info) THEN
			
			WRITE(*,*) "     SDSS DR7 chisq (lnlike) = ", chisqSDSS, chisqSDSS/2.0
			WRITE(*,*) "     WiggleZ chisq (lnlike)  = ", chisqWiggleZ, chisqWiggleZ/2.0
			WRITE(*,*) "     6dFGS chisq (lnlike)    = ", chisq6dFGS, chisq6dFGS/2.0
			WRITE(*,*) "     chisq_bao_ver1 (lnlike)  = ", de_chisq_bao_ver1,  de_chisq_bao_ver1/2.0
			WRITE(*,*)
		ENDIF
		
	END FUNCTION de_chisq_bao_ver1

  !----------------------------------------------------------------
  ! chisq function for new bao, refer to arXiv:1312.4877 and                           !bossli
  ! 1401.0358(6dFGS + BOSS DR7 + BOSS DR11 + Imp. WiggleZ)
	DOUBLE PRECISION FUNCTION de_chisq_bao_ver3()
		DOUBLE PRECISION :: baodiff(8), covinv(8,8)

		baodiff(1) = de_rstodv(0.106d0) - 0.336d0
		baodiff(2) = 1.0d0 / de_rstodv(0.35d0) - 8.88d0
		baodiff(3) = 1.0d0 / de_rstodv(0.32d0) - 8.25d0   !BOSS DR11 DV(0.32)/rd 
		baodiff(4) = de_DAtord(0.57d0) - 1421d0           !DA(0.57)*(rd_fid/rd)
		baodiff(5) = de_Hrd(0.57d0) - 96.8d0              !H(0.57)*(rd/rd_fid).
		baodiff(6) = de_DVtors(0.44d0) - 1716.4d0         !
		baodiff(7) = de_DVtors(0.60d0) - 2220.8d0         !Imp. WiggleZ DV(z)*(rs_fid/rs)
		baodiff(8) = de_DVtors(0.73d0) - 2516.1d0         !

		covinv = 0.0d0
		covinv(1,1) = 4444.4d0

		covinv(2,2) = 34.602d0

		covinv(3,3) = 39.0625d0
		covinv(4,4) = 0.003523712451964d0
		covinv(5,5) = 0.121927762078074d0

		covinv(6,6) = 2.17898878d-4
		covinv(7,7) = 1.70712004d-4
		covinv(8,8) = 1.65283175d-4

		covinv(4,5) = -0.011172240969447d0
		covinv(5,4) = covinv(4,5)

		covinv(6,7) = -11163.3221d0
		covinv(7,6) = covinv(6,7)

		covinv(6,8) = 4.6982851d-4
		covinv(8,6) = covinv(6,8)
		
		covinv(7,8) = -7.1847155d-4
		covinv(8,7) = covinv(7,8)

		de_chisq_bao_ver3 = dot_product(baodiff,matmul(covinv,baodiff))
		IF(pr_chisq_info) THEN
			WRITE(*,*) "        de_rstodv(0.10) = ", de_rstodv(0.1d0)		
			WRITE(*,*) "        de_rstodv(0.35) = ", de_rstodv(0.35d0)
			WRITE(*,*) "        de_rstodv(0.32) = ", de_rstodv(0.32d0)			
			WRITE(*,*) "        de_DAtord(0.57d0) = ", de_DAtord(0.57d0)			
			WRITE(*,*) "        de_Hrd(0.57d0) = ", de_Hrd(0.57d0)
			WRITE(*,*) "        de_DVtors(0.44d0) = ", de_DVtors(0.44d0)
			WRITE(*,*) "        de_DVtors(0.60d0) = ", de_DVtors(0.60d0)
			WRITE(*,*) "        de_DVtors(0.73d0) = ", de_DVtors(0.73d0)			
			WRITE(*,*) "     chisq_bao_ver3 (lnlike) = ", de_chisq_bao_ver3, de_chisq_bao_ver3/2.0
			WRITE(*,*)
		ENDIF
	END FUNCTION de_chisq_bao_ver3
  ! BOSS DR11 anisotropy
  !  DV = (1264 ± 25 Mpc)(rd /rd,fid ) at z = 0.32 (isotropic constraint at z=0.32).
  !  DA = (1421 ± 20 Mpc)(rd /rd,fid ) and H = (96.8 ± 3.4 km/s/Mpc)(rd,fid /rd) at z = 0.57, 
  !    with a correlation coefficient between DA and H of 0.539
  !  we shall use rd,fid = 153.19 for EH98 (using 149.28 for camb)
	DOUBLE PRECISION FUNCTION de_chisq_dr11()
		DOUBLE PRECISION :: covinv(2,2), baodiff(2), rd_fid = 153.19d0
		baodiff(1) = de_DAtord(0.57d0)*rd_fid - 1421.0d0           !DA(0.57)*(rd_fid/rd)
		baodiff(2) = de_Hrd(0.57d0)*de_const_c/rd_fid - 96.8d0              !H(0.57)*(rd/rd_fid).
		covinv(1,1) = 0.003523712451964d0
		covinv(2,2) = 0.121927762078074d0
		covinv(1,2) = -0.011172240969447d0
		covinv(2,1) = covinv(1,2)
		de_chisq_dr11 = dot_product(baodiff,matmul(covinv,baodiff)) &
			+ ((1.0d0/de_rstodv(0.32d0)-8.25d0)/0.16d0)**2.0 ! in fact 1/rstodv is dv/rd
		if(pr_chisq_info) print *, "     chisq_dr11 (lnlike) = ", de_chisq_dr11, de_chisq_dr11/2.0
	end function de_chisq_dr11
  ! Improved WiggleZ 
  !  DV * (rs_fiducial / rs) = 1716±83 Mpc, 2221±101 Mpc, 2516±86 Mpc (68% CL) at z = 0.44, 0.6, 0.73.
  !  we shall use rs,fid = 152.3 rather than 148.6
	DOUBLE PRECISION FUNCTION de_chisq_impwig()
		DOUBLE PRECISION :: covinv(3,3), baodiff(3), rs_fid = 152.3d0
		baodiff(1) = de_DVtors(0.44d0)*rs_fid - 1716.4d0         !
		baodiff(2) = de_DVtors(0.60d0)*rs_fid - 2220.8d0         !Imp. WiggleZ DV(z)*(rs_fid/rs)
		baodiff(3) = de_DVtors(0.73d0)*rs_fid - 2516.1d0         !            
		! cov, diag
		covinv(1,1) = 2.17898878d-4
		covinv(2,2) = 1.70712004d-4
		covinv(3,3) = 1.65283175d-4
		! cov, nodiag: 1,2
		covinv(1,2) = -1.11633221d-4
		covinv(2,1) = covinv(1,2)
		! cov, nodiag: 1,3
		covinv(1,3) = 0.46982851d-4
		covinv(3,1) = covinv(1,3)
		! cov, nodiag: 2,3
		covinv(2,3) = -0.71847155d-4
		covinv(3,2) = covinv(2,3)
		de_chisq_impwig = dot_product(baodiff,matmul(covinv,baodiff)) 
		if(pr_chisq_info) print *, "     chisq_impwig (lnlike) = ", de_chisq_impwig, de_chisq_impwig/2.0
	end function de_chisq_impwig

!#####################################################
!#####################################################
  ! H chisq
!#####################################################
!#####################################################
  !----------------------------------------------------------------
  ! chisq FUNCTION for H
  !----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION de_chisq_h_Riess()
		de_chisq_h_Riess    = ((de_CP%h-0.738d0)/0.024d0)**2.0
		IF(pr_chisq_info) THEN
			WRITE(*,*) "     chisq_h_Riess (lnlike) = ", de_chisq_h_Riess, de_chisq_h_Riess/2.0
		ENDIF
	END FUNCTION de_chisq_h_Riess
	

  !----------------------------------------------------------------
  ! chisq FUNCTION for H
  !----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION de_chisq_h_Carnegie()
		de_chisq_h_Carnegie    = ((de_CP%h-0.743d0)/0.021d0)**2.0
		IF(pr_chisq_info) THEN
			WRITE(*,*) "     chisq_h_Carnegie (lnlike) = ", de_chisq_h_Carnegie, de_chisq_h_Carnegie/2.0
		ENDIF
	END FUNCTION de_chisq_h_Carnegie
END MODULE de_chisqs
