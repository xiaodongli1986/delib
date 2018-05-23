

!-----------------------------------------------------------------------
!------------------ BEGINNING OF THE module Union2 ---------------------
!-----------------------------------------------------------------------
module de_chisqs

use de_model_init
use de_chisqs_JLA

IMPLICIT NONE
        
        private
        public :: de_chisq_all
        
        !supernovae chisq functions
        public :: de_chisq_snls3, de_chisq_union2p1, de_chisq_union2p1zc, de_chisq_union2p1zcs, &
                de_chisq_g06, de_chisq_jla, de_chisq_bao_ver3, de_chisq_bao_desi, de_chisq_ap_desi_wcdm, de_chisq_ap_desi_cpl, &
                de_chisq_bao_plc_desi_wcdm, de_chisq_bao_plc_desi_cpl
        
        !bao chisq functions
        public :: de_chisq_6dFGS!, de_chisq_2dFGS
        public :: de_chisq_SDSSDR7_MGS
        public :: de_chisq_SDSSDR7_old, de_chisq_SDSSDR7_new, de_chisq_BOSSDR9, de_chisq_dr11
        public :: de_chisq_WigZ_rsDV, de_chisq_WigZ_A, de_chisq_impwig
        public :: de_chisq_bao_ver1, de_chisq_bao_ver2
        public :: de_chisq_sdssdr12_tomobao_zhao,  de_chisq_sdssdr12_tomobao_wang, de_chisq_sdssdr12_wang, de_chisq_SDSSDR14_QSO, de_chisq_SDSSDR12_Alam
        public :: de_chisq_sdssdr12_APinlcdm

        !CMB chisq functions
        public :: de_chisq_wmap7, de_chisq_wmap9, de_chisq_planck, de_chisq_planck3yr_flat
        
        !H0 chisq functions
        public :: de_chisq_h_Riess, de_chisq_h_Carnegie


!Some useful constants
        double precision,   parameter  ::  logzero = 1.0d10
        double precision,   parameter  ::  zfacsq = 25.0d0/(LOG(10.0d0))**2

!Some snls parameters
        double precision,   parameter  ::  de_alphatol = 1.0d-10, betatol = 1.0d-10
        double precision,   parameter  ::  scriptmcut = 10.0d0
        double precision,   parameter  ::  intrinsicsq(4) = (/ 0.0675d0**2, 0.1133d0**2, 0.0815d0**2, 0.0989d0**2 /)
        double precision,   parameter  ::  pecz = 0.0005d0
        character,          parameter  ::  uplo = 'U' !For LAPACK

!Supernova data type
        type, private :: supernova
                character(LEN=20) :: name  !The name of the SN
                double precision  :: zhel, zcmb    !The heliocentric and CMB frame redshifts
                double precision  :: z_var         !The variance of the redshift
                double precision  :: mag           !The K-corrected peak magnitude
                double precision  :: mag_var       !The variance of mag
                double precision  :: stretch       !The light-curve fit stretch parameter
                double precision  :: stretch_var   !The variance in the stretch
                double precision  :: colour        !The colour of the SN
                double precision  :: colour_var    !The variance of colour
                double precision  :: thirdvar      !Third variable for scripm split
                double precision  :: thirdvar_var  !Variance in thirdvar
                double precision  :: cov_mag_stretch !Covariance between mag and stretch
                double precision  :: cov_mag_colour  !Covariance between mag and colour
                double precision  :: cov_stretch_colour !Covariance between stretch and colour
                integer           :: dataset       !Subset identifier if subset dependent intrinsic disp is used
        end type supernova

!---------------------------------------------------
!!! SNLS3
!Number of sn
        integer, parameter   ::   nsnls3  =  472

!Supernova data
        type( supernova )    ::   sndata(nsnls3)

!The covarariance matrix 
        double precision, private     ::   m_covmat(nsnls3,nsnls3), s_covmat(nsnls3,nsnls3)
        double precision, private     ::   c_covmat(nsnls3,nsnls3), m_s_covmat(nsnls3,nsnls3)
        double precision, private     ::   m_c_covmat(nsnls3,nsnls3), s_c_covmat(nsnls3,nsnls3)

!The inverse covariance matrix
        double precision, private     ::   invcovmat(nsnls3,nsnls3)

!Some useful arrays
        double precision, private     ::   pre_vars(nsnls3), A1(nsnls3), A2(nsnls3)
        double precision, private     ::   snls3_lumdists(nsnls3), diffmag(nsnls3), invvars(nsnls3)
  
!Some useful variables
        logical, public     :: snls_read       =  .FALSE.
        logical, public     :: snls_prepped    =  .FALSE.
        logical                    :: snls_cov_check  =  .FALSE.
        logical, private    :: first_inversion =  .TRUE.
        logical, private    :: use_ab_prev     =  .TRUE.
        double precision    :: alpha_prev, beta_prev

!---------------------------------------------------
!!! Gold06
        integer, private, parameter :: g06num = 292
        double precision :: g06z(g06num), g06mu(g06num), g06dmu(g06num)
        integer :: g06mark(g06num) ! 1 means gold; 0 means silver
        integer :: g06goodbad(g06num) ! 1 means good; 0 means bad
        logical, public :: g06_read = .false.                
                                


!---------------------------------------------------
!!! Union2.1
!Union2.1 data
        integer, private, parameter           :: union2p1_num  =  580
        double precision, private         :: union2p1_z(union2p1_num), union2p1_moduli(union2p1_num), &
                union2p1_modulierr(union2p1_num), union2p1_plow(union2p1_num), union2p1_sumninv
        double precision, private         :: union2p1_Ninv(union2p1_num,union2p1_Num)

!Union2 data
        integer, private, parameter           :: union2_num  =  558
        double precision, private         :: union2_z(union2_num), union2_moduli(union2_num), &
                union2_modulierr(union2_num), union2_plow(union2_num), union2_sumninv
        double precision, private         :: union2_Ninv(union2_num,union2_Num)

! BAO related
        double precision, private :: sdssdr12_tomobao_zhao_z(18), &
                sdssdr12_tomobao_zhao_obs(18), sdssdr12_tomobao_zhao_invcov(18,18)
        logical, private :: sdssdr12_tomobao_zhao_inited = .FALSE.
        double precision, private :: sdssdr12_tomobao_wang_z(18), &
                sdssdr12_tomobao_wang_obs(18), sdssdr12_tomobao_wang_invcov(18,18)
        logical, private :: sdssdr12_tomobao_wang_inited = .FALSE.
        double precision :: APlcdm_oms(999), APlcdm_chisqs(999)
        logical :: APlcdm_readin = .false.

!Some useful variables
        logical, public                    :: union2p1_syscovmat  =  .TRUE.
        logical, public             :: union2p1_inited     =  .FALSE.
        double precision            :: union2p1_lumdists(union2p1_num), union2p1_diffs(union2p1_num)
        integer, private            :: count_chisq     =   0

        logical, public                    :: union2_syscovmat  =  .TRUE.
        logical, public             :: union2_inited     =  .FALSE.
        double precision            :: union2_lumdists(union2p1_num), union2_diffs(union2_num)



  CONTAINS


  !------------------------------------------
  ! The total chisq. 
  !------------------------------------------
        double precision function de_chisq_all()
        
                if(de_CP%Odm0<0.0d0 .or. de_CP%Odm0>1.0d0 &
                        .or. de_CP%h<0.5d0 .or. de_CP%h>1.0d0) then
                        de_chisq_all = logzero
                        RETURN
                endif
                
                CALL de_Init()

                de_chisq_all =  de_chisq_union2p1() + de_chisq_snls3() &
                         + de_chisq_wmap7() + de_chisq_h_Carnegie() + de_chisq_h_Riess() &
                        + de_chisq_bao_ver1() + de_chisq_bao_ver2()
                !+ de_chisq_union2p1() + chisq_bao_old()
                If(de_chisq_all .ge. logzero) de_chisq_all = logzero
                count_chisq = count_chisq + 1
        end function de_chisq_all

!#####################################################
!#####################################################
  ! Supernovae chisq
!#####################################################
!#####################################################

  ! Gold06
  ! Gold06 read_in
        SUBROUTINE g06_init()
                integer :: i, i_gold, i_silver
                character(LEN=100) :: tmpstr, tmpstr1
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
        end SUBROUTINE g06_init
  ! Gold06 Chisq
        double precision function de_chisq_g06(ignoresilver, nobadsilver)
                ! Dummy 
                logical, INTENT(IN) :: ignoresilver 
                logical, INTENT(IN), optional :: nobadsilver
                ! Local
                integer :: i, num_bin
                double precision :: z,mu,dmu,lumdist,muth,diffmu,dmufac, A, B, C

                if(g06_read.eqv..FALSE.) then
                        CALL g06_init()
                endif

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
                endDO

                de_chisq_g06 = A - B*B / C
        end function de_chisq_g06


  ! SNLS3
  ! SNLS3 read_in
          !---------------------------------------------------------------
          ! READ in the SNLS data. Including 1 data FILE, 6 covmats.     
          ! The global variable snls_READ will be setted to .TRUE.
          !---------------------------------------------------------------
        SUBROUTINE read_snls_dataset   
                integer  :: I
                double precision :: dz, dm, ds, dc, dt

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
                endDO
                CLOSE(70); CLOSE(71); CLOSE(72); CLOSE(73); CLOSE(74); CLOSE(75); CLOSE(76);

                snls_read = .TRUE.
                first_inversion = .TRUE.
                snls_prepped = .FALSE.
                write(*,*) "snls read_in complete"   
        end SUBROUTINE read_snls_dataset
  ! SNLS3 preparation
          !----------------------------------------------------------------
          ! Prepares the data for fitting by pre-calculating the parts of  
          !  the errors that can be done ahead of time.                   
          ! READ_snls_dataset() must have been CALLed before CALLing this       
          ! The global variable snls_prepped will be setted to .TRUE.
          !----------------------------------------------------------------
        SUBROUTINE snls_prep()
                integer :: i

                pre_vars = sndata%mag_var + intrinsicsq(sndata%dataset+1)
                DO I=i,nsnls3
                pre_vars(i) = pre_vars(i) + &
                        zfacsq * ( sndata(i)%z_var + pecz**2 ) * &
                          ( (1.0d0 + sndata(i)%zcmb)/&
                        (sndata(i)%zcmb*(1.0d0 + 0.5d0*sndata(i)%zcmb)) )**2
                endDO

                DO i=1, nsnls3
                        if (sndata(i)%thirdvar .LE. scriptmcut ) then
                        A1(i) = 1.0d0;  A2(i) = 0.0d0
                        ELSE
                        A1(i) = 0.0d0;  A2(i) = 1.0d0
                        endif
                endDO
!                write(*,*) "test pre_vars 2: ", pre_vars(1), pre_vars(100), pre_vars(300)
                write(*,*) "snls preparation complete"   
                snls_prepped = .TRUE.
                first_inversion = .TRUE.
        end SUBROUTINE snls_prep 
  ! SNLS3 inv covmat
        SUBROUTINE inv_cov_mat( alpha, beta, STATUS )     
                double precision  ::  alpha, beta
                integer :: STATUS, i, j
                double precision  ::  alphasq, betasq, alphabeta 

!                if (use_ab_prev .EQ. .TRUE. .AND. .NOT. first_inversion .AND. (ABS(alpha-alpha_prev) .LT. alphatol) .AND. &
!                        ( ABS(beta-beta_prev) .LT. betatol )) then           !Previous invcovmatrix is close enough
!                        STATUS = 0
!                        RETURN
!                endif

                alphasq = alpha * alpha
                betasq = beta * beta
                alphabeta = alpha * beta

!                Build the covariance matrix    
                 invcovmat = m_covmat &
                        + alphasq * s_covmat &
                        + betasq * c_covmat &
                        + 2.0d0 * alpha * m_s_covmat &  
                        - 2.0d0 * beta * m_c_covmat &
                        - 2.0d0 * alphabeta * s_c_covmat

                if(snls_cov_check) then
                        write(*,*) "SNLS covmat check:"
                        write(*,*) m_covmat(1,1), m_covmat(235,256), m_covmat(nsnls3, nsnls3)                        
                        write(*,*) s_covmat(1,1), s_covmat(235,256), s_covmat(nsnls3, nsnls3)
                        write(*,*) c_covmat(1,1), c_covmat(235,256), c_covmat(nsnls3, nsnls3)
                        write(*,*) m_s_covmat(1,1), m_s_covmat(235,256), m_s_covmat(nsnls3, nsnls3)        
                        write(*,*) m_c_covmat(1,1), m_c_covmat(235,256), m_c_covmat(nsnls3, nsnls3)        
                        write(*,*) s_c_covmat(1,1), s_c_covmat(235,256), s_c_covmat(nsnls3, nsnls3)        
                        write(*,*) invcovmat(1,1), invcovmat(235,256), invcovmat(nsnls3, nsnls3)
                endif
                
!                Update the diagonal terms
                DO I=1, nsnls3
                        invcovmat(I,I) = invcovmat(I,I) + pre_vars(I) &
                                + alphasq * sndata(I)%stretch_var &
                                + betasq  * sndata(I)%colour_var &
                                + 2.0d0 * alpha * sndata(I)%cov_mag_stretch &
                                - 2.0d0 * beta * sndata(I)%cov_mag_colour &
                                - 2.0d0 * alphabeta * sndata(I)%cov_stretch_colour
                endDO
                
                print *, 'EXIT!'; stop
                !CALL DPOTRF(uplo,nsnls3,invcovmat,nsnls3,STATUS)
                print *, 'EXIT!'; stop
                !CALL DPOTRI(uplo,nsnls3,invcovmat,nsnls3,STATUS)
        
                first_inversion = .FALSE.
                alpha_prev = alpha
                beta_prev  = beta
        end SUBROUTINE inv_cov_mat
  ! SNLS3 chisq
        double precision function de_chisq_snls3()
                double precision :: alpha, alphasq, beta, betasq, alphabeta
                double precision :: ogamma, Omegar, Omegade1, Omegade2
                double precision :: now_broader
                double precision :: estimated_scriptm, wtval
                double precision :: zhel, zcmb, chisq
                double precision :: amarg_A, amarg_B, amarg_C, amarg_D, amarg_E, amarg_F, tempG
                integer :: i, num_bin, STATUS

                if(snls_read .EQV. .FALSE.)        CALL read_snls_dataset()
                if(snls_prepped .EQV. .FALSE. )  CALL snls_prep()

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
                endDO


!                Calculate estimated_scriptm and diffmag
                invvars = 1.0d0 / ( pre_vars + alphasq * sndata%stretch_var &
                        + betasq * sndata%colour_var &
                        + 2.0d0 * alpha * sndata%cov_mag_stretch &
                        - 2.0d0 * beta * sndata%cov_mag_colour &
                        - 2.0d0 * alphabeta * sndata%cov_stretch_colour )
!                write(*,*) alpha, beta, alphasq, betasq, alphabeta

                wtval = SUM( invvars )
                estimated_scriptm= SUM( (sndata%mag - snls3_lumdists)*invvars ) / wtval
                diffmag = sndata%mag - snls3_lumdists + alpha*( sndata%stretch - 1.0d0 ) &
                        - beta * sndata%colour - estimated_scriptm 

                CALL inv_cov_mat(alpha,beta,STATUS)

!                Now find the amarg_ parameters
!                We re-use the invvars variable to hold the intermediate product
!                which is sort of naughty
!                invvars = V^-1 * diffmag (invvars = 1.0*invcovmat*diffmag+0*invvars)
                print *, 'EXIT!'; stop
                !CALL DSYMV(uplo,nsnls3,1.0d0,invcovmat,nsnls3,diffmag,1,0.0d0,invvars,1)

                amarg_A = DOT_PRODUCT( diffmag, invvars ) ! diffmag*V^-1*diffmag
                amarg_B = DOT_PRODUCT( invvars, A1 ) !diffmag*V^-1*A1
                amarg_C = DOT_PRODUCT( invvars, A2 ) !diffmag*V^-1*A2

!                Be naughty again and stick V^-1 * A1 in invvars
                print *, 'EXIT!'; stop
                !CALL DSYMV(uplo,nsnls3,1.0d0,invcovmat,nsnls3,A1,1,0.0d0,invvars,1)
                amarg_D = DOT_PRODUCT( invvars, A2 ) !A2*V^-1*A1
                amarg_E = DOT_PRODUCT( invvars, A1 ) !A1*V^-1*A1

!                now V^-1 * A2
                                print *, 'EXIT!'; stop
                !CALL DSYMV(uplo,nsnls3,1.0d0,invcovmat,nsnls3,A2,1,0.0d0,invvars,1)
                amarg_F = DOT_PRODUCT( invvars, A2 ) !A2*V^-1*A2
                tempG = amarg_F - amarg_D*amarg_D/amarg_E; !g/e

!                Marginalized chisq
                chisq = amarg_A + LOG( amarg_E*de_inv_twopi ) + &
                        LOG( tempG * de_inv_twopi ) - amarg_C*amarg_C/tempG - &
                        amarg_B*amarg_B*amarg_F / ( amarg_E*tempG ) + &
                        2.0d0*amarg_B*amarg_C*amarg_D/(amarg_E*tempG )

                de_chisq_snls3 = chisq 
                if(pr_chisq_info) then
                        write(*,*) "     chisq_snls3 (lnlike) = ", de_chisq_snls3, de_chisq_snls3/2.0
                        write(*,*)
                endif
        end function de_chisq_snls3
  ! Union2.1
  ! Union2.1 Read in
        SUBROUTINE union2p1_init()
                character(LEN=20) :: union2p1_name
                integer :: i

                OPEN(UNIT=101, FILE=trim(adjustl(de_data_path))//'sn_z_mu_dmu_plow_union2.1.txt')
                DO i = 1, union2p1_num
                        READ(101,*) union2p1_name, union2p1_z(i), union2p1_moduli(i), union2p1_modulierr(i), union2p1_plow(i)
                endDO
                CLOSE(101)

                if(union2p1_syscovmat) then
                        OPEN(UNIT=102, FILE=trim(adjustl(de_data_path))//'sn_wmat_sys_union2.1.txt')
                        write(*,*) "We are using sys in SN data."
                        ELSE
                        OPEN(UNIT=102, FILE=trim(adjustl(de_data_path))//'sn_wmat_nosys_union2.1.txt')
                endif

                DO i = 1, union2p1_num
                        READ(102,*) union2p1_ninv(i,1:union2p1_num)
                endDO
                CLOSE(102)
!                union2p1_sumninv = sum(union2p1_ninv)
        end SUBROUTINE union2p1_init
  ! Union2.1 chisq
        double precision function de_chisq_union2p1()
                ! Local
                integer :: i, num_bin
                
                if(union2p1_inited .EQV. .FALSE.) then
                        CALL union2p1_init()
                        write(*,*) "Union2.1 read in completed"
                        union2p1_inited = .TRUE.
                endif

                DO i = 1, union2p1_num
                        num_bin = Ceiling(union2p1_z(i)*16.0d0)
                        union2p1_lumdists(i) = de_Simpson(de_inv_e, 0.0d0, union2p1_z(i), num_bin) 
                        union2p1_lumdists(i) = de_fk(union2p1_lumdists(i))/ de_CP%H0
                        union2p1_lumdists(i) = 5.0 * LOG10((1.0d0+union2p1_z(i))*union2p1_lumdists(i)) + 25        
                        union2p1_diffs(i) = union2p1_lumdists(i) - union2p1_moduli(i)
                endDO

!                write(*,'(<560>(f10.5,1x))') union2p1_diffs
                de_chisq_union2p1 = dot_product(union2p1_diffs,matmul(union2p1_ninv,union2p1_diffs)) 

                if(pr_chisq_info) then
                        write(*,*) "     chisq_Union2.1 = ", de_chisq_union2p1
                        write(*,*)
                endif
        end function de_chisq_union2p1


        double precision function de_chisq_union2p1zcs(zcut_max, SNIaM)
                ! Dummy
                double precision, INTENT(IN) :: zcut_max, SNIaM
                ! Local
                integer :: i, num_bin
                
                if(union2p1_inited .EQV. .FALSE.) then
                        CALL union2p1_init()
                        write(*,*) "Union2.1 read in completed"
                        union2p1_inited = .TRUE.
                endif

                DO i = 1, union2p1_num
                        num_bin = Ceiling(union2p1_z(i)*16.0d0)
                        union2p1_lumdists(i) = de_Simpson(de_inv_e, 0.0d0, union2p1_z(i), num_bin) 
                        union2p1_lumdists(i) = de_fk(union2p1_lumdists(i))
                        union2p1_lumdists(i) = 5.0 * LOG10((1.0d0+union2p1_z(i))*union2p1_lumdists(i)) + 25 + SNIaM
                        union2p1_diffs(i) = union2p1_lumdists(i) - union2p1_moduli(i)
                        if(union2p1_z(i) > zcut_max) then
                                union2p1_diffs(i) = 0.0
                        endif
                endDO
                de_chisq_union2p1zcs = dot_product(union2p1_diffs,matmul(union2p1_ninv,union2p1_diffs)) 

        end function de_chisq_union2p1zcs
        ! chisq of Union2.1, with maximal cut of redshift
        ! Actually it is very demanding; I do a search in H direction
        double precision function de_chisq_union2p1zc(zcut_max)
                ! Dummy
                double precision, INTENT(IN) :: zcut_max
                ! Local
                double precision :: H1, H2, Hmid, chisqA, chisqB, chisqmid, deltaH = 0.0001

                H1 = -100; H2 = 100.0;
                do while(abs(H2-H1).ge.3.0*deltaH)
                        Hmid = (H1+H2)/2.0
                        chisqA = de_chisq_union2p1zcs(zcut_max, Hmid)
                        chisqB = de_chisq_union2p1zcs(zcut_max, Hmid+max((H2-H1)*0.01,deltaH))
!                        print *, real(H1), real(H2), real(chisqA), real(chisqB)
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
        double precision function de_chisq_wmap7()
                double precision :: zstar, R, lA
                double precision :: g1, g2, DAzstar, rszstar
                double precision :: zstar_ML, R_ML, lA_ML, DVecCMB(3), CovCMB(3,3)

                CovCMB(1,1)=2.305d0;   CovCMB(1,2)=29.698d0;     CovCMB(1,3)=-1.333d0;
                CovCMB(2,1)=29.698d0;  CovCMB(2,2)=6825.27d0;    CovCMB(2,3)=-113.180d0;
                CovCMB(3,1)=-1.333d0;  CovCMB(3,2)=-113.180d0;   CovCMB(3,3)=3.414d0;
                lA_ML=302.09d0; R_ML=1.725d0; zstar_ML=1091.3d0;
   
                DVecCMB   = (/de_CP%lA-lA_ML, de_CP%R-R_ML, de_CP%zstar-zstar_ML/)
                de_chisq_wmap7 = DOT_PRODUCT(DVecCMB,matmul(CovCMB,DVecCMB))
                
                if(pr_chisq_info) then
                        write(*,*) "     chisq_wmap7 (lnlike) = ", de_chisq_wmap7, de_chisq_wmap7/2.0
                        write(*,*)
                endif
        end function de_chisq_wmap7
  ! planck 1yr
        double precision function de_chisq_planck()   !de_chisq_planck   !de_chisq_wmap7()
                double precision :: lA, R, Ob0hsq    !zstar
                double precision :: g1, g2, DAzstar, rszstar
                double precision :: lA_ML, R_ML, Ob0hsq_ML, DVecCMB(3), CovCMB(3,3)   !zstar_ML

                lA     = de_CP%lA
                R      = de_CP%R 
                Ob0hsq = de_CP%Ob0hsq

                CovCMB(1,1)=43.0179647667078;   CovCMB(1,2)=-366.771828690038;     CovCMB(1,3)=2972.52785899757;   !-----------------
                CovCMB(2,1)=-366.771828690038;  CovCMB(2,2)=24872.6578215152;      CovCMB(2,3)=446498.465930772;   ! planck's covmat
                CovCMB(3,1)=2972.52785899758;   CovCMB(3,2)=446498.465930772;      CovCMB(3,3)=21554701.3253821;   ! and mean value
                lA_ML=301.57d0; R_ML=1.7407d0; Ob0hsq_ML=0.02228d0                                                 !----------------- 

                DVecCMB         = (/lA-lA_ML, R-R_ML, Ob0hsq-Ob0hsq_ML/)
                !DVecCMB         = (/lA-lA_ML, 0.0d0, 0.0d0 /) !R-R_ML, Ob0hsq-Ob0hsq_ML/)
                de_chisq_planck = DOT_PRODUCT(DVecCMB,matmul(CovCMB,DVecCMB))

                if(pr_chisq_info) then
                        write(*,*) "     chisq_planck (lnlike) = ", de_chisq_planck, de_chisq_planck/2.0
                        write(*,*)
                endif

                !write(*,*) "la,r,ob0hsq", lA, R, Ob0hsq    !bossli
        end function de_chisq_planck   !de_chisq_planck   !de_chisq_wmap7
  ! planck 3-yr
        double precision function de_chisq_planck3yr_flat()   ! using 1509.02198; flat prior
                double precision :: lA, R, Ob0hsq, ns, siglA, sigR, sigOb0hsq, signs    !zstar
                double precision :: g1, g2, DAzstar, rszstar
                double precision :: lA_ML, R_ML, Ob0hsq_ML, ns_ML, p(4), cov(4,4),invcov(4,4), sigs(4)   !zstar_ML
                integer :: i, j


                siglA = 0.090; sigR = 0.0048; sigOb0hsq = 0.00016; signs = 0.0048
                sigs = (/ siglA, sigR, sigOb0hsq, signs /)
                cov(1,1) = 1.0; cov(1,2) = 0.3996; cov(1,3) = -0.3181; cov(1,4) = -0.3004;
                                cov(2,2) = 1.0;    cov(2,3) = -0.6891; cov(2,4) = -0.7677;
                                                   cov(3,3) = 1.0;     cov(3,4) = 0.5152;
                                                                       cov(4,4) = 1.0
                cov(2,1) = cov(1,2); 
                cov(3,1) = cov(1,3); cov(3,2) = cov(2,3);
                cov(4,1) = cov(1,4); cov(4,2) = cov(2,4); cov(4,3) = cov(3,4);
                do i = 1,4
                  do j = 1,4
                    cov(i,j) = cov(i,j) * sigs(i) * sigs(j)
                  enddo
                enddo
                call de_nizhen(cov,invcov,4)

                
                lA_ML=301.77d0; R_ML=1.7482d0; Ob0hsq_ML=0.02226d0; ns_ML = 0.9653d0;                                                 !-----------------

                lA     = de_CP%lA
                R      = de_CP%R 
                Ob0hsq = de_CP%Ob0hsq
                ns     = de_CP%ns

                p  = (/lA-lA_ML, R-R_ML, Ob0hsq-Ob0hsq_ML, ns-ns_ML/)
                de_chisq_planck3yr_flat = DOT_PRODUCT(p,matmul(invcov,p))

                if(pr_chisq_info) then
                        write(*,*) "     chisq_planck (lnlike) = ", de_chisq_planck3yr_flat, de_chisq_planck3yr_flat/2.0
                        write(*,'(A,f8.2,f9.2,f8.5,f6.3)') "       la,r,ob0hsq,ns      = ", real(lA), real(R), real(Ob0hsq), real(ns)    !bossli
                        write(*,'(A,f8.2,f9.2,f8.5,f6.3)') "       la,r,ob0hsq,ns (ML) = ", real(lA_ML), real(R_ML), real(Ob0hsq_ML), real(ns_ML)    !bossli
                        write(*,*)
                endif

        end function de_chisq_planck3yr_flat
  ! WMAP9
        double precision function de_chisq_wmap9()   !de_chisq_planck   !de_chisq_wmap7()
                double precision :: lA, R, Ob0hsq    !zstar
                double precision :: g1, g2, DAzstar, rszstar
                double precision :: lA_ML, R_ML, Ob0hsq_ML, DVecCMB(3), CovCMB(3,3)   !zstar_ML

                lA     = de_CP%lA
                R      = de_CP%R 
                Ob0hsq = de_CP%Ob0hsq

                CovCMB(1,1)=3.68712321402858;   CovCMB(1,2)=-14.1725878878498;     CovCMB(1,3)=2566.01640514695;   !-----------------
                CovCMB(2,1)=-14.1725878878498;  CovCMB(2,2)=5179.04970797377;      CovCMB(2,3)=73212.4449017310;   ! wmap9's covmat
                CovCMB(3,1)=2566.01640514695;   CovCMB(3,2)=73212.4449017310;      CovCMB(3,3)=6692540.05113778;   ! and mean value
                lA_ML=302.02d0; R_ML=1.7327d0; Ob0hsq_ML=0.02260d0                                                 !----------------- 

                DVecCMB         = (/lA-lA_ML, R-R_ML, Ob0hsq-Ob0hsq_ML/)
                de_chisq_wmap9 = DOT_PRODUCT(DVecCMB,matmul(CovCMB,DVecCMB))

                if(pr_chisq_info) then
                        write(*,*) "     chisq_wmap9 (lnlike) = ", de_chisq_wmap9, de_chisq_wmap9/2.0
                        write(*,*)
                endif

                !write(*,*) "la,r,ob0hsq", lA, R, Ob0hsq    !bossli
        end function de_chisq_wmap9   !de_chisq_planck   !de_chisq_wmap7

!#####################################################
!#####################################################
  ! BAO chisq
!#####################################################
!#####################################################
  ! 6dFGS rstodv(0.106) = 0.336 +/- 0.015 arXiv:1106.3366
        double precision function de_chisq_6dFGS()
                de_chisq_6dFGS = ((de_rstodv(0.106d0) - 0.336d0)/0.015d0)**2.0 
                if(pr_chisq_info) then
                        write(*,*) "        de_rstodv(0.106) = ", &
                                de_rstodv(0.106d0)
                        write(*,*) "     chisq_6dFGS (lnlike) = ", &
                                de_chisq_6dFGS, de_chisq_6dFGS/2.0
                        write(*,*)
                endif
        end function de_chisq_6dFGS
  ! SDSS DR7 MGS dv(0.15) = 664 +/- 25 * (rd/rd_fid) arXiv:1409.3242
        double precision function de_chisq_SDSSDR7_MGS()
                double precision :: rdfid = 148.69d0, x, rdrescale
                rdrescale = 1.027369826 !1.0253413514035492 !
                ! We shall use our own formular to compute rd in their fiducial cosmology, and compare with their rdfid
                ! If our rdfid is 1.027 times larger/smaller then their, ...
                !print *, 'rd= ', de_CP%rszd
                !print *, 'rd/rdfid = ', de_CP%rszd/rdfid
                x = 1.0/(de_rstodv(0.15d0)/rdrescale) * rdfid
                de_chisq_SDSSDR7_MGS = (x-664.0d0)**2.0 / 25.0**2.0 
                if(pr_chisq_info) then
                        write(*,*) "        de_rstodv(0.15) = ", &
                                de_rstodv(0.15d0)
                        write(*,*) "     chisq_SDSSDR7_MGS (lnlike) = ", &
                                de_chisq_SDSSDR7_MGS, de_chisq_SDSSDR7_MGS/2.0
                        write(*,*)
                endif
        end function de_chisq_SDSSDR7_MGS

  ! SDSS DR12 3bin Alam arXiv:1607.03155    APLike
          double precision function de_chisq_SDSSDR12_Alam()
                double precision :: baodiff(6), covinv(6,6)
                double precision :: rd_fid_camb = 147.78d0, rd_fid, rdrescale = 1.02411520606716 !1.0253413514035492 !
                ! We shall use our own formular to compute rd in their fiducial cosmology, and compare with their rdfid
                ! If our rdfid is 1.027 times larger/smaller then their, ...
                !print *, 'rd= ', de_CP%rszd*de_const_c ! derive_rdfid_ratio
                !print *, 'rd/rdfid = ', de_CP%rszd/rdfid*de_const_c ! derive_rdfid_ratio
                rd_fid = rd_fid_camb * rdrescale
!                print *, 'de_DAtord(0.38d0), rd_fid =', de_DAtord(0.38d0), rd_fid
                baodiff(1) = de_DAtord(0.38d0)*(1+0.38d0)*rd_fid - 1518.36d0           !DA(0.38)*(rd_fid/rd)
                baodiff(2) = de_Hrd(0.38d0)*de_const_c/rd_fid - 81.5095d0              !H(0.38)*(rd/rd_fid)
                baodiff(3) = de_DAtord(0.51d0)*(1+0.51d0)*rd_fid - 1977.44d0           !DA(0.51)*(rd_fid/rd)
                baodiff(4) = de_Hrd(0.51d0)*de_const_c/rd_fid - 90.4474d0              !H(0.51)*(rd/rd_fid)
                baodiff(5) = de_DAtord(0.61d0)*(1+0.61d0)*rd_fid - 2283.18d0           !DA(0.61)*(rd_fid/rd)
                baodiff(6) = de_Hrd(0.61d0)*de_const_c/rd_fid - 97.2556d0              !H(0.61)*(rd/rd_fid)
                covinv(1,1) = 0.00278698d0
                covinv(2,1) = -0.00663325d0
                covinv(3,1) = -0.0012451d0
                covinv(4,1) = 0.00347492d0
                covinv(5,1) = 0.000155427d0
                covinv(6,1) = -0.000573549d0
                covinv(1,2) = covinv(2,1)
                covinv(2,2) = 0.378681d0
                covinv(3,2) = 0.00203742d0
                covinv(4,2) = -0.188209d0
                covinv(5,2) = -0.000660168d0
                covinv(6,2) = 0.0186241d0
                covinv(1,3) = covinv(3,1)
                covinv(2,3) = covinv(3,2)
                covinv(3,3) = 0.00258017d0
                covinv(4,3) = -0.00711667d0
                covinv(5,3) = -0.000919586d0
                covinv(6,3) = 0.0032798d0
                covinv(1,4) = covinv(4,1)
                covinv(2,4) = covinv(4,2)
                covinv(3,4) = covinv(4,3)
                covinv(4,4) = 0.48902d0
                covinv(5,4) = 0.00224051d0
                covinv(6,4) = -0.206989d0
                covinv(1,5) = covinv(5,1)
                covinv(2,5) = covinv(5,2)
                covinv(3,5) = covinv(5,3)
                covinv(4,5) = covinv(5,4)
                covinv(5,5) = 0.00141172d0
                covinv(6,5) = -0.00485036d0
                covinv(1,6) = covinv(6,1)
                covinv(2,6) = covinv(6,2)
                covinv(3,6) = covinv(6,3)
                covinv(4,6) = covinv(6,4)
                covinv(5,6) = covinv(6,5)
                covinv(6,6) = 0.342161d0
                de_chisq_SDSSDR12_Alam = dot_product(baodiff,matmul(covinv,baodiff))
                if(pr_chisq_info) then                
                        write(*,*) "        de_DAtord(0.38) = ", de_DAtord(0.38d0)                        
                        write(*,*) "        de_Hrd(0.38) = ", de_Hrd(0.38d0)
                        write(*,*) "        de_DAtord(0.51) = ", de_DAtord(0.51d0)                        
                        write(*,*) "        de_Hrd(0.51) = ", de_Hrd(0.51d0)
                        write(*,*) "        de_DAtord(0.61) = ", de_DAtord(0.61d0)                        
                        write(*,*) "        de_Hrd(0.61) = ", de_Hrd(0.61d0)
                        write(*,*) "        baodiff = ", real(baodiff)
                        write(*,*) "     chisq_SDSSDR12_Alam (lnlike) = ", de_chisq_SDSSDR12_Alam, de_chisq_SDSSDR12_Alam/2.0
                        write(*,*)
                endif
        end function de_chisq_SDSSDR12_Alam

  ! BOSS DR12 tomographic BAO, Zhao et al. 
        double precision function de_chisq_SDSSDR12_tomoBAO_Zhao()
                integer :: i
                double precision :: theory(18), diff(18)
                double precision :: rdrescale = 1.02533817318468d0
                ! fiducial cosmo: {Ω M , Ω b , Ω K , h, σ 8 } = {0.307115, 0.0480, 0, 0.6777, 0.8288}
                !   rszd (from camb) = 147.73
                !   rszd (fitting formular) =    151.473208324573     
                !   rszd (fitting formular) / rszd(from camb) =    1.02533817318468
                if(.not. sdssdr12_tomobao_zhao_inited) then
                        open(unit=9987,file=trim(adjustl(de_data_path))//&
                                '/SDSSDR12_tomoBAO_Zhao/DR12_tomoBAO_Zhao_measurement.txt', action='read')
                        do i = 1, 18
                                read(9987,*) sdssdr12_tomobao_zhao_z(i), sdssdr12_tomobao_zhao_obs(i)
                        enddo
                        close(9987)

                        open(unit=9988,file=trim(adjustl(de_data_path))//&
                                '/SDSSDR12_tomoBAO_Zhao/DR12_tomoBAO_Zhao_invcov.txt', action='read')
                        do i = 1, 18
                                read(9988,*) sdssdr12_tomobao_zhao_invcov(i,1:18)
                        enddo
                        close(9988)
                        write(*,*) 'sdss dr12 tomogarphic bao, Zhao et al. read in complete'
                        if(pr_chisq_info.and..false.) then
                                write(*,*) '     Read in SDSS DR12 tomoBAO Zhao ...'
                                do i = 1, 18
                                        print *, real(sdssdr12_tomobao_zhao_z(i)), real(sdssdr12_tomobao_zhao_obs(i))
                                enddo
                                do i = 1, 18
                                        print *, real(sdssdr12_tomobao_zhao_invcov(i,:))
                                enddo
                        endif
                        sdssdr12_tomobao_zhao_inited=.true.
                endif
                !PRINT *, 'You need to write the definition of chisq!!!'; STOP
                DO i = 1, 9
                        theory(i) = de_DAtord(sdssdr12_tomobao_zhao_z(i))*rdrescale ! Our rd is larger than their rd, we need to rescale 
                        theory(9+i) = de_Hrd(sdssdr12_tomobao_zhao_z(9+i)) * de_const_c / rdrescale
                endDO
                diff = theory - sdssdr12_tomobao_zhao_obs
                  de_chisq_SDSSDR12_tomoBAO_Zhao = dot_product(diff,matmul(sdssdr12_tomobao_zhao_invcov,diff))
                if(pr_chisq_info) then
                        do i = 1, 9
                                print *, 'z, DA/rd theory, obs = ', &
                                        real(sdssdr12_tomobao_zhao_z(i)), real(theory(i)), real(sdssdr12_tomobao_zhao_obs(i))
                        enddo
                        do i = 10, 18
                                print *, 'z, H*rd theory, obs = ', &
                                        real(sdssdr12_tomobao_zhao_z(i)), real(theory(i)), real(sdssdr12_tomobao_zhao_obs(i))
                        enddo
                        write(*,*) "     de_chisq_SDSSDRDR12_tomoBAO_Zhao (lnlike) = ", &
                                de_chisq_SDSSDR12_tomoBAO_Zhao, de_chisq_SDSSDR12_tomoBAO_Zhao/2.0
                endif
        end function de_chisq_SDSSDR12_tomoBAO_Zhao


  ! BOSS DR12 tomographic BAO, Wang et al. 
        double precision function de_chisq_SDSSDR12_tomoBAO_Wang()
                integer :: i
                double precision :: theory(18), diff(18)
                double precision :: rdrescale = 1.02533817318468d0
                ! fiducial cosmo: {Ω M , Ω b , Ω K , h, σ 8 } = {0.307115, 0.0480, 0, 0.6777, 0.8288}
                !   rszd (from camb) = 147.73
                !   rszd (fitting formular) =    151.473208324573     
                !   rszd (fitting formular) / rszd(from camb) =    1.02533817318468
                if(.not. sdssdr12_tomobao_wang_inited) then
                        open(unit=19987,file=trim(adjustl(de_data_path))//&
                                '/SDSSDR12_tomoBAO_Wang/tomographic_BAO.txt', action='read')
                        do i = 1, 18
                                read(19987,*) sdssdr12_tomobao_wang_z(i), sdssdr12_tomobao_wang_obs(i)
                        enddo
                        close(19987)

                        open(unit=19988,file=trim(adjustl(de_data_path))//&
                                '/SDSSDR12_tomoBAO_Wang/tomographic_BAO_invcov.txt', action='read')
                        do i = 1, 18
                                read(19988,*) sdssdr12_tomobao_wang_invcov(i,1:18)
                        enddo
                        close(19988)
                        write(*,*) 'sdss dr12 tomogarphic bao, Wang et al. read in complete'
                        if(pr_chisq_info.and..false. .or. .true.) then
                                write(*,*) '     Read in SDSS DR12 tomoBAO Wang ...'
                                do i = 1, 18
                                        print *, real(sdssdr12_tomobao_wang_z(i)), real(sdssdr12_tomobao_wang_obs(i))
                                enddo
                                do i = 1, 18
                                        print *, real(sdssdr12_tomobao_wang_invcov(i,:))
                                enddo
                        endif
                        sdssdr12_tomobao_wang_inited=.true.
                endif
                !PRINT *, 'You need to write the definition of chisq!!!'; STOP
                DO i = 1, 9
                        theory(2*i-1) = de_DAtord(sdssdr12_tomobao_wang_z(i))*rdrescale ! Our rd is larger than their rd, we need to rescale 
                        theory(2*i) = de_Hrd(sdssdr12_tomobao_wang_z(9+i)) * de_const_c / rdrescale
                endDO
                diff = theory - sdssdr12_tomobao_wang_obs
                  de_chisq_SDSSDR12_tomoBAO_Wang = dot_product(diff,matmul(sdssdr12_tomobao_wang_invcov,diff))
                if(pr_chisq_info) then
                        do i = 1, 9
                                print *, 'z, DA/rd theory, obs = ', &
                                        real(sdssdr12_tomobao_wang_z(i)), real(theory(i)), real(sdssdr12_tomobao_wang_obs(i))
                        enddo
                        do i = 10, 18
                                print *, 'z, H*rd theory, obs = ', &
                                        real(sdssdr12_tomobao_wang_z(i)), real(theory(i)), real(sdssdr12_tomobao_wang_obs(i))
                        enddo
                        write(*,*) "     de_chisq_SDSSDRDR12_tomo_Wang (lnlike) = ", &
                                de_chisq_SDSSDR12_tomoBAO_Wang, de_chisq_SDSSDR12_tomoBAO_Wang/2.0
                endif
        end function de_chisq_SDSSDR12_tomoBAO_Wang
  ! BOSS DR12 9bin Wang arXiv:1607.03154    APLike
          double precision function de_chisq_SDSSDR12_Wang()
                double precision :: baodiff(18), covinv(18,18)
!                double precision :: zsrd = 151.528756316527, rdfid = 147.74d0
                double precision :: rdrescale = 1.02564475644055d0
                integer :: i, j
                baodiff(1) = de_DAtord(0.31d0)*rdrescale - 6.29d0;    baodiff(2) = de_Hrd(0.31d0)*de_const_c/rdrescale - 11550d0;
                baodiff(3) = de_DAtord(0.36d0)*rdrescale - 7.09d0;    baodiff(4) = de_Hrd(0.36d0)*de_const_c/rdrescale - 11810d0;
                baodiff(5) = de_DAtord(0.40d0)*rdrescale - 7.70d0;    baodiff(6) = de_Hrd(0.40d0)*de_const_c/rdrescale - 12120d0;
                baodiff(7) = de_DAtord(0.44d0)*rdrescale - 8.20d0;    baodiff(8) = de_Hrd(0.44d0)*de_const_c/rdrescale - 12530d0;            
                baodiff(9) = de_DAtord(0.48d0)*rdrescale - 8.64d0;    baodiff(10) = de_Hrd(0.48d0)*de_const_c/rdrescale - 12970d0;
                baodiff(11) = de_DAtord(0.52d0)*rdrescale - 8.90d0;   baodiff(12) = de_Hrd(0.52d0)*de_const_c/rdrescale - 13940d0;
                baodiff(13) = de_DAtord(0.56d0)*rdrescale - 9.16d0;   baodiff(14) = de_Hrd(0.56d0)*de_const_c/rdrescale - 13790d0;
                baodiff(15) = de_DAtord(0.59d0)*rdrescale - 9.45d0;   baodiff(16) = de_Hrd(0.59d0)*de_const_c/rdrescale - 14550d0;
                baodiff(17) = de_DAtord(0.64d0)*rdrescale - 9.62d0;   baodiff(18) = de_Hrd(0.64d0)*de_const_c/rdrescale - 14600d0;
                open(234,file=trim(adjustl(de_data_path))//'bao_DR12_Wang_invcov.txt',status='old',action='read')
                  do i=1,18
                    read(234,*) covinv(i,1:18)
                  end do
                close(234)
                de_chisq_SDSSDR12_Wang = dot_product(baodiff,matmul(covinv,baodiff))
                if(pr_chisq_info) then                
                        write(*,*) "     chisq_SDSSDR12_Wang (lnlike) = ", de_chisq_SDSSDR12_Wang, de_chisq_SDSSDR12_Wang/2.0
                        write(*,*)
                endif
        end function de_chisq_SDSSDR12_Wang

  ! BAO like for DESI
        subroutine de_chisq_bao_desi(Rs, DAs, Hs, set_fid, chisqiso, chisqaniso)
                double precision, intent(inout) :: Rs(:), DAs(:), Hs(:)
                double precision, intent(out) :: chisqiso, chisqaniso
                logical, intent(in) :: set_fid
                integer, parameter :: nfile = 10
                character(len=10000) :: nowfiles(nfile), nowfile, tmpstr
                integer :: i,ifile,nline
                double precision :: zs(1000)=0.0, siglnRs(1000)=0.0, siglnDAs(1000)=0.0, &
                        siglnHs(1000)=0.0, chisqisos(1000)=0.0, chisqanisos(1000)=0.0, &
                        Rs_th(1000), DAs_th(1000), Hs_th(1000), &
                        sig1, sig2, det, dR, dDA, dH, coeff, dtheta(2), invcov(2,2)

                ! correlation coefficient between DA & H
                coeff = 0.41

                nowfiles(1) = '/home/xiaodongli/software/delib/data/DESI14kBGS.dat'
                nowfiles(2) = '/home/xiaodongli/software/delib/data/DESI14k.dat'
                nowfiles(3) = '/home/xiaodongli/software/delib/data/DESI14kQSO.dat'

                nline = 0
                chisqisos=0.; chisqanisos=0.; 
                do ifile = 1, 3   
                        nowfile = nowfiles(ifile)
                        !write(*, '(A,A)') ' open file for read: ', trim(adjustl(nowfile))
                        open(unit=132894,file=nowfile)
                        read(132894,*) tmpstr
                        do while(.true.)
                                read(132894,*,end=1909) zs(nline+1), siglnRs(nline+1), siglnDAs(nline+1), siglnHs(nline+1)
                                siglnRs(nline+1)  = siglnRs(nline+1)  /100.
                                siglnDAs(nline+1) = siglnDAs(nline+1) /100.
                                siglnHs(nline+1)  = siglnHs(nline+1)  /100.
                                nline = nline + 1
                                cycle
1909                                exit
                        enddo
                        close(132894); 
                enddo
                if(set_fid) then
                        do i = 1, nline
                                Rs(i)  = de_DV(zs(i)) /de_CP%rszd
                                DAs(i) = de_DA(zs(i)) /de_CP%rszd
                                Hs(i)  = de_CP%H0 / de_inv_e(zs(i)) * de_CP%rszd
                        enddo
                        chisqiso = 0.0; chisqaniso = 0.0
                        return
                else
                        chisqisos = 0.0; chisqanisos = 0.0
                        do i = 1, nline
                                Rs_th(i)  = de_DV(zs(i)) /de_CP%rszd
                                DAs_th(i) = de_DA(zs(i)) /de_CP%rszd
                                Hs_th(i)  = de_CP%H0 / de_inv_e(zs(i)) *de_CP%rszd 
                                ! isotropic bao chisq
                                chisqisos(i) = ( (Rs_th(i) - Rs(i)) / (siglnRs(i)*Rs(i)) ) ** 2.0
                                ! anisotropic bao chisq
                                sig1 = siglnDAs(i) * DAs(i)
                                sig2 = siglnHs(i) * Hs(i)
                                dtheta(1) = DAs_th(i) - DAs(i)
                                dtheta(2) = Hs_th(i) - Hs(i)
                                det = (1-coeff**2.0) * (sig1**2.0 * sig2**2.0)
                                invcov(1,1) = sig2**2. / det
                                invcov(1,2) = -coeff*sig1*sig2/det
                                invcov(2,2) = sig1**2. / det
                                invcov(2,1) = invcov(1,2)
                                chisqanisos(i) = dot_product(dtheta,matmul(invcov,dtheta))
                        !        write(*,'(f5.2,1x, 4f10.5)')  zs(i), Rs(i), Rs_th(i), siglnRs(i), chisqisos(i)
                        enddo
                        chisqiso = sum(chisqisos)
                        chisqaniso = sum(chisqanisos)
                        if(pr_chisq_info) then
                           write(*,'(A,2f8.2)') '(de_chisq_bao_desi) iso/aniso = ', real(chisqiso), real(chisqaniso)
                        endif
                endif
        end subroutine de_chisq_bao_desi

  ! AP like for wCDM, using desi data
        double precision function de_chisq_ap_desi_wcdm(om, w)
                double precision, intent(in) :: om, w
                double precision :: invcov(2,2), ombf = 0.3121d0, wbf = -1.0d0, diffp(2)
                invcov(1,1) = 55553.54
                invcov(1,2) = -1558.22
                invcov(2,1) = invcov(1,2)
                invcov(2,2) = 2418.82
                diffp(1) = om - ombf
                diffp(2) = w  - wbf
                de_chisq_ap_desi_wcdm = dot_product(diffp, matmul(invcov,diffp))
         end function de_chisq_ap_desi_wcdm
  ! wCDM, planck + desi 3-bin aniso bao
        double precision function de_chisq_bao_plc_desi_wcdm(om, w)
                double precision, intent(in) :: om, w
                double precision :: invcov(2,2), ombf = 0.3121d0, wbf = -1.0d0, diffp(2)
!127435.3725561 ,  -19618.00704179],
!        [ -19618.00704179,    4252.48088164
                invcov(1,1) = 127435.37
                invcov(1,2) = -19618.01
                invcov(2,1) = invcov(1,2)
                invcov(2,2) = 4252.48
                diffp(1) = om - ombf
                diffp(2) = w  - wbf
                de_chisq_bao_plc_desi_wcdm = dot_product(diffp, matmul(invcov,diffp))
         end function de_chisq_bao_plc_desi_wcdm
  ! AP like for CPL, using desi data
        double precision function de_chisq_ap_desi_cpl(om, w0, wa)
                double precision, intent(in) :: om, w0, wa
                double precision :: invcov(3,3), ombf = 0.3121, wbf = -1.0, wabf = 0.0, diffp(3)
                integer :: i,j
                invcov(1,1) = 49775.74
                invcov(1,2) = 178.03
                invcov(1,3) = 2877.09
                invcov(2,2) = 2027.35
                invcov(2,3) = 368.48
                invcov(3,3) = 241.63
                do i = 1,3; do j = i+1, 3;
                  invcov(j,i) = invcov(i,j)
                enddo; enddo
                diffp(1) = om - ombf
                diffp(2) = w0 - wbf
                diffp(3) = wa - wabf
                de_chisq_ap_desi_cpl = dot_product(diffp, matmul(invcov,diffp))
         end function de_chisq_ap_desi_cpl
  ! CPL, planck + desi 3-bin aniso bao
        double precision function de_chisq_bao_plc_desi_cpl(om, w0, wa)
                double precision, intent(in) :: om, w0, wa
                double precision :: invcov(3,3), ombf = 0.3121, wbf = -1.0, wabf = 0.0, diffp(3)
                integer :: i,j
!129173.38927712  -19756.55366779   -2606.85526009]
! [ -19756.55366779    4246.24380324     793.19795982]
! [  -2606.85526009     793.19795982     189.19577691
                invcov(1,1) = 129173.39
                invcov(1,2) = -19756.55
                invcov(1,3) = -2606.86
                invcov(2,2) = 4246.24
                invcov(2,3) = 793.20
                invcov(3,3) = 189.20
                do i = 1,3; do j = i+1, 3;
                  invcov(j,i) = invcov(i,j)
                enddo; enddo
                !print *, invcov
                diffp(1) = om - ombf
                diffp(2) = w0 - wbf
                diffp(3) = wa - wabf
                de_chisq_bao_plc_desi_cpl = dot_product(diffp, matmul(invcov,diffp))
         end function de_chisq_bao_plc_desi_cpl
![[ 129173.38927712  -19756.55366779   -2606.85526009]
! [ -19756.55366779    4246.24380324     793.19795982]
! [  -2606.85526009     793.19795982     189.19577691]]


 

  ! AP like for Lambda CDM, using SDSS DR12
        double precision function de_chisq_sdssdr12_APinlcdm(om,filestr)
                character(len=*), INTENT(IN), OPTIONAL :: filestr
                double precision :: om, om1, om2, y1, y2
                integer :: i, i1, i2
                character(len=1000) :: inputfile
                if (.not.present(filestr)) then
                        inputfile = 'sdssdr12_AP_lcdm_chisq.mucut0.97.s20to25.table'
                else
                        inputfile = filestr
                endif
                if (.not. APlcdm_readin) then
                        open(unit=3489,file=trim(adjustl(de_data_path))//trim(adjustl(inputfile)))
                        do i = 1, 999
                                read(3489,*) APlcdm_oms(i), APlcdm_chisqs(i)
                        enddo
                        close(3489)
                        APlcdm_readin = .true.
                endif
                i1 = floor(om/0.001d0); i2 = i1+1
                om1 = APlcdm_oms(i1); om2 = APlcdm_oms(i2)
                y1  = APlcdm_chisqs(i1);  y2  = APlcdm_chisqs(i2)
                de_chisq_sdssdr12_APinlcdm = y1 + (y2-y1) * (om-om1) / (om2-om1)
                !print *, om, om1, om2
                !print *, y1, y2, de_chisq_sdssdr12_APinlcdm, y1+y2-2*de_chisq_sdssdr12_APinlcdm
        end function de_chisq_sdssdr12_APinlcdm



  ! SDSS DR14 QSO dv(1.52) = 3843 +/- 147 * (rd/rd_fid) arXiv:1705.06373
        double precision function de_chisq_SDSSDR14_QSO()
                double precision :: rdfid = 147.78d0, x, rdrescale
                rdrescale = 1.02411520606716
                ! rszd =    151.343745152605     
                ! rszd / rszd(from camb) =    1.02411520606716
                !x = 1.0/(de_rstodv(1.52d0)/rdrescale) * rdfid
                x = 1.0/(de_rstodv(1.52d0)) * 151.343745152605d0
                de_chisq_SDSSDR14_QSO = (x-3843.0d0)**2.0 / 147.0**2.0 
                if(pr_chisq_info) then
                        write(*,*) "        de_rstodv(1.52) = ", &
                                de_rstodv(1.52d0)
                        write(*,*) "     chisq_SDSSDR14_QSO (lnlike) = ", &
                                de_chisq_SDSSDR14_QSO, de_chisq_SDSSDR14_QSO/2.0
                        write(*,*)
                endif
        end function de_chisq_SDSSDR14_QSO


!####
  ! BOSS DR9 dvtors(0.57) = 13.67 +/- 0.22 arXiv:1203.6594
          double precision function de_chisq_BOSSDR9()
                  double precision :: rstoDV0p57
                  
                  rstoDV0p57 = de_rstodv(0.57d0)
                  de_chisq_BOSSDR9 = ((1.0/rstoDV0p57 - 13.67d0)/0.22d0)**2.0
                  
                          
                  if(pr_chisq_info) then
                          write(*,*) "        de_rstodv(0.57) = ", &
                                  rstoDV0p57
                          write(*,*) "     chisq_BOSSDR9 (lnlike) = ", &
                                  de_chisq_BOSSDR9, de_chisq_BOSSDR9/2.0
                          write(*,*)
                  endif
          end function de_chisq_BOSSDR9
  ! SDSS DR7 reconstructed dvtors(0.35) = 8.88 +/- 0.17 arXiv:1202.0090
        double precision function de_chisq_SDSSDR7_new()
                de_chisq_SDSSDR7_new = ((1.0d0/de_rstodv(0.35d0)-8.88d0)/0.17d0)**2.0 
                if(pr_chisq_info) then
                        write(*,*) "        de_rstodv(0.35) = ", &
                                de_rstodv(0.35d0)
                        write(*,*) "     chisq_SDSSDR7_new (lnlike) = ", &
                                de_chisq_SDSSDR7_new, de_chisq_SDSSDR7_new/2.0
                        write(*,*)
                endif
        end function de_chisq_SDSSDR7_new  
  ! SDSS DR7 arXiv:0907.1660 see wmap9 paper arXiv:1212.5226
          double precision function de_chisq_SDSSDR7_old()
                  double precision :: baodiff(2), covinv(2,2)
                  
                  baodiff(1) = de_rstodv(0.2d0)-0.1905d0
                  baodiff(2) = de_rstodv(0.35d0)-0.1097d0
                  
                  covinv(1,1) =  30124.0d0
                covinv(2,2) =  86977.0d0 
                
                covinv(1,2) = -17227.0d0
                  covinv(2,1) = covinv(1,2)
                  
                  de_chisq_SDSSDR7_old = dot_product(baodiff,matmul(covinv,baodiff))
                  
                  if(pr_chisq_info) then
                        write(*,*) "        de_rstodv(0.2) = ", &
                                de_rstodv(0.2d0)                        
                        write(*,*) "        de_rstodv(0.35) = ", &
                                de_rstodv(0.35d0)
                        write(*,*) "     chisq_SDSSDR7_old (lnlike) = ", &
                                de_chisq_SDSSDR7_old, de_chisq_SDSSDR7_old/2.0
                        write(*,*)                        
                endif
        end function de_chisq_SDSSDR7_old
  ! WiggleZ rstoDV see wmap9 paper arXiv:1212.5226
          double precision function de_chisq_WigZ_rsDV()
                  double precision :: baodiff(3), covinv(3,3)

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
                  
                if(pr_chisq_info) then
                        write(*,*) "        de_rstodv(0.44) = ", &
                                de_rstodv(0.44d0)                        
                        write(*,*) "        de_rstodv(0.60) = ", &
                                de_rstodv(0.6d0)
                        write(*,*) "        de_rstodv(0.73) = ", &
                                de_rstodv(0.73d0)                        
                        write(*,*) "     chisq_WiggleZ_rstoDV (lnlike) = ", &
                                de_chisq_WigZ_rsDV, de_chisq_WigZ_rsDV/2.0
                        write(*,*)                        
                endif
        end function de_chisq_WigZ_rsDV
  ! WiggleZ A  see wmap9 paper arXiv:1212.5226
          double precision function de_chisq_WigZ_A()
                  double precision :: baodiff(3), covinv(3,3)

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
                  
                if(pr_chisq_info) then
                        write(*,*) "        de_A_bao(0.44) = ", &
                                de_rstodv(0.44d0)                        
                        write(*,*) "        de_A_bao(0.60) = ", &
                                de_rstodv(0.6d0)
                        write(*,*) "        de_A_bao(0.73) = ", &
                                de_rstodv(0.73d0)                        
                        write(*,*) "     chisq_WiggleZ_A (lnlike) = ", &
                                de_chisq_WigZ_A, de_chisq_WigZ_A/2.0
                        write(*,*)                        
                endif
        end function de_chisq_WigZ_A
        
  !----------------------------------------------------------------
  ! used in wmap9 paper arXiv:1212.5226
  !  (6dFGS + WiggleZ rsDV + BOSS)
        double precision function de_chisq_bao_ver2()
                double precision :: baodiff(6), covinv(6,6)

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
                if(pr_chisq_info) then
                        write(*,*) "        de_rstodv(0.10) = ", de_rstodv(0.1d0)                
                        write(*,*) "        de_rstodv(0.35) = ", de_rstodv(0.35d0)
                        write(*,*) "        de_rstodv(0.57) = ", de_rstodv(0.57d0)                        
                        write(*,*) "        de_rstodv(0.44) = ", de_rstodv(0.44d0)                        
                        write(*,*) "        de_rstodv(0.60) = ", de_rstodv(0.6d0)
                        write(*,*) "        de_rstodv(0.73) = ", de_rstodv(0.73d0)                        
                        write(*,*) "     chisq_bao_ver2 (lnlike) = ", de_chisq_bao_ver2, de_chisq_bao_ver2/2.0
                        write(*,*)
                endif
        end function de_chisq_bao_ver2
        
  !----------------------------------------------------------------
  ! SDSS DR7 + WiggleZ A + 6dFGS
        double precision function de_chisq_bao_ver1()
                double precision :: z1, z2, z3, z4, z5, z6
                double precision :: p1(2), p2(3)
                 double precision :: cov1(2,2)
                double precision :: cov2(3,3)
                double precision :: chisqSDSS, chisqWiggleZ, chisq6dFGS
                integer :: i
                
                cov1(1,1) = 30124.0d0
                cov1(1,2) = -17227.0d0
                cov1(2,1) = -17227.0d0
                cov1(2,2) = 86977.0d0
                cov2(1,1) = 1040.3d0
                cov2(1,2) = -807.5d0
                cov2(1,3) = 336.8d0
                cov2(2,1) = -807.5d0
                cov2(2,2) = 3720.3d0
                cov2(2,3) = -1551.9d0
                cov2(3,1) = 336.8d0
                cov2(3,2) = -1551.9d0
                cov2(3,3) = 2914.9d0
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
                
                if(pr_chisq_info) then
                        
                        write(*,*) "     SDSS DR7 chisq (lnlike) = ", chisqSDSS, chisqSDSS/2.0
                        write(*,*) "     WiggleZ chisq (lnlike)  = ", chisqWiggleZ, chisqWiggleZ/2.0
                        write(*,*) "     6dFGS chisq (lnlike)    = ", chisq6dFGS, chisq6dFGS/2.0
                        write(*,*) "     chisq_bao_ver1 (lnlike)  = ", de_chisq_bao_ver1,  de_chisq_bao_ver1/2.0
                        write(*,*)
                endif
                
        end function de_chisq_bao_ver1

  !----------------------------------------------------------------
  ! chisq function for new bao, refer to arXiv:1312.4877 and                           !bossli
  ! 1401.0358(6dFGS + BOSS DR7 + BOSS DR11 + Imp. WiggleZ)
        double precision function de_chisq_bao_ver3()
                double precision :: baodiff(8), covinv(8,8)

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
                if(pr_chisq_info) then
                        write(*,*) "        de_rstodv(0.10) = ", de_rstodv(0.1d0)                
                        write(*,*) "        de_rstodv(0.35) = ", de_rstodv(0.35d0)
                        write(*,*) "        de_rstodv(0.32) = ", de_rstodv(0.32d0)                        
                        write(*,*) "        de_DAtord(0.57d0) = ", de_DAtord(0.57d0)                        
                        write(*,*) "        de_Hrd(0.57d0) = ", de_Hrd(0.57d0)
                        write(*,*) "        de_DVtors(0.44d0) = ", de_DVtors(0.44d0)
                        write(*,*) "        de_DVtors(0.60d0) = ", de_DVtors(0.60d0)
                        write(*,*) "        de_DVtors(0.73d0) = ", de_DVtors(0.73d0)                        
                        write(*,*) "     chisq_bao_ver3 (lnlike) = ", de_chisq_bao_ver3, de_chisq_bao_ver3/2.0
                        write(*,*)
                endif
        end function de_chisq_bao_ver3
  ! BOSS DR11 anisotropy
  !  DV = (1264 ± 25 Mpc)(rd /rd,fid ) at z = 0.32 (isotropic constraint at z=0.32).
  !  DA = (1421 ± 20 Mpc)(rd /rd,fid ) and H = (96.8 ± 3.4 km/s/Mpc)(rd,fid /rd) at z = 0.57, 
  !    with a correlation coefficient between DA and H of 0.539
  !  we shall use rd,fid = 153.19 for EH98 (using 149.28 for camb)
        double precision function de_chisq_dr11()
                double precision :: covinv(2,2), baodiff(2), rd_fid = 153.19d0
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
        double precision function de_chisq_impwig()
                double precision :: covinv(3,3), baodiff(3), rs_fid = 152.3d0
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
  ! chisq function for H
  !----------------------------------------------------------------
        double precision function de_chisq_h_Riess()
                de_chisq_h_Riess    = ((de_CP%h-0.738d0)/0.024d0)**2.0
                if(pr_chisq_info) then
                        write(*,*) "     chisq_h_Riess (lnlike) = ", de_chisq_h_Riess, de_chisq_h_Riess/2.0
                endif
        end function de_chisq_h_Riess
        

  !----------------------------------------------------------------
  ! chisq function for H
  !----------------------------------------------------------------
        double precision function de_chisq_h_Carnegie()
                de_chisq_h_Carnegie    = ((de_CP%h-0.743d0)/0.021d0)**2.0
                if(pr_chisq_info) then
                        write(*,*) "     chisq_h_Carnegie (lnlike) = ", de_chisq_h_Carnegie, de_chisq_h_Carnegie/2.0
                endif
        end function de_chisq_h_Carnegie
end module de_chisqs
