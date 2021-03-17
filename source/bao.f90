    !BAO, f_sigma8 and other measurements at single redshifts, with correlations

    !AL/JH Oct 2012: encorporate DR9 data into something close to new cosmomc format
    !Dec 2013: merged in DR11 patch (Antonio J. Cuesta, for the BOSS collaboration)
    !Sept 2014: refactored structure using new "bao_dataset[NAME] = file.dataset" syntax
    !           added refactoed DR7 MGS code (thanks Lado Samushia)
    !Jan 2014 added bao_HZ_rs, f_sigma8, F_AP, and more general .dat format
    !Feb 2016: TC 1) H(z)rs in km units,
    !             2) DR12 BAO Likelihoods from alpha_par and alpha_plel
    !             3) consensus likelihood distribution for LOWZ and CMASS
    !             4) RSD DR12 results included with QPM and MD-PATCHY inverse
    !             covariance matrices
    !             5) corrected DR11 rs_rescale factor (leads to ~0.1% bias
    !             correction at most)
    !Oct 2016: add support for dataset files no measurement_type, instead
    !              specified for each point in data file
    !PL Apr 2019: added Lyman alpha BAO from eBOSS DR12 and DR14
    !Dec 2019 EM and ADM for eBOSS collaboration

    module bao
    use MatrixUtils
    use settings
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    use IniObjects
    implicit none
    private

    character(LEN=Ini_Enumeration_Len), parameter :: measurement_types(11) = &
        [character(Ini_Enumeration_Len)::'Az','DV_over_rs','rs_over_DV','DA_over_rs', &
        'F_AP', 'f_sigma8','bao_Hz_rs','bao_Hz_rs_103','dilation','DM_over_rs', &
        'DH_over_rs']

    integer, parameter :: bao_Az =1, bao_DV_over_rs = 2, bao_rs_over_DV = 3, bao_DA_over_rs = 4, &
        F_AP= 5, f_sigma8=6, bao_Hz_rs = 7, bao_Hz_rs_103 = 8, dilation = 9, bao_DM_over_rs = 10, &
        bao_DH_over_rs = 11

    type, extends(TCosmoCalcLikelihood) :: TBAOLikelihood
        integer :: num_bao ! total number of points used
        integer, allocatable :: type_bao(:) !one of the constants defined above
        real(mcp) :: rs_rescale = 1._mcp !if not generated with numerical CAMB rs
        real(mcp), allocatable, dimension(:) :: bao_z, bao_obs, bao_err
        real(mcp), allocatable, dimension(:,:) :: bao_invcov
        real(mcp) :: Hrd_fid,DA_rd_fid
    contains
    procedure :: LogLike => BAO_LnLike
    procedure :: ReadIni => BAO_ReadIni
    procedure :: InitProbDist => BAO_InitProbDist
    procedure, private :: Acoustic
    procedure, private :: Get_rs_drag
    end type TBAOLikelihood

    Type, extends(TBAOLikelihood) :: DR1xLikelihood
        real(mcp), allocatable, dimension(:) :: alpha_perp_file,alpha_plel_file
        real(mcp), allocatable ::  prob_file(:,:)
        integer alpha_npoints_perp, alpha_npoints_plel
    contains
    procedure :: LogLike => DR1x_loglike
    procedure :: InitProbDist => DR1x_InitProbDist
    end type
    
    Type, extends(TBAOLikelihood) :: BAO_DR1xLikelihood
        real(mcp), allocatable, dimension(:) :: DM_file,DH_file
        real(mcp), allocatable ::  prob_file(:,:)
        integer DM_npoints, DH_npoints 
    contains
    procedure :: LogLike => BAO_DR1x_loglike
    procedure :: InitProbDist => BAO_DR1x_InitProbDist
    end type

    Type, extends(TBAOLikelihood) :: BAORSD_DR1xLikelihood
        real(mcp), allocatable, dimension(:) :: DM_file,DH_file,fs8_file
        real(mcp), allocatable ::  prob_file(:,:,:)
        integer DM_npoints, DH_npoints, fs8_npoints
    contains
    procedure :: LogLike => BAORSD_DR1x_loglike
    procedure :: InitProbDist => BAORSD_DR1x_InitProbDist
    end type

    Type, extends(TBAOLikelihood) :: RSD_DR1xLikelihood
        real(mcp), allocatable, dimension(:) :: fs8_file
        real(mcp), allocatable ::  prob_file(:)
        integer fs8_npoints
    contains
    procedure :: LogLike => RSD_DR1x_loglike
    procedure :: InitProbDist => RSD_DR1x_InitProbDist
    end type
    
    Type, extends(TBAOLikelihood) :: BAOiso_DR1xLikelihood
        real(mcp), allocatable, dimension(:) :: DV_over_rs_file
        real(mcp), allocatable ::  prob_file(:)
        integer DV_over_rs_npoints
    contains
    procedure :: LogLike => BAOiso_DR1x_loglike
    procedure :: InitProbDist => BAOiso_DR1x_InitProbDist
    end type


    Type, extends(TBAOLikelihood) :: MGSLikelihood
        real(mcp), allocatable :: alpha_prob(:)
    contains
    procedure :: LogLike => BAO_MGS_loglike
    procedure :: InitProbDist => BAO_MGS_InitProbDist
    end type

    real(mcp) :: BAO_fixed_rs = -1._mcp

    public BAOLikelihood_Add
    contains

    subroutine BAOLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    class(TBAOLikelihood), pointer :: this
    integer i
    Type(TSettingIni) :: DataSets, OverrideSettings

    if (.not. Ini%Read_Logical('use_BAO',.false.)) return

    call Ini%TagValuesForName('bao_dataset', DataSets, filename=.true.)
    if (DataSets%Count==0) call MpiStop('Use_BAO but no bao_dataset[NAMETAG] defined')

    if (Ini%Haskey('BAO_fixed_rs')) then
        BAO_fixed_rs= Ini%Read_Double('BAO_fixed_rs',-1._mcp)
    end if

    do i= 1, DataSets%Count
       call Ini%SettingValuesForTagName('bao_dataset',DataSets%Name(i),OverrideSettings)
       if (Datasets%Name(i)=='MGS') then
          allocate(MGSLikelihood::this)
       else if (Datasets%Name(i)=='DR11CMASS') then
          allocate(DR1xLikelihood::this)
       else if (Datasets%Name(i)=='DR12CMASS') then
          allocate(DR1xLikelihood::this)
       else if (Datasets%Name(i)=='DR12LYA') then
          allocate(DR1xLikelihood::this)
       else if (Datasets%Name(i)=='DR14LYA') then
          allocate(DR1xLikelihood::this)
       else if (Datasets%Name(i)=='DR16LYAUTO') then
          allocate(BAO_DR1xLikelihood::this)
       else if (Datasets%Name(i)=='DR16LYxQSO') then
          allocate(BAO_DR1xLikelihood::this)
       else if (Datasets%Name(i)=='DR16ELGbao') then
          allocate(BAOiso_DR1xLikelihood::this)
       else if (Datasets%Name(i)=='DR16ELGrsd') then
          allocate(RSD_DR1xLikelihood::this)
       else if (Datasets%Name(i)=='DR16ELGbaoFS') then
          allocate(BAORSD_DR1xLikelihood::this)
       else
          allocate(TBAOLikelihood::this)
       end if
       call this%ReadDatasetFile(Datasets%Value(i),OverrideSettings)
       this%tag = Datasets%Name(i)
       this%LikelihoodType = 'BAO'
       this%needs_background_functions = .true.
       call LikeList%Add(this)
    end do
    
    if (Feedback>1) write(*,*) 'read BAO data sets'

    end subroutine BAOLikelihood_Add

    subroutine BAO_ReadIni(this, Ini)
    class(TBAOLikelihood) this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: bao_measurement, bao_measurements_file
    integer i,iopb
    real (mcp) :: rd_fid,H_fid,DA_fid
    Type(TTextFile) :: F
    logical :: hastype, haserror
    integer :: status
    character(LEN=Ini_Enumeration_Len) :: tp
    character(LEN=:), allocatable :: InLine

    if (Feedback > 0 .and. MpiRank==0) write (*,*) 'reading BAO data set: '//trim(this%name)
    this%num_bao = Ini%Read_Int('num_bao',1)

    call Ini%Read('rs_rescale',this%rs_rescale)

    allocate(this%bao_z(this%num_bao))
    allocate(this%bao_obs(this%num_bao))
    allocate(this%bao_err(this%num_bao))

    hasType = Ini%HasKey('measurement_type')
    if (hasType) then
        call Ini%Read_Enumeration_List('measurement_type',measurement_types, this%type_bao, nvalues = this%num_bao)
    else
        allocate(this%type_bao(this%num_bao))
    end if

    if ((Ini%HasKey('Hrd_fid')).and.(Ini%HasKey('DA_rd_fid'))) then
        this%Hrd_fid = Ini%Read_Double('Hrd_fid')
        this%DA_rd_fid = Ini%Read_Double('DA_rd_fid')
    else if ((Ini%HasKey('rd_fid')).and.(Ini%HasKey('H_fid')).and.(Ini%HasKey('DA_fid'))) then
        rd_fid = Ini%Read_Double('rd_fid')
        H_fid = Ini%Read_Double('H_fid')
        DA_fid = Ini%Read_Double('DA_fid')
        this%Hrd_fid = H_fid*rd_fid
        this%DA_rd_fid = DA_fid/rd_fid
    end if

    if (Ini%HasKey('zeff')) then
        bao_measurement  = Ini%Read_String('zeff')
        read (bao_measurement,*, iostat=status) this%bao_z
        if (status/=0) then !assume just one redshift
            this%bao_z = Ini%Read_Double('zeff')
        end if
        if (this%type_bao(1) /= dilation) then
            bao_measurement  = Ini%Read_String('bao_measurement')
            if (this%num_bao>1) then
                read (bao_measurement,*) this%bao_obs(:)
            else
                read (bao_measurement,*) this%bao_obs(1),this%bao_err(1)
            end if
        end if
    else
        bao_measurements_file = Ini%ReadRelativeFileName('bao_measurements_file')
        haserror = Ini%Read_Logical('bao_measurements_file_has_error',.true.)
        call F%Open(bao_measurements_file)
        do i=1,this%num_bao
            if (F%ReadLineSkipEmptyAndComments(InLine)) then
                if (haserror) then
                    if (hasType) then
                        read (InLine,*, iostat=iopb) this%bao_z(i),this%bao_obs(i),this%bao_err(i)
                    else
                        read (InLine,*, iostat=iopb) this%bao_z(i),this%bao_obs(i),this%bao_err(i), tp
                        this%type_bao(i) = Ini%EnumerationValue(tp, measurement_types)
                    end if
                else
                    if (hasType) then
                        read (InLine,*, iostat=iopb) this%bao_z(i),this%bao_obs(i)
                    else
                        read (InLine,*, iostat=iopb) this%bao_z(i),this%bao_obs(i), tp
                        this%type_bao(i) = Ini%EnumerationValue(tp, measurement_types)
                    end if
                end if
                if (iopb /= 0) call MpiStop('BAO_ReadIni: Error reading bao_measurements_file: ' &
                    //trim(this%name))
            else
                call MpiStop('BAO_ReadIni: Missing line in bao_measurements_file: ' &
                    //trim(this%name))
            end if
        end do
        call F%Close()
    end if
    if (any(this%bao_z< 0.0001)) call MpiStop('Error reading BAO measurements')

    this%needs_powerspectra =  any(this%type_bao == f_sigma8)

    if (this%needs_powerspectra) then
        this%needs_exact_z =  .true.
        allocate(this%exact_z(this%num_bao))
        this%exact_z = this%bao_z
    end if


    call this%InitProbDist(Ini)

    end subroutine BAO_ReadIni

    subroutine BAO_InitProbDist(this, Ini)
    class(TBAOLikelihood) this
    class(TSettingIni) :: Ini
    integer i
    character(LEN=:), allocatable ::bao_invcov_file

    allocate(this%bao_invcov(this%num_bao,this%num_bao))
    this%bao_invcov=0

    if (Ini%HasKey('bao_invcov_file')) then
        bao_invcov_file  = Ini%ReadRelativeFileName('bao_invcov_file')
        call File%ReadTextMatrix(bao_invcov_file, this%bao_invcov)
    else if (Ini%HasKey('bao_cov_file')) then
        bao_invcov_file  = Ini%ReadRelativeFileName('bao_cov_file')
        call File%ReadTextMatrix(bao_invcov_file, this%bao_invcov)
        call Matrix_Inverse(this%bao_invcov)
    else
        do i=1,this%num_bao
            !diagonal, or actually just 1..
            this%bao_invcov(i,i) = 1/this%bao_err(i)**2
        end do
    end if

    end subroutine BAO_InitProbDist


    real(mcp) function get_rs_drag(this,Theory)
    class(TBAOLikelihood) :: this
    Class(TCosmoTheoryPredictions), target :: Theory

    if (BAO_fixed_rs>0) then
        !this is just for use for e.g. BAO 'only' constraints
        get_rs_drag =  BAO_fixed_rs
    else
        get_rs_drag =  Theory%derived_parameters( derived_rdrag )
    end if

    end function

    function Acoustic(this,CMB,z)
    class(TBAOLikelihood) :: this
    class(CMBParams) CMB
    real(mcp) Acoustic
    real(mcp), intent(IN) :: z
    real(mcp) omh2,ckm,omegam,h

    omegam = 1.d0 - CMB%omv - CMB%omk
    h = CMB%h0/100
    ckm = const_c/1e3_mcp !JD c in km/s

    omh2 = omegam*h**2.d0
    Acoustic = 100*this%Calculator%BAO_D_v(z)*sqrt(omh2)/(ckm*z)
    end function Acoustic

    function BAO_LnLike(this, CMB, Theory, DataParams)
    Class(TBAOLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    integer j
    real(mcp) BAO_LnLike, rs, z, f
    real(mcp)  :: BAO_theory(this%num_bao)

    rs = this%get_rs_drag(Theory) * this%rs_rescale

    do j=1, this%num_bao
        z= this%bao_z(j)

        select case(this%type_bao(j))
        case (bao_DV_over_rs)
            BAO_theory(j) = this%Calculator%BAO_D_v(z)/rs
        case (bao_Hz_rs)
            BAO_theory(j) = this%Calculator%Hofz_Hunit(z)*rs
        case (bao_Hz_rs_103)
            BAO_theory(j) = this%Calculator%Hofz_Hunit(z)*rs*1.0d-3
        case (bao_rs_over_DV)
            BAO_theory(j) = rs/this%Calculator%BAO_D_v(z)
        case (bao_Az)
            BAO_theory(j) = this%Acoustic(CMB,z)
        case (bao_DA_over_rs)
            BAO_theory(j) = this%Calculator%AngularDiameterDistance(z)/rs
        case (bao_DM_over_rs)
            BAO_theory(j) = (1+z)*this%Calculator%AngularDiameterDistance(z)/rs
        case (F_AP)
            BAO_theory(j) = (1+z)*this%Calculator%AngularDiameterDistance(z)* &
                this%Calculator%Hofz(z)
        case (bao_DH_over_rs)
            BAO_theory(j) = 1./this%Calculator%Hofz(z)/rs
        case (f_sigma8) 
            BAO_theory(j) = Theory%growth_z%Value(z)
         case default
            call MpiStop('BAO_LnLike: Unknown type_bao')
        end select
    end do

    BAO_theory = BAO_theory - this%bao_obs
    BAO_LnLike = Matrix_QuadForm(this%bao_invcov,BAO_theory) / 2

    if(feedback>1) write(*,*) trim(this%name)//' BAO likelihood = ', BAO_LnLike

    end function BAO_LnLike


    !!!DR11/DR12 CMASS/LOWZ/LYA

    subroutine DR1x_InitProbDist(this, Ini)
    class(DR1xLikelihood) this
    class(TSettingIni) :: Ini
    real(mcp) :: tmp0,tmp1,tmp2
    integer ios,ii,jj
    Type(TTExtFile) F
    integer :: alpha_npoints_plel, alpha_npoints_perp
    logical is_square

    is_square = Ini%Read_Logical('is_square',.true.)
    alpha_npoints_perp = Ini%Read_Int('alpha_npoints')
    if (is_square .eqv. .true.) then
        alpha_npoints_plel = alpha_npoints_perp
    else
        alpha_npoints_plel = Ini%Read_Int('alpha_npoints_plel')
    end if

    allocate(this%alpha_perp_file(alpha_npoints_perp),this%alpha_plel_file(alpha_npoints_plel))
    allocate(this%prob_file(alpha_npoints_perp,alpha_npoints_plel))

    call F%Open(Ini%ReadRelativeFileName('prob_dist'))
    do ii=1, alpha_npoints_perp
        do jj=1, alpha_npoints_plel
            read (F%unit,*,iostat=ios) tmp0,tmp1,tmp2
            if (ios /= 0) call MpiStop('Error reading BAO file')
            this%alpha_perp_file(ii)   = tmp0
            this%alpha_plel_file(jj)   = tmp1
            this%prob_file(ii,jj)      = tmp2
        end do
    end do
    call F%Close()

    !Normalize distribution (so that the peak value is 1.0)

    this%prob_file=this%prob_file/ maxval(this%prob_file)
    this%alpha_npoints_perp = alpha_npoints_perp
    this%alpha_npoints_plel = alpha_npoints_plel

    end subroutine DR1x_InitProbDist

    subroutine BAO_DR1x_InitProbDist(this, Ini)
    class(BAO_DR1xLikelihood) this
    class(TSettingIni) :: Ini
    real(mcp) :: tmp0,tmp1,tmp2
    integer ios,ii,jj
    Type(TTExtFile) F
    integer :: DH_npoints, DM_npoints
    logical is_square
    is_square = Ini%Read_Logical('is_square',.true.)
    DM_npoints = Ini%Read_Int('DM_npoints')
    if (is_square .eqv. .true.) then
        DH_npoints = DM_npoints
    else
        DH_npoints = Ini%Read_Int('DH_npoints')
    end if

    allocate(this%DH_file(DM_npoints),this%DM_file(DH_npoints))
    allocate(this%prob_file(DM_npoints,DH_npoints))

    call F%Open(Ini%ReadRelativeFileName('prob_dist'))
    read (F%unit,*,iostat=ios) tmp0,tmp1,tmp2
    do ii=1, DM_npoints
        do jj=1, DH_npoints
            read (F%unit,*,iostat=ios) tmp0,tmp1,tmp2
            if (ios /= 0) call MpiStop('Error reading BAO file')
            this%DM_file(ii)           = tmp0
            this%DH_file(jj)           = tmp1
            this%prob_file(ii,jj)      = tmp2
        end do
    end do
    call F%Close()
    !Normalize distribution (so that the peak value is 1.0)

    this%prob_file=this%prob_file/ maxval(this%prob_file)
    this%DM_npoints = DM_npoints
    this%DH_npoints = DH_npoints

    end subroutine BAO_DR1x_InitProbDist

    function DR1x_loglike(this, CMB, Theory, DataParams)
    Class(DR1xLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real (mcp) z, DR1x_loglike, alpha_perp, alpha_plel, prob
    integer i, j, ii,jj
    real(mcp) rsdrag_theory

    z = this%bao_z(1)
    rsdrag_theory = this%get_rs_drag(Theory)
    alpha_perp=(this%Calculator%AngularDiameterDistance(z)/rsdrag_theory)/(this%DA_rd_fid)!CMASS/LOWZ/LYA
    alpha_plel=(this%Hrd_fid)/((this%Calculator%Hofz_Hunit(z))*rsdrag_theory)!CMASS/LOWZ/LYA
    if ((alpha_perp < this%alpha_perp_file(1)).or.(alpha_perp > this%alpha_perp_file(this%alpha_npoints_perp-1)).or. &
        &   (alpha_plel < this%alpha_plel_file(1)).or.(alpha_plel > this%alpha_plel_file(this%alpha_npoints_plel-1))) then
        DR1x_loglike = logZero
    else
        do i=1,this%alpha_npoints_perp
            if (alpha_perp - this%alpha_perp_file(i) .le. 0) then
                ii = i-1
                exit
            end if
        end do
        do j=1,this%alpha_npoints_plel
            if (alpha_plel - this%alpha_plel_file(j) .le. 0) then
                jj = j-1
                exit
            end if
        end do
        prob=(1./((this%alpha_perp_file(ii+1)-this%alpha_perp_file(ii))*(this%alpha_plel_file(jj+1)-this%alpha_plel_file(jj))))*  &
            &       (this%prob_file(ii,jj)*(this%alpha_perp_file(ii+1)-alpha_perp)*(this%alpha_plel_file(jj+1)-alpha_plel) &
            &       -this%prob_file(ii+1,jj)*(this%alpha_perp_file(ii)-alpha_perp)*(this%alpha_plel_file(jj+1)-alpha_plel) &
            &       -this%prob_file(ii,jj+1)*(this%alpha_perp_file(ii+1)-alpha_perp)*(this%alpha_plel_file(jj)-alpha_plel) &
            &       +this%prob_file(ii+1,jj+1)*(this%alpha_perp_file(ii)-alpha_perp)*(this%alpha_plel_file(jj)-alpha_plel))
        if  (prob > 0) then
            DR1x_loglike = -log( prob )
        else
            DR1x_loglike = logZero
        endif
    endif
    end function DR1x_loglike

    !!!! SDSS DR7 main galaxy sample http://arxiv.org/abs/1409.3242
    !Adapted from code by Lado Samushia
   function BAO_DR1x_loglike(this, CMB, Theory, DataParams)
    Class(BAO_DR1xLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real (mcp) z, BAO_DR1x_loglike, DM, DH, prob
    integer i, j, ii,jj
    real(mcp) rsdrag_theory

    z = this%bao_z(1)
    rsdrag_theory = this%get_rs_drag(Theory)
    DM =(1+z)*this%Calculator%AngularDiameterDistance(z)/rsdrag_theory
    DH = 1./this%Calculator%Hofz(z)/rsdrag_theory
    if ((DM < this%DM_file(1)).or.(DM > this%DM_file(this%DM_npoints-1)).or. &
        &   (DH < this%DH_file(1)).or.(DH > this%DH_file(this%DH_npoints-1))) then
        BAO_DR1x_loglike = logZero
    else
        do i=1,this%DM_npoints
            if (DM - this%DM_file(i) .le. 0) then
                ii = i-1
                exit
            end if
        end do
        do j=1,this%DH_npoints
            if (DH - this%DH_file(j) .le. 0) then
                jj = j-1
                exit
            end if
        end do
        prob=(1./((this%DM_file(ii+1)-this%DM_file(ii))*(this%DH_file(jj+1)-this%DH_file(jj))))*&
            &(this%prob_file(ii,jj)*(this%DM_file(ii+1)-DM)*(this%DH_file(jj+1)-DH)&
            &-this%prob_file(ii+1,jj)*(this%DM_file(ii)-DM)*(this%DH_file(jj+1)-DH)&
            &-this%prob_file(ii,jj+1)*(this%DM_file(ii+1)-DM)*(this%DH_file(jj)-DH)&
            &+this%prob_file(ii+1,jj+1)*(this%DM_file(ii)-DM)*(this%DH_file(jj)-DH))
        if  (prob > 0) then
            BAO_DR1x_loglike = -log( prob )
        else
            BAO_DR1x_loglike = logZero
        endif
    endif
    end function BAO_DR1x_loglike
 
    subroutine BAO_MGS_InitProbDist(this, Ini)
    class(MGSLikelihood) this
    class(TSettingIni) :: Ini

    call File%LoadTxt(Ini%ReadRelativeFileName('prob_dist'), this%alpha_prob)

    end subroutine BAO_MGS_InitProbDist


    function BAO_MGS_loglike(this, CMB, Theory, DataParams)
    Class(MGSLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) BAO_MGS_loglike
    real (mcp) alphamgs, chi2
    real(mcp),parameter :: rsfidmgs = 148.69_mcp, DVfidmgs = 638.9518_mcp
    real(mcp), parameter :: alpha_min=0.8005_mcp, alpha_max = 1.1985_mcp
    integer ii

    alphamgs =   this%Calculator%BAO_D_v(this%bao_z(1))/this%get_rs_drag(Theory) / (DVfidmgs / rsfidmgs)
    if ((alphamgs > alpha_max).or.(alphamgs < alpha_min)) then
        BAO_MGS_loglike = logZero
    else
        ii = 1+floor((alphamgs - alpha_min)/0.001)
        chi2 = (this%alpha_prob(ii) + this%alpha_prob(ii+1))/2.0
        BAO_MGS_loglike = chi2/2.0
    endif

    end function BAO_MGS_loglike

    subroutine BAORSD_DR1x_InitProbDist(this, Ini)
    class(BAORSD_DR1xLikelihood) this
    class(TSettingIni) :: Ini
    real(mcp) :: tmp0,tmp1,tmp2,tmp3
    integer ios,ii,jj,kk
    Type(TTExtFile) F
    integer :: DM_npoints, DH_npoints, fs8_npoints
    DM_npoints = Ini%Read_Int('DM_npoints')
    DH_npoints = Ini%Read_Int('DH_npoints')
    fs8_npoints = Ini%Read_Int('fs8_npoints')

    allocate(this%DM_file(DM_npoints),this%DH_file(DH_npoints),this%fs8_file(fs8_npoints))
    allocate(this%prob_file(DM_npoints,DH_npoints,fs8_npoints))

    call F%Open(Ini%ReadRelativeFileName('prob_dist'))
    do ii=1, DM_npoints
        do jj=1, DH_npoints
            do kk=1, fs8_npoints
                read (F%unit,*,iostat=ios) tmp0,tmp1,tmp2,tmp3
                if (ios /= 0) call MpiStop('Error reading RSD file')
                this%DM_file(ii)         = tmp0
                this%DH_file(jj)          = tmp1
                this%fs8_file(kk)        = tmp2
                this%prob_file(ii,jj,kk) = tmp3
            end do
        end do
    end do
    call F%Close()

    !Normalize distribution (so that the peak value is 1.0)

    this%prob_file=this%prob_file/ maxval(this%prob_file)
    this%DM_npoints = DM_npoints
    this%DH_npoints = DH_npoints
    this%fs8_npoints = fs8_npoints

    end subroutine BAORSD_DR1x_InitProbDist

    function BAORSD_DR1x_loglike(this, CMB, Theory, DataParams)
    Class(BAORSD_DR1xLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real (mcp) z, BAORSD_DR1x_loglike, DM, DH, fs8, prob, dx, dy, dz
    integer i, j, k, ii, jj, kk
    real(mcp) rsdrag_theory

    z = this%bao_z(1)
    rsdrag_theory = this%get_rs_drag(Theory)
    DM = (1+z)*this%Calculator%AngularDiameterDistance(z)/rsdrag_theory
    DH = 1./this%Calculator%Hofz(z)/rsdrag_theory
    fs8 = Theory%growth_z%Value(z)
    if ((DM < this%DM_file(1)).or.(DM > this%DM_file(this%DM_npoints-1)).or. &
        & (DH < this%DH_file(1)).or.(DH > this%DH_file(this%DH_npoints-1)).or. &
        & (fs8 < this%fs8_file(1)).or.(fs8 > this%fs8_file(this%fs8_npoints-1))) then
        BAORSD_DR1x_loglike = logZero
    else
        do i=1,this%DM_npoints
            if (DM - this%DM_file(i) .le. 0) then
                ii = i-1
                exit
            end if
        end do
        do j=1,this%DH_npoints
            if (DH - this%DH_file(j) .le. 0) then
                jj = j-1
                exit
            end if
        end do
        do k=1,this%fs8_npoints
            if (fs8 - this%fs8_file(k) .le. 0) then
                kk = k-1
                exit
            end if
        end do
        dx = (DM-this%DM_file(ii))/(this%DM_file(ii+1)-this%DM_file(ii))
        dy = (DH-this%DH_file(jj))/(this%DH_file(jj+1)-this%DH_file(jj))
        dz = (fs8-this%fs8_file(kk))/(this%fs8_file(kk+1)-this%fs8_file(kk))
        prob = this%prob_file(ii,jj,kk)*(1-dx)*(1-dy)*(1-dz) + this%prob_file(ii+1,jj,kk)*dx*(1-dy)*(1-dz) &
            & + this%prob_file(ii,jj+1,kk)*(1-dx)*dy*(1-dz) + this%prob_file(ii,jj,kk+1)*(1-dx)*(1-dy)*dz &
            & + this%prob_file(ii+1,jj+1,kk)*dx*dy*(1-dz) + this%prob_file(ii,jj+1,kk+1)*(1-dx)*dy*dz &
            & + this%prob_file(ii+1,jj,kk+1)*dx*(1-dy)*dz + this%prob_file(ii+1,jj+1,kk+1)*dx*dy*dz
        if  (prob > 0) then
            BAORSD_DR1x_loglike = -log( prob )
        else
            BAORSD_DR1x_loglike = logZero
        endif
    endif
    end function BAORSD_DR1x_loglike

    subroutine BAOiso_DR1x_InitProbDist(this, Ini)
    class(BAOiso_DR1xLikelihood) this
    class(TSettingIni) :: Ini
    real(mcp) :: tmp0,tmp1
    integer ios,ii,jj,kk
    Type(TTExtFile) F
    integer :: DV_over_rs_npoints

    DV_over_rs_npoints = Ini%Read_Int('DV_over_rs_npoints')

    allocate(this%DV_over_rs_file(DV_over_rs_npoints))
    allocate(this%prob_file(DV_over_rs_npoints))

    call F%Open(Ini%ReadRelativeFileName('prob_dist'))
    do ii=1, DV_over_rs_npoints
        read (F%unit,*,iostat=ios) tmp0,tmp1
        if (ios /= 0) call MpiStop('Error reading RSD file')
        this%DV_over_rs_file(ii)         = tmp0
        this%prob_file(ii) = tmp1
    end do
    call F%Close()

    !Normalize distribution (so that the peak value is 1.0)

    this%prob_file=this%prob_file/ maxval(this%prob_file)
    this%DV_over_rs_npoints = DV_over_rs_npoints

    end subroutine BAOiso_DR1x_InitProbDist

    function BAOiso_DR1x_loglike(this, CMB, Theory, DataParams)
    Class(BAOiso_DR1xLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real (mcp) z, BAOiso_DR1x_loglike, DV_over_rs, prob
    integer i, ii
    real(mcp) rsdrag_theory

    z = this%bao_z(1)
    rsdrag_theory = this%get_rs_drag(Theory)
    DV_over_rs= this%Calculator%BAO_D_v(z)/rsdrag_theory
    if ((DV_over_rs < this%DV_over_rs_file(1)).or.&
        & (DV_over_rs > this%DV_over_rs_file(this%DV_over_rs_npoints-1))) then
        BAOiso_DR1x_loglike = logZero
    else
        do i=1,this%DV_over_rs_npoints
            if (DV_over_rs - this%DV_over_rs_file(i) .le. 0) then
                ii = i-1
                exit
            end if
        end do
        prob=(1./((this%DV_over_rs_file(ii+1)-this%DV_over_rs_file(ii))))*&
            &(this%prob_file(ii)*(this%DV_over_rs_file(ii+1)-DV_over_rs)&
            &-this%prob_file(ii+1)*(this%DV_over_rs_file(ii)-DV_over_rs))
        if  (prob > 0) then
            BAOiso_DR1x_loglike = -log( prob )
        else
            BAOiso_DR1x_loglike = logZero
        endif
    endif
    end function BAOiso_DR1x_loglike
   
    subroutine RSD_DR1x_InitProbDist(this, Ini)
    class(RSD_DR1xLikelihood) this
    class(TSettingIni) :: Ini
    real(mcp) :: tmp0,tmp1
    integer ios,ii,jj,kk
    Type(TTExtFile) F
    integer :: fs8_npoints

    fs8_npoints = Ini%Read_Int('fs8_npoints')

    allocate(this%fs8_file(fs8_npoints))
    allocate(this%prob_file(fs8_npoints))

    call F%Open(Ini%ReadRelativeFileName('prob_dist'))
    do ii=1, fs8_npoints
        read (F%unit,*,iostat=ios) tmp0,tmp1
        if (ios /= 0) call MpiStop('Error reading RSD file')
        this%fs8_file(ii) = tmp0
        this%prob_file(ii) = tmp1
    end do
    call F%Close()

    !Normalize distribution (so that the peak value is 1.0)

    this%prob_file=this%prob_file/ maxval(this%prob_file)
    this%fs8_npoints = fs8_npoints

    end subroutine RSD_DR1x_InitProbDist
 
    function RSD_DR1x_loglike(this, CMB, Theory, DataParams)
    Class(RSD_DR1xLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) :: z, RSD_DR1x_loglike, fs8, prob
    integer i, ii

    z = this%bao_z(1)
    fs8 = Theory%growth_z%Value(z)
    if ((fs8 < this%fs8_file(1)).or.&
        & (fs8 > this%fs8_file(this%fs8_npoints-1))) then
        RSD_DR1x_loglike = logZero
    else
        do i=1,this%fs8_npoints
            if (fs8 - this%fs8_file(i) .le. 0) then
                ii = i-1
                exit
            end if
        end do
        prob=(1./((this%fs8_file(ii+1)-this%fs8_file(ii))))*&
            &(this%prob_file(ii)*(this%fs8_file(ii+1)-fs8)&
            &-this%prob_file(ii+1)*(this%fs8_file(ii)-fs8))
        if  (prob > 0) then
            RSD_DR1x_loglike = -log( prob )
        else
            RSD_DR1x_loglike = logZero
        endif
    endif
    end function RSD_DR1x_loglike
    end module bao
