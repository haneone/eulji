!*******************************************************************************
! {kNN_MI}
!
! Input:	SNP genotype, Multivariate Quantitative Trait phenotype samples
!       
! Output:	kNN Mutual Information
!
! version:
!   kNN_MI_v1: based on KDE_MI, implemented kNN entropy estimation, averaged over all possible k
!   kNN_MI_v2:
!   kNN_MI_v3:
!
! Note: 
!   1. v1: p-value estimation has not been implemented, yet.
!   2. v2: same as v1. adapted only for the reading of 'survMDR' simulation data, which has the following structure
!           - 20 SNPs, 2 Phenotypes(far to the right)
!           - 400 samples for each set, 100 sets, 5 models, 2 MAFs => 400x100x5x2=400,000 lines in a single file
!           - line number at the first column
!   3. v3: based on v1.
!           - implemented p-value estimation by permutation, where the null dist is constructed by max values.
!
! Command Line Arguments:   0: program name
!                           1: input file name
!                           2: signal for header
!                     --->  3: phenotype position('R' fixed) ---> set number
!                           4: the number of samples
!                           5: the number of phenotypes
!                           6: the number of SNPs
!                           7: interaction dimension(upto 'max_dim_set')
!                           8: output filename
!                           9: the number of top rankers to write('0' for all)
!                           10: signal for appending('1' to append)
!                           11: signal for permutation('1' to permute)
!                           [12]: appending file name, if if_append = 1.
!                           [12 or 13]: the number of resampling; 12 if if_append /= 1, 13 if if_append = 1.
!*******************************************************************************

!*******************************************************************************
! module stuffs
!*******************************************************************************

module stuffs

implicit none

save

!-------------------------------------------------------------------------------
! General numbers
!-------------------------------------------------------------------------------
    real, parameter :: epsilon = 1.0E-7             ! small number
    real, parameter :: LargeNum = 1.0E+7            ! large number
    real, parameter :: pi = 3.14159265358979        ! pi
    real, parameter :: a_half = 0.5                 ! 0.5
    real, parameter :: undefined_value = -1.0       ! to denote the undefined element
    
!-------------------------------------------------------------------------------
! Parameters
!-------------------------------------------------------------------------------
	integer, parameter :: n_gt = 3				    ! number of genotypes; 0(AA), 1(Aa), 2(aa) - minor allele additive
	integer, parameter :: max_dim_set = 3		    ! upper limit of n_dim of this program(extendable)
                                                    ! to change the max dim, change max_dim_set and also add for that dim to "get_n_total_pair", "get_kNN_MI", "write_tops", "write_tops_permuted"
	integer, parameter :: n_top_set = 100		    ! set # of top association SNP pairs
    integer, parameter :: infile_unit = 1
    integer, parameter :: outfile_unit = 9
    integer, parameter :: appendfile_unit = 11
    
!-------------------------------------------------------------------------------
! Parameters to input
!-------------------------------------------------------------------------------
    character(len=1) :: position_phenotype          ! position of the outcome column('L/R' for Leftmost/Rightmost)
    integer :: if_permute                           ! perform permutation if set to 1
    integer :: if_append                            ! signal for appending ('1' for appending)
    integer :: if_header                            ! signal for header ('1' if header exists in the input file)
	integer :: n_dim		                        ! interaction dimension to estimate
	integer :: n_sample                             ! # of samples
  	integer :: n_SNP                                ! # SNPs in this microarray data  
    integer :: n_top_write                          ! # of top combinations to write out
    integer :: n_resample_set                       ! # of resampling for permutation (with replacement)
    integer :: n_phenotype                         ! multivariate phenotype
    integer :: n_set                               ! set number to estimate in the input 'survMDR' file
    integer :: neg_seed                            ! seed number for function 'ran2'
    
!-------------------------------------------------------------------------------
! Data holding arrays
!-------------------------------------------------------------------------------
	integer(kind=1), allocatable, dimension(:,:) :: genotype_array_save     ! (:,:)=(n_sample,n_SNP), save the original genotype data
    real, allocatable, dimension(:,:) :: phenotype_array_save                 ! (:,:)=(n_sample,n_phenotype), save multivariate quantitative phenotype
    real, allocatable, dimension(:,:) :: log_rho_table_save                       ! (:,:)=(n_sample,n_sample), estimate all possible distance bet the multivariate phenotype vectors and save
    real :: H_X_kNN                                                             ! H(P): overall entropy estimated solely with multivariate phenotype values
    
!-------------------------------------------------------------------------------
! Data structure holding the estimated results
!-------------------------------------------------------------------------------
    type kNN_sort_and_write
        integer, dimension(max_dim_set) :: SNP      ! SNP combination
        real :: H_X_kNN                                 ! H(P)
        real :: H_XGY                               ! H(P|G)
        real :: kNN_MI                              ! mutual information; decision making variable
        integer :: n_sample_valid                   ! may change SNP by SNP; cannot be global unlike n_sample
        real :: std_kNN_MI                          ! standardized kNN_MI after permutative iteration
        real :: p_value                             ! multiple comparison considered
    end type kNN_sort_and_write
    
    type (kNN_sort_and_write), allocatable, dimension(:) :: kNN_top_ranks   ! (:)=(n_top_write)
    
    real, allocatable, dimension(:) :: null_kNN_MI_by_max       ! (:)=(n_resample_set), holding max kNN_MI from each iteration among all SNP combination
    real :: null_kNN_MI_mean_by_max                 ! mean value of the array of null_kNN_MI_by_max
    real :: null_kNN_MI_stdev_by_max                ! standard deviation of the array of null_kNN_MI_by_max
    
!-------------------------------------------------------------------------------
! Common variables
!-------------------------------------------------------------------------------
	character(len=180) :: infile	 	            ! data filename
	character(len=180) :: outfile	                ! output filename
    character(len=180) :: appendingfile             ! appending filename
    
Contains
!-------------------------------------------------------------------------------
! Functions
!	integer function n_pick: returns a random integer from 1 to n_upto
!-------------------------------------------------------------------------------
	integer function n_pick(n_upto)

	implicit none

	integer, intent(in) :: n_upto
	real :: r_num

	call random_number(r_num)
	n_pick = int(r_num*float(n_upto)) + 1

	end function n_pick

end module stuffs

!*******************************************************************************
! program kNN_MI
!*******************************************************************************
! Input: 
! Output: 
!-------------------------------------------------------------------------------

program	kNN_MI

use stuffs

implicit none

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
call initial_steps

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
call open_and_read

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
call main_steps

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
call cleanups

print *, 'All done!'

end program kNN_MI

!*******************************************************************************
! subroutine initial_steps
!*******************************************************************************
! Accept inputs
!-------------------------------------------------------------------------------

subroutine initial_steps

use stuffs

implicit none

! Command line argument variables
!-------------------------------------------------------------------------------
character(len=180) :: comm_arg_buff    ! command line argument bufer

real :: ran2    ! external function; random number generator bet 0 and 1, exclusive of both ends
integer, parameter :: dummy_int = 10000   ! seed generation purpose using 'n_pick'
real :: r_dummy     ! random number generator warming up purpose
integer :: i       ! counter

! Get data file name
write(*,*) 'Input data filename = '
! read(*,'(A)') infile
call getarg(1, comm_arg_buff)
read(comm_arg_buff,*) infile

! Get the signal for the header in input file
write(*,*) 'Signal for the header in input file(''1'' if header exist) = '
! read(*,'(A)') if_header
call getarg(2, comm_arg_buff)
read(comm_arg_buff,*) if_header

! Get the set number to estimate in the input 'survMDR' file
write(*,*) 'Input the set number to estimate in the input ''survMDR'' file = '
! read(*,'(A)') n_set
call getarg(3, comm_arg_buff)
read(comm_arg_buff,*) n_set

!! Get the signal for the position of outcome(trait) in input file
!write(*,*) 'Signal for the position of outcome in input file(''L/R'' for Leftmost/Rightmost column) = '
!! read(*,'(A)') position_phenotype
!call getarg(3, comm_arg_buff)
!read(comm_arg_buff,*) position_phenotype

! Get the number of samples
write(*,*) 'Input the number of samples = '
! read(*,'(A)') n_sample
call getarg(4, comm_arg_buff)
read(comm_arg_buff,*) n_sample

! Get the adaptive kernel sensitivity parameter, alpha
write(*,*) 'Input the number of phenotypes = '
! read(*,'(A)') n_phenotype
call getarg(5, comm_arg_buff)
read(comm_arg_buff,*) n_phenotype

! Get the number of SNPs
write(*,*) 'Input the number of SNPs = '
! read(*,'(A)') n_SNP
call getarg(6, comm_arg_buff)
read(comm_arg_buff,*) n_SNP

! Get n_dim to estimate
write(*,*) 'Input n_dim = '
! read(*,'(A)') n_dim
call getarg(7, comm_arg_buff)
read(comm_arg_buff,*) n_dim

! Get the output filename
write(*,*) 'Output filename = '
! read(*,'(A)') outfile
call getarg(8, comm_arg_buff)
read(comm_arg_buff,*) outfile

! Get the number of output top ranks
write(*,*) 'Number of top ranks to write(''0'' to include all) = '
! read(*,'(A)') n_top_write
call getarg(9, comm_arg_buff)
read(comm_arg_buff,*) n_top_write

! Get the signal for appending
write(*,*) 'Signal for appending(''1'' to append) = '
! read(*,'(A)') if_append
call getarg(10, comm_arg_buff)
read(comm_arg_buff,*) if_append

! Get the signal for permutation
write(*,*) 'Signal for permutation(''1'' to permute) = '
! read(*,'(A)') if_permute
call getarg(11, comm_arg_buff)
read(comm_arg_buff,*) if_permute

! Get the appending filename
if (if_append == 1) then
    write(*,*) 'Appending filename = '
!    read(*,'(A)') appendingfile
    call getarg(12, comm_arg_buff)
    read(comm_arg_buff,*) appendingfile
endif

! Get the number of resampling
if (if_permute == 1) then
    write(*,*) 'Number of resampling = '
!    read(*,'(A)') n_resample_set
    if (if_append == 1) then
        call getarg(13, comm_arg_buff)
    else
        call getarg(12, comm_arg_buff)
    endif
    read(comm_arg_buff,*) n_resample_set
endif

!! warming up the random number generator
neg_seed = - n_pick(dummy_int)      ! seed should be negative integer
r_dummy = ran2(neg_seed)
do i = 1, 1000                      ! any number of iteration that is properly large
    r_dummy = ran2(neg_seed)
enddo

end subroutine initial_steps

!*******************************************************************************
! subroutine open_and_read
!*******************************************************************************
! Open the input and output files
! Read in the data file and store in the sample_array
!-------------------------------------------------------------------------------

subroutine open_and_read

use stuffs

implicit none

integer(kind=1), allocatable, dimension(:,:) :: genotype_array      ! (:,:)=(n_sample,n_SNP),holding genotype part of input data
real, allocatable, dimension(:,:) :: phenotype_array                  ! (:)=(n_sample,n_phenotype), holding phenotype part of the input data

integer, parameter :: n_sample_in_a_set = 400         ! each set has 400 samples

integer :: n_total_pair
integer :: i, j, k, kk
integer :: iend, istat	                                            ! file handling status; success if 0, fail otherwise
integer, dimension(n_SNP) :: dummy                                  ! dummy to skip interger type data
real,dimension(n_phenotype) :: QT_dummy                                 ! dummy to skip real type data
character(len=1) :: line_dummy                                              ! dummy to skip the line number of 'survMDR' data
integer :: n_advance                                   ! determine to the number of lines to reach the intended set

! open the appending file if requested
if (if_append == 1) then
    open(unit=appendfile_unit, file=appendingfile, status='old', action='write', access='append', iostat=istat)
    if (istat /= 0) stop '...Error opening appending file...' 
endif

! open the input data file
open(unit=infile_unit, file=infile, status='old', action='read', iostat=istat)
if (istat /= 0) stop '...Error opening input file...' 

! open the output file
open(unit=outfile_unit,file=outfile,status='new',action='write',iostat=istat)
if (istat /= 0) stop '...Error opening outfile ...'

! allocate genotype_array that will hold the input data
allocate (genotype_array(n_sample,n_SNP), stat=istat)
if (istat /= 0) stop '...Error allocating genotype_array ...'

! allocate phenotype_array that will hold the quantitative trait phenotype values
allocate (phenotype_array(n_sample,n_phenotype), stat=istat)
if (istat /= 0) stop '...Error allocating phenotype_array ...'

! reading phenotype part of the input data
if (if_header == 1) then
    read(infile_unit,*)                                                   ! skip the header line
endif    

! advance to the intended set; each set has 400 samples
n_advance = (n_set-1)*n_sample_in_a_set
do i = 1, n_advance
    read(infile_unit,*)
enddo

! set position_phenotype = 'R' only for the 'survMDR' data
position_phenotype = 'R'
if (position_phenotype == 'L') then
    do i = 1, n_sample
        read(1,*,iostat=iend) line_dummy, (phenotype_array(i,kk), kk=1, n_phenotype)                ! read in the phenotype part only, sample by sample
        if (iend > 0) stop 'Error reading input file ...'       ! iostat [> 0; error] [=0; no error, no eof] [<0; eof]
    enddo
elseif (position_phenotype == 'R') then
    do i = 1, n_sample
        read(1,*,iostat=iend) line_dummy, (dummy(k), k=1, n_SNP), (phenotype_array(i,kk), kk=1, n_phenotype)    ! read in the phenotype part only, sample by sample
        if (iend > 0) stop 'Error reading input file ...'                   ! iostat [> 0; error] [=0; no error, no eof] [<0; eof]
    enddo
else
    stop '... Wrong phenotype position specified ...'
endif

! reading genotype part of the input data
rewind(unit=infile_unit)                                       ! rewind the input file

if (if_header == 1) then
    read(1,*)                                                   ! skip the header line
endif    

n_advance = (n_set-1)*n_sample_in_a_set                         ! advance to the intended set; each set has 400 samples
do i = 1, n_advance
    read(infile_unit,*)
enddo

if (position_phenotype == 'L') then
    do i = 1, n_sample
        read(1,*,iostat=iend) line_dummy, (QT_dummy(kk), kk=1, n_phenotype), (genotype_array(i,k), k=1, n_SNP)   ! read in the genotype part only, sample by sample
        if (iend > 0) stop 'Error reading input file ...'                   ! iostat [> 0; error] [=0; no error, no eof] [<0; eof]
    enddo
elseif (position_phenotype == 'R') then
    do i = 1, n_sample
        read(1,*,iostat=iend) line_dummy, (genotype_array(i,k), k=1, n_SNP) ! read in the genotype part only, sample by sample
        if (iend > 0) stop 'Error reading input file ...'       ! iostat [> 0; error] [=0; no error, no eof] [<0; eof]
    enddo
else
    stop '... Wrong phenotype position specified ...'
endif

! save the original data array for permutation, by the ascending order of phenotype values
allocate (genotype_array_save(n_sample,n_SNP), stat=istat)      ! to save the original genotype data for permutation
if (istat /= 0) stop '...Error allocating genotype_array_save ...'
allocate (phenotype_array_save(n_sample,n_phenotype), stat=istat)           ! to save the original phenotype data for sorting
if (istat /= 0) stop '...Error allocating phenotype_array_save ...'
allocate (log_rho_table_save(n_sample,n_sample), stat=istat)                ! holding the index of the sorting result
if (istat /= 0) stop '...Error allocating log_rho_table_save ...'

phenotype_array_save = phenotype_array          ! copy to saving array
genotype_array_save = genotype_array            ! copy to saving array

! pre-process and save all possible distance bet multivarite phenotype vectors
call prepare_log_rho_table_save

! close data file
close(infile_unit)
write(*,*) '=== Data file reading is completed! ==='

! deallocate data holding arrays
deallocate(phenotype_array,stat=istat)
if (istat /= 0) stop 'Error deallocating phenotype_array ...'
deallocate(genotype_array,stat=istat)
if (istat /= 0) stop 'Error deallocating genotype_array ...'

! allocate and initialize top_ranks
if (n_top_write == 0) then
    if (n_dim == 1) then
        n_top_write = n_SNP
    else
        call get_n_total_pair(n_dim,n_SNP,n_total_pair)
        n_top_write = n_total_pair
    endif
endif

! allocate and initialize the estimation holding structure
allocate (kNN_top_ranks(n_top_write), stat=istat)
if (istat /= 0) stop '...Error allocating top_ranks ...'
do i = 1, n_top_write
    kNN_top_ranks(i)%SNP = 0                        ! SNP=0 means unassigned
    kNN_top_ranks(i)%H_X_kNN = 0.0                      ! H_X_kNN=0.0 means nothing there
    kNN_top_ranks(i)%H_XGY = 0.0                    ! H_XGY=0.0 means no info from G=y
    kNN_top_ranks(i)%kNN_MI = -LargeNum             ! if not properly assigned, a ridiculous value of -LargeNum will appear 
    kNN_top_ranks(i)%n_sample_valid = -LargeNum     ! if not properly assigned, a ridiculous value of -LargeNum will appear 
    kNN_top_ranks(i)%std_kNN_MI = -LargeNum         ! if not properly assigned, a ridiculous value of -LargeNum will appear 
    kNN_top_ranks(i)%p_value = 0.0                  ! initialize to zero because p-value is essentially a count
enddo    

! allocate and initialize the null dist. making array
if (if_permute == 1) then
    allocate(null_kNN_MI_by_max(n_resample_set), stat=istat)
    if (istat /= 0) stop '...Error allocating null_kNN_MI_by_max ...'
    do i = 1, n_resample_set
        null_kNN_MI_by_max(i) = -LargeNum           ! if not properly assigned, a ridiculous value of -LargeNum will appear
    enddo    
endif

end subroutine open_and_read

!*******************************************************************************
! subroutine cleanups
!*******************************************************************************

subroutine cleanups

use stuffs

implicit none

integer :: istat

! deallocate
deallocate(kNN_top_ranks,stat=istat)
if (istat /= 0) stop 'Error deallocating kNN_top_ranks ...'
deallocate(genotype_array_save,stat=istat)
if (istat /= 0) stop 'Error deallocating genotype_array_save ...'
deallocate(phenotype_array_save,stat=istat)
if (istat /= 0) stop 'Error deallocating phenotype_array_save ...'
deallocate(log_rho_table_save,stat=istat)
if (istat /= 0) stop 'Error deallocating log_rho_table_save ...'

! close output file
close(outfile_unit)
write(*,*) '=== Output file writing is completed! ==='
    
end subroutine cleanups

!*******************************************************************************
! subroutine get_n_total_pair
!*******************************************************************************
!     call get_n_total_pair(n_dim,n_SNP,n_total_pair)
!-------------------------------------------------------------------------------
subroutine get_n_total_pair(n_d,n_s,n_t)

use stuffs

implicit none

integer, intent(in) :: n_d  ! pairing dimension
integer, intent(in) :: n_s  ! number of SNP
integer, intent(out) :: n_t ! total number of SNP pair
integer :: i, j, k

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------
n_t = 0
if (n_d == 1) then
    n_t = n_s
elseif (n_d == 2) then
    do i = 1, n_s - 1
        do j = i+1, n_s
            n_t = n_t + 1
        enddo
    enddo
elseif (n_d == 3) then
    do i = 1, n_s - 2
        do j = i+1, n_s - 1
            do k = j+1, n_s
                n_t = n_t + 1
            enddo
        enddo
    enddo
endif

! print *,'Total number of pair is ',n_t

end subroutine get_n_total_pair

!*******************************************************************************
! subroutine prepare_log_rho_table_save
!*******************************************************************************

subroutine prepare_log_rho_table_save

use stuffs

implicit none

real :: temp
integer :: i, k, d

!initialize the rho table
log_rho_table_save = undefined_value     ! denote if -1.0, the element has not been defined

do i = 1, n_sample
    do k = 1, n_sample
        temp = 0.0
        if (i /= k) then
            do d = 1, n_phenotype
                temp = temp + (phenotype_array_save(i,d) - phenotype_array_save(k,d))**2.0
            enddo
        endif
        
        if (temp > 0.0) then
            log_rho_table_save(i,k) = log(temp)
        else
            log_rho_table_save(i,k) = 0.0
        endif
        
    enddo
enddo

end subroutine prepare_log_rho_table_save

!*******************************************************************************
! subroutine main_steps
!*******************************************************************************

subroutine main_steps

use stuffs
use omp_lib

implicit none

integer(kind=1), allocatable, dimension(:,:) :: genotype_array          ! (:,:)=(n_sample,n_SNP), holds genotype values
real, allocatable, dimension(:,:) :: phenotype_array          ! (:,:)=(n_sample,n_phenotype) holds phenotype values
real, allocatable, dimension(:,:) :: log_rho_table           ! (:,:)=(n_sample,n_sample) holds the all possible distances between the multivariate vectors

real :: temp_null_kNN_MI_max                        ! temporarily hold the max kNN_MI from each permutation
integer :: flag_resample                            ! tells if the resampled data is used or not (0/1) = (original/resampled)
integer, parameter :: iter_mon = 100                ! prompt the progress at every 100 iteration
integer :: istat
integer :: start_time, stop_time, count_rate        ! for system_clock
integer :: i, k, kk

! start timer
call SYSTEM_CLOCK(start_time, count_rate)

! allocate genotype, phenotype arrays and log_rho_table
allocate (genotype_array(n_sample,n_SNP), stat=istat)
if (istat /= 0) stop '...Error allocating genotype_array ...'
allocate (phenotype_array(n_sample,n_phenotype), stat=istat)
if (istat /= 0) stop '...Error allocating phenotype_array ...'
allocate (log_rho_table(n_sample,n_sample), stat=istat)
if (istat /= 0) stop '...Error allocating log_rho_table ...'

! prepare the data arrays and flag
flag_resample = 0                                   ! to signal to "get_kNN_MI" -> "estimate_kNN_MI" to carry out "sort_tops_kNN_MI"
genotype_array = genotype_array_save
phenotype_array = phenotype_array_save
log_rho_table = log_rho_table_save

! calculate H(P), which will be the same throughout the program
call cal_H_X_kNN(n_sample, n_phenotype)         
write(*,*) '... finished evaluating H(P) ...'

! calculate H(P|G) for every SNP pair and obtain MI, then sort
call get_kNN_MI(genotype_array, log_rho_table, flag_resample, temp_null_kNN_MI_max)   ! Here, temp_null_kNN_MI_max is only a dummy
write(*,*) '... finished evaluating kNN_MI for the original data'

! deallocate data holding arrays
deallocate(phenotype_array,stat=istat)              ! this array will be replaced by phenotype_array_permute in the following process
if (istat /= 0) stop 'Error deallocating phenotype_array ...'
deallocate(log_rho_table,stat=istat)
if (istat /= 0) stop 'Error deallocating log_rho_table ...'
deallocate(genotype_array,stat=istat)               ! this array will be replaced by permuted genotype array in the following process
if (istat /= 0) stop 'Error deallocating genotype_array ...'

! perform permuatation if requested
if (if_permute == 1) then

    flag_resample = 1
    
!    call omp_set_num_threads(3)                    ! using up all the available threads will degrade the performance
    call omp_set_dynamic(.TRUE.)
    
!$OMP parallel private(genotype_array,log_rho_table, temp_null_kNN_MI_max) shared(null_kNN_MI_by_max)

! allocate genotype_array and log_rho_table
    allocate (genotype_array(n_sample,n_SNP), stat=istat)
    if (istat /= 0) stop '...Error allocating genotype_array ...'
    allocate (log_rho_table(n_sample,n_sample), stat=istat)
    if (istat /= 0) stop '...Error allocating log_rho_table ...'

! prepare the data arrays and flag
    genotype_array = genotype_array_save
    log_rho_table = log_rho_table_save

!$OMP do
    do i = 1, n_resample_set
        call prepare_permuted_table(log_rho_table)
        call get_kNN_MI(genotype_array, log_rho_table, flag_resample, temp_null_kNN_MI_max)
        null_kNN_MI_by_max(i) = temp_null_kNN_MI_max    ! the max null kNN_MI in this iteration among all the possible SNP combinations
    
        if (mod(i,iter_mon)==0) then
            write(*,*) 'i_iter = ',i
        endif
    enddo
!$OMP end do

! deallocate data holding arrays
    deallocate(log_rho_table,stat=istat)
    if (istat /= 0) stop 'Error deallocating log_rho_table ...'
    deallocate(genotype_array,stat=istat)               ! this array will be replaced by permuted genotype array in the following process
    if (istat /= 0) stop 'Error deallocating genotype_array ...'

!$OMP end parallel

! get p-value by null distribution made of the array of null_kNN_MI_by_max
    call get_p_value_by_max             

endif   ! if (if_permute == 1)

! write out the result
if (if_permute == 1) then
    call write_tops_permuted
else
    call write_tops
endif

! end timer
call system_clock(stop_time, count_rate)
write(*,*) 'Elapsed time = ', float(stop_time - start_time) / float(count_rate), ' s'

end subroutine main_steps

!*******************************************************************************

subroutine cal_H_X_kNN(n,n_ph)
!===============================================================================
! Given Xi vector in d-space with n elements,
!	estimate the entropy, H(P), nonparametrically with kernel density estimation
!   refer to the LabNote p.3 dated 20190529
!===============================================================================

use stuffs

implicit none

integer, intent(in) :: n                            ! number of samples
integer, intent(in) :: n_ph                              ! number of phenotype

real :: digamma_ftn                    ! external function calculating digamma function
real :: term_1, term_2, term_3         ! terms consitituting H(P)
integer :: i, k                                     ! counters

! calculating term_1
term_1 = 0.0
do i = 1, n
    do k = 1, n     !log_rho_table(i,i) = 0, therefore would not contribute to the result anyway.
        if (i /= k) then
            term_1 = term_1 + log_rho_table_save(i,k)
        endif
    enddo
enddo
term_1 = term_1 * float(n_ph) / float(n*(n-1))

! calculation of term_2
term_2 = log(float(n-1))

! calculating term_3
term_3 = 0.0
do k = 1, n-1
    term_3 = term_3 + digamma_ftn(k)
enddo
term_3 = term_3 / float(n-1)

! finally, H(P)
!   log(Vd) term is omitted because it will be cancelled out in mutual information
H_X_kNN = term_1 + term_2 - term_3      

end subroutine cal_H_X_kNN

!*******************************************************************************
! real function digamma_ftn
!*******************************************************************************

real function digamma_ftn(k)

implicit none

integer, intent(in) :: k

real, parameter :: gammaEM = 0.577215664901532

integer :: i

if (k == 1) then
    digamma_ftn = -gammaEM
elseif (k >= 2) then
    digamma_ftn = 0.0
    do i = 1, k-1
        digamma_ftn = digamma_ftn + 1.0/float(i)
    enddo
    digamma_ftn = -gammaEM + digamma_ftn
else
    stop '... invalid k while calculating digamma function...'
endif

end function digamma_ftn

!*******************************************************************************

subroutine cal_H_XGY_kNN(Si,Gi,log_rho_table,n,d,H_XGY)
!===============================================================================
! estimate the conditional entropy, H(P|G)
!   refer to LabNote p.3. dated 20190529
!===============================================================================

use stuffs

implicit none

integer, intent(in) :: n                            ! number of samples w/ valid QT or genotype data
integer, dimension(n), intent(in) :: Si                ! Si: sample_Id_array corresponding to GT_array
integer, dimension(n), intent(in) :: Gi             ! Gi: GT_array containing genotype data
real, dimension(n_sample,n_sample), intent(in) :: log_rho_table
integer, intent(in) :: d                            ! genotype combination dimension
real, intent(out) :: H_XGY                          ! conditional entropy, H(P|G)

integer, allocatable, dimension(:) :: n_y           ! # of samples w/ genotype Y=y
integer, allocatable, dimension(:,:) :: Si_y           ! (:,:)=(0:n_gt_d-1,n_y_max); subset of Sample_Id_array, Si

integer, parameter :: n_y_limit = 3                 ! consider only the genotype at least w/ this number of data points

real :: digamma_ftn     ! digamma function

integer :: n_gt_d                                   ! number of genotype combination for d-dim: = 3**d
integer :: n_y_g                                    ! temporarily holding n_y(i_g)
integer :: istat                                    ! indicator for allocation
integer :: i_n, i_g, n_y_count                      ! counters
integer :: n_y_sum, n_y_max                         ! temporary summing variables

real :: term_1, term_2, term_3                      ! holding each term for H(P|G)

integer :: i, k                                     ! counters

! Determine the number of (joint) genotypes
n_gt_d = n_gt**d

! Allocate n_y, # of data points for each (joint) genotype
allocate(n_y(0:n_gt_d-1), stat=istat)
if (istat /= 0) stop '... Error allocating n_y...'

! Get n_y, # of data points for each (joint) genotype
n_y_max = 0

! determine n_y(i_g)
do i_g = 0, n_gt_d-1

    n_y_sum = 0

    do i_n = 1, n                                   ! count the number of samples w/ genotype Y=y or Y=i_g
        if (Gi(i_n) == i_g) then
            n_y_sum = n_y_sum + 1
        endif
    enddo
    
    n_y(i_g) = n_y_sum                              ! # of samples w/ Y=y ie. size of the subset Xi_y
    
    if (n_y(i_g) > n_y_max) then
        n_y_max = n_y(i_g)                          ! need this for allocation of Xi_y
    endif

enddo

! Allocate Si_y, holding array of valid sample Id
allocate(Si_y(0:n_gt_d-1,n_y_max), stat=istat)
if (istat /= 0) stop '...Error allocating Si_y...'

! Assign Si_y array value
do i_g = 0, n_gt_d-1
    n_y_count = 0
    do i_n = 1, n
        if (Gi(i_n) == i_g) then
            n_y_count = n_y_count + 1
            Si_y(i_g,n_y_count) = Si(i_n)           ! make the subset Si_y ie. Si w/ Y=y or Y=i_g
        endif
    enddo
    if (n_y_count /= n_y(i_g)) stop '...n_y_count conflict...'
enddo

! Calculating term_1
term_1 = 0.0
do i_g = 0, n_gt_d - 1
    n_y_g = n_y(i_g)
    if (n_y_g >= n_y_limit) then
        do i = 1, n_y_g
            do k = 1, n_y_g
                term_1 = term_1 + log_rho_table(Si_y(i_g,i),Si_y(i_g,k)) / float(n_y_g - 1)
            enddo
        enddo
    endif
enddo
term_1 = term_1 * float(n_phenotype) / float(n)

! Calculating term_2
term_2 = 0.0
do i_g = 0, n_gt_d - 1
    n_y_g = n_y(i_g)
    if (n_y_g >= n_y_limit) then
        term_2 = term_2 + float(n_y_g) * log(float(n_y_g - 1))
    endif
enddo
term_2 = term_2 / float(n)

! Calculating term_3
term_3 = 0.0
do i_g = 0, n_gt_d - 1
    n_y_g = n_y(i_g)
    if (n_y_g >= n_y_limit) then
        do k = 1, n_y_g - 1
            term_3 = term_3 + digamma_ftn(k) * float(n_y_g) / float(n_y_g - 1)
        enddo
    endif
enddo
term_3 = term_3 / float(n)

!cy
!print *, term_1, term_2, term_3

! Now calculating H(P|G)
! log(Vd) is omitted because it will be cancelled out in mutual information
H_XGY = term_1 + term_2 - term_3

!=== deallocation
deallocate(n_y,stat=istat)
deallocate(Si_y,stat=istat)
if (istat /= 0) stop 'Error deallocating either n_y, Si_y...'    

end subroutine cal_H_XGY_kNN

!*******************************************************************************
! subroutine get_kNN_MI
!*******************************************************************************

subroutine get_kNN_MI(genotype_array, log_rho_table, flag_resample, null_kNN_MI_max)

use stuffs

implicit none

integer(kind=1), dimension(n_sample,n_SNP), intent(in) :: genotype_array
real, dimension(n_sample,n_sample), intent(in) :: log_rho_table
integer, intent(in) :: flag_resample
real, intent(out) :: null_kNN_MI_max

real, parameter :: unlikely_negative = -10000.0

real :: old_kNN_MI_max, new_kNN_MI_max
integer :: j, k, i_pick
integer :: istat
integer :: i1, i2, i3
integer, dimension(n_dim) :: ii
integer, dimension(n_dim) :: SNP_pick
integer :: put_gt_combi
integer :: gt_count
integer :: n_sample_valid
integer, dimension(n_dim) :: picked_g
integer :: g1, g2, g3
integer :: n_total_pair
integer, dimension(n_sample) :: GT_array            ! GT=genotype as in multifactor
integer, dimension(n_sample) :: sample_Id_array       ! record valid sample Id corresponding to GT_array

! initialize
old_kNN_MI_max = unlikely_negative

! get the total number of pair
call get_n_total_pair(n_dim,n_SNP,n_total_pair)

!-------------------------------------------------------------------------------
! n_dim = 1, 2, 3 for now
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! n_dim = 1
!-------------------------------------------------------------------------------
if_D1: if (n_dim == 1) then
    
! pick all the possible SNP     
    do_pick_D1: do i1 = 1, n_SNP
        ii(1) = i1	! TT; NOT elegant...
	    do k = 1, n_dim
		    SNP_pick(k) = ii(k)
        enddo
        
!-------------------------------------------------------------------------------
! Core calculation
!-------------------------------------------------------------------------------
! Put the data in the PT_array for a ptcular SNP combi, GT_array
        n_sample_valid = 0                          ! initialize the number of valid samples for this SNP pick

        do j = 1, n_sample                          ! for every sample for this picked SNP
            do i_pick = 1, n_dim
                picked_g(i_pick) = genotype_array(j,SNP_pick(i_pick))    ! data format: phenotype genotype (ascending order in phenotype)
            enddo
            
            if (picked_g(1) >= 0) then              ! valid genotype = 0, 1, 2 (invalid genotype such as missing data is encode as -9)
                n_sample_valid = n_sample_valid + 1 ! increment of number of valid sample
                gt_count = 0                        ! initializing the gt_combi assignment
                put_gt_combi_D1: do g1 = 0, n_gt-1
                    if (picked_g(1) == g1) then
                        put_gt_combi = gt_count     ! proper number of gt_combi found
                        exit put_gt_combi_D1        ! no need to find further
                    else
                        gt_count = gt_count + 1     ! examine the next number for proper gt_combi assignment
                    endif                        
                enddo put_gt_combi_D1                
                        
                sample_Id_array(n_sample_valid) = j     ! put valid sample Id corresponding to GT_array
                GT_array(n_sample_valid) = put_gt_combi             ! put_gt_combi = 0,1,2
            endif
        enddo

! estimate the GxP association for this ptcular SNP combi    
        call estimate_kNN_MI(flag_resample,n_sample_valid,SNP_pick,sample_Id_array,GT_array,log_rho_table,old_kNN_MI_max,new_kNN_MI_max)
        old_kNN_MI_max = new_kNN_MI_max
        
    enddo do_pick_D1
    
endif if_D1    

!-------------------------------------------------------------------------------
! n_dim = 2
!-------------------------------------------------------------------------------
if_D2: if (n_dim == 2) then

! pick all the possible SNP     
    do_pick_D2: do i1 = 1, n_SNP-1
    	ii(1) = i1	! TT
        do i2 = i1+1, n_SNP
		    ii(2) = i2	! TT
            do k = 1, n_dim
			    SNP_pick(k) = ii(k)
            enddo

!-------------------------------------------------------------------------------
! Core calculation
!-------------------------------------------------------------------------------
! Put the data in the PT_array for a ptcular SNP combi, GT_array
            n_sample_valid = 0                      ! initialize the number of valid samples for this SNP pick

            do j = 1, n_sample                      ! for every sample for this picked SNP
                do i_pick = 1, n_dim
                    picked_g(i_pick) = genotype_array(j,SNP_pick(i_pick))   ! data format: phenotype genotype (ascending order in phenotype)
                enddo

                if ((picked_g(1) >= 0).AND.(picked_g(2) >=0)) then          ! valid genotype = 0, 1, 2 (invalid genotype such as missing data is encode as -9)
                    n_sample_valid = n_sample_valid + 1                     ! increment of number of valid sample
                    gt_count = 0                                            ! initializing the gt_combi assignment
                    put_gt_combi_D2: do g1 = 0, n_gt-1
                        do g2 = 0, n_gt-1
                            if ((picked_g(1) == g1).AND.(picked_g(2) == g2)) then
                                put_gt_combi = gt_count                     ! proper number of gt_combi found
                                exit put_gt_combi_D2                        ! no need to find further
                            else
                                gt_count = gt_count + 1                     ! examine the next number for proper gt_combi assignment
                            endif            
                        enddo
                    enddo put_gt_combi_D2
                        
                    sample_Id_array(n_sample_valid) = j             ! put valid sample Id corresponding to GT_array
                    GT_array(n_sample_valid) = put_gt_combi                 ! put_gt_combi = 0,1,2,3,4,5,6,7,8

                endif
            enddo

! estimate the GxP association for this ptcular SNP combi    
            call estimate_kNN_MI(flag_resample,n_sample_valid,SNP_pick,sample_Id_array,GT_array,log_rho_table,old_kNN_MI_max,new_kNN_MI_max)
            old_kNN_MI_max = new_kNN_MI_max

        enddo
    enddo do_pick_D2
    
endif if_D2

!-------------------------------------------------------------------------------
! n_dim = 3
!-------------------------------------------------------------------------------
if_D3: if (n_dim == 3) then

! pick all the possible SNP     
    do_pick_D3: do i1 = 1, n_SNP-2
    	ii(1) = i1	! TT
        do i2 = i1+1, n_SNP-1
		    ii(2) = i2	! TT
            do i3 = i2+1, n_SNP
		        ii(3) = i3	! TT
                do k = 1, n_dim
			        SNP_pick(k) = ii(k)
                enddo

!-------------------------------------------------------------------------------
! Core calculation
!-------------------------------------------------------------------------------
! Put the data in the PT_array for a ptcular SNP combi, GT_array
                n_sample_valid = 0                  ! initialize the number of valid samples for this SNP pick

                do j = 1, n_sample                  ! for every sample for this picked SNP
                    do i_pick = 1, n_dim
                        picked_g(i_pick) = genotype_array(j,SNP_pick(i_pick))   ! data format: phenotype genotype (ascending order in phenotype)
                    enddo

                    if ((picked_g(1) >= 0).AND.(picked_g(2) >=0).AND.(picked_g(3) >=0)) then  ! valid genotype = 0, 1, 2 (invalid genotype such as missing data is encode as -9)
                        n_sample_valid = n_sample_valid + 1                     ! increment of number of valid sample
                        gt_count = 0                                            ! initializing the gt_combi assignment
                        put_gt_combi_D3: do g1 = 0, n_gt-1
                            do g2 = 0, n_gt-1
                                do g3 = 0, n_gt-1
                                    if ((picked_g(1) == g1).AND.(picked_g(2) == g2).AND.(picked_g(3) == g3)) then
                                        put_gt_combi = gt_count                 ! proper number of gt_combi found
                                        exit put_gt_combi_D3                    ! no need to find further
                                    else
                                        gt_count = gt_count + 1                 ! examine the next number for proper gt_combi assignment
                                    endif
                                enddo
                            enddo
                        enddo put_gt_combi_D3
                        
                        sample_Id_array(n_sample_valid) = j         ! put valid sample Id corresponding to GT_array
                        GT_array(n_sample_valid) = put_gt_combi                 ! put_gt_combi = 0,1,2 ..., 3x3x3-1

                    endif
                enddo

! estimate the GxP association for this ptcular SNP combi    
                call estimate_kNN_MI(flag_resample,n_sample_valid,SNP_pick,sample_Id_array,GT_array,log_rho_table,old_kNN_MI_max,new_kNN_MI_max)
                old_kNN_MI_max = new_kNN_MI_max

            enddo
        enddo
    enddo do_pick_D3
    
endif if_D3

! Got the null_kNN_MI_max from "estimate_kNN_MI" for this permutated dataset.
! If this came from the process with if_permute == 0, 
! then the info was given to "sort_tops_kNN_MI" at the end of "estimate_kNN_MI", and the following line means nothing.
null_kNN_MI_max = old_kNN_MI_max    ! immediately after each 'estimate_kNN_MI', new_kNN_MI_max is put into old_kNN_MI_max

end subroutine get_kNN_MI

!*******************************************************************************
! subroutine estimate_kNN_MI
!*******************************************************************************

subroutine estimate_kNN_MI(flag_resample,n_sample_valid,SNP_pick,sample_Id_array,GT_array,log_rho_table,old_kNN_MI_max,new_kNN_MI_max)

use stuffs

implicit none

integer, intent(in) :: flag_resample
integer, intent(in) :: n_sample_valid
integer, dimension(n_dim), intent(in) :: SNP_pick
integer, dimension(n_sample_valid), intent(in) :: sample_Id_array
integer, dimension(n_sample_valid), intent(in) :: GT_array
real, dimension(n_sample,n_sample), intent(in) :: log_rho_table
real, intent(in) :: old_kNN_MI_max
real, intent(out) :: new_kNN_MI_max

real :: H_XGY                                       ! H(P|G)
real :: kNN_MI                                      ! MI = H(P) - H(P|G) = H_X - H_XGY
integer :: i, j, k
integer :: istat

! Calculate H(P|G) = H(X|Y)
call cal_H_XGY_kNN(sample_Id_array, GT_array, log_rho_table, n_sample_valid, n_dim, H_XGY)

! Calculate kNN_MI = MI 
kNN_MI = H_X_kNN - H_XGY                                ! H_X_kNN is defined globally, calculated already using "cal_H_X_kNN"

! sort/store(flag_resample == 0) or add another kNN_MI to null dist(flag_resample == 1)
if (flag_resample == 0) then                        ! no permutation
    call sort_tops_kNN_MI(n_top_write,n_dim,SNP_pick,H_XGY,kNN_MI,n_sample_valid)
else                                                ! do permutation
    if (kNN_MI > old_kNN_MI_max) then
        new_kNN_MI_max = kNN_MI                     ! update the candidate for the null_kNN_MI_max
    else
        new_kNN_MI_max = old_kNN_MI_max
    endif        
endif

end subroutine estimate_kNN_MI

!*******************************************************************************
! subroutine sort_tops_kNN_MI
!*******************************************************************************
subroutine sort_tops_kNN_MI(n_write,n,SNP,H_XGY,kNN_MI,n_sample_valid)

use stuffs

implicit none

integer, intent(in) :: n_write                                  ! == n_top_write; number of top rankers to write
integer, intent(in) :: n                                        ! == n_dim; number of SNPs involved in the association
integer, dimension(n), intent(in) :: SNP                        ! store the SNP numbers involved
real, intent(in) :: H_XGY, kNN_MI
integer, intent(in) :: n_sample_valid
integer :: i, j                                                 ! counters

!-------------------------------------------------------------------------------
! find out the rank and squeeze into the proper place
!-------------------------------------------------------------------------------
do i = 1, n_write
    if (kNN_MI >= kNN_top_ranks(i)%kNN_MI) then    ! sorting criteria of kNN_MI used
        do j = n_write, i+1, -1   ! demote by one to make place for the newly found top ranker
            kNN_top_ranks(j) = kNN_top_ranks(j-1)
        enddo
! put the new top ranker in top_ranks record; other components are defined after permutation
        kNN_top_ranks(i)%SNP = SNP
        kNN_top_ranks(i)%H_X_kNN = H_X_kNN
        kNN_top_ranks(i)%H_XGY = H_XGY
        kNN_top_ranks(i)%kNN_MI = kNN_MI
        kNN_top_ranks(i)%n_sample_valid = n_sample_valid
        
        exit    ! vital; stop comparing and return control to the calling process
        
    endif
enddo    

end subroutine sort_tops_kNN_MI

!*******************************************************************************
! subroutine write_tops
!*******************************************************************************
! 
!-------------------------------------------------------------------------------
subroutine write_tops

use stuffs

implicit none

integer :: i    ! counter

!-------------------------------------------------------------------------------
! use stuffs: n_dim, top_ranks, n_top_write
!-------------------------------------------------------------------------------
if (n_dim == 1) then
    do i = 1, n_top_write
        write(9,'(2I15,3F30.15,I15)') i,                     &
                            kNN_top_ranks(i)%SNP(1),        &
                            kNN_top_ranks(i)%H_X_kNN,           &
                            kNN_top_ranks(i)%H_XGY,         &
                            kNN_top_ranks(i)%kNN_MI,        &
                            kNN_top_ranks(i)%n_sample_valid
    enddo
    
    ! write appending file
    if (if_append == 1) then
        ! append 1st ranked result
        write(11,'(3I15,3F30.15,I15)') n_sample, n_SNP,      &  
                            kNN_top_ranks(1)%SNP(1),        &
                            kNN_top_ranks(1)%H_X_kNN,           &
                            kNN_top_ranks(1)%H_XGY,         &
                            kNN_top_ranks(1)%kNN_MI,        &
                            kNN_top_ranks(1)%n_sample_valid
    endif

elseif (n_dim == 2) then
    do i = 1, n_top_write
        write(9,'(3I15,3F30.15,I15)') i,                     &
                            kNN_top_ranks(i)%SNP(1),        &
                            kNN_top_ranks(i)%SNP(2),        &
                            kNN_top_ranks(i)%H_X_kNN,           &
                            kNN_top_ranks(i)%H_XGY,         &
                            kNN_top_ranks(i)%kNN_MI,        &
                            kNN_top_ranks(i)%n_sample_valid
    enddo
    
    ! write appending file
    if (if_append == 1) then
        ! append 1st ranked result
        write(11,'(4I15,3F30.15,I15)') n_sample, n_SNP,      & 
                            kNN_top_ranks(1)%SNP(1),        &
                            kNN_top_ranks(1)%SNP(2),        &
                            kNN_top_ranks(1)%H_X_kNN,           &
                            kNN_top_ranks(1)%H_XGY,         &
                            kNN_top_ranks(1)%kNN_MI,        &
                            kNN_top_ranks(1)%n_sample_valid
    endif

elseif (n_dim == 3) then
    do i = 1, n_top_write
        write(9,'(4I15,3F30.15,I15)') i,                     &
                            kNN_top_ranks(i)%SNP(1),        &
                            kNN_top_ranks(i)%SNP(2),        &
                            kNN_top_ranks(i)%SNP(3),        &
                            kNN_top_ranks(i)%H_X_kNN,           &
                            kNN_top_ranks(i)%H_XGY,         &
                            kNN_top_ranks(i)%kNN_MI,        &
                            kNN_top_ranks(i)%n_sample_valid
    enddo
    
    ! write appending file
    if (if_append == 1) then
        ! append 1st ranked result
        write(11,'(5I15,3F30.15,I15)') n_sample, n_SNP,      &
                            kNN_top_ranks(1)%SNP(1),        &
                            kNN_top_ranks(1)%SNP(2),        &
                            kNN_top_ranks(1)%SNP(3),        &
                            kNN_top_ranks(1)%H_X_kNN,           &
                            kNN_top_ranks(1)%H_XGY,         &
                            kNN_top_ranks(1)%kNN_MI,        &
                            kNN_top_ranks(1)%n_sample_valid
    endif
    
endif

end subroutine write_tops

!*******************************************************************************
! subroutine prepare_permuted_table
!*******************************************************************************

subroutine prepare_permuted_table(log_rho_table)

use stuffs

implicit none

real, dimension(n_sample,n_sample), intent(out) :: log_rho_table

real :: ran2       ! external function; generate random number bet 0 and 1 exclusive of both ends
integer :: i, k
integer :: r1, r2

! permute the log_rho_table using the save original data in log_rho_table_sav
do i = 1, n_sample
    do k = 1, n_sample
        r1 = int(ran2(neg_seed)*float(n_sample)) + 1            ! generate random integer bet 1 and n_sample
        r2 = int(ran2(neg_seed)*float(n_sample)) + 1
        log_rho_table(i,k) = log_rho_table_save(r1,r2)      ! permutation w/ replacement
    enddo
enddo

end subroutine prepare_permuted_table

!*******************************************************************************
! subroutine get_p_value_by_max
!*******************************************************************************

!-------------------------------------------------------------------------------
subroutine get_p_value_by_max

use stuffs

implicit none

integer :: i, k
real :: temp_sum_1, temp_sum_2, temp_dev    ! used for corrected two-pass algorithm
real :: null_mean, null_stdev

! calculate the mean of the null distribution
null_mean = 0.0
do k = 1, n_resample_set
    null_mean = null_mean + null_kNN_MI_by_max(k)
enddo
null_mean = null_mean / float(n_resample_set)
null_kNN_MI_mean_by_max = null_mean
   
! calculate the standard deviation of the null dist.
! use the corrected two-pass algorithm to enhance the accuracy
temp_sum_1 = 0.0
temp_sum_2 = 0.0
do k = 1, n_resample_set
    temp_sum_1 = temp_sum_1 + (null_kNN_MI_by_max(k) - null_mean)**2
    temp_sum_2 = temp_sum_2 + (null_kNN_MI_by_max(k) - null_mean)
enddo
temp_sum_2 = temp_sum_2*temp_sum_2 / float(n_resample_set)
temp_dev = sqrt((temp_sum_1 - temp_sum_2) / float(n_resample_set - 1))
null_stdev = temp_dev
null_kNN_MI_stdev_by_max = null_stdev

! For the top_ranks
do i = 1, n_top_write
    
! get std_kNN_MI
    kNN_top_ranks(i)%std_kNN_MI = (kNN_top_ranks(i)%kNN_MI - null_mean) / null_stdev

! calculate the significance(p-value)
! count the number of events from null distribution that exceed the observed value from the top sorted values
    do k = 1, n_resample_set
        if (null_kNN_MI_by_max(k) >= kNN_top_ranks(i)%kNN_MI) then
            kNN_top_ranks(i)%p_value = kNN_top_ranks(i)%p_value + 1.0
        endif
    enddo
    kNN_top_ranks(i)%p_value = kNN_top_ranks(i)%p_value + 1.0   ! if all possible permutation is performed, at least one(i.e. self) should satisfy ">=" condition
    kNN_top_ranks(i)%p_value = kNN_top_ranks(i)%p_value / float(n_resample_set)

enddo

end subroutine get_p_value_by_max

!*******************************************************************************
! subroutine write_tops_permuted
!*******************************************************************************
! 
!-------------------------------------------------------------------------------
subroutine write_tops_permuted

use stuffs

implicit none

integer :: i    ! counter

!-------------------------------------------------------------------------------
! use stuffs: n_dim, top_ranks, n_top_write
!-------------------------------------------------------------------------------
if (n_dim == 1) then
    do i = 1, n_top_write
        write(9,'(2I15,7F30.15,I15)') i,                    &
                            kNN_top_ranks(i)%SNP(1),        &
                            kNN_top_ranks(i)%H_X_kNN,           &
                            kNN_top_ranks(i)%H_XGY,         &
                            kNN_top_ranks(i)%kNN_MI,        &
                            null_kNN_MI_mean_by_max,        &
                            null_kNN_MI_stdev_by_max,       &
                            kNN_top_ranks(i)%std_kNN_MI,    &
                            kNN_top_ranks(i)%p_value,       &
                            kNN_top_ranks(i)%n_sample_valid
    enddo
    
    ! write appending file
    if (if_append == 1) then
        ! append 1st ranked result
        write(11,'(3I15,7F30.15,I15)') n_sample, n_SNP,     &
                            kNN_top_ranks(1)%SNP(1),        &
                            kNN_top_ranks(1)%H_X_kNN,           &
                            kNN_top_ranks(1)%H_XGY,         &
                            kNN_top_ranks(1)%kNN_MI,        &
                            null_kNN_MI_mean_by_max,        &  
                            null_kNN_MI_stdev_by_max,       &
                            kNN_top_ranks(1)%std_kNN_MI,    &
                            kNN_top_ranks(1)%p_value,       &
                            kNN_top_ranks(1)%n_sample_valid
    endif

elseif (n_dim == 2) then
    do i = 1, n_top_write
        write(9,'(3I15,7F30.15,I15)') i,                    &
                            kNN_top_ranks(i)%SNP(1),        &
                            kNN_top_ranks(i)%SNP(2),        &
                            kNN_top_ranks(i)%H_X_kNN,           &
                            kNN_top_ranks(i)%H_XGY,         &
                            kNN_top_ranks(i)%kNN_MI,        &
                            null_kNN_MI_mean_by_max,        &
                            null_kNN_MI_stdev_by_max,       &
                            kNN_top_ranks(i)%std_kNN_MI,    &
                            kNN_top_ranks(i)%p_value,       &
                            kNN_top_ranks(i)%n_sample_valid
    enddo
    
    ! write appending file
    if (if_append == 1) then
        ! append 1st ranked result
        write(11,'(4I15,7F30.15,I15)') n_sample, n_SNP,     & 
                            kNN_top_ranks(1)%SNP(1),        &
                            kNN_top_ranks(1)%SNP(2),        &
                            kNN_top_ranks(1)%H_X_kNN,           &
                            kNN_top_ranks(1)%H_XGY,         &
                            kNN_top_ranks(1)%kNN_MI,        &
                            null_kNN_MI_mean_by_max,        &
                            null_kNN_MI_stdev_by_max,       &
                            kNN_top_ranks(1)%std_kNN_MI,    &
                            kNN_top_ranks(1)%p_value,       &
                            kNN_top_ranks(1)%n_sample_valid
    endif

elseif (n_dim == 3) then
    do i = 1, n_top_write
        write(9,'(4I15,7F30.15,I15)') i,                    &
                            kNN_top_ranks(i)%SNP(1),        &
                            kNN_top_ranks(i)%SNP(2),        &
                            kNN_top_ranks(i)%SNP(3),        &
                            kNN_top_ranks(i)%H_X_kNN,           &
                            kNN_top_ranks(i)%H_XGY,         &
                            kNN_top_ranks(i)%kNN_MI,        &
                            null_kNN_MI_mean_by_max,        &
                            null_kNN_MI_stdev_by_max,       &
                            kNN_top_ranks(i)%std_kNN_MI,    &
                            kNN_top_ranks(i)%p_value,       &
                            kNN_top_ranks(i)%n_sample_valid
    enddo
    
    ! write appending file
    if (if_append == 1) then
        ! append 1st ranked result
        write(11,'(5I15,7F30.15,I15)') n_sample, n_SNP,     &
                            kNN_top_ranks(1)%SNP(1),        &
                            kNN_top_ranks(1)%SNP(2),        &
                            kNN_top_ranks(1)%SNP(3),        &
                            kNN_top_ranks(1)%H_X_kNN,           &
                            kNN_top_ranks(1)%H_XGY,         &
                            kNN_top_ranks(1)%kNN_MI,        &
                            null_kNN_MI_mean_by_max,        &
                            null_kNN_MI_stdev_by_max,       &
                            kNN_top_ranks(1)%std_kNN_MI,    &
                            kNN_top_ranks(1)%p_value,       &
                            kNN_top_ranks(1)%n_sample_valid
    endif
    
endif

end subroutine write_tops_permuted

!*******************************************************************************
! function ran2
!       input: idum - any negative integer
!       returns: uniform deviates between 0 and 1, exclusive of both end points
!*******************************************************************************

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,            &
                 IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,    &
                 NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software *!%91a9I.

!*******************************************************************************
! subroutine indexx
!*******************************************************************************
!     call indexx(n,arr,indx) : arr(indx(j)) should be in ascending order
!-------------------------------------------------------------------------------

      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END subroutine indexx
      
! C  (C) Copr. 1986-92 Numerical Recipes Software *!%91a9I.

!*******************************************************************************