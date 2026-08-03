! mod_perturb_param.f90
! =====================================
!  Module for Parameter Perturbations
! ======================================

module mod_perturb_param
  use common_pom_var
!  use MYNNF_lev25_2012
  implicit none
!  private

!  public :: init_perturb_alp2, update_perturb_alp2


contains

  ! ===================================================================
  !  Initialization: Set initial parameter fields (Called at initialize.f90 once)
  ! ===================================================================
  subroutine init_perturb_alp2(im, jm, irestart, nens, iens)
    integer, intent(in) :: im, jm
    logical, intent(in) :: irestart
    integer, intent(in) :: nens, iens


    real(kind = r_size) :: assigned_val

    integer :: iyr, imon, iday   !compute with julday_ymd

    real(kind = r_size), parameter :: alp2_default = 0.53d0

    ! シード用変数の宣言を追加
    integer :: seed_size
    integer, allocatable :: seed_array(:)

    if (.not. allocated(alp2)) allocate(alp2(im, jm))


    ! --- 乱数シードの初期化（ここで1回だけ実行！） ---
    call random_seed(size=seed_size)
    allocate(seed_array(seed_size))
    
    ! 年、月、アンサンブルメンバ番号から一意なシードを生成
    ! 例: iens=3, iyr=2012, imon=8 の場合 -> 3000000 + 201200 + 8 = 3201208
    call julian_ymd(int(julday_start+time),iyr,imon,iday)
    seed_array(:) = iens * 1000000 + iyr * 100 + imon
    
    call random_seed(put=seed_array)
    deallocate(seed_array)

    ! Do nothing for restart runs to preserve LETKF-updated values
    if (irestart) return

    ! Set initial values for cold start
    if (i_init_alp2 == 1) then
       ! Pattern 1: Linspace (Equally spaced spread)
       if (nens > 1) then
          assigned_val = alp2_min + &
               (alp2_max - alp2_min) * real(iens - 1, kind=r_size) / real(nens - 1, kind=r_size)
       else
          assigned_val = alp2_default
       end if
       alp2(:,:) = assigned_val

    else if (i_init_alp2 == 0) then
       ! Pattern 0: Common default value for all members
       alp2(:,:) = alp2_default
    else
       alp2(:,:) = alp2_default 
    end if

  end subroutine init_perturb_alp2

  ! ===================================================================
  !  Update: Time evolution and Physical Constraints
  ! ===================================================================
  subroutine update_perturb_alp2(im, jm, nens)
    integer, intent(in) :: im, jm
    integer, intent(in) :: nens

    integer :: i, j
    real(kind = r_size) :: rand_val
    real(kind = r_size), allocatable :: rand_field(:,:)
    
    real(kind = r_size), intent(in) :: pert_amp_alp2  ! perturbation amplitude


    ! --- 1. Add time-varying noise (Modes 1 and 2) ---
    ! 毎ステップ、単純に random_number を引くだけで済むように！
    if (i_update_alp2 == 1 .or. i_update_alp2 == 2) then
       
       if (i_update_alp2 == 2) then
          ! Mode 2: 2D spatial random noise
          allocate(rand_field(im, jm))
          call random_number(rand_field)
          do j = 1, jm
             do i = 1, im
                alp2(i,j) = alp2(i,j) + (rand_field(i,j) - 0.5d0) * 2.0d0 * pert_amp_alp2
             end do
          end do
          deallocate(rand_field)
       else
          ! Mode 1: 1D uniform random noise
          call random_number(rand_val)
          alp2(:,:) = alp2(:,:) + (rand_val - 0.5d0) * 2.0d0 * pert_amp_alp2
       end if
    end if

    ! --- 2. Reflecting Boundary ---
    ! Reflect values inward if they exceed boundaries to maintain natural variance
    do j = 1, jm
       do i = 1, im
          ! Reflection at the upper boundary
          if (alp2(i,j) > alp2_max) then
             alp2(i,j) = 2.0d0 * alp2_max - alp2(i,j)
          end if
          
          ! Reflection at the lower boundary
          if (alp2(i,j) < alp2_min) then
             alp2(i,j) = 2.0d0 * alp2_min - alp2(i,j)
          end if
          
          ! Failsafe (prevents values from escaping bounds due to extreme cases)
          alp2(i,j) = max(alp2_min, min(alp2_max, alp2(i,j)))
       end do
    end do

  end subroutine update_perturb_alp2

end module mod_perturb_param
