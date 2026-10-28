! mod_perturb_param.f90
! =====================================
!  Module for Parameter Perturbations
! ======================================

module mod_perturb_param
  use common_pom_var
  implicit none
  private

  public :: init_perturb_alp2, update_perturb_alp2

contains

  ! ===================================================================
  !  Initialization: Set initial parameter fields (Called at initialize.f90)
  ! ===================================================================
  subroutine init_perturb_alp2(im, jm, nens, iens)
    integer, intent(in) :: im, jm
    integer, intent(in) :: nens, iens

    real(kind = r_size) :: assigned_val
    integer :: iyr, imon, iday

    real(kind = r_size), parameter :: alp2_default = 0.53d0

    integer :: seed_size
    integer, allocatable :: seed_array(:)

    if (.not. allocated(alp2)) allocate(alp2(im, jm))

    ! --- 乱数シードの初期化（メンバ一意） ---
    call random_seed(size=seed_size)
    allocate(seed_array(seed_size))

    call julian_ymd(int(julday_start+time), iyr, imon, iday)
    seed_array(:) = iens * 1000000 + iyr * 100 + imon

    call random_seed(put=seed_array)
    deallocate(seed_array)

    if (i_init_alp2 == 0) then
        ! Pattern 0: Common default value
        alp2(:,:) = alp2_default

    else if (i_init_alp2 == 1) then
        ! Pattern 1: Linspace (Equally spaced spread)
        if (nens > 1) then
           assigned_val = alp2_min + &
                (alp2_max - alp2_min) * real(iens - 1, kind=r_size) / real(nens - 1, kind=r_size)
        else
           assigned_val = alp2_default
        end if
        alp2(:,:) = assigned_val

    else
        alp2(:,:) = alp2_default
    end if

  end subroutine init_perturb_alp2

  subroutine init_perturb_horcon(im, jm, nens, iens)
    integer, intent(in) :: im, jm
    integer, intent(in) :: nens, iens

    real(kind = r_size) :: assigned_val
    integer :: iyr, imon, iday

    real(kind = r_size), parameter :: horcon_default = 0.53d0

    integer :: seed_size
    integer, allocatable :: seed_array(:)

    if (.not. allocated(alp2)) allocate(horcon(im, jm))

    ! --- 乱数シードの初期化（メンバ一意） ---
    call random_seed(size=seed_size)
    allocate(seed_array(seed_size))

    call julian_ymd(int(julday_start+time), iyr, imon, iday)
    seed_array(:) = iens * 1000000 + iyr * 100 + imon

    call random_seed(put=seed_array)
    deallocate(seed_array)

    if (i_init_horcon == 0) then
        ! Pattern 0: Common default value
        horcon(:,:) = horcon_default

    else if (i_init_horcon == 1) then
        ! Pattern 1: Linspace (Equally spaced spread)
        if (nens > 1) then
           assigned_val = horcon_min + &
                (horcon_max - horcon_min) * real(iens - 1, kind=r_size) / real(nens - 1, kind=r_size)
        else
           assigned_val = horcon_default
        end if
        horcon(:,:) = assigned_val

    else
        horcon(:,:) = horcon_default
    end if

  end subroutine init_perturb_horcon

  ! ===================================================================
  !  Update: Time evolution and Constraints
  ! ===================================================================
  subroutine update_perturb_alp2(im, jm)
    integer, intent(in) :: im, jm

    integer :: i, j
    real(kind = r_size) :: rand_val
    real(kind = r_size), allocatable :: rand_field(:,:)

    ! --- 1. Add time-varying noise (前ステップの値を維持してノイズ加算) ---
    if (i_update_alp2 == 1) then
        ! Mode 1: 1D uniform random noise (領域一律)
        call random_number(rand_val)
        alp2(:,:) = alp2(:,:) + (rand_val - 0.5d0) * 2.0d0 * pert_amp_alp2

    else if (i_update_alp2 == 2) then
        ! Mode 2: 2D spatial random noise
        allocate(rand_field(im, jm))
        call random_number(rand_field)
        do j = 1, jm
           do i = 1, im
              alp2(i,j) = alp2(i,j) + (rand_field(i,j) - 0.5d0) * 2.0d0 * pert_amp_alp2
           end do
        end do
        deallocate(rand_field)
    end if

    ! --- 2. Reflecting Boundary & Physical Constraints ---
    do j = 1, jm
       do i = 1, im
          ! フェールセーフ（安全クランプ）
          alp2(i,j) = max(alp2_min, min(alp2_max, alp2(i,j)))
       end do
    end do

  end subroutine update_perturb_alp2

  subroutine update_perturb_horcon(im, jm)
    integer, intent(in) :: im, jm

    integer :: i, j
    real(kind = r_size) :: rand_val
    real(kind = r_size), allocatable :: rand_field(:,:)

    ! --- 1. Add time-varying noise (前ステップの値を維持してノイズ加算) ---
    if (i_update_horcon == 1) then
        ! Mode 1: 1D uniform random noise (領域一律)
        call random_number(rand_val)
        horcon(:,:) = horcon(:,:) + (rand_val - 0.5d0) * 2.0d0 * pert_amp_horcon

    else if (i_update_horcon == 2) then
        ! Mode 2: 2D spatial random noise
        allocate(rand_field(im, jm))
        call random_number(rand_field)
        do j = 1, jm
           do i = 1, im
              horcon(i,j) = horcon(i,j) + (rand_field(i,j) - 0.5d0) * 2.0d0 * pert_amp_horcon
           end do
        end do
        deallocate(rand_field)
    end if

    ! --- 2. Reflecting Boundary & Physical Constraints ---
    do j = 1, jm
       do i = 1, im
          ! フェールセーフ（安全クランプ）
          horcon(i,j) = max(horcon_min, min(horcon_max, horcon(i,j)))
       end do
    end do

  end subroutine update_perturb_horcon
end module mod_perturb_param
