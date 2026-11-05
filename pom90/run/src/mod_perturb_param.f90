module mod_perturb_param
  use common_pom_var, only: r_size
  implicit none
  private

  public :: init_perturb_2d, update_perturb_2d

contains

  subroutine init_perturb_2d(field, im, jm, nens, iens, &
                             default_val, min_val, max_val, init_mode)
    integer, intent(in) :: im, jm, nens, iens, init_mode
    real(kind=r_size), intent(in) :: default_val, min_val, max_val
    real(kind=r_size), intent(inout), allocatable :: field(:,:)

    real(kind=r_size) :: assigned_val

    if (.not. allocated(field)) allocate(field(im, jm))

    if (init_mode == 1) then
       if (nens > 1) then
          assigned_val = min_val + (max_val - min_val) * &
               real(iens - 1, kind=r_size) / real(nens - 1, kind=r_size)
       else
          assigned_val = default_val
       end if
       field(:,:) = assigned_val
    else
       ! init_mode == 0 または未定義のモード
       field(:,:) = default_val
    end if

    ! どの初期化モードでも許容範囲内にする
    field(:,:) = max(min_val, min(max_val, field(:,:)))

  end subroutine init_perturb_2d


  subroutine update_perturb_2d(field, im, jm, update_mode, &
                               min_val, max_val, amp)
    integer, intent(in) :: im, jm, update_mode
    real(kind=r_size), intent(in) :: min_val, max_val, amp
    real(kind=r_size), intent(inout) :: field(im, jm)

    integer :: i, j
    real(kind=r_size) :: rand_val
    real(kind=r_size), allocatable :: rand_field(:,:)

    if (update_mode == 1) then
       ! 全格子点に同じノイズを加える
       call random_number(rand_val)
       field(:,:) = field(:,:) + (2.0_r_size * rand_val - 1.0_r_size) * amp

    else if (update_mode == 2) then
       ! 格子点ごとに異なるノイズを加える
       allocate(rand_field(im, jm))
       call random_number(rand_field)

       do j = 1, jm
          do i = 1, im
             field(i,j) = field(i,j) + &
                  (2.0_r_size * rand_field(i,j) - 1.0_r_size) * amp
          end do
       end do

       deallocate(rand_field)
    end if
    ! update_mode == 0 または未定義のモードではノイズを加えない

    ! ノイズを加えない場合も含めて、毎回許容範囲内にする
    field(:,:) = max(min_val, min(max_val, field(:,:)))

  end subroutine update_perturb_2d

end module mod_perturb_param
