!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Operator for precip fields to be fed to CoMorph precipramp.
!> @details Defines a coarsening operation that should roughly preserve ppnrate scaling
!!          of entrainment given coarsened precip fields.
!!          This kernel only works for the lowest-order W3 and Wtheta spaces and assumes that
!!          resdep_precipramp is true.

module sci_inv_precipramp_kernel_mod

use constants_mod,           only: i_def, r_double, r_single
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   GH_READ, GH_WRITE,         &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_2, &
                                   GH_COARSE, GH_FINE, CELL_COLUMN

implicit none

private

type, public, extends(kernel_type) :: inv_precipramp_kernel_type
   private
   type(arg_type) :: meta_args(4) = (/                                   &
        arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1, &
                                              mesh_arg=GH_COARSE),       &
        arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2, &
                                              mesh_arg=GH_FINE  ),       &
        arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1, &
                                              mesh_arg=GH_COARSE),       &
        arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2, &
                                              mesh_arg=GH_FINE  )        &
        /)
  integer :: operates_on = CELL_COLUMN
end type inv_precipramp_kernel_type

public :: inv_precipramp_kernel_code

  ! Generic interface for real32 and real64 types
  interface inv_precipramp_kernel_code
    module procedure  &
      inv_precipramp_code_r_single, &
      inv_precipramp_code_r_double
  end interface

contains

  !> @brief Thein restriction operation for scalar fields.
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_y  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] coarse_field             Coarse grid field to write to
  !> @param[in]     fine_field               The fine grid field to restrict
  !> @param[in]     delta_coarse             Grid-spacing for the coarse grid
  !> @param[in]     delta_fine               Grid spacing for the fine grid
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid

  ! R_SINGLE PRECISION
  ! ==================
  subroutine inv_precipramp_code_r_single(                          &
                                           nlayers,                 &
                                           cell_map,                &
                                           ncell_fine_per_coarse_x, &
                                           ncell_fine_per_coarse_y, &
                                           ncell_fine,              &
                                           coarse_field,            &
                                           fine_field,              &
                                           delta_coarse,            &
                                           delta_fine,              &
                                           undf_coarse,             &
                                           map_coarse,              &
                                           ndf,                     &
                                           undf_fine,               &
                                           map_fine)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in)    :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)
    integer(kind=i_def), intent(in)    :: ncell_fine
    integer(kind=i_def), intent(in)    :: ndf
    integer(kind=i_def), intent(in)    :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in)    :: map_coarse(ndf)
    integer(kind=i_def), intent(in)    :: undf_fine, undf_coarse
    real(kind=r_single), intent(inout) :: coarse_field(undf_coarse)
    real(kind=r_single), intent(in)    :: fine_field(undf_fine)
    real(kind=r_single), intent(in)    :: delta_coarse(undf_coarse)
    real(kind=r_single), intent(in)    :: delta_fine(undf_fine)

    integer(kind=i_def) :: df, k, x_idx, y_idx, top_df
    real(kind=r_single) :: coarse_value(nlayers-1+ndf), denom
    ! Hardcoded variables needed for inverting precipramp scaling
    real(kind=r_single) :: par_radius_knob, par_radius_knob_max, par_radius_ppn_max
    real(kind=r_single) :: ppnrate, dxfac

    ! Assume lowest order W3 or Wtheta space
    df = 1
    ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
    top_df = nlayers - 2 + ndf
    ! Initialise variables
    coarse_value = 0.0_r_single
    denom = 1.0_r_single/real(ncell_fine_per_coarse_x*ncell_fine_per_coarse_y, kind=r_single)
    ! Assign hardcoded values for comorph precipramp
    par_radius_knob = 0.45_r_single
    par_radius_knob_max = 2.0_r_single
    par_radius_ppn_max = 3.5e-4_r_single

    ! Build up 1D array of new coarse values for this column
    do k = 0, top_df
      do y_idx = 1, ncell_fine_per_coarse_y
        do x_idx = 1, ncell_fine_per_coarse_x
          !coarse_value(k+1) =  min(coarse_value(k+1), &
          !  fine_field(map_fine(df,cell_map(x_idx,y_idx))+k))
          ! dxfac = delta_fine(map_fine(df,cell_map(x_idx,y_idx))+k) / 50.0e3_r_single
          dxfac = 1.0_r_single
          ppnrate = dxfac * fine_field(map_fine(df,cell_map(x_idx,y_idx))+k)
          coarse_value(k+1) = coarse_value(k+1) + (1.0_r_single/(1.0_r_single + &
                              ((par_radius_knob_max/par_radius_knob) - 1.0_r_single) * &
                              (ppnrate / par_radius_ppn_max))) * denom
        end do
      end do
    end do

    ! Copy over values into coarse field and invert it back to ppnrate
    do k = 0, top_df
      ! dxfac = delta_coarse(map_coarse(df) + k) / 50.0e3_r_single
      dxfac = 1.0_r_single
      coarse_field(map_coarse(df) + k) = par_radius_ppn_max * ((1.0_r_single/coarse_value(k+1)) - 1.0_r_single) / &
                                         (((par_radius_knob_max/par_radius_knob) - 1.0_r_single) * dxfac)
    end do

  end subroutine inv_precipramp_code_r_single

  ! R_DOUBLE PRECISION
  ! ==================
  subroutine inv_precipramp_code_r_double(              &
                                           nlayers,                 &
                                           cell_map,                &
                                           ncell_fine_per_coarse_x, &
                                           ncell_fine_per_coarse_y, &
                                           ncell_fine,              &
                                           coarse_field,            &
                                           fine_field,              &
                                           delta_coarse,            &
                                           delta_fine,              &
                                           undf_coarse,             &
                                           map_coarse,              &
                                           ndf,                     &
                                           undf_fine,               &
                                           map_fine)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in)    :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)
    integer(kind=i_def), intent(in)    :: ncell_fine
    integer(kind=i_def), intent(in)    :: ndf
    integer(kind=i_def), intent(in)    :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in)    :: map_coarse(ndf)
    integer(kind=i_def), intent(in)    :: undf_fine, undf_coarse
    real(kind=r_double), intent(inout) :: coarse_field(undf_coarse)
    real(kind=r_double), intent(in)    :: fine_field(undf_fine)
    real(kind=r_double), intent(in)    :: delta_coarse(undf_coarse)
    real(kind=r_double), intent(in)    :: delta_fine(undf_fine)

    integer(kind=i_def) :: df, k, x_idx, y_idx, top_df
    real(kind=r_double) :: coarse_value(nlayers-1+ndf), denom
    ! Hardcoded variables needed for inverting precipramp scaling
    real(kind=r_double) :: par_radius_knob, par_radius_knob_max, par_radius_ppn_max
    real(kind=r_double) :: ppnrate, dxfac

    ! Assume lowest order W3 or Wtheta space
    df = 1
    ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
    top_df = nlayers - 2 + ndf
    ! Initialise variables
    coarse_value = 0.0_r_double
    denom = 1.0_r_double/real(ncell_fine_per_coarse_x*ncell_fine_per_coarse_y, kind=r_double)
    ! Assign hardcoded values for comorph precipramp
    par_radius_knob = 0.45_r_double
    par_radius_knob_max = 2.0_r_double
    par_radius_ppn_max = 3.5e-4_r_double

    ! Build up 1D array of new coarse values for this column
    do k = 0, top_df
      do y_idx = 1, ncell_fine_per_coarse_y
        do x_idx = 1, ncell_fine_per_coarse_x
          !coarse_value(k+1) =  min(coarse_value(k+1), &
          !  fine_field(map_fine(df,cell_map(x_idx,y_idx))+k))
          ! dxfac = delta_fine(map_fine(df,cell_map(x_idx,y_idx))+k) / 50.0e3_r_double
          dxfac = 1.0_r_double
          ppnrate = dxfac * fine_field(map_fine(df,cell_map(x_idx,y_idx))+k)
          coarse_value(k+1) = coarse_value(k+1) + (1.0_r_double/(1.0_r_double + &
                              ((par_radius_knob_max/par_radius_knob) - 1.0_r_double) * &
                              (ppnrate / par_radius_ppn_max))) * denom
        end do
      end do
    end do

    ! Copy over values into coarse field and invert it back to ppnrate
    do k = 0, top_df
      ! dxfac = delta_coarse(map_coarse(df) + k) / 50.0e3_r_double
      dxfac = 1.0_r_double
      coarse_field(map_coarse(df) + k) = par_radius_ppn_max * ((1.0_r_double/coarse_value(k+1)) - 1.0_r_double) / &
                                         (((par_radius_knob_max/par_radius_knob) - 1.0_r_double) * dxfac)
    end do

  end subroutine inv_precipramp_code_r_double


end module sci_inv_precipramp_kernel_mod
