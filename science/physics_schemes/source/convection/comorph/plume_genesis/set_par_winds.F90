! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module set_par_winds_mod

implicit none

contains

! Subroutine to set parcel initial fields that are assumed equal in
! all the sub-grid regions:
! - winds
! - tracers
! - parcel radius
subroutine set_par_winds( n_points, n_points_super, n_fields_tot,              &
                          cmpr_init, k, l_tracer, l_down,                      &
                          fields_k, turb_pert_k, par_radius_k,                 &
                          par_gen_par, par_gen_mean, par_gen_core )

use comorph_constants_mod, only: real_cvprec, one, par_gen_core_fac,           &
                                 n_tracers, name_length, l_par_core,           &
                                 i_check_bad_values_cmpr, i_check_bad_none
use fields_type_mod, only: i_wind_u, i_wind_w, i_tracers,                      &
                           i_temperature, i_q_vap, i_qc_first, i_qc_last,      &
                           field_positive, field_names
use parcel_type_mod, only: n_par, i_radius, i_edge_virt_temp
use cmpr_type_mod, only: cmpr_type

use calc_virt_temp_mod, only: calc_virt_temp
use check_bad_values_mod, only: check_bad_values_cmpr
use raise_error_mod, only: raise_warning, raise_fatal

implicit none

! Number of points
integer, intent(in) :: n_points

! Size of input super-arrays (maybe larger than n_points,
! as arrays are dimensioned with the max size they could need
! on any level and reused for all levels)
integer, intent(in) :: n_points_super

! Total number of fields in the input super-arrays,
! including tracers if they are used
integer, intent(in) :: n_fields_tot

! Indices of points being processed, for making error messages more informative
type(cmpr_type), intent(in) :: cmpr_init
integer, intent(in) :: k

! Flag for whether tracers are in use
logical, intent(in) :: l_tracer

! Flag for updraft vs downdraft
logical, intent(in) :: l_down

! Primary model-fields at level k
real(kind=real_cvprec), intent(in) :: fields_k                                 &
                                      ( n_points_super, n_fields_tot )

! Turbulence-based perturbations to u,v,w,Tl,qt, at level k
real(kind=real_cvprec), intent(in) :: turb_pert_k                              &
                                      ( n_points, i_wind_u:i_q_vap )

! Parcel initial radius at level k
real(kind=real_cvprec), intent(in) :: par_radius_k(n_points)

! Updraft or downdraft initiaing properties:
! Super-array containing parcel radius and edge virtual temperature
real(kind=real_cvprec), intent(in out) :: par_gen_par                          &
                                          ( n_points, n_par )
! In-parcel means of primary fields
real(kind=real_cvprec), intent(in out) :: par_gen_mean                         &
                                          ( n_points, n_fields_tot )
! Parcel core primary fields
real(kind=real_cvprec), intent(in out) :: par_gen_core                         &
                                          ( n_points, n_fields_tot )

! Factor scaling the turbulence-based perturbations
real(kind=real_cvprec) :: factor

! String containing info to print in error messages
character(len=64) :: call_string
character(len=name_length) :: field_name

! Loop counters
integer :: ic, i_field

character(len=*), parameter :: routinename = "SET_PAR_WINDS"

call raise_warning(routinename, &
 "Start set_par_winds, n_points = "//trim(adjustl(str(n_points)))//", n_par = "//trim(adjustl(str(n_par))))
call raise_warning(routinename, &
 "Parcel radius copy: i_radius = "//trim(adjustl(str(i_radius))))
write(*,*) 'Test write statement'

do ic = 1, n_points
  ! Copy parcel radius into the parcel
  par_gen_par(ic,i_radius) = par_radius_k(ic)
end do

call raise_warning(routinename, &
 "Call write virt temp")
! Set initial edge virtual temperature to the environment value at k
call calc_virt_temp( n_points, n_points_super,                                 &
                     fields_k(:,i_temperature),                                &
                     fields_k(:,i_q_vap),                                      &
                     fields_k(:,i_qc_first:i_qc_last),                         &
                     par_gen_par(:,i_edge_virt_temp) )

! Set sign of perturbation to add; sign is reversed for downdrafts
if ( l_down ) then
  factor = -one
else
  factor = one
end if
call raise_warning(routinename, &
  "Set in-parcel mean winds")
! Set in-parcel mean winds
do i_field = i_wind_u, i_wind_w
  do ic = 1, n_points
    par_gen_mean(ic,i_field) = fields_k(ic,i_field)                            &
                             + factor * turb_pert_k(ic,i_field)
  end do
end do

! Parcel core has perturbations scaled up by par_gen_core_fac
factor = factor * par_gen_core_fac
call raise_warning(routinename, &
  "Set core winds")
! Set parcel core winds
do i_field = i_wind_u, i_wind_w
  do ic = 1, n_points
    par_gen_core(ic,i_field) = fields_k(ic,i_field)                            &
                             + factor * turb_pert_k(ic,i_field)
  end do
end do
call raise_warning(routinename, &
  "Copy tracer values")
! For now, copy grid-mean tracer values into the parcel
! (planning to add a turbulence-based perturbation to the
!  tracer fields if l_turb_par_gen, but not yet implemented).
if ( l_tracer .and. n_tracers > 0 ) then
  do i_field = i_tracers(1), i_tracers(n_tracers)
    do ic = 1, n_points
      par_gen_mean(ic,i_field) = fields_k(ic,i_field)
    end do
  end do
  if ( l_par_core ) then
    do i_field = i_tracers(1), i_tracers(n_tracers)
      do ic = 1, n_points
        par_gen_core(ic,i_field) = fields_k(ic,i_field)
      end do
    end do
  end if
end if

if ( i_check_bad_values_cmpr > i_check_bad_none ) then
  ! Check outputs for bad values (NaN, Inf etc).
  call raise_warning(routinename, "name_length assigned to call_string is "//trim(adjustl(str(name_length))))
  !call raise_fatal( routinename, "Force flushing before setting call_string" ) 
  if ( l_down ) then
    call_string = "On output from set_par_winds; dndraft"
  else
    call_string = "On output from set_par_winds; updraft"
  end if
  
  do i_field = i_wind_u, i_wind_w
    call raise_warning(routinename, "Check parcel mean winds")
    field_name = "par_gen_mean_" // trim(adjustl(field_names(i_field)))
    call check_bad_values_cmpr( cmpr_init, k, par_gen_mean(:,i_field),         &
                                call_string, field_name,                       &
                                field_positive(i_field) )
    if ( l_par_core ) then
      call raise_warning(routinename, "Check core winds")
      field_name = "par_gen_core_" // trim(adjustl(field_names(i_field)))
      call check_bad_values_cmpr( cmpr_init, k, par_gen_core(:,i_field),       &
                                  call_string, field_name,                     &
                                  field_positive(i_field) )
    end if
  end do
  if ( l_tracer .and. n_tracers > 0 ) then
    call raise_warning(routinename, 'Check tracer values')
    do i_field = i_tracers(1), i_tracers(n_tracers)
      field_name = "par_gen_mean_" // trim(adjustl(field_names(i_field)))
      call check_bad_values_cmpr( cmpr_init, k, par_gen_mean(:,i_field),       &
                                  call_string, field_name,                     &
                                  field_positive(i_field) )
      if ( l_par_core ) then
        field_name = "par_gen_core_" // trim(adjustl(field_names(i_field)))
        call check_bad_values_cmpr( cmpr_init, k, par_gen_core(:,i_field),     &
                                    call_string, field_name,                   &
                                    field_positive(i_field) )
      end if
    end do
  end if

end if  ! ( i_check_bad_values_cmpr > i_check_bad_none )


return
end subroutine set_par_winds

! Helper function for casting integers to string
pure function str(i) result(s)
    integer, intent(in) :: i
    character(:), allocatable :: s

    character(32) :: tmp

    write(tmp, '(I0)') i
    s = trim(tmp)
end function

end module set_par_winds_mod
