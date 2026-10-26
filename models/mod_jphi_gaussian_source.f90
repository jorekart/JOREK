!> Tilted toroidally-symmetric Gaussian current source.
!>
!> Defines a circular (single-width) 2D Gaussian current distribution J_phi(R,Z),
!> purely in the toroidal direction (J = (0,0,J_phi)), then evaluates the current
!> density that results from rigidly tilting that distribution about the cartesian
!> y-axis (through R=0, Z=0) by a fixed angle. The amplitude ramps up linearly in
!> time from 0 at t=0 to jgs_amplitude at t=jgs_ramp_time, then stays constant
!> (jgs_ramp_time<=0 means no ramp, full amplitude immediately). The input
!> parameters (jgs_amplitude, jgs_R0, jgs_Z0, jgs_width, jgs_tilt_angle,
!> jgs_ramp_time) are declared in phys_module.f90.
module mod_jphi_gaussian_source

  use phys_module, only: jgs_amplitude, jgs_R0, jgs_Z0, jgs_width, jgs_tilt_angle, jgs_ramp_time, t_now
  use mod_coordinate_transforms, only: cylindrical_to_cartesian, cartesian_to_cylindrical, &
                                        vector_cylindrical_to_cartesian, vector_cartesian_to_cylindrical, &
                                        rotate_vectors_cart_y
  implicit none
  private

  public :: jphi_gaussian_axisym
  public :: gaussian_source_to_cartesian
  public :: rotate_gaussian_source_cartesian
  public :: cartesian_source_to_cylindrical
  public :: jphi_gaussian_rotated_source

contains

  !> rotate a single 3-vector about the cartesian y-axis by a given angle.
  !> rotate_vectors_cart_y operates on (3,n_v) batches of vectors; this wraps
  !> it for the single-vector case used throughout this module.
  function rotate_cart_y_vec3(angle,v) result(v_rot)
    real*8, intent(in) :: angle
    real*8, intent(in) :: v(3)
    real*8              :: v_rot(3)
    real*8              :: tmp(3,1)

    tmp(:,1) = v
    call rotate_vectors_cart_y(1,angle,tmp)
    v_rot = tmp(:,1)
  end function rotate_cart_y_vec3

  !> Step 1: toroidally symmetric, circular (single-width) 2D Gaussian J_phi(R,Z).
  pure function jphi_gaussian_axisym(R,Z) result(Jphi)
    real*8, intent(in) :: R,Z !< position where the source is evaluated
    real*8             :: Jphi

    Jphi = jgs_amplitude * exp( -0.5d0*((R-jgs_R0)**2 + (Z-jgs_Z0)**2) / jgs_width**2 )
  end function jphi_gaussian_axisym

  !> Linear ramp factor for the source amplitude: rises linearly from 0 at t=0
  !> to 1 at t=jgs_ramp_time, then stays at 1. jgs_ramp_time<=0 means no ramp
  !> (always fully on).
  pure function jgs_ramp_scale(t) result(scale)
    real*8, intent(in) :: t
    real*8             :: scale

    if (jgs_ramp_time .le. 0.d0) then
      scale = 1.d0
    else
      scale = min(1.d0, max(0.d0,t) / jgs_ramp_time)
    endif
  end function jgs_ramp_scale

  !> time derivative of jgs_ramp_scale.
  pure function jgs_ramp_scale_rate(t) result(rate)
    real*8, intent(in) :: t
    real*8             :: rate

    if (jgs_ramp_time .le. 0.d0 .or. t .lt. 0.d0 .or. t .ge. jgs_ramp_time) then
      rate = 0.d0
    else
      rate = 1.d0 / jgs_ramp_time
    endif
  end function jgs_ramp_scale_rate

  !> Step 2: express a point of the axisymmetric source, and its purely-toroidal
  !> current vector there, in cartesian (x,y,z) coordinates.
  subroutine gaussian_source_to_cartesian(R,Z,phi,xyz,Jxyz)
    real*8, intent(in)  :: R,Z,phi !< position in JOREK (R,Z,phi) coordinates
    real*8, intent(out) :: xyz(3)  !< position in cartesian (x,y,z) coordinates
    real*8, intent(out) :: Jxyz(3) !< current density in cartesian (Jx,Jy,Jz) coordinates
    real*8              :: Jphi

    Jphi = jphi_gaussian_axisym(R,Z)
    xyz  = cylindrical_to_cartesian( (/R,Z,phi/) )
    Jxyz = vector_cylindrical_to_cartesian(phi, (/0.d0,0.d0,Jphi/))
  end subroutine gaussian_source_to_cartesian

  !> Step 3: rigidly rotate a cartesian position and its current vector about the
  !> cartesian y-axis (through R=0, Z=0) by a given angle.
  subroutine rotate_gaussian_source_cartesian(angle,xyz,Jxyz)
    real*8, intent(in)    :: angle   !< rotation angle about the cartesian y-axis
    real*8, intent(inout) :: xyz(3)  !< position in cartesian coordinates
    real*8, intent(inout) :: Jxyz(3) !< current density in cartesian coordinates

    xyz  = rotate_cart_y_vec3(angle,xyz)
    Jxyz = rotate_cart_y_vec3(angle,Jxyz)
  end subroutine rotate_gaussian_source_cartesian

  !> Step 4: transform a rotated cartesian position and current vector back to
  !> JOREK (R,Z,phi) coordinates.
  subroutine cartesian_source_to_cylindrical(xyz,Jxyz,R,Z,phi,J_R,J_Z,J_phi)
    real*8, intent(in)  :: xyz(3)  !< position in cartesian coordinates
    real*8, intent(in)  :: Jxyz(3) !< current density in cartesian coordinates
    real*8, intent(out) :: R,Z,phi
    real*8, intent(out) :: J_R,J_Z,J_phi
    real*8              :: cyl(3), Jcyl(3)

    cyl = cartesian_to_cylindrical(xyz)
    R = cyl(1); Z = cyl(2); phi = cyl(3)
    Jcyl = vector_cartesian_to_cylindrical(phi, Jxyz)
    J_R = Jcyl(1); J_Z = Jcyl(2); J_phi = Jcyl(3)
  end subroutine cartesian_source_to_cylindrical

  !> Evaluate the tilted, time-ramped Gaussian current source at an arbitrary
  !> query point (R,Z,phi), returning the resulting current density components
  !> J_R, J_Z, J_phi there, along with their time derivatives dJ_R_dt, dJ_Z_dt,
  !> dJ_phi_dt (from the linear amplitude ramp only; not currently used by any
  !> equation, but exported for future use e.g. a term involving dJ/dt).
  !> This is the routine to call from a matrix assembly / Gauss-point loop.
  !>
  !> Because the tilt displaces the source, evaluating it at an arbitrary query
  !> point (rather than only at the points the forward rotation happens to land on)
  !> requires locating where that query point came from before the tilt: the query
  !> point is mapped backward (rotated by -jgs_tilt_angle, i.e. steps 2+3 run in
  !> reverse on the position only) to find the corresponding point on the original,
  !> untilted axisymmetric source; the source is evaluated there (step 1); and the
  !> resulting current vector is rotated forward by +jgs_tilt_angle (step 3 on the
  !> vector) before being expressed back in (J_R,J_Z,J_phi) components at the query
  !> point's actual phi (step 4). This correctly reproduces "take the axisymmetric
  !> source, tilt the whole object, then ask what the tilted object's current
  !> density is at this point" for any point, not just the image of a chosen input.
  !>
  !> The amplitude ramp only rescales the magnitude of the (fixed-shape, fixed-tilt)
  !> source in time, so its time derivative carries through the same (linear)
  !> coordinate transforms as the value itself: evaluate the untilted profile's
  !> value and rate once, then push both through the identical rotate+convert steps.
  subroutine jphi_gaussian_rotated_source(R,Z,phi,J_R,J_Z,J_phi,dJ_R_dt,dJ_Z_dt,dJ_phi_dt)
    real*8, intent(in)  :: R,Z,phi                !< query point in JOREK coordinates
    real*8, intent(out) :: J_R,J_Z,J_phi          !< current density at the query point
    real*8, intent(out) :: dJ_R_dt,dJ_Z_dt,dJ_phi_dt !< time derivative of the above (ramp only)

    real*8 :: xyz0(3), Jxyz(3), dJxyz_dt(3)
    real*8 :: cyl0(3), R0, Z0, phi0, Jphi0_full, scale, rate
    real*8 :: Jcyl(3), dJcyl_dt(3)

    ! locate the pre-image of the query point under the tilt
    xyz0 = cylindrical_to_cartesian( (/R,Z,phi/) )
    xyz0 = rotate_cart_y_vec3(-jgs_tilt_angle,xyz0)
    cyl0 = cartesian_to_cylindrical(xyz0)
    R0 = cyl0(1); Z0 = cyl0(2); phi0 = cyl0(3)

    ! evaluate the original (full-amplitude) axisymmetric source there, apply
    ! the time ramp, and carry both the value and its time derivative forward
    ! through the same tilt
    Jphi0_full = jphi_gaussian_axisym(R0,Z0)
    scale      = jgs_ramp_scale(t_now)
    rate       = jgs_ramp_scale_rate(t_now)

    Jxyz     = vector_cylindrical_to_cartesian(phi0, (/0.d0,0.d0,Jphi0_full*scale/))
    dJxyz_dt = vector_cylindrical_to_cartesian(phi0, (/0.d0,0.d0,Jphi0_full*rate /))

    Jxyz     = rotate_cart_y_vec3(jgs_tilt_angle,Jxyz)
    dJxyz_dt = rotate_cart_y_vec3(jgs_tilt_angle,dJxyz_dt)

    Jcyl     = vector_cartesian_to_cylindrical(phi, Jxyz)
    dJcyl_dt = vector_cartesian_to_cylindrical(phi, dJxyz_dt)

    J_R = Jcyl(1); J_Z = Jcyl(2); J_phi = Jcyl(3)
    dJ_R_dt = dJcyl_dt(1); dJ_Z_dt = dJcyl_dt(2); dJ_phi_dt = dJcyl_dt(3)
  end subroutine jphi_gaussian_rotated_source

end module mod_jphi_gaussian_source
