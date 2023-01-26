!======================================================================================================
! The module contains subroutines transforming the angles in local reference frame to the angles in
! B-field reference frame and back.
!======================================================================================================
module AngTransformation
contains
  subroutine AngTrans_Atm2Bfield(theta,fi,theta_Bz,theta_B,fi_B)
  implicit none
  real*8,intent(in)::theta,fi,theta_Bz
  real*8::theta_B,fi_B
  real*8::x,y,z,xb,yb,zb
    if(theta_Bz.eq.0.d0)then
      theta_B=theta
    else
      call Polar2Decart(theta,fi,x,y,z)
      call VecRotation3d(x,y,z,2,-theta_Bz,xb,yb,zb)
      call Decart2Polar(xb,yb,zb,theta_B,fi_B)
    end if
  return
  end subroutine AngTrans_Atm2Bfield

  subroutine AngTrans_Bfield2Atm(theta_B,fi_B,theta_Bz,theta,fi)
  implicit none
  real*8,intent(in)::theta_B,fi_B,theta_Bz
  real*8::theta,fi
  real*8::x,y,z,xb,yb,zb
    if(theta_Bz.eq.0.d0)then
      theta=theta_B
    else
      call Polar2Decart(theta_B,fi_B,xb,yb,zb)
      call VecRotation3d(xb,yb,zb,2,theta_Bz,x,y,z)
      call Decart2Polar(x,y,z,theta,fi)
    end if
  return
  end subroutine AngTrans_Bfield2Atm
end module AngTransformation



!=================================================================================================================!
! Angle transformation from B-field reference frame to observer's reference frame
! The function returns latitude in observer's RF:
! i.e. the angle between direction to the observer and direction to the point on the NS surface
!=================================================================================================================!
real*8 function B2obs(ksi,theta_B,fi_B)
implicit none
real*8,intent(in)::ksi,theta_B,fi_B
  !B2obs=acos(cos(theta_B)*cos(ksi)-sin(fi_B)*sin(theta_B)*cos(ksi))
  B2obs=acos(cos(theta_B)*cos(ksi)-sin(fi_B)*sin(theta_B)*sin(ksi))
return
end function B2obs
!=================================================================================================================!



!=========================================================================
! Azimutal angle in observer's RF.
!=========================================================================
real*8 function B2obs_fi(ksi,theta_B,fi_B)
implicit none
real*8,intent(in)::ksi,theta_B,fi_B
real*8::B2obs                !==function==!
real*8::theta_0,sin_fi0,cos_fi0
  theta_0 = B2obs(ksi,theta_B,fi_B)
  if( sin(theta_0).eq.0.d0 )then
    B2obs_fi=fi_B !0.d0
  else
    sin_fi0=(sin(fi_B)*sin(theta_B)*cos(ksi)+cos(theta_B)*sin(ksi))/sin(theta_0)
    cos_fi0=cos(fi_B)*sin(theta_B)/sin(theta_0)
    B2obs_fi=acos(max(-1.d0,min(cos_fi0,1.d0)))
    if( sin_fi0.lt.0.d0 )then
      B2obs_fi=-B2obs_fi
    end if
  end if
return
end function B2obs_fi
!=========================================================================



!=================================================================================================================!
! Angle transformation from observer's reference frame to B-field reference frame
! The function returns the latitude in B-field RF:
! i.e. angle between the magnetic field direction and direction to the point on the NS surface.
! SBCh.
!=================================================================================================================!
real*8 function obs2B(ksi,theta0,fi0)
implicit none
real*8,intent(in)::ksi,theta0,fi0
  obs2B=acos(sin(ksi)*sin(theta0)*sin(fi0)+cos(ksi)*cos(theta0))
return
end function obs2B
!=================================================================================================================!
