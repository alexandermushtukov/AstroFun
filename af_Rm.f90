!================================================================================
! The function calculates the magnetosperic radius using the standard formula.
!================================================================================
real*8 function Rm_norm8(B12,L37,Lambda,m,R6)
implicit none
real*8,intent(in)::B12,L37,Lambda,m,R6
  Rm_norm8=2.5d0*Lambda*B12**(4.d0/7.d0)/L37**(2.d0/7.d0)*(m/1.4d0)**(1.d0/7.d0)*R6**(10.d0/7.d0)
return
end function Rm_norm8



!================================================================================
! The function calculates the magnetosperic radius using the standard formula.
!================================================================================
real*8 function Rm_norm8_(B12,dotm,Lambda,m,R6)
implicit none
real*8,intent(in)::B12,dotm,Lambda,m,R6
real*8::L37
  L37=26.666d0*dotm*m**2/R6
  Rm_norm8_=2.5d0*Lambda*B12**(4.d0/7.d0)/L37**(2.d0/7.d0)*(m/1.4d0)**(1.d0/7.d0)*R6**(10.d0/7.d0)
return
end function Rm_norm8_


!============================================================================
!============================================================================
!subroutine Rmag_SC(dotm0,m,R6,B12,Lambda,eps_w,Rmag8,dotm_Rm)

subroutine Rmag_SC(dotm0,m,R6,B12,Lambda,eps_w,Rmag8,dotm_Rm)
implicit none
real*8,intent(in)::dotm0,m,R6,B12
real*8::eps_w,C_r_to_r8,C_r8_to_r,Rmag8,dotm_Rm
real*8::Rm_norm8_,r1,r2,dotm,Lambda
integer::i
  C_r_to_r8=m*9.d-3 !m/3.333d2
  C_r8_to_r=1.d0/C_r_to_r8
  i=1
  dotm=dotm0
  r1=Rm_norm8_(B12,dotm0,Lambda,m,R6)*C_r8_to_r
  do while(i.le.15)
    !write(*,*)r1,dotm
    dotm=dotm_r(r1,dotm0,eps_w)  !==dotM at Rm
    r1=Rm_norm8_(B12,dotm,Lambda,m,R6)*C_r8_to_r  !==new Rm
    i=i+1
  end do
  Rmag8=r1*C_r_to_r8
  dotm_Rm=dotm
return
contains
  real*8 function dotm_r(r,dotm0,eps_w)
  implicit none
  real*8,intent(in)::r,dotm0,eps_w
  real*8::m_in
  real*8::r_sp  !==function==!
    m_in=dotm_in(dotm0,eps_w)
    if(r.lt.r_sp(dotm0,eps_w))then
      dotm_r=m_in+(dotm0-m_in)*r/r_sp(dotm0,eps_w)
    else
      dotm_r=dotm0
    end if
  return
  end function dotm_r

  real*8 function dotm_in(dotm0,eps_w)
  implicit none
  real*8,intent(in)::dotm0,eps_w
  real*8::a
    a=eps_w*(0.83d0-0.25d0*eps_w)
    dotm_in=dotm0*(1.d0-a)/(1.d0-a/(0.4d0*dotm0))
  return
  end function dotm_in
end subroutine Rmag_SC


!========================================================================================
! The function gives an approximate spherization radius in units if gravitational radius.
! mdot0 - the initial mass accretion rate in units of Eddington mass accretion rate.
! eps_w - the fraction of heat, which is going into the wind launhing.
!========================================================================================
real*8 function r_sp(mdot0,eps_w)
implicit none
real*8,intent(in)::mdot0,eps_w
  r_sp=mdot0*(1.34d0-0.4d0*eps_w+0.1d0*eps_w**2-(1.1d0-0.7d0*eps_w)*mdot0**(-2./3.))
return
end function r_sp
