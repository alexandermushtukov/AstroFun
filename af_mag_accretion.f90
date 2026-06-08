!=======================================================================================================
! The module contains functions and subroutines to calculate free-free absorption in a strong B-field
! according to Meszaros 1992.
!=======================================================================================================
module AstroFun_magnetic_accretion
contains
  !=======================================================================================
  ! The function calculates the magnetosperic radius as a function of accretion luminosity
  ! using the standard formula.
  !=======================================================================================
  real*8 function Rm_norm8(B12,L37,Lambda,m,R6)
  implicit none
  real*8,intent(in)::B12,L37,Lambda,m,R6
    Rm_norm8=2.5d0*Lambda*B12**(4.d0/7)/L37**(2.d0/7)*(m/1.4d0)**(1.d0/7)*R6**(10.d0/7)
  return
  end function Rm_norm8

  !=======================================================================================
  ! The function calculates the magnetosperic radius as a function of mass accretion rate
  ! using the standard formula.
  ! dotm - mass accretion rate in Eddington units.
  !=======================================================================================
  real*8 function Rm_norm8_(B12,dotm,Lambda,m,R6)
  implicit none
  real*8,intent(in)::B12,dotm,Lambda,m,R6
  real*8::L37
    L37=26.666d0*dotm*m**2/R6
    Rm_norm8_=2.5d0*Lambda*B12**(4.d0/7)/L37**(2.d0/7)*(m/1.4d0)**(1.d0/7)*R6**(10.d0/7)
  return
  end function Rm_norm8_

  !========================================================================================
  ! The function calculates magnetic field strength [1.e12 G] at the NS surface from known
  ! inner disc radius, Lambda-parameter, accretion luminosity, NS mass and radius.
  !========================================================================================
  real*8 function Rm2B12(Rm8,Lambda,L37,m,R10)
  implicit none
  real*8,intent(in)::Rm8,Lambda,L37,m,R10
    Rm2B12=( Rm8/2.5d0/Lambda*L37**(2./7)/(m/1.4)**(1./7)/R10**(10./7) )**(7./4)
  return
  end function Rm2B12

  !=========================================================================================
  ! The function calculates the corotational radius in units of 1.e8 cm.
  ! m - NS mass in units of Solar masses, P - spin period [s].
  !=========================================================================================
  real*8 function R_corot_8(m,P)
  implicit none
  real*8,intent(in)::m,P
    R_corot_8=1.68d0*(m/1.4d0)**(1.d0/3)*P*(2.d0/3)
  return
  end function R_corot_8

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


  !====================================================================================================================
  ! The subroutine calculates the mass accretion rate at R_m accounting for mass losses due to the winds from the disc.
  !   dotm0 - mass accretion rate from a companion star in units of Eddington mass accretion rates.
  !   m, R6, B12 - NS mass, radius and surface magnetic field strength
  !   Rmag8 - inner disc radius
  !   dotm_Rm - mass accrewtion rate at the inner disc radius.
  !====================================================================================================================
  subroutine Rmag_SC(dotm0,m,R6,B12,Lambda,eps_w,Rmag8,dotm_Rm)
  implicit none
  real*8,intent(in)::dotm0,m,R6,B12
  real*8::eps_w,C_r_to_r8,C_r8_to_r,Rmag8,dotm_Rm
  real*8::r1,r2,dotm,Lambda
  integer::i
    !Lambda=1.d0
    !eps_w=0.5d0          !==efficientcy of outflow==!
    C_r_to_r8=m*9.d-3 !m/3.333d2
    C_r8_to_r=1.d0/C_r_to_r8
    i=1
    dotm=dotm0
    r1=Rm_norm8_(B12,dotm0,Lambda,m,R6)*C_r8_to_r
    do while(i.le.15)
      dotm=dotm_r(r1,dotm0,eps_w)                   !==mass accretion rate at Rm
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

end module AstroFun_magnetic_accretion



!================================================================================================
! The subroutine calculats the velosity profile in accretion envelope accounting for mass
! outflow from accretion disc.
!   mdot - initial mass accretion rate in units of Eddngton mass accretion rates.
!   mas(i,1) - lambda, mas(i,2) - beta, mas(i,3) - tau.
!================================================================================================
subroutine VelosityProf_2(Lambda,eps_w,mdot,B12,P,m,mas,n,tau_min,t,Rmag8,L39)
use AstroFun_magnetic_accretion
implicit none
real*8::pi=3.141592653589793d0
real*8,intent(in)::Lambda,eps_w,mdot,B12,P,m
real*8::mas(n,3)  !== the array contains velosity and tau
integer,intent(in)::n
real*8::L39,omega,Rm,H2Rm,omega_k,T_d,beta_0
real*8::l,dl,r,beta,cos_l,xi,tau,f1,f2,ff,g1,g2,gg,P_tot,P_mag,l_max
integer::i
real*8::t, tau_min,R6,Rmag8,mdot_Rm

  R6=1.d0  !==NS radius==!

  omega=2*pi/P
  call Rmag_SC(mdot,m,R6,B12,Lambda,eps_w,Rmag8,mdot_Rm)
  !write(*,*)"Rmag8 ",Rmag8

  L39=0.138d0*m*mdot_Rm
  Rm=Rmag8*1.d8

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  H2Rm = min(0.1d0*L39/m/(Rm/1.d8),0.7d0)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  omega_k=sqrt(1.33d26*m/Rm**3)
  T_d=0.25d0*L39**(13./28.)/B12**(3./7.)
  !beta_0=sqrt(2.d0*T_d/9.d5)
  !beta_0=Rm/3.d10*sqrt(3.d0/2.d0*(omega_k-omega)**2)
  beta_0=0.056d0/sqrt(Lambda)*m**(3./7.)*L39**(1./7.)/B12**(2./7.)*abs(1.d0-omega/omega_k)
  !beta_0=sqrt(m*3.d5/Rm)

  !write(*,*)"# ",beta_0,Rm,omega_k,omega,atan(0.5d0*H2Rm)

  tau_min = 1.d10
  t = 0.d0
  beta = beta_0 !/1.1d0
  !dl=0.0001d0
  l = atan(0.5d0*H2Rm)
  !write(*,*)l,H2Rm
  !read(*,*)
  l_max = acos(sqrt(1.d6/Rm))
  dl = (l_max-l)/(n-1)
  r = Rm*cos(l)*cos(l)
  i = 1
  do while(i.le.n)
    cos_l=cos(l)
    xi=atan(1.d0/tan(l)/2)
    tau=0.12d0*L39**(9.d0/7)/B12**(4.d0/7)/beta/Lambda**(1.5d0)/m**(3.d0/14)/cos_l**3
    !f1=cos(xi)*(1.d0-5.d0*L39/m/tau)  !==with radiation pressure==!
    f1=cos(xi)                         !==without radiation pressure==!
    f2=2.6d-3*omega**2*Lambda**3/m**(4.d0/7)*B12**(12.d0/7)/L39**(6.d0/7)*cos_l**7 * cos(pi+xi-l)
    ff=2.1d-3*m**(6.d0/7)*sqrt(4.d0-3*cos_l**2)*L39**(2.d0/7)/Lambda/B12**(4.d0/7)/cos_l**3
    g1=sin(xi)*(1.d0-5.d0/m/tau)
    g2=2.6d-3*omega**2*Lambda**3/m**(4.d0/7)*B12**(12.d0/7)/L39**(6.d0/7)*cos_l**7 * sin(pi+xi-l)
    gg=1.1d10*tau*m*L39**(4.d0/7)/Lambda**2/m**(2.d0/7)/B12**(8.d0/7)/cos_l**4
    P_tot=gg*(g1+g2)   !==why "-"?
    P_mag=1.63d13/B12**(10.d0/7)*L39**(12.d0/7)/cos_l**12

    beta=beta+ (dl/beta)*ff*(f1+f2)
    mas(i,1)=l
    mas(i,2)=beta
    mas(i,3)=tau
    if(tau_min.gt.tau)then
     tau_min=tau
    end if
    if(beta.le.0.d0)then
      r=13.d0
    end if
    l=l+dl
    r=Rm*cos(l)*cos(l)
    t=t+dl*Rm*cos(l)*sqrt(4.d0-3.d0*cos(l)*cos(l))/(beta*3.d10)
    i=i+1
  end do
return
end subroutine VelosityProf_2

