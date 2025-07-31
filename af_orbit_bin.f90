!=======================================================================================================================
! The subroutine calculates the variations in orbit due to the mass transfer in a close binary system.
! Input parameters:
!    m - NS mass in units of Solar masses
!    m_c - mass of a companion star
!    P_d - orbital period in a system in days
!    dotM20 - mass accretion rate from a companion star in units of 1.e20 g/s
!    B12,R6 - magnetic field strength at NS surface and NS radius
!    det=1 - accretion through the Lagrangian point, det=2 - wind accretion
! Output parameters:
!    dot_a - major semi-axis derivative [cm/s]
!    a - major semi-axis [cm]
!    r_lc - NS Roshe lobe radius [cm]
!    Rmag8 - inner disc radius [1.e8 cm]
!    dotM20_delta - the expected mass loss rate from the companion calculated under assumption of fixed radius of a star.
!=======================================================================================================================
subroutine b_orbit_dynamics(dot_a,a,r_lc,Rmag8,dotM20_delta,m,m_c,P_d,dotM20,B12,R6,det)
use AstroFun_magnetic_accretion
implicit none
real*8,intent(in)::m,m_c,P_d,dotM20,B12,R6
integer,intent(in)::det
real*8::q,beta,alpha,dotm0,Rmag8,dotm_Rm,dot_a,a,r_lc,dotM20_delta
real*8::rho_L1,T3_L1
real*8::pi=3.141592653589793d0
real*8::R_Roshe_lobe  !==function==!
  q=m/m_c
  if(det.eq.1)then
    !==Lagrangian point==!
    dotm0=dotM20/0.01875d0
    call Rmag_SC(dotm0,m,R6,B12,0.5d0,0.d0,Rmag8,dotm_Rm)
    beta=dotm_Rm/dotm0
    alpha=1.d0/(1.d0+q)**2        !==alpha for the case, when mass losses from the NS surface (winds from inner ragions of accretion disc)==!
  else
    !==wind accretion==!
    beta=0.d0
    alpha=q**2/(1.d0+q)**2           !==alpha for the case, when mass losses from companion star==!
  end if
  dot_a=1.47d-2*(m+m_c)**(1.d0/3)/m_c*P_d**(2.d0/3)*dotM20&
        *( 2*(1.d0-beta/q) - (1.d0-beta)/(1.d0+q) - 2*alpha*(1.d0-beta)*(1.d0+q)/q )
  a=2.94d11*(m+m_c)**(1.d0/3)*P_d**(2.d0/3)
  r_lc=R_Roshe_lobe(m_c,m,a)

  T3_L1=3.d0                                           !==expected tamperature of the stellar surface==!
  rho_L1=2.2d-5*dotM20/T3_L1**(3./2)/P_d**2            !==expected mass density at L1==!
  dotM20_delta=4*pi*r_lc**2*abs(dot_a)*rho_L1/1.d20
return
end subroutine b_orbit_dynamics


subroutine test_b_orbit()
implicit none
real*8::m,m_c,P_d,dotM20,B12,R6,dot_a,a,Rmag8,dotM20_delta,r_lc
real*8::R_Roshe_lobe   !==function==!
  201 format (8(es11.4,"   "))
  m=1.4d0
  m_c=5.d0
  P_d=23.d0
  B12=1.d0
  R6=1.d0
  dotM20=1.d-4
  write(*,*)
  do while(dotM20.le.1.d1)
    call b_orbit_dynamics(dot_a,a,r_lc,Rmag8,dotM20_delta,m,m_c,P_d,dotM20,B12,R6,2)
    write(*,201)dotM20,dot_a,a,r_lc,Rmag8*1.d8,dotM20_delta
    dotM20=dotM20*1.1d0
  end do
return
end subroutine test_b_orbit
!=============================================================================================================


!==================================================================================
! The function calculates the radius of a Roche lobe around mass m_c.
! m, m_c - masses
! a - separation between the components in a binary.
!==================================================================================
real*8 function R_Roshe_lobe(m_c,m,a)
implicit none
real*8,intent(in)::m,m_c,a
real*8::q
  q=m/m_c
  R_Roshe_lobe=0.49d0/q**(2.d0/3)/( 0.6d0/q**(2.d0/3)+log(1.d0+1.d0/q**(1.d0/3)) ) * a
return
end function R_Roshe_lobe
!==================================================================================




!=============================================================================================================
! The subroutine calculated the fraction (!) of mass accretion rate from the wind in a binary system.
! m, m_c - NS mass and mass of companion in units of Solar mass.
! P_orb_d - orbital parion in days.
! Circular orbit is assumed.
!=============================================================================================================
real*8 function M_dot_frac_wind(m,m_c,P_orb_d)
implicit none
real*8,intent(in)::m,m_c,P_orb_d
real*8::q,a_10,v_fi_9,v_wind_9,tan_b,ksi
real*8::pi=3.141592653589793d0
real*8::alpha_2,alpha_3,R_c_10,v_inf_9     !==parameters, see Lipunov (76.II)==!
  R_c_10=3.d1   !==radius of a companion star==!
  alpha_2=1.d0
  alpha_3=1.d0
  ksi=1.d0          !==parameter==!
  q=m/m_c
  a_10 = 2.94d1*(m+m_c)**(1.d0/3)*P_orb_d**(2.d0/3)        !==orbital separation==!
  v_fi_9=1.15d0*(m_c/(a_10*100))                           !==Keplerian velocity [cm/s]==!
  !v_inf_9=1.15d0*sqrt(m_c/(R_c_10*100)) *sqrt(2.d0)       !==parabolic velocity at the companion star surface==!
  v_inf_9=7.d-2                 !==typical value from El Mellah paper==!
  v_wind_9=v_inf_9*(1.d0 - (R_c_10/a_10)**alpha_2)**alpha_3
  tan_b=v_fi_9/v_wind_9
  M_dot_frac_wind = ksi*tan_b**4/pi/(1.d0+tan_b**2)**1.5d0 *q**2*(1.d0+q)**2
return
end function M_dot_frac_wind


subroutine test_M_dot_frac_wind()
implicit none
real*8::m,m_c,P_orb_d,res
real*8::M_dot_frac_wind  !==function==!
  m=1.4d0
  m_c=10.d0
  P_orb_d=77.5d0
  do while(P_orb_d.le.10.d0)
    res=M_dot_frac_wind(m,m_c,P_orb_d)
    write(*,*)P_orb_d,res
    P_orb_d=P_orb_d+0.3d0
  end do
return
end subroutine test_M_dot_frac_wind
!========================================================================================================


!====================================================================
! Subroutine simulates motion of particle in a binary.
! RF rotates together with a binary.
!====================================================================
subroutine trace_particle_in_binary()
implicit none
real*8::pi=3.141592653589793d0
real*8::P_day,m_1,m_2,a8,a8_1,a8_2,r8(3),v6(3),a4(3),r8_1(3),r8_2(3)
real*8::dt2,t2
real*8::a_cen4(3),a_cor4(3),ag4_1(3),ag4_2(3),delta_r1(3),delta_r2(3),a_tot4(3)
real*8::geom3d_length,delta_3d  !== functions ==!
real*8::d8_1,d8_2,RL8_1,RL8_2,q,r_st8_1,r_st8_2,random,theta,fi,v6_ini,v6_rf(3),v6_lab(3),acc,t_print,dt_print
  !== parameters ==!
  P_day = 10.d0   !== orbital period ==!
  m_1 = 2.d0     !== mass of 1st star ==!
  m_2 = 2.d0     !== mass of 2nd star ==!
  q = m_1/m_2
  r_st8_1 = 1000.   !== radius of 1st star in [1.e8 cm] ==!
  r_st8_2 = 1000.   !== radius of 2nd star in [1.e8 cm] ==!
  v6_ini = 90.d0     !== initial wind velocity [1.e6 cm/s] ==!
  !================!

  t_print = 0.d0; dt_print = 20.d0

  a8 = 2.9d3 * m_1**(1./3) * (1.d0+m_2/m_1)**(1./3) * P_day**(2./3)  !== separation b/w stars ==!
  a8_1 = m_2/(m_1+m_2)*a8
  a8_2 = a8 - a8_1
  r8_1(1) = +a8_1; r8_1(2:3)=0.d0    !== coordinates of 1st star ==!
  r8_2(1) = -a8_2; r8_2(2:3)=0.d0    !== coordinates of 2nd star ==!
  RL8_1 = a8 * 0.49*q**(2./3)/( 0.6*q**(2./3)+log( 1. + q**(1./3) ) )
  q = 1./q
  RL8_2 = a8 * 0.49*q**(2./3)/( 0.6*q**(2./3)+log( 1. + q**(1./3) ) )
  !write(*,*)m_1,m_2,a8,RL8_1,RL8_2
  !read(*,*)

  !== initial parameters of particle ==!
  ! r8(1:2) = 0.d0; r8(3) = a8    !== initial coordinate of a partile ==!
  ! v6(1:3) = 0.d0                !== initial velocity of a particle ==!
  !== random place of particle start from the 2nd star ==!
  call RANDOM_NUMBER(random); theta=acos(1.d0-random)
  call RANDOM_NUMBER(random); fi = 2*pi*random
  theta = 0.d0
  fi = 0.d0
  r8(1) = sin(theta)*cos(fi);  r8(2) = sin(theta)*sin(fi);  r8(3) = cos(theta)
  v6(1:3) = v6_ini * r8(1:3) / geom3d_length(r8)   !== particle is emitted along normal to the stellar atmosphere ==!
  !== correction for position of a star ==!
  r8(1:3) = (1.02d0*r_st8_2) * r8(1:3) + r8_2(1:3) !== we start particle a little bit above stellar surface ==!
  !====================================!


  d8_1 = delta_3d(r8,r8_1)
  d8_2 = delta_3d(r8,r8_2)

  dt2 = 1.d0
  t2 = 0.d0
  do while(t2 .le. 2.d5)
    delta_r1(1:3) = r8_1(1:3) - r8(1:3)
    ag4_1(1:3) = 1.328d6 * m_1 * delta_r1(1:3)/ ( geom3d_length(delta_r1) )**3
    delta_r2(1:3) = r8_2(1:3) - r8(1:3)
    ag4_2(1:3) = 1.328d6 * m_2 * delta_r2(1:3)/ ( geom3d_length(delta_r2) )**3

    !== centrifugal acceleration ==!
    a_cen4(1:2) = 5.3d-5 / P_day**2 * r8(1:2)
    a_cen4(3) = 0.d0
    !==============================!

    !== coriolis acceleration ==!
    a_cor4(1) = +2 * 7.27d-3 / P_day * v6(2)
    a_cor4(2) = -2 * 7.27d-3 / P_day * v6(1)
    a_cor4(3) = 0.d0
    !============================!
    a_tot4(1:3) = ag4_1(1:3) + ag4_2(1:3) + a_cen4(1:3) + a_cor4(1:3)
    acc = geom3d_length(a_tot4)
    dt2 = 0.1 * sqrt(10.d0 / acc)
    r8(1:3) = r8(1:3) + dt2*( v6(1:3)+dt2*a_tot4(1:3)/2 )
    v6(1:3) = v6(1:3) + dt2 * a_tot4(1:3)     !== note: it is velocity in rotating RF ==!
    d8_1 = delta_3d(r8,r8_1)
    d8_2 = delta_3d(r8,r8_2)
    if( (d8_1.lt.r_st8_1).or.(d8_2.lt.r_st8_2) )then
      exit
    end if
    !== component of velocity due to RF rotation ==!
    v6_rf(1) = +7.27d-3/P_day * r8(1)
    v6_rf(2) = +7.27d-3/P_day * r8(2)
    v6_rf(3) = 0.d0
    v6_lab(1:3) = v6(1:3) - v6_rf(1:3)
    if( t2.gt.t_print )then
      write(*,'(10(ES13.6,"  "))') t2, d8_1/RL8_1, d8_2/RL8_2, geom3d_length(v6), geom3d_length(v6_rf), geom3d_length(v6_lab), geom3d_length(a_tot4)
      t_print = t_print + dt_print
    end if
    !write(*,*) v6(1:3),dt2*a_tot4(1:3)
    !write(*,'(10(ES13.6,"  "))') t2,v6(1:3),v6_rf(1:3),v6_lab(1:3)
    !read(*,*)
    t2 = t2 + dt2
  end do

return
end subroutine trace_particle_in_binary
