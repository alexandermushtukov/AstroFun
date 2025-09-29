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



!=========================================================================================================
!   stat(1) - goes to the infinity
!   stat(2) - absorbed by 1sr star
!   stat(3) - obsorbed by 2nd star
!   stat(4) - form disc aroun binary
!=========================================================================================================
subroutine binary_wind_accreretion()
implicit none
real*8::stat(4),m_1,m_2,P_day,v6_ini
integer::n_w
  201 format (8(es11.4,"   "))
  n_w = int(8.e3)
  m_1 = 2.d0     !== NS mass ==!
  m_2 = 8.d0     !== companion mass: NGC 7793 P13 - 10.-20.; M51 ULX-7 >8. ==!
  P_day = 5.d0  !5.7d0   !== NGC 5907 ULX-1: 5.66; NGC 7793 P13 - 65; M51 ULX7 - 2  ==!
  v6_ini = 150.d0
  write(*,*)"# P_day, v6_ini, outflow frac, abs frac"
  do while(v6_ini.le.250.d0)
    call trace_particle_in_binary(stat,m_1,m_2,P_day,v6_ini,n_w)
    write(*,201)P_day,v6_ini,stat(1:2)
    v6_ini = v6_ini*1.1
  end do
return
end subroutine binary_wind_accreretion


!====================================================================
! Subroutine simulates motion of particle in a binary.
! RF rotates together with a binary.
! Terminates when:
!   det_res=0 : particle is sufficiently far AND unbound
!   det_res=1 : collision with star 1
!   det_res=2 : collision with star 2
!   det_res=3 : time > 10 orbital periods
!====================================================================
subroutine trace_particle_in_binary(stat,m_1,m_2,P_day,v6_ini,n_w)
implicit none
real*8,intent(in)::m_1,m_2,P_day,v6_ini
integer,intent(in):: n_w
integer::det_res, iw
real*8::pi=3.141592653589793d0
real*8::a8,a8_1,a8_2,r8(3),v6(3),a4(3),r8_1(3),r8_2(3)
real*8::dt2,t2,t2_max,omega_cu
real*8::a_cen4(3),a_cor4(3),ag4_1(3),ag4_2(3),delta_r1(3),delta_r2(3),a_tot4(3)
real*8::geom3d_length,delta_3d  !== functions ==!
real*8::d8_1,d8_2,RL8_1,RL8_2,q,r_st8_1,r_st8_2,random,theta,fi,v6_rf(3),v6_lab(3),acc,t_print,dt_print,omega
real*8::E_spec, r_norm
real*8::stat(4)
  !== parameters ==!
  omega = 2.d0 * pi / (P_day * 86400.d0)   !== [rad/s], P_day в сутках ==!
  q = m_1/m_2
  r_st8_1 = 1000.    !== radius of 1st star in [1.e8 cm] ==!
  r_st8_2 = 1000.    !== radius of 2nd star in [1.e8 cm] ==!
  !================!

  a8 = 2.9d3 * m_1**(1./3) * (1.d0+m_2/m_1)**(1./3) * P_day**(2./3)  !== separation b/w stars ==!
  !write(*,*)"# ",a8,r_st8_1,r_st8_2
  a8_1 = m_2/(m_1+m_2)*a8
  a8_2 = a8 - a8_1
  r8_1(1) = +a8_1; r8_1(2:3)=0.d0    !== coordinates of 1st star ==!
  r8_2(1) = -a8_2; r8_2(2:3)=0.d0    !== coordinates of 2nd star ==!
  RL8_1 = a8 * 0.49*q**(2./3)/( 0.6*q**(2./3)+log( 1. + q**(1./3) ) )
  q = 1./q
  RL8_2 = a8 * 0.49*q**(2./3)/( 0.6*q**(2./3)+log( 1. + q**(1./3) ) )

  !== time limit: 10 orbital periods in code units (uses same omega constant as below) ==!
  omega_cu = 7.2722d-3 / P_day
  t2_max   = 10.d0 * 2.d0*pi / omega_cu

  stat(1:4)=0.d0  !== stat(1) - goes to the infinity, stat(2) - absorbed by 1sr star, stat(3) - obsorbed by 2nd star, stat(4) - form disc aroun binary ==!
  !== cycle over particles ==!
  do iw = 1, n_w

    t_print = 0.d0; dt_print = 40.d0
    det_res = -1

    !== initial parameters of particle ==!
    !== random place of particle start from the 2nd star ==!
    call RANDOM_NUMBER(random); theta=acos(1.d0-random)
    call RANDOM_NUMBER(random); fi = 2*pi*random
    r8(1) = sin(theta)*cos(fi);  r8(2) = sin(theta)*sin(fi);  r8(3) = cos(theta)
    v6(1:3) = v6_ini * r8(1:3) / geom3d_length(r8)   !== particle is emitted along normal to the stellar atmosphere ==!
    !== correction for position of a star ==!
    r8(1:3) = (1.02d0*r_st8_2) * r8(1:3) + r8_2(1:3)  !== we start particle a little bit above stellar surface ==!
    !====================================!

    d8_1 = delta_3d(r8,r8_1)
    d8_2 = delta_3d(r8,r8_2)

    dt2 = 1.d0
    t2 = 0.d0
    do
      delta_r1(1:3) = r8_1(1:3) - r8(1:3)
      ag4_1(1:3) = 1.328d6 * m_1 * delta_r1(1:3)/ ( geom3d_length(delta_r1) )**3
      delta_r2(1:3) = r8_2(1:3) - r8(1:3)
      ag4_2(1:3) = 1.328d6 * m_2 * delta_r2(1:3)/ ( geom3d_length(delta_r2) )**3

      !== centrifugal acceleration ==!
      a_cen4(1:2) = 5.2885d-5 / P_day**2 * r8(1:2)
      a_cen4(3) = 0.d0
      !==============================!

      !== coriolis acceleration ==!
      a_cor4(1) = +2 * 7.2722d-3 / P_day * v6(2)
      a_cor4(2) = -2 * 7.2722d-3 / P_day * v6(1)
      a_cor4(3) = 0.d0
      !============================!

      a_tot4(1:3) = ag4_1(1:3) + ag4_2(1:3) + a_cen4(1:3) + a_cor4(1:3)
      acc = geom3d_length(a_tot4)

      !== adaptive timestep ==!
      dt2 = 0.05d0 * sqrt(10.d0 / acc)

      !== advance (position + velocity in rotating RF) ==!
      r8(1:3) = r8(1:3) + dt2*( v6(1:3) + dt2*a_tot4(1:3)/2.d0 )
      v6(1:3) = v6(1:3) + dt2 * a_tot4(1:3)

      !== distances after move ==!
      d8_1 = delta_3d(r8,r8_1)
      d8_2 = delta_3d(r8,r8_2)
      r_norm = geom3d_length(r8)

      !== component of velocity due to RF rotation & lab-frame velocity ==!
      v6_rf(1) = -7.27d-3 / P_day * r8(2)
      v6_rf(2) = +7.27d-3 / P_day * r8(1)
      v6_rf(3) = 0.d0
      v6_lab(1:3) = v6(1:3) + v6_rf(1:3)

      !== check collision ==!
      !if( d8_1.lt.r_st8_1 )then
      if( d8_1 .lt. RL8_1 )then
        det_res = 1
        exit
      end if
      if( d8_2 .lt. r_st8_2 )then
        det_res = 2
        exit
      end if

      !== check "sufficiently far and unbound" ==!
      E_spec = 0.5d0*geom3d_length(v6_lab)**2 - 1.328d6*( m_1/d8_1 + m_2/d8_2 )
      if( (r_norm.gt.40.d0*a8) .and. (E_spec.gt.0.d0) )then
        det_res = 0
        exit
      end if

      !== diagnostics ==!
      if( t2.gt.t_print )then
        !write(*,'(10(ES13.6,"  "))') t2, r_norm/a8, d8_1/RL8_1, d8_2/RL8_2, &
        !                             geom3d_length(v6), geom3d_length(v6_rf), geom3d_length(v6_lab), geom3d_length(a_tot4)
        t_print = t_print + dt_print
      end if

      t2 = t2 + dt2

      !== time limit (>10 orbits) ==!
      if( t2.gt.t2_max )then
        det_res = 3
        exit
      end if
    end do
    stat(det_res+1)=stat(det_res+1)+1.d0
  end do  !== iw loop ==!
  stat(1:4) = stat(1:4)/SUM(stat)
return
end subroutine trace_particle_in_binary
