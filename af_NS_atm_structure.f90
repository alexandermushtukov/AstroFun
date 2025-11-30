!==========================================================================================================!
! ...
! mas_x_rho_tau_2(n,1) - x
! mas_x_rho_tau_2(n,2) - rho
! mas_x_rho_tau_2(n,3) - tau
! mas_x_rho_tau_2(n,4) - eps24
! mas_x_rho_tau_2(n,5) - T_keV (?)
!==========================================================================================================
subroutine test_acc_atm_structure()
implicit none
real*8,allocatable::mas_x_rho_tau_2(:,:),F_tau(:,:)
real*8::rho_min,rho_max,g14,T_keV,dot_m_6,ln_Lambda,m,R6,mas_tau_TkeV(100,2),x_scale,drhodx
integer::i,n
  n = 1000
  allocate( mas_x_rho_tau_2(n,5),F_tau(n,2) )

  m = 1.4d0; R6 = 1.d0
  dot_m_6 = 1.1d0
  g14 = 1.328d0*m/R6**2
  T_keV = 10.d0
  ln_Lambda = 40.d0        !== Coulomb logarithm ==!

  i=1
  do while(i.le.100)
    mas_tau_TkeV(i,1) = 0.1d0*i;  mas_tau_TkeV(i,2) = T_keV
    i=i+1
  end do

  rho_min = 0.7d-4*dot_m_6
  rho_max = 1.d+1
  !call acc_atm_structure_2(mas_x_rho_tau_2,n,rho_min,rho_max,x_scale,g14,T_keV,dot_m_6,ln_Lambda,F_tau)
  !call acc_atm_structure_3(mas_x_rho_tau_2,n,rho_min,rho_max,x_scale,g14,mas_tau_TkeV,100,dot_m_6,ln_Lambda,F_tau)
   call acc_atm_structure_4(mas_x_rho_tau_2,n,1.d-2,5.d2,x_scale, &
                              g14,mas_tau_TkeV,100,dot_m_6,ln_Lambda,rho_min,F_tau)
  i = 2
  do while(i.le.n-1)
    drhodx = ( mas_x_rho_tau_2(i,3) - mas_x_rho_tau_2(i-1,3) )/( mas_x_rho_tau_2(i,1) - mas_x_rho_tau_2(i-1,1) )
    write(*,*)mas_x_rho_tau_2(i,1:3),drhodx;! read(*,*)
    !write(*,*)i,F_tau(i,1:2)
    i=i+10
  end do
  write(*,*)
return
end subroutine test_acc_atm_structure
!================================================================================



!============================================================================================
! The subroutine calculates structure of NS atmosphere influenced by accretion process.
! Parameters:
!    mas_x_rho_2 - array that contains dependence of mass dencity on coordinate in [cm]
!    n - number of lines in the array
!    rho_min,rho_max - the dencity interval that is under consideration
! mas_x_rho_tau_2(n,1) - x
! mas_x_rho_tau_2(n,2) - rho
! mas_x_rho_tau_2(n,3) - tau
! mas_x_rho_tau_2(n,4) - eps24
! mas_x_rho_tau_2(n,5) - T_keV (?)
! Add: possibility of temperature profile - the current version assumes fixed temperature.
!============================================================================================
subroutine acc_atm_structure_2(mas_x_rho_tau_2,n,rho_min,rho_max,x_scale,g14,T_keV,dot_m_6,ln_Lambda,F_tau)
implicit none
real*8::mas_x_rho_tau(n,5),mas_x_rho_tau_2(n,5),x_scale
integer,intent(in)::n
real*8,intent(in)::rho_min,rho_max,g14,T_keV,dot_m_6,ln_Lambda
real*8::rho,rho2,dx,x,drho,drho2,tau2,x_min,x_max
integer::i
real*8::beta_ff,beta,Z,dbeta_dx,dTkeV_dx,T_keV_,T_keV_0,eps24,F_tau(n,2)

  beta_ff = 0.5d0
  Z = 1.d0
  dTkeV_dx = 0.d0    !== initially ==!
  T_keV_ = 1.d0; T_keV_0 = 1.d0   !== "fake" temperature ==!
 
  drho = (rho_max-rho_min)/n
  rho = rho_min - drho
  dx = 1.d-3
  x = 0.d0
  i=1
  tau2 = 0.d0     !== local optical depth ==!
  beta = beta_ff
  rho2 = 1.d-4
  do while(rho.lt.rho_max)
    dbeta_dx = 3.24d-4*rho2*Z**2/beta**3*ln_Lambda
    beta = max(0.d0,beta - dx*dbeta_dx)
    if( beta.gt.0.d0 )then
      drho2 = ( 0.5*0.1*g14/T_keV*rho2 + 31.15*dot_m_6*0.5/T_keV*dbeta_dx )*dx   !== 0.5 because of H-plasma ==!
      eps24 = 9.e2*dot_m_6*beta*dbeta_dx
      T_keV_ = ( eps24/1.75d0/rho2**2 )**2!/100000
      dTkeV_dx = ( T_keV_0 - T_keV_ )/dx
    else
      drho2 = ( 0.5*0.1*g14/T_keV*rho2 + 0.d0 )*dx
    end if
    !write(*,*)"v2:", x, dx, rho2, drho2, dbeta_dx, beta!; read(*,*)
    rho2 = rho2 + drho2
    tau2 = tau2 + 0.34*rho2*dx
    T_keV_0 = T_keV_
    222 if( (rho2.ge.(rho+drho)) )then
      if( i.le.n )then
        mas_x_rho_tau(i,1) = x + dx*( rho+drho - (rho2-drho2) )/drho2
        mas_x_rho_tau(i,2) = rho+drho
        mas_x_rho_tau(i,3) = tau2
        if(beta.ne.0.d0)then
          mas_x_rho_tau(i,4) = eps24 !9.e2*dot_m_6*beta*dbeta_dx   !== [1.e24 erg/cm^3/s]  ==!
          mas_x_rho_tau(i,5) = T_keV_
        else
          mas_x_rho_tau(i,4) = 0.d0
          mas_x_rho_tau(i,5) = 0.d0
        end if
      end if
      rho = rho + drho
      i=i+1
      goto 222
    else
      goto 333
    end if
    333 x = x+dx
  end do

  F_tau(1,1)=0.d0; F_tau(1,2)=0.d0;
  i=2
  do while( i.le.n )
    F_tau(i,1)=mas_x_rho_tau(i,3)-mas_x_rho_tau(1,3)
    F_tau(i,2)=F_tau(i-1,2) + mas_x_rho_tau(i,4)*(mas_x_rho_tau(i,1)-mas_x_rho_tau(i-1,1))
    i=i+1
  end do
  F_tau(1:n,2)=F_tau(1:n,2)/F_tau(n,2)   !== the commulative function of sources disctribution ==!

  !== inverting array ==!
  x_min = mas_x_rho_tau(i,1); x_max = mas_x_rho_tau(n,1)
  i=1
  do while(i.le.n)
    mas_x_rho_tau_2(i,1) = x_max - mas_x_rho_tau(n-i+1,1)
    mas_x_rho_tau_2(i,2) = mas_x_rho_tau(n-i+1,2)
    mas_x_rho_tau_2(i,3) = mas_x_rho_tau(n-i+1,3)
    mas_x_rho_tau_2(i,4) = mas_x_rho_tau(n-i+1,4)
    mas_x_rho_tau_2(i,5) = mas_x_rho_tau(n-i+1,5)
    i=i+1
  end do
  x_scale = x_max - x_min
return
end subroutine acc_atm_structure_2
!===============================================================================


!============================================================================================
! The subroutine calculates structure of NS atmosphere influenced by accretion process.
! Parameters:
!    mas_x_rho_2 - array that contains dependence of mass dencity on coordinate in [cm]
!    n - number of lines in the array
!    rho_min,rho_max - the dencity interval that is under consideration
! mas_x_rho_tau_2(n,1) - x
! mas_x_rho_tau_2(n,2) - rho
! mas_x_rho_tau_2(n,3) - tau
! mas_x_rho_tau_2(n,4) - eps24
! mas_x_rho_tau_2(n,5) - T_keV (?)
! This version of a code assumes some variations of temperature with optical depth.
!============================================================================================
subroutine acc_atm_structure_3(mas_x_rho_tau_2,n,rho_min,rho_max,x_scale,g14,mas_tau_TkeV,n_mas,dot_m_6,ln_Lambda,F_tau)
implicit none
real*8::mas_x_rho_tau(n,5),mas_x_rho_tau_2(n,5),x_scale
integer,intent(in)::n,n_mas
real*8,intent(in)::rho_min,rho_max,g14,dot_m_6,ln_Lambda,mas_tau_TkeV(n_mas,2)
real*8::rho,rho2,dx,x,drho,drho2,tau2,x_min,x_max,T_keV
integer::i
real*8::beta_ff,beta,Z,dbeta_dx,dTkeV_dx,T_keV_,T_keV_0,eps24,F_tau(n,2)
real*8::find_H_ordered_inc  !== function ==!

  beta_ff = 0.5d0
  Z = 1.d0
  dTkeV_dx = 0.d0    !== initially ==!
  T_keV_ = 1.d0; T_keV_0 = 1.d0   !== "fake" temperature ==!
 
  drho = (rho_max-rho_min)/n
  rho = rho_min !- drho
  dx = 1.d-3
  x = 0.d0
  i=1
  tau2 = 0.d0     !== local optical depth ==!
  beta = beta_ff
  rho2 = rho_min !1.d-4
  !write(*,*)rho_min,rho_max; read(*,*)
  do while(rho.lt.rho_max)
    !write(*,*)beta,rho2; read(*,*)
    dbeta_dx = 3.24d-4*rho2*Z**2/beta**3*ln_Lambda
    beta = max(0.d0,beta - dx*dbeta_dx)
    if( beta.gt.0.d0 )then
      T_keV = find_H_ordered_inc(tau2,mas_tau_TkeV,n_mas)
      drho2 = ( 0.5*0.1*g14/T_keV*rho2 + 31.15*dot_m_6*0.5/T_keV*dbeta_dx )*dx   !== 0.5 because of H-plasma ==!
      eps24 = 9.e2*dot_m_6*beta*dbeta_dx
      T_keV_ = ( eps24/1.75d0/rho2**2 )**2!/100000
      dTkeV_dx = ( T_keV_0 - T_keV_ )/dx
    else
      drho2 = ( 0.5*0.1*g14/T_keV*rho2 + 0.d0 )*dx
    end if

    !write(*,*)"v2:", x, dx, rho2, drho2, dbeta_dx, beta!; read(*,*)
    rho2 = rho2 + drho2

    tau2 = tau2 + 0.34*rho2*dx
    T_keV_0 = T_keV_
    222 if( (rho2.ge.(rho+drho)) )then
      if( i.le.n )then
        mas_x_rho_tau(i,1) = x + dx*( rho+drho - (rho2-drho2) )/drho2
        mas_x_rho_tau(i,2) = rho+drho
        mas_x_rho_tau(i,3) = tau2
        if(beta.ne.0.d0)then
          mas_x_rho_tau(i,4) = eps24 !9.e2*dot_m_6*beta*dbeta_dx   !== [1.e24 erg/cm^3/s]  ==!
          mas_x_rho_tau(i,5) = T_keV_
        else
          mas_x_rho_tau(i,4) = 0.d0
          mas_x_rho_tau(i,5) = 0.d0
        end if
      end if
      rho = rho + drho
      i=i+1
      goto 222
    else
      goto 333
    end if
    333 x = x+dx
  end do

  F_tau(1,1)=0.d0; F_tau(1,2)=0.d0;
  i=2
  do while( i.le.n )
    F_tau(i,1)=mas_x_rho_tau(i,3)-mas_x_rho_tau(1,3)
    F_tau(i,2)=F_tau(i-1,2) + mas_x_rho_tau(i,4)*(mas_x_rho_tau(i,1)-mas_x_rho_tau(i-1,1))
    i=i+1
  end do
  F_tau(1:n,2)=F_tau(1:n,2)/F_tau(n,2)   !== the commulative function of sources disctribution ==!

  !== non(!) inverting array ==!
  x_min = mas_x_rho_tau(i,1); x_max = mas_x_rho_tau(n,1)
  i=1
  do while(i.le.n)
    mas_x_rho_tau_2(i,1) = mas_x_rho_tau(i,1)
    mas_x_rho_tau_2(i,2) = mas_x_rho_tau(i,2)
    mas_x_rho_tau_2(i,3) = mas_x_rho_tau(i,3)
    mas_x_rho_tau_2(i,4) = mas_x_rho_tau(i,4)
    mas_x_rho_tau_2(i,5) = mas_x_rho_tau(i,5)
    i=i+1
  end do
  x_scale = x_max - x_min
return
end subroutine acc_atm_structure_3
!===============================================================================


!============================================================================================
! OUTPUT GRID: tau targets are [0, logspace(tau_min..tau_max)]:
!   mas_x_rho_tau_2(1,:) at tau=0
!   mas_x_rho_tau_2(2..n,:) log-spaced between tau_min and tau_max
!
! Requires tau_min > 0.
!============================================================================================
subroutine acc_atm_structure_4(mas_x_rho_tau_2,n,tau_min,tau_max,x_scale, &
                              g14,mas_tau_TkeV,n_mas,dot_m_6,ln_Lambda,rho0,F_tau)
  implicit none

  integer, intent(in)  :: n, n_mas
  real*8,  intent(in)  :: tau_min, tau_max, g14, dot_m_6, ln_Lambda, rho0
  real*8,  intent(in)  :: mas_tau_TkeV(n_mas,2)

  real*8,  intent(out) :: mas_x_rho_tau_2(n,5), x_scale
  real*8,  intent(out) :: F_tau(n,2)

  ! local
  integer :: i, steps, max_steps, k
  real*8  :: dx, x, x_old
  real*8  :: tau2, tau_old, target_tau, f
  real*8  :: rho2, rho_old, drho2
  real*8  :: beta, beta_ff, dbeta_dx, Z
  real*8  :: T_keV, T_keV_, T_old
  real*8  :: eps24, eps_old
  real*8  :: rho_floor, beta_floor, T_floor
  real*8  :: tau_targets(n), log_ratio, expo

  real*8  :: find_H_ordered_inc
  external :: find_H_ordered_inc

  rho_floor  = 1.d-30
  beta_floor = 1.d-30
  T_floor    = 1.d-30

  beta_ff = 0.5d0
  Z       = 1.d0

  dx   = 1.d-3
  x    = 0.d0
  tau2 = 0.d0

  rho2 = max(rho0, rho_floor)
  beta = beta_ff

  eps24  = 0.d0
  T_keV_ = 0.d0

  !---- build tau grid: tau_targets(1)=0, tau_targets(2..n)=logspace(tau_min..tau_max)
  tau_targets(1) = 0.d0

  if (n .eq. 1) then
    ! nothing else
  else if (n .eq. 2) then
    tau_targets(2) = tau_max
  else
    ! require tau_min>0 and tau_max>tau_min
    log_ratio = log(tau_max / tau_min)
    do k = 2, n
      expo = dble(k-2) / dble(n-2)   ! runs 0..1
      tau_targets(k) = tau_min * exp( expo * log_ratio )
    end do
    tau_targets(n) = tau_max
  end if

  !---- first row exactly at tau=0
  mas_x_rho_tau_2(1,1) = 0.d0
  mas_x_rho_tau_2(1,2) = rho2
  mas_x_rho_tau_2(1,3) = 0.d0
  mas_x_rho_tau_2(1,4) = 0.d0
  mas_x_rho_tau_2(1,5) = 0.d0

  i = 2
  target_tau = tau_targets(i)

  steps = 0
  max_steps = 200000000  ! safety cap

  do while (i .le. n)
    if (steps .ge. max_steps) exit
    steps = steps + 1

    ! save old state for interpolation
    x_old   = x
    tau_old = tau2
    rho_old = rho2
    eps_old = eps24
    T_old   = T_keV_

    ! temperature from tau-table at current tau (before step)
    T_keV = find_H_ordered_inc(tau_old, mas_tau_TkeV, n_mas)
    if (T_keV .le. T_floor) T_keV = T_floor

    ! evolve beta
    if (beta .gt. 0.d0) then
      dbeta_dx = 3.24d-4 * rho2 * Z**2 / max(beta,beta_floor)**3 * ln_Lambda
      beta = max(0.d0, beta - dx*dbeta_dx)
    else
      dbeta_dx = 0.d0
      beta = 0.d0
    end if

    ! evolve rho2
    if (beta .gt. 0.d0) then
      drho2 = ( 0.5d0*0.1d0*g14/T_keV*rho2 + 31.15d0*dot_m_6*0.5d0/T_keV*dbeta_dx ) * dx
    else
      drho2 = ( 0.5d0*0.1d0*g14/T_keV*rho2 ) * dx
    end if
    rho2 = max(rho_floor, rho2 + drho2)

    ! local source + proxy temperature
    if (beta .gt. 0.d0) then
      eps24  = 9.d2 * dot_m_6 * beta * dbeta_dx
      T_keV_ = ( eps24 / 1.75d0 / max(rho2,rho_floor)**2 )**2
    else
      eps24  = 0.d0
      T_keV_ = 0.d0
    end if

    ! evolve optical depth and x
    tau2 = tau2 + 0.34d0 * rho2 * dx
    x    = x + dx

    ! write any crossed tau targets (may cross several in one dx step)
    do while ( (tau2 .ge. target_tau) .and. (i .le. n) )
      if (tau2 .gt. tau_old) then
        f = (target_tau - tau_old) / (tau2 - tau_old)
      else
        f = 0.d0
      end if

      mas_x_rho_tau_2(i,1) = x_old + f*dx
      mas_x_rho_tau_2(i,2) = rho_old + f*(rho2 - rho_old)
      mas_x_rho_tau_2(i,3) = target_tau
      mas_x_rho_tau_2(i,4) = eps_old + f*(eps24 - eps_old)
      mas_x_rho_tau_2(i,5) = T_old   + f*(T_keV_ - T_old)

      i = i + 1
      if (i .le. n) then
        target_tau = tau_targets(i)
      end if
    end do
  end do

  ! pad if stopped early
  if (i .le. n) then
    do while (i .le. n)
      mas_x_rho_tau_2(i,1) = x
      mas_x_rho_tau_2(i,2) = rho2
      mas_x_rho_tau_2(i,3) = tau_targets(i)
      mas_x_rho_tau_2(i,4) = eps24
      mas_x_rho_tau_2(i,5) = T_keV_
      i = i + 1
    end do
  end if

  x_scale = mas_x_rho_tau_2(n,1) - mas_x_rho_tau_2(1,1)

  ! cumulative source distribution vs tau (normalized)
  F_tau(1,1) = 0.d0
  F_tau(1,2) = 0.d0
  do i = 2, n
    F_tau(i,1) = mas_x_rho_tau_2(i,3) - mas_x_rho_tau_2(1,3)
    F_tau(i,2) = F_tau(i-1,2) + mas_x_rho_tau_2(i,4) * (mas_x_rho_tau_2(i,1) - mas_x_rho_tau_2(i-1,1))
  end do
  if (F_tau(n,2) .gt. 0.d0) F_tau(1:n,2) = F_tau(1:n,2) / F_tau(n,2)

  return
end subroutine acc_atm_structure_4
!============================================================================================



!====================================================================================================
!====================================================================================================
subroutine acc_atm_source_F_distrib(F_tau_eps,mas_rho,n,rho_min,rho_max,m,R6,dot_m_6,mas_tau_TkeV,n_mas,ln_Lambda)
implicit none
real*8::F_tau_eps(n,2),mas_rho(n_mas)
integer,intent(in)::n,n_mas
real*8,intent(in)::rho_min,rho_max,m,R6,dot_m_6,mas_tau_TkeV(n_mas,2),ln_Lambda
real*8::g14,mas_x_rho_tau_2(n,5),x_scale,mas_tau_rho(n,2)
real*8::find_H_ordered_inc !== function ==!
integer::i
  g14 = 1.328d0*m/R6**2
  call acc_atm_structure_3(mas_x_rho_tau_2,n,rho_min,rho_max,x_scale,g14,mas_tau_TkeV,n_mas,dot_m_6,ln_Lambda,F_tau_eps)
  mas_tau_rho(1:n,1)=mas_x_rho_tau_2(1:n,3)
  mas_tau_rho(1:n,2)=mas_x_rho_tau_2(1:n,2)
  i=1
  do while(i.le.n_mas)
    mas_rho(i) = find_H_ordered_inc(mas_tau_TkeV(i,1),mas_tau_rho,n)
    i=i+1
  end do
return
end subroutine acc_atm_source_F_distrib
!====================================================================================================
