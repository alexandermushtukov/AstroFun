module xray_disk_mc_module
  !! Monte-Carlo X-ray disk irradiation kernel ported from the provided Python script.
  !!
  !! Public API:
  !!   call xray_disk_mc_profile(..., profile, ierr)
  !!
  !! Output profile(i,1:7):
  !!   1  R_cm
  !!   2  F_viscous_erg_s_cm2
  !!   3  F_direct_irradiation_erg_s_cm2
  !!   4  F_scattered_irradiation_erg_s_cm2
  !!   5  T_viscous_K
  !!   6  T_direct_irradiation_K
  !!   7  T_scattered_irradiation_K
  !!
  !! No plotting, no file I/O, no terminal output.

  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  private

  public :: xray_disk_mc_profile,test_xray_disk_mc_profile

  integer, parameter :: dp = real64

  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: mec2_kev = 511.0_dp
  real(dp), parameter :: sigma_sb_cgs = 5.670374419e-5_dp
  real(dp), parameter :: g_cgs = 6.67430e-8_dp
  real(dp), parameter :: m_sun_g = 1.98847e33_dp

  integer, parameter :: max_steps_default = 1000000
  integer, parameter :: max_internal_scatters_default = 100000
  integer, parameter :: planck_nmax = 800

  real(dp), save :: planck_cdf(planck_nmax)
  logical, save :: planck_ready = .false.

  type :: photon_result
     logical :: disk_absorbed = .false.
     integer :: disk_component = 0     ! 0 none, 1 direct, 2 scattered
     real(dp) :: disk_rho_model = 0.0_dp
     real(dp) :: launch_weight = 0.0_dp
     real(dp) :: e0_kev = 0.0_dp
     real(dp) :: efinal_kev = 0.0_dp
  end type photon_result

contains

  subroutine xray_disk_mc_profile(N, L39, B12, T_keV, nbins, profile, ierr, &
       r_g, R_ns, R_ns_km, tau_crit, r_stop_factor, step_coarse, tol, advance_factor, &
       L_x_erg_s, alpha_visc, m_msun, disk_rho_max_factor, disk_rho_min_factor, &
       forced_thin_scattering, min_packet_weight, seed, use_two_sided_disk_area)
    !! One-call production subroutine.
    !! Caller allocates profile(nbins,7). All units match the Python script.

    integer, intent(in) :: N, nbins
    real(dp), intent(in) :: L39, B12, T_keV
    real(dp), intent(out) :: profile(nbins, 7)
    integer, intent(out) :: ierr

    real(dp), intent(in), optional :: r_g, R_ns, R_ns_km, tau_crit, r_stop_factor
    real(dp), intent(in), optional :: step_coarse, tol, advance_factor
    real(dp), intent(in), optional :: L_x_erg_s, alpha_visc, m_msun
    real(dp), intent(in), optional :: disk_rho_max_factor, disk_rho_min_factor
    real(dp), intent(in), optional :: min_packet_weight
    logical, intent(in), optional :: forced_thin_scattering, use_two_sided_disk_area
    integer, intent(in), optional :: seed

    real(dp) :: rg, Rns, Rnskm, taucrit, rstopfac, stepc, tolerance, advfac
    real(dp) :: Lx, alpha, mstar, rho_max_fac, rho_min_fac, minw
    logical :: forced, two_sided

    real(dp) :: Rm, rho_min_model, rho_max_model
    real(dp) :: rho_edges(nbins+1), rho_centers(nbins)
    real(dp) :: R_edges_cm(nbins+1), R_centers_cm(nbins), dR_cm(nbins)
    real(dp) :: direct_bins(nbins), scattered_bins(nbins)
    real(dp) :: total_wE0, total_wEdirect, total_wEscattered
    real(dp) :: direct_dL(nbins), scattered_dL(nbins)
    real(dp) :: F_direct(nbins), F_scattered(nbins), T_direct(nbins), T_scattered(nbins)
    real(dp) :: F_visc(nbins), T_visc(nbins)
    real(dp) :: area(nbins), Rns_cm, Rinner_cm
    real(dp) :: rawsum, Ldir, Lscat
    integer :: i, ibin
    type(photon_result) :: res

    ierr = 0
    profile = 0.0_dp

    if (N <= 0 .or. nbins <= 0 .or. L39 <= 0.0_dp .or. B12 <= 0.0_dp .or. T_keV <= 0.0_dp) then
       ierr = 1
       return
    end if

    rg = opt_real(r_g, 0.42_dp)
    Rns = opt_real(R_ns, 1.0_dp)
    Rnskm = opt_real(R_ns_km, 10.0_dp)
    taucrit = opt_real(tau_crit, 0.5_dp)
    rstopfac = opt_real(r_stop_factor, 100.0_dp)
    stepc = opt_real(step_coarse, 0.5_dp)
    tolerance = opt_real(tol, 1.0e-5_dp)
    advfac = opt_real(advance_factor, 50.0_dp)
    Lx = opt_real(L_x_erg_s, L39 * 1.0e39_dp)
    alpha = opt_real(alpha_visc, 0.1_dp)
    mstar = opt_real(m_msun, 1.4_dp)
    rho_max_fac = opt_real(disk_rho_max_factor, 100.0_dp)
    rho_min_fac = opt_real(disk_rho_min_factor, 1.0_dp)
    minw = opt_real(min_packet_weight, 1.0e-12_dp)
    forced = opt_logical(forced_thin_scattering, .false.)
    two_sided = opt_logical(use_two_sided_disk_area, .true.)

    call init_planck_cdf()
    if (present(seed)) call seed_rng(seed)

    Rm = Rm_from_LB(L39, B12)
    rho_min_model = max(rho_min_fac * Rm, Rm, 1.0e-12_dp)
    rho_max_model = rho_max_fac * Rm
    if (rho_max_model <= rho_min_model) rho_max_model = rho_min_model * (1.0_dp + 1.0e-6_dp)

    do i = 1, nbins + 1
       rho_edges(i) = 10.0_dp ** (log10(rho_min_model) + real(i-1,dp) * &
            (log10(rho_max_model) - log10(rho_min_model)) / real(nbins,dp))
    end do
    do i = 1, nbins
       rho_centers(i) = sqrt(rho_edges(i) * rho_edges(i+1))
    end do

    direct_bins = 0.0_dp
    scattered_bins = 0.0_dp
    total_wE0 = 0.0_dp
    total_wEdirect = 0.0_dp
    total_wEscattered = 0.0_dp

    !$omp parallel default(shared) private(i,res,ibin) &
    !$omp reduction(+:total_wE0,total_wEdirect,total_wEscattered,direct_bins,scattered_bins)
    !$omp do schedule(dynamic)
    do i = 1, N
       call propagate_one_photon(L39, B12, T_keV, rg, Rns, Rnskm, taucrit, rstopfac, &
            stepc, tolerance, advfac, Lx, alpha, mstar, rho_max_fac, forced, minw, res)

       total_wE0 = total_wE0 + res%launch_weight * res%e0_kev

       if (res%disk_absorbed) then
          ibin = locate_bin(rho_edges, nbins, res%disk_rho_model)
          if (res%disk_component == 1) then
             total_wEdirect = total_wEdirect + res%launch_weight * res%efinal_kev
             if (ibin > 0) direct_bins(ibin) = direct_bins(ibin) + res%launch_weight * res%efinal_kev
          else if (res%disk_component == 2) then
             total_wEscattered = total_wEscattered + res%launch_weight * res%efinal_kev
             if (ibin > 0) scattered_bins(ibin) = scattered_bins(ibin) + res%launch_weight * res%efinal_kev
          end if
       end if
    end do
    !$omp end do
    !$omp end parallel

    Rns_cm = Rnskm * 1.0e5_dp
    R_edges_cm = rho_edges * Rns_cm
    R_centers_cm = rho_centers * Rns_cm
    dR_cm = R_edges_cm(2:nbins+1) - R_edges_cm(1:nbins)

    if (two_sided) then
       area = 4.0_dp * pi * R_centers_cm * dR_cm
    else
       area = 2.0_dp * pi * R_centers_cm * dR_cm
    end if

    direct_dL = 0.0_dp
    scattered_dL = 0.0_dp
    if (total_wE0 > 0.0_dp) then
       rawsum = sum(direct_bins)
       Ldir = Lx * total_wEdirect / total_wE0
       if (rawsum > 0.0_dp) direct_dL = Ldir * direct_bins / rawsum

       rawsum = sum(scattered_bins)
       Lscat = Lx * total_wEscattered / total_wE0
       if (rawsum > 0.0_dp) scattered_dL = Lscat * scattered_bins / rawsum
    end if

    F_direct = 0.0_dp
    F_scattered = 0.0_dp
    T_direct = 0.0_dp
    T_scattered = 0.0_dp
    do i = 1, nbins
       if (area(i) > 0.0_dp) then
          F_direct(i) = direct_dL(i) / area(i)
          F_scattered(i) = scattered_dL(i) / area(i)
       end if
       if (F_direct(i) > 0.0_dp) T_direct(i) = (F_direct(i) / sigma_sb_cgs) ** 0.25_dp
       if (F_scattered(i) > 0.0_dp) T_scattered(i) = (F_scattered(i) / sigma_sb_cgs) ** 0.25_dp
    end do

    Rinner_cm = Rm * Rns_cm
    call viscous_disk_profile(R_centers_cm, nbins, Lx, Rnskm, mstar, Rinner_cm, F_visc, T_visc)

    do i = 1, nbins
       profile(i,1) = R_centers_cm(i)
       profile(i,2) = F_visc(i)
       profile(i,3) = F_direct(i)
       profile(i,4) = F_scattered(i)
       profile(i,5) = T_visc(i)
       profile(i,6) = T_direct(i)
       profile(i,7) = T_scattered(i)
    end do
  end subroutine xray_disk_mc_profile

  subroutine propagate_one_photon(L39, B12, T_keV, r_g, R_ns, R_ns_km, tau_crit, r_stop_factor, &
       step_coarse, tol, advance_factor, Lx, alpha, mstar, disk_rho_max_factor, forced, minw, res)
    real(dp), intent(in) :: L39, B12, T_keV, r_g, R_ns, R_ns_km, tau_crit, r_stop_factor
    real(dp), intent(in) :: step_coarse, tol, advance_factor, Lx, alpha, mstar
    real(dp), intent(in) :: disk_rho_max_factor, minw
    logical, intent(in) :: forced
    type(photon_result), intent(out) :: res

    real(dp) :: Rm, r_stop, disk_rho_max_model, advance_dist
    real(dp) :: pos(3), normal_out(3), dir(3), hit(3), nnew(3), layer_normal(3)
    real(dp) :: E0, E, weight, mu_launch, t_stop, t_ns, t_m, t_disk, t_event, tau, pscat
    real(dp) :: disk_rho, side_sign, E_before, e_dep, offset_sign
    integer :: step, event, n_scat, n_internal, entry_sign, exit_side, status
    logical :: ok_stop, ok_ns, ok_m, ok_disk, did_scatter

    res = photon_result()
    Rm = Rm_from_LB(L39, B12)
    r_stop = r_stop_factor * Rm
    disk_rho_max_model = disk_rho_max_factor * Rm
    advance_dist = max(advance_factor * tol, 1.0e-7_dp)

    call sample_emission_point_on_polar_ring(Rm, R_ns, pos)
    normal_out = normalize(pos)
    E0 = sample_thermal_energy(T_keV)
    E = E0
    call sample_uniform_outward_hemisphere_direction(normal_out, dir, mu_launch)
    weight = mu_launch
    n_scat = 0

    do step = 1, max_steps_default
       if (weight < minw) exit

       ok_stop = ray_sphere_first_intersection(pos, dir, r_stop, advance_dist, t_stop)
       ok_ns   = ray_sphere_first_intersection(pos, dir, R_ns,   advance_dist, t_ns)
       ok_m    = ray_dipole_first_intersection(pos, dir, Rm, step_coarse, tol, advance_dist, merge(t_stop, huge(1.0_dp), ok_stop), t_m)
       ok_disk = ray_flared_disk_first_intersection(pos, dir, Rm, Lx, R_ns_km, alpha, mstar, &
            advance_dist, merge(t_stop, huge(1.0_dp), ok_stop), step_coarse, tol, disk_rho_max_model, t_disk, disk_rho, side_sign)

       if (.not. ok_stop .and. .not. ok_ns .and. .not. ok_m .and. .not. ok_disk) exit

       call choose_event(ok_stop,t_stop, ok_ns,t_ns, ok_m,t_m, ok_disk,t_disk, t_event, event)
       hit = pos + t_event * dir

       select case(event)
       case(1) ! outer boundary
          exit

       case(2) ! disk
          res%disk_absorbed = .true.
          if (n_scat == 0) then
             res%disk_component = 1
          else
             res%disk_component = 2
          end if
          res%disk_rho_model = disk_rho
          exit

       case(3) ! neutron-star reflection
          pos = R_ns * normalize(hit)
          call sample_uniform_outward_hemisphere_direction(normalize(pos), dir, mu_launch)
          weight = weight * mu_launch
          pos = pos + advance_dist * dir

       case(4) ! magnetosphere
          tau = tau_at_point(hit, Rm, L39, B12, r_g, R_ns)
          if (tau <= tau_crit) then
             pscat = p_scatter_from_tau(tau)
             did_scatter = .true.
             if (.not. forced) then
                did_scatter = (urand() < pscat)
             else
                weight = weight * pscat
                did_scatter = (pscat > 0.0_dp .and. weight >= minw)
             end if

             if (did_scatter) then
                E_before = E
                call isotropic_random_direction(nnew)
                E = compton_scatter_energy(E, dot3(dir, nnew))
                dir = normalize(nnew)
                n_scat = n_scat + 1
                pos = hit + advance_dist * dir
             else
                pos = hit + advance_dist * dir
             end if
          else
             layer_normal = dipole_surface_normal(hit, Rm)
             call transport_through_layer(tau, dir, E, layer_normal, nnew, E, n_internal, e_dep, entry_sign, exit_side, status)
             dir = normalize(nnew)
             n_scat = n_scat + n_internal
             if (exit_side == 2) then
                offset_sign = real(entry_sign,dp)
             else
                offset_sign = -real(entry_sign,dp)
             end if
             pos = hit + offset_sign * advance_dist * layer_normal
             if (status /= 0) exit
          end if
       end select
    end do

    res%launch_weight = weight
    res%e0_kev = E0
    res%efinal_kev = E
  end subroutine propagate_one_photon

  subroutine choose_event(os,ts,on,tn,om,tm,od,td,t,event)
    logical, intent(in) :: os,on,om,od
    real(dp), intent(in) :: ts,tn,tm,td
    real(dp), intent(out) :: t
    integer, intent(out) :: event
    t = huge(1.0_dp); event = 0
    if (os .and. ts < t) then; t = ts; event = 1; end if
    if (od .and. td < t) then; t = td; event = 2; end if
    if (on .and. tn < t) then; t = tn; event = 3; end if
    if (om .and. tm < t) then; t = tm; event = 4; end if
  end subroutine choose_event

  function locate_bin(edges, nbins, x) result(idx)
    integer, intent(in) :: nbins
    real(dp), intent(in) :: edges(nbins+1), x
    integer :: idx, i
    idx = 0
    do i = 1, nbins
       if ((x >= edges(i) .and. x < edges(i+1)) .or. (i == nbins .and. x == edges(i+1))) then
          idx = i
          return
       end if
    end do
  end function locate_bin

  function opt_real(x, default) result(y)
    real(dp), intent(in), optional :: x
    real(dp), intent(in) :: default
    real(dp) :: y
    if (present(x)) then; y = x; else; y = default; end if
  end function opt_real

  function opt_logical(x, default) result(y)
    logical, intent(in), optional :: x
    logical, intent(in) :: default
    logical :: y
    if (present(x)) then; y = x; else; y = default; end if
  end function opt_logical

  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
       put(i) = seed + 104729 * (i - 1)
    end do
    call random_seed(put=put)
    deallocate(put)
  end subroutine seed_rng

  function urand() result(u)
    real(dp) :: u
    call random_number(u)
    u = max(u, 1.0e-300_dp)
  end function urand

  function dot3(a,b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c
    c = sum(a*b)
  end function dot3

  function norm3(a) result(n)
    real(dp), intent(in) :: a(3)
    real(dp) :: n
    n = sqrt(sum(a*a))
  end function norm3

  function normalize(a) result(b)
    real(dp), intent(in) :: a(3)
    real(dp) :: b(3), n
    n = norm3(a)
    if (n <= 0.0_dp) then
       b = [1.0_dp, 0.0_dp, 0.0_dp]
    else
       b = a / n
    end if
  end function normalize

  function cross3(a,b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
  end function cross3

  subroutine isotropic_random_direction(n)
    real(dp), intent(out) :: n(3)
    real(dp) :: u, v, mu, phi, st
    u = urand(); v = urand()
    mu = 1.0_dp - 2.0_dp*u
    phi = 2.0_dp*pi*v
    st = sqrt(max(0.0_dp, 1.0_dp - mu*mu))
    n = [st*cos(phi), st*sin(phi), mu]
  end subroutine isotropic_random_direction

  subroutine orthonormal_basis(normal_out, e1, e2, e3)
    real(dp), intent(in) :: normal_out(3)
    real(dp), intent(out) :: e1(3), e2(3), e3(3)
    real(dp) :: ref(3)
    e3 = normalize(normal_out)
    if (abs(e3(3)) < 0.9_dp) then
       ref = [0.0_dp,0.0_dp,1.0_dp]
    else
       ref = [1.0_dp,0.0_dp,0.0_dp]
    end if
    e1 = cross3(ref, e3)
    if (norm3(e1) < 1.0e-14_dp) e1 = cross3([0.0_dp,1.0_dp,0.0_dp], e3)
    e1 = normalize(e1)
    e2 = normalize(cross3(e3,e1))
  end subroutine orthonormal_basis

  subroutine sample_uniform_outward_hemisphere_direction(normal_out, direction, mu_check)
    real(dp), intent(in) :: normal_out(3)
    real(dp), intent(out) :: direction(3), mu_check
    real(dp) :: e1(3), e2(3), e3(3), mu, phi, st
    call orthonormal_basis(normal_out, e1, e2, e3)
    mu = urand()
    phi = 2.0_dp*pi*urand()
    st = sqrt(max(0.0_dp, 1.0_dp - mu*mu))
    direction = st*cos(phi)*e1 + st*sin(phi)*e2 + mu*e3
    direction = normalize(direction)
    mu_check = max(0.0_dp, dot3(direction, normalize(normal_out)))
  end subroutine sample_uniform_outward_hemisphere_direction

  function F_dipole_xyz(x,y,z,Rm) result(f)
    real(dp), intent(in) :: x,y,z,Rm
    real(dp) :: f, r2, r
    r2 = x*x + y*y + z*z
    r = sqrt(r2)
    if (r == 0.0_dp) then
       f = huge(1.0_dp)
    else
       f = r*r*r - Rm*r2 + Rm*z*z
    end if
  end function F_dipole_xyz

  function F_along_ray(t,p,n,Rm) result(f)
    real(dp), intent(in) :: t,p(3),n(3),Rm
    real(dp) :: f, q(3)
    q = p + t*n
    f = F_dipole_xyz(q(1),q(2),q(3),Rm)
  end function F_along_ray

  function adaptive_step(p,n,t,step_min) result(dt)
    real(dp), intent(in) :: p(3), n(3), t, step_min
    real(dp) :: dt, q(3), r
    q = p + t*n
    r = norm3(q)
    dt = min(max(5.0_dp, step_min), max(step_min, 0.005_dp * max(r, 1.0_dp)))
  end function adaptive_step

  function bisect_dipole_root(p,n,Rm,a,b,fa,fb,tol) result(root)
    real(dp), intent(in) :: p(3), n(3), Rm, a, b, fa, fb, tol
    real(dp) :: root, lo, hi, flo, mid, fmid
    integer :: i
    lo = a; hi = b; flo = fa
    do i = 1, 200
       mid = 0.5_dp * (lo + hi)
       fmid = F_along_ray(mid,p,n,Rm)
       if ((hi-lo) <= tol .or. fmid == 0.0_dp) exit
       if (flo*fmid <= 0.0_dp) then
          hi = mid
       else
          lo = mid; flo = fmid
       end if
    end do
    root = 0.5_dp * (lo + hi)
  end function bisect_dipole_root

  function ray_dipole_first_intersection(p0,n_dir,Rm,step_coarse,tol,t_min,t_max,t_root) result(ok)
    real(dp), intent(in) :: p0(3), n_dir(3), Rm, step_coarse, tol, t_min, t_max
    real(dp), intent(out) :: t_root
    logical :: ok
    real(dp) :: n(3), tprev, fprev, t, fcur, dt, root
    n = normalize(n_dir)
    ok = .false.; t_root = 0.0_dp
    tprev = max(0.0_dp, t_min)
    fprev = F_along_ray(tprev,p0,n,Rm)
    if (abs(fprev) < 1.0e-14_dp) then
       tprev = tprev + max(1.0e-7_dp, 0.5_dp*tol)
       fprev = F_along_ray(tprev,p0,n,Rm)
    end if
    do while (tprev < t_max - 1.0e-15_dp)
       dt = adaptive_step(p0,n,tprev,step_coarse)
       t = min(tprev + dt, t_max)
       fcur = F_along_ray(t,p0,n,Rm)
       if (fprev == 0.0_dp .or. fcur == 0.0_dp .or. fprev*fcur < 0.0_dp) then
          if (fprev == 0.0_dp) then
             root = tprev
          else if (fcur == 0.0_dp) then
             root = t
          else
             root = bisect_dipole_root(p0,n,Rm,tprev,t,fprev,fcur,tol)
          end if
          if (root >= t_min .and. root > 0.0_dp) then
             t_root = root; ok = .true.; return
          end if
       end if
       tprev = t; fprev = fcur
    end do
  end function ray_dipole_first_intersection

  function ray_sphere_first_intersection(p0,n_dir,R,t_min,t_hit) result(ok)
    real(dp), intent(in) :: p0(3), n_dir(3), R, t_min
    real(dp), intent(out) :: t_hit
    logical :: ok
    real(dp) :: n(3), b, c, disc, s, t1, t2
    n = normalize(n_dir)
    b = 2.0_dp * dot3(p0,n)
    c = dot3(p0,p0) - R*R
    disc = b*b - 4.0_dp*c
    ok = .false.; t_hit = 0.0_dp
    if (disc < 0.0_dp) return
    s = sqrt(disc)
    t1 = (-b - s) / 2.0_dp
    t2 = (-b + s) / 2.0_dp
    if (t1 > t_min .and. t2 > t_min) then
       t_hit = min(t1,t2); ok = .true.
    else if (t1 > t_min) then
       t_hit = t1; ok = .true.
    else if (t2 > t_min) then
       t_hit = t2; ok = .true.
    end if
  end function ray_sphere_first_intersection

  function disk_H_coeff(Lx,Rnskm,alpha,mstar) result(c)
    real(dp), intent(in) :: Lx, Rnskm, alpha, mstar
    real(dp) :: c, L37, Rns_cm, R6
    L37 = Lx / 1.0e37_dp
    Rns_cm = Rnskm * 1.0e5_dp
    R6 = Rns_cm / 1.0e6_dp
    c = 0.08_dp * alpha**(-0.1_dp) * L37**(3.0_dp/20.0_dp) * &
        mstar**(-21.0_dp/40.0_dp) * R6**(3.0_dp/20.0_dp) * &
        (Rns_cm/1.0e8_dp)**(1.0_dp/8.0_dp)
  end function disk_H_coeff

  function disk_H_model(rho,Lx,Rnskm,alpha,mstar) result(H)
    real(dp), intent(in) :: rho,Lx,Rnskm,alpha,mstar
    real(dp) :: H
    if (rho <= 0.0_dp) then
       H = 0.0_dp
    else
       H = disk_H_coeff(Lx,Rnskm,alpha,mstar) * rho**(9.0_dp/8.0_dp)
    end if
  end function disk_H_model

  function disk_g(t,p,n,side,Lx,Rnskm,alpha,mstar) result(g)
    real(dp), intent(in) :: t,p(3),n(3),side,Lx,Rnskm,alpha,mstar
    real(dp) :: g, q(3), rho, H
    q = p + t*n
    rho = sqrt(q(1)*q(1) + q(2)*q(2))
    H = disk_H_model(rho,Lx,Rnskm,alpha,mstar)
    g = q(3) - side*H
  end function disk_g

  function bisect_disk_root(p,n,side,a,b,fa,fb,Lx,Rnskm,alpha,mstar,tol) result(root)
    real(dp), intent(in) :: p(3),n(3),side,a,b,fa,fb,Lx,Rnskm,alpha,mstar,tol
    real(dp) :: root, lo, hi, flo, mid, fmid
    integer :: i
    lo=a; hi=b; flo=fa
    do i=1,100
       mid = 0.5_dp*(lo+hi)
       fmid = disk_g(mid,p,n,side,Lx,Rnskm,alpha,mstar)
       if ((hi-lo) <= tol .or. fmid == 0.0_dp) exit
       if (flo*fmid <= 0.0_dp) then
          hi=mid
       else
          lo=mid; flo=fmid
       end if
    end do
    root = 0.5_dp*(lo+hi)
  end function bisect_disk_root

  function ray_flared_disk_first_intersection(p0,n_dir,Rm,Lx,Rnskm,alpha,mstar,t_min,t_max,step_coarse,tol,rho_max_model,t_hit,rho_hit,side_hit) result(ok)
    real(dp), intent(in) :: p0(3), n_dir(3), Rm, Lx, Rnskm, alpha, mstar, t_min, t_max, step_coarse, tol, rho_max_model
    real(dp), intent(out) :: t_hit, rho_hit, side_hit
    logical :: ok
    real(dp) :: n(3), side, tprev, fprev, t, fcur, root, q(3), rho, dt
    real(dp) :: best_t, best_rho, best_side
    integer :: iside
    ok = .false.; best_t = huge(1.0_dp); best_rho = 0.0_dp; best_side = 1.0_dp
    n = normalize(n_dir)
    do iside = 1, 2
       side = merge(1.0_dp, -1.0_dp, iside == 1)
       tprev = max(0.0_dp, t_min)
       fprev = disk_g(tprev,p0,n,side,Lx,Rnskm,alpha,mstar)
       if (abs(fprev) < 1.0e-14_dp) then
          tprev = tprev + max(tol, 1.0e-7_dp)
          fprev = disk_g(tprev,p0,n,side,Lx,Rnskm,alpha,mstar)
       end if
       do while (tprev < t_max - 1.0e-15_dp)
          dt = adaptive_step(p0,n,tprev,step_coarse)
          t = min(tprev + dt, t_max)
          fcur = disk_g(t,p0,n,side,Lx,Rnskm,alpha,mstar)
          if (fprev == 0.0_dp .or. fcur == 0.0_dp .or. fprev*fcur < 0.0_dp) then
             if (fprev == 0.0_dp) then
                root = tprev
             else if (fcur == 0.0_dp) then
                root = t
             else
                root = bisect_disk_root(p0,n,side,tprev,t,fprev,fcur,Lx,Rnskm,alpha,mstar,tol)
             end if
             if (root > t_min) then
                q = p0 + root*n
                rho = sqrt(q(1)*q(1) + q(2)*q(2))
                if (rho + 1.0e-12_dp >= Rm .and. rho <= rho_max_model + 1.0e-12_dp) then
                   if (root < best_t) then
                      best_t = root; best_rho = rho; best_side = side; ok = .true.
                   end if
                   exit
                end if
             end if
          end if
          tprev = t; fprev = fcur
       end do
    end do
    t_hit = best_t; rho_hit = best_rho; side_hit = best_side
  end function ray_flared_disk_first_intersection

  function Rm_from_LB(L39,B12) result(Rm)
    real(dp), intent(in) :: L39, B12
    real(dp) :: Rm
    Rm = 35.0_dp * 1.4_dp**(1.0_dp/7.0_dp) * B12**(4.0_dp/7.0_dp) * L39**(-2.0_dp/7.0_dp)
  end function Rm_from_LB

  function beta_from_r(r,r_g) result(beta)
    real(dp), intent(in) :: r, r_g
    real(dp) :: beta
    if (r <= 0.0_dp) then
       beta = 0.0_dp
    else
       beta = min(sqrt(r_g/r), 1.0_dp)
    end if
  end function beta_from_r

  function theta_ring_from_Rm(Rm,Rns) result(theta)
    real(dp), intent(in) :: Rm, Rns
    real(dp) :: theta
    theta = sqrt(max(0.0_dp, Rns/Rm))
  end function theta_ring_from_Rm

  function tau_at_point(hit,Rm,L39,B12,r_g,Rns) result(tau)
    real(dp), intent(in) :: hit(3), Rm, L39, B12, r_g, Rns
    real(dp) :: tau, r, beta, cl, cl0, mu, pref
    r = norm3(hit)
    if (r <= 0.0_dp) then
       tau = huge(1.0_dp); return
    end if
    beta = max(beta_from_r(r,r_g), 1.0e-12_dp)
    mu = hit(3) / r
    cl = sqrt(max(0.0_dp, 1.0_dp - mu*mu))
    cl0 = sin(theta_ring_from_Rm(Rm,Rns))
    pref = 70.0_dp * L39**(6.0_dp/7.0_dp) * B12**(2.0_dp/7.0_dp) / beta
    tau = pref * (cl0 / max(cl, 1.0e-12_dp))**3
  end function tau_at_point

  function p_scatter_from_tau(tau) result(p)
    real(dp), intent(in) :: tau
    real(dp) :: p
    if (tau <= 0.0_dp) then
       p = 0.0_dp
    else
       p = 1.0_dp - exp(-tau)
    end if
  end function p_scatter_from_tau

  subroutine sample_emission_point_on_polar_ring(Rm,Rns,pos)
    real(dp), intent(in) :: Rm, Rns
    real(dp), intent(out) :: pos(3)
    real(dp) :: theta0, phi, zsign, st, ct
    theta0 = theta_ring_from_Rm(Rm,Rns)
    phi = 2.0_dp*pi*urand()
    zsign = merge(1.0_dp, -1.0_dp, urand() < 0.5_dp)
    st = sin(theta0); ct = cos(theta0)
    pos = [Rns*st*cos(phi), Rns*st*sin(phi), zsign*Rns*ct]
  end subroutine sample_emission_point_on_polar_ring

  subroutine init_planck_cdf()
    integer :: i
    real(dp) :: s
    if (planck_ready) return
    s = 0.0_dp
    do i = 1, planck_nmax
       s = s + 1.0_dp / real(i,dp)**3
       planck_cdf(i) = s
    end do
    planck_cdf = planck_cdf / planck_cdf(planck_nmax)
    planck_ready = .true.
  end subroutine init_planck_cdf

  function sample_thermal_energy(T_keV) result(E)
    real(dp), intent(in) :: T_keV
    real(dp) :: E, u, x
    integer :: n
    u = urand()
    n = 1
    do while (n < planck_nmax .and. planck_cdf(n) < u)
       n = n + 1
    end do
    x = gamma_k3_scale(1.0_dp / real(n,dp))
    E = T_keV * x
  end function sample_thermal_energy

  function gamma_k3_scale(scale) result(x)
    !! Gamma(shape=3, scale=scale), equivalent to sum of 3 exponentials.
    real(dp), intent(in) :: scale
    real(dp) :: x
    x = -scale * log(urand()*urand()*urand())
  end function gamma_k3_scale

  function compton_scatter_energy(E,cospsi) result(Eout)
    real(dp), intent(in) :: E, cospsi
    real(dp) :: Eout, c
    c = max(-1.0_dp, min(1.0_dp, cospsi))
    Eout = E / (1.0_dp + (E/mec2_kev)*(1.0_dp-c))
  end function compton_scatter_energy

  function dipole_surface_normal(hit,Rm) result(n)
    real(dp), intent(in) :: hit(3), Rm
    real(dp) :: n(3), r
    r = norm3(hit)
    n = [hit(1)*(3.0_dp*r - 2.0_dp*Rm), hit(2)*(3.0_dp*r - 2.0_dp*Rm), 3.0_dp*r*hit(3)]
    n = normalize(n)
  end function dipole_surface_normal

  subroutine transport_through_layer(tau0, direction_in, E_in, normal, n_out, E_out, n_internal, e_dep, entry_sign, exit_side, status)
    real(dp), intent(in) :: tau0, direction_in(3), E_in, normal(3)
    real(dp), intent(out) :: n_out(3), E_out, e_dep
    integer, intent(out) :: n_internal, entry_sign, exit_side, status

    real(dp) :: n(3), m(3), E, mu, u, tau_to_boundary, delta_tau, nnew(3)
    logical :: entered_from_low
    integer :: i, exit_boundary

    status = 0
    n = normalize(direction_in)
    m = normalize(normal)
    E = E_in
    mu = dot3(n,m)
    if (abs(mu) < 1.0e-14_dp) mu = merge(1.0e-14_dp, -1.0e-14_dp, urand() < 0.5_dp)
    entry_sign = merge(1, -1, mu > 0.0_dp)
    entered_from_low = (mu > 0.0_dp)
    u = merge(0.0_dp, tau0, mu > 0.0_dp)
    n_internal = 0

    do i = 1, max_internal_scatters_default
       mu = dot3(n,m)
       if (abs(mu) < 1.0e-14_dp) mu = merge(1.0e-14_dp, -1.0e-14_dp, urand() < 0.5_dp)
       if (mu > 0.0_dp) then
          tau_to_boundary = max((tau0-u)/mu, 0.0_dp)
       else
          tau_to_boundary = max(u/(-mu), 0.0_dp)
       end if
       delta_tau = -log(urand())
       if (delta_tau >= tau_to_boundary) then
          exit_boundary = merge(2, 1, mu > 0.0_dp) ! 2 high, 1 low
          if (entered_from_low) then
             exit_side = merge(1, 2, exit_boundary == 1) ! 1 same, 2 opposite
          else
             exit_side = merge(1, 2, exit_boundary == 2)
          end if
          n_out = n; E_out = E; e_dep = E_in - E
          return
       end if
       u = min(max(u + mu*delta_tau, 0.0_dp), tau0)
       call isotropic_random_direction(nnew)
       E = compton_scatter_energy(E, dot3(n,nnew))
       n = normalize(nnew)
       n_internal = n_internal + 1
    end do

    status = 1
    exit_side = 1
    n_out = n; E_out = E; e_dep = E_in - E
  end subroutine transport_through_layer

  subroutine viscous_disk_profile(R_cm, n, Lx, Rnskm, mstar, Rinner, F, T)
    integer, intent(in) :: n
    real(dp), intent(in) :: R_cm(n), Lx, Rnskm, mstar, Rinner
    real(dp), intent(out) :: F(n), T(n)
    real(dp) :: Rns_cm, M_g, mdot, boundary
    integer :: i
    Rns_cm = Rnskm * 1.0e5_dp
    M_g = mstar * m_sun_g
    mdot = Lx * Rns_cm / (g_cgs * M_g)
    F = 0.0_dp; T = 0.0_dp
    do i = 1, n
       if (R_cm(i) > 0.0_dp .and. R_cm(i) >= Rinner) then
          boundary = max(1.0_dp - sqrt(Rinner / R_cm(i)), 0.0_dp)
          F(i) = 3.0_dp * g_cgs * M_g * mdot / (8.0_dp*pi*R_cm(i)**3) * boundary
          if (F(i) > 0.0_dp) T(i) = (F(i)/sigma_sb_cgs)**0.25_dp
       end if
    end do
  end subroutine viscous_disk_profile


  subroutine test_xray_disk_mc_profile()
  implicit none
  integer, parameter :: nbins = 40
  integer :: ierr,i
  real(8) :: profile(nbins,7)
    call xray_disk_mc_profile( &
        N       = int(5.e5), &
        L39     = 0.001d0, &
        B12     = 3.0d0, &
        T_keV   = 5.0d0, &
        nbins   = nbins, &
        profile = profile, disk_rho_max_factor = 100.d0, &
        ierr    = ierr )
    if (ierr /= 0) stop 'xray_disk_mc_profile failed'
    open(unit=10, file='./res/xray_disk_profile_L1e36_B3e12', status='replace')
    do i = 1, nbins
      write(10,'(7E20.10)') profile(i,1), profile(i,2), profile(i,3), profile(i,4), &
                            profile(i,5), profile(i,6), profile(i,7)
    end do
    close(10)
  return
  end subroutine test_xray_disk_mc_profile
end module xray_disk_mc_module




!=================================================================================
!
!
!=================================================================================
subroutine test_AC()
use omp_lib
use iso_fortran_env, only: int64
implicit none
integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: G        = 6.67430d-8
  real(dp), parameter :: c_light  = 2.99792458d10
  real(dp), parameter :: sigma_T  = 6.6524587321d-25
  real(dp), parameter :: m_p      = 1.67262192369d-24
  real(dp), parameter :: k_B      = 1.380649d-16
  real(dp), parameter :: a_rad    = 7.5657d-15
  real(dp), parameter :: M_sun    = 1.98847d33
  real(dp), parameter :: erg_keV  = 1.602176634d-9
  real(dp), parameter :: pi       = 3.1415926535897932384626433832795d0

  real(dp), parameter :: M_NS = 1.4d0 * M_sun
  real(dp), parameter :: R_NS = 1.0d6

  real(dp), parameter :: R_COL = 1.0d4
  real(dp), parameter :: H_COL = 2.0d5

  real(dp), parameter :: MDOT_COL = 7.0d16

  integer, parameter :: NR = 40
  integer, parameter :: NZ = 80

  integer, parameter :: N_PACKETS = 5000
  integer, parameter :: N_ITER    = 50

  real(dp), parameter :: OPACITY_SCALE = 1.0d0

  integer, parameter :: MAX_SCATTERS = 500000
  integer, parameter :: MAX_STEPS    = 300000

  real(dp), parameter :: V_MIN = 1.0d6
  real(dp), parameter :: VELOCITY_RELAX = 0.03d0

  integer, parameter :: N_STREAMS = 8

  real(dp) :: r_edges(0:NR), z_edges(0:NZ)
  real(dp) :: r_cent(NR), z_cent(NZ)
  real(dp) :: dr, dz

  real(dp) :: cell_vol(NR,NZ)
  real(dp) :: ring_area(NR)
  real(dp) :: side_area_zbin

  real(dp) :: R2D(NR,NZ), Z2D(NR,NZ)
  real(dp) :: g_grid(NR,NZ)
  real(dp) :: v_ff_profile(NZ)
  real(dp) :: v_ff_grid(NR,NZ)

  real(dp) :: v_grid(NR,NZ), v_new(NR,NZ), v_target(NR,NZ)
  real(dp) :: rho_grid(NR,NZ), ne_grid(NR,NZ), alpha_grid(NR,NZ)
  real(dp) :: dv_dz_grid(NR,NZ)

  real(dp) :: q_grid(NR,NZ), lum_cells(NR,NZ)
  real(dp) :: u_accum(NR,NZ), Fr_accum(NR,NZ), Fz_accum(NR,NZ)
  real(dp) :: u_rad(NR,NZ), Fr(NR,NZ), Fz(NR,NZ)
  real(dp) :: T_K(NR,NZ), T_keV(NR,NZ)
  real(dp) :: P_gas(NR,NZ), P_rad(NR,NZ)
  real(dp) :: dPgas_dz(NR,NZ), dPrad_dz(NR,NZ)
  real(dp) :: dv_dz_suggested(NR,NZ)

  real(dp), allocatable :: u_thr(:,:,:), Fr_thr(:,:,:), Fz_thr(:,:,:)
  real(dp), allocatable :: visit_thr(:,:,:), side_thr(:,:)

  real(dp) :: cell_cdf(NR*NZ)
  real(dp) :: base_ring_cdf(NR)
  real(dp) :: base_lum_rings(NR)

  real(dp) :: L_volume, L_base, L_total
  real(dp) :: packet_luminosity, p_base
  real(dp) :: tau_radial_mean, tau_vertical_mean

  real(dp) :: area_col
  real(dp) :: rho_base(NR), v_base(NR), mdot_ring(NR)

  integer :: i, j, k, iteration
  integer :: i_vprof

  integer(int64) :: seeds(N_STREAMS)

  integer :: scatter_hist(0:MAX_SCATTERS+1)
  integer, allocatable :: scatter_thr(:,:)

  real(dp) :: dv_dz_residual(NR,NZ)
  real(dp) :: iter_metric
  real(dp), parameter :: RES_EPS = 1.0d-30
  real(dp) :: dv_raw, dv_lim

  call execute_command_line("mkdir -p ./res")

  call omp_set_num_threads(N_STREAMS)

  scatter_hist = 0

  do i = 0, NR
     r_edges(i) = R_COL * real(i,dp) / real(NR,dp)
  end do

  do j = 0, NZ
     z_edges(j) = H_COL * real(j,dp) / real(NZ,dp)
  end do

  do i = 1, NR
     r_cent(i) = 0.5d0 * (r_edges(i-1) + r_edges(i))
  end do

  do j = 1, NZ
     z_cent(j) = 0.5d0 * (z_edges(j-1) + z_edges(j))
  end do

  dr = r_edges(1) - r_edges(0)
  dz = z_edges(1) - z_edges(0)

  do i = 1, NR
     ring_area(i) = pi * (r_edges(i)**2 - r_edges(i-1)**2)
  end do

  do i = 1, NR
     do j = 1, NZ
        cell_vol(i,j) = ring_area(i) * dz
        R2D(i,j) = r_cent(i)
        Z2D(i,j) = z_cent(j)
     end do
  end do

  side_area_zbin = 2.0d0 * pi * R_COL * dz

  i_vprof = 1
  do i = 1, NR
     if (abs(r_cent(i) - R_COL/3.0d0) < abs(r_cent(i_vprof) - R_COL/3.0d0)) then
        i_vprof = i
     end if
  end do

  do j = 1, NZ
     v_ff_profile(j) = v_freefall(z_cent(j))
  end do

  do i = 1, NR
     do j = 1, NZ
        g_grid(i,j) = gravity(z_cent(j))
        v_ff_grid(i,j) = v_ff_profile(j)
        v_grid(i,j) = v_ff_grid(i,j)
     end do
  end do

  do i = 1, N_STREAMS
     seeds(i) = int(12345 + 104729*i, int64)
  end do

  allocate(u_thr(NR,NZ,N_STREAMS))
  allocate(Fr_thr(NR,NZ,N_STREAMS))
  allocate(Fz_thr(NR,NZ,N_STREAMS))
  allocate(visit_thr(NR,NZ,N_STREAMS))
  allocate(side_thr(NZ,N_STREAMS))
  allocate(scatter_thr(0:MAX_SCATTERS+1,N_STREAMS))

  dv_dz_residual = 0.0d0
  call write_velocity_profile(0)
  do iteration = 1, N_ITER

     write(*,*)
     write(*,*) "================================================"
     write(*,'(A,I4,A,I4)') "ITERATION ", iteration, " / ", N_ITER
     write(*,*) "================================================"

     area_col = pi * R_COL**2

     do i = 1, NR
        do j = 1, NZ
           rho_grid(i,j)   = MDOT_COL / (area_col * v_grid(i,j))
           ne_grid(i,j)    = rho_grid(i,j) / m_p
           alpha_grid(i,j) = OPACITY_SCALE * ne_grid(i,j) * sigma_T
        end do
     end do

     tau_radial_mean   = sum(alpha_grid) / real(NR*NZ,dp) * R_COL
     tau_vertical_mean = sum(alpha_grid) / real(NR*NZ,dp) * H_COL

     write(*,'(A,ES12.4)') "Mean radial optical depth   = ", tau_radial_mean
     write(*,'(A,ES12.4)') "Mean vertical optical depth = ", tau_vertical_mean

     call gradient_z(v_grid, dv_dz_grid, dz)

     do i = 1, NR
        do j = 1, NZ
           q_grid(i,j) = rho_grid(i,j) * v_grid(i,j) * &
                         (g_grid(i,j) + v_grid(i,j) * dv_dz_grid(i,j))
           if (q_grid(i,j) < 0.0d0) q_grid(i,j) = 0.0d0
           lum_cells(i,j) = q_grid(i,j) * cell_vol(i,j)
        end do
     end do

     L_volume = sum(lum_cells)

     do i = 1, NR
        v_base(i) = v_grid(i,1)
        rho_base(i) = rho_grid(i,1)
        mdot_ring(i) = rho_base(i) * v_base(i) * ring_area(i)
        base_lum_rings(i) = 0.5d0 * mdot_ring(i) * v_base(i)**2
     end do

     L_base = sum(base_lum_rings)
     L_total = L_volume + L_base

     write(*,'(A,ES12.4)') "L_volume = ", L_volume
     write(*,'(A,ES12.4)') "L_base   = ", L_base
     write(*,'(A,ES12.4)') "L_total  = ", L_total

     if (L_total <= 0.0d0) stop "ERROR: total luminosity is non-positive."

     call build_cdf_2d(lum_cells, NR, NZ, cell_cdf)
     call build_cdf_1d(base_lum_rings, NR, base_ring_cdf)

     p_base = L_base / L_total
     packet_luminosity = L_total / real(N_PACKETS,dp)

     u_thr = 0.0d0
     Fr_thr = 0.0d0
     Fz_thr = 0.0d0
     visit_thr = 0.0d0
     side_thr = 0.0d0
     scatter_thr = 0

     !$omp parallel default(shared) private(k)
     call mc_parallel_loop()
     !$omp end parallel

     u_accum = sum(u_thr, dim=3)
     Fr_accum = sum(Fr_thr, dim=3)
     Fz_accum = sum(Fz_thr, dim=3)

     do i = 1, NR
        do j = 1, NZ
           u_rad(i,j) = u_accum(i,j) / cell_vol(i,j)
           Fr(i,j)    = Fr_accum(i,j) / cell_vol(i,j)
           Fz(i,j)    = Fz_accum(i,j) / cell_vol(i,j)
        end do
     end do

     Fr(1,:) = 0.0d0

     do i = 1, NR
        do j = 1, NZ
           T_K(i,j)   = max((u_rad(i,j)/a_rad)**0.25d0, 1.0d0)
           T_keV(i,j) = k_B * T_K(i,j) / erg_keV
           P_gas(i,j) = 2.0d0 * ne_grid(i,j) * k_B * T_K(i,j)
           P_rad(i,j) = u_rad(i,j) / 3.0d0
        end do
     end do

     call gradient_z(P_gas, dPgas_dz, dz)
     call gradient_z(P_rad, dPrad_dz, dz)

     do i = 1, NR
        do j = 1, NZ
           dv_dz_suggested(i,j) = - ( g_grid(i,j) + dPgas_dz(i,j)/rho_grid(i,j) &
                                      + dPrad_dz(i,j)/rho_grid(i,j) ) &
                                  / max(v_grid(i,j), 1.0d-30)
        end do
     end do

    ! ========================================================
    ! Iteration residual:
    ! current velocity gradient versus target velocity gradient
    ! ========================================================
    do i = 1, NR
      do j = 1, NZ
        dv_dz_residual(i,j) = &
             (dv_dz_grid(i,j) - dv_dz_suggested(i,j)) &
             / max(abs(dv_dz_suggested(i,j)), RES_EPS)
       end do
    end do
    iter_metric = sqrt(sum(dv_dz_residual**2) / real(NR*NZ, dp))
    write(*,'(A,ES12.4)') "Velocity-gradient residual metric = ", iter_metric

     v_target = 0.0d0
     v_target(:,NZ) = v_ff_grid(:,NZ)

     do j = NZ-1, 1, -1
        do i = 1, NR
           v_target(i,j) = v_target(i,j+1) - dv_dz_suggested(i,j+1) * dz
        end do
     end do

     do i = 1, NR
        do j = 1, NZ
           v_target(i,j) = max(v_target(i,j), V_MIN)
           v_target(i,j) = min(v_target(i,j), v_ff_grid(i,j))

    !========================================================
    ! Under-relaxation with limited velocity correction
    !========================================================
    dv_raw = VELOCITY_RELAX * (v_target(i,j) - v_grid(i,j))
    ! Limit change to 5 percent per iteration
    dv_lim = 0.10d0 * v_grid(i,j)
    if (dv_raw > dv_lim) then
      dv_raw = dv_lim
    end if
    if (dv_raw < -dv_lim) then
      dv_raw = -dv_lim
    end if
    v_new(i,j) = v_grid(i,j) + dv_raw
        end do
     end do

     call smooth_z(v_new)

     v_grid = v_new
     call write_velocity_profile(iteration)
     scatter_hist = scatter_hist + sum(scatter_thr, dim=2)
  end do
  call write_outputs()
  deallocate(u_thr, Fr_thr, Fz_thr, visit_thr, side_thr, scatter_thr)

contains
  !==================================================================
  !==================================================================
  subroutine write_velocity_profile(iter)
  integer, intent(in) :: iter
  integer :: jj, unit
  character(len=256) :: fname
    write(fname,'("./res/AC_gif/velocity_profile_iter_",I4.4,".dat")') iter
    open(newunit=unit, file=fname, status="replace", action="write")
    write(unit,'(A)') "# z  v_center  v_r_one_third  v_freefall  residual_r_one_third"
    do jj = 1, NZ
       write(unit,'(5ES24.15)') z_cent(jj), &
                              v_grid(1,jj), &
                              v_grid(i_vprof,jj), &
                              v_ff_profile(jj), &
                              dv_dz_residual(i_vprof,jj)
    end do
    close(unit)
  end subroutine write_velocity_profile

  real(dp) function gravity(z)
    real(dp), intent(in) :: z
    gravity = G * M_NS / (R_NS + z)**2
  end function gravity

  real(dp) function v_freefall(z)
    real(dp), intent(in) :: z
    v_freefall = sqrt(2.0d0 * G * M_NS / (R_NS + z))
  end function v_freefall

  real(dp) function rng_uniform(state)
    integer(int64), intent(inout) :: state
    integer(int64) :: x
    x = state
    x = ieor(x, ishft(x, 13))
    x = ieor(x, ishft(x, -7))
    x = ieor(x, ishft(x, 17))
    state = x
    rng_uniform = real(iand(x, huge(1_int64)), dp) / real(huge(1_int64), dp)
    if (rng_uniform <= 0.0d0) rng_uniform = 1.0d-16
    if (rng_uniform >= 1.0d0) rng_uniform = 1.0d0 - 1.0d-16
  end function rng_uniform

  subroutine isotropic_direction_rng(state, dir)
    integer(int64), intent(inout) :: state
    real(dp), intent(out) :: dir(3)
    real(dp) :: mu, phi, st
    mu = 2.0d0*rng_uniform(state) - 1.0d0
    phi = 2.0d0*pi*rng_uniform(state)
    st = sqrt(max(0.0d0, 1.0d0 - mu*mu))
    dir(1) = st*cos(phi)
    dir(2) = st*sin(phi)
    dir(3) = mu
  end subroutine isotropic_direction_rng

  subroutine upward_direction_rng(state, dir)
    integer(int64), intent(inout) :: state
    real(dp), intent(out) :: dir(3)
    real(dp) :: mu, phi, st
    mu = rng_uniform(state)
    phi = 2.0d0*pi*rng_uniform(state)
    st = sqrt(max(0.0d0, 1.0d0 - mu*mu))
    dir(1) = st*cos(phi)
    dir(2) = st*sin(phi)
    dir(3) = mu
  end subroutine upward_direction_rng

  subroutine build_cdf_2d(arr, n1, n2, cdf)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: arr(n1,n2)
    real(dp), intent(out) :: cdf(n1*n2)
    integer :: ii, jj, kk
    real(dp) :: total
    total = sum(arr)
    kk = 0
    do ii = 1, n1
       do jj = 1, n2
          kk = kk + 1
          if (kk == 1) then
             cdf(kk) = arr(ii,jj) / total
          else
             cdf(kk) = cdf(kk-1) + arr(ii,jj) / total
          end if
       end do
    end do
    cdf(n1*n2) = 1.0d0
  end subroutine build_cdf_2d

  subroutine build_cdf_1d(arr, n, cdf)
    integer, intent(in) :: n
    real(dp), intent(in) :: arr(n)
    real(dp), intent(out) :: cdf(n)
    integer :: ii
    real(dp) :: total
    total = sum(arr)
    do ii = 1, n
       if (ii == 1) then
          cdf(ii) = arr(ii) / total
       else
          cdf(ii) = cdf(ii-1) + arr(ii) / total
       end if
    end do
    cdf(n) = 1.0d0
  end subroutine build_cdf_1d

  integer function sample_cdf(cdf, n, u)
    integer, intent(in) :: n
    real(dp), intent(in) :: cdf(n), u
    integer :: lo, hi, mid
    lo = 1
    hi = n
    do while (lo < hi)
       mid = (lo + hi) / 2
       if (u <= cdf(mid)) then
          hi = mid
       else
          lo = mid + 1
       end if
    end do
    sample_cdf = lo
  end function sample_cdf

  subroutine sample_volume_source(state, pos)
    integer(int64), intent(inout) :: state
    real(dp), intent(out) :: pos(3)
    integer :: kk, ii, jj
    real(dp) :: r1, r2, rr, phi

    kk = sample_cdf(cell_cdf, NR*NZ, rng_uniform(state))
    ii = (kk-1)/NZ + 1
    jj = mod(kk-1,NZ) + 1

    r1 = r_edges(ii-1)
    r2 = r_edges(ii)

    rr = sqrt(r1*r1 + rng_uniform(state)*(r2*r2 - r1*r1))
    phi = 2.0d0*pi*rng_uniform(state)

    pos(1) = rr*cos(phi)
    pos(2) = rr*sin(phi)
    pos(3) = z_edges(jj-1) + rng_uniform(state)*dz
  end subroutine sample_volume_source

  subroutine sample_base_source(state, pos)
    integer(int64), intent(inout) :: state
    real(dp), intent(out) :: pos(3)
    integer :: ii
    real(dp) :: r1, r2, rr, phi

    ii = sample_cdf(base_ring_cdf, NR, rng_uniform(state))

    r1 = r_edges(ii-1)
    r2 = r_edges(ii)

    rr = sqrt(r1*r1 + rng_uniform(state)*(r2*r2 - r1*r1))
    phi = 2.0d0*pi*rng_uniform(state)

    pos(1) = rr*cos(phi)
    pos(2) = rr*sin(phi)
    pos(3) = 1.0d-5
  end subroutine sample_base_source

  subroutine get_cell_indices(pos, ii, jj, inside)
    real(dp), intent(in) :: pos(3)
    integer, intent(out) :: ii, jj
    logical, intent(out) :: inside
    real(dp) :: rr

    rr = sqrt(pos(1)**2 + pos(2)**2)

    if (rr < 0.0d0 .or. rr >= R_COL .or. pos(3) < 0.0d0 .or. pos(3) >= H_COL) then
       inside = .false.
       ii = -1
       jj = -1
       return
    end if

    ii = min(int(rr/dr) + 1, NR)
    jj = min(int(pos(3)/dz) + 1, NZ)
    inside = .true.
  end subroutine get_cell_indices

  real(dp) function distance_to_z_boundary(pos, dir)
    real(dp), intent(in) :: pos(3), dir(3)
    if (dir(3) > 0.0d0) then
       distance_to_z_boundary = (H_COL - pos(3)) / dir(3)
    else if (dir(3) < 0.0d0) then
       distance_to_z_boundary = -pos(3) / dir(3)
    else
       distance_to_z_boundary = huge(1.0d0)
    end if
  end function distance_to_z_boundary

  real(dp) function distance_to_outer_cylinder(pos, dir)
    real(dp), intent(in) :: pos(3), dir(3)
    real(dp) :: A, B, Cq, D, sqrtD, s1, s2, smin

    A = dir(1)**2 + dir(2)**2
    if (A < 1.0d-14) then
       distance_to_outer_cylinder = huge(1.0d0)
       return
    end if

    B = 2.0d0 * (pos(1)*dir(1) + pos(2)*dir(2))
    Cq = pos(1)**2 + pos(2)**2 - R_COL**2
    D = B**2 - 4.0d0*A*Cq

    if (D <= 0.0d0) then
       distance_to_outer_cylinder = huge(1.0d0)
       return
    end if

    sqrtD = sqrt(D)
    s1 = (-B + sqrtD) / (2.0d0*A)
    s2 = (-B - sqrtD) / (2.0d0*A)

    smin = huge(1.0d0)
    if (s1 > 1.0d-10) smin = min(smin, s1)
    if (s2 > 1.0d-10) smin = min(smin, s2)

    distance_to_outer_cylinder = smin
  end function distance_to_outer_cylinder

  subroutine gradient_z(arr, grad, dzloc)
    real(dp), intent(in) :: arr(NR,NZ), dzloc
    real(dp), intent(out) :: grad(NR,NZ)
    integer :: ii, jj

    do ii = 1, NR
       grad(ii,1) = (arr(ii,2) - arr(ii,1)) / dzloc
       do jj = 2, NZ-1
          grad(ii,jj) = (arr(ii,jj+1) - arr(ii,jj-1)) / (2.0d0*dzloc)
       end do
       grad(ii,NZ) = (arr(ii,NZ) - arr(ii,NZ-1)) / dzloc
    end do
  end subroutine gradient_z

  subroutine smooth_z(arr)
    real(dp), intent(inout) :: arr(NR,NZ)
    real(dp) :: tmp(NR,NZ)
    integer :: ii, jj
    tmp = arr
    do ii = 1, NR
       do jj = 2, NZ-1
          arr(ii,jj) = 0.25d0*tmp(ii,jj-1) + 0.5d0*tmp(ii,jj) + 0.25d0*tmp(ii,jj+1)
       end do
    end do
  end subroutine smooth_z

  subroutine mc_parallel_loop()
    integer :: tid, n, ii, jj, jwall
    integer :: scatter_count, step_count
    integer(int64) :: state
    real(dp) :: pos(3), dir(3), mid(3), er(3)
    real(dp) :: tau_to_scatter, alpha, s_scatter, s_z, s_side, s_escape, s_next
    real(dp) :: r_mid, r_now, z_now, mu_r, mu_z
    logical :: inside

    tid = omp_get_thread_num() + 1
    state = seeds(tid)

!$omp do schedule(dynamic)
    do n = 1, N_PACKETS

       if (rng_uniform(state) < p_base) then
          call sample_base_source(state, pos)
          call upward_direction_rng(state, dir)
       else
          call sample_volume_source(state, pos)
          call isotropic_direction_rng(state, dir)
       end if

       tau_to_scatter = -log(rng_uniform(state))

       scatter_count = 0
       step_count = 0

       do

          step_count = step_count + 1
          if (step_count > MAX_STEPS) exit

          call get_cell_indices(pos, ii, jj, inside)
          if (.not. inside) exit

          visit_thr(ii,jj,tid) = visit_thr(ii,jj,tid) + 1.0d0

          alpha = alpha_grid(ii,jj)

          if (alpha > 0.0d0) then
             s_scatter = tau_to_scatter / alpha
          else
             s_scatter = huge(1.0d0)
          end if

          s_z    = distance_to_z_boundary(pos, dir)
          s_side = distance_to_outer_cylinder(pos, dir)

          s_escape = min(s_z, s_side)
          s_next   = min(s_scatter, s_escape)

          if (.not. ieee_is_finite_local(s_next)) exit

          if (s_next < 1.0d-12) then
             pos = pos + 1.0d-5 * dir
             cycle
          end if

          mid = pos + 0.5d0*s_next*dir
          r_mid = sqrt(mid(1)**2 + mid(2)**2)

          if (r_mid > 0.5d0*dr) then
             er(1) = mid(1) / r_mid
             er(2) = mid(2) / r_mid
             er(3) = 0.0d0
             mu_r = dir(1)*er(1) + dir(2)*er(2)
          else
             mu_r = 0.0d0
          end if

          mu_z = dir(3)

          u_thr(ii,jj,tid)  = u_thr(ii,jj,tid)  + packet_luminosity*s_next/c_light
          Fr_thr(ii,jj,tid) = Fr_thr(ii,jj,tid) + packet_luminosity*mu_r*s_next
          Fz_thr(ii,jj,tid) = Fz_thr(ii,jj,tid) + packet_luminosity*mu_z*s_next

          pos = pos + s_next*dir

          if (s_escape <= s_scatter) then

             r_now = sqrt(pos(1)**2 + pos(2)**2)
             z_now = pos(3)

             if (r_now >= R_COL*(1.0d0 - 1.0d-8)) then

                if (z_now >= 0.0d0 .and. z_now < H_COL) then
                   jwall = min(int(z_now/dz) + 1, NZ)
                   side_thr(jwall,tid) = side_thr(jwall,tid) + packet_luminosity
                end if

                exit

             else if (z_now >= H_COL*(1.0d0 - 1.0d-8)) then

                exit

             else if (z_now <= 1.0d-10) then

                dir(3) = abs(dir(3))
                pos(3) = 1.0d-5
                cycle

             end if

          end if

          scatter_count = scatter_count + 1

          if (scatter_count > MAX_SCATTERS) exit

          call isotropic_direction_rng(state, dir)
          tau_to_scatter = -log(rng_uniform(state))

       end do

       if (scatter_count > MAX_SCATTERS) then
          scatter_thr(MAX_SCATTERS+1,tid) = scatter_thr(MAX_SCATTERS+1,tid) + 1
       else
          scatter_thr(scatter_count,tid) = scatter_thr(scatter_count,tid) + 1
       end if

    end do
!$omp end do

    seeds(tid) = state

  end subroutine mc_parallel_loop

  logical function ieee_is_finite_local(x)
    real(dp), intent(in) :: x
    ieee_is_finite_local = (abs(x) < huge(1.0d0))
  end function ieee_is_finite_local

  subroutine write_outputs()
    integer :: ii, jj, unit
    real(dp) :: F_side_esc_z(NZ)
    real(dp) :: visit_final(NR,NZ)

    F_side_esc_z = sum(side_thr, dim=2) / side_area_zbin
    visit_final = sum(visit_thr, dim=3)

    open(newunit=unit, file="./res/text_AC", status="replace", action="write")
    do jj = 1, NZ
      write(unit,'(6ES24.15)') z_cent(jj), v_grid(1,jj), v_grid(i_vprof,jj), &
                            dv_dz_residual(i_vprof,jj), &
                            dv_dz_grid(i_vprof,jj) , dv_dz_suggested(i_vprof,jj)
    end do
    close(unit)

    open(newunit=unit, file="./res/test_AC2d", status="replace", action="write")
    do ii = 1, NR
       do jj = 1, NZ
          write(unit,'(3ES24.15)') r_cent(ii), z_cent(jj), v_grid(ii,jj)
       end do
       write(unit,*)
    end do
    close(unit)

    open(newunit=unit, file="./res/test_AC_maps", status="replace", action="write")
    do ii = 1, NR
       do jj = 1, NZ
          write(unit,'(9ES24.15)') r_cent(ii), z_cent(jj), q_grid(ii,jj), &
               u_rad(ii,jj), Fr(ii,jj), Fz(ii,jj), T_keV(ii,jj), &
               P_gas(ii,jj), P_rad(ii,jj)
       end do
       write(unit,*)
    end do
    close(unit)

    open(newunit=unit, file="./res/test_AC_side_flux", status="replace", action="write")
    do jj = 1, NZ
       write(unit,'(2ES24.15)') z_cent(jj), F_side_esc_z(jj)
    end do
    close(unit)

    open(newunit=unit, file="./res/test_AC_scatter_hist", status="replace", action="write")
    do ii = 0, MAX_SCATTERS+1
       write(unit,'(I10,1X,I16)') ii, scatter_hist(ii)
    end do
    close(unit)

    open(newunit=unit, file="./res/test_AC_visit", status="replace", action="write")
    do ii = 1, NR
       do jj = 1, NZ
          write(unit,'(3ES24.15)') r_cent(ii), z_cent(jj), visit_final(ii,jj)
       end do
       write(unit,*)
    end do
    close(unit)

    write(*,*)
    write(*,*) "Output written to:"
    write(*,*) "./res/text_AC"
    write(*,*) "./res/test_AC2d"
    write(*,*) "./res/test_AC_maps"
    write(*,*) "./res/test_AC_side_flux"
    write(*,*) "./res/test_AC_scatter_hist"
    write(*,*) "./res/test_AC_visit"

  end subroutine write_outputs

end subroutine test_AC
