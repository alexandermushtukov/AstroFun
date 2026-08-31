!=======================================================================================
! Monte-Carlo X-ray disk irradiation kernel ported from the provided Python script.
! Public API:
!   call xray_disk_mc_profile(..., profile, ierr)
! Output profile(i,1:7):
!   1  R_cm
!   2  F_viscous_erg_s_cm2
!   3  F_direct_irradiation_erg_s_cm2
!   4  F_scattered_irradiation_erg_s_cm2
!   5  T_viscous_K
!   6  T_direct_irradiation_K
!   7  T_scattered_irradiation_K
! No plotting, no file I/O, no terminal output.
!=======================================================================================
module xray_disk_mc_module
use, intrinsic :: iso_fortran_env, only: real64
use omp_lib
implicit none
private
public :: xray_disk_mc_profile, test_xray_disk_mc_profile
integer, parameter :: dp = real64
real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
real(dp), parameter :: mec2_kev = 511.0_dp
real(dp), parameter :: sigma_sb_cgs = 5.670374419e-5_dp
real(dp), parameter :: g_cgs = 6.67430e-8_dp
real(dp), parameter :: m_sun_g = 1.98847e33_dp
real(dp), parameter :: m_p_cgs = 1.67262192369e-24_dp
real(dp), parameter :: sigma_t_cgs = 6.6524587321e-25_dp
real(dp), parameter :: kappa_t_h = sigma_t_cgs / m_p_cgs
integer, parameter :: max_steps_default = 1000000
integer, parameter :: max_internal_scatters_default = 100000
integer, parameter :: planck_nmax = 800
integer(kind=8), parameter :: rng_modulus = 2147483647_8
integer(kind=8), parameter :: rng_multiplier = 16807_8

real(dp), save :: planck_cdf(planck_nmax)
logical, save :: planck_ready = .false.

type :: photon_result
   logical :: disk_absorbed = .false.
   integer :: disk_component = 0
   real(dp) :: disk_rho_model = 0.0_dp
   real(dp) :: disk_phi = 0.0_dp
   real(dp) :: disk_mu = 0.0_dp
   real(dp) :: launch_weight = 0.0_dp
   real(dp) :: packet_weight = 0.0_dp
   real(dp) :: e0_kev = 0.0_dp
   real(dp) :: efinal_kev = 0.0_dp
end type photon_result

type :: curtain_geometry
   real(dp) :: m_ns_g = 0.0_dp
   real(dp) :: r_ns_cm = 0.0_dp
   real(dp) :: rm_model = 0.0_dp
   real(dp) :: rm_cm = 0.0_dp
   real(dp) :: dotm_g_s = 0.0_dp
   real(dp) :: omega = 0.0_dp
   real(dp) :: alpha = 0.0_dp
   real(dp) :: phase = 0.0_dp
   real(dp) :: delta_lambda_disc = 0.0_dp
   real(dp) :: v_initial_cms = 0.0_dp
   real(dp) :: velocity_floor_cms = 0.0_dp
   real(dp) :: accretion_coverage_fraction = 0.0_dp
   integer :: nphi = 0
   integer :: nstream = 0
   integer :: n_valid = 0
   real(dp) :: magnetic_to_global(3,3) = 0.0_dp
   real(dp) :: global_to_magnetic(3,3) = 0.0_dp
   logical, allocatable :: coverage_exists(:,:)
   logical, allocatable :: coverage_reaches_surface(:,:)
   real(dp), allocatable :: coverage_phi(:)
   real(dp), allocatable :: coverage_lambda_start(:,:)
   real(dp), allocatable :: coverage_lambda_stop(:,:)
   real(dp), allocatable :: coverage_lambda_surface(:,:)
   integer, allocatable :: valid_idx(:)
   integer, allocatable :: valid_hem(:)
end type curtain_geometry

contains

subroutine xray_disk_mc_profile(N, L39, B12, T_keV, nbins, profile, ierr, &
   R_ns, R_ns_km, r_stop_factor, step_coarse, tol, advance_factor, &
   L_x_erg_s, alpha_visc, m_msun, disk_rho_max_factor, disk_rho_min_factor, &
   forced_thin_scattering, min_packet_weight, seed, use_two_sided_disk_area, &
   dipole_tilt_rad, phase_cycles, H_over_Rm, spin_period_s, v_initial_cms, &
   velocity_floor_cms, coverage_azimuth_bins, coverage_points_along_stream, tau_crit)
integer, intent(in) :: N, nbins
real(dp), intent(in) :: L39, B12, T_keV
real(dp), intent(out) :: profile(nbins,9)
integer, intent(out) :: ierr
real(dp), intent(in), optional :: R_ns, R_ns_km, r_stop_factor, step_coarse, tol, advance_factor
real(dp), intent(in), optional :: L_x_erg_s, alpha_visc, m_msun
real(dp), intent(in), optional :: disk_rho_max_factor, disk_rho_min_factor
real(dp), intent(in), optional :: min_packet_weight
logical, intent(in), optional :: forced_thin_scattering, use_two_sided_disk_area
integer, intent(in), optional :: seed
real(dp), intent(in), optional :: dipole_tilt_rad, phase_cycles, H_over_Rm, spin_period_s
real(dp), intent(in), optional :: v_initial_cms, velocity_floor_cms, tau_crit
integer, intent(in), optional :: coverage_azimuth_bins, coverage_points_along_stream
real(dp) :: Rns, Rnskm, rstopfac, stepc, tolerance, advfac
real(dp) :: Lx, alpha_disk, mstar, rho_max_fac, rho_min_fac, minw, taucrit
real(dp) :: tilt, phasecyc, hoverrm, pspin, vinit, vfloor
integer :: cov_nphi, cov_nstream
logical :: forced, two_sided
real(dp) :: Rm, rho_min_model, rho_max_model
real(dp) :: rho_edges(nbins+1), rho_centers(nbins)
real(dp) :: R_edges_cm(nbins+1), R_centers_cm(nbins), dR_cm(nbins)
real(dp) :: direct_bins(nbins), scattered_bins(nbins)
real(dp) :: direct_mu_bins(nbins), scattered_mu_bins(nbins)
real(dp) :: mean_mu_direct(nbins), mean_mu_scattered(nbins)
real(dp) :: total_wE0, total_wEdirect, total_wEscattered
real(dp) :: direct_dL(nbins), scattered_dL(nbins)
real(dp) :: F_direct(nbins), F_scattered(nbins), T_direct(nbins), T_scattered(nbins)
real(dp) :: F_visc(nbins), T_visc(nbins)
real(dp) :: area(nbins), Rns_cm, Rinner_cm
real(dp) :: rawsum, Ldir, Lscat
integer :: i, ibin, base_seed
integer(kind=8) :: state
integer :: nthreads, tid
real(dp), allocatable :: direct_thr(:,:), scattered_thr(:,:), direct_mu_thr(:,:), scattered_mu_thr(:,:)
real(dp), allocatable :: total0_thr(:), totaldir_thr(:), totalscat_thr(:)
type(photon_result) :: res
type(curtain_geometry) :: geom

profile = 0.0_dp
ierr = 0
if (N <= 0 .or. nbins <= 0 .or. L39 <= 0.0_dp .or. B12 <= 0.0_dp .or. T_keV <= 0.0_dp) then
   ierr = 1
   return
end if

Rns = opt_real(R_ns, 1.0_dp)
Rnskm = opt_real(R_ns_km, 10.0_dp)
rstopfac = opt_real(r_stop_factor, 100.0_dp)
stepc = opt_real(step_coarse, 0.5_dp)
tolerance = opt_real(tol, 1.0e-5_dp)
advfac = opt_real(advance_factor, 50.0_dp)
Lx = opt_real(L_x_erg_s, L39*1.0e39_dp)
alpha_disk = opt_real(alpha_visc, 0.1_dp)
mstar = opt_real(m_msun, 1.4_dp)
rho_max_fac = opt_real(disk_rho_max_factor, 100.0_dp)
rho_min_fac = opt_real(disk_rho_min_factor, 1.0_dp)
minw = opt_real(min_packet_weight, 1.0e-12_dp)
forced = opt_logical(forced_thin_scattering, .false.)
two_sided = opt_logical(use_two_sided_disk_area, .true.)
tilt = opt_real(dipole_tilt_rad, 0.0_dp)
phasecyc = opt_real(phase_cycles, 0.0_dp)
hoverrm = opt_real(H_over_Rm, 0.06_dp)
pspin = opt_real(spin_period_s, 5000.0_dp)
vinit = opt_real(v_initial_cms, 2.0e7_dp)
vfloor = opt_real(velocity_floor_cms, 1.0e5_dp)
taucrit = opt_real(tau_crit, 0.5_dp)
cov_nphi = opt_int(coverage_azimuth_bins, 720)
cov_nstream = opt_int(coverage_points_along_stream, 500)
base_seed = 1357911
if (present(seed)) base_seed = seed

if (Rns <= 0.0_dp .or. Rnskm <= 0.0_dp .or. rho_max_fac <= 0.0_dp .or. &
    pspin <= 0.0_dp .or. hoverrm < 0.0_dp .or. vinit < 0.0_dp .or. vfloor <= 0.0_dp .or. &
    cov_nphi < 8 .or. cov_nstream < 2 .or. taucrit <= 0.0_dp) then
   ierr = 2
   return
end if

call init_planck_cdf()
Rm = Rm_from_LB(L39,B12)
call init_curtain_geometry(geom, Rm, Lx, Rnskm, mstar, tilt, phasecyc, hoverrm, pspin, &
                           vinit, vfloor, cov_nphi, cov_nstream, ierr)
if (ierr /= 0) then
   call destroy_curtain_geometry(geom)
   return
end if

rho_min_model = max(rho_min_fac*Rm, Rm, 1.0e-12_dp)
rho_max_model = rho_max_fac*Rm
if (rho_max_model <= rho_min_model) rho_max_model = rho_min_model*(1.0_dp + 1.0e-6_dp)

do i = 1, nbins+1
   rho_edges(i) = 10.0_dp**(log10(rho_min_model) + real(i-1,dp)* &
        (log10(rho_max_model)-log10(rho_min_model))/real(nbins,dp))
end do
do i = 1, nbins
   rho_centers(i) = sqrt(rho_edges(i)*rho_edges(i+1))
end do

nthreads = omp_get_max_threads()
allocate(direct_thr(nbins,nthreads), scattered_thr(nbins,nthreads))
allocate(direct_mu_thr(nbins,nthreads), scattered_mu_thr(nbins,nthreads))
allocate(total0_thr(nthreads), totaldir_thr(nthreads), totalscat_thr(nthreads))
direct_thr = 0.0_dp
scattered_thr = 0.0_dp
direct_mu_thr = 0.0_dp
scattered_mu_thr = 0.0_dp
total0_thr = 0.0_dp
totaldir_thr = 0.0_dp
totalscat_thr = 0.0_dp

!$omp parallel default(shared) private(i,tid,res,ibin,state)
tid = omp_get_thread_num()+1
!$omp do schedule(dynamic,10)
do i = 1, N
   state = photon_seed(base_seed,i)
   call propagate_one_photon(state, geom, L39, B12, T_keV, Rns, Rnskm, rstopfac, &
        stepc, tolerance, advfac, Lx, alpha_disk, mstar, rho_max_fac, forced, minw, taucrit, res)

   total0_thr(tid) = total0_thr(tid) + res%launch_weight*res%e0_kev
   if (res%disk_absorbed) then
      ibin = locate_bin(rho_edges,nbins,res%disk_rho_model)
      if (res%disk_component == 1) then
         totaldir_thr(tid) = totaldir_thr(tid) + res%packet_weight*res%efinal_kev
         if (ibin > 0) then
            direct_thr(ibin,tid) = direct_thr(ibin,tid) + res%packet_weight*res%efinal_kev
            direct_mu_thr(ibin,tid) = direct_mu_thr(ibin,tid) + &
                 res%packet_weight*res%efinal_kev*res%disk_mu
         end if
      else if (res%disk_component == 2) then
         totalscat_thr(tid) = totalscat_thr(tid) + res%packet_weight*res%efinal_kev
         if (ibin > 0) then
            scattered_thr(ibin,tid) = scattered_thr(ibin,tid) + res%packet_weight*res%efinal_kev
            scattered_mu_thr(ibin,tid) = scattered_mu_thr(ibin,tid) + &
                 res%packet_weight*res%efinal_kev*res%disk_mu
         end if
      end if
   end if
end do
!$omp end do
!$omp end parallel

direct_bins = sum(direct_thr,dim=2)
scattered_bins = sum(scattered_thr,dim=2)
direct_mu_bins = sum(direct_mu_thr,dim=2)
scattered_mu_bins = sum(scattered_mu_thr,dim=2)
total_wE0 = sum(total0_thr)
total_wEdirect = sum(totaldir_thr)
total_wEscattered = sum(totalscat_thr)

deallocate(direct_thr,scattered_thr,direct_mu_thr,scattered_mu_thr)
deallocate(total0_thr,totaldir_thr,totalscat_thr)

Rns_cm = Rnskm*1.0e5_dp
R_edges_cm = rho_edges*Rns_cm
R_centers_cm = rho_centers*Rns_cm
dR_cm = R_edges_cm(2:nbins+1)-R_edges_cm(1:nbins)
if (two_sided) then
   area = 4.0_dp*pi*R_centers_cm*dR_cm
else
   area = 2.0_dp*pi*R_centers_cm*dR_cm
end if

direct_dL = 0.0_dp
scattered_dL = 0.0_dp
if (total_wE0 > 0.0_dp) then
   rawsum = sum(direct_bins)
   Ldir = Lx*total_wEdirect/total_wE0
   if (rawsum > 0.0_dp) direct_dL = Ldir*direct_bins/rawsum
   rawsum = sum(scattered_bins)
   Lscat = Lx*total_wEscattered/total_wE0
   if (rawsum > 0.0_dp) scattered_dL = Lscat*scattered_bins/rawsum
end if

mean_mu_direct = 0.0_dp
mean_mu_scattered = 0.0_dp
do i = 1, nbins
   if (direct_bins(i) > 0.0_dp) mean_mu_direct(i) = direct_mu_bins(i)/direct_bins(i)
   if (scattered_bins(i) > 0.0_dp) mean_mu_scattered(i) = scattered_mu_bins(i)/scattered_bins(i)
end do

F_direct = 0.0_dp
F_scattered = 0.0_dp
T_direct = 0.0_dp
T_scattered = 0.0_dp
do i = 1, nbins
   if (area(i) > 0.0_dp) then
      F_direct(i) = direct_dL(i)/area(i)
      F_scattered(i) = scattered_dL(i)/area(i)
   end if
   if (F_direct(i) > 0.0_dp) T_direct(i) = (F_direct(i)/sigma_sb_cgs)**0.25_dp
   if (F_scattered(i) > 0.0_dp) T_scattered(i) = (F_scattered(i)/sigma_sb_cgs)**0.25_dp
end do

Rinner_cm = Rm*Rns_cm
call viscous_disk_profile(R_centers_cm,nbins,Lx,Rnskm,mstar,Rinner_cm,F_visc,T_visc)

do i = 1, nbins
   profile(i,1) = R_centers_cm(i)
   profile(i,2) = F_visc(i)
   profile(i,3) = F_direct(i)
   profile(i,4) = F_scattered(i)
   profile(i,5) = T_visc(i)
   profile(i,6) = T_direct(i)
   profile(i,7) = T_scattered(i)
   profile(i,8) = mean_mu_direct(i)
   profile(i,9) = mean_mu_scattered(i)
end do

call destroy_curtain_geometry(geom)
end subroutine xray_disk_mc_profile

subroutine init_curtain_geometry(g, Rm, Lx, Rnskm, mstar, tilt, phase_cycles, hoverrm, &
                                 pspin, vinit, vfloor, nphi, nstream, ierr)
type(curtain_geometry), intent(inout) :: g
real(dp), intent(in) :: Rm,Lx,Rnskm,mstar,tilt,phase_cycles,hoverrm,pspin,vinit,vfloor
integer, intent(in) :: nphi,nstream
integer, intent(out) :: ierr
integer :: k,row,hem
real(dp) :: phi,ls,lstop,lsurf
logical :: ex,re

ierr = 0
g%m_ns_g = mstar*m_sun_g
g%r_ns_cm = Rnskm*1.0e5_dp
g%rm_model = Rm
g%rm_cm = Rm*g%r_ns_cm
g%dotm_g_s = Lx*g%r_ns_cm/(g_cgs*g%m_ns_g)
g%omega = 2.0_dp*pi/pspin
g%alpha = tilt
g%phase = 2.0_dp*pi*phase_cycles
g%delta_lambda_disc = hoverrm
g%v_initial_cms = vinit
g%velocity_floor_cms = vfloor
g%nphi = nphi
g%nstream = nstream
call build_rotations(g)
allocate(g%coverage_exists(2,nphi),g%coverage_reaches_surface(2,nphi))
allocate(g%coverage_phi(nphi))
allocate(g%coverage_lambda_start(2,nphi),g%coverage_lambda_stop(2,nphi))
allocate(g%coverage_lambda_surface(2,nphi))
allocate(g%valid_idx(2*nphi),g%valid_hem(2*nphi))
g%n_valid = 0

do k = 1, nphi
   phi = 2.0_dp*pi*real(k-1,dp)/real(nphi,dp)
   g%coverage_phi(k) = phi
   do row = 1, 2
      if (row == 1) then
         hem = -1
      else
         hem = 1
      end if
      call branch_extent(g,phi,hem,ex,re,ls,lstop,lsurf)
      g%coverage_exists(row,k) = ex
      g%coverage_reaches_surface(row,k) = re
      g%coverage_lambda_start(row,k) = ls
      g%coverage_lambda_stop(row,k) = lstop
      g%coverage_lambda_surface(row,k) = lsurf
      if (re) then
         g%n_valid = g%n_valid+1
         g%valid_idx(g%n_valid) = k
         g%valid_hem(g%n_valid) = hem
      end if
   end do
end do
if (g%n_valid <= 0) then
   ierr = 3
   return
end if
g%accretion_coverage_fraction = real(g%n_valid,dp)/real(2*nphi,dp)
end subroutine init_curtain_geometry

subroutine destroy_curtain_geometry(g)
type(curtain_geometry), intent(inout) :: g
if (allocated(g%coverage_exists)) deallocate(g%coverage_exists)
if (allocated(g%coverage_reaches_surface)) deallocate(g%coverage_reaches_surface)
if (allocated(g%coverage_phi)) deallocate(g%coverage_phi)
if (allocated(g%coverage_lambda_start)) deallocate(g%coverage_lambda_start)
if (allocated(g%coverage_lambda_stop)) deallocate(g%coverage_lambda_stop)
if (allocated(g%coverage_lambda_surface)) deallocate(g%coverage_lambda_surface)
if (allocated(g%valid_idx)) deallocate(g%valid_idx)
if (allocated(g%valid_hem)) deallocate(g%valid_hem)
end subroutine destroy_curtain_geometry

subroutine build_rotations(g)
type(curtain_geometry), intent(inout) :: g
real(dp) :: Rx(3,3),Rz(3,3),ca,sa,cp,sp
ca = cos(g%alpha); sa = sin(g%alpha)
cp = cos(g%phase); sp = sin(g%phase)
Rx = 0.0_dp
Rx(1,1)=1.0_dp; Rx(2,2)=ca; Rx(2,3)=sa; Rx(3,2)=-sa; Rx(3,3)=ca
Rz = 0.0_dp
Rz(1,1)=cp; Rz(1,2)=-sp; Rz(2,1)=sp; Rz(2,2)=cp; Rz(3,3)=1.0_dp
g%magnetic_to_global = matmul(Rz,Rx)
g%global_to_magnetic = transpose(g%magnetic_to_global)
end subroutine build_rotations

subroutine global_to_mag(g,vg,vm)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: vg(3)
real(dp), intent(out) :: vm(3)
vm = matmul(g%global_to_magnetic,vg)
end subroutine global_to_mag

subroutine mag_to_global(g,vm,vg)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: vm(3)
real(dp), intent(out) :: vg(3)
vg = matmul(g%magnetic_to_global,vm)
end subroutine mag_to_global

subroutine magnetic_spherical(p,r,phi,lamb)
real(dp), intent(in) :: p(3)
real(dp), intent(out) :: r,phi,lamb
r = norm3(p)
if (r <= 0.0_dp) then
   phi=0.0_dp; lamb=0.0_dp; return
end if
phi = atan2(p(2),p(1))
if (phi < 0.0_dp) phi = phi+2.0_dp*pi
lamb = asin(max(-1.0_dp,min(1.0_dp,p(3)/r)))
end subroutine magnetic_spherical

function disc_latitude(g,phi) result(lamb)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: phi
real(dp) :: lamb,ee_m(3),ez_m(3),ee_g(3),ez_g(3)
ee_m = [cos(phi),sin(phi),0.0_dp]
ez_m = [0.0_dp,0.0_dp,1.0_dp]
call mag_to_global(g,ee_m,ee_g)
call mag_to_global(g,ez_m,ez_g)
lamb = atan2(-ee_g(3),ez_g(3))
if (lamb > 0.5_dp*pi) lamb = lamb-pi
if (lamb < -0.5_dp*pi) lamb = lamb+pi
end function disc_latitude

function field_R0(g,phi) result(R0)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: phi
real(dp) :: R0,c
c = cos(disc_latitude(g,phi))
if (abs(c) < 1.0e-12_dp) then
   R0 = huge(1.0_dp)
else
   R0 = g%rm_model/(c*c)
end if
end function field_R0

function surface_latitude_abs(g,phi) result(lamb)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: phi
real(dp) :: lamb,R0,q
R0 = field_R0(g,phi)
q = 1.0_dp/R0
if (q <= 0.0_dp .or. q > 1.0_dp) then
   lamb = huge(1.0_dp)
else
   lamb = acos(sqrt(q))
end if
end function surface_latitude_abs

function branch_start(g,phi,hem) result(lamb)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: phi
integer, intent(in) :: hem
real(dp) :: lamb
if (hem > 0) then
   lamb = disc_latitude(g,phi)+g%delta_lambda_disc
else
   lamb = disc_latitude(g,phi)-g%delta_lambda_disc
end if
end function branch_start

function branch_surface(g,phi,hem) result(lamb)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: phi
integer, intent(in) :: hem
real(dp) :: lamb,q
q = surface_latitude_abs(g,phi)
if (hem > 0) then
   lamb = q
else
   lamb = -q
end if
end function branch_surface

function branch_geometric(g,phi,hem) result(ok)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: phi
integer, intent(in) :: hem
logical :: ok
real(dp) :: a,b
a = branch_start(g,phi,hem)
b = branch_surface(g,phi,hem)
if (abs(b) >= 0.5_dp*huge(1.0_dp)) then
   ok = .false.
else if (hem > 0) then
   ok = (a < b)
else
   ok = (a > b)
end if
end function branch_geometric

subroutine magnetic_position(g,phi,lamb,p)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: phi,lamb
real(dp), intent(out) :: p(3)
real(dp) :: pm(3),R0,r,c
R0 = field_R0(g,phi)
c = cos(lamb)
r = R0*c*c
pm = [r*c*cos(phi),r*c*sin(phi),r*sin(lamb)]
call mag_to_global(g,pm,p)
end subroutine magnetic_position

function effective_potential(g,p) result(phi_eff)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: p(3)
real(dp) :: phi_eff,r_cm,varpi_cm
r_cm = norm3(p)*g%r_ns_cm
if (r_cm <= 0.0_dp) then
   phi_eff = -huge(1.0_dp); return
end if
varpi_cm = sqrt(p(1)*p(1)+p(2)*p(2))*g%r_ns_cm
phi_eff = -g_cgs*g%m_ns_g/r_cm - 0.5_dp*g%omega*g%omega*varpi_cm*varpi_cm
end function effective_potential

function ballistic_v2(g,phi,lamb,hem) result(v2)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: phi,lamb
integer, intent(in) :: hem
real(dp) :: v2,p0(3),p(3),ls
if (.not. branch_geometric(g,phi,hem)) then
   v2 = -1.0_dp; return
end if
ls = branch_start(g,phi,hem)
call magnetic_position(g,phi,ls,p0)
call magnetic_position(g,phi,lamb,p)
v2 = g%v_initial_cms*g%v_initial_cms + 2.0_dp*(effective_potential(g,p0)-effective_potential(g,p))
end function ballistic_v2

function refine_barrier(g,phi,hem,a0,b0) result(root)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: phi,a0,b0
integer, intent(in) :: hem
real(dp) :: root,a,b,m,fm
integer :: k
a=a0; b=b0
do k=1,60
   m=0.5_dp*(a+b)
   fm=ballistic_v2(g,phi,m,hem)
   if (fm > 0.0_dp) then
      a=m
   else
      b=m
   end if
end do
root=0.5_dp*(a+b)
end function refine_barrier

subroutine branch_extent(g,phi,hem,exists,reaches,ls,lstop,lsurf)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: phi
integer, intent(in) :: hem
logical, intent(out) :: exists,reaches
real(dp), intent(out) :: ls,lstop,lsurf
real(dp) :: prev,cur,v2
integer :: k
exists=.false.; reaches=.false.; ls=0.0_dp; lstop=0.0_dp; lsurf=0.0_dp
if (.not. branch_geometric(g,phi,hem)) return
ls=branch_start(g,phi,hem)
lsurf=branch_surface(g,phi,hem)
prev=ls
v2=ballistic_v2(g,phi,prev,hem)
if (v2 <= 0.0_dp) return
exists=.true.
do k=2,g%nstream
   cur=ls+(lsurf-ls)*real(k-1,dp)/real(g%nstream-1,dp)
   v2=ballistic_v2(g,phi,cur,hem)
   if (v2 <= 0.0_dp) then
      lstop=refine_barrier(g,phi,hem,prev,cur)
      return
   end if
   prev=cur
end do
lstop=lsurf
reaches=.true.
end subroutine branch_extent

function nearest_coverage(g,phi) result(k)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: phi
integer :: k
real(dp) :: q,dphi
q = modulo(phi,2.0_dp*pi)
dphi = 2.0_dp*pi/real(g%nphi,dp)
k = int(floor(q/dphi+0.5_dp))+1
if (k > g%nphi) k=1
end function nearest_coverage

function classify_branch(g,phi,lamb) result(hem)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: phi,lamb
integer :: hem,k
real(dp) :: a,b
hem=0
k=nearest_coverage(g,phi)
if (g%coverage_exists(2,k)) then
   a=g%coverage_lambda_start(2,k); b=g%coverage_lambda_stop(2,k)
   if (lamb >= a .and. lamb <= b) then
      hem=1; return
   end if
end if
if (g%coverage_exists(1,k)) then
   a=g%coverage_lambda_start(1,k); b=g%coverage_lambda_stop(1,k)
   if (lamb >= b .and. lamb <= a) hem=-1
end if
end function classify_branch

function surface_density(g,phi,lamb,hem) result(sigma)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: phi,lamb
integer, intent(in) :: hem
real(dp) :: sigma,v2,v,R0cm,c
v2=ballistic_v2(g,phi,lamb,hem)
if (v2 <= 0.0_dp) then
   sigma=-1.0_dp; return
end if
v=max(sqrt(v2),g%velocity_floor_cms)
R0cm=field_R0(g,phi)*g%r_ns_cm
c=max(abs(cos(lamb)),1.0e-12_dp)
sigma=(g%dotm_g_s/(4.0_dp*pi*g%accretion_coverage_fraction))/(v*R0cm*c**3)
end function surface_density

function tau_thomson_normal(g,p) result(tau)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: p(3)
real(dp) :: tau,pm(3),r,phi,lamb,s
integer :: hem
call global_to_mag(g,p,pm)
call magnetic_spherical(pm,r,phi,lamb)
hem=classify_branch(g,phi,lamb)
if (hem == 0) then
   tau=-1.0_dp; return
end if
s=surface_density(g,phi,lamb,hem)
if (s <= 0.0_dp) then
   tau=-1.0_dp
else
   tau=kappa_t_h*s
end if
end function tau_thomson_normal

function curtain_residual(g,p) result(f)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: p(3)
real(dp) :: f,pm(3),r,phi,lamb,R0
call global_to_mag(g,p,pm)
call magnetic_spherical(pm,r,phi,lamb)
if (r <= 0.0_dp) then
   f=huge(1.0_dp); return
end if
R0=field_R0(g,phi)
f=r-R0*cos(lamb)**2
end function curtain_residual

function loaded_curtain_point(g,p) result(ok)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: p(3)
logical :: ok
real(dp) :: pm(3),r,phi,lamb
call global_to_mag(g,p,pm)
call magnetic_spherical(pm,r,phi,lamb)
ok=(classify_branch(g,phi,lamb) /= 0)
end function loaded_curtain_point

subroutine curtain_normal(g,p,n)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: p(3)
real(dp), intent(out) :: n(3)
real(dp), parameter :: h=1.0e-4_dp
real(dp) :: pp(3),pm(3)
integer :: k
do k=1,3
   pp=p; pm=p
   pp(k)=pp(k)+h; pm(k)=pm(k)-h
   n(k)=(curtain_residual(g,pp)-curtain_residual(g,pm))/(2.0_dp*h)
end do
if (norm3(n) < 1.0e-14_dp) n=p
n=normalize(n)
end subroutine curtain_normal

function residual_ray(g,t,p,n) result(f)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: t,p(3),n(3)
real(dp) :: f,q(3)
q=p+t*n
f=curtain_residual(g,q)
end function residual_ray

function bisect_curtain_root(g,p,n,a0,b0,fa0,fb0,tol) result(root)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: p(3),n(3),a0,b0,fa0,fb0,tol
real(dp) :: root,a,b,fa,m,fm
integer :: k
a=a0; b=b0; fa=fa0
do k=1,200
   m=0.5_dp*(a+b)
   fm=residual_ray(g,m,p,n)
   if ((b-a) <= tol .or. fm == 0.0_dp) then
      root=m; return
   end if
   if (fa*fm <= 0.0_dp) then
      b=m
   else
      a=m; fa=fm
   end if
end do
root=0.5_dp*(a+b)
end function bisect_curtain_root

function ray_curtain_first_intersection(g,p0,n_dir,t_min,t_max,step_min,tol,t_hit) result(ok)
type(curtain_geometry), intent(in) :: g
real(dp), intent(in) :: p0(3),n_dir(3),t_min,t_max,step_min,tol
real(dp), intent(out) :: t_hit
logical :: ok
real(dp) :: n(3),tprev,tcur,fprev,fcur,step,root,rhere,q(3),hit(3),root_skip
n=normalize(n_dir)
ok=.false.; t_hit=0.0_dp
root_skip=max(5.0e-5_dp,5.0_dp*tol)
tprev=max(t_min,0.0_dp)
fprev=residual_ray(g,tprev,p0,n)
if (abs(fprev) < 1.0e-12_dp) then
   tprev=tprev+max(root_skip,tol)
   if (tprev >= t_max) return
   fprev=residual_ray(g,tprev,p0,n)
end if
do while (tprev < t_max)
   q=p0+tprev*n
   rhere=norm3(q)
   step=max(step_min,0.02_dp*rhere)
   step=min(step,10.0_dp)
   tcur=min(tprev+step,t_max)
   fcur=residual_ray(g,tcur,p0,n)
   root=-1.0_dp
   if (fprev == 0.0_dp) then
      root=tprev
   else if (fcur == 0.0_dp) then
      root=tcur
   else if (fprev*fcur < 0.0_dp) then
      root=bisect_curtain_root(g,p0,n,tprev,tcur,fprev,fcur,tol)
   end if
   if (root > t_min) then
      hit=p0+root*n
      if (loaded_curtain_point(g,hit)) then
         t_hit=root; ok=.true.; return
      end if
      tprev=root+max(root_skip,step_min*1.0e-3_dp)
      if (tprev >= t_max) return
      fprev=residual_ray(g,tprev,p0,n)
   else
      if (tcur >= t_max) return
      tprev=tcur; fprev=fcur
   end if
end do
end function ray_curtain_first_intersection

subroutine sample_footprint(g,state,p)
type(curtain_geometry), intent(in) :: g
integer(kind=8), intent(inout) :: state
real(dp), intent(out) :: p(3)
integer :: q,k,row,hem
real(dp) :: phi,lamb
q=1+int(urand(state)*real(g%n_valid,dp))
if (q > g%n_valid) q=g%n_valid
k=g%valid_idx(q); hem=g%valid_hem(q)
if (hem < 0) then
   row=1
else
   row=2
end if
phi=g%coverage_phi(k)
lamb=g%coverage_lambda_surface(row,k)
call magnetic_position(g,phi,lamb,p)
p=normalize(p)
end subroutine sample_footprint

subroutine propagate_one_photon(state, geom, L39, B12, T_keV, R_ns, R_ns_km, r_stop_factor, &
     step_coarse, tol, advance_factor, Lx, alpha_disk, mstar, disk_rho_max_factor, forced, minw, taucrit, res)
integer(kind=8), intent(inout) :: state
type(curtain_geometry), intent(in) :: geom
real(dp), intent(in) :: L39,B12,T_keV,R_ns,R_ns_km,r_stop_factor,step_coarse,tol,advance_factor
real(dp), intent(in) :: Lx,alpha_disk,mstar,disk_rho_max_factor,minw,taucrit
logical, intent(in) :: forced
type(photon_result), intent(out) :: res
real(dp) :: Rm,r_stop,disk_rho_max_model,advance_dist
real(dp) :: pos(3),normal_out(3),dir(3),hit(3),nnew(3),layer_normal(3)
real(dp) :: E0,E,weight,initial_weight,mu_launch,t_stop,t_ns,t_m,t_disk,t_event,tau,pscat
real(dp) :: disk_rho,side_sign,E_before,E_after,e_dep,offset_sign,cospsi,mu_n
integer :: step,event,n_scat,n_internal,entry_sign,exit_side,status
logical :: ok_stop,ok_ns,ok_m,ok_disk,did_scatter

res=photon_result()
Rm=geom%rm_model
r_stop=r_stop_factor*Rm
disk_rho_max_model=disk_rho_max_factor*Rm
advance_dist=max(advance_factor*tol,1.0e-7_dp)
call sample_footprint(geom,state,pos)
normal_out=normalize(pos)
E0=sample_thermal_energy(T_keV,state)
E=E0
call sample_uniform_outward_hemisphere_direction(normal_out,state,dir,mu_launch)
weight=mu_launch
initial_weight=weight
n_scat=0

do step=1,max_steps_default
   if (weight < minw) exit
   ok_stop=ray_sphere_first_intersection(pos,dir,r_stop,advance_dist,t_stop)
   ok_ns=ray_sphere_first_intersection(pos,dir,R_ns,advance_dist,t_ns)
   ok_disk=ray_flared_disk_first_intersection(pos,dir,Rm,Lx,R_ns_km,alpha_disk,mstar, &
        advance_dist,merge(t_stop,huge(1.0_dp),ok_stop),step_coarse,tol,disk_rho_max_model, &
        t_disk,disk_rho,side_sign)
   ok_m=ray_curtain_first_intersection(geom,pos,dir,advance_dist, &
        merge(t_stop,huge(1.0_dp),ok_stop),step_coarse,tol,t_m)

   if (.not.ok_stop .and. .not.ok_ns .and. .not.ok_m .and. .not.ok_disk) exit
   call choose_event(ok_stop,t_stop,ok_ns,t_ns,ok_m,t_m,ok_disk,t_disk,t_event,event)
   hit=pos+t_event*dir

   select case(event)
   case(1)
      exit
   case(2)
      res%disk_absorbed=.true.
      if (n_scat == 0) then
         res%disk_component=1
      else
         res%disk_component=2
      end if
      res%disk_rho_model=disk_rho
      res%disk_phi=atan2(hit(2),hit(1))
      if (res%disk_phi < 0.0_dp) res%disk_phi=res%disk_phi+2.0_dp*pi
      res%disk_mu=disk_surface_mu(hit,dir,side_sign,Lx,R_ns_km,alpha_disk,mstar)
      exit
   case(3)
      pos=R_ns*normalize(hit)
      call sample_uniform_outward_hemisphere_direction(normalize(pos),state,dir,mu_launch)
      weight=weight*mu_launch
      pos=pos+advance_dist*dir
   case(4)
      tau=tau_thomson_normal(geom,hit)
      if (tau <= 0.0_dp) then
         pos=hit+advance_dist*dir
         cycle
      end if
      call curtain_normal(geom,hit,layer_normal)
      mu_n=max(abs(dot3(dir,layer_normal)),1.0e-12_dp)
      if (tau <= taucrit) then
         pscat=1.0_dp-exp(-tau/mu_n)
         if (.not.forced) then
            did_scatter=(urand(state) < pscat)
         else
            weight=weight*pscat
            did_scatter=(pscat > 0.0_dp .and. weight >= minw)
         end if
         if (did_scatter) then
            E_before=E
            call sample_klein_nishina_direction(E_before,dir,state,nnew,cospsi)
            E=compton_scatter_energy(E_before,cospsi)
            dir=nnew
            n_scat=n_scat+1
         end if
         pos=hit+advance_dist*dir
      else
         call transport_through_layer(tau,dir,E,layer_normal,state,nnew,E_after,n_internal,e_dep, &
                                      entry_sign,exit_side,status)
         E=E_after
         dir=normalize(nnew)
         n_scat=n_scat+n_internal
         if (exit_side == 2) then
            offset_sign=real(entry_sign,dp)
         else
            offset_sign=-real(entry_sign,dp)
         end if
         pos=hit+offset_sign*advance_dist*layer_normal
         if (status /= 0) exit
      end if
   end select
end do
res%launch_weight=initial_weight
res%packet_weight=weight
res%e0_kev=E0
res%efinal_kev=E
end subroutine propagate_one_photon

subroutine transport_through_layer(tau0,direction_in,E_in,normal,state,n_out,E_out,n_internal,e_dep, &
                                   entry_sign,exit_side,status)
real(dp), intent(in) :: tau0,direction_in(3),E_in,normal(3)
integer(kind=8), intent(inout) :: state
real(dp), intent(out) :: n_out(3),E_out,e_dep
integer, intent(out) :: n_internal,entry_sign,exit_side,status
real(dp) :: n(3),m(3),E,mu,u,tau_to_boundary,delta_tau,nnew(3),cospsi,E_before
logical :: entered_from_low
integer :: i,exit_boundary
status=0
n=normalize(direction_in); m=normalize(normal); E=E_in
mu=dot3(n,m)
if (abs(mu) < 1.0e-14_dp) mu=merge(1.0e-14_dp,-1.0e-14_dp,urand(state)<0.5_dp)
entry_sign=merge(1,-1,mu>0.0_dp)
entered_from_low=(mu>0.0_dp)
u=merge(0.0_dp,tau0,mu>0.0_dp)
n_internal=0

do i=1,max_internal_scatters_default
   mu=dot3(n,m)
   if (abs(mu) < 1.0e-14_dp) mu=merge(1.0e-14_dp,-1.0e-14_dp,urand(state)<0.5_dp)
   if (mu > 0.0_dp) then
      tau_to_boundary=max((tau0-u)/mu,0.0_dp)
   else
      tau_to_boundary=max(u/(-mu),0.0_dp)
   end if
   delta_tau=-log(urand(state))
   if (delta_tau >= tau_to_boundary) then
      exit_boundary=merge(2,1,mu>0.0_dp)
      if (entered_from_low) then
         exit_side=merge(1,2,exit_boundary==1)
      else
         exit_side=merge(1,2,exit_boundary==2)
      end if
      n_out=n; E_out=E; e_dep=E_in-E
      return
   end if
   u=min(max(u+mu*delta_tau,0.0_dp),tau0)
   E_before=E
   call sample_klein_nishina_direction(E_before,n,state,nnew,cospsi)
   E=compton_scatter_energy(E_before,cospsi)
   n=nnew
   n_internal=n_internal+1
end do
status=1
exit_side=1
n_out=n; E_out=E; e_dep=E_in-E
end subroutine transport_through_layer

subroutine sample_klein_nishina_direction(E_keV,n_in,state,n_out,cospsi)
real(dp), intent(in) :: E_keV,n_in(3)
integer(kind=8), intent(inout) :: state
real(dp), intent(out) :: n_out(3),cospsi
real(dp) :: nin(3),e1(3),e2(3),ref(3),mu,phi,st,eps,q,pdf
nin=normalize(n_in)
eps=max(E_keV,0.0_dp)/mec2_kev
if (abs(nin(3)) < 0.9_dp) then
   ref=[0.0_dp,0.0_dp,1.0_dp]
else
   ref=[1.0_dp,0.0_dp,0.0_dp]
end if
e1=normalize(cross3(ref,nin))
e2=normalize(cross3(nin,e1))
do
   mu=2.0_dp*urand(state)-1.0_dp
   q=1.0_dp/(1.0_dp+eps*(1.0_dp-mu))
   pdf=0.5_dp*q*q*(q+1.0_dp/q-(1.0_dp-mu*mu))
   if (urand(state) <= pdf) exit
end do
phi=2.0_dp*pi*urand(state)
st=sqrt(max(0.0_dp,1.0_dp-mu*mu))
n_out=mu*nin+st*cos(phi)*e1+st*sin(phi)*e2
n_out=normalize(n_out)
cospsi=mu
end subroutine sample_klein_nishina_direction

subroutine choose_event(os,ts,on,tn,om,tm,od,td,t,event)
logical, intent(in) :: os,on,om,od
real(dp), intent(in) :: ts,tn,tm,td
real(dp), intent(out) :: t
integer, intent(out) :: event
t=huge(1.0_dp); event=0
if (os .and. ts<t) then; t=ts; event=1; end if
if (od .and. td<t) then; t=td; event=2; end if
if (on .and. tn<t) then; t=tn; event=3; end if
if (om .and. tm<t) then; t=tm; event=4; end if
end subroutine choose_event

function locate_bin(edges,nbins,x) result(idx)
integer, intent(in) :: nbins
real(dp), intent(in) :: edges(nbins+1),x
integer :: idx,i
idx=0
do i=1,nbins
   if ((x>=edges(i) .and. x<edges(i+1)) .or. (i==nbins .and. x==edges(i+1))) then
      idx=i; return
   end if
end do
end function locate_bin

function opt_real(x,default) result(y)
real(dp), intent(in), optional :: x
real(dp), intent(in) :: default
real(dp) :: y
if (present(x)) then; y=x; else; y=default; end if
end function opt_real

function opt_int(x,default) result(y)
integer, intent(in), optional :: x
integer, intent(in) :: default
integer :: y
if (present(x)) then; y=x; else; y=default; end if
end function opt_int

function opt_logical(x,default) result(y)
logical, intent(in), optional :: x
logical, intent(in) :: default
logical :: y
if (present(x)) then; y=x; else; y=default; end if
end function opt_logical

function photon_seed(seed,iphot) result(state)
integer, intent(in) :: seed,iphot
integer(kind=8) :: state
state=mod(int(seed,8)+104729_8*int(iphot,8),rng_modulus-1_8)+1_8
end function photon_seed

function urand(state) result(u)
integer(kind=8), intent(inout) :: state
real(dp) :: u
state=mod(rng_multiplier*state,rng_modulus)
if (state <= 0_8) state=state+rng_modulus-1_8
u=real(state,dp)/real(rng_modulus,dp)
u=max(1.0e-15_dp,min(1.0_dp-1.0e-15_dp,u))
end function urand

function dot3(a,b) result(c)
real(dp), intent(in) :: a(3),b(3)
real(dp) :: c
c=sum(a*b)
end function dot3

function norm3(a) result(n)
real(dp), intent(in) :: a(3)
real(dp) :: n
n=sqrt(sum(a*a))
end function norm3

function normalize(a) result(b)
real(dp), intent(in) :: a(3)
real(dp) :: b(3),n
n=norm3(a)
if (n<=0.0_dp) then
   b=[1.0_dp,0.0_dp,0.0_dp]
else
   b=a/n
end if
end function normalize

function cross3(a,b) result(c)
real(dp), intent(in) :: a(3),b(3)
real(dp) :: c(3)
c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
end function cross3

subroutine orthonormal_basis(normal_out,e1,e2,e3)
real(dp), intent(in) :: normal_out(3)
real(dp), intent(out) :: e1(3),e2(3),e3(3)
real(dp) :: ref(3)
e3=normalize(normal_out)
if (abs(e3(3))<0.9_dp) then
   ref=[0.0_dp,0.0_dp,1.0_dp]
else
   ref=[1.0_dp,0.0_dp,0.0_dp]
end if
e1=cross3(ref,e3)
if (norm3(e1)<1.0e-14_dp) e1=cross3([0.0_dp,1.0_dp,0.0_dp],e3)
e1=normalize(e1)
e2=normalize(cross3(e3,e1))
end subroutine orthonormal_basis

subroutine sample_uniform_outward_hemisphere_direction(normal_out,state,direction,mu_check)
real(dp), intent(in) :: normal_out(3)
integer(kind=8), intent(inout) :: state
real(dp), intent(out) :: direction(3),mu_check
real(dp) :: e1(3),e2(3),e3(3),mu,phi,st
call orthonormal_basis(normal_out,e1,e2,e3)
mu=urand(state)
phi=2.0_dp*pi*urand(state)
st=sqrt(max(0.0_dp,1.0_dp-mu*mu))
direction=st*cos(phi)*e1+st*sin(phi)*e2+mu*e3
direction=normalize(direction)
mu_check=max(0.0_dp,dot3(direction,normalize(normal_out)))
end subroutine sample_uniform_outward_hemisphere_direction

function ray_sphere_first_intersection(p0,n_dir,R,t_min,t_hit) result(ok)
real(dp), intent(in) :: p0(3),n_dir(3),R,t_min
real(dp), intent(out) :: t_hit
logical :: ok
real(dp) :: n(3),b,c,disc,s,t1,t2
n=normalize(n_dir)
b=2.0_dp*dot3(p0,n)
c=dot3(p0,p0)-R*R
disc=b*b-4.0_dp*c
ok=.false.; t_hit=0.0_dp
if (disc<0.0_dp) return
s=sqrt(disc); t1=(-b-s)/2.0_dp; t2=(-b+s)/2.0_dp
if (t1>t_min .and. t2>t_min) then
   t_hit=min(t1,t2); ok=.true.
else if (t1>t_min) then
   t_hit=t1; ok=.true.
else if (t2>t_min) then
   t_hit=t2; ok=.true.
end if
end function ray_sphere_first_intersection

function disk_H_coeff(Lx,Rnskm,alpha,mstar) result(c)
real(dp), intent(in) :: Lx,Rnskm,alpha,mstar
real(dp) :: c,L37,Rns_cm,R6
L37=Lx/1.0e37_dp
Rns_cm=Rnskm*1.0e5_dp
R6=Rns_cm/1.0e6_dp
c=0.08_dp*alpha**(-0.1_dp)*L37**(3.0_dp/20.0_dp)*mstar**(-21.0_dp/40.0_dp)* &
  R6**(3.0_dp/20.0_dp)*(Rns_cm/1.0e8_dp)**(1.0_dp/8.0_dp)
end function disk_H_coeff

function disk_H_model(rho,Lx,Rnskm,alpha,mstar) result(H)
real(dp), intent(in) :: rho,Lx,Rnskm,alpha,mstar
real(dp) :: H
if (rho<=0.0_dp) then
   H=0.0_dp
else
   H=disk_H_coeff(Lx,Rnskm,alpha,mstar)*rho**(9.0_dp/8.0_dp)
end if
end function disk_H_model

function disk_g(t,p,n,side,Lx,Rnskm,alpha,mstar) result(g)
real(dp), intent(in) :: t,p(3),n(3),side,Lx,Rnskm,alpha,mstar
real(dp) :: g,q(3),rho,H
q=p+t*n
rho=sqrt(q(1)*q(1)+q(2)*q(2))
H=disk_H_model(rho,Lx,Rnskm,alpha,mstar)
g=q(3)-side*H
end function disk_g

function disk_surface_mu(pos,dir,side,Lx,Rnskm,alpha,mstar) result(mu)
real(dp), intent(in) :: pos(3),dir(3),side,Lx,Rnskm,alpha,mstar
real(dp) :: mu,rho,H,dHdr,n_surf(3)
rho=sqrt(pos(1)*pos(1)+pos(2)*pos(2))
if (rho<=0.0_dp) then
   mu=abs(dir(3)); return
end if
H=disk_H_model(rho,Lx,Rnskm,alpha,mstar)
dHdr=(9.0_dp/8.0_dp)*H/rho
n_surf(1)=-side*dHdr*pos(1)/rho
n_surf(2)=-side*dHdr*pos(2)/rho
n_surf(3)=side
n_surf=normalize(n_surf)
mu=abs(dot3(normalize(dir),n_surf))
mu=max(0.0_dp,min(1.0_dp,mu))
end function disk_surface_mu

function bisect_disk_root(p,n,side,a,b,fa,fb,Lx,Rnskm,alpha,mstar,tol) result(root)
real(dp), intent(in) :: p(3),n(3),side,a,b,fa,fb,Lx,Rnskm,alpha,mstar,tol
real(dp) :: root,lo,hi,flo,mid,fmid
integer :: i
lo=a; hi=b; flo=fa
do i=1,100
   mid=0.5_dp*(lo+hi)
   fmid=disk_g(mid,p,n,side,Lx,Rnskm,alpha,mstar)
   if ((hi-lo)<=tol .or. fmid==0.0_dp) exit
   if (flo*fmid<=0.0_dp) then
      hi=mid
   else
      lo=mid; flo=fmid
   end if
end do
root=0.5_dp*(lo+hi)
end function bisect_disk_root

function adaptive_step(p,n,t,step_min) result(dt)
real(dp), intent(in) :: p(3),n(3),t,step_min
real(dp) :: dt,q(3),r
q=p+t*n; r=norm3(q)
dt=min(max(5.0_dp,step_min),max(step_min,0.005_dp*max(r,1.0_dp)))
end function adaptive_step

function ray_flared_disk_first_intersection(p0,n_dir,Rm,Lx,Rnskm,alpha,mstar,t_min,t_max, &
                                            step_coarse,tol,rho_max_model,t_hit,rho_hit,side_hit) result(ok)
real(dp), intent(in) :: p0(3),n_dir(3),Rm,Lx,Rnskm,alpha,mstar,t_min,t_max,step_coarse,tol,rho_max_model
real(dp), intent(out) :: t_hit,rho_hit,side_hit
logical :: ok
real(dp) :: n(3),side,tprev,fprev,t,fcur,root,q(3),rho,dt,best_t,best_rho,best_side
integer :: iside
ok=.false.; best_t=huge(1.0_dp); best_rho=0.0_dp; best_side=1.0_dp
n=normalize(n_dir)
do iside=1,2
   side=merge(1.0_dp,-1.0_dp,iside==1)
   tprev=max(0.0_dp,t_min)
   fprev=disk_g(tprev,p0,n,side,Lx,Rnskm,alpha,mstar)
   if (abs(fprev)<1.0e-14_dp) then
      tprev=tprev+max(tol,1.0e-7_dp)
      fprev=disk_g(tprev,p0,n,side,Lx,Rnskm,alpha,mstar)
   end if
   do while (tprev<t_max-1.0e-15_dp)
      dt=adaptive_step(p0,n,tprev,step_coarse)
      t=min(tprev+dt,t_max)
      fcur=disk_g(t,p0,n,side,Lx,Rnskm,alpha,mstar)
      if (fprev==0.0_dp .or. fcur==0.0_dp .or. fprev*fcur<0.0_dp) then
         if (fprev==0.0_dp) then
            root=tprev
         else if (fcur==0.0_dp) then
            root=t
         else
            root=bisect_disk_root(p0,n,side,tprev,t,fprev,fcur,Lx,Rnskm,alpha,mstar,tol)
         end if
         if (root>t_min) then
            q=p0+root*n
            rho=sqrt(q(1)*q(1)+q(2)*q(2))
            if (rho+1.0e-12_dp>=Rm .and. rho<=rho_max_model+1.0e-12_dp) then
               if (root<best_t) then
                  best_t=root; best_rho=rho; best_side=side; ok=.true.
               end if
               exit
            end if
         end if
      end if
      tprev=t; fprev=fcur
   end do
end do
t_hit=best_t; rho_hit=best_rho; side_hit=best_side
end function ray_flared_disk_first_intersection

function Rm_from_LB(L39,B12) result(Rm)
real(dp), intent(in) :: L39,B12
real(dp) :: Rm
Rm=35.0_dp*1.4_dp**(1.0_dp/7.0_dp)*B12**(4.0_dp/7.0_dp)*L39**(-2.0_dp/7.0_dp)
end function Rm_from_LB

subroutine init_planck_cdf()
integer :: i
real(dp) :: s
if (planck_ready) return
s=0.0_dp
do i=1,planck_nmax
   s=s+1.0_dp/real(i,dp)**3
   planck_cdf(i)=s
end do
planck_cdf=planck_cdf/planck_cdf(planck_nmax)
planck_ready=.true.
end subroutine init_planck_cdf

function sample_thermal_energy(T_keV,state) result(E)
real(dp), intent(in) :: T_keV
integer(kind=8), intent(inout) :: state
real(dp) :: E,u,x
integer :: n
u=urand(state)
n=1
do while (n<planck_nmax .and. planck_cdf(n)<u)
   n=n+1
end do
x=gamma_k3_scale(1.0_dp/real(n,dp),state)
E=T_keV*x
end function sample_thermal_energy

function gamma_k3_scale(scale,state) result(x)
real(dp), intent(in) :: scale
integer(kind=8), intent(inout) :: state
real(dp) :: x
x=-scale*log(urand(state)*urand(state)*urand(state))
end function gamma_k3_scale

function compton_scatter_energy(E,cospsi) result(Eout)
real(dp), intent(in) :: E,cospsi
real(dp) :: Eout,c
c=max(-1.0_dp,min(1.0_dp,cospsi))
Eout=E/(1.0_dp+(E/mec2_kev)*(1.0_dp-c))
end function compton_scatter_energy

subroutine viscous_disk_profile(R_cm,n,Lx,Rnskm,mstar,Rinner,F,T)
integer, intent(in) :: n
real(dp), intent(in) :: R_cm(n),Lx,Rnskm,mstar,Rinner
real(dp), intent(out) :: F(n),T(n)
real(dp) :: Rns_cm,M_g,mdot,boundary
integer :: i
Rns_cm=Rnskm*1.0e5_dp
M_g=mstar*m_sun_g
mdot=Lx*Rns_cm/(g_cgs*M_g)
F=0.0_dp; T=0.0_dp
do i=1,n
   if (R_cm(i)>0.0_dp .and. R_cm(i)>=Rinner) then
      boundary=max(1.0_dp-sqrt(Rinner/R_cm(i)),0.0_dp)
      F(i)=3.0_dp*g_cgs*M_g*mdot/(8.0_dp*pi*R_cm(i)**3)*boundary
      if (F(i)>0.0_dp) T(i)=(F(i)/sigma_sb_cgs)**0.25_dp
   end if
end do
end subroutine viscous_disk_profile

subroutine test_xray_disk_mc_profile()
implicit none
integer, parameter :: nbins=40
integer :: ierr,i
real(dp) :: profile(nbins,9)
call xray_disk_mc_profile( &
   N=100000, L39=0.04_dp, B12=3.0_dp, T_keV=5.0_dp, nbins=nbins, &
   profile=profile, ierr=ierr, disk_rho_max_factor=50.0_dp, &
   dipole_tilt_rad=50.0_dp*pi/180.0_dp, &
   phase_cycles=0.0_dp, &
   H_over_Rm=0.06_dp, spin_period_s=5.0_dp, v_initial_cms=2.0e7_dp, &
   velocity_floor_cms=1.0e5_dp, coverage_azimuth_bins=720, &
   coverage_points_along_stream=500)
if (ierr/=0) stop 'xray_disk_mc_profile failed'
open(unit=10,file='./res/res_new    ',status='replace')
do i=1,nbins
   write(10,'(9E20.10)') profile(i,1),profile(i,2),profile(i,3),profile(i,4), &
        profile(i,5),profile(i,6),profile(i,7),profile(i,8),profile(i,9)
end do
close(10)
end subroutine test_xray_disk_mc_profile

end module xray_disk_mc_module





