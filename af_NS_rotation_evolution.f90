!===============================================================================
! ns_population_torque_models.f90
! OpenMP-ready, final-only population evolution of accreting magnetised neutron stars in
! Be/X-ray binaries.
! IMPORTANT FOR YOUR MAIN CODE:
!   The module name is intentionally kept as
!       ns_population_torque_models
!   and the zero-argument wrapper is still
!       call run_ns_population()
!   so the calling convention does not need to change.
!
! Parallelisation:
!   The population loop over stars is parallelised with OpenMP. Each thread
!   generates one accretion history, integrates one star, stores only the final
!   population summary, and immediately deallocates the history.
!
! Compilation with OpenMP:
!   gfortran -O3 -fopenmp -c ns_population_torque_models.f90 -o obj/ns_population_torque_models.o
!
! Linking must also include -fopenmp:
!   gfortran -O3 -fopenmp -o main obj/ns_population_torque_models.o obj/c.o ...
!
! Run with, for example:
!   OMP_NUM_THREADS=8 ./main
! Memory:
!   This version does not allocate orbit-by-orbit track arrays. It stores only
!   initial/final values and summary diagnostics for each star.
! Output:
!   All non-data text lines begin with "#", so stdout can be used directly as a
!   data stream for plotting scripts. Table data rows do not start with "#".
!===============================================================================
module ns_population_torque_models
implicit none
private
public :: run_ns_population
public :: run_population
public :: NSParams
public :: AccretionHistory
public :: PopulationResult
public :: TORQUE_BA2021
public :: TORQUE_LINEAR_FASTNESS
public :: TORQUE_SATURATED_FASTNESS
public :: TORQUE_TANH_FASTNESS
public :: OBL_ACCRETOR_ONLY
public :: OBL_POSITIVE_INTERACTION
public :: OBL_SAME_AS_SPIN

integer, parameter :: dp = selected_real_kind(15, 300)

real(dp), parameter :: pi   = 3.1415926535897932384626433832795_dp
real(dp), parameter :: Gcgs = 6.67430e-8_dp
real(dp), parameter :: ccgs = 2.99792458e10_dp
real(dp), parameter :: Msun = 1.98847e33_dp
real(dp), parameter :: yr   = 365.25_dp * 24.0_dp * 3600.0_dp
real(dp), parameter :: Mdot_Edd_solar = 1.6e17_dp

integer, parameter :: TORQUE_BA2021               = 1
integer, parameter :: TORQUE_LINEAR_FASTNESS      = 2
integer, parameter :: TORQUE_SATURATED_FASTNESS   = 3
integer, parameter :: TORQUE_TANH_FASTNESS        = 4

integer, parameter :: OBL_ACCRETOR_ONLY           = 1
integer, parameter :: OBL_POSITIVE_INTERACTION    = 2
integer, parameter :: OBL_SAME_AS_SPIN            = 3

  type :: NSParams
     ! Neutron star
     real(dp) :: M0 = 1.5_dp * Msun
     real(dp) :: R  = 1.25e6_dp
     real(dp) :: I_per_Msun = 1.0e45_dp
     real(dp) :: xi = 0.5_dp

     ! Initial NS parameters
     real(dp) :: P0 = 100.0_dp
     real(dp) :: B0 = 1.0e13_dp
     integer  :: B_convention = 2
     ! B_convention = 1: equatorial field, mu = B R^3
     ! B_convention = 2: polar field,      mu = B R^3 / 2

     ! Initial geometry; overwritten star-by-star in population runs
     real(dp) :: alpha0_deg = 60.0_dp
     real(dp) :: chi0_deg   = 60.0_dp

     ! Fallback/default accretion rate.
     ! In Be-outburst mode this is used as the mean Type-I outburst Mdot unless
     ! Mdot_typeI_mean is explicitly set positive.
     real(dp) :: Mdot0 = 1.0e15_dp
     real(dp) :: eta   = 0.99_dp

     !---------------------------------------------------------------------------
     ! Be/XRB outburst model
     !---------------------------------------------------------------------------

     logical  :: use_be_outbursts = .true.

     ! Orbital period and Type-I outburst duration.
     real(dp) :: Porb_yr = 0.27_dp
     real(dp) :: outburst_duration_yr = 0.03_dp

     ! Probability that a regular Type-I outburst occurs on a given orbit.
     ! 1.0 means a Type-I outburst every orbit.
     real(dp) :: typeI_probability_per_orbit = 1.0_dp

     ! Quiescent accretion between Type-I outbursts.
     ! Set to 0 for no accretion torque outside outbursts.
     real(dp) :: Mdot_quiescent = 0.0_dp

     ! Correlated angular-momentum direction of captured material.
     ! Orbital angular momentum is the z-axis.
     ! OU process is defined in the small-tilt plane:
     !   u = (u_x, u_y)
     !   j = normalize( u_x, u_y, 1 )
     real(dp) :: beta_mean_deg = 15.0_dp
     real(dp) :: beta_rms_deg  = 10.0_dp
     real(dp) :: corr_time_orbits = 20.0_dp

     ! If true, the azimuth of the mean Be-disk tilt is randomized for each star.
     logical  :: random_mean_tilt_azimuth = .true.
     real(dp) :: beta_mean_phi_deg = 0.0_dp

     ! Outburst Mdot model.
     ! log(Mdot/Mdot_typeI_mean) follows an OU process.
     real(dp) :: Mdot_typeI_mean = -1.0_dp
     real(dp) :: Mdot_typeI_sigma_log = 0.5_dp
     real(dp) :: Mdot_corr_time_orbits = 5.0_dp

     ! Optional rare giant outbursts.
     logical  :: include_typeII = .false.
     real(dp) :: typeII_probability_per_orbit = 1.0e-3_dp
     real(dp) :: typeII_boost = 30.0_dp
     real(dp) :: typeII_duration_orbits = 5.0_dp

     ! Numerical substeps inside active/quiescent parts of each orbit.
     integer  :: active_substeps = 10
     integer  :: quiescent_substeps = 1

     !---------------------------------------------------------------------------
     ! Torque prescription
     !---------------------------------------------------------------------------

     integer  :: torque_model = TORQUE_SATURATED_FASTNESS
     real(dp) :: omega_c = 1.0_dp
     real(dp) :: tanh_width = 0.25_dp

     ! For fastness torque models this is usually false, to avoid double-counting
     ! spin-down. It may be used for BA2021-like runs.
     logical  :: include_BA_magnetospheric_braking = .false.

     ! Geometry/obliquity torque scaling
     integer  :: obliquity_torque_mode = OBL_ACCRETOR_ONLY
     real(dp) :: interaction_decay_power = 2.0_dp

     ! Physics switches
     logical :: include_accretion_obliquity = .true.
     logical :: include_alpha_alignment = .true.
     logical :: include_pulsar_torque = .true.
     logical :: include_mass_growth = .false.
     logical :: include_magnetic_burial = .false.
     logical :: include_parfrey_enhancement = .false.

     ! Angle clipping
     real(dp) :: chi_floor_deg   = 1.0e-5_dp
     real(dp) :: chi_ceiling_deg = 89.99999_dp
  end type NSParams


  type :: AccretionHistory
     integer :: n_orbits = 0
     real(dp) :: Porb_yr = 0.0_dp
     real(dp) :: outburst_duration_yr = 0.0_dp

     ! jvec(:,n): direction of captured angular momentum in orbit n.
     real(dp), allocatable :: jvec(:,:)      ! (3,n_orbits)

     ! Mdot during the active part of orbit n.
     real(dp), allocatable :: mdot(:)

     ! Type-I and Type-II activity flags.
     logical, allocatable :: is_typeI(:)
     logical, allocatable :: is_typeII(:)

     ! Small-tilt coordinates used to generate jvec.
     real(dp), allocatable :: u_tilt(:,:)    ! (2,n_orbits)
  end type AccretionHistory


  type :: TorqueTerms
     real(dp) :: d(3)
     real(dp) :: cos_alpha
     real(dp) :: alpha
     real(dp) :: M
     real(dp) :: I
     real(dp) :: Mdot
     real(dp) :: mu
     real(dp) :: rm
     real(dp) :: rco
     real(dp) :: rlc
     real(dp) :: omega_s
     logical  :: accretor
     logical  :: inside_lc
     logical  :: outburst_active
     integer  :: orbit_index
     real(dp) :: N0
     real(dp) :: n_fastness
     real(dp) :: interaction_factor
     real(dp) :: N_spin
     real(dp) :: Nmag_BA
     real(dp) :: Kpsr
     real(dp) :: N_geom
  end type TorqueTerms


  type :: PopulationResult
     integer :: n_stars = 0
     integer :: n_steps = 0   ! number of stored orbit-boundary steps

     real(dp), allocatable :: alpha0_deg(:), chi0_deg(:), P0_s(:)
     real(dp), allocatable :: alpha_final_deg(:), chi_final_deg(:), P_final_s(:)
     real(dp), allocatable :: delta_alpha_deg(:), delta_chi_deg(:), delta_P_s(:)

     real(dp), allocatable :: mean_q_alpha(:)
     real(dp), allocatable :: mean_mdot(:)
     real(dp), allocatable :: final_rm_over_rco(:), final_omega_s(:)
     real(dp), allocatable :: final_n_fastness(:), final_interaction_factor(:)

     ! Full tracks are deliberately not allocated in this final-only version.
     ! These components remain for compatibility with old output routines.
     real(dp), allocatable :: t_yr(:)
     real(dp), allocatable :: P_track(:,:), alpha_track(:,:), chi_track(:,:)
     real(dp), allocatable :: rm_over_rco_track(:,:), omega_s_track(:,:)
     real(dp), allocatable :: q_alpha_track(:,:), mdot_track(:,:)
  end type PopulationResult

contains
!===============================================================================
! Public wrapper
!===============================================================================
  subroutine run_ns_population()
  implicit none
  type(NSParams) :: p
  type(PopulationResult) :: res

  integer :: n_stars
  integer :: seed
  real(dp) :: Lx
  real(dp) :: t_end_yr

    !---------------------------------------------------------------------------
    ! User configuration block
    !---------------------------------------------------------------------------
    n_stars = 200
    seed = 1

    Lx = 1.0e35_dp

    p%P0 = 100.0_dp !1.0_dp
    p%B0 = 1.0e13_dp
    p%Mdot0 = mdot_from_lx(Lx, p%M0, p%R, 1.0_dp)
    p%eta = 0.99_dp

    ! Be/XRB accretion history
    p%use_be_outbursts = .true.

    p%Porb_yr = 0.27_dp
    p%outburst_duration_yr = 0.03_dp
    p%Mdot_quiescent = 0.0_dp
    p%typeI_probability_per_orbit = 0.3_dp

    p%beta_mean_deg = 15.0_dp
    p%beta_rms_deg  = 10.0_dp
    p%corr_time_orbits = 20.0_dp
    p%random_mean_tilt_azimuth = .true.

    ! If negative, the code uses p%Mdot0.
    p%Mdot_typeI_mean = -1.0_dp
    p%Mdot_typeI_sigma_log = 0.5_dp
    p%Mdot_corr_time_orbits = 5.0_dp

    p%include_typeII = .false.
    p%typeII_probability_per_orbit = 1.0e-3_dp
    p%typeII_boost = 30.0_dp
    p%typeII_duration_orbits = 5.0_dp

    p%active_substeps = 10
    p%quiescent_substeps = 1

    !===================================
    ! Torque model
    !   TORQUE_BA2021
    !   TORQUE_LINEAR_FASTNESS
    !   TORQUE_SATURATED_FASTNESS
    !   TORQUE_TANH_FASTNESS
    !===================================!
    p%torque_model = TORQUE_SATURATED_FASTNESS

    !   OBL_ACCRETOR_ONLY
    !   OBL_POSITIVE_INTERACTION
    !   OBL_SAME_AS_SPIN
    p%obliquity_torque_mode = OBL_ACCRETOR_ONLY

    p%omega_c = 1.0_dp
    p%tanh_width = 0.25_dp
    p%include_BA_magnetospheric_braking = .true.

    p%include_pulsar_torque = .true.
    p%include_accretion_obliquity = .true.
    p%include_alpha_alignment = .true.

    p%include_mass_growth = .false.
    p%include_magnetic_burial = .false.
    p%include_parfrey_enhancement = .false.

    t_end_yr = 4.0e5_dp

    call run_population(p, n_stars, t_end_yr, seed, res)

    call print_configuration(p, n_stars, t_end_yr, res%n_steps)
    call print_population_table(res)
    call print_spin_statistics(res)
    call print_histograms(res, 12, 9)

    ! Optional output for plotting later:
    ! call write_population_ascii("ns_be_population_summary.dat", res)
    ! Do not call in final-only mode:
    ! call write_tracks_ascii("ns_be_population_tracks.dat", res)
  end subroutine run_ns_population


!===============================================================================
! Population driver
!===============================================================================
  subroutine run_population(p_in, n_stars, t_end_yr, seed, res)
    implicit none

    type(NSParams), intent(in) :: p_in
    integer, intent(in) :: n_stars
    real(dp), intent(in) :: t_end_yr
    integer, intent(in) :: seed
    type(PopulationResult), intent(out) :: res

    type(NSParams) :: p
    type(AccretionHistory) :: hist

    integer :: i, n_orbits
    real(dp) :: alpha0, chi0
    real(dp) :: qsum

    if (p_in%use_be_outbursts) then
       n_orbits = ceiling_to_int(t_end_yr / p_in%Porb_yr)
       if (n_orbits < 1) n_orbits = 1
    else
       n_orbits = ceiling_to_int(t_end_yr / max(p_in%Porb_yr, 1.0e-10_dp))
       if (n_orbits < 1) n_orbits = 1
    end if

    ! Final-only storage. The number of orbit steps is not used for allocation.
    call allocate_population_result(res, n_stars, 0)

    !---------------------------------------------------------------------------
    ! Parallel block:
    ! Each star is independent. Its accretion history is generated, used, and
    ! deallocated inside the same iteration. Therefore histories for all stars are
    ! never stored simultaneously.
    !
    ! The standard Fortran random_number/random_seed generator is global. To avoid
    ! OpenMP races, all random generation for a given star is done inside a
    ! critical section. Integration itself is still parallel.
    !---------------------------------------------------------------------------

    !$omp parallel do default(shared) private(i,p,hist,alpha0,chi0,qsum) schedule(dynamic)
    do i = 1, n_stars

       p = p_in

       !$omp critical(rng_history_generation)
       call init_random_seed(seed + 1000003*i)

       call generate_accretion_history(p, t_end_yr, hist)

       alpha0 = draw_alpha_isotropic_deg(0.0_dp, 180.0_dp)
       chi0   = draw_chi_isotropic_folded_deg(0.0_dp, 90.0_dp)
       !$omp end critical(rng_history_generation)

       p%alpha0_deg = alpha0
       p%chi0_deg = chi0

       call integrate_one_star(p, hist, i, res, qsum)

       call deallocate_history(hist)

    end do
    !$omp end parallel do

  end subroutine run_population


  subroutine allocate_population_result(res, n_stars, n_steps)
    implicit none

    type(PopulationResult), intent(out) :: res
    integer, intent(in) :: n_stars, n_steps

    res%n_stars = n_stars
    res%n_steps = 0

    allocate(res%alpha0_deg(n_stars), res%chi0_deg(n_stars), res%P0_s(n_stars))
    allocate(res%alpha_final_deg(n_stars), res%chi_final_deg(n_stars), res%P_final_s(n_stars))
    allocate(res%delta_alpha_deg(n_stars), res%delta_chi_deg(n_stars), res%delta_P_s(n_stars))
    allocate(res%mean_q_alpha(n_stars), res%mean_mdot(n_stars))
    allocate(res%final_rm_over_rco(n_stars), res%final_omega_s(n_stars))
    allocate(res%final_n_fastness(n_stars), res%final_interaction_factor(n_stars))

    res%alpha0_deg = 0.0_dp
    res%chi0_deg = 0.0_dp
    res%P0_s = 0.0_dp
    res%alpha_final_deg = 0.0_dp
    res%chi_final_deg = 0.0_dp
    res%P_final_s = 0.0_dp
    res%delta_alpha_deg = 0.0_dp
    res%delta_chi_deg = 0.0_dp
    res%delta_P_s = 0.0_dp
    res%mean_q_alpha = 0.0_dp
    res%mean_mdot = 0.0_dp
    res%final_rm_over_rco = 0.0_dp
    res%final_omega_s = 0.0_dp
    res%final_n_fastness = 0.0_dp
    res%final_interaction_factor = 0.0_dp

  end subroutine allocate_population_result


!===============================================================================
! Be/XRB accretion history
!===============================================================================

  subroutine generate_accretion_history(p, t_end_yr, hist)
    implicit none

    type(NSParams), intent(in) :: p
    real(dp), intent(in) :: t_end_yr
    type(AccretionHistory), intent(out) :: hist

    integer :: n, n_orb, n_typeII_left
    real(dp) :: rho_j, rho_m
    real(dp) :: sigma_component, sigma_log
    real(dp) :: beta_mean, phi_mean
    real(dp) :: mu_u(2), u(2), z1, z2
    real(dp) :: logm, zlog, mdot_mean
    real(dp) :: j(3)

    if (p%use_be_outbursts) then
       n_orb = ceiling_to_int(t_end_yr / p%Porb_yr)
    else
       n_orb = ceiling_to_int(t_end_yr / max(p%Porb_yr, 1.0e-10_dp))
    end if
    if (n_orb < 1) n_orb = 1

    hist%n_orbits = n_orb
    hist%Porb_yr = p%Porb_yr
    hist%outburst_duration_yr = p%outburst_duration_yr

    allocate(hist%jvec(3,n_orb))
    allocate(hist%mdot(n_orb))
    allocate(hist%is_typeI(n_orb))
    allocate(hist%is_typeII(n_orb))
    allocate(hist%u_tilt(2,n_orb))

    hist%is_typeI = .false.
    hist%is_typeII = .false.

    if (.not. p%use_be_outbursts) then
       do n = 1, n_orb
          hist%jvec(:,n) = (/0.0_dp, 0.0_dp, 1.0_dp/)
          hist%mdot(n) = p%Mdot0
          hist%is_typeI(n) = .true.
          hist%is_typeII(n) = .false.
          hist%u_tilt(:,n) = 0.0_dp
       end do
       return
    end if

    beta_mean = deg_to_rad(p%beta_mean_deg)

    if (p%random_mean_tilt_azimuth) then
       phi_mean = 2.0_dp*pi*randu()
    else
       phi_mean = deg_to_rad(p%beta_mean_phi_deg)
    end if

    mu_u(1) = beta_mean * cos(phi_mean)
    mu_u(2) = beta_mean * sin(phi_mean)

    ! If beta_rms is the rms amplitude of the 2D tilt, each component has
    ! beta_rms/sqrt(2).
    sigma_component = deg_to_rad(p%beta_rms_deg) / sqrt(2.0_dp)
    rho_j = exp(-1.0_dp / max(p%corr_time_orbits, 1.0e-12_dp))

    call random_normal_pair(z1, z2)
    u(1) = mu_u(1) + sigma_component*z1
    u(2) = mu_u(2) + sigma_component*z2

    if (p%Mdot_typeI_mean > 0.0_dp) then
       mdot_mean = p%Mdot_typeI_mean
    else
       mdot_mean = p%Mdot0
    end if

    sigma_log = p%Mdot_typeI_sigma_log
    rho_m = exp(-1.0_dp / max(p%Mdot_corr_time_orbits, 1.0e-12_dp))

    call random_normal_pair(zlog, z2)
    logm = sigma_log*zlog

    n_typeII_left = 0

    do n = 1, n_orb
       if (n > 1) then
          call random_normal_pair(z1, z2)
          u = mu_u + rho_j*(u - mu_u) + sigma_component*sqrt(max(0.0_dp,1.0_dp-rho_j**2))*(/z1,z2/)

          call random_normal_pair(zlog, z2)
          logm = rho_m*logm + sigma_log*sqrt(max(0.0_dp,1.0_dp-rho_m**2))*zlog
       end if

       hist%u_tilt(:,n) = u

       j = (/u(1), u(2), 1.0_dp/)
       call normalize_vector(j)
       hist%jvec(:,n) = j

       hist%mdot(n) = mdot_mean * exp(logm)

       if (randu() < min(1.0_dp, max(0.0_dp, p%typeI_probability_per_orbit))) then
          hist%is_typeI(n) = .true.
       else
          hist%is_typeI(n) = .false.
       end if

       if (p%include_typeII) then
          if (n_typeII_left > 0) then
             hist%is_typeII(n) = .true.
             hist%is_typeI(n) = .true.
             hist%mdot(n) = hist%mdot(n) * p%typeII_boost
             n_typeII_left = n_typeII_left - 1
          else
             if (randu() < p%typeII_probability_per_orbit) then
                hist%is_typeII(n) = .true.
                hist%is_typeI(n) = .true.
                hist%mdot(n) = hist%mdot(n) * p%typeII_boost
                n_typeII_left = max(0, nint(p%typeII_duration_orbits) - 1)
             end if
          end if
       end if
    end do

  end subroutine generate_accretion_history


  subroutine deallocate_history(hist)
    implicit none
    type(AccretionHistory), intent(inout) :: hist

    if (allocated(hist%jvec)) deallocate(hist%jvec)
    if (allocated(hist%mdot)) deallocate(hist%mdot)
    if (allocated(hist%is_typeI)) deallocate(hist%is_typeI)
    if (allocated(hist%is_typeII)) deallocate(hist%is_typeII)
    if (allocated(hist%u_tilt)) deallocate(hist%u_tilt)

    hist%n_orbits = 0
  end subroutine deallocate_history


  subroutine accretion_input(p, hist, t, d, Mdot, active, orbit_index)
    implicit none

    type(NSParams), intent(in) :: p
    type(AccretionHistory), intent(in) :: hist
    real(dp), intent(in) :: t
    real(dp), intent(out) :: d(3)
    real(dp), intent(out) :: Mdot
    logical, intent(out) :: active
    integer, intent(out) :: orbit_index

    real(dp) :: t_yr_now, phase_yr
    integer :: n

    if (.not. p%use_be_outbursts) then
       d = (/0.0_dp, 0.0_dp, 1.0_dp/)
       Mdot = p%Mdot0
       active = .true.
       orbit_index = 1
       return
    end if

    t_yr_now = t / yr

    n = int(t_yr_now / p%Porb_yr) + 1
    n = max(1, min(hist%n_orbits, n))

    phase_yr = t_yr_now - real(n-1,dp)*p%Porb_yr

    d = hist%jvec(:,n)
    orbit_index = n

    if (hist%is_typeI(n) .and. phase_yr >= 0.0_dp .and. phase_yr < p%outburst_duration_yr) then
       active = .true.
       Mdot = hist%mdot(n)
    else
       active = .false.
       Mdot = p%Mdot_quiescent
    end if
  end subroutine accretion_input


!===============================================================================
! One-star integration
!===============================================================================
  subroutine integrate_one_star(p, hist, star_id, res, qsum_out)
    implicit none

    type(NSParams), intent(in) :: p
    type(AccretionHistory), intent(in) :: hist
    integer, intent(in) :: star_id
    type(PopulationResult), intent(inout) :: res
    real(dp), intent(out) :: qsum_out

    integer :: n
    real(dp) :: y(6)
    real(dp) :: t
    real(dp) :: Omega, chi, Macc
    real(dp) :: s(3)
    real(dp) :: q_alpha
    type(TorqueTerms) :: tt

    call initial_state(p, hist, y)

    t = 0.0_dp
    qsum_out = 0.0_dp

    ! Initial diagnostics.
    Omega = y(1)
    s = y(2:4)
    chi = y(5)
    Macc = y(6)

    call compute_torque_terms(p, hist, t, Omega, s, chi, Macc, tt)

    q_alpha = sin(tt%alpha)**2 * cos(tt%alpha)
    qsum_out = qsum_out + q_alpha

    res%alpha0_deg(star_id) = rad_to_deg(tt%alpha)
    res%chi0_deg(star_id)   = rad_to_deg(fold_chi_to_0_90(chi))
    res%P0_s(star_id)       = 2.0_dp*pi / Omega

    ! Main evolution loop. Nothing is stored at intermediate orbits.
    do n = 1, hist%n_orbits
       call integrate_one_orbit(p, hist, n, t, y)

       Omega = y(1)
       s = y(2:4)
       chi = y(5)
       Macc = y(6)

       call compute_torque_terms(p, hist, t, Omega, s, chi, Macc, tt)

       q_alpha = sin(tt%alpha)**2 * cos(tt%alpha)
       qsum_out = qsum_out + q_alpha
    end do

    ! Final diagnostics.
    Omega = y(1)
    s = y(2:4)
    chi = y(5)
    Macc = y(6)

    call compute_torque_terms(p, hist, t, Omega, s, chi, Macc, tt)

    res%alpha_final_deg(star_id) = rad_to_deg(tt%alpha)
    res%chi_final_deg(star_id)   = rad_to_deg(fold_chi_to_0_90(chi))
    res%P_final_s(star_id)       = 2.0_dp*pi / Omega

    res%delta_alpha_deg(star_id) = res%alpha_final_deg(star_id) - res%alpha0_deg(star_id)
    res%delta_chi_deg(star_id)   = res%chi_final_deg(star_id)   - res%chi0_deg(star_id)
    res%delta_P_s(star_id)       = res%P_final_s(star_id)       - res%P0_s(star_id)

    res%mean_q_alpha(star_id) = qsum_out / real(hist%n_orbits + 1, dp)

    if (hist%n_orbits > 0) then
       res%mean_mdot(star_id) = sum(hist%mdot) / real(hist%n_orbits, dp)
    else
       res%mean_mdot(star_id) = 0.0_dp
    end if

    res%final_rm_over_rco(star_id) = tt%rm / tt%rco
    res%final_omega_s(star_id) = tt%omega_s
    res%final_n_fastness(star_id) = tt%n_fastness
    res%final_interaction_factor(star_id) = tt%interaction_factor

  end subroutine integrate_one_star


  subroutine integrate_one_orbit(p, hist, n, t, y)
    implicit none

    type(NSParams), intent(in) :: p
    type(AccretionHistory), intent(in) :: hist
    integer, intent(in) :: n
    real(dp), intent(inout) :: t
    real(dp), intent(inout) :: y(6)

    integer :: k, nact, nq
    real(dp) :: dt_active, dt_quiet
    real(dp) :: active_time, quiet_time

    if (.not. p%use_be_outbursts) then
       nact = max(1, p%active_substeps)
       dt_active = p%Porb_yr * yr / real(nact, dp)
       do k = 1, nact
          call rk4_step(p, hist, t, y, dt_active)
          t = t + dt_active
       end do
       return
    end if

    active_time = min(p%outburst_duration_yr, p%Porb_yr) * yr
    quiet_time  = max(0.0_dp, p%Porb_yr - p%outburst_duration_yr) * yr

    nact = max(1, p%active_substeps)
    dt_active = active_time / real(nact, dp)

    do k = 1, nact
       call rk4_step(p, hist, t, y, dt_active)
       t = t + dt_active
    end do

    if (quiet_time > 0.0_dp) then
       nq = max(1, p%quiescent_substeps)
       dt_quiet = quiet_time / real(nq, dp)
       do k = 1, nq
          call rk4_step(p, hist, t, y, dt_quiet)
          t = t + dt_quiet
       end do
    end if
  end subroutine integrate_one_orbit


  subroutine initial_state(p, hist, y)
    implicit none

    type(NSParams), intent(in) :: p
    type(AccretionHistory), intent(in) :: hist
    real(dp), intent(out) :: y(6)

    real(dp) :: s0(3)
    real(dp) :: d0(3)

    if (p%use_be_outbursts) then
       d0 = hist%jvec(:,1)
    else
       d0 = (/0.0_dp, 0.0_dp, 1.0_dp/)
    end if

    call initial_spin_vector_from_alpha(p%alpha0_deg, d0, s0)

    y(1) = 2.0_dp*pi / p%P0
    y(2:4) = s0
    y(5) = deg_to_rad(p%chi0_deg)
    y(6) = 0.0_dp

    call physicalize_state(p, y)
  end subroutine initial_state


  subroutine rk4_step(p, hist, t, y, dt)
    implicit none

    type(NSParams), intent(in) :: p
    type(AccretionHistory), intent(in) :: hist
    real(dp), intent(in) :: t, dt
    real(dp), intent(inout) :: y(6)

    real(dp) :: k1(6), k2(6), k3(6), k4(6)
    real(dp) :: yt(6)
    type(TorqueTerms) :: tt

    if (dt <= 0.0_dp) return

    call rhs(p, hist, t, y, k1, tt)

    yt = y + 0.5_dp*dt*k1
    call physicalize_state(p, yt)
    call rhs(p, hist, t + 0.5_dp*dt, yt, k2, tt)

    yt = y + 0.5_dp*dt*k2
    call physicalize_state(p, yt)
    call rhs(p, hist, t + 0.5_dp*dt, yt, k3, tt)

    yt = y + dt*k3
    call physicalize_state(p, yt)
    call rhs(p, hist, t + dt, yt, k4, tt)

    y = y + dt*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
    call physicalize_state(p, y)
  end subroutine rk4_step


  subroutine rhs(p, hist, t, y, dydt, tt)
    implicit none

    type(NSParams), intent(in) :: p
    type(AccretionHistory), intent(in) :: hist
    real(dp), intent(in) :: t
    real(dp), intent(in) :: y(6)
    real(dp), intent(out) :: dydt(6)
    type(TorqueTerms), intent(out) :: tt

    real(dp) :: yy(6)
    real(dp) :: Omega, chi, Macc
    real(dp) :: s(3)
    real(dp) :: dMacc_dt, dI_dt
    real(dp) :: dsdt(3)
    real(dp) :: A, Nacc_chi

    yy = y
    call physicalize_state(p, yy)

    Omega = yy(1)
    s = yy(2:4)
    chi = yy(5)
    Macc = yy(6)

    call compute_torque_terms(p, hist, t, Omega, s, chi, Macc, tt)

    if (tt%outburst_active .and. tt%accretor .and. p%include_mass_growth) then
       dMacc_dt = tt%Mdot
    else
       dMacc_dt = 0.0_dp
    end if

    if (p%include_mass_growth) then
       dI_dt = p%I_per_Msun / Msun * dMacc_dt
    else
       dI_dt = 0.0_dp
    end if

    dydt(1) = (tt%N_spin + tt%Nmag_BA + tt%Kpsr*(1.0_dp + sin(chi)**2) - dI_dt*Omega) / tt%I

    if (p%include_alpha_alignment) then
       dsdt = (tt%N_geom / (tt%I * Omega)) * (tt%d - tt%cos_alpha * s)
    else
       dsdt = 0.0_dp
    end if

    dydt(2:4) = dsdt

    if (p%include_accretion_obliquity) then
       A = modulation_A(p%eta, tt%alpha, chi)
       Nacc_chi = p%eta * A * tt%N_geom &
                  * sin(tt%alpha)**2 * tt%cos_alpha &
                  * sin(chi) * cos(chi)
    else
       Nacc_chi = 0.0_dp
    end if

    dydt(5) = (Nacc_chi + tt%Kpsr*sin(chi)*cos(chi)) / (tt%I * Omega)
    dydt(6) = dMacc_dt
  end subroutine rhs


  subroutine physicalize_state(p, y)
    implicit none
    type(NSParams), intent(in) :: p
    real(dp), intent(inout) :: y(6)

    real(dp) :: s_norm
    real(dp) :: chi_floor, chi_ceiling

    y(1) = max(y(1), 1.0e-30_dp)

    s_norm = sqrt(sum(y(2:4)**2))
    if (s_norm <= 0.0_dp) then
       y(2:4) = (/0.0_dp, 0.0_dp, 1.0_dp/)
    else
       y(2:4) = y(2:4) / s_norm
    end if

    chi_floor = deg_to_rad(p%chi_floor_deg)
    chi_ceiling = deg_to_rad(p%chi_ceiling_deg)
    y(5) = min(max(y(5), chi_floor), chi_ceiling)

    if (.not. p%include_mass_growth) y(6) = 0.0_dp
    if (p%include_mass_growth) y(6) = max(y(6), 0.0_dp)
  end subroutine physicalize_state


!===============================================================================
! Torque terms
!===============================================================================
  subroutine compute_torque_terms(p, hist, t, Omega, s, chi, Macc, tt)
    implicit none

    type(NSParams), intent(in) :: p
    type(AccretionHistory), intent(in) :: hist
    real(dp), intent(in) :: t, Omega, s(3), chi, Macc
    type(TorqueTerms), intent(out) :: tt

    logical :: active
    integer :: orbit_index

    call accretion_input(p, hist, t, tt%d, tt%Mdot, active, orbit_index)

    tt%outburst_active = active
    tt%orbit_index = orbit_index

    tt%cos_alpha = max(-1.0_dp, min(1.0_dp, dot_product(s, tt%d)))
    tt%alpha = acos(tt%cos_alpha)

    if (p%include_mass_growth) then
       tt%M = p%M0 + Macc
    else
       tt%M = p%M0
    end if

    tt%I = inertia(tt%M, p)
    tt%mu = compute_mu(Macc, max(tt%Mdot, 1.0e-300_dp), p)

    if (tt%Mdot > 0.0_dp) then
       tt%rm = magnetosphere_radius(tt%mu, tt%M, tt%Mdot, p%xi)
       tt%N0 = tt%Mdot * sqrt(Gcgs * tt%M * tt%rm)
    else
       tt%rm = 1.0e99_dp
       tt%N0 = 0.0_dp
    end if

    tt%rco = corotation_radius(tt%M, Omega)
    tt%rlc = light_cylinder_radius(Omega)

    tt%omega_s = fastness_parameter(tt%rm, tt%rco)

    tt%accretor = (tt%Mdot > 0.0_dp .and. tt%rm < tt%rco)
    tt%inside_lc = (tt%rm < tt%rlc)

    tt%n_fastness = torque_model_factor(p, tt%omega_s)
    tt%interaction_factor = interaction_factor_model(p, tt%omega_s, tt%accretor)

    if (p%torque_model == TORQUE_BA2021) then
       if (tt%accretor) then
          tt%N_spin = tt%N0 * tt%cos_alpha
       else
          tt%N_spin = 0.0_dp
       end if
    else
       if (tt%Mdot > 0.0_dp) then
          tt%N_spin = tt%N0 * tt%n_fastness * tt%cos_alpha
       else
          tt%N_spin = 0.0_dp
       end if
    end if

    tt%Nmag_BA = 0.0_dp
    if (p%include_BA_magnetospheric_braking .and. tt%inside_lc) then
       tt%Nmag_BA = -tt%mu**2 / (3.0_dp * tt%rco**3)
    end if

    tt%Kpsr = 0.0_dp
    if (p%include_pulsar_torque) then
       tt%Kpsr = -tt%mu**2 / tt%rlc**3
       if (p%include_parfrey_enhancement .and. tt%accretor) then
          tt%Kpsr = tt%Kpsr * (tt%rlc / max(tt%rm, 1.0e-300_dp))**2
       end if
    end if

    tt%N_geom = tt%N0 * tt%interaction_factor
  end subroutine compute_torque_terms


  real(dp) function torque_model_factor(p, omega_s) result(n)
    implicit none
    type(NSParams), intent(in) :: p
    real(dp), intent(in) :: omega_s
    real(dp) :: x

    x = omega_s / p%omega_c

    select case (p%torque_model)
    case (TORQUE_BA2021)
       n = 1.0_dp

    case (TORQUE_LINEAR_FASTNESS)
       n = 1.0_dp - x

    case (TORQUE_SATURATED_FASTNESS)
       n = (1.0_dp - x) / (1.0_dp + abs(x))

    case (TORQUE_TANH_FASTNESS)
       n = tanh((1.0_dp - x) / max(p%tanh_width, 1.0e-12_dp))

    case default
       write(*,*) '# ERROR: Unknown torque_model = ', p%torque_model
       stop
    end select
  end function torque_model_factor


  real(dp) function interaction_factor_model(p, omega_s, accretor) result(f)
    implicit none
    type(NSParams), intent(in) :: p
    real(dp), intent(in) :: omega_s
    logical, intent(in) :: accretor
    real(dp) :: x

    select case (p%obliquity_torque_mode)
    case (OBL_ACCRETOR_ONLY)
       if (accretor) then
          f = 1.0_dp
       else
          f = 0.0_dp
       end if

    case (OBL_SAME_AS_SPIN)
       if (p%torque_model == TORQUE_BA2021) then
          if (accretor) then
             f = 1.0_dp
          else
             f = 0.0_dp
          end if
       else
          f = torque_model_factor(p, omega_s)
       end if

    case (OBL_POSITIVE_INTERACTION)
       x = omega_s / p%omega_c
       if (x <= 1.0_dp) then
          f = 1.0_dp
       else
          f = 1.0_dp / (1.0_dp + (x - 1.0_dp)**p%interaction_decay_power)
       end if

    case default
       write(*,*) '# ERROR: Unknown obliquity_torque_mode = ', p%obliquity_torque_mode
       stop
    end select
  end function interaction_factor_model


!===============================================================================
! Diagnostics and storage
!===============================================================================
  subroutine store_diagnostics(p, hist, t, y, star_id, k, res, q_alpha, tt)
    implicit none

    type(NSParams), intent(in) :: p
    type(AccretionHistory), intent(in) :: hist
    real(dp), intent(in) :: t, y(6)
    integer, intent(in) :: star_id, k
    type(PopulationResult), intent(inout) :: res
    real(dp), intent(out) :: q_alpha
    type(TorqueTerms), intent(out) :: tt

    real(dp) :: yy(6)
    real(dp) :: Omega, chi, Macc
    real(dp) :: s(3)

    yy = y
    call physicalize_state(p, yy)

    Omega = yy(1)
    s = yy(2:4)
    chi = yy(5)
    Macc = yy(6)

    call compute_torque_terms(p, hist, t, Omega, s, chi, Macc, tt)

    q_alpha = sin(tt%alpha)**2 * cos(tt%alpha)

    ! Avoid accidental use in final-only mode.
    if (.not. allocated(res%P_track)) then
       write(*,*) '# ERROR: store_diagnostics called, but full tracks are not allocated.'
       stop
    end if

    ! Avoid a harmless but formal OpenMP race: only one star writes the common
    ! time grid.
    if (star_id == 1) res%t_yr(k) = t / yr

    res%P_track(k, star_id) = 2.0_dp*pi / Omega
    res%alpha_track(k, star_id) = rad_to_deg(tt%alpha)
    res%chi_track(k, star_id) = rad_to_deg(fold_chi_to_0_90(chi))
    res%rm_over_rco_track(k, star_id) = tt%rm / tt%rco
    res%omega_s_track(k, star_id) = tt%omega_s
    res%q_alpha_track(k, star_id) = q_alpha
    res%mdot_track(k, star_id) = tt%Mdot
  end subroutine store_diagnostics


!===============================================================================
! Geometry
!===============================================================================
  subroutine initial_spin_vector_from_alpha(alpha0_deg, d0, s0)
    implicit none
    real(dp), intent(in) :: alpha0_deg
    real(dp), intent(in) :: d0(3)
    real(dp), intent(out) :: s0(3)

    real(dp) :: z(3), e(3)
    real(dp) :: alpha0

    z = (/0.0_dp, 0.0_dp, 1.0_dp/)

    e = z - dot_product(z, d0) * d0
    if (sqrt(sum(e**2)) < 1.0e-12_dp) then
       e = (/1.0_dp, 0.0_dp, 0.0_dp/)
    end if
    call normalize_vector(e)

    alpha0 = deg_to_rad(alpha0_deg)
    s0 = cos(alpha0)*d0 + sin(alpha0)*e
    call normalize_vector(s0)
  end subroutine initial_spin_vector_from_alpha


  subroutine normalize_vector(v)
    implicit none
    real(dp), intent(inout) :: v(3)
    real(dp) :: n

    n = sqrt(sum(v**2))
    if (n <= 0.0_dp) then
       v = (/0.0_dp, 0.0_dp, 1.0_dp/)
    else
       v = v / n
    end if
  end subroutine normalize_vector


  real(dp) function angle_between(u, v) result(a)
    implicit none
    real(dp), intent(in) :: u(3), v(3)
    real(dp) :: uu(3), vv(3), dotuv

    uu = u
    vv = v
    call normalize_vector(uu)
    call normalize_vector(vv)

    dotuv = max(-1.0_dp, min(1.0_dp, dot_product(uu, vv)))
    a = acos(dotuv)
  end function angle_between


!===============================================================================
! Physical functions
!===============================================================================

  real(dp) function inertia(M, p) result(I)
    implicit none
    real(dp), intent(in) :: M
    type(NSParams), intent(in) :: p
    I = p%I_per_Msun * (M / Msun)
  end function inertia


  real(dp) function magnetosphere_radius(mu, M, Mdot, xi) result(rm)
    implicit none
    real(dp), intent(in) :: mu, M, Mdot, xi
    rm = xi * (mu**4 / (2.0_dp * Gcgs * M * max(Mdot,1.0e-300_dp)**2))**(1.0_dp/7.0_dp)
  end function magnetosphere_radius


  real(dp) function corotation_radius(M, Omega) result(rco)
    implicit none
    real(dp), intent(in) :: M, Omega
    rco = (Gcgs * M / Omega**2)**(1.0_dp/3.0_dp)
  end function corotation_radius


  real(dp) function light_cylinder_radius(Omega) result(rlc)
    implicit none
    real(dp), intent(in) :: Omega
    rlc = ccgs / Omega
  end function light_cylinder_radius


  real(dp) function fastness_parameter(rm, rco) result(omega_s)
    implicit none
    real(dp), intent(in) :: rm, rco
    omega_s = (rm / rco)**1.5_dp
  end function fastness_parameter


  real(dp) function modulation_A(eta, alpha, chi) result(A)
    implicit none
    real(dp), intent(in) :: eta, alpha, chi
    real(dp) :: denom

    denom = 1.0_dp - 0.5_dp*eta * &
       (sin(chi)**2 * sin(alpha)**2 + 2.0_dp*cos(chi)**2 * cos(alpha)**2)

    A = 1.0_dp / max(denom, 1.0e-14_dp)
  end function modulation_A


  real(dp) function compute_mu(Macc, Mdot, p) result(mu)
    implicit none
    real(dp), intent(in) :: Macc, Mdot
    type(NSParams), intent(in) :: p

    real(dp) :: mu0, Mdot1, mu0_30, R_12p5, Mscale

    if (p%B_convention == 1) then
       mu0 = p%B0 * p%R**3
    else
       mu0 = 0.5_dp * p%B0 * p%R**3
    end if

    mu = mu0

    if (p%include_magnetic_burial) then
       Mdot1 = max(Mdot / Mdot_Edd_solar, 1.0e-30_dp)
       mu0_30 = mu0 / 1.0e30_dp
       R_12p5 = p%R / 1.25e6_dp
       Mscale = 1.1e-5_dp * Mdot1**(1.0_dp/7.0_dp) &
                * mu0_30**(3.0_dp/14.0_dp) * R_12p5**3 * Msun
       mu = mu0 * (1.0_dp + max(Macc,0.0_dp)/Mscale)**(-14.0_dp/11.0_dp)
    end if
  end function compute_mu


  real(dp) function mdot_from_lx(Lx, M, R, efficiency_factor) result(Mdot)
    implicit none
    real(dp), intent(in) :: Lx, M, R, efficiency_factor
    Mdot = Lx * R / (efficiency_factor * Gcgs * M)
  end function mdot_from_lx


  real(dp) function lx_from_mdot(Mdot, M, R, efficiency_factor) result(Lx)
    implicit none
    real(dp), intent(in) :: Mdot, M, R, efficiency_factor
    Lx = efficiency_factor * Gcgs * M * Mdot / R
  end function lx_from_mdot


!===============================================================================
! Random numbers
!===============================================================================

  subroutine init_random_seed(seed)
    implicit none
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: seed_array(:)

    call random_seed(size=n)
    allocate(seed_array(n))

    do i = 1, n
       seed_array(i) = seed + 37*i + 101*i*i
    end do

    call random_seed(put=seed_array)
    deallocate(seed_array)
  end subroutine init_random_seed


  real(dp) function randu() result(x)
    implicit none
    call random_number(x)
  end function randu


  subroutine random_normal_pair(z1, z2)
    implicit none
    real(dp), intent(out) :: z1, z2
    real(dp) :: u1, u2

    u1 = max(randu(), 1.0e-300_dp)
    u2 = randu()

    z1 = sqrt(-2.0_dp*log(u1)) * cos(2.0_dp*pi*u2)
    z2 = sqrt(-2.0_dp*log(u1)) * sin(2.0_dp*pi*u2)
  end subroutine random_normal_pair


  real(dp) function draw_alpha_isotropic_deg(alpha_min_deg, alpha_max_deg) result(alpha_deg)
    implicit none
    real(dp), intent(in) :: alpha_min_deg, alpha_max_deg

    real(dp) :: amin, amax, cmin, cmax, c

    amin = deg_to_rad(alpha_min_deg)
    amax = deg_to_rad(alpha_max_deg)

    cmin = cos(amax)
    cmax = cos(amin)

    c = cmin + (cmax - cmin) * randu()
    c = max(-1.0_dp, min(1.0_dp, c))

    alpha_deg = rad_to_deg(acos(c))
  end function draw_alpha_isotropic_deg


  real(dp) function draw_chi_isotropic_folded_deg(chi_min_deg, chi_max_deg) result(chi_deg)
    implicit none
    real(dp), intent(in) :: chi_min_deg, chi_max_deg

    real(dp) :: cmin, cmax, c

    cmin = cos(deg_to_rad(chi_max_deg))
    cmax = cos(deg_to_rad(chi_min_deg))

    c = cmin + (cmax - cmin) * randu()
    c = max(-1.0_dp, min(1.0_dp, c))

    chi_deg = rad_to_deg(acos(c))
  end function draw_chi_isotropic_folded_deg


!===============================================================================
! Printing/output
!===============================================================================

  subroutine print_configuration(p, n_stars, t_end_yr, n_steps)
    implicit none

    type(NSParams), intent(in) :: p
    integer, intent(in) :: n_stars, n_steps
    real(dp), intent(in) :: t_end_yr

    write(*,*)
    write(*,*) '#=============================================================='
    write(*,*) '# Be/XRB NS population run configuration'
    write(*,*) '# Final-only mode: full tracks are not stored'
    write(*,*) '#=============================================================='
    write(*,'(A,I8)')      '# N stars                   = ', n_stars
    write(*,'(A,I8)')      '# N stored orbit steps       = ', n_steps
    write(*,'(A,ES12.4)')  '# t_end [yr]                = ', t_end_yr
    write(*,'(A,ES12.4)')  '# P0 [s]                    = ', p%P0
    write(*,'(A,ES12.4)')  '# B0 [G]                    = ', p%B0
    write(*,'(A,ES12.4)')  '# Mdot0/Type-I mean [g/s]   = ', p%Mdot0
    write(*,'(A,ES12.4)')  '# eta                       = ', p%eta
    write(*,'(A,I8)')      '# torque_model              = ', p%torque_model
    write(*,'(A,I8)')      '# obliquity_torque_mode     = ', p%obliquity_torque_mode
    write(*,'(A,L8)')      '# use_be_outbursts          = ', p%use_be_outbursts
    write(*,'(A,ES12.4)')  '# Porb [yr]                 = ', p%Porb_yr
    write(*,'(A,ES12.4)')  '# outburst duration [yr]    = ', p%outburst_duration_yr
    write(*,'(A,ES12.4)')  '# Type-I probability/orbit  = ', p%typeI_probability_per_orbit
    write(*,'(A,ES12.4)')  '# beta_mean [deg]           = ', p%beta_mean_deg
    write(*,'(A,ES12.4)')  '# beta_rms [deg]            = ', p%beta_rms_deg
    write(*,'(A,ES12.4)')  '# corr_time [orbits]        = ', p%corr_time_orbits
    write(*,'(A,ES12.4)')  '# Mdot sigma log            = ', p%Mdot_typeI_sigma_log
    write(*,'(A,L8)')      '# include Type II           = ', p%include_typeII
    write(*,*) '#=============================================================='
    write(*,*)
  end subroutine print_configuration


  subroutine print_population_table(res)
    implicit none

    type(PopulationResult), intent(in) :: res
    integer :: i

    write(*,*)
    write(*,*) '# Population table: initial and final parameters'
    write(*,*) '#---------------------------------------------------------------------------------------------'
    write(*,'(A)') '# id   alpha0   alphaf     chi0     chif       P0          Pf       rm/rco_f   omega_f   <Qalpha>'
    write(*,*) '#---------------------------------------------------------------------------------------------'

    do i = 1, res%n_stars
       write(*,'(I3,2X,F8.3,2X,F8.3,2X,F8.3,2X,F8.3,2X,ES10.3,2X,ES10.3,2X,F9.4,2X,F9.4,2X,F9.4)') &
          i, res%alpha0_deg(i), res%alpha_final_deg(i), &
          res%chi0_deg(i), res%chi_final_deg(i), &
          res%P0_s(i), res%P_final_s(i), &
          res%final_rm_over_rco(i), res%final_omega_s(i), res%mean_q_alpha(i)
    end do

    write(*,*) '#---------------------------------------------------------------------------------------------'
    write(*,*)
  end subroutine print_population_table


  subroutine print_spin_statistics(res)
  implicit none
  type(PopulationResult), intent(in) :: res
  integer :: i, n
  real(dp) :: mean_P, sigma_P
  real(dp) :: mean_logP, sigma_logP
  real(dp) :: x, dx
  real(dp) :: Pmin, Pmax
    n = res%n_stars
    if (n <= 0) then
      write(*,*) '# WARNING: no stars in population; spin statistics skipped.'
      return
    end if
    mean_P = 0.0_dp
    mean_logP = 0.0_dp

    Pmin = huge(1.0_dp)
    Pmax = -huge(1.0_dp)
    do i = 1, n
      x = res%P_final_s(i)
      mean_P = mean_P + x
      mean_logP = mean_logP + log10(max(x, 1.0e-300_dp))
      Pmin = min(Pmin, x)
      Pmax = max(Pmax, x)
    end do
    mean_P = mean_P / real(n, dp)
    mean_logP = mean_logP / real(n, dp)
    sigma_P = 0.0_dp
    sigma_logP = 0.0_dp
    if (n > 1) then
      do i = 1, n
        x = res%P_final_s(i)
        dx = x - mean_P
        sigma_P = sigma_P + dx*dx
        dx = log10(max(x, 1.0e-300_dp)) - mean_logP
        sigma_logP = sigma_logP + dx*dx
      end do
      sigma_P = sqrt(sigma_P / real(n - 1, dp))
      sigma_logP = sqrt(sigma_logP / real(n - 1, dp))
    else
      sigma_P = 0.0_dp
      sigma_logP = 0.0_dp
    end if
   write(*,*)
   write(*,*) '# Spin-period statistics at final time'
   write(*,*) '#--------------------------------------------------------------'
   write(*,'(A,ES14.6)') '# <P_final> [s]              = ', mean_P
   write(*,'(A,ES14.6)') '# sigma(P_final) [s]         = ', sigma_P
   write(*,'(A,ES14.6)') '# min(P_final) [s]           = ', Pmin
   write(*,'(A,ES14.6)') '# max(P_final) [s]           = ', Pmax
   write(*,'(A,F10.5)')  '# <log10 P_final>            = ', mean_logP
   write(*,'(A,F10.5)')  '# sigma(log10 P_final) [dex] = ', sigma_logP
   write(*,*) '#--------------------------------------------------------------'
   write(*,*)
end subroutine print_spin_statistics


  subroutine print_histograms(res, nbins_alpha, nbins_chi)
    implicit none

    type(PopulationResult), intent(in) :: res
    integer, intent(in) :: nbins_alpha, nbins_chi

    integer, allocatable :: h_alpha0(:), h_alphaf(:), h_chi0(:), h_chif(:)
    integer :: i
    real(dp) :: lo, hi

    allocate(h_alpha0(nbins_alpha), h_alphaf(nbins_alpha))
    allocate(h_chi0(nbins_chi), h_chif(nbins_chi))

    call make_histogram(res%alpha0_deg, res%n_stars, 0.0_dp, 180.0_dp, nbins_alpha, h_alpha0)
    call make_histogram(res%alpha_final_deg, res%n_stars, 0.0_dp, 180.0_dp, nbins_alpha, h_alphaf)
    call make_histogram(res%chi0_deg, res%n_stars, 0.0_dp, 90.0_dp, nbins_chi, h_chi0)
    call make_histogram(res%chi_final_deg, res%n_stars, 0.0_dp, 90.0_dp, nbins_chi, h_chif)

    write(*,*)
    write(*,*) '# Histogram: alpha initial/final'
    write(*,*) '#---------------------------------------------'
    write(*,'(A)') '# bin_low  bin_high   N_initial   N_final'
    do i = 1, nbins_alpha
       lo = 180.0_dp * real(i-1, dp) / real(nbins_alpha, dp)
       hi = 180.0_dp * real(i, dp) / real(nbins_alpha, dp)
       write(*,'(F8.2,2X,F8.2,2X,I8,2X,I8)') lo, hi, h_alpha0(i), h_alphaf(i)
    end do

    write(*,*)
    write(*,*) '# Histogram: chi initial/final'
    write(*,*) '#---------------------------------------------'
    write(*,'(A)') '# bin_low  bin_high   N_initial   N_final'
    do i = 1, nbins_chi
       lo = 90.0_dp * real(i-1, dp) / real(nbins_chi, dp)
       hi = 90.0_dp * real(i, dp) / real(nbins_chi, dp)
       write(*,'(F8.2,2X,F8.2,2X,I8,2X,I8)') lo, hi, h_chi0(i), h_chif(i)
    end do
    write(*,*)

    deallocate(h_alpha0, h_alphaf, h_chi0, h_chif)
  end subroutine print_histograms


  subroutine make_histogram(x, n, xmin, xmax, nbins, h)
    implicit none

    integer, intent(in) :: n, nbins
    real(dp), intent(in) :: x(n), xmin, xmax
    integer, intent(out) :: h(nbins)

    integer :: i, b
    real(dp) :: dx

    h = 0
    dx = (xmax - xmin) / real(nbins, dp)

    do i = 1, n
       if (x(i) < xmin .or. x(i) > xmax) cycle

       if (x(i) == xmax) then
          b = nbins
       else
          b = int((x(i) - xmin) / dx) + 1
       end if

       if (b >= 1 .and. b <= nbins) h(b) = h(b) + 1
    end do
  end subroutine make_histogram


  subroutine write_population_ascii(filename, res)
    implicit none

    character(len=*), intent(in) :: filename
    type(PopulationResult), intent(in) :: res
    integer :: i, u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,'(A)') '# id alpha0 alphaf chi0 chif P0 Pf rm_over_rco_f omega_f mean_q_alpha mean_mdot'
    do i = 1, res%n_stars
       write(u,'(I8,1X,10(ES20.10,1X))') i, &
          res%alpha0_deg(i), res%alpha_final_deg(i), &
          res%chi0_deg(i), res%chi_final_deg(i), &
          res%P0_s(i), res%P_final_s(i), &
          res%final_rm_over_rco(i), res%final_omega_s(i), &
          res%mean_q_alpha(i), res%mean_mdot(i)
    end do
    close(u)
  end subroutine write_population_ascii


  subroutine write_tracks_ascii(filename, res)
    implicit none

    character(len=*), intent(in) :: filename
    type(PopulationResult), intent(in) :: res
    integer :: i, k, u

    if (.not. allocated(res%P_track)) then
       write(*,*) '# ERROR: write_tracks_ascii called, but full tracks are not allocated in final-only mode.'
       stop
    end if

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,'(A)') '# star_id t_yr P alpha chi rm_over_rco omega_s q_alpha mdot'
    do i = 1, res%n_stars
       do k = 0, res%n_steps
          write(u,'(I8,1X,8(ES20.10,1X))') i, res%t_yr(k), &
             res%P_track(k,i), res%alpha_track(k,i), res%chi_track(k,i), &
             res%rm_over_rco_track(k,i), res%omega_s_track(k,i), &
             res%q_alpha_track(k,i), res%mdot_track(k,i)
       end do
    end do
    close(u)
  end subroutine write_tracks_ascii


!===============================================================================
! Math/helpers
!===============================================================================

  integer function ceiling_to_int(x) result(n)
    implicit none
    real(dp), intent(in) :: x
    n = int(x)
    if (real(n,dp) < x) n = n + 1
  end function ceiling_to_int


  real(dp) function deg_to_rad(x) result(y)
    implicit none
    real(dp), intent(in) :: x
    y = x*pi/180.0_dp
  end function deg_to_rad


  real(dp) function rad_to_deg(x) result(y)
    implicit none
    real(dp), intent(in) :: x
    y = x*180.0_dp/pi
  end function rad_to_deg


  real(dp) function fold_chi_to_0_90(chi) result(chif)
    implicit none
    real(dp), intent(in) :: chi
    real(dp) :: x

    x = modulo(chi, pi)
    chif = min(x, pi - x)
  end function fold_chi_to_0_90

end module ns_population_torque_models

