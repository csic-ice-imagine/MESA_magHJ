!-----------------------------------------------------------------------
! Modifiable file.
!
! The files included below are located in the MESA folder /include/.
!-----------------------------------------------------------------------
module run_star_extras

   use star_lib
   use star_def
   use const_def
   use math_lib

   implicit none

   include "test_suite_extras_def.inc"

   contains

      include "test_suite_extras.inc"

!-----------------------------------------------------------------------
! Point to the new heat source subroutine.
!-----------------------------------------------------------------------
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         ! Point to the new heat source subroutine (also in this file).
         s% other_energy => default_other_energy
      end subroutine extras_controls

      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup

      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve

      ! Returns either keep_going, retry or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
      end function extras_check_model

      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns

      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns

      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns

      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns

      ! Returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step

!-----------------------------------------------------------------------
! Joule heating.
! Needed variables: R, M and L (for B) and conductivity profile (sigma).
!
! Notes:
! - Follow the instructions provided in the README file located in the
! MESA folder /star/other/.
! - The module const_def is defined in /const/public/const_def.f90 and
! contains standard definitions of constants.
! - The derived type star_info (of the variable s) is defined in the
! file /star_data/public/star_data_def.f90 (as well as the subroutine
! star_ptr). Most declarations of the derived type are in the included
! file star_data.inc in the same folder.
! - Surface irradiation from the host star is controlled by the input
! variables irradiation_flux and column_depth_for_irradiation (defined
! in inlist_evolution). The input variables are processed by the file
! /star/private/ctrls_io.f90 where the "controls" namelist is defined.
!-----------------------------------------------------------------------
      subroutine default_other_energy(id, ierr)
         ! use const_def, only: Rsun ! Solar radius (in cm).
         use const_def, only: Lsun,Msun,Rsun,m_jupiter,r_jupiter, &
                              m_earth,r_earth,dp,pi
         ! Note: The "only" attribute may be more relevant to naming conflicts.

         ! Default definitions.
         implicit none
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ! Custom (local) definitions.
         integer :: k,ndim,ntot
         real(dp) :: age,density,pressure,temperature,entropy,gravity
         real(dp) :: radius,rr,mass,dmass,luminosity
         real(dp) :: radius_p,mass_p,luminosity_p ! Planetary quantities (in cgs).
         real(dp) :: radius_j,mass_j ! Jupiter radius and mass (in cgs).
         real(dp) :: B_dynamo,B_dipole
         real(dp) :: c_cgs ! Speed of light.

         real(dp) :: eps_irrad,eps_joule,eps_joule1,eps_joule2
         real(dp) :: lum_irrad,lum_joule,lum_joule1,lum_joule2
         real(dp) :: lum_joule2_above,lum_joule2_below
         real(dp) :: lum_joule_above,lum_joule_below ! Old, remove!

         real(dp) :: irradiation,Teff,Teq
         real(dp) :: log_sigma,sigma2_amp,log_v,v_atm,v_amp
         real(dp) :: sigma1_log,log_sigma_eye,sigma1_eye,log_sigma_tsang,sigma1_tsang
         real(dp) :: length_scale_input
         real(dp) :: current1,joule_heating1,length_scale1,sigma1
         real(dp) :: current2,joule_heating2,length_scale2,sigma2
         real(dp) :: Bphi2_l,Bphi2_n,J2_l,J2_n,Q2_l,Q2_n,exp_fn ! Linear & non-linear cases.
         real(dp) :: c1,c2,c3,c4,c5,aa,bb,cc,rj_tj20,mu_0 ! Analytic fit for sigma.
         real(dp) :: inertia_moment,radius_1

         ! External table for sigma1.
         integer :: k_tab
         integer :: ndim_tab = 501
         real(dp), dimension(501) :: p_tab,r_tab,sigma_tab
         real(dp) :: dummy,sigma1_int
         character(1) :: dummy_string

         ! Needed for the radiative-convective boundary (RCB).
         integer :: n1_rcb,n3_rcb,n4_rcb
         real(dp) :: r1_old,r3_old,r4_old,r1_rcb,r3_rcb,r4_rcb ! Radii.
         real(dp) :: p1_old,p3_old,p4_old,p1_rcb,p3_rcb,p4_rcb ! Pressures.
         real(dp) :: p_surface,p12_boundary,delta_rcb

         ! Internally used call counters.
         ! Variables with the SAVE attribute are initialized in the declaration.
         integer, save :: n_call = 0
         logical, save :: first_call = .TRUE.
         logical, save :: open_file_1 = .TRUE.
         logical, save :: write_file_1 = .FALSE.
         logical, save :: open_file_2 = .TRUE.
         logical, save :: write_file_2 = .FALSE.
         logical, save :: open_file_3 = .TRUE.
         logical, save :: write_file_3 = .FALSE.

         ! Speed of light (in cgs).
         c_cgs = 2.99792458d10

         ! Jupiter mass and radius (in cgs).
         mass_j = m_jupiter
         radius_j = r_jupiter

         ! Print to screen.
         !write(*,*)"Defined constants:"
         !write(*,*)Lsun,Msun,Rsun
         !write(*,*)m_jupiter,r_jupiter,m_earth,r_earth
         !write(*,*)dp,pi
         !stop

         !--------------------------------------------------------------
         ! Initialize the pointer s (of derived type star_info).
         ! Before the call to star_ptr the pointer s is uninitialized.
         ! Note: Assigning uninitialized pointers to variables can cause
         ! segmentation fault.
         !--------------------------------------------------------------
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! After the call to star_ptr the pointer s is initialized.
         ndim = s% nz ! Number of radial zones. (Variable!)
         ntot = size(s% r(:)) ! Total radial array size. (Fixed!)

         age = s% star_age ! Current age [in yrs].

         Teff = s% Teff ! Effective temperature of the planet.

         ! Planetary luminosity, mass and radius (in cgs).
         ! The surface values are at index 1.
         luminosity_p = s% L(1)
         mass_p = s% m(1)
         radius_p = s% r(1)

         !--------------------------------------------------------------
         ! Load external table for sigma1.
         ! (Should be done only once in a run. To be cleaned up!)
         !--------------------------------------------------------------
         open(21,file = "in/conductivity_table.txt")
         read(21,*)dummy_string
         read(21,*)dummy_string
         read(21,*)dummy_string
         read(21,*)dummy_string
         ! There are multiple profiles of sigma in the input file.
         ! Decide which one to read here.
         do k_tab = 1,ndim_tab
            read(21,*)dummy,r_tab(k_tab),p_tab(k_tab),sigma_tab(k_tab)
         enddo
         close(21)

!         ! Check the input file.
!         open(22,file = "in/check.txt")
!         write(22,*)"# File output as a check. (See the source code for details.)"
!         do k = 1,ndim
!            pressure = s% P(k)
!            !-----------------------------------------------------------
!            ! Conductivity from an external table.
!            !-----------------------------------------------------------
!            sigma1_int = sigma_tab(ndim_tab)
!            if(p_tab(1).gt.pressure)then
!               sigma1_int = sigma_tab(1)
!               goto 10
!            endif
!            do k_tab = 2,ndim_tab
!               if(p_tab(k_tab).gt.pressure)then
!                  sigma1_int = sigma_tab(k_tab) - (sigma_tab(k_tab) - sigma_tab(k_tab-1))* &
!                     (p_tab(k_tab) - pressure) / (p_tab(k_tab) - p_tab(k_tab-1))
!                  goto 10
!               endif
!            enddo
!10          continue
!            write(22,*)pressure,sigma1_int
!         enddo
!         close(22)
!         stop

         !--------------------------------------------------------------
         ! Read the input parameter length_scale_input.
         !--------------------------------------------------------------
         open(22,file = "in/additional_input.txt")
         read(22,*)dummy_string
         read(22,*)length_scale_input
         close(22)

         !--------------------------------------------------------------
         ! Magnetic field strength from scaling laws.
         !--------------------------------------------------------------
         ! Dynamo magnetic field [in Gauss].
         ! (Reiners & Christensen 2010, eq.1.)
         B_dynamo = 4.8d3*((mass_p/Msun)*(luminosity_p/Lsun)**2d0/(radius_p/Rsun)**7d0)**(1d0/6d0)
         ! Dipole magnetic field strength at the equator [in Gauss].
         ! (Reiners & Christensen 2010, eq.2.)
         B_dipole = (B_dynamo/(2d0*sqrt(2d0)))*(1d0 - 0.17d0*mass_j/mass_p)**3d0

         !--------------------------------------------------------------
         ! Surface irradiation (dayside flux from the input).
         !--------------------------------------------------------------
         irradiation = s% irradiation_flux
         Teq = (irradiation/4d0/5.67d-5)**0.25d0

         lum_irrad = 0d0
         lum_joule = 0d0
         lum_joule_below = 0d0 ! Below the RCB.
         lum_joule_above = 0d0 ! Above the RCB.

         !--------------------------------------------------------------
         ! Determine the pressure at the RCB (linear interpolation).
         !--------------------------------------------------------------
         ! Define the surface pressure (at outermost cell).
         p_surface = s% P(1)
         ! Initial values.
         r1_old = radius_p
         r3_old = radius_p
         r4_old = radius_p
         p1_old = p_surface
         p3_old = p_surface
         p4_old = p_surface
         r1_rcb = (s% conv_mx1_top_r) * Rsun
         r3_rcb = (s% conv_mx2_top_r) * Rsun
         r4_rcb = (s% conv_mx2_bot_r) * Rsun
         p1_rcb = 0d0
         p3_rcb = 0d0
         p4_rcb = 0d0
         n1_rcb = 0
         n3_rcb = 0
         n4_rcb = 0

         ! Loop over radial index. (The radial dimension ndim is variable.)
         ! From outside inwards (k=1 corresponds to the surface).
         do k = 1,ndim
            radius = s% r(k)
            pressure = s% P(k) ! Total pressure (P = P_rad + P_gas).

            ! Main convection region (extending from r1_rcb down to the center).
            if(r1_old.gt.r1_rcb.and.radius.le.r1_rcb)then
               p1_rcb = p1_old + (pressure - p1_old) * (r1_rcb - r1_old) / (radius - r1_old)
               n1_rcb = n1_rcb + 1 ! Should always be 1.
               !write(*,*)"RCB",r1_rcb,r1_old,radius
               !write(*,*)"RCB",p1_rcb,p1_old,pressure
            endif
            r1_old = radius
            p1_old = pressure

            ! Secondary convection region (from r3_rcb down to r4_rcb).
            if(r3_rcb.ne.0d0.and.r3_old.gt.r3_rcb.and.radius.le.r3_rcb)then
               p3_rcb = p3_old + (pressure - p3_old) * (r3_rcb - r3_old) / (radius - r3_old)
               n3_rcb = n3_rcb + 1 ! Should always be 1 unless it is zero.
            endif
            r3_old = radius
            p3_old = pressure

            if(r4_rcb.ne.0d0.and.r4_old.gt.r4_rcb.and.radius.le.r4_rcb)then
               p4_rcb = p4_old + (pressure - p4_old) * (r4_rcb - r4_old) / (radius - r4_old)
               n4_rcb = n4_rcb + 1 ! Should always be 1 unless it is zero.
            endif
            r4_old = radius
            p4_old = pressure
         enddo

         !--------------------------------------------------------------
         ! Boundary between the two Ohmic heating terms (J1 and J2).
         !--------------------------------------------------------------
         delta_rcb = 0d0
         ! p12_boundary = p1_rcb*(1d0 + delta_rcb)
         ! p12_boundary = 8d8 ! Old Ohmic heating 1 range.
         p12_boundary = 1d9 ! Cut-off at 10^3 bar.
         ! p12_boundary = 1d11 ! Cut-off at 10^5 bar.

         !--------------------------------------------------------------
         ! Radial profiles for the new heating terms.
         ! (Multiple output files at specific ages.)
         !--------------------------------------------------------------
         open(1,file = "LOGS/radial_profile_last.txt")
         write(1,'(a1,5a18)')"#","Age:","Teq:","Mass:","R_RCB:","P_RCB:"
         write(1,'(a1,5e18.9)')"#",age,Teq,mass_p,r1_rcb,p1_rcb
         write(1,'(a1,a6,25a18)')"#","INDEX","RADIUS [RJ]","RADIUS [RP]", &
         "DENSITY","PRESSURE","TEMPERATURE","ENTROPY","GRAVITY", &
         "MASS","DMASS","CHI_RHO","CHI_T","GAMMA_1","GAMMA_3", &
         "SIGMA1","SIGMA2","EPS_JOULE1","EPS_JOULE2","V_ATM", &
         "SIGMA1_LOG","SIGMA1_EYE","SIGMA1_TSANG","SIGMA1_INT", &
         "LUMINOSITY","EXTRA_HEAT","IRRAD_HEAT"
         write(1,'(a1,a6,25a18)')"#","COL 1","2","3","4","5","6","7","8", &
         "9","10","11","12","13","14","15","16","17","18","19","20","21", &
         "22","23","24","25","26"

         ! Snapshot at (near) 100 Myr.
         if(open_file_1)then
            if(age >= 1d8)then
               open(2,file = "LOGS/radial_profile_0.1Gyr.txt")
               write(2,'(a1,5a18)')"#","Age:","Teq:","Mass:","R_RCB:","P_RCB:"
               write(2,'(a1,5e18.9)')"#",age,Teq,mass_p,r1_rcb,p1_rcb
               write(2,'(a1,a6,25a18)')"#","INDEX","RADIUS [RJ]","RADIUS [RP]", &
               "DENSITY","PRESSURE","TEMPERATURE","ENTROPY","GRAVITY", &
               "MASS","DMASS","CHI_RHO","CHI_T","GAMMA_1","GAMMA_3", &
               "SIGMA1","SIGMA2","EPS_JOULE1","EPS_JOULE2","V_ATM", &
               "SIGMA1_LOG","SIGMA1_EYE","SIGMA1_TSANG","SIGMA1_INT", &
               "LUMINOSITY","EXTRA_HEAT","IRRAD_HEAT"
               write(2,'(a1,a6,25a18)')"#","COL 1","2","3","4","5","6","7","8", &
               "9","10","11","12","13","14","15","16","17","18","19","20","21", &
               "22","23","24","25","26"

               open_file_1 = .FALSE.
               write_file_1 = .TRUE.
            endif
         endif

         ! Snapshot at (near) 1 Gyr.
         if(open_file_2)then
            if(age >= 1d9)then
               open(3,file = "LOGS/radial_profile_1Gyr.txt")
               write(3,'(a1,5a18)')"#","Age:","Teq:","Mass:","R_RCB:","P_RCB:"
               write(3,'(a1,5e18.9)')"#",age,Teq,mass_p,r1_rcb,p1_rcb
               write(3,'(a1,a6,25a18)')"#","INDEX","RADIUS [RJ]","RADIUS [RP]", &
               "DENSITY","PRESSURE","TEMPERATURE","ENTROPY","GRAVITY", &
               "MASS","DMASS","CHI_RHO","CHI_T","GAMMA_1","GAMMA_3", &
               "SIGMA1","SIGMA2","EPS_JOULE1","EPS_JOULE2","V_ATM", &
               "SIGMA1_LOG","SIGMA1_EYE","SIGMA1_TSANG","SIGMA1_INT", &
               "LUMINOSITY","EXTRA_HEAT","IRRAD_HEAT"
               write(3,'(a1,a6,25a18)')"#","COL 1","2","3","4","5","6","7","8", &
               "9","10","11","12","13","14","15","16","17","18","19","20","21", &
               "22","23","24","25","26"

               open_file_2 = .FALSE.
               write_file_2 = .TRUE.
            endif
         endif

         ! Snapshot at (near) 5 Gyr.
         if(open_file_3)then
            if(age >= 5d9)then
               open(4,file = "LOGS/radial_profile_5Gyr.txt")
               write(4,'(a1,5a18)')"#","Age:","Teq:","Mass:","R_RCB:","P_RCB:"
               write(4,'(a1,5e18.9)')"#",age,Teq,mass_p,r1_rcb,p1_rcb
               write(4,'(a1,a6,25a18)')"#","INDEX","RADIUS [RJ]","RADIUS [RP]", &
               "DENSITY","PRESSURE","TEMPERATURE","ENTROPY","GRAVITY", &
               "MASS","DMASS","CHI_RHO","CHI_T","GAMMA_1","GAMMA_3", &
               "SIGMA1","SIGMA2","EPS_JOULE1","EPS_JOULE2","V_ATM", &
               "SIGMA1_LOG","SIGMA1_EYE","SIGMA1_TSANG","SIGMA1_INT", &
               "LUMINOSITY","EXTRA_HEAT","IRRAD_HEAT"
               write(4,'(a1,a6,25a18)')"#","COL 1","2","3","4","5","6","7","8", &
               "9","10","11","12","13","14","15","16","17","18","19","20","21", &
               "22","23","24","25","26"

               open_file_3 = .FALSE.
               write_file_3 = .TRUE.
            endif
         endif

         !--------------------------------------------------------------
         ! Main loop for Ohmic heating terms.
         !--------------------------------------------------------------
         ! Manual moment of inertia calculation.
         ! (This was needed by Albert, but not anymore.)
         ! Moment of inertia is s% i_rot.
         ! (However, this is calculated only when rotation_flag is true.)
         inertia_moment = 0d0

         ! Loop over radial index. (The radial dimension ndim is variable.)
         ! From outside inwards (k=1 corresponds to the surface).
         do k = 1,ndim
            ! All quantities are in cgs units.
            radius = s% r(k)
            rr = radius/radius_p ! Dimensionless radial coordinate.
            mass = s% m(k) ! Mass enclosed [in g].
            dmass = s% dm(k) ! Baryonic mass of cell k.
            density = s% rho(k)
            ! There is also s% entropy(k), which seems to be something different.
            entropy = exp(s% lnS(k)) ! Log of S (specific entropy).
            ! Gravitational acceleration: g = standard_cgrav*mass/radius**2.
            gravity = s% grav(k)
            pressure = s% P(k) ! Total pressure (P = P_rad + P_gas).
            temperature = s% T(k)
            luminosity = s% L(k) ! For output purpose only.

            !-----------------------------------------------------------
            ! Manual moment of inertia calculation. (Was needed by Albert.)
            !-----------------------------------------------------------
            if(k.ne.ndim)then
               radius_1 = s% r(k+1)
               inertia_moment = inertia_moment + (8d0*pi/3d0) * density * radius**4 * (radius - radius_1)
            endif

            !-----------------------------------------------------------
            ! Joule 1:
            ! Ohmic heating in the interior.
            !-----------------------------------------------------------
            current1 = 0d0
            joule_heating1 = 0d0
            sigma1 = 0d0 ! In cgs units. (1 S/m = 8.98755 x 10^9 1/s.)
            sigma1_log = 0d0
            sigma1_eye = 0d0
            sigma1_tsang = 0d0
            sigma1_int = 0d0

            !-----------------------------------------------------------
            ! Electrical conductivity profile.
            !-----------------------------------------------------------
            ! (1) Logarithmic profile.
            !-----------------------------------------------------------
            ! (1a) Linear in radius (on a log-linear scale).
            ! log_sigma = - 30d0*rr + 31.2d0 ! Note: 1 S/cm = 100 S/m.
            ! sigma1_log = 8.98755d9*10d0**log_sigma
            ! (2a) Linear in pressure (on a log-log scale).
            sigma1_log = 10d0**((5d0/6d0)*log10(pressure) + 31d0/6d0)
            !-----------------------------------------------------------
            ! (2) Fit by eye.
            !-----------------------------------------------------------
            ! The radial range needs to be restricted to 0.935 for this case.
            if(rr.lt.0.935d0)then
               log_sigma_eye = 2d0 + 6.02d0*(0.935d0 - rr)**0.16d0
               sigma1_eye = 8.98755d9*10d0**log_sigma_eye
            endif
            !-----------------------------------------------------------
            ! (3) Fit of Tsang & Jones (2020) (eq.14 and fig.2).
            !-----------------------------------------------------------
            ! Radius of Jupiter used by Tsang & Jones (2020), eq.9 [in m].
            rj_tj20 = 6.9894d7
            ! Permeability of free space.
            ! mu_0 = 1.25663706143d-6
            mu_0 = 4d0*pi*1d-7
            ! Coefficients.
            c1 = -4.279d-6
            c2 = 274.9d0
            c3 = -2.544d-8
            c4 = 1.801d0
            c5 = 20.28d0
            aa = 1d0
            bb = c1*rr*rj_tj20 + c2 + c3*rr*rj_tj20 + c4
            cc = (c1*rr*rj_tj20 + c2)*(c3*rr*rj_tj20 + c4) - c5
            log_sigma_tsang = (-bb + sqrt(bb**2 - 4d0*aa*cc))/(2d0*aa)
            sigma1_tsang = 8.98755d9/mu_0/exp(log_sigma_tsang)
            !-----------------------------------------------------------
            ! (4) Conductivity from an external table.
            !-----------------------------------------------------------
            ! This is not an optimal implementation (but rather, somewhat of a brute force).
            sigma1_int = sigma_tab(ndim_tab)
            if(p_tab(1).gt.pressure)then
               sigma1_int = sigma_tab(1)
               goto 100
            endif
            do k_tab = 2,ndim_tab
               if(p_tab(k_tab).gt.pressure)then
                  sigma1_int = sigma_tab(k_tab) - (sigma_tab(k_tab) - sigma_tab(k_tab-1))* &
                     (p_tab(k_tab) - pressure) / (p_tab(k_tab) - p_tab(k_tab-1))
                  goto 100
               endif
            enddo
100         continue

            ! Assign a sigma profile.
            sigma1 = sigma1_log
            ! sigma1 = sigma1_int

            !-----------------------------------------------------------
            ! Ohmic heating in the interior.
            !-----------------------------------------------------------
            ! Current from Ampere's law: J = c/(4 pi) x curl B.
            ! length_scale1 = radius_p/1d3 ! Length scale of curl B [in cm].
            length_scale1 = radius_p*length_scale_input
            current1 = (c_cgs/4d0/pi) * B_dynamo / length_scale1

            if(pressure.ge.p12_boundary)then ! Always below the pressure cutoff line.
            ! if(rr.ge.0.84d0.and.rr.le.0.94d0)then ! Old radial range.
               ! Joule heating: Q = J^2/sigma.
               joule_heating1 = current1**2 / sigma1
            endif

            !-----------------------------------------------------------
            ! Joule 2:
            ! Ohmic heating near the surface (atmospheric circulation).
            ! Pressure range for Joule heating.
            ! Linear and non-linear regimes as functions of temperature.
            !-----------------------------------------------------------
            current2 = 0d0
            J2_l = 0d0
            J2_n = 0d0
            joule_heating2 = 0d0
            Q2_l = 0d0
            Q2_n = 0d0
            ! Temperature dependent amplitude.
            ! Normalize the amplitude to 1 for T = 1000K.
            !sigma2_amp = 10d0**(9d0 - (4d0 - Teq/1d3)**2)
            ! Reference temperature for the amplitude 1400 K.
            ! Conductivity as a function of equilibrium temperature (Teq = 0 is problematic).
            ! sigma2_amp =  3d-1 * 10d0**(9d0 - (4d0 - 1250d0/1d3)**2) &
            !              * (Teq/1400d0)**0.75d0 * exp(25188d0*(1d0/1400d0 - 1d0/Teq))
            ! Conductivity as a function of temperature profile in the interior.
            sigma2_amp =  3d-1 * 10d0**(9d0 - (4d0 - 1250d0/1d3)**2) &
                          * (temperature/1400d0)**0.75d0 * exp(25188d0*(1d0/1400d0 - 1d0/temperature))
            sigma2 = 0d0 ! In cgs units. (1 S/m = 9x10^9 1/s, approximately.)
            ! Assume max. amplitude of v = 4x10^5 cm/s for T = 4000K.
            ! v is a function of Teq and pressure. (Vanishes for Teq = 0.)
            v_amp = Teq/1d3
            v_atm = 0d0

            ! The numeric surface pressure is just below 10^5 Ba (i.e. 0.1 bar).
            if(pressure.lt.p12_boundary)then
            ! if(pressure.ge.1d3.and.pressure.lt.1d9)then ! Pressure range.

               ! Electrical conductivity profile (in units of s^-1).
               ! Linear logarithmic profile (cf. fig.3 of Batygin et al. 2011).
               ! Note, 1 bar = 10^5 Pa (Pascal) [SI] = 10^6 Ba (barye) [cgs].
               log_sigma = 7d0/6d0*log10(pressure/1d6) - 15d0/6d0 + 9d0
               ! Multiply by some extra amplitude.
               sigma2 = sigma2_amp * 10d0**log_sigma

               ! Wind speed of 1 km/s = 10^5 cm/s at the top (at 0.001 bar).
               ! (A1) Linear on log-log scale.
               ! log_v = -5d0/9d0*log10(pressure/1d6) + 10d0/3d0
               ! v_atm = 8d0*10d0**log_v
               ! (A2) Old version for linear velocity profile.
               ! log_v = -2d0/7d0*log10(pressure/1d6) - 1d0/7d0
               ! v_atm = 4d3*10d0**log_v
               ! (B) Hyperbolic profile.
               ! Limiting value for P -> 0 is v_amp*1d5.
               ! For Jupiter approx. Teq = 100K, so we assume max. v = 1d4 cm/s.
               v_atm = v_amp * 5d4 * (1d0 + tanh((-0.5d0 - log10(pressure/1d6))/0.5d0))

               ! Current from Ohm's law: J = sigma v x B / c.
               ! Linear regime.
               Bphi2_l = B_dipole ! Assume magnetic Reynolds number of 1.
               J2_l = sigma2 * v_atm * Bphi2_l/c_cgs
               ! Non-linear regime.
               Bphi2_n = v_atm * (4d0 * pi * density)**0.5d0
               length_scale2 = radius_p/1d2 ! Length scale of curl B [in cm].
               J2_n = c_cgs * Bphi2_n / (4d0 * pi * length_scale2)

               ! Joule heating: Q = J^2/sigma.
               ! joule_heating2 = current2**2 / sigma2
               Q2_l = J2_l**2 / sigma2
               Q2_n = J2_n**2 / sigma2
               exp_fn = exp((Q2_l - Q2_n) / (Q2_l + Q2_n))
               joule_heating2 = exp(1d0) / (1d0 / exp_fn / Q2_l + exp_fn / Q2_n)
            endif

            !-----------------------------------------------------------
            ! Contribution to luminosity due to Joule heating.
            !-----------------------------------------------------------
            eps_joule1 = joule_heating1/density
            eps_joule2 = joule_heating2/density
            ! eps_joule = 0d0
            eps_joule = eps_joule1
            ! eps_joule = eps_joule2
            ! eps_joule = eps_joule1 + eps_joule2

            ! The same amount of heating anywhere below the RCB gives the same evolution.
            ! eps_joule = 0d0
            ! if(rr.le.0.5d0.and.lum_joule.le.4d28)then
            !    eps_joule = 0.05d0
            ! endif

            ! Luminosity (i.e. power) due to irradiation.
            eps_irrad = s% irradiation_heat(k) ! In erg/g/sec.
            lum_irrad = lum_irrad + eps_irrad*dmass

            ! Luminosity due to Joule heating (integrated over mass).
            lum_joule = lum_joule + eps_joule*dmass

            ! Joule heating above and below the RCB (as defined by region-1).
            if(radius.gt.r1_rcb)then
               lum_joule_above = lum_joule_above + eps_joule*dmass
            else
               lum_joule_below = lum_joule_below + eps_joule*dmass
            endif

            !-----------------------------------------------------------
            ! Update the component of s with the extra heat source.
            !-----------------------------------------------------------
            s% extra_heat(k) = eps_joule ! In erg/g/sec.

            !-----------------------------------------------------------
            ! Write radial output.
            !-----------------------------------------------------------

            !-----------------------------------------------------------
            ! Radial profiles for various quantities (at the last step).
            ! Pressure is in cgs units (1 Ba, barye = 1e-6 bar).
            !-----------------------------------------------------------
            write(1,'(i7,30e18.9)')k,radius/radius_j,rr,density,pressure, &
            temperature,entropy,gravity,mass,dmass, &
            s% chiRho(k),s% chiT(k),s% gamma1(k),s% gamma3(k), &
            sigma1,sigma2,eps_joule1,eps_joule2,v_atm, &
            sigma1_log,sigma1_eye,sigma1_tsang,sigma1_int, &
            luminosity,s% extra_heat(k),s% irradiation_heat(k)

            ! Note: Heat equation terms.
            ! eps_heat = eps_nuc - non_nuc_neu + extra_heat + irradiation_heat
            ! Here eps_nuc is the sum for all reactions (power per unit mass),
            ! and non_nuc_neu is for the non-nuclear reaction neutrino losses.
            ! The heat terms are not summed until after this subroutine.
            ! (So, the total term s% eps_heat is not updated at this point!)
            ! Other terms: s% eps_heat(k),s% eps_nuc(k),s% non_nuc_neu(k)

            if(write_file_1)then
               write(2,'(i7,30e18.9)')k,radius/radius_j,rr,density,pressure, &
               temperature,entropy,gravity,mass,dmass, &
               s% chiRho(k),s% chiT(k),s% gamma1(k),s% gamma3(k), &
               sigma1,sigma2,eps_joule1,eps_joule2,v_atm, &
               sigma1_log,sigma1_eye,sigma1_tsang,sigma1_int, &
               luminosity,s% extra_heat(k),s% irradiation_heat(k)
            endif

            if(write_file_2)then
               write(3,'(i7,30e18.9)')k,radius/radius_j,rr,density,pressure, &
               temperature,entropy,gravity,mass,dmass, &
               s% chiRho(k),s% chiT(k),s% gamma1(k),s% gamma3(k), &
               sigma1,sigma2,eps_joule1,eps_joule2,v_atm, &
               sigma1_log,sigma1_eye,sigma1_tsang,sigma1_int, &
               luminosity,s% extra_heat(k),s% irradiation_heat(k)
            endif

            if(write_file_3)then
               write(4,'(i7,30e18.9)')k,radius/radius_j,rr,density,pressure, &
               temperature,entropy,gravity,mass,dmass, &
               s% chiRho(k),s% chiT(k),s% gamma1(k),s% gamma3(k), &
               sigma1,sigma2,eps_joule1,eps_joule2,v_atm, &
               sigma1_log,sigma1_eye,sigma1_tsang,sigma1_int, &
               luminosity,s% extra_heat(k),s% irradiation_heat(k)
            endif
         enddo
         ! End of loop over radial index.

         ! Close the opened files (for the radial profiles).
         close(1)

         if(write_file_1)then
            write_file_1 = .FALSE.
            close(2)
         endif

         if(write_file_2)then
            write_file_2 = .FALSE.
            close(3)
         endif

         if(write_file_3)then
            write_file_3 = .FALSE.
            close(4)
         endif

         !--------------------------------------------------------------
         ! Output of quantitites as a function of age.
         !--------------------------------------------------------------
         if(n_call.eq.0)then
            open(11,file = "LOGS/evolution.txt")
            write(11,'(a22,e16.7)')"# SOLAR RADIUS (RSUN):",Rsun
            write(11,'(a22,e16.7)')"# SOLAR MASS (MSUN):  ",Msun
            write(11,'(a1,a15,8a16)')"#","AGE","RADIUS","MASS","T_EFF","T_EQ", &
            "LUMINOSITY","B_DYNAMO","B_DIPOLE","INERTIA"
            write(11,'(a1,a15,8a16)')"#","[YEARS]","[RSUN]","[MSUN]","[K]","[K]", &
            "[LSUN]","[GAUSS]","[GAUSS]","[CGS]"
            write(11,'(a1,a15,8a16)')"#","I","II","III","IV","V","VI","VII","VIII","IX"
         else
            ! In old versions of Fortran ACCESS was accepted for appending files.
            ! However, in new versions POSITION should be used instead.
            open(11,file = "LOGS/evolution.txt",position = "append")
            write(11,'(10e16.7)')age,radius_p/Rsun,mass_p/Msun,Teff,Teq, &
            luminosity_p,B_dynamo,B_dipole,inertia_moment
         endif
         close(11)

         ! Further output.
         if(n_call.eq.0)then
            open(12,file = "LOGS/evolution_extra.txt")
            write(12,'(a1,a15,2a16,a32)')"#","AGE","RADIUS","MASS","HEAT TERMS..."
         else
            open(12,file = "LOGS/evolution_extra.txt",position = "append")
            write(12,'(15e16.7,2i8)')age,radius_p/Rsun,mass_p/Msun, &
            luminosity_p/Lsun,lum_irrad,lum_joule, &
            lum_joule_above,lum_joule_below, &
            s% total_extra_heating,s% total_irradiation_heating, &
            s% cumulative_extra_heating,s% cumulative_irradiation_heating, &
            s% cumulative_extra_heating_old,s% cumulative_irradiation_heating_old, &
            s% r(ndim)/Rsun,ndim,ntot
         endif
         close(12)

         ! Convective regions.
         if(n_call.eq.0)then
            open(13,file = "LOGS/evolution_conv.txt")
            write(13,'(a22,e16.7)')"# SOLAR RADIUS (RSUN):",Rsun
            write(13,'(a22,e16.7)')"# SOLAR MASS (MSUN):  ",Msun
            write(13,'(a1,a15,2a16,a32)')"#","AGE (YRS)","RADIUS (RSUN)", &
            "MASS (MSUN)","CONVECTIVE REGIONS..."

            open(14,file = "LOGS/evolution_conv2.txt")
            write(14,'(a22,e16.7)')"# SOLAR RADIUS (RSUN):",Rsun
            write(14,'(a22,e16.7)')"# SOLAR MASS (MSUN):  ",Msun
            write(14,'(a1,a15,2a16,a48)')"#","AGE (YRS)","RADIUS (RSUN)", &
            "MASS (MSUN)","CONVECTIVE REGIONS... (See code for columns.)"
         else
            open(13,file = "LOGS/evolution_conv.txt",position = "append")
            write(13,'(11e16.7)')age,radius_p/Rsun,mass_p/Msun, &
            s% conv_mx1_top_r,s% conv_mx1_bot_r, &
            s% conv_mx2_top_r,s% conv_mx2_bot_r, &
            s% mx1_top_r,s% mx1_bot_r,s% mx2_top_r,s% mx2_bot_r

            open(14,file = "LOGS/evolution_conv2.txt",position = "append")
            write(14,'(5e16.7,3(2e16.7,i8))')age,radius_p/Rsun,mass_p/Msun, &
            p_surface,p12_boundary,r1_rcb/Rsun,p1_rcb,n1_rcb, &
            r3_rcb/Rsun,p3_rcb,n3_rcb,r4_rcb/Rsun,p4_rcb,n4_rcb
         endif
         close(13)
         close(14)

         ! Update saved counters.
         n_call = n_call + 1
         if(first_call) first_call = .FALSE.

         ! Note: Central (core) quantites.
         ! L_center and M_center are defined at the inner radius R_center.
         ! These are added to the values of L, M and R at ndim (innermost index).
         ! At ndim, the values for ndim+1 are replaced by those of "center":
         ! i.e., m(ndim) = dm(ndim) + M_center (and similarly for L & R).
         ! At any other k: m(k) = m(k+1) + dm(k).
         ! Beyond ndim is not used: m(ndim+1) = 0.
         ! The luminosity at k is: L(k) = L(k+1) + dm(k)*eps.
         ! Some output checks:
         ! write(*,'(3e16.7)')s% L_center,s% M_center,s% R_center
         ! write(*,'(3e16.7)')10d0*m_earth,s% M_center
         ! write(*,'(3e16.7)')s% m(ndim),s% dm(ndim) + s% M_center

         ! write(*,*)"Disable Ohmic heating in inlist_create and/or hard stop!"
         ! stop

         return

      end subroutine default_other_energy

end module run_star_extras
