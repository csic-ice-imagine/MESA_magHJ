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
         type (star_info), pointer :: s     ! Collects basically all the variables from MESA main program
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
         use const_def, only: Lsun,Msun,Rsun,m_jupiter,r_jupiter, &
                              boltz_sigma,clight,me,dp,pi 

         ! Default definitions.
         implicit none
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ! Constants
         real(dp), parameter :: m_u = 1.660539d-24  ! Atomic mass unit (in g).
         real(dp), parameter :: e_e = 4.8032d-10 ! Elementary charge in cgs units (g^0.5 cm^1.5 s^-1).
         real(dp), parameter :: p_dynamo = 1d6   ! Pressure at which we consider the dynamo to act, i.e. with metallic H (bar)

         ! Conversion factors
         real(dp), parameter :: sigma_SI_to_cgs = clight**2d0 / 1d11 ! Conversion factor for conductivity from S/m (SI) to cgs (s^-1).
         real(dp), parameter :: p_cgs_to_bar = 1d-6 ! Pressure (1 barye = 1e-6 bar).

         ! Custom (local) definitions.
         integer :: k,ndim
         real(dp) :: age,density,pressure,temperature,entropy,gravity
         real(dp) :: radius,rr,mass,dmass,luminosity       ! Variables of the structure equations
         real(dp) :: radius_p,luminosity_p                 ! Planetary radius (cm) quantities (in cgs).
         real(dp) :: r_dynamo,B_dynamo,B_dip_surface                ! Inferred values of B dynamo and B_dip_surface
         real(dp) :: luminosity_int ! Magnetic field and internal luminosity just above the dynamo region (1e6 bar). CHECK PROFILES OF LUMINOSITIES IN THE OUTPUT!
         real(dp), save :: X_frac, Y_frac, Z_frac, mass_p ! Composition and mass (fixed)

         ! Irradiation and Ohmic terms
         real(dp) :: eps_irrad,eps_joule ! Specific heat from irradiation and Joule effect (erg/g/s)
         real(dp) :: lum_irrad,lum_joule,lum_joule_above,lum_joule_below ! Integrated luminosity for irradiation, Joule, Joule above/below RCB (erg/s)
         real(dp), save :: irradiation,column_depth_for_irradiation,Teq
         real(dp) :: x_e,sv_e,f_time_ohm

         ! Parameters read from the input (see below for their description)
         real(dp), save :: J0, pcmin, pcmax, s0
         real(dp), save :: a_K, mu_mean, age_ohmic_on, age_ohmic_on_full
         real(dp) :: joule_heating, sigma, window

         ! Pressures and radii of the main (1) and superficial (2, if any)
         ! radiative-convective boundaries (RCB).
         real(dp) :: r_top_rcb1,r_top_rcb2,r_bot_rcb2
         real(dp) :: p_top_rcb1,p_top_rcb2,p_bot_rcb2
         real(dp) :: r_top_rcb1_old,r_top_rcb2_old,r_bot_rcb2_old
         real(dp) :: p_top_rcb1_old,p_top_rcb2_old,p_bot_rcb2_old
         real(dp) :: p_surface
         ! index for the presence(1)/absence(0) of the convective region
         integer :: conv_sup

         ! Variables for the outputs
         logical, save :: first_call = .TRUE.
         logical, save :: write_profile = .FALSE.
         integer, save :: ia = 1
         integer :: iaa
         ! Maximum age of the output and number of outputs
         ! (set like this to have the same f3.1 format)
         real(dp), parameter :: max_age_out = 9.9d9
         integer, parameter :: nage = 99
         real(dp), dimension(nage), save :: age_output
         character(len = 7), parameter :: AGE_STR_FORMAT = "(f3.1)"
         character(len = 3) :: str
         character(len = 3) :: tmp

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
         ndim = s% nz ! Number of radial zones (variable).
         age = s% star_age ! Current age [yrs]

         if (first_call) then

           do iaa = 1,nage
             age_output(iaa) = iaa * (max_age_out) / nage
           end do
    
           mass_p = s% m(1)   ! Total mass of the planet (outermost value, constant in time)
           irradiation = s% irradiation_flux      ! Stellar irradiation flux [erg cm^-2 s^-1]
           column_depth_for_irradiation = s% column_depth_for_irradiation  ! Column depth for irradiation deposition [cm^2 g^-1]
           ! Composition (constant in time and space)
           X_frac = s% X(1)   ! Mass fraction of neutral Hydrogen
           Y_frac = s% Y(1)   ! Mass fraction of neutral Helium
           Z_frac = s% Z(1)   ! Mass fraction of neutral heavier elements

           ! Derived quantities
           ! Mean molecular weight, assuming for heavy elements Z/A=1/2 
           ! and negligible H-dissociation and negligible ionization
           mu_mean = 1d0/(0.5d0*X_frac + 0.25d0*Y_frac + Z_frac*0.5d0)

           ! Equilibrium temperature (K) assuming perfect redistribution
           Teq = (irradiation/4d0/boltz_sigma)**0.25d0

           ! Read Ohmic input parameters of the models.
           open(1,file = "in/input_ohmic.txt")
           read(1,*) J0        ! Current density (statamp/cm^2) for a sigma=1 S/m - Order of magnitude: 1e3
           read(1,*) pcmin     ! Minimum pressure (bar) at which the Ohmic dissipation is considered - Std: 1
           read(1,*) pcmax     ! Maximum pressure (bar) at which the Ohmic dissipation is considered - Std: 1e4
           read(1,*) s0        ! Width of the tanh of the deposition region, i.e., how smooth is the beginning and end of the region (log(p)) - std: 0.25
           read(1,*) age_ohmic_on       ! Age at which we gradually start switching on the Ohmic term - std: 1d5
           read(1,*) age_ohmic_on_full  ! Age from which the Ohmic term is fully on - std: 1e6
           read(1,*) a_K       ! Potassium abundance - Solar value: 1e-7
           close(1)

           ! Convert pressures from bar to cgs
           pcmin = pcmin/p_cgs_to_bar
           pcmax = pcmax/p_cgs_to_bar
         endif

         ! Planetary luminosity, mass and radius (in cgs).
         ! The surface values are at index 1.
         luminosity_p = s% L(1)
         radius_p = s% r(1)

         ! Initialize variables
         joule_heating = 0d0
         sigma = 0d0
         lum_irrad = 0d0
         lum_joule = 0d0
         lum_joule_below = 0d0
         lum_joule_above = 0d0
         luminosity_int = 0d0

         !--------------------------------------------------------------
         ! Determine the pressure at the RCB (linear interpolation).
         ! Here we track the radius and pressure delimiting the main 
         ! and superficial (if any) convective regions
         ! The main region bottom is the deepest layer by definition
         !--------------------------------------------------------------
         ! Define the surface pressure (at outermost cell).
         p_surface = s% P(1)
         ! Initial values.
         r_top_rcb1_old = radius_p
         r_top_rcb2_old = radius_p
         r_bot_rcb2_old = radius_p
         p_top_rcb1_old = p_surface
         p_top_rcb2_old = p_surface
         p_bot_rcb2_old = p_surface
         r_top_rcb1 = (s% conv_mx1_top_r) * Rsun   ! top of the main convective region, converted to cm
         r_top_rcb2 = (s% conv_mx2_top_r) * Rsun   ! top of the possible second convective region (cm)
         r_bot_rcb2 = (s% conv_mx2_bot_r) * Rsun   ! bottom of the possible second convective region (cm)
         p_top_rcb1 = 0d0   ! corresponding pressures
         p_top_rcb2 = 0d0
         p_bot_rcb2 = 0d0
         conv_sup = 0     ! 1 if there is a secondary convective region

         ! Loop to calculate the pressures of the convective regions
         ! Loop over radial index. (The radial dimension ndim is variable.)
         ! From outside inwards (k=1 corresponds to the surface).
         do k = 1,ndim
            radius = s% r(k)
            pressure = s% P(k) ! Total pressure (P = P_rad + P_gas).
            luminosity = s% L(k)

            ! Determine the luminosity and radius of the top of the dynamo region
            if(pressure.le.p_dynamo/p_cgs_to_bar)then
               luminosity_int = luminosity
               r_dynamo = radius
            endif

            ! main convection region (extending from r_top_rcb1 down to the center).
            if(r_top_rcb1_old.gt.r_top_rcb1.and.radius.le.r_top_rcb1)then
               p_top_rcb1 = p_top_rcb1_old + (pressure - p_top_rcb1_old) * (r_top_rcb1 - r_top_rcb1_old) / (radius - r_top_rcb1_old)
            endif
            r_top_rcb1_old = radius
            p_top_rcb1_old = pressure

            ! Secondary convection region (from r_top_rcb2 down to r_bot_rcb2).
            if(r_top_rcb2.ne.0d0.and.r_top_rcb2_old.gt.r_top_rcb2.and.radius.le.r_top_rcb2)then
               p_top_rcb2 = p_top_rcb2_old + (pressure - p_top_rcb2_old) * (r_top_rcb2 - r_top_rcb2_old) / (radius - r_top_rcb2_old)
               conv_sup = 1 ! It's one if there is a secondary convective region
            endif
            r_top_rcb2_old = radius
            p_top_rcb2_old = pressure

            if(r_bot_rcb2.ne.0d0.and.r_bot_rcb2_old.gt.r_bot_rcb2.and.radius.le.r_bot_rcb2)then
               p_bot_rcb2 = p_bot_rcb2_old + (pressure - p_bot_rcb2_old) * (r_bot_rcb2 - r_bot_rcb2_old) / (radius - r_bot_rcb2_old)
            endif
            r_bot_rcb2_old = radius
            p_bot_rcb2_old = pressure

         enddo

         !--------------------------------------------------------------
         ! Magnetic field strength from scaling laws.
         !--------------------------------------------------------------

         ! Dynamo magnetic field [in Gauss] from Reiners & Christensen 2010, eq. 1
         ! We use the internal lumonisity at the top of the dynamo region
         B_dynamo = 4.8d3*((mass_p/Msun)*(luminosity_int/Lsun)**2d0/(radius_p/Rsun)**7d0)**(1d0/6d0)

         ! Extrapolation of the dipolar component to the surface of the planet,
         ! assuming a dynamo surface at p_dynamo
         B_dip_surface = B_dynamo*(r_dynamo/radius_p)**3d0

         !--------------------------------------------------------------
         ! Radial profiles header
         !--------------------------------------------------------------
         ! Writing the header
         if(ia .le. nage .and. age >= age_output(ia) .and. age <= max_age_out)then
               write(tmp, AGE_STR_FORMAT) age/1d9
               str = trim(tmp)
               open(2,file = "LOGS/radial_profile_"//str//"Gyr.txt")

               ! Header with general information               
               write(2,'(a1,6a18)')"#","Age","Teq","Mass","Radius","R_rcb","P_rcb"
               write(2,'(a1,6a18)')"#","[yr]","[K]","[Mj]","[Rj]","[Rj]","[bar]"
               write(2,'(a1,6es18.9)')"#",age,Teq,mass_p/m_jupiter, &
               radius_p/r_jupiter,r_top_rcb1/r_jupiter,p_top_rcb1*p_cgs_to_bar
               
               write(2,'(a1,a6,25a18)')"#","1.Index","2.Radius","3.Radius", &
               "4.Density","5.Pressure","6.Temperature","7.Entropy","8.Gravity", &
               "9.Mass","10.dm","11.Chi_rho","12.Chi_T","13.Gamma_1","14.Gamma_3", &
               "15.Lum_int","16.Spec_Irr_heat",&
               "17.Spec_Ohm_heat","18.Vol_Ohm_heat","19.Conductivity"
               
               write(2,'(a1,a6,25a18)')"#","","[Rj]","[Rp]", &
               "[g cm^-3]","[bar]","[K]","Entropy","[cm s^-2]", &
               "[Mj]","[Mj]","Chi_rho","Chi_T","Gamma_1","Gamma_3", &
               "[erg s^-1]","[erg s^-1 g^-1]","[erg s^-1 g^-1]",&
               "[erg s^-1 cm^-3]","[s^-1]"

               write(2,'(a1,a6,25a18)')"#","COL 1","2","3","4","5","6","7","8", &
               "9","10","11","12","13","14","15","16","17","18","19"      

               write_profile = .TRUE.
         endif

         ! Defining F(t): gradual increase of the Ohmic term to avoid
         ! unrealistic higher values of J at early age, where sigma
         ! turns out to be very high due to the temperature
         f_time_ohm = 1d0
         if (age .lt. age_ohmic_on_full .and. age .gt. age_ohmic_on) then
            f_time_ohm = (age - age_ohmic_on)/(age_ohmic_on_full-age_ohmic_on)
         elseif (age .lt. age_ohmic_on_full) then
            f_time_ohm = 0d0
         endif

         !--------------------------------------------------------------
         ! main loop for Ohmic heating terms.
         !--------------------------------------------------------------
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
            ! Electrical conductivity profile sigma (s^-1).
            ! Recipe used in Perna et al. 2010 and following works
            !-----------------------------------------------------------
            
            ! Fraction of electrons assuming only contributions from Potassium, Perna et al. 2010 Eq.1
            x_e = 6.47d-13 * (a_K * 1d7)**0.5d0 * (temperature / 1d3)**0.75d0 &
                  * (2.4d15 / (density / mu_mean / m_u))**0.5d0 &
                  * exp(-25188d0 / temperature) / 1.15d-11

            ! Eq.2 (See eq.14 of Draine et al. 1983 for the correct units!)
            sv_e = 1d-15 * (128d0 * boltzm * temperature / (9d0 * pi * me))**0.5d0
            ! Conductivity (eq.3) (s^-1)
            sigma = x_e * e_e**2 / me / sv_e


            ! Region in which we consider the deposition of heat

            ! Joule heating term: Q = J^2/sigma * F(t) * window
            ! J = J0*(sigma/sigma_ref), where we take 1S/m as sigma_ref = 1 S/m = sigma_SI_to_cgs
            ! J and J0 are in statamp/cm^2 = (g^1/2 cm^-1/2 s^-2), 1 statamp = 10/clight A

            joule_heating = (J0/sigma_SI_to_cgs)**2d0*sigma

            ! Multiplying by a smoothened radial window function and by the activation function F(t)
            window = 0.25d0*(1+dtanh((dlog(pressure)-dlog(pcmin))/s0))*(1+dtanh((dlog(pcmax)-dlog(pressure))/s0))
            joule_heating = joule_heating * window * f_time_ohm

            ! Specific irradiation and Joule heating rates (erg g^-1 s^-1).
            eps_irrad = s% irradiation_heat(k)
            eps_joule = joule_heating/density

            ! Update the component of s with the extra heat source.
            s% extra_heat(k) = eps_joule

            ! Luminosity due to irradiation and Joule heating (erg s^-1)
            ! obtained by integrating specific heating rates in mass
            lum_irrad = lum_irrad + eps_irrad*dmass
            lum_joule = lum_joule + eps_joule*dmass

            ! Joule heating above and below the main RCB.
            if(radius.gt.r_top_rcb1)then
               lum_joule_above = lum_joule_above + eps_joule*dmass
            else
               lum_joule_below = lum_joule_below + eps_joule*dmass
            endif

            !-----------------------------------------------------------
            ! Write radial output.
            !-----------------------------------------------------------

            ! Note: Heat equation terms.
            ! eps_heat = eps_nuc - non_nuc_neu + extra_heat + irradiation_heat
            ! Here eps_nuc is the sum for all reactions (power per unit mass),
            ! and non_nuc_neu is for the non-nuclear reaction neutrino losses.
            ! The heat terms are not summed until after this subroutine.
            ! (So, the total term s% eps_heat is not updated at this point!)
            ! Other terms: s% eps_heat(k),s% eps_nuc(k),s% non_nuc_neu(k)
        	   if(write_profile)then
               write(2,'(i7,30es18.9)')k,radius/r_jupiter,rr,density,pressure*p_cgs_to_bar, &
               temperature,entropy,gravity,mass,dmass, &
               s% chiRho(k),s% chiT(k),s% gamma1(k),s% gamma3(k), &
               luminosity,s% irradiation_heat(k),s% extra_heat(k),joule_heating, &
               sigma
            endif
          
         enddo
         ! End of loop over radial index.

         ! Close the opened files (for the radial profiles).
         if(write_profile)then
            write_profile = .FALSE.
            close(2)
            ia = ia + 1
         endif

         !--------------------------------------------------------------
         ! Output of quantitites as a function of age.
         !--------------------------------------------------------------
         if(first_call)then

            open(11,file = "LOGS/evolution.txt")
            write(11,'(a20,2f8.2)') "# Mass [Mj], Teq [K]: ",mass_p/m_jupiter,Teq
            write(11,'(a1,18a16)')"#","1.Age","2.Radius","3.L_surf","4.Teff",&
            "5.L_int","6.R_dynamo","7.B_dynamo","8.B_dip_surf",&
            "9.L_irr","10.L_joule","11.L_joule>RCB","12.L_joule<RCB"

            write(11,'(a1,18a16)')"#","[yr]","[Rj]","[erg s^-1]","[K]",&
            "[erg s^-1]","[Rj]","[G]","[G]","[erg s^-1]","[erg s^-1]",&
            "[erg s^-1]","[erg s^-1]"

            open(12,file = "LOGS/evolution_conv.txt")
            write(12,'(a1,12a16)')"#","1.Age","2.Radius",&
            "3.R_bot main RCB","4.P_bot main RCB","5.R_top main RCB","6.P_top main RCB", &
            "7.R_bot Sh.RCB","8.P_bot Sh.RCB","9.R_top Sh.RCB","10.P_top.Sh.RCB", &
            "11.Surf P","12.Sh.RCB?"

            write(12,'(a1,12a16)')"#","[yr]","[Rj]","[Rj]","[bar]", &
            "[Rj]","[bar]","[Rj]","[bar]","[Rj]","[bar]","[bar]","[1=yes]"

         else
            open(11,file = "LOGS/evolution.txt",position = "append")
            write(11,'(18es16.7)')age,radius_p/r_jupiter, &
            luminosity_p,s% Teff,luminosity_int,r_dynamo/r_jupiter,&
            b_dynamo,b_dip_surface,lum_irrad,lum_joule,lum_joule_above,lum_joule_below

            open(12,file = "LOGS/evolution_conv.txt",position = "append")
            write(12,'(11es16.7,i8)') age,radius_p/r_jupiter, &
            (s% r(ndim))/r_jupiter,(s% P(ndim))*p_cgs_to_bar, &
            r_top_rcb1/r_jupiter,p_top_rcb1*p_cgs_to_bar, &
            r_bot_rcb2/r_jupiter,p_bot_rcb2*p_cgs_to_bar, &
            r_top_rcb2/r_jupiter,p_top_rcb2*p_cgs_to_bar, &
            p_surface*p_cgs_to_bar,conv_sup

         endif

         close(11)
         close(12)
         if(first_call) first_call = .FALSE.

         return

      end subroutine default_other_energy

end module run_star_extras
