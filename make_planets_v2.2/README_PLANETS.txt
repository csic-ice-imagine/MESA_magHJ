# Make Planets v2.2 Notes

# Python scripts.
- Main script for sequential runs: python_run.py
- Subprograms needed by the main script: python_subs.py
- Script for adding # to comments in output (for gnuplot): python_add_hash.py

# Other scripts.
- A cleaning script (invoked by python_run.py, or can be run separately).
- Various MESA scripts (ck, clean, mk, re, rn and rn1).

# Input files.
There are three inlist files which are run consecutively (see the rn script):
- inlist_create creates a planet;
- inlist_core adds a core;
- inlist_evolve follows the subsequent evolution.
The file inlist is overwritten by the script rn.

The irradiation (due to the host star) is added only for the evolution part.
(Not for the creation and core parts of the code.)

# Output customization.
A customized history_columns_file, defining the columns in the history.data
output file (under the LOGS folder), has to be defined in all three inlists.
The history_columns_file parameter is part of the star_job namelist. These
are defined in the main MESA folder, in a file in the star/private folder.

The sample history_columns_file provided includes further explanations of
various output quantitites.

For sequential runs, all output in LOGS is copied into a custom folder out/
created by the script python_runs.py (only when there is more than one run).

# Definitions of constants.
Important constants are defined in ${MESA_DIR}/const/public/const_def.f90, including
Msun, m_jupiter and m_earth:

standard_cgrav = 6.67430d-8 ! Gravitational constant G (g^-1 cm^3 s^-2).

! Definitions according to IAU 2015 Resolution B3.
! Standard gravitational parameters = G*M (units cm^3 s^-2).
mu_sun = 1.3271244d26
mu_earth = 3.986004d20
mu_jupiter = 1.2668653d23

Msun = mu_sun / standard_cgrav ! Solar mass (g).
Rsun = 6.957d10 ! Solar radius (cm).
Lsun = 3.828d33 ! Solar luminosity (erg s^-1).

m_earth = mu_earth/standard_cgrav ! Earth mass (g).
r_earth = 6.3781d8 ! Earth equatorial radius (cm).
r_earth_polar = 6.3568d8 ! Earth polar radius (cm).

m_jupiter = mu_jupiter/standard_cgrav ! Jupiter mass (g).
r_jupiter = 7.1492d9 ! Jupiter equatorial radius (cm).
r_jupiter_polar = 6.6854d9 ! Jupiter polar radius (cm).

# Mixing. (Summary from the file history_columns.list for output customization.)
! Mixing regions.
! mx1 refers to the largest (by mass) convective region. mx2 is the 2nd largest.
! conv_mx1_top and conv_mx1_bot are the region where mixing_type == convective_mixing.
! mx1_top and mx1_bot are the extent of all kinds of mixing, convective and other.

! Values in mass (m/Mstar): conv_mx1_top, conv_mx1_bot, mx1_top, mx1_bot, etc.
! Values in radius (in Rsun units): conv_mx1_top_r, conv_mx1_bot_r, etc.

! Mixing types are: no_mixing, convective_mixing, overshoot_mixing,
! semiconvective_mixing and salt_finger_mixing (numbered from 0 to 4).

# Other notes.
! Avoid using the file unit 6. After closing with the statement close(6) some of
! the terminal output might be directed to a file fort.6 instead.
