# Files & folders (from the official MESA distribution).
# (https://docs.mesastar.org)
mesa-r15140/		- Main MESA folder.
mesasdk/		- MESA software development kit.

# Python scripts for plotting.
# (https://github.com/wmwolf/py_mesa_reader)
py_mesa_reader-master/	- MESA_READER (for plotting).

# Check what version is currently specified in the paths.
echo ${MESA_DIR}

# Main MESA folder.
The principal (core) folder for MESA is ${MESA_DIR}/star/.
The folder star/ contains defaults/ and test_suites/, among other things.
Additional modules (opacities, eos, etc.) are in their own separate folders.

# Make and clean.
The make script is typically included in a file named mk (run by ./mk).
This points to the makefile under the folder make/.
This, in turn, points to the default makefile in the original MESA folder.

# To clean and compile typically run (see scripts for more details):
./clean
./mk

# To run MESA invoke the script rn by:
./rn

# To restart from a "photo" (i.e. binary dump):
./re x137

# PGPLOT.
MESA uses pgplot to plot in real time the evolution of profiles.
inlist_pgstar is read at each timestep, i.e. plot changes can be seen live.

# Input.
Input parameters are listed in "inlist" files (as Fortran namelists).
Inlists are defined in separate files in the folder ${MESA_DIR}/star/private/.
(See star_job and ctrls as examples.)
Default values are defined in the folder ${MESA_DIR}/star/defaults/.

# Customized output.
Optional history and profile columns (for the output files in LOGS/) can be
specified through the parameters history_columns_file and profile_columns_file.
Default history and profile files are provided in the folder star/defaults/.
A customized history_columns_file can be specified in the inlist, under
star_job. (See make_planets as an example.)

# Definitions of constants.
Important constants are defined in ${MESA_DIR}/const/public/const_def.f90,
including Msun, m_jupiter and m_earth.

# User-specified routines.
See the readme file under ${MESA_DIR}/star/other/ for further instructions.
