#!/usr/bin/env python
from pathlib import Path
import math
import python_subs as my
import numpy as np
import os
import shutil
import sys

#-------------------------------------------------------------------------------
# Use the recommended indentation of 4 spaces (as per PEP 8).
# Check scripts using the flake8 command (flake8 --ignore=E265,E501 file).
#-------------------------------------------------------------------------------

# Constants (in cgs).
msun = 1.9884099e33			# Solar mass.
rsun = 6.9598e10			# Solar radius.
mjup = 1.8981246e30			# Jupiter mass (IAU).
rjup = 7.1492e9				# Jupiter radius.
mearth = 5.9721679e27			# Earth mass.
sigma = 5.67e-5				# Stefan-Boltzmann constant.
au = 1.496e13				# Astronomical unit.

# Parameters. (Some of these may be overwritten in the loops below!)
mp = 1.0				# Planet mass in mjup.
rp = 4.0				# Inital planet radius in rjup.
mcore = 10.0				# Planet core mass in mearth.
mcore_array = [0.0, 10.0, 100.0]	# Core mass as an array of values.
rhocore = 10.0				# Core density.
z = 0.02                                # Metallicity of both planet and star.
y = 0.24                                # He fraction of both planet and star.
maxage = 1.1e10                         # Ending age (go slightly over 1e10).

# Parameters related to irradiation. (These are usually redefined below!)
rs = 1.0				# Stellar radius in rsun.
Teff = 5800				# Stellar Teff (needed for the flux calculation).
orb_sep = 0.05				# Orbital sepration in AU (needed for the flux).
irrad_col = 250.0			# Column depth for depositing stellar radiation.
flux_dayside = sigma * Teff**4 * (rs * rsun / orb_sep / au)**2		# Flux hitting planet's dayside.
Teq = (flux_dayside / 4.0 / sigma)**0.25    # Calculated only for output purposes.

# Flags to skip steps.
do_create_planet = True
do_put_in_core = True
do_evolve_planet = True

# Some output files.
f = open('logfile', 'w')
file_model = open('models.txt', 'w')

# Number of models to run (see loops below).
imax = 4
jmax = 5

#-------------------------------------------------------------------------------
# Prompt user before executing run.
#-------------------------------------------------------------------------------
print()
print("Sequential run for", imax * jmax, "models (imax*jmax).")
print("Warning: Will delete previous runs! (Will run script_clean.sh!)")
print("Warning: Will recompile code! (Will run mk!)")
print()
answer = input("Do you want to continue? [y/n]: ")
print()
if answer.lower().startswith("y"):
    print("Starting...")
    os.system("./script_clean.sh")
    print("Compiling...")
    print()
    os.system("./mk")
    print()
elif answer.lower() != "y":
    print("Aborted...")
    print()
    exit()

#-------------------------------------------------------------------------------
# Start of loop for sequential runs.
# Note that python does not include the upper bound of the range.
#-------------------------------------------------------------------------------
for i in range(1, imax + 1, 1):

    # Varying planetary mass.
    if imax == 1:
        mp = 1.0
        #mp = 8.0
    else:
        # Logarithmic.
        mp = math.pow(2.0, 3.0 * (i - 1) / (imax - 1))
        # Linear.
        #mp = 1.0 + (8.0 - 1.0) * (i - 1) / (imax - 1)

    # Varying core mass.
    #mcore = 0.0 + (20.0 - 0.0)*(i-1)/(imax-1)
    #mp = 1.0
    #mcore = mcore_array[i-1]

    # Planet mass without the core.
    # Normally, the core mass should be subtracted from the planet mass.
    # Warning: However, the final planet mass does not add up to mp_wo_core + mcore!
    #mp_wo_core = mp - mcore * mearth / mjup
    # Perhaps better to use just the planet mass without the core mass.
    mp_wo_core = mp

    for j in range(1, jmax + 1, 1):
        # Varying column depth (for depositing stellar radiation as heat).
        #irrad_col = 300.0 + (2300.0 - 300.0)*(i-1)/(imax-1)
        #Teq = 1250.0

        # Varying stellar effective temperature.
        #Teff = 7000.0 + (12000.0 - 7000.0)*(i-1)/(imax-1)

        # Varying planetary equilibrium temperature.
        # This is not physical, but rather a limiting case when no internal heat sources are present.
        if jmax == 1:
            #Teq = 0.0
            Teq = 1500.0
        else:
            #Teq = 0.0 + (2500.0 - 0.0) * (j - 1) / (jmax - 1)
            Teq = 0.0 + (2000.0 - 0.0) * (j - 1) / (jmax - 1)

        # Flux hitting planet's dayside.
        # As a function of stellar Teff.
        #flux_dayside = sigma*Teff**4 * (rs*rsun/orb_sep/au)**2
        # As a function of planetary Teq.
        flux_dayside = 4.0 * sigma * Teq**4

        # Calculated only for output purposes.
        # Alternatively: Teq = Teff*(rs*rsun/2.0/orb_sep/au)**0.5
        #Teq = (flux_dayside/4.0/sigma)**0.25

        # New version (same name, overwriting files for each model).
        createmodel = "planet_1_create.mod"
        coremodel = "planet_2_core.mod"
        evolvemodel = "planet_3_evolve.mod"

        # Old version (custom name for each model).
        #createmodel = "planet_create_" + str(mp_wo_core)[0:6] + "_MJ_" + str(rp) + "_RJ.mod"
        #coremodel = "planet_core_" + str(mp)[0:6] + "_MJ_" + str(mcore)[0:6] + "_ME_" + str(rp) + "_RJ.mod"
        #evolvemodel = "planet_evolve_" + str(mp)[0:6] + "_MJ_" + str(mcore)[0:6] + "_ME_" + str(rp) + "_RJ.mod"

        # Print to the screen or into a file.
        count = (i - 1) * jmax + j
        # Format the count variable to occupy two spaces.
        count_formatted = "{:2d}".format(count)
        # Redirect print from the screen into a file (with the file option at the end).
        #print("Model:",count_formatted,"Indices (i,j):",i,j,"Parameters (Mass,Mcore,Teq):",mp,mcore,Teq)
        print("Model:", count_formatted, "Indices (i,j):", i, j, "Parameters (Mass,Mcore,Teq):", mp, mcore, Teq, file=file_model)

        #print("Parameters: (Mass, Core, Teff, Teq:",mp,mcore,Teff,Teq,round(Teq,2))
        #print("Files:",createmodel,coremodel,evolvemodel)
        #print()
        #my.print_parameters(mp,rp,mcore,rhocore,mp_wo_core,irrad_col,flux_dayside,Teq,y,z,maxage)

        #-----------------------------------------------------------------------
        # Uncomment to skip the rest of the script. (Useful for a dry run.)
        #-----------------------------------------------------------------------
        # Uncommenting will cause an error since the folder out/ part is inside the loop.
        # But this error can be ignored (won't happen when it's commented out).
        #continue

        #-----------------------------------------------------------------------
        # Create planet.
        #-----------------------------------------------------------------------
        if do_create_planet:
            #inlist1 = "inlist_create_" + str(mp)[0:6] + "_MJ_" + str(mcore)[0:6] + "_ME_" + str(rp) + "_RJ"
            inlist1 = "inlist_create_" + str(mp)[0:6] + "_MJ_" + str(Teq)[0:6] + "_Teq"
            run_time = my.create_planet(mp_wo_core, rp, y, z, inlist1, createmodel)

        success = True

        # Why is all the rest needed? Not all seem necessary.
        if not os.path.exists(createmodel):
            #print("Heyoooooooooouuu!!! (createmodel)",os.path.exists(createmodel))
            success = False

        k = open('LOGS/history.data', 'r')
        for line in k.readlines():
            pass
        last_temp = line
        k.close()
        last = last_temp.split()
        #print("Final model number in create =",last[0])
        #print("Is last[0]==1000?",last[0]=="1000")
        if last[0] == "1000":
            success = False
        # Write into the logfile.
        outstring = '%6.3f\t%6.3f\t%6.3f\t%s\n' % (mp, rp, run_time, success)
        f.write(outstring)

        # Note: Continue skips the rest and goes directly to the next loop!
        #if not success:
            #print("Aboo...")
            #continue

        #-----------------------------------------------------------------------
        # Put in core.
        #-----------------------------------------------------------------------
        if do_put_in_core:
            if mcore > 0.0:
                #inlist2 = "inlist_core_" + str(mp)[0:6] + "_MJ_" + str(mcore)[0:6] + "_ME_" + str(rp) + "_RJ"
                inlist2 = "inlist_core_" + str(mp)[0:6] + "_MJ_" + str(Teq)[0:6] + "_Teq"
                run_time = my.put_core_in_planet(mcore, rhocore, inlist2, createmodel, coremodel)
            else:
                shutil.copyfile(createmodel, coremodel)

        #if not os.path.exists(coremodel):
            #print("Heyoooooooooouuu!!! (coremodel)",os.path.exists(coremodel))
            #continue

        # Delete unnecessary prints.
        # Skip the rest and go to the next for loop.
        #print("Skibbidy...")
        #continue
        #quit("Quitting...")

        #-----------------------------------------------------------------------
        # Evolve planet.
        #-----------------------------------------------------------------------
        if do_evolve_planet:
            #inlist3 = "inlist_evolve_" + str(mp)[0:6] + "_MJ_" + str(mcore)[0:6] + "_ME_" + str(rp) + "_RJ"
            inlist3 = "inlist_evolve_" + str(mp)[0:6] + "_MJ_" + str(Teq)[0:6] + "_Teq"
            run_time = my.evolve_planet(irrad_col, flux_dayside, maxage, inlist3, coremodel, evolvemodel)

        #if not os.path.exists(evolvemodel):
            #print("Heyoooooooooouuu!!! (evolvemodel)",os.path.exists(evolvemodel))
            #continue

        # Copy the LOGS folder (deleting old ones).
        if imax * jmax > 1:

            # Construct the directory path (including the parent directory).
            #dirpath = Path("LOGS_" + str(mp)[0:6] + "_MJ")
            dirpath = Path("out/LOGS_" + str(count))

            # Check if the parent directory exists, and create it if it doesn't.
            if not dirpath.parent.exists():
                dirpath.parent.mkdir(parents=True)

            # Remove the directory if it exists.
            if dirpath.exists() and dirpath.is_dir():
                shutil.rmtree(dirpath)

            # Copy the contents of "LOGS" directory to the destination directory
            shutil.copytree("LOGS", dirpath)

#-------------------------------------------------------------------------------
# End of main loop.
# Close the output files.
#-------------------------------------------------------------------------------
f.close()
file_model.close()

#-------------------------------------------------------------------------------
# Display the contents of the file models.txt which lists all the models.
#-------------------------------------------------------------------------------
filename = 'models.txt'

# Open the file in read mode.
with open(filename, 'r') as file:
    # Read the contents of the file.
    file_contents = file.read()
# The file is automatically closed after the block (no need to manually close it).

# Print the contents to the screen.
print()
print("This run contains the following models: (see models.txt)")
print(file_contents)

# Copy the file into the output directory (only for sequential runs).
if imax * jmax > 1:
    shutil.copy(filename, "out/")

#-------------------------------------------------------------------------------
# Done.
#-------------------------------------------------------------------------------
print("Python Script: All done...")
