#!/usr/bin/env python
from pathlib import Path
import math
import python_subs as my
import numpy as np
import os
import shutil
import sys


def PyRun(index,mass,mcore,rhocore,y,z,rpin,Teq,irrad_col,J0,pcmax,pcmin,s0,ohm_on,ohm_full,ak):

	# Common evolutionary parameters
	maxage = 1.1e10	  # Ending age [yr].
	# Flux dayside [erg cm^-2 s^-1].
	sigma_sb = 5.670374419e-5		# Stefan-Boltzmann constant [erg cm^-2 s^-1 K^-4].
	flux_dayside = 4.0 * sigma_sb * Teq**4

	#-------------------------------------------------------------------------------
	# Use the recommended indentation of 4 spaces (as per PEP 8).
	# Check scripts using the flake8 command (flake8 --ignore=E265,E501 filename).
	#-------------------------------------------------------------------------------

	# Function for cleaning all contents of a directory (both files and subfolders).
	def empty_directory(directory):
		for item in os.listdir(directory):
			item_path = os.path.join(directory, item)
			try:
				if os.path.isfile(item_path):
					os.remove(item_path)
				elif os.path.isdir(item_path):
					shutil.rmtree(item_path)
			except Exception as e:
				print(f"Failed to remove {item_path}. Reason: {e}")

	# Flags to skip steps.
	do_create_planet = True
	do_put_in_core = True
	do_evolve_planet = True

	#-------------------------------------------------------------------------------
	# Prompt user before executing run.
	#-------------------------------------------------------------------------------
	print("Starting...")
	os.system("./script_clean.sh")
	print("Compiling...")
	print()
	os.system("./mk")
	print()

	# Some additional output files.
	file_log = open('logfile_'+str(index)+'.txt', 'w')

    # Same name, overwriting files for each model.
	createmodel = "planet_1_create.mod"
	coremodel = "planet_2_core.mod"
	evolvemodel = "planet_3_evolve.mod"

    # Uncomment to skip the rest of the script. (Useful for a dry run.)
    # continue

    #-----------------------------------------------------------------------
    # Create planet.
    #-----------------------------------------------------------------------
	if do_create_planet:
		inlist1 = "inlist_create_" + str(index)
		run_time = my.create_planet(mass, rpin, y, z, inlist1, createmodel)

		success = True

		if not os.path.exists(createmodel):
			success = False

		k = open('LOGS/history.data', 'r')
		for line in k.readlines():
			pass
		last_temp = line
		k.close()
		last = last_temp.split()
		if last[0] == "1000":
			success = False
		# Write into the logfile.
		outstring = '%6.3f\t%s\n' % (run_time, success)
		file_log.write(outstring)

	    # Note: Continue skips the rest and goes directly to the next loop!
	    #if not success:
	        #print("Error (create planet)...")
	        #continue

	    # Put in core.
		if do_put_in_core:
			if mcore > 0.0:
				inlist2 = "inlist_core_" + str(index)
				run_time = my.put_core_in_planet(mcore, rhocore, inlist2, createmodel, coremodel)
			else:
				shutil.copyfile(createmodel, coremodel)

        # Skip the rest and go to the next for loop.
	    #continue

	    #-----------------------------------------------------------------------
	    # Evolve planet.
	    #-----------------------------------------------------------------------
		if do_evolve_planet:
			inlist3 = "inlist_evolve_" + str(index)

			# Call the evolution run with the input variables
			run_time = my.evolve_planet(irrad_col, flux_dayside, maxage, inlist3, coremodel, evolvemodel, J0, pcmin, pcmax, s0, ohm_on, ohm_full, ak)

				#continue


	#-------------------------------------------------------------------------------
	# Close the output files.
	#-------------------------------------------------------------------------------
	file_log.close()



if __name__=="__main__":

	# out_path = "/home/user/MESA/work_2.1/mp_v2.2/parameter_explor/"
	out_path = "/home/daniele/codes/MESA/version_r15140/make_planets_v2.2/Output_M07-03/"
	if not os.path.exists(out_path):
		os.mkdir(out_path)
		print("Creating the output directory: ",out_path)

	# Main planetary parameters.
	mass_arr = [0.7,0.3] # Planetary mass [Jupiter mass].
	y_arr = [0.24]		# Mass fraction of Helium of the planet envelope.
	z_arr = [0.02]		# Mass fraction of heavier elements of the planet envelope.
	rpin = 3.5	   	# Initial planet radius [Rj].

	# Irradiation parameters.
	Teq_arr = np.arange(1000,2251,1000) 	# Equilibrium temperature [K].
	Teq_arr = [0.,1500.,2000.]
	irrad_col_arr = [250.]					# Column depth for irradiation, Sigma* [cm^2 g^-1].

	# Core parameters
	mcore_arr = [0.,5.] 	# Core masses [Earth masses].
	rhocore_arr = [10.]	# Core density [g cm^-3].

	# Ohmic heating parameters
	J0_arr = np.linspace(4000,8000,4) # Amplitude of current, J = J0*(sigma[S/m])  J and J0 in [statamp cm^-2]
	J0_arr = [0]
	pcmax_arr = [5e4]    # Maximum pressure [bar] (pressure below which the Ohmic term is gradually forced to zero)
	pcmin_arr = [1e-8]  # Minimum pressure [bar] (pressure below which the Ohmic term is gradually forced to zero)
	s0_arr = [0.25]	     # Width of the tanh in log(p[bar]) of the gradual decrease to zero for the radial profile of Ohmic term
	ohm_on_arr = [1e5]	 # Age at which the Ohmic heat is turned on linearly (avoid early times)
	ohm_full_arr = [1e6] # Age at which the Ohmic heat reaches its full value
	ak_arr = [1e-7]		 # Potassium mass abundance (for the conductivity)

	file_model = open(out_path+'model_list.txt', 'w')
	file_crash = open(out_path+'model_crashed_list.txt', 'w')
	print("index      M[Mj]   Mc[Me] rhoc[gcc]  Rp[Rj] Teq[K] S*[cm2/g]      Y       Z       J0[cgs] pcmax-pcmin     s0  ohm_on/full[yr]    aK",file=file_model)
	print("index (see model_list.txt for the parameters)",file=file_crash)
	file_model.close()
	file_crash.close()

	index = 0

	for aa in range(len(pcmax_arr)):
		pcmax = pcmax_arr[aa]
		for bb in range(len(pcmin_arr)):
			pcmin = pcmin_arr[bb]
			for cc in range(len(y_arr)):
				y = y_arr[cc]
				for dd in range(len(z_arr)):
					z = z_arr[dd]
					for ee in range(len(mass_arr)):
						mass = mass_arr[ee]
						for ff in range(len(mcore_arr)):
							mcore = mcore_arr[ff]
							for mm in range(len(rhocore_arr)):
								rhocore = rhocore_arr[mm]
								for gg in range(len(Teq_arr)):
									Teq= Teq_arr[gg]
									for nn in range(len(irrad_col_arr)):
										irrad_col = irrad_col_arr[nn]
										for hh in range(len(J0_arr)):
											J0 = J0_arr[hh]
											for ii in range(len(s0_arr)):
												s0 = s0_arr[ii]
												for jj in range(len(ohm_on_arr)):
													ohm_on = ohm_on_arr[jj]
													for kk in range(len(ohm_full_arr)):
														ohm_full = ohm_full_arr[kk]
														for ll in range(len(ak_arr)):
															ak = ak_arr[ll]

															# Print in the list of models the input parameters
															index = index + 1
															file_model = open(out_path+'model_list.txt', 'a')
															print("{:5d}".format(index),"{:8.1f}".format(mass),"{:8d}".format(math.floor(mcore)),"{:8d}".format(math.floor(rhocore)),"{:8d}".format(math.floor(rpin)),"{:8d}".format(math.floor(Teq)),"{:8d}".format(math.floor(irrad_col)),"{:8.2f}".format(y),"{:8.2f}".format(z),"{:8d}".format(math.floor(J0)),"{:8.0e}".format(pcmax).replace("+0",""),"{:8.0e}".format(pcmin).replace("+0",""),"{:8.2f}".format(s0),"{:8.0e}".format(ohm_on).replace("+0",""),"{:8.0e}".format(ohm_full).replace("+0",""),"{:8.0e}".format(ak).replace("-0","-"),file=file_model)
															file_model.close()

															PyRun(index,mass,mcore,rhocore,y,z,rpin,Teq,irrad_col,J0,pcmax,pcmin,s0,ohm_on,ohm_full,ak)
																																								# Name of the run
															name = str(math.floor(index))

															print("="*20)
															print("Planet: index="+"{:3d}".format(index),"M="+"{:4.1f}".format(mass),"Mc="+"{:3d}".format(math.floor(mcore)),"rhoc="+"{:2d}".format(math.floor(rhocore)),"Rpin="+"{:3d}".format(math.floor(rpin))," Teq="+"{:4d}".format(math.floor(Teq)),"Sigma="+"{:4d}".format(math.floor(irrad_col)),"Y="+"{:4.2f}".format(y),"Z="+"{:4.2f}".format(z),"J0="+"{:4d}".format(math.floor(J0)),"pmax="+"{:.0e}".format(pcmax).replace("+0",""),"pmin="+"{:.0e}".format(pcmin).replace("+0",""),"s0="+"{:4.2f}".format(s0),"ohmon="+"{:.0e}".format(ohm_on).replace("+0",""),"ohmfull="+"{:.0e}".format(ohm_full).replace("+0","")," ak="+"{:.0e}".format(ak).replace("-0","-"))

															print("="*20)
															if os.path.isdir("LOGS") and os.path.isfile("LOGS/evolution.txt"):
																new_path = out_path+"LOGS_"+name
																shutil.move("LOGS",new_path)

																evo_file = out_path+"evolution_"+name+".txt"
																shutil.copy(new_path+"/evolution.txt",evo_file)

																if (os.path.exists("logfile_"+str(index)+".txt")):
																	shutil.move("logfile_"+str(index)+".txt",out_path+"logfile_"+str(index)+".txt")
																if (os.path.exists("inlist_create_"+str(index))):
																	shutil.move("inlist_create_"+str(index),out_path+"inlist_create_"+str(index))
																if (os.path.exists("inlist_core_"+str(index))):
																	shutil.move("inlist_core_"+str(index),out_path+"inlist_core_"+str(index))
																if (os.path.exists("inlist_evolve_"+str(index))):
																	shutil.move("inlist_evolve_"+str(index),out_path+"inlist_evolve_"+str(index))
																if (os.path.exists("in/input_ohmic.txt")):
																	shutil.move("in/input_ohmic.txt",out_path+"inlist_ohmic_"+str(index))

																print("="*20)
																print("The output folder moved to "+new_path)
															else:
																file_crash = open(out_path+'model_crashed_list.txt', 'a')
																print("Simulation crashes!\n")
																print("Model index="+"{:3d}".format(index),file=file_crash)
																file_crash.close()
	

	# Done.
	print("Python Script: All done!")
