#!/usr/bin/env python
import numpy as np
import os
import shutil
import time

#---------------------------------------------------------------------------
# Subroutines needed by the main python_run.py script.
#---------------------------------------------------------------------------

# Constants.
msun = 1.9884099e33
rsun = 6.9598e10
mjup = 1.8981246e30
rjup = 7.1492e9
mearth = 5.9721679e27
au = 1.496e13

#---------------------------------------------------------------------------
# Create the initial planet without the core.
#---------------------------------------------------------------------------
def create_planet(mp_wo_core,rp,y,z,inlist1,createmodel):
    # Keep track of run time (see further below).
    start_time = time.time()

    print()
    print("Create initial planet...")

    # Read and write inlist.
    f = open('inlist_create','r')
    g = f.read()
    f.close()

    g = g.replace("<<m>>",str(mp_wo_core*mjup))
    g = g.replace("<<r>>",str(rp*rjup))
    g = g.replace("<<z>>",str(z))
    g = g.replace("<<y>>",str(y))
    g = g.replace("<<ritefile>>",'"' + createmodel + '"')

    h = open(inlist1,'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist1,"inlist")

    # Run.
    # print("*** Create: Uncomment the execution!")
    os.system('./star')

    # Keep track of run time (see above).
    run_time = time.time() - start_time
    print("Run time for create_planets in seconds:",run_time)
    return run_time

#---------------------------------------------------------------------------
# Put in the core.
#---------------------------------------------------------------------------
def put_core_in_planet(mcore,rhocore,inlist2,createmodel,coremodel):
    start_time = time.time()

    print()
    print("Put in core...")

    # Read and write inlist.
    f = open('inlist_core','r')
    g = f.read()
    f.close()

    g = g.replace("<<new_core_mass>>",str(mcore*mearth/msun))
    g = g.replace("<<core_avg_rho>>",str(rhocore))
    g = g.replace("<<loadfile>>",'"' + createmodel + '"')
    g = g.replace("<<ritefile>>",'"' + coremodel + '"')

    h = open(inlist2,'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist2,"inlist")

    # Run.
    # print("*** Core: Uncomment the execution!")
    os.system('./star')

    run_time = time.time() - start_time
    print("Run time to put in core in seconds:",run_time)
    return run_time

#---------------------------------------------------------------------------
# Evolve to the desired starting age.
#---------------------------------------------------------------------------
def evolve_planet(irrad_col,flux_dayside,maxage,inlist3,coremodel,evolvemodel):
    start_time = time.time()

    print()
    print("Evolve planet...")

    # Read and write inlist.
    f = open('inlist_evolve','r')
    g = f.read()
    f.close()

    g = g.replace("<<irrad_col>>",str(irrad_col))
    g = g.replace("<<flux_dayside>>",str(flux_dayside))
    g = g.replace("<<maxage>>",str(maxage))
    g = g.replace("<<loadfile>>",'"' + coremodel + '"')
    g = g.replace("<<ritefile>>",'"' + evolvemodel + '"')

    h = open(inlist3,'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist3,"inlist")

    # Run.
    # print("*** Evolve: Uncomment the execution!")
    os.system('./star')

    run_time = time.time() - start_time
    print("Run time to evolve in seconds:",run_time)
    return run_time

#---------------------------------------------------------------------------
# Print stuff.
#---------------------------------------------------------------------------
def print_parameters(mp,rp,mcore,rhocore,mp_wo_core,irrad_col,flux_dayside,Teq,y,z,maxage):
    print('######################################################')
    print('Parameters:')
    print('mp/mj =',mp)
    print('rp/rj =',rp)
    print('mcore/me =',mcore)
    print('mcore/msun =',mcore*mearth/msun)
    print('rhocore/cgs =',rhocore)
    print('(mp-mcore)/mj =',mp_wo_core)
    print('irrad_col =',irrad_col)
    print('flux_dayside/1.e9 =',flux_dayside/1.e9)
    print('Teq =',Teq)
    print('Z =',z)
    print('Y =',y)
    print('evolve to age/Myr =',maxage/1.e6)
    print('######################################################')
    return
