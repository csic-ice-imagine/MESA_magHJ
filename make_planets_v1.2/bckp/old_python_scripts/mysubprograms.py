#!/usr/bin/env python

import numpy as np
import os
import shutil
import time

msun = 1.9892e33
rsun = 6.9598e10
mjup = 1.8986e30
rjup = 6.9911e9
mearth = 5.97e27
sigma=5.67e-5
au = 1.496e13

# make the initial planet without the core

def create_planet(mp_wo_core,rp,y,z,inlist1,createmodel):
    start_time = time.time()
    print("create initial planet")
    f = open('inlist_create','r')
    g = f.read()
    f.close()
    g = g.replace("<<m>>",str(mp_wo_core*mjup))
    g = g.replace("<<r>>",str(rp*rjup))
    g = g.replace("<<z>>",str(z))
    g = g.replace("<<y>>",str(y))
    g = g.replace("<<smwtfname>>", '"' + createmodel + '"')
    h = open(inlist1,'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist1,"inlist")
    os.system('./star')
    run_time = time.time() - start_time
    print("run time for create_planets in sec=",run_time)
    return run_time


# put in the core if me > 0.0

def put_core_in_planet(mcore,rhocore,inlist2,createmodel,coremodel):
    start_time = time.time()
    print("put core in planet")
    f = open('inlist_core','r')
    g = f.read()
    f.close()
    g = g.replace("<<loadfile>>",'"' + createmodel + '"')
    g = g.replace("<<smwtfname>>", '"' + coremodel + '"')
    g = g.replace("<<new_core_mass>>",str(mcore*mearth/msun))
    g = g.replace("<<core_avg_rho>>",str(rhocore))
    h = open(inlist2,'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist2,"inlist")
    os.system('./star')
    run_time = time.time() - start_time
    print("run time to put in core in sec=",run_time)
    return run_time


# evolve to the desired starting age for binary_rlo

def evolve_planet(irrad_col,flux_dayside,maxage,inlist3,coremodel,evolvemodel):
    start_time = time.time()
    print("evolve planet")
    f = open('inlist_evolve','r')
    g = f.read()
    f.close()
    g = g.replace("<<loadfile>>",'"' + coremodel + '"')
    g = g.replace("<<smwtfname>>", '"' + evolvemodel + '"')
    g = g.replace("<<irrad_col>>", str(irrad_col) )
    g = g.replace("<<flux_dayside>>", str(flux_dayside) )
    g = g.replace("<<maxage>>",str(maxage))
    h = open(inlist3,'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist3,"inlist")
    os.system('./star')
    run_time = time.time() - start_time
    print("run time to evolve in sec=",run_time)
    return run_time


# print stuff
def print_parameters(mp,rp,mcore,rhocore,mp_wo_core,irrad_col,flux_dayside,Teq,y,z,maxage):
    print('######################################################')
    print('parameters:')
    print('mp/mj=',mp)
    print('rp/rj=',rp)
    print('mcore/me=',mcore)
    print('mcore/msun=',mcore*mearth/msun)
    print('rhocore/cgs=',rhocore)
    print('(mp-mcore)/mj=',mp_wo_core)
    print('irrad_col=',irrad_col)
    print('flux_dayside/1.e9=',flux_dayside/1.e9)
    print('Teq=',Teq)
    print('z=',z)
    print('y=',y)
    print('evolve to age/Myr=',maxage/1.e6)
    print('######################################################')
    return
