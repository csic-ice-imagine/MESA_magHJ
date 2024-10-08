
Hi,

Here is my standard setup. I've attached an example of a Mtot=1MJ,
Mcore=10Mearth planet, which is evolved for 10Gyr,
and is irradiated by mesa's surface irradiation.

In order,

1) inlist_create_1.0_MJ_10.0_ME_2.0_RJ

makes an initial model with R=2RJ using create_initial_model.  things seem
to fail less often if you set initial_model_relax_num_steps = 0.

2) inlist_core_1.0_MJ_10.0_ME_2.0_RJ

put in the core of 10Mearth, rho=10g/cc and eps_core=0

3) inlist_evolve_1.0_MJ_10.0_ME_2.0_RJ

evolve for 10Gyr (or Teff_lower_limit = 100.0, for un-irradiated models).
irradiation put in using simple method by
        column_depth_for_irradiation = 300.0 ! cm^2/g
        irradiation_flux = 555501654.562 ! erg/cm^2/s ! day side flux!

You're probably not interested, but I've also included a python script
runscript.py (with functions in mysubprograms.py) that I use to generate
a grid of models from template inlists inlist_create, inlist_core,
inlist_evolve. Presently this script is set up to generate a grid of models
with different masses and the same starting radius (as well as other parameters).
the script does search and replace on template inlists.

Cheers,

Phil


