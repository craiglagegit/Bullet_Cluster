# Bullet_Cluster
Code used to simulate the Bullet Cluster collision
Craig Lage 15-Feb-2018

UC Davis

************************************************************************************

This is a brief description of the code I used to generate the simulations described in  Lage, C. and Farrar, G. "Constrained Simulation of the Bullet Cluster". Astrophysical Journal, V787, P144, 2014. Available on the arXiv at http://arxiv.org/abs/1312.0959v1, and included in the docs directory.

First, note that I did not intend for others to use the code as it stands, so it is not very readable, but I will do my best to describe how it works.  With some study, you should be able to reconstruct what I have done from these notes.

(1) The code directory contains the code used, and the directory labeled bsearch_matrixn110 contains the best fitting simulation, with the actual enzo files in bsearch_matrixn110/batchfiles/batch3/ddfiles/run4.

(2) The overall running of the various pieces of code was managed by the bsearch_matrixn110/script.py file.  This is a fairly complex piece of code that was designed to run many copies of enzo in parallel with different initial conditions, then evaluate the fit between the simulation and the observations, then launch a new round of simulations with modified initial conditions.  I did this many times to try and find the best fit initial conditions.  If you are trying to run these simulations, you probably won't want to use the script.py code, but by studying it you can see the sequence that is run to perform a complete simulation.  The sequence is roughly as follows:

     (A) The garfields code was run to generate a Gaussian random magnetic field.  Typically this was only run once, then the field was scaled later based on the gas density (see step H).  It only needed to be re-run if a change in the field coherence length was desired. The code was controlled using the code/garfields/garFields_par file, and it outputs a file called magField.dat which you should find in bsearch_matrixn110.  Its use will be described later.

     (B) The smile program was then run to generate the DM halos.  I modified the code slightly, but you should find the code I ran in code/TriaxialHalo_Smile203.  The running is managed by the files smile.ini, smileb.script, smilebbatch, and smilebscript.csh, and the corresponding files with smileb replaced by smilem.  Throughout this, I use the letter "b" to refer to the smaller or "bullet" cluster, and the letter "m" to refer to the larger or "main" cluster.  These steps generates the files nbodym.dat and nbodyb.dat, which are text files giving the DM particle positions and velocities.

     (C) The next several steps are controlled by the transcript.csh script, as follows.

     (D) The code/TriaxialHalo_SGas8GTemp/TriaxialHalo is custom code that I wrote to decorate the DM halos with gas in a way which is in hydrostatic equilibrium.  The running of the code is managed by the triaxialb.cfg file.  It takes the nbodyb.dat, potentialb.dat, and triaxialb.cfg as inputs, and generates the bsearch_matrixn110/batchfiles/batch3/ddfiles/run4/bullet.dat file, which is a binary file containing the DM and gas particles in Gadget format. Of course, the same is done with m instead of b, generating main.dat.

     (E) The code code/CombineGalaxies_NewIP/CombineGalaxies is then run, taking the main.dat and bullet.dat files as inputs and generating a combined Gadget file collision.dat with the two halos in the proper initial locations and with teh proper initial velocities.

     (F) The small code code/pysubs/mod_grid_locations_big.py is run to modify the bsearch_matrixn110/batchfiles/batch3/ddfiles/run4/AMRTest.enzo file to center the enzo initial grids on the densest parts of the two clusters.  The adaptive mesh refinement will update this as the code runs, but it is good to start with the best possible grid locations.  Also, the AMRTest.enzo file controls the Gadget to Enzo translation, which is upcoming.

     (G) The small code code/pysubs/make_cool_tfudge.py is run to generate the bsearch_matrixn110/batchfiles/batch3/ddfiles/run4/cool_rates.infile, which adjust the radiative cooling based on the gas metallicity.  This is one of the modifications that I made to the Enzo code.

     (H) The code code/pysubs/translate_gadget_to_enzo_DEF_big.py then takes the files collision.dat, AMRTest.enzo, magField.dat as input and generates the binary Enzo initial conditions files.  These are in bsearch_matrixn110/batchfiles/batch3/ddfiles/run4 and have names like Grid*, Particle*, InternalEnergy* and TotalEnergy*.

     (I) Finally the Enzo code itself can be run.  The code is in code/enzo-2.0/src/enzo_MHD_Visc14 and is controlled by the bsearch_matrixn110/batchfiles/batch3/ddfiles/run4/AMRTest.enzo file. I have left the output directory fro the t=0 time step in the directory, whwere it is called DD0000.  A typical run generates about 100 of these directories, so about 100GB of output, which was too big to load up to github.

     (J) I typically evaluated the output using the code yt.  Using this code is a whole lesson in itself, so I won't try to describe it here, but I have left in a few of the simple yt scripts.

I hope this helps.  If you are trying to come up to speed on these simulations, the easiest way might be to work backwards and try running Enzo with the initial condition files I've provided, then try generating your own initial conditions files using the scripts I've written as starting points.

Don't hesitate to contact me at cslage@ucdavis.edu if I can be of assistance.


***ADDENDUM***  After I wrote the above and tried to push the files to github, I found that the following files exceeded github's 100MB file size limit, so I had to delete them.  You should be able to recreate these files by following the steps I outlined above, but if you need these files, contact me and we'll find a way to get them to you.

nbodym.dat
main.dat
collision.dat
The entire directory DD0000

***SECOND ADDENDUM***  After I completed the two papers in the docs directory, and before I completed my PhD, I wrote the code in the code/halo_builder directory.  This code dispenses with many of the steps in the complicated sequence outlined above.  You still need to run smile to generate the nbody.dat files, but then the halo_builder code does all of the remaining steps, including translating the input halos, adding the gas, and generating the magnetic field. It then outputs the Enzo initial condition files directly, without the intermediate step of generating Gadget files.  It is also capable of combining more than two clusters.  The running of the code is controlled by the .cfg files in code/halo_builder.  If you are embarking on running a new set of simulations using Enzo, you might consider using this code as a starting point.
