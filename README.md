# Resonance contributions to inclusive electron-proton scattering observables.

You can run with either of the following two options:
1) gfortran main.f -fcoarray=single
2) ./compile.sh

For option 2), before being able to (permanently) run it, you will need to turn the file into an executable (one single time):
chmod 755 compile.sh

When running main.f, the following output files are created in the folder Output:
1) sing_res_i.dat, where i runs from 1 to 19, are tables of the contributions of each of the separate 19 resonances to inclusive observables.
2) sing_res.dat and sing_res_interf.dat are tables for the incoherent and coherent sum of resonance contributions to the inclusive observables, respectively.
