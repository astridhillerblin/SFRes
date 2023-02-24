# Resonance contributions to inclusive electron-proton scattering observables.

Details about the physics and formalism behind the codes can be found in our works:
"Nucleon resonance contributions to unpolarised inclusive electron scattering", A. N. Hiller Blin et al., Phys. Rev. C 100 (2019) 035201, arXiv:1904.08016 [hep-ph]
"Resonant contributions to inclusive nucleon structure functions from exclusive meson electroproduction data", A. N. Hiller Blin et al., Phys. Rev. C 104 (2021) 025201, arXiv:2105.05834 [hep-ph]
"Nucleon resonance contributions to unpolarised inclusive electron scattering", A. N. Hiller Blin, V. I. Mokeev, and W. Melnitchouk, arXiv:2212.11952 [hep-ph].

## Fortran code
### Includes complete and speedy uncertainty propagation
The code can be run with either of the following two options:
1) gfortran main.f -fcoarray=single
2) ./compile.sh

For option 2), before being able to (permanently) run it, you will need to turn the file into an executable (one single time):
chmod 755 compile.sh

When running the code, the following steps happen in order:
1) The terminal prompts the user to enter a Q^2 value in GeV^2 at which output files should be generated. Please enter the value of choice.
2) The terminal asks whether truncated moments should be generated (y) or not (n). These observables are integrated over energies and given as functions of Q^2 in the range 1 to 3.5 GeV^2. Therefore, it is not necessary to re-generate these files every time, since they will remain the same for any choice in 1). Due to the integrations, the running time is slightly longer and thus it might be of advantage to skip this step most of the time.
3) The final question is whether the uncertainty propagation should be generated (y or n). This involves generating 10000 Monte-Carlo samples, for which reason a long runtime is expected. Therefore, while for final results one should generate the uncertainties, for tests one should consider skipping this step.

The following output files are then created in the folder Output:

A) sing_res_i.dat, where i is the index denoting a single resonance, running from 1 to 19. The order is:
1) N(1440)1/2+
2) N(1520)3/2-
3) N(1535)1/2-
4) N(1650)1/2-
5) N(1675)5/2-
6) N(1680)5/2+
7) N(1710)1/2+
8) N(1720)3/2+
9) Couplings currently set to 0, allows for new resonance in future
10) Couplings currently set to 0, allows for new resonance in future
11) Couplings currently set to 0, allows for new resonance in future
12) Delta(1232)3/2+
13) Delta(1620)1/2-
14) Delta(1700)3/2-
15) Couplings currently set to 0, allows for new resonance in future
16) Couplings currently set to 0, allows for new resonance in future
17) Couplings currently set to 0, allows for new resonance in future
18) Couplings currently set to 0, allows for new resonance in future
19) N'(1720)3/2+

B) sum_res.dat and sum_res_interf.dat, for the incoherent and coherent sum of resonances, respectively.

If the user chose the truncated moments to be generated, then:

C) Trunc_Reg1.dat, Trunc_Reg2.dat, Trunc_Reg3.dat, and Trunc_RegTot.dat are created, giving the truncated moments in the 1st, 2nd, and 3rd resonance regions, as well as over the entire resonance regime, respectively.

If the user chose to generate the uncertainty band propagation, then the following steps happen:
D) res_samp.dat (incoherent) and resinterf_samp.dat (coherent) are generated. They store all the MC-sampled observables for as many events as in the MC.

E) From the files generated in D), finally the mean and standard deviation values are computed for the observables and stored in res_meanstd.dat (incoherent) and resinterf_meanstd.dat (coherent).

Each generated file is a table of observables vs the proton-virtual-photon center-of-mass energy W (or vs Q^2 in the case of the truncated moments). Please refer to the README file in the Output/ folder for further explanation of these tables.

## Python code
### Simple visualization, but no uncertainty propagation

The Python folder contains equivalent Python (.py) and JupyterNotebook (.ipynb) codes that may be helpful to more quickly visualize observables - however, no uncertainty propagration is included, as it would be too slow for the needed MC generation. All observables are created from the central values of the electrocouplings.
