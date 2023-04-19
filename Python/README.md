
# Python code for simple visualization -- no uncertainty propagation

The Python (.py) and JupyterNotebook (.ipynb) codes are equivalent, they are just there to allow you to choose your prefered format.

## Running the Python code
Simply run ipython main.py on the terminal.

When running the code, the following steps happen in order:
1) The terminal prompts the user to enter a Q^2 value in GeV^2 at which the output should be generated. Please enter the value of choice.
2) The terminal prompts the user to enter a W value in GeV at which the output should be generated. Please enter the value of choice.
3) The terminal asks whether interference effects should be included: choose 1 or 0 for yes or no.
4) The terminal prompts the user to choose which of the 12 resonances to include.
If return is pressed (or the wrong input format is inserted), then the default 1,1,1,1,1,1,1,1,1,1,1,1 is taken: all resonances are included.
Otherwise, if 12 comma-separated 1s and 0s are inserted, the respective resonances are switched on and off. The order of resonances is:
N(1440)1/2+,N(1520)3/2-,N(1535)1/2-,N(1650)1/2-,N(1675)5/2-,N(1680)5/2+,N(1710)1/2+,N(1720)3/2+,Delta(1232)3/2+,Delta(1620)1/2-,Delta(1700)3/2-,N'(1720)3/2+

The output is a list of all observables at the chosen values of $Q^2$ and $W$:
$F_1$, $F_2$, $F_L$, $g_1$, $g_2$, $H_{1/2}$, $H_{3/2}$, $A_1$, $A_2$, $\sigma_T~[\mu b]$, $\sigma_L~[\mu b]$, $d\sigma/(dQ^2dW)~[nb/GeV^3]$
