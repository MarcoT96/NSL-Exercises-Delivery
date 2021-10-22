## Exercise 04

<br>
The folder contains four subfolders related to the four points of the exercise.<br>

In each of these folders there is a makefile: if deemed appropriate download the C++ files
and code the `make` command to generate the executable <strong>'MolDyn_NVE.exe'</strong>; otherwise,
folders already contain simulation results saved to files typically in .dat format.<br>
In turn, each of these sub-folders contains three more, related to the three phases of the simulated Argon system: *solid*, *liquid* and *gas*. <br>
In particular, I would like to highlight the following issues:

- folder `04.2` contains only the code updated with the implementation of the **blocking method**, I did not do any simulation here.
- folder `04.3` contains an additional fourth sub-folder, called `dt_reduced`, in which I saved the results obtained by decreasing by a factor of 10 the MD integration time step *dt*; this reduction has led to heavier simulations and consequently files too large to be uploaded to `GitHub`. In particular, `dt_reduced/solid`, `dt_reduced/liquid` and `dt_reduced/gas` lack files related to the instantaneous values of the observables (output_epot.dat, output_ekin.dat, output_etot.dat and output_temp.dat). If you need these files for any reason, please contact me. However, all C++ codes are able to generate these files.
- folder `Optional Exercise` contains the data that I will use for the comparison between MD and MC Argon simulations, and include the calculation of pressure and g(r) in the best equilibrium conditions that I managed to obtain for the three phases.

The comment on the results is contained in the jupyter-notebook named <strong>'Numerical Exercises 4.ipynb'</strong>.<br>
This notebook does not appear in full form when viewed through `GitHub`: download the jupyter file and
open it locally to see <em>colors</em>, <em>titles</em>, <em>equations</em> and <em>graphs</em> correctly.
