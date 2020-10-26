# Modified-Gravity

The codes available here can be used to compute some theoretical quantities in modified gravity.

The folder contains the following codes:

-MG_Power_Spectrum.c: This code computes the linear matter power spectrum in general relativity, f(R), and symmetron theories;

-params_mg.dat: This is an example of the input file for the code above;

-Spherical_Collapse_phy_dc.c: This code computes delta_c (the threshold for the halo spherical collapse) for a generic initial profile in general relativity and f(R);

-Spherical_Collapse_phy_dct.c: This code computes delta_t (the overdensity at the turnaround) for a generic initial profile in general relativity and f(R).

The complete list of input options for each code is given by the code when started without any input.

The code MG_Power_Spectrum.c needs the CAMB to work and the codes Spherical_Collapse_phy_dc.c and Spherical_Collapse_phy_dt.c need the library GSL.

If you use the code MG_Power_Spectrum.c, please cite the paper: https://arxiv.org/abs/1609.02544.

If you use the codes Spherical_Collapse_phy_dc.c and/or Spherical_Collapse_phy_dt.c, please cite at least one of the papers: https://arxiv.org/abs/1805.09918 and https://arxiv.org/abs/1809.10321.

Please, contact me if you have any question: rodrigo.voivodic@usp.br.
