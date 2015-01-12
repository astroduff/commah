# Python script to create mass accretion rates and NFW concentrations for haloes of any mass at and redshift for any flat LCDM cosmology based on Correa et al 2015 a,b,c
# This process can take a few seconds therefore a range of standard cosmologies have been calculated (using "creategrid" to loop over "run") with a simple "loadval" routine
# enabling fast access by interpolation over the pickle files. Example commands are in "runexample" and plot commands in "plotexample".
# If user wants a cosmology outside of the ones provided (WMAP1, WMAP3, WMAP5, WMAP7, WMAP9, Planck13) then specify a dictionary cosmology and pass to creategrid
# then name this output file something memorable that will then be access in future by the loadval routine.
