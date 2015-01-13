import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm

from commah import *

def runcommand(cosmology='WMAP5'):
  """ Example interface commands """

  ## Return the WMAP5 cosmology concentration predicted for z=0 range of masses
  M = [1e8, 1e9, 1e10]
  z = 0.
  print "Concentrations for haloes of mass ", M, " at z=",z
  print loadval(cosmology = cosmology, z=z, M=M, val = 'c')

  ## Return the WMAP5 cosmology concentration predicted for MW mass (2e12 Msol) across redshift
  M = 2e12
  z = [0.,0.5,1.,1.5,2.,2.5]
  print "Concentrations for haloes of mass ", M, " at z=",z
  print loadval(cosmology = cosmology, z=z, M=M, val = 'c')

  ## Return the WMAP5 cosmology concentration and rarity of high-z cluster
  M = 2e14
  z = 6.
  print "Concentrations for haloes of mass ", M, " at z=",z
  print loadval(cosmology = cosmology, z=z, M=M, val = 'c')
  print "Fluctuation sigma of haloes of mass ", M, " at z=",z
  print loadval(cosmology = cosmology, z=z, M=M, val = 'sig0')
  print "Rarity for haloes of mass ", M, " at z=",z
  print loadval(cosmology = cosmology, z=z, M=M, val = 'nu')    

  ## Return the WMAP5 cosmology accretion rate prediction for haloes at range of redshift and mass
  M = [1e8, 1e9, 1e10]
  z = [0.,0.5,1.,1.5,2.,2.5]
  print "Concentrations for haloes of mass ", M, " at z=",z
  print loadval(cosmology = cosmology, z=z, M=M, val = 'dMdt')

  ## Return the WMAP5 cosmology Halo Mass History for haloes with M(z=0) = 1e8
  M = [1e8]
  z = [0.,0.5,1.,1.5,2.,2.5]
  print "Halo Mass History for z=0 mass of ", M, " across z=",z
  print loadval(cosmology = cosmology, z=z, M=M, val = 'Mz')

  ## Return the WMAP5 cosmology formation redshifts for haloes at range of redshift and mass
  M = [1e8, 1e9, 1e10]
  z = [0.]
  print "Formation Redshifts for haloes of mass ", M, " at z=",z
  print loadval(cosmology = cosmology, z=z, M=M, val = 'zf')

  return "Done"

def plotcommand(filename='Full_WMAP5_COM.pkl', plotout=None):
  """ Example ways to interrogate the dataset and plot the commah output """

  ## Check that the file type is indeed a pickle
  if filename.rfind('.pkl') >= 0:
    print filename + " is a pickle file"
    ## Load the file

    with open(filename, 'rb') as fin:
      data = pkl.load(fin)
      dataset = data.get("dataset",[])
      Mhalo = data.get("Mhalo",[])
      Redshift = data.get("Redshift",[])

  if filename.rfind('_COM') >= 0:
    com = True
    mah = True

  if filename.rfind('_MAH') >= 0:
    com = False
    mah = True

  ## Plot specifics
  plt.rc('text', usetex=True)

  ## Plot the Concentration Mass Relation
  if com:
  ## Plot the c-M relation as a function of redshift
    xval = 'M'
    xarray = 10.**(np.arange(1.,15.,1.))
    yval = 'c'

    ## Interpolate the output from commah
    interp = RectBivariateSpline(Mhalo, Redshift, dataset[yval])
 
    ## Specify the redshift range
    zarray = np.arange(0.,5.,0.5)

    xtitle = r"Halo Mass (h$^{-1}$ M$_{sol}$)"
    ytitle = r"Concentration"
    linelabel = "z="

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(xtitle)            
    ax.set_ylabel(ytitle)   
    colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

    for zind, zval in enumerate(zarray):
      ## Interpolate to the desired output values
      yarray = interp(xarray,zarray[zind]).flatten()

      ## Plot each line in turn with different colour   
      ax.plot(xarray, yarray, label=linelabel+str(zval), color=colors[zind],)
      ## Overplot the D08 predictions in black
      ax.plot(xarray, cduffy(zval, xarray), color="black")

    ax.set_xscale('log')
    ax.set_yscale('log')

    leg = ax.legend(loc=1)
    leg.get_frame().set_alpha(0) # this will make the box totally transparent
    leg.get_frame().set_edgecolor('white')
    for label in leg.get_texts():                
      label.set_fontsize('small') # the font size   
    for label in leg.get_lines():
      label.set_linewidth(4)  # the legend line width
  
    if plotout:
      fig.tight_layout(pad=0.2)
      print "Plotting to ",plotout+"_CM_relation.png"
      fig.savefig(plotout+"_CM_relation.png", dpi=fig.dpi*5) 
    else:
      plt.show()


  ## Plot the c-z relation as a function of mass
    xval = 'z'
    xarray = 10.**(np.arange(0.,1.,0.01)) - 1.
    yval = 'c'

    ## Interpolate the output from commah
    interp = RectBivariateSpline(Mhalo, Redshift, dataset[yval])
 
    ## Specify the mass range
    zarray = np.arange(6.,15.,2.)

    xtitle = r"Redshift"
    ytitle = r"NFW Concentration"
    linelabel = r"log M$_{0}$(h$^{-1}$ M$_{sol}$)="

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(xtitle)            
    ax.set_ylabel(ytitle)   
    colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

    for zind, zval in enumerate(zarray):
      ## Interpolate to the desired output values
      yarray = interp(10.**zval,xarray).flatten()

      ## Plot each line in turn with different colour   
      ax.plot(xarray, yarray, label=linelabel+"{0:.1f}".format( zval ), color=colors[zind],)

    leg = ax.legend(loc=1)
    leg.get_frame().set_alpha(0) # this will make the box totally transparent
    leg.get_frame().set_edgecolor('white')
    for label in leg.get_texts():                
      label.set_fontsize('small') # the font size   
    for label in leg.get_lines():
      label.set_linewidth(4)  # the legend line width
  
    if plotout:
      fig.tight_layout(pad=0.2)
      print "Plotting to ",plotout+"_Cz_relation.png"
      fig.savefig(plotout+"_Cz_relation.png", dpi=fig.dpi*5) 
    else:
      plt.show()

  ## Plot the zf-z relation as a function of redshift
    xval = 'z'
    xarray = 10.**(np.arange(0.,1.,0.01)) - 1.
    yval = 'zf'

    ## Interpolate the output from commah
    interp = RectBivariateSpline(Mhalo, Redshift, dataset[yval])
 
    ## Specify the mass range
    zarray = np.arange(6.,15.,2.)

    xtitle = r"Redshift"
    ytitle = r"Formation Redshift"
    linelabel = r"log M$_0$(h$^{-1}$ M$_{sol}$)="

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(xtitle)            
    ax.set_ylabel(ytitle)   
    colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

    for zind, zval in enumerate(zarray):
      ## Interpolate to the desired output values
      yarray = interp(10.**zval,xarray).flatten()

      ## Plot each line in turn with different colour   
      ax.plot(xarray, yarray, label=linelabel+"{0:.1f}".format( zval ), color=colors[zind],)

    leg = ax.legend(loc=2)
    leg.get_frame().set_alpha(0) # this will make the box totally transparent
    leg.get_frame().set_edgecolor('white')
    for label in leg.get_texts():                
      label.set_fontsize('small') # the font size   
    for label in leg.get_lines():
      label.set_linewidth(4)  # the legend line width
  
    if plotout:
      fig.tight_layout(pad=0.2)
      print "Plotting to ",plotout+"_zfz_relation.png"
      fig.savefig(plotout+"_zfz_relation.png", dpi=fig.dpi*5) 
    else:
      plt.show()      

## Plot the dMdt-M relation as a function of redshift
  xval = 'M'
  xarray = 10.**(np.arange(1.,15.,1.))
  yval = 'dMdt'

  ## Interpolate the output from commah
  interp = RectBivariateSpline(Mhalo, Redshift, dataset[yval])

  ## Specify the redshift range
  zarray = np.arange(0.,5.,0.5)

  xtitle = r"Halo Mass h$^{-1}$ M$_{sol}$"
  ytitle = r"Accretion Rate h$^{-1}$ M$_{sol}$ Gyr$^{-1}$"
  linelabel = "z="

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_xlabel(xtitle)            
  ax.set_ylabel(ytitle)   
  colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

  for zind, zval in enumerate(zarray):
    ## Interpolate to the desired output values
    yarray = interp(xarray,zval).flatten()

    ## Plot each line in turn with different colour   
    ax.plot(xarray, yarray, label=linelabel+str(zval), color=colors[zind],)

  ax.set_xscale('log')
  ax.set_yscale('log')

  leg = ax.legend(loc=2)
  leg.get_frame().set_alpha(0) # this will make the box totally transparent
  leg.get_frame().set_edgecolor('white')
  for label in leg.get_texts():                
    label.set_fontsize('small') # the font size   
  for label in leg.get_lines():
    label.set_linewidth(4)  # the legend line width

  if plotout:
    fig.tight_layout(pad=0.2)
    print "Plotting to ",plotout+"_MAH_M_relation.png"
    fig.savefig(plotout+"_MAH_M_relation.png", dpi=fig.dpi*5) 
  else:
    plt.show()


## Plot the (dM/M)dt-M relation as a function of redshift
  xval = 'M'
  xarray = 10.**(np.arange(1.,15.,1.))
  yval = 'dMdt'

  ## Interpolate the output from commah
  interp = RectBivariateSpline(Mhalo, Redshift, dataset[yval])

  ## Specify the redshift range
  zarray = np.arange(0.,5.,0.5)

  xtitle = r"Halo Mass h$^{-1}$ M$_{sol}$"
  ytitle = r"Specific Accretion Rate Gyr$^{-1}$"
  linelabel = "z="

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_xlabel(xtitle)            
  ax.set_ylabel(ytitle)   
  colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

  for zind, zval in enumerate(zarray):
    ## Interpolate to the desired output values
    yarray = interp(xarray,zval).flatten()/xarray

    ## Plot each line in turn with different colour   
    ax.plot(xarray, yarray, label=linelabel+str(zval), color=colors[zind],)

  ax.set_xscale('log')
  ax.set_yscale('log')

  leg = ax.legend(loc=1)
  leg.get_frame().set_alpha(0) # this will make the box totally transparent
  leg.get_frame().set_edgecolor('white')
  for label in leg.get_texts():                
    label.set_fontsize('small') # the font size   
  for label in leg.get_lines():
    label.set_linewidth(4)  # the legend line width

  if plotout:
    fig.tight_layout(pad=0.2)
    print "Plotting to ",plotout+"_MAH_M_relation.png"
    fig.savefig(plotout+"_MAH_M_relation.png", dpi=fig.dpi*5) 
  else:
    plt.show()



## Plot the Mz-z relation as a function of mass
  xval = 'z'
  xarray = 10.**(np.arange(0.,1.,0.01)) - 1.
  yval = 'Mz'

  ## Interpolate the output from commah
  interp = RectBivariateSpline(Mhalo, Redshift, dataset[yval])

  ## Specify the mass range
  zarray = np.arange(6.,15.,2.)

  xtitle = r"Redshift"
  ytitle = r"M(z) (h$^{-1}$ M$_{sol}$)"
  linelabel = r"log M$_{0}$(h$^{-1}$ M$_{sol}$)="

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_xlabel(xtitle)            
  ax.set_ylabel(ytitle)   
  colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

  for zind, zval in enumerate(zarray):
    ## Interpolate to the desired output values
    yarray = interp(10.**zval,xarray).flatten()

    ## Plot each line in turn with different colour   
    ax.plot(xarray, yarray, label=linelabel+"{0:.1f}".format( zval ), color=colors[zind],)

  ax.set_yscale('log')

  leg = ax.legend(loc=1)
  leg.get_frame().set_alpha(0) # this will make the box totally transparent
  leg.get_frame().set_edgecolor('white')
  for label in leg.get_texts():                
    label.set_fontsize('small') # the font size   
  for label in leg.get_lines():
    label.set_linewidth(4)  # the legend line width

  if plotout:
    fig.tight_layout(pad=0.2)
    print "Plotting to ",plotout+"_Mzz_relation.png"
    fig.savefig(plotout+"_Mzz_relation.png", dpi=fig.dpi*5) 
  else:
    plt.show()


## Plot the Mz/M0-z relation as a function of mass
  xval = 'z'
  xarray = 10.**(np.arange(0.,1.,0.01)) - 1.
  yval = 'Mz'

  ## Interpolate the output from commah
  interp = RectBivariateSpline(Mhalo, Redshift, dataset[yval])

  ## Specify the mass range
  zarray = np.arange(6.,15.,2.)

  xtitle = r"Redshift"
  ytitle = r"M(z)/M$_{0}$"
  linelabel = r"log M$_{0}$(h$^{-1}$ M$_{sol}$)="

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_xlabel(xtitle)            
  ax.set_ylabel(ytitle)   
  colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

  for zind, zval in enumerate(zarray):
    ## Interpolate to the desired output values
    yarray = interp(10.**zval,xarray).flatten()

    ## Plot each line in turn with different colour   
    ax.plot(xarray, yarray/10.**zval, label=linelabel+"{0:.1f}".format( zval ), color=colors[zind],)

  ax.set_xscale('log')
  ax.set_yscale('log')

  leg = ax.legend(loc=3)
  leg.get_frame().set_alpha(0) # this will make the box totally transparent
  leg.get_frame().set_edgecolor('white')
  for label in leg.get_texts():                
    label.set_fontsize('small') # the font size   
  for label in leg.get_lines():
    label.set_linewidth(4)  # the legend line width

  if plotout:
    fig.tight_layout(pad=0.2)
    print "Plotting to ",plotout+"_MzM0z_relation.png"
    fig.savefig(plotout+"_MzM0z_relation.png", dpi=fig.dpi*5) 
  else:
    plt.show()

  return "Done"