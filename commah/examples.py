import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import commah


def runcommand(cosmology='WMAP5'):
    """ Example interface commands """

    # Return the WMAP5 cosmology concentration predicted for
    # z=0 range of masses
    Mi = [1e8, 1e9, 1e10]
    zi = 0.
    print "Concentrations for haloes of mass ", Mi, " at z=", zi
    output = commah.run(cosmology=cosmology, zi=zi, Mi=Mi)

    print output['c'].flatten()

    # Return the WMAP5 cosmology concentration predicted for
    # z=0 range of masses AND cosmological parameters
    Mi = [1e8, 1e9, 1e10]
    zi = 0.
    print "Concentrations for haloes of mass ", Mi, " at z=", zi
    output, cosmo = commah.run(cosmology=cosmology, zi=zi, Mi=Mi,
                               retcosmo=True)

    print output['c'].flatten()
    print cosmo

    # Return the WMAP5 cosmology concentration predicted for MW
    # mass (2e12 Msol) across redshift
    Mi = 2e12
    z = [0., 0.5, 1., 1.5, 2., 2.5]
    output = commah.run(cosmology=cosmology, zi=0., Mi=Mi, z=z)
    for zval in z:
        print "M(z=0)=", Mi, " has c(z=", zval, ")= ",\
              output[output['z'] == zval]['c'].flatten()

    # Return the WMAP5 cosmology concentration predicted for MW
    # mass (2e12 Msol) across redshift
    Mi = 2e12
    zi = [0., 0.5, 1., 1.5, 2., 2.5]
    output = commah.run(cosmology=cosmology, zi=zi, Mi=Mi)
    for zval in zi:
        print "M(z=", zval, ")=", Mi, " has concentration ",\
              output[(output['zi'] == zval) &
                     (output['z'] == zval)]['c'].flatten()

    # Return the WMAP5 cosmology concentration and
    # rarity of high-z cluster
    Mi = 2e14
    zi = 6.
    output = commah.run(cosmology=cosmology, zi=zi, Mi=Mi)
    print "Concentrations for haloes of mass ", Mi, " at z=", zi
    print output['c'].flatten()
    print "Mass variance sigma of haloes of mass ", Mi, " at z=", zi
    print output['sig'].flatten()
    print "Fluctuation for haloes of mass ", Mi, " at z=", zi
    print output['nu'].flatten()

    # Return the WMAP5 cosmology accretion rate prediction
    # for haloes at range of redshift and mass
    Mi = [1e8, 1e9, 1e10]
    zi = [0.]
    z = [0., 0.5, 1., 1.5, 2., 2.5]
    output = commah.run(cosmology=cosmology, zi=zi, Mi=Mi, z=z)
    for Mval in Mi:
        print "dM/dt for halo of mass ", Mval, " at z=", zi,\
              " across redshift ", z, " is: "
        print output[output['Mi'] == Mval]['dMdt'].flatten()

    # Return the WMAP5 cosmology Halo Mass History for haloes with M(z=0) = 1e8
    M = [1e8]
    z = [0., 0.5, 1., 1.5, 2., 2.5]
    print "Halo Mass History for z=0 mass of ", M, " across z=", z
    output = commah.run(cosmology=cosmology, zi=0., Mi=M, z=z)
    print output['Mz'].flatten()

    # Return the WMAP5 cosmology formation redshifts for haloes at
    # range of redshift and mass
    M = [1e8, 1e9, 1e10]
    z = [0.]
    print "Formation Redshifts for haloes of mass ", M, " at z=", z
    output = commah.run(cosmology=cosmology, zi=0., Mi=M, z=z)
    for Mval in M:
        print output[output['Mi'] == Mval]['zf'].flatten()

    return "Done"


def plotcommand(cosmology='WMAP5', plotname=None):
    """ Example ways to interrogate the dataset and plot the commah output """

    # Plot the c-M relation as a functon of redshift
    xarray = 10.**(np.arange(1., 15., 0.2))
    yval = 'c'

    # Specify the redshift range
    zarray = np.arange(0., 5., 0.5)

    xtitle = r"Halo Mass (M$_{sol}$)"
    ytitle = r"Concentration"
    linelabel = "z="

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    plt.ylim([2., 30.])

    colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

    for zind, zval in enumerate(zarray):
        output = commah.run(cosmology=cosmology, zi=zval, Mi=xarray)

        # Access the column yval from the data file
        yarray = output[yval].flatten()

        # Plot each line in turn with different colour
        ax.plot(xarray, yarray, label=linelabel+str(zval), color=colors[zind])
        # Overplot the D08 predictions in black
        ax.plot(xarray, commah.cduffy(zval, xarray), color="black")

    ax.set_xscale('log')
    ax.set_yscale('log')

    leg = ax.legend(loc=1)
    # Make box totally transparent
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    for label in leg.get_texts():
        label.set_fontsize('small')  # the font size
    for label in leg.get_lines():
        label.set_linewidth(4)   # the legend line width

    if plotname:
        fig.tight_layout(pad=0.2)
        print "Plotting to ", plotname+"_CM_relation.png"
        fig.savefig(plotname+"_CM_relation.png", dpi=fig.dpi*5)
    else:
        plt.show()

    # Plot the c-z relation as a function of mass (so always Mz=M0)
    xarray = 10.**(np.arange(0., 1., 0.05)) - 1.
    yval = 'c'

    # Specify the mass range
    zarray = 10.**np.arange(6., 14., 2.)

    xtitle = r"Redshift"
    ytitle = r"NFW Concentration"
    linelabel = r"log$_{10}$ M$_{z}$(M$_{sol}$)="

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

    for zind, zval in enumerate(zarray):
        output = commah.run(cosmology=cosmology, zi=xarray, Mi=zval)

        # Access the column yval from the data file
        yarray = output[yval].flatten()

        # Plot each line in turn with different colours
        ax.plot(xarray, yarray,
                label=linelabel+"{0:.1f}".format(np.log10(zval)),
                color=colors[zind],)

    leg = ax.legend(loc=1)
    # Make box totally transparent
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    for label in leg.get_texts():
        label.set_fontsize('small')  # the font size
    for label in leg.get_lines():
        label.set_linewidth(4)   # the legend line width

    if plotname:
        fig.tight_layout(pad=0.2)
        print "Plotting to ", plotname+"_Cz_relation.png"
        fig.savefig(plotname+"_Cz_relation.png", dpi=fig.dpi*5)
    else:
        plt.show()

    # Plot the zf-z relation for different masses (so always Mz=M0)
    xarray = 10.**(np.arange(0., 1., 0.05)) - 1.
    yval = 'zf'

    # Specify the mass range
    zarray = 10.**np.arange(6., 14., 2.)

    xtitle = r"Redshift"
    ytitle = r"Formation Redshift"
    linelabel = r"log$_{10}$ M$_{z}$(M$_{sol}$)="

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

    for zind, zval in enumerate(zarray):
        output = commah.run(cosmology=cosmology, zi=xarray, Mi=zval)

        yarray = output[yval].flatten()

        # Plot each line in turn with different colour
        ax.plot(xarray, yarray,
                label=linelabel+"{0:.1f}".format(np.log10(zval)),
                color=colors[zind],)

    leg = ax.legend(loc=2)
    # Make box totally transparent
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    for label in leg.get_texts():
        label.set_fontsize('small')  # the font size
    for label in leg.get_lines():
        label.set_linewidth(4)   # the legend line width

    if plotname:
        fig.tight_layout(pad=0.2)
        print "Plotting to ", plotname+"_zfz_relation.png"
        fig.savefig(plotname+"_zfz_relation.png", dpi=fig.dpi*5)
    else:
        plt.show()

    # Plot the dM/dt-z relation for different masses (so always Mz=M0)
    xarray = 10.**(np.arange(0., 1., 0.05)) - 1.
    yval = 'dMdt'

    # Specify the mass range
    zarray = 10.**np.arange(10., 14., 0.5)

    xtitle = r"log$_{10}$ (1+z)"
    ytitle = r"log$_{10}$ Accretion Rate M$_{sol}$ yr$^{-1}$"
    linelabel = r"log$_{10}$ M$_z$(M$_{sol}$)="

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

    cosmo = commah.getcosmo(cosmology)
    for zind, zval in enumerate(zarray):
        output = commah.run(cosmology=cosmology, zi=xarray, Mi=zval,
                            com=False, mah=True)

        yarray = output[yval].flatten()

        # Plot each line in turn with different colour
        ax.plot(np.log10(xarray+1.), np.log10(yarray),
                label=linelabel+"{0:.1f}".format(np.log10(zval)),
                color=colors[zind],)

        # Plot the semi-analytic approximate formula from Correa et al 2015b
        semianalytic_approx = 71.6 * (zval / 1e12) * (cosmo['h'] / 0.7) *\
            (-0.24 + 0.75 * (xarray + 1.)) * np.sqrt(
            cosmo['omega_M_0'] * (xarray + 1.)**3. + cosmo['omega_lambda_0'])

        ax.plot(np.log10(xarray + 1.), np.log10(semianalytic_approx),
                color='black')

    leg = ax.legend(loc=2)
    # Make box totally transparent
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    for label in leg.get_texts():
        label.set_fontsize('small')  # the font size
    for label in leg.get_lines():
        label.set_linewidth(4)   # the legend line width

    if plotname:
        fig.tight_layout(pad=0.2)
        print "Plotting to ", plotname+"_dMdtz_relation.png"
        fig.savefig(plotname+"_dMdtz_relation.png", dpi=fig.dpi*5)
    else:
        plt.show()

    # Plot the dMdt-M relation as a function of redshift
    xarray = 10.**(np.arange(10., 14., 0.5))
    yval = 'dMdt'

    # Specify the redshift range
    zarray = np.arange(0., 5., 0.5)

    xtitle = r"Halo Mass M$_{sol}$"
    ytitle = r"Accretion Rate M$_{sol}$ yr$^{-1}$"
    linelabel = "z="

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

    for zind, zval in enumerate(zarray):
        output = commah.run(cosmology=cosmology, zi=zval, Mi=xarray,
                            com=False, mah=True)

        yarray = output[yval].flatten()

        # Plot each line in turn with different colour
        ax.plot(xarray, yarray, label=linelabel+str(zval),
                color=colors[zind],)

    ax.set_xscale('log')
    ax.set_yscale('log')

    leg = ax.legend(loc=2)
    # Make box totally transparent
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    for label in leg.get_texts():
        label.set_fontsize('small')  # the font size
    for label in leg.get_lines():
        label.set_linewidth(4)   # the legend line width

    if plotname:
        fig.tight_layout(pad=0.2)
        print "Plotting to ", plotname+"_MAH_M_relation.png"
        fig.savefig(plotname+"_MAH_M_relation.png", dpi=fig.dpi*5)
    else:
        plt.show()

    # Plot the (dM/M)dt-M relation as a function of redshift
    xarray = 10.**(np.arange(10., 14., 0.5))
    yval = 'dMdt'

    # Specify the redshift range
    zarray = np.arange(0., 5., 0.5)

    xtitle = r"Halo Mass M$_{sol}$"
    ytitle = r"Specific Accretion Rate yr$^{-1}$"
    linelabel = "z="

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

    for zind, zval in enumerate(zarray):
        output = commah.run(cosmology=cosmology, zi=zval, Mi=xarray,
                            mah=True, com=False)

        yarray = output[yval].flatten()

        # Plot each line in turn with different colour
        ax.plot(xarray, yarray/xarray, label=linelabel+str(zval),
                color=colors[zind],)

    ax.set_xscale('log')
    ax.set_yscale('log')

    leg = ax.legend(loc=1)
    # Make box totally transparent
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    for label in leg.get_texts():
        label.set_fontsize('small')  # the font size
    for label in leg.get_lines():
        label.set_linewidth(4)   # the legend line width

    if plotname:
        fig.tight_layout(pad=0.2)
        print "Plotting to ", plotname+"_specificMAH_M_relation.png"
        fig.savefig(plotname+"_specificMAH_M_relation.png", dpi=fig.dpi*5)
    else:
        plt.show()

    # Plot the Mz-z relation as a function of mass
    # (so mass is decreasing to zero as z-> inf)
    xarray = 10.**(np.arange(0., 1., 0.05)) - 1.
    yval = 'Mz'

    # Specify the mass range
    zarray = 10.**np.arange(10., 14., 0.5)

    xtitle = r"Redshift"
    ytitle = r"M(z) (M$_{sol}$)"
    linelabel = r"log$_{10}$ M$_{0}$(M$_{sol}$)="

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

    for zind, zval in enumerate(zarray):
        output = commah.run(cosmology=cosmology, zi=0., Mi=zval, z=xarray)

        yarray = output[yval].flatten()

        # Plot each line in turn with different colour
        ax.plot(xarray, yarray,
                label=linelabel+"{0:.1f}".format(np.log10(zval)),
                color=colors[zind],)

    ax.set_yscale('log')

    leg = ax.legend(loc=1)
    # Make box totally transparent
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    for label in leg.get_texts():
        label.set_fontsize('small')  # the font size
    for label in leg.get_lines():
        label.set_linewidth(4)   # the legend line width

    if plotname:
        fig.tight_layout(pad=0.2)
        print "Plotting to ", plotname+"_Mzz_relation.png"
        fig.savefig(plotname+"_Mzz_relation.png", dpi=fig.dpi*5)
    else:
        plt.show()

    # Plot the Mz/M0-z relation as a function of mass
    xarray = 10.**(np.arange(0., 1., 0.02)) - 1.
    yval = 'Mz'

    # Specify the mass range
    zarray = 10.**np.arange(10., 14., 0.5)

    xtitle = r"Redshift"
    ytitle = r"log$_{10}$ M(z)/M$_{0}$"
    linelabel = r"log$_{10}$ M$_{0}$(M$_{sol}$)="

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

    for zind, zval in enumerate(zarray):
        output = commah.run(cosmology=cosmology, zi=0., Mi=zval, z=xarray)

        yarray = output[yval].flatten()

        # Plot each line in turn with different colour
        ax.plot(xarray, np.log10(yarray/zval),
                label=linelabel+"{0:.1f}".format(np.log10(zval)),
                color=colors[zind],)

    leg = ax.legend(loc=3)
    # Make box totally transparent
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    for label in leg.get_texts():
        label.set_fontsize('small')  # the font size
    for label in leg.get_lines():
        label.set_linewidth(4)   # the legend line width

    if plotname:
        fig.tight_layout(pad=0.2)
        print "Plotting to ", plotname+"_MzM0z_relation.png"
        fig.savefig(plotname+"_MzM0z_relation.png", dpi=fig.dpi*5)
    else:
        plt.show()

    return "Done"
