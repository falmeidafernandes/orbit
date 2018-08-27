#!/ust/bin/env

################################################################################
__author__  = "Felipe Almeida-Fernandes"
__emais__   = ["felipe11@astro.ufrj.br", "felipefer42@gmail.com"]
__version__ = 1
__date__    = "2018-08-16"

#\todo expand this docstring
"""
The goal of this script is to perform orbital integrations using galpy.
"""

# Import the necessary modules #################################################
import_error_msg = ("could not import library {0}. You may try to install it " 
                    "using:\n$ pip install {0}")

import sys, os, warnings, argparse
import collections

# Galpy is the module that deals with:
# (i) the convertion between observables and the quantities used for integration
# (ii) the orbital integration routines
# (iii) the potential in which the orbits are integrated
try:
    import galpy
except ImportError:
    raise ImportError(import_error_msg.format("galpy"))
else:
    from galpy.orbit import Orbit
    from galpy.potential import MWPotential2014 as MW14

# astropy is used for:
# (i) especification of units for galpy 
# (ii) conversion between coordinate
try:
    import astropy
except ImportError:
    raise ImportError(import_error_msg.format("astropy"))
else:
    from astropy import units as u
    from astropy.coordinates import SkyCoord

# numpy is used for arrays and load/save files handling
try:
    import numpy as np
except ImportError:
    raise ImportError(import_error_msg.format("numpy"))

# pandas is used for data frames and load/save files handling
try:
    import pandas as pd
except ImportError:
    raise ImportError(import_error_msg.format("pandas"))

# matplotlib is used to produce the plots
bool_plot = None
try:
    from matplotlib import pyplot as plt
except ImportError:
    warnings.warn(import_error_msg.format("matplotlib")+
                  "\nThe code will continue but no plots will be made")
    bool_plot = False

def ensure_dir(savedir):
    if not os.path.exists(savedir):
        os.makedirs(savedir)
################################################################################

# PARSE ########################################################################

parser = argparse.ArgumentParser(description='Integrate galactic orbit')
parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="count")

# PARSE OBSERVABLES
parser.add_argument("--ra", help="right ascention in degrees",
                    type=float)

parser.add_argument("--dec", help="declination in degrees",
                    type=float)

parser.add_argument("--d", help="distance in kpc",
                    type=float)

parser.add_argument("--pmra", help="right ascention proper motion in mas/yr",
                    type=float)

parser.add_argument("--pmde", help="declination proper motion in mas/yr",
                    type=float)

parser.add_argument("--rv", help="radial (line of sight) velocity in km/s",
                    type=float)

# PARSE INTEGRATOR OPTIONS
parser.add_argument("-t", "--time", help="Integration time in Gyr",
                    type=float, default = 1)

parser.add_argument("-N", "--Nsteps", help="Number of steps to integrate",
                    type=int, default=10000)

parser.add_argument("--ro", help="Solar galactocentric distance in kpc",
                    type=float, default = 8.)

parser.add_argument("--vo", help="LSR velocity in km/s",
                    type=float, default = 220.)

parser.add_argument("--solar", help="Solar motion",
                    default='schoenrich')

parser.add_argument("--pot", help="Potential used for the integration",
                    default='MW14')

# PARSE SAVE FILE OPTIONS
parser.add_argument("-s", "--savedir", type=str, default="results_default",
                    help="directory where the results will be saved")

parser.add_argument("--obs", action="count",
                    help=("save orbit's observables ra[deg] dec[deg], "
                          "d[kpc], pmra[mas/yr], pmde[mas/yr], rv[km/s]"))

parser.add_argument("--xyz", action="count",
                    help="save orbit's x, y, z, r positions in kpc")

parser.add_argument("--rphiz", action="count",
                    help="save orbit's r [kpc], phi [degrees] and z [kpc]")

parser.add_argument("--lb", help="save orbits l, b positions in degrees", 
                    action="count")

parser.add_argument("--pmlb", help="save orbits pml, pmb positions in mas/yr", 
                    action="count")

parser.add_argument("--vxvyvz", action="count",
                    help="save orbit's vx, vy, vz velocities in km/s")

parser.add_argument("--vrvphivz", action="count",
                    help="save orbit's vr, vphi, vz velocities in km/s")

parser.add_argument("--UVW", action="count",
                    help="save orbit's U, V, W velocities in km/s")

parser.add_argument("-a", "--all", action="count",
                    help="compute all available info")

parser.add_argument("-n", "--notsave", action="count",
                    help="will not save any file")

parser.add_argument("--plt", action="count",
                    help="plot available info")

# PARSE ARGS
args = parser.parse_args()

verbose = args.verbose

################################################################################

# Prepare inputs and outputs ###################################################

if verbose: print "Preparing input observables."

# Observables
ra, dec, d = args.ra, args.dec, args.d
pmra, pmde, rv = args.pmra, args.pmde, args.rv

vxvv = [ra, dec, d, pmra, pmde, rv]

if verbose: print "Preparing integration input parameters."

# Integrator params
ro, vo = args.ro, args.vo
N = args.Nsteps
tf = args.time
ts = np.linspace(0, tf, N)*u.Gyr

if args.solar in ['schoenrich']:
    solarmotion = args.solar
else:
    raise ValueError("Unsupported Solar Motion")

if args.pot in ['MW14', 'MWPotential2014']:
    pot = MW14
else:
    raise ValueError("Unsupported Potential")

if verbose: print "Preparing output data frame."

o_data = collections.OrderedDict()
o_data['ts'] = ts

if args.all:
    cols = ['ra', 'dec', 'd', 'pmra', 'pmde', 'rv', 
            'x', 'y', 'z', 'r', 'phi',
            'vx', 'vy', 'vz', 'vr', 'vphi',
            'U','V','W',
            'l', 'b', 'pml', 'pmb']
    
    for key in cols: o_data[key] = np.nan

else:
    
    if args.obs:
        cols = ['ra', 'dec', 'd', 'pmra', 'pmde', 'rv']
        for key in cols: o_data[key] = np.nan
    
    if args.xyz:
        cols = ['x', 'y', 'z']
        for key in cols: o_data[key] = np.nan
    
    if args.rphiz:
        cols = ['r', 'phi', 'z']
        for key in cols: o_data[key] = np.nan
    
    if args.vxvyvz: 
        cols = ['vx', 'vy', 'vz']
        for key in cols: o_data[key] = np.nan
    
    if args.vrvphivz:
        cols = ['vr', 'vphi', 'vz']
        for key in cols: o_data[key] = np.nan
    
    if args.vrvphivz:
        cols = ['U', 'V', 'W']
        for key in cols: o_data[key] = np.nan
    
    if args.lb: 
        cols = ['l', 'b']
        for key in cols: o_data[key] = np.nan
    
    if args.pmlb: 
        cols = ['pml', 'pmb']
        for key in cols: o_data[key] = np.nan

################################################################################

# Integrate orbit ##############################################################

if verbose: print "Initializing orbit"
o = Orbit(vxvv=vxvv, radec = True, ro = ro, vo = vo, solarmotion = solarmotion)

if verbose: print "Integrating orbit"
o.integrate(ts, pot, method = 'leapfrog')

if verbose: print "Saving results to data frame."
if 'ra' in o_data.keys(): o_data['ra'] = o.ra(ts)
if 'dec' in o_data.keys(): o_data['dec'] = o.dec(ts)
if 'd' in o_data.keys(): o_data['d'] = o.dist(ts)
if 'pmra' in o_data.keys(): o_data['pmra'] = o.pmra(ts)
if 'pmde' in o_data.keys(): o_data['pmde'] = o.pmdec(ts)
if 'rv' in o_data.keys(): o_data['rv'] = o.vlos(ts)

if 'x' in o_data.keys(): o_data['x'] = o.x(ts)
if 'y' in o_data.keys(): o_data['y'] = o.y(ts)
if 'z' in o_data.keys(): o_data['z'] = o.z(ts)
if 'r' in o_data.keys(): o_data['r'] = o.R(ts)
if 'phi' in o_data.keys(): o_data['phi'] = o.phi(ts)

if 'vx' in o_data.keys(): o_data['vx'] = o.vx(ts)
if 'vy' in o_data.keys(): o_data['vy'] = o.vy(ts)
if 'vz' in o_data.keys(): o_data['vz'] = o.vz(ts)
if 'vr' in o_data.keys(): o_data['vr'] = o.vR(ts)
if 'vphi' in o_data.keys(): o_data['vphi'] = o.vT(ts)

if 'U' in o_data.keys(): o_data['U'] = o.U(ts)
if 'V' in o_data.keys(): o_data['V'] = o.V(ts)
if 'W' in o_data.keys(): o_data['W'] = o.W(ts)

if 'l' in o_data.keys(): o_data['l'] = o.ll(ts)
if 'b' in o_data.keys(): o_data['b'] = o.bb(ts)
if 'pml' in o_data.keys(): o_data['pml'] = o.pmll(ts)
if 'pmb' in o_data.keys(): o_data['pmb'] = o.pmbb(ts)

################################################################################

# SAVE DATA ####################################################################

savefiles = not args.notsave

# Prepare savefile
savedir = args.savedir
ensure_dir(savedir)

dframe = pd.DataFrame(o_data)

if savefiles:
    # Save input parameters
    if verbose: "Saving input parameters to file '{0}/o_inputs.tsv'.".format(savedir)
    with open(savedir+"/o_inputs.tsv", "w") as f:
        f.write("Orbit calculated for input parameters:\n")
        f.write("RA   = {: 10.5f} degrees\n".format(ra))
        f.write("DEC  = {: 10.5f} degrees\n".format(dec))
        f.write("DIST = {: 10.5f} kpc\n".format(d))
        f.write("PMRA = {: 10.5f} mas/yr\n".format(pmra))
        f.write("PMDE = {: 10.5f} mas/yr\n".format(pmde))
        f.write("RV   = {: 10.5f} km/s\n\n".format(rv))
        f.write("Orbit input parameters:\n")
        f.write("Ro = {: 10.5f} kpc\n".format(ro))
        f.write("Vo = {: 10.5f} km/s\n".format(vo))
        f.write("Integration time: {: 10.5f} Gyr\n".format(tf))
        f.write("Steps           : {0}\n".format(N))
        f.write("Solar Motion    : {0}\n".format(solarmotion))
        f.write("Potential       : {0}\n\n".format(pot))
        f.write("Calculated observables at t = 0\n")
        f.write("U = {: 10.5f} km/s\n".format(o.U(0)[0]))
        f.write("V = {: 10.5f} km/s\n".format(o.V(0)[0]))
        f.write("W = {: 10.5f} km/s\n".format(o.W(0)[0]))
        f.write("l = {: 10.5f} degrees\n".format(o.ll(0)[0]))
        f.write("b = {: 10.5f} degrees\n".format(o.bb(0)[0]))
        f.write("pml = {: 10.5f} mas/yr\n".format(o.pmll(0)[0]))
        f.write("pmb = {: 10.5f} mas/yr".format(o.pmbb(0)[0]))
    
    # Save integration resuls
    if verbose: "Saving integration time steps results to file '{0}/o_results.tsv'.".format(savedir)
    with open(savedir+"/o_results.tsv", "w") as f:
        # Save header
        f.write("#       ts")
        for key in dframe.keys()[1:]:
            f.write("{:>11.10s}".format(key))
        f.write('\n')
        
        np.savetxt(f, np.array(dframe), delimiter = " ", fmt = "% 10.5f")
    
    # Save orbital parameters
    if verbose: "Saving orbital parameters to file '{0}/o_params.tsv'.".format(savedir)
    with open(savedir+"/o_params.tsv", 'w') as f:
        f.write("# Orbit integrated using galpy\n")
        f.write("{:>10.10s} ".format('e'))
        f.write("{:>10.10s} ".format('zmax'))
        f.write("{:>10.10s} ".format('rap'))
        f.write("{:>10.10s}\n".format('rperi'))
        f.write("{:10.6f} ".format(o.e()))
        f.write("{:10.6f} ".format(o.zmax()))
        f.write("{:10.6f} ".format(o.rap()))
        f.write("{:10.6f}".format(o.rperi()))

################################################################################

# PLOT RESULTS #################################################################

if bool_plot is not False:
    if args.plt: bool_plot = True

from galpy_plot import save_plot, plot_xy, plot_xz, plot_yz, plot_xyz, plot_lb

if bool_plot:
    savedir = args.savedir
    ensure_dir(savedir+"/plots")
    
    if args.all or args.xyz:
        if verbose: "Plotting xyz orbits"
        
        plot_xy(o_data, 'x', 'y')
        save_plot(savedir+'/plots/xy.png')
        
        plot_xz(o_data, 'x', 'z')
        save_plot(savedir+'/plots/xz.png')
        
        plot_xz(o_data, 'x', 'z', lim = False)
        save_plot(savedir+'/plots/xz_2.png')
        
        plot_yz(o_data, 'y', 'z')
        save_plot(savedir+"/plots/yz.png")
        
        plot_yz(o_data, 'y', 'z', lim = False)
        save_plot(savedir+"/plots/yz_2.png")
        
        plot_xyz(o_data, lim = True)
        save_plot(savedir+"/plots/xyz.png")
        
        plot_xyz(o_data, lim = False)
        save_plot(savedir+"/plots/xyz_2.png")
        
        plot_xyz(o_data, lim = 2)
        save_plot(savedir+"/plots/xyz_3.png")
    
    if args.all or args.lb:
        plot_lb(o_data['l'], o_data['b'])
        save_plot(savedir+"/plots/lb.png")
        
        if tf >= 0.1:
            plot_lb(o_data['l'][ts < 0.1*u.Gyr], o_data['b'][ts < 0.1*u.Gyr])
            save_plot(savedir+"/plots/lb_2.png")
        
        elif tf <= -0.1:
            plot_lb(o_data['l'][ts > -0.1*u.Gyr], o_data['b'][ts > -0.1*u.Gyr])
            save_plot(savedir+"/plots/lb_2.png")

################################################################################

# PRINT RESULTS ################################################################

if verbose:
    print "\nOrbit calculated for input parameters:"
    print "RA   = {: 10.5f} degrees".format(ra)
    print "DEC  = {: 10.5f} degrees".format(dec)
    print "DIST = {: 10.5f} kpc".format(d)
    print "PMRA = {: 10.5f} mas/yr".format(pmra)
    print "PMDE = {: 10.5f} mas/yr".format(pmde)
    print "RV   = {: 10.5f} km/s".format(rv)
    print ""
    print "Orbit input parameters:"
    print "Ro               = {: 10.5f} kpc".format(ro)
    print "Vo               = {: 10.5f} km/s".format(vo)
    print "Integration time = {: 10.5f} Gyr".format(tf)
    print "Steps            = {: n}".format(N)
    print "Solar Motion     : {0}".format(solarmotion)
    print "Potential        : {0}".format(pot)

if verbose > 1:
    print "\nCalculated observables at t = 0"
    print "U = {: 10.5f} km/s".format(o.U(0)[0])
    print "V = {: 10.5f} km/s".format(o.V(0)[0])
    print "W = {: 10.5f} km/s".format(o.W(0)[0])
    print "l = {: 10.5f} degrees".format(o.ll(0)[0])
    print "b = {: 10.5f} degrees".format(o.bb(0)[0])
    print "pml = {: 10.5f} mas/yr".format(o.pmll(0)[0])
    print "pmb = {: 10.5f} mas/yr".format(o.pmbb(0)[0])
    
print "\nDerived orbital parameters:"
print "eccentricity = {: 10.5f}".format(o.e())
print "zmax         = {: 10.5f} kpc".format(o.zmax())
print "rap          = {: 10.5f} kpc".format(o.rap())
print "rperi        = {: 10.5f} kpc".format(o.rperi())

################################################################################

