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

import sys, warnings, argparse

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
    raise ImportError(import_error_msg.format("astropy"))

# matplotlib is used to produce the plots
bool_plot = None
try:
    from matplotlib import pyplot as plt
except ImportError:
    warnings.warn(import_error_msg.format("matplotlib")+
                  "\nThe code will continue but no plots will be made")
    bool_plot = False

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
parser.add_argument("-s", "--savefile", type=str, default="orbit_result.tsv",
                    help="file where the values will be stored")

parser.add_argument("--obs", action="count",
                    help=("save orbit's observables ra[deg] dec[deg], "
                          "d[kpc], pmra[mas/yr], pmde[mas/yr], rv[km/s]"))

parser.add_argument("--xyz", action="count",
                    help="save orbit's x, y, z, r positions in kpc")

parser.add_argument("--rphi", action="count",
                    help="save orbit's r [kpc] and phi [degrees]")

parser.add_argument("--lb", help="save orbits l, b positions in degrees", 
                    action="count")

parser.add_argument("--vrvphi", action="count",
                    help="save orbit's vr, vphi, vz velocities in km/s")

# PARSE PLOT OPTIONS
parser.add_argument("-plt_xyz", help="if true, plot x,y and x,z results", 
                    action="store_true")

parser.add_argument("-plt_lb", help="if true, plot l,b results", 
                    action="store_true")

# PARSE ARGS
args = parser.parse_args()

# Observables
ra, dec, d = args.ra, args.dec, args.d
pmra, pmde, rv = args.pmra, args.pmde, args.rv

vxvv = [ra, dec, d, pmra, pmde, rv]

# Integrator params
ro, vo = args.ro, args.vo
N = args.Nsteps
ts = np.linspace(0, args.time, N)*u.Gyr

if args.solar in ['schoenrich']:
    solar = args.solar
else:
    raise ValueError("Unsupported Solar Motion")

if args.pot in ['MW14', 'MWPotential2014']:
    pot = MW14
else:
    raise ValueError("Unsupported Potential")

# Preparing save options
Ncols = 1 # Always include a column for the time
save_header = "t         " + save_header

if args.obs:
    Ncols+=6
    save_header += "ra        dec       d         pmra      pmde      rv        "

if args.xyz:
    Ncols+=3
    save_header += "x         y         z         "

if args.rphi:
    Ncols+=2
    save_header += "r         phi       "

if args.vrvphivz: 
    Ncols+=3
    save_header += "vr        vphi      vz        "

if args.lb: 
    Ncols+=2
    save_header += "l         b         "

# Prepare savefile
savefile = args.savefile

# Prepare array to store the calculated data
save_array = np.empty((N, Ncols))
save_array[:] = np.nan

################################################################################

# Integrate orbit ##############################################################

o = Orbit(vxvv=vxvv, radec = True, ro = ro, vo = vo, solarmotion = solarmotion)
o.integrate(ts, pot, method = 'leapfrog')

save_array[:,0] = ts
i = 1
if args.obs:
    save_array[:,i] = o.ra(ts); i+=1
    save_array[:,i] = o.dec(ts); i+=1
    save_array[:,i] = o.dist(ts); i+=1
    save_array[:,i] = o.pmra(ts); i+=1
    save_array[:,i] = o.pmdec(ts); i+=1
    save_array[:,i] = o.vlos(ts); i+=1

if args.xyz:
    save_array[:,i] = o.x(ts); i+=1
    save_array[:,i] = o.y(ts); i+=1
    save_array[:,i] = o.z(ts); i+=1

if args.rphi:
    save_array[:,i] = o.R(ts); i+=1
    save_array[:,i] = o.phi(ts); i+=1

if args.vrvphivz:
    save_array[:,i] = o.R(ts); i+=1
    save_array[:,i] = o.phi(ts); i+=1


################################################################################


    
print args
print vxvv
print ro, vo, solar, pot, ts
print Ncols
print save_header

