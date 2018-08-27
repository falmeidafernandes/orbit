#!/ust/bin/env

################################################################################
__author__  = "Felipe Almeida-Fernandes"
__emais__   = ["felipe11@astro.ufrj.br", "felipefer42@gmail.com"]
__version__ = 1
__date__    = "2018-08-26"

# Import the necessary modules #################################################

import_error_msg = ("could not import library {0}. You may try to install it " 
                    "using:\n$ pip install {0}")

try:
    import numpy as np
except ImportError:
    raise ImportError(import_error_msg.format("numpy"))

try:
    from matplotlib import pyplot as plt
except ImportError:
    raise ImportError(import_error_msg.format("matplotlib"))
else:
    from mpl_toolkits.mplot3d import Axes3D

################################################################################

def save_plot(savefile):
    plt.savefig(savefile)
    plt.close()
    plt.clf()

# PLOT xyz #####################################################################

def plot_x_y_z(o_data, x1, x2, lim = True):
        plt.plot(o_data[x1], o_data[x2])
        
        if lim:
            ax_min = np.min((o_data[x1].min(), o_data[x2].min()))
            ax_max = np.max((o_data[x1].max(), o_data[x2].max()))
            
            DX = ax_max-ax_min
            ax_min = ax_min - 0.1*DX
            ax_max = ax_max + 0.1*DX
            
            plt.gca().set_xlim([ax_min, ax_max])
            plt.gca().set_ylim([ax_min, ax_max])
        
        plt.gca().set_xlabel(x1+' (kpc)')
        plt.gca().set_ylabel(x2+' (kpc)')

plot_xy = plot_x_y_z
plot_xz = plot_x_y_z
plot_yz = plot_x_y_z

def plot_xyz(o_data, lim = 2):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    ax.plot(o_data['x'], o_data['y'], o_data['z'])
    
    if lim == 1:
        ax_min = np.min((o_data['x'].min(), o_data['y'].min()))
        ax_max = np.max((o_data['x'].max(), o_data['y'].max()))
    elif lim == 2:
        ax_min = np.min((o_data['x'].min(), o_data['y'].min(), o_data['z'].min()))
        ax_max = np.max((o_data['x'].max(), o_data['y'].max(), o_data['z'].max()))
    
    if lim:
        DX = ax_max - ax_min
        ax_min = ax_min - 0.1*DX
        ax_max = ax_max + 0.1*DX
        
        ax.set_xlim([ax_min, ax_max])
        ax.set_ylim([ax_min, ax_max])
    
    if lim == 2:
        ax.set_zlim([ax_min, ax_max])
    
    ax.set_xlabel('x (kpc)')
    ax.set_ylabel('y (kpc)')
    ax.set_zlabel('z (kpc)')
    
################################################################################

# PLOT COORDS ##################################################################

radconv = np.pi/180.


def plot_grid(alpha1 = 0.2, alpha2 = 0.1, N = 1000, **kargs):
    coord1_grid = np.linspace(-180, 180, N)*radconv
    coord2_grid = np.linspace(-90, 90, N)*radconv
    
    for coord1, coord2 in zip([-120, -60, 0, 60, 120], [-60, -30, 0, 30, 60]):
        plt.plot([coord1*radconv]*N, coord2_grid, alpha = alpha1, **kargs)
        plt.plot(coord1_grid, [coord2*radconv]*N, alpha = alpha1, **kargs)
    
    for coord1, coord2 in zip([-150, -90, -30, 30,  90, 150], [-75, -45, -15, 15, 45, 75]):
        plt.plot([coord1*radconv]*N, coord2_grid, alpha = alpha2, **kargs)
        plt.plot(coord1_grid, [coord2*radconv]*N, alpha = alpha2, **kargs)



def plot_coord(coord1, coord2, zorder = 1, **kargs):
    Npoints = len(coord1)
    
    # Find boundary breaks at -180 and 180
    boundary_breaks = []
    
    # crossings for 180
    for j in range(1, Npoints):
        delta_max = 10 # only consider as break with delta < delta_max
        # crossings for 180+
        if (coord1[j-1] <= 180) and (coord1[j] > 180):
            delta = abs(coord1[j-1] - coord1[j])
            if delta < delta_max:
                boundary_breaks.append(j)
        # crossings for 180-
        elif (coord1[j-1] >= 180) and (coord1[j] < 180):
            delta = abs(coord1[j-1] - coord1[j])
            if delta < delta_max:
                boundary_breaks.append(j)
    
    boundary_breaks.append(Npoints)
    
    # fix coord1 > 180
    coord1[coord1 >= 180] = coord1[coord1 >= 180] - 360
    
    # plots
    j0 = 0
    for k in range(len(boundary_breaks)):
        jf = boundary_breaks[k]
        plt.plot(coord1[j0:jf]*radconv, coord2[j0:jf]*radconv, zorder = zorder, **kargs)
        j0 = jf   

################################################################################

# PLOT lb ######################################################################

def plot_lb(l, b):
    plt.figure(figsize = (12,6.3))
    plt.subplot(111, projection="aitoff")
    
    plot_grid(alpha1 = 0.2, alpha2 = 0.2, color = '#666666', lw = 1, zorder = 0)
    plot_coord(l, b, alpha = 0.5, ls = '-', lw = 2, color = "#1f77b4")
    
    plt.tick_params(labelbottom=False, labelleft=False)
    plt.subplots_adjust(bottom = 0.01, top = 0.99, right = 0.99, left = 0.01)

################################################################################
