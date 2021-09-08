# -*- coding: utf-8 -*-
from PyQt5 import QtGui
import numpy
import math
import scipy.optimize

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot
package = {'text.latex.preamble' : [r'\usepackage{textcomp}']}
matplotlib.pyplot.rcParams.update(package)
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{textcomp}"]
matplotlib.pyplot.rc('font',family='serif')
matplotlib.pyplot.rc('text',usetex=True)

from mpl_toolkits import mplot3d

import scipy.optimize,scipy.stats
import heterogeneous
import PFC
import pressure


###############################################################################
def example1():
    """ It creates a plot showing the behavior of :math:`P_{0.5}`, the acoustic pressure that induces vaporization, 
    as a function of the droplet radius for (a number of) one droplet.
    """
    # Experimental data of droplet radius and vaporization pressure at p = 1/2 
    R_w =  numpy.arange(0.05e-6,30e-6,0.01e-6)
    N = len(R_w)
    #print(R)
    n = 1    
    T = 293.15 # kelvin
    f = 1.1e6 # Hz
    PFH = PFC.PFC('PFH')

    # interfacial tension around the nuclei
    g_lw = 0.030
    g_gw = 0.010
    g_gl = 0.022
    surface = heterogeneous.convexSurface(PFH,T,f)
 
    P_r = numpy.zeros(N)
    P_l = numpy.zeros(N)
    for i in range(0,N):
        args = (R_w[i],n,g_lw,g_gw,g_gl)
        P_root = scipy.optimize.root(surface.func_P,1e6,args=args,method='hybr')
        P_ls = scipy.optimize.least_squares(surface.func_P,1e6,args=args,method='lm')
        P_r[i] = P_root.x[0]
        P_l[i] = P_ls.x[0]
    
    marker_size = 12
    markerwidth = 2
    width = 1
    cap = 3
    
    fig = matplotlib.pyplot.figure(num=0,figsize=(7,5),facecolor='w',edgecolor='w')
    ax = fig.add_subplot(111)
    matplotlib.pyplot.subplots_adjust(left=0.175,right=0.97,bottom=0.170,top=0.970,hspace = 0.12,wspace = 0)
    matplotlib.pyplot.tick_params(axis='both',which='both',direction ='in',bottom=True,top=True,left=True,right=True)
    
    matplotlib.pyplot.scatter(R_w*1e6, P_r*1e-6)
    #matplotlib.pyplot.scatter(R_w*1e6, P_l*1e-6)
        
    matplotlib.pyplot.xlabel(r'$R$ ($\mu$m)',fontsize = 28)
    matplotlib.pyplot.ylabel(r'$P_{0.5}$ (MPa)',fontsize = 28)
    matplotlib.pyplot.xticks(fontsize = 24)
    matplotlib.pyplot.yticks(fontsize = 24)
    matplotlib.pyplot.xlim(0,30)
    #matplotlib.pyplot.ylim(2,2.2)
    #ax.set_yticks([0,1,2,3,4])
    #ax.xaxis.set_minor_locator(MultipleLocator(1))
    #ax.yaxis.set_minor_locator(MultipleLocator(0.2))
    matplotlib.pyplot.savefig('../fig_P_vs_R_1droplet.png',dpi=600)
    return 1
###############################################################################

###############################################################################
def example2():
    """  It creates a plot showing the behavior of :math:`P_{0.5}`, the acoustic pressure that induces vaporization, 
    as a function of number of droplets for a constant radius
    """
    # Experimental data of droplet radius and vaporization pressure at p = 1/2 
    n =  numpy.arange(1e2,3e5,1e2)
    N = len(n)
    R = 20e-6
    T = 293.15 # kelvin
    f = 1.1e6 # Hz
    PFH = PFC.PFC('PFH')

    # interfacial tension around the nuclei
    g_lw = 0.030
    g_gw = 0.010
    g_gl = 0.022
    surface =   heterogeneous.convexSurface(PFH,T,f)
 
    P_r = numpy.zeros(N)
    #P_l = numpy.zeros(N)
    for i in range(0,N):
        args = (R,n[i],g_lw,g_gw,g_gl)
        P_root = scipy.optimize.root(surface.func_P,1e6,args=args,method='hybr')
        #P_ls = scipy.optimize.least_squares(surface.func_P,1e6,args=args,method='lm')
        P_r[i] = P_root.x[0]
        #P_l[i] = P_ls.x[0]
    
    marker_size = 12
    markerwidth = 2
    width = 1
    cap = 3
    
    fig = matplotlib.pyplot.figure(num=0,figsize=(7,5),facecolor='w',edgecolor='w')
    ax = fig.add_subplot(111)
    matplotlib.pyplot.subplots_adjust(left=0.175,right=0.97,bottom=0.170,top=0.970,hspace = 0.12,wspace = 0)
    matplotlib.pyplot.tick_params(axis='both',which='both',direction ='in',bottom=True,top=True,left=True,right=True)
    
    matplotlib.pyplot.scatter(n*1e-5, P_r*1e-6)
    #matplotlib.pyplot.scatter(n, P_l*1e-6)
        
    matplotlib.pyplot.xlabel(r'$n$ ($10^{5}$ droplets)',fontsize = 28)
    matplotlib.pyplot.ylabel(r'$P_{0.5}$ (MPa)',fontsize = 28)
    matplotlib.pyplot.xticks(fontsize = 24)
    matplotlib.pyplot.yticks(fontsize = 24)
    matplotlib.pyplot.xlim(0,3)
    #matplotlib.pyplot.ylim(0,4)
    #ax.set_yticks([0,1,2,3,4])
    #ax.xaxis.set_minor_locator(MultipleLocator(1))
    #ax.yaxis.set_minor_locator(MultipleLocator(0.2))
    matplotlib.pyplot.savefig('../fig_PvsR.png',dpi=600)

    return 1
###############################################################################


###############################################################################
def example3():
    """  It derive the values of surface tension from the experimental values of :math:`P_{0.5}`, the acoustic pressure that induces vaporization
    made from a single and multi core droplets.
    Then identical values of surface tensions, the :math:`P_{0.5}` values are derived for the two cases.
    """
    # The external radius of the single and multi-core droplets is
    R = 20e-6
    T = 293.15 # kelvin
    f = 1.1e6 # Hz
    PFH = PFC.PFC('PFH')

    # interfacial tension around the nuclei
    g_lw = 0.025
    g_gw = 0.020
    g_gl = 0.022
    surface =   heterogeneous.convexSurface(PFH,T,f)

    P_ex_single = 1.5e6
    Rw_single = 10.5e-6  
    n_single = 1 

    P_ex_multi = 2.2e6    
    Rw_multi = 0.215e-6 
    n_multi = math.pow(R,3)/math.pow(Rw_multi,3)*0.4
    
    args = (P_ex_single,Rw_single,n_single)
    P_single = scipy.optimize.root(surface.func_3g,(g_lw,g_gw,g_gl),args=args,method='hybr')
    P_ls_single = scipy.optimize.least_squares(surface.func_3g,(g_lw,g_gw,g_gl),args=args,method='lm')
    print(P_single.x,P_ls_single.x)
    
    args = (P_ex_multi,Rw_multi,n_multi)
    P_multi = scipy.optimize.root(surface.func_3g,(g_lw,g_gw,g_gl),args=args,method='hybr')
    P_ls_multi = scipy.optimize.least_squares(surface.func_3g,(g_lw,g_gw,g_gl),args=args,method='lm')
    print(P_multi.x,P_ls_multi.x)
    
    # The equation was solved for all pressure to find common values of g_lw, g_gw and g_gl
    g = numpy.array([0.020,0.015,0.020])
    sol = scipy.optimize.least_squares(surface.func3,g,method='lm')
    print(sol.x)
    
    args = (Rw_single,n_single,P_ls_single.x[0],P_ls_single.x[1],P_ls_single.x[2])
    P_s = scipy.optimize.root(surface.func_P,1e6,args=args,method='lm')
    args = (Rw_multi,n_multi,P_ls_multi.x[0],P_ls_multi.x[1],P_ls_single.x[2])
    P_m = scipy.optimize.root(surface.func_P,1e6,args=args,method='lm')   
    print(P_s.x[0]*1e-6,P_m.x[0]*1e-6)    
    return 1
###############################################################################


###############################################################################
def example4():
    """  It creates a plot showing the behavior of :math:`P_{0.5}`, the acoustic pressure that induces vaporization,
    as a function of volume fraction of water inside a droplet.
    """
    R = 0.215e-6
    T = 293.15 # kelvin
    f = 1.1e6 # Hz
    PFH = PFC.PFC('PFH')
    
    # interfacial tension around the nuclei
    g_lw = 0.025
    g_gw = 0.020
    g_gl = 0.022
    surface =   heterogeneous.convexSurface(PFH,T,f)


    phi = numpy.array([0.05,0.4,0.6,0.8])
    N = len(phi)
    Rw = 0.215e-6
    P_r = numpy.zeros(N)

    for i in range(0,N):
        n = math.pow(R,3)/math.pow(Rw,3)*phi[i]
        print(n)
        args = (Rw,n,g_lw,g_gw,g_gl)
        P_multi = scipy.optimize.root(surface.func_P,1e6,args=args,method='hybr')
        P_r[i] = P_multi.x[0]
    
    marker_size = 12
    markerwidth = 2
    width = 1
    cap = 3
    
    fig = matplotlib.pyplot.figure(num=0,figsize=(7,5),facecolor='w',edgecolor='w')
    ax = fig.add_subplot(111)
    matplotlib.pyplot.subplots_adjust(left=0.175,right=0.97,bottom=0.170,top=0.970,hspace = 0.12,wspace = 0)
    matplotlib.pyplot.tick_params(axis='both',which='both',direction ='in',bottom=True,top=True,left=True,right=True)
    
    matplotlib.pyplot.scatter(phi, P_r*1e-6)
        
    matplotlib.pyplot.xlabel(r'$\phi_{w}$',fontsize = 28)
    matplotlib.pyplot.ylabel(r'$P_{0.5}$ (MPa)',fontsize = 28)
    matplotlib.pyplot.xticks(fontsize = 24)
    matplotlib.pyplot.yticks(fontsize = 24)
    matplotlib.pyplot.xlim(0,1)
    #matplotlib.pyplot.ylim(0,4)
    #ax.set_yticks([0,1,2,3,4])
    #ax.xaxis.set_minor_locator(MultipleLocator(1))
    #ax.yaxis.set_minor_locator(MultipleLocator(0.2))
    matplotlib.pyplot.savefig('../fig_multi_phiw.png',dpi=600)

        
    return 1
###############################################################################

#example1()
#example2()
#example3()
#example4()
#matplotlib.pyplot.show()

