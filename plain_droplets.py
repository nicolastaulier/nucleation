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
def plot_func():
    """ It plots a surface plot of the solution of the equation :math:`W^{*}-W_{0.5}` as a function of two
    interfacial tensions
    """
    R =  numpy.array([0.1,2.5,5,20,30])*1e-6
    P = numpy.array([2.49,2.39,2.35,2.08,1.34])*1e6
    n = numpy.array([2.4e10,1.5e6,1.9e5,2.9e3,8.8e2]) 
    #R = 0.1e-6
    #P = 2.49e6
    #n = 1.4e10
    
    T = 293.15 # kelvin
    freq = 1.1e6 # Hz
    PFH = PFC.PFC('PFH')
    
    concave =   heterogeneous.concaveSurface(PFH,T,freq)
    marker_size = 12
    markerwidth = 2
    width = 1
    cap = 3
    
    fig = matplotlib.pyplot.figure(num=0,figsize=(7,5),facecolor='w',edgecolor='w')
    ax = fig.add_subplot(projection='3d')
    
    g_gl = numpy.arange(0.001,0.05,0.0001,dtype=float)
    g_lw = numpy.arange(0.002,0.05,0.0001,dtype=float)
    g_gw = 0.010
    
    Z = numpy.zeros((len(g_gl),len(g_lw)),dtype=float)
    for l in range(0,len(P)):
        i = 0
        for gl in g_gl:
            j = 0
            for diff in g_lw:
                g = (diff,g_gw,gl)
                Z[i,j] =  concave.func_3g_b(g,P[l],R[l],n[l])
                j +=1
            i+=1
        X,Y = numpy.meshgrid(g_lw,g_gl)    
        ax.plot_surface(X,Y,Z,rstride=10)
    
    ax.set_xlabel(r'$\gamma_{lw}$')
    ax.set_xlim(0,0.05)
    ax.set_ylabel(r'$\gamma_{gl}$')
    ax.set_ylim(0,0.05)
    ax.set_zlabel(r'$W^{*} - W_{0.5}$')
    ax.set_zlim(-0.25e-17,0.2e-17)
    ax.set_title(r'Plain PFH droplets')    
    matplotlib.pyplot.savefig('../Surf.png',dpi=600)
    return 1
###############################################################################

###############################################################################
def example1():
    """ It plot the behavior of  the acoustic pressure of vaporization as a function of droplet radius for one droplet.
    """
    # Experimental data of droplet radius and vaporization pressure at p = 1/2 
    R =  numpy.arange(0.05e-6,30e-6,0.1e-6)
    N = len(R)
    #print(R)
    n = 1    
    T = 293.15 # kelvin
    f = 1.1e6 # Hz
    PFH = PFC.PFC('PFH')

    # interfacial tension around the nuclei
    g_lw = 0.025
    g_gw =  0.012
    g_gl = 0.011
    surface =  heterogeneous.concaveSurface(PFH,T,f)
 
    P_r = numpy.zeros(N)
    P_l = numpy.zeros(N)
    for i in range(0,N):
        args = (R[i],n,g_lw,g_gw,g_gl)
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
    matplotlib.pyplot.subplots_adjust(left=0.150,right=0.97,bottom=0.170,top=0.970,hspace = 0.12,wspace = 0)
    matplotlib.pyplot.tick_params(axis='both',which='both',direction ='in',bottom=True,top=True,left=True,right=True)
    
    matplotlib.pyplot.scatter(R*1e6, P_r*1e-6)
    #matplotlib.pyplot.scatter(R*1e6, P_l*1e-6)
        
    matplotlib.pyplot.xlabel(r'$R$ ($\mu$m)',fontsize = 28)
    matplotlib.pyplot.ylabel(r'$P_{0.5}$ (MPa)',fontsize = 28)
    matplotlib.pyplot.xticks(fontsize = 24)
    matplotlib.pyplot.yticks(fontsize = 24)
    matplotlib.pyplot.xlim(0,30)
    #matplotlib.pyplot.ylim(0,4)
    #ax.set_yticks([0,1,2,3,4])
    #ax.xaxis.set_minor_locator(MultipleLocator(1))
    #ax.yaxis.set_minor_locator(MultipleLocator(0.2))
    matplotlib.pyplot.savefig('../fig_PvsR_1droplet.png',dpi=600)
    return 1
###############################################################################

###############################################################################
def example2():
    """ It plot the behavior of  the acoustic pressure of vaporization as a function of number of droplets 
    for a constant radius
    """
    # Experimental data of droplet radius and vaporization pressure at p = 1/2
    n =  numpy.arange(1e2,3e5,1e2)
    N = len(n)
    R = 40e-6
    T = 293.15 # kelvin
    f = 1.1e6 # Hz
    PFH = PFC.PFC('PFH')

    # interfacial tension around the nuclei
    g_lw = 0.025
    g_gw =  0.012
    g_gl = 0.011
    surface =   heterogeneous.concaveSurface(PFH,T,f)
 
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
    """ It calculates the surface tensions of each of our experimental acoustic pressure,
    then calculate the acoustic pressure from mean surface tension and plot the resulting pressures.
    """
    # Experimental data of droplet radius and vaporization pressure at p = 1/2 
    R =  numpy.array([0.1,2.5,5,20,30])*1e-6
    P = numpy.array([2.49,2.39,2.35,2.08,1.34])*1e6
    err_P = numpy.array([0.1,0.05,0.1,0.1,0.1])*1e6
    n = numpy.array([2.4e11,1.5e6,1.9e5,2.9e3,8.8e2])    
    T = 293.15 # kelvin
    f = 1.1e6 # Hz
    PFH = PFC.PFC('PFH')
    N = len(P)

    # interfacial tension around the nuclei
    g_lw = 0.040
    g_gw =  0.010
    g_gl = 0.020
    surface =   heterogeneous.concaveSurface(PFH,T,f)
 
    # The equation was solved for all pressure to find common values of g_lw, g_gw and g_gl
    g = numpy.array([0.025,0.020,0.015])
    bounds = ([0.005,0.005,0.005],[0.05,0.05,0.050])
    sol = scipy.optimize.least_squares(surface.func3, g, method='trf', bounds=bounds)
    print(sol.x)
    # Using these values, the pressures are calculated
    solP = numpy.zeros(N)
    for i in range(0,N):
        args = (R[i],n[i],sol.x[0],sol.x[1],sol.x[2])
        P_root = scipy.optimize.root(surface.func_P,1e6,args=args,method='hybr')
        P_ls = scipy.optimize.least_squares(surface.func_P,1e6,args=args,method='lm')
        print(round(P_root.x[0]*1e-6,2),round(P_ls.x[0]*1e-6,2))
        solP[i] = P_ls.x[0]

    # The equations were solved separately for each pressure, 
    # then a derived values for g_diff and g_gl at each pressure
    P_r = numpy.zeros(N)
    for i in range(0,N):
        args = (P[i],R[i],n[i])
        g = numpy.array([g_lw,g_gw,g_gl])
        sol2 = scipy.optimize.root(surface.func_3g,g,args=args,method='hybr')
        sol3 = scipy.optimize.least_squares(surface.func_3g,g,args=args,method='lm')
        print('gamma ',sol2.x, sol3.x)
        args = (R[i],n[i],sol2.x[0],sol2.x[1],sol2.x[2])  
        P_root = scipy.optimize.root(surface.func_P,1e6,args=args,method='hybr')
        print(P_root.x[0])
        P_r[i] =P_root.x[0]
        args = (R[i],n[i],sol3.x[0],sol3.x[1],sol3.x[2])
        P_ls = scipy.optimize.least_squares(surface.func_P,1e6,args=args,method='lm')
        #print(P_ls)
        print('P',round(P_root.x[0]*1e-6,1),round(P_ls.x[0]*1e-6,1))
       
    
    marker_size = 12
    markerwidth = 2
    width = 1
    cap = 3
    
    fig = matplotlib.pyplot.figure(num=0,figsize=(7,5),facecolor='w',edgecolor='w')
    ax = fig.add_subplot(111)
    matplotlib.pyplot.subplots_adjust(left=0.150,right=0.97,bottom=0.170,top=0.970,hspace = 0.12,wspace = 0)
    matplotlib.pyplot.tick_params(axis='both',which='both',direction ='in',bottom=True,top=True,left=True,right=True)
    
    matplotlib.pyplot.errorbar(R*1e6, P*1e-6, yerr = err_P*1e-6, marker='o', markeredgecolor='k', markerfacecolor='none', markersize=marker_size ,linestyle='', fillstyle='full', markeredgewidth=markerwidth, clip_on=False, ecolor='k', capsize=cap, elinewidth=width)
    
    
    matplotlib.pyplot.scatter(R*1e6, solP*1e-6)
    matplotlib.pyplot.scatter(R*1e6, P_r*1e-6)
    
    #matplotlib.pyplot.annotate(text = r'$2.41 - 6.86 \times 10^{-5} \times R^{2.84}$', xy = (10,3.2), fontsize = 28)
    #(chisq,p) = scipy.stats.chisquare(P,fit(R,*popt))
    #print('chisquare = ' + str(chisq) + ' p-value = ' + str(p))
    #print('Pearson',scipy.stats.pearsonr(R,P))
    
    matplotlib.pyplot.xlabel(r'$R$ ($\mu$m)',fontsize = 28)
    matplotlib.pyplot.ylabel(r'$P_{0.5}$ (MPa)',fontsize = 28)
    matplotlib.pyplot.xticks(fontsize = 24)
    matplotlib.pyplot.yticks(fontsize = 24)
    #matplotlib.pyplot.xlim(0,35)
    #matplotlib.pyplot.ylim(0,4)
    #ax.set_yticks([0,1,2,3,4])
    #ax.xaxis.set_minor_locator(MultipleLocator(1))
    #ax.yaxis.set_minor_locator(MultipleLocator(0.2))
    matplotlib.pyplot.savefig('../fig_PvsR.png',dpi=600)

    return 1
###############################################################################

#plot_func()
#example1()
#example2()
#example3()
#matplotlib.pyplot.show()

