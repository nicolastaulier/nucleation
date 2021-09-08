import math
import numpy
import scipy.constants


###############################################################################
def surfaceTension(R,g_0,delta):
    r""" It provides the effective surface tension :math:`\gamma`  according to Tolman.
        
     .. math::
        \gamma = \frac{\gamma_0}{1+\frac{2\delta}{R}}
        
        
    See  The effect of droplet size on surface tension.  J. Chem. Physics 17 (1949) 333-337.
    doi:10.1063/1.1747247
    
    :param R: Radius of a droplet :math:`R`, in m.
    :type command: float
        
    :param g_0: Surface tension :math:`\gamma_{0}` of a flat interface, in N/m.
    :type command: float

    :param delta: Tolman length :math:`\delta`, in m.
    :type command: float
        
    :return: surface tension :math:`\gamma` of a curve surface of radius R.
    :rtype: float
    
    """
    value = g_0/(1+2*delta/R)
    return value
###############################################################################

###############################################################################
def surfaceTensionRc(P,g_0,delta):
    r""" It provides the effective surface tension :math:`\gamma`  for :math:`r = r^{*}`.
        
    .. math::
        \gamma = \gamma_0 - \delta P
        
    
    :param P: Pressure :math:`P`, in Pa. 
    :type command: float
    
    :param g_0: Surface tension :math:`\gamma_{0}` of a flat interface, in N/m.
    :type command: float

    :param delta: Tolman length :math:`\delta`, in m.
    :type command: float
        
    :return: Surface tension :math:`\gamma` of a curved surface of radius :math:`r = r^{*}`, in N/m.
    :rtype: float
    """
    value = g_0 - delta*P
    return value
###############################################################################

###############################################################################
def Gamma_0(g,rho,M):
    r""" It provides the nucleation rate :math:`\Gamma_{0}`
    
    .. math:: 
            \Gamma_{0} = N_A \rho \sqrt{\frac{2\gamma}{\pi M}}


    :param g: Surface tension :math:`\gamma` of a curved interface, in N/m.
    :type command: float

    :param rho: Density :math:`\rho` of the medium, in kg/m.
    :type command: float

    :param M: Mass :math:`M` of one molecule of the medium, in kg.
    :type command: float
        
    :return: :math:`\Gamma_{0}`
    :rtype: float        
        
    """
    value =  scipy.constants.N_A*rho*numpy.sqrt(2*g/(math.pi*M))
    return value
###############################################################################

###############################################################################
def P_05(T,G0,R,f,n,P,g):
    r""" It derive the pressure vaporisation :math:`P_{p}` at probability :math:`p = 0.5`.
    
    .. math::
            P_{0.5} = \sqrt{\frac{16\pi \gamma^{3}}{3k_{B}T}\frac{1}{\ln{\frac{\Gamma_{0}\pi n4R^{3}}{6f\ln{2}}}}}
    
    where :math:`k_{B}` is the Boltzmann constant, :math:`T` is the temperature, :math:`f` is the acoustic frequency,
    :math:`R` is the radius of the droplet and :math:`n` the number of droplets, :math:`\gamma` is the interface tension
    at the interface between the inner and outer space of the droplet and :math:`Gamma_{0}` is the nucleation rate.
    
    
    :param T: Temperature :math:`T`, in K
    :type command: float

    :param G0:  Nucleation rate :math:`\Gamma_{0}`,  
    :type command: float

    :param R:  Radius :math:`R` in meter of a droplet, in m
    :type command: float

    :param f: Acoustic frequency :math:`f`, in Hz
    :type command: float    

    :param n: Number of droplets :math:`n`
    :type command: float    
    
    :param g: Surface tensions :math:`\gamma` of a curved interface, in N/m
    :type command: float
        
    :return: Acoustic pressure :math:`P_{0.5}`, in Pa
    :rtype: float     
       
    """
    value = numpy.sqrt((16*math.pi*numpy.power(g,3))/(3*scipy.constants.k*T*numpy.log((G0*math.pi*n*4*numpy.power(R,3))/(6*f*math.log(2)))))
    return value
###############################################################################

###############################################################################
def r_c(g_0,delta,P):
    r""" It derives the critical radius of a nucleus from the Tolman equation
    
    .. math::     r^{*} = \frac{2 \gamma_{0}}{P^{*}} - 2 \delta
    
    :param g_0: surface tension of a flat interface :math:`\gamma_{0}`, in N/m
    :type command: class

    :param delta: Tolman length :math:`\delta`, in meter
    :type command: float

    :param P: Critical pressure  :math:`P_{*}` (in Pa) at which vaporization occurs
    :type command: float
        
    :return: :math:`r^{*}`
    :rtype: float     
    
    """
    value = 2*gamma_0/P - 2*delta
    return value
###############################################################################

###############################################################################
def delta(g,g_0,P):
    r""" It derives the Tolman length from the Tolman equation
    
    .. math::
        \delta = \frac{\gamma_{0}-\gamma}{P} 
    
    :param g_0: surface tension of a flat interface :math:`\gamma_{0}` in N/m
    :type command: class

    :param g: surface tension of a curved interface in N/m
    :type command: float

    :param P: Critical pressure at which vaporization occurs
    :type command: float
        
    :return: :math:`\delta`
    :rtype: float     
    
    """
    value = (g_0 - g)/P
    return value
###############################################################################

