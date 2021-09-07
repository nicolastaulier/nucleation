import math
import numpy
import scipy.constants

###############################################################################
def P_out(P,R,g):
    r""" It derives the pressure out of a droplet from the pressure inside a droplet
    
    .. math::
        P_{out} = P_{in} - \frac{2\gamma}{R}
        
    :param P: Pressure :math:`P_{in}` inside the droplet, in Pa.
    :type command: float

    :param R: Pressure :math:`R` of the droplet, in m.
    :type command: float

    :param g: Interfacial tension :math:`\gamma` of the interface between the inside and outside space of the droplet, in N/m.
    :type command: float
        
    :return: :math:`P_{out}` outside the droplet, in Pa.
    :rtype: float        
                
    """
    value = P - 2*g/R
    return value
###############################################################################

###############################################################################
def P_in(P,R,g):
    r""" It derives the pressure inside droplet from the pressure outside a droplet

    .. math::
        P_{in} = P_{out} +  \frac{2\gamma}{R}

    :param P: Pressure :math:`P_{out}` outside the droplet, in Pa.
    :type command: float

    :param R: Pressure :math:`R` of the droplet, in m.
    :type command: float

    :param g: Interfacial tension :math:`\gamma` of the interface between the inside and outside space of the droplet, in N/m.
    :type command: float
        
    :return: :math:`P_{in}` inside the droplet, in Pa.
    :rtype: float        
                
    """
    value = P + 2*g/R
    return value
###############################################################################



