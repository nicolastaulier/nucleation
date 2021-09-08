import math
import numpy
import scipy.constants


import pressure

###############################################################################
class curvedSurface():

    def __init__(self,PFC,T,freq):
        r""" It set the values of the interface tension between PFC liquid/water (:math:`\gamma_{lw}`),
        PFC gas/water (:math:`\gamma_{gw}`) and PFC gas/ PFC liquid (:math:`\gamma_{gl}`).
        The class PFC is asked.

        :param T: The temperature :math:`T`, in K.
        :type command: float

        :param freq: The acoustic frequency :math:`f` of the acoustic wave, in Hz.
        :type command: float

        :param PFC: Type of perfluorocarbone.
        :type command: string
        
        """
        self.PFC = PFC
        self.T = T
        self.freq = freq
    
    def Pi_0(self,m,g_lw):
        r""" It provides the surface nucleation rate :math:`\Pi_{0}` (see https://doi.org/10.1002/aic.690210502):
    
        .. math:: 
                \Pi_{0} = (N_A \rho)^{2/3} \left( \frac{1+m}{2} \right) \sqrt{\frac{2\gamma_{lw}}{\pi M}}
                
        where         
        * :math:`N_A` is the Avogadro number,
        * :math:`\rho` is the density of the liquid in which nucleation occurs,
        * :math:`\gamma_{lw}` is the interface tension of the surface where nucleation occurs,
        * :math:`M` is the mass of a molecule of the liquid 
        (that is :math:`M_{w}/N_{A}`, where :math:`M_{w}` is the molecular weigh of the molecule),
        * finally  :math:`m` is the ratio :math:`\frac{\gamma_{lw}-\gamma_{gw}}{\gamma_{gl}} = \cos\theta`.

        :param m:  It is the ratio :math:`m = \frac{\gamma_{lw}-\gamma_{gw}}{\gamma_{gl}} = \cos\theta`.
        :type command: float
        
        :param g_lw: Interfacial tension :math:`\gamma_{lw}` between the PFC liquid and the water inside the droplet, in N/m. 
        :type command: float    
        
        :return: The surface nucleation rate :math:`\Pi_{0}`
        :rtype: float        
        
        """
        value =  math.pow(scipy.constants.N_A*self.PFC.rho,2/3) * (1 + m)*0.5*math.sqrt(2*g_lw/(self.PFC.Mw/scipy.constants.N_A*math.pi))
        return value
    
    def Wdiff(self,P_out,R,n,g_lw,g_gw,g_gl):
        r""" It derives the value of the difference :math:`(W^{*}-W_{0.5})` that is equal to:
    
        .. math::
            \frac{16\pi\gamma^{3}_{gl}}{3P_{0.5}^{2}} f(m,x) - k_{B}T \ln\left(\frac{\Pi_{0}Af}{2\ln(2)}\right)
    
        where
        
        *  :math:`P_{0.5}` is the acoustic pressure,
        * :math:`k_{B}` is the Boltzmann constant, 
        * :math:`f` is the acoustic frequency,
        * :math:`T` is the temperature,
        * :math:`A` is the droplet surface (:math:`= 4\pi R^{2}` where :math:`R` is the droplet radius),
        * :math:`\gamma_{gl}`, is the interfacial tension between gaseous PFC and liquid PFC,
        * :math:`\gamma_{lw}`, is the interfacial tension between liquid PFC liquid and water,
        * :math:`\gamma_{gw}`, is the interfacial tension between gaseous PFC liquid and water,
        * :math:`m = \frac{\gamma_{lw}-\gamma_{gw}}{\gamma_{gl}}`
        * :math:`\Pi_{0}` is the surface nucleation rate.
       
        :param P_out: Acoustic pressure :math:`P_{0.5}`, in Pa.
        :type command: float        
    
        :param R: Radius :math:`R` in meter of a droplet on which surface nucleation occurs, in m
        :type command: float

        :param n: Number :math:`n` of droplets. 
        :type command: float    
        
        :param g_lw: Interfacial tension :math:`\gamma_{lw}` between the PFC liquid and the water inside the droplet, in N/m. 
        :type command: float    
        
        :param g_gw: Interfacial tension :math:`\gamma_{gw}` between the PFC gas inside the nuclei and the water droplets, in N/m. 
        :type command: float    

        :param g_gl: Interfacial tension :math:`\gamma_{lw}` between the PFC gas inside the nuclei and the PFC liquid, in N/m.
        :type command: float
        
        :return:  the difference :math:`(W^{*}-W_{0.5})`
        :rtype: float     
       
        """
        P = pressure.P_in(P_out,R,g_lw)
        m = (g_lw-g_gw)/g_gl
        rc = 2*g_gl/P
        x = R/rc
        A = 4*math.pi*numpy.power(R,2)
        value = 16*math.pi*numpy.power(g_gl,3)*self.f(x,m)/(3*math.pow(P,2)) - scipy.constants.k*self.T*numpy.log((self.Pi_0(m,g_lw)*n*self.freq*A)/(2*math.log(2)))
        return value
    
    def f(self,x,m):
        """ To be define in a subclass
        """
        print('function f is not defined')
        return False

    def func_P(self,P_out,R,n,g_lw, g_gw,g_gl):
        r""" This function is used to find the root of the equation :math:`(W^{*}-W_{0.5})`  derived by the function Wdiff.
        The unknown is P_out, and the other parameters R,n,g_gl,g_diff are given as args.
        
        :param P_out: Acoustic pressure :math:`P_{0.5}`, in Pa.
        :type command: float        
    
        :param R: Radius :math:`R` in meter of a droplet on which surface nucleation occurs, in m
        :type command: float

        :param n: Number :math:`n` of droplets. 
        :type command: float    
        
        :param g_lw: Interfacial tension :math:`\gamma_{lw}` between the PFC liquid and the water inside the droplet, in N/m. 
        :type command: float    
        
        :param g_gw: Interfacial tension :math:`\gamma_{gw}` between the PFC gas inside the nuclei and the water droplets, in N/m. 
        :type command: float    

        :param g_gl: Interfacial tension :math:`\gamma_{lw}` between the PFC gas inside the nuclei and the PFC liquid, in N/m.
        :type command: float
        
        :return:  the difference :math:`(W^{*}-W_{0.5})`
        :rtype: float           
        """
        value = self.Wdiff(P_out,R,n,g_lw, g_gw,g_gl)
        return value

    def func_3g(self,g,P,R,n):
        r""" This function is used to find the root of the equation :math:`(W^{*}-W_{0.5})`  derived by the function Wdiff.
        The unknown is the vector g = [g_lw,g_gw,g_gl], and the other parameters P,R,n are given as args.
        
        :param P: Acoustic pressure :math:`P_{0.5}`, in Pa.
        :type command: float        
    
        :param R: Radius :math:`R` in meter of a droplet on which surface nucleation occurs, in m
        :type command: float

        :param n: Number :math:`n` of droplets. 
        :type command: float    
        
        :param g: Vector of interfacial tensions [:math:`\gamma_{lw}`, :math:`\gamma_{gw}`, :math:`\gamma_{lw}`], in N/m.
        :type command: float
        
        :return:  the vector [:math:`(W^{*}-W_{0.5})`,0,0].
        :rtype: float           
        """
        g_lw = g[0]
        g_gw = g[1]
        g_gl = g[2]
        return numpy.array([self.Wdiff(P,R,n,g_lw,g_gw,g_gl),0,0])    

    def func_3g_b(self,g,P,R,n):
        r""" This function is used to find the root of the equation :math:`(W^{*}-W_{0.5})`  derived by the function Wdiff.
        The unknown is the vector g = [g_lw,g_gw,g_gl], and the other parameters P,R,n are given as args.
        
        :param P: Acoustic pressure :math:`P_{0.5}`, in Pa.
        :type command: float        
    
        :param R: Radius :math:`R` in meter of a droplet on which surface nucleation occurs, in m
        :type command: float

        :param n: Number :math:`n` of droplets. 
        :type command: float    
        
        :param g: Vector of interfacial tensions [:math:`\gamma_{lw}`, :math:`\gamma_{gw}`, :math:`\gamma_{lw}`], in N/m.
        :type command: float
        
        :return:  the difference :math:`(W^{*}-W_{0.5})`.
        :rtype: float           
        """
        g_lw = g[0]
        g_gw = g[1]
        g_gl = g[2]
        return self.Wdiff(P,R,n,g_lw,g_gw,g_gl)
###############################################################################


###############################################################################
class convexSurface(curvedSurface):
       
    def v(self,r,R,m):
        r""" It derives the volume :math:`v`  of a gas nuclei:
        
        .. math::
            v = \frac{\pi}{3} \left[r^{3}(2 - 3 \cos\psi +  \cos^{3}\psi) - R_{w}^{3} (2 - 3\cos\phi + \cos^{3}\phi)  \right]
            
        Where :math:`R_{w}` and :math:`r` are the radius of the droplet surface at which the nuclei occurs
        and of the surface of the nuclei volume not in contact with the droplet, respectively,
        the angles :math:`\phi` and :math:`\psi` are defined by         

        .. math::
            \cos\psi = \frac{mR_{w}-r}{d}
            
        .. math::
            \cos\phi = \frac{R_{w}-mr}{d}
            
        where :math:`d` is the distance of the center of the two curved surfaces (droplet and nuclei) and
        :math:`m =\cos\theta` where :math:`\theta` is the contact angle between the nuclei and droplet surface.
        
        :param r: Radius :math:`r` of the nuclei, in m.
        :type command: float

        :param R: Radius :math:`R` of the droplet, in m.
        :type command: float

        :param m: :math:`m = \frac{\gamma_{lw}-\gamma_{gw}}{\gamma_{gl}}` ,without units.
        :type command: float
        
        :return: Volume :math:`v` of a nuclei, in m\ :sup:`3`.
        :rtype: float     
        
        """
        d = numpy.sqrt(numpy.power(R,2) + numpy.power(r,2) - 2*R*r*m)
        cos_phi = (m*R-r)/d
        cos_psi = (R-m*r)/d
        value = math.pi/3*(numpy.power(r,3)*(2-3*cos_psi+numpy.power(cos_psi,3)) - numpy.power(R,3)*(2-3*cos_phi+numpy.power(cos_phi,3)))
        return value
    
    def f(self,x,m):
        r""" It derives the function :math:`f` (see https://doi.org/10.1063/1.3146810):
        
        .. math::
            f(m,x) = \frac{1}{2}\left\lbrace 1- \left( \frac{mx-1}{g} \right)^{3} + x^{3} \left[2-3\left(\frac{x-m}{g}\right) + \left(\frac{x-m}{g}\right)^{3}\right] +3mx^{2}\left(\frac{x-m}{g}-1\right)\right\rbrace
        
        where :math:`g = \sqrt{1 + x^{2}-2mx}`.
        
        :param x: Parameter :math:`x = \frac{R}{r^{*}}`, without unit.
        :type command: x

        :param m: Parameter :math:`m = \frac{\gamma_{lw}-\gamma_{gw}}{\gamma_{gl}}`, without unit.
        :type command: float
        
        :return: :math:`f`, without unit.
        :rtype: float           
        """
        g = numpy.sqrt(1 + numpy.power(x,2) - 2*m*x)
        t  =  (x - m)/g
        value = 0.5*(1 - numpy.power((m*x - 1)/g,3) + numpy.power(x,3)*(2 - 3*t + numpy.power(t,3) ) + 3*m*numpy.power(x,2)*(t-1))
        return value
    
    def func3(self,g):
        P_ex_single = 1.5e6
        Rw_single = 10.5e-6  
        n_single = 1 

        P_ex_multi = 2.2e6    
        Rw_multi = 0.215e-6 
        n_multi = math.pow(20e-6,3)/math.pow(Rw_multi,3)*0.4
        R =  numpy.array([Rw_single,Rw_multi])      
        P = numpy.array([P_ex_single,P_ex_multi])
        n = numpy.array([n_single,n_multi])
        g_lw = g[0]
        g_gw = g[1]
        g_gl = g[2]
        N = 2
        value = numpy.zeros(N+1)
        for i in range(0,N):
            value[i] = self.Wdiff(P[i],R[i],n[i],g_lw,g_gw,g_gl)    
        return value
    
###############################################################################

###############################################################################
class concaveSurface(curvedSurface):
    
    def v(self,r,R):
        r""" It derives the volume :math:`v`  of a gas nuclei 
        
        .. math::
            v = \frac{\pi}{3} \left[ R^{3} (2 - 3\cos\phi + \cos^{3}\phi) + r^{3}(2 - 3 \cos\psi +  \cos^{3}\psi) \right]
            
        Where :math:`R` and :math:`r` are the radius of the droplet and  of the nuclei, respectively, and 

        .. math::
            \cos\psi = \frac{Rm+r}{d}
            
        .. math::
            \cos\phi = \frac{R+rm}{d}        
        
        :param r: Radius :math:`r` of the nuclei, in m.
        :type command: float

        :param R: Radius :math:`R` of the droplet, in m.
        :type command: float
        
        :return: Volume :math:`v` of a nuclei, in m\ :sup:`3`.
        :rtype: float     
        
        """
        d = numpy.sqrt(numpy.power(R,2) + numpy.power(r,2)+ 2*R*r*self.m)
        m_psi = (R*self.m + r)/d
        m_phi = (R + r*self.m)/d  
        value = math.pi/3*(numpy.power(R,3)*(2 - 3*m_phi + numpy.power(m_phi,3)) + numpy.power(r,3)*(2 - 3*m_psi + numpy.power(m_psi,3)))
        return value
    
    def f(self,x,m):
        r""" It derives the function :math:`f`
        
        .. math:: 
            f(m,x) = \frac{1}{2}\left\lbrace 1- \left( \frac{1+mx}{g} \right)^{3} - x^{3} \left[2-3\left(\frac{x+m}{g}\right) + \left(\frac{x+m}{g}\right)^{3}\right] -3mx^{2}\left(1-\frac{x+m}{g}\right)\right\rbrace 
        
        
        with :math:`g = \sqrt{1 + x^{2} + 2mx}`

        :param x: Parameter :math:`x = \frac{R}{r^{*}}`, without unit.
        :type command: x

        :param m: Parameter :math:`m = \frac{\gamma_{lw}-\gamma_{gw}}{\gamma_{gl}}`, without unit.
        :type command: float
        
        :return: :math:`f`, without unit.
        :rtype: float   
        
        """
        g = numpy.sqrt(1 + numpy.power(x,2) + 2*m*x)
        t  =  (x + m)/g
        value = 0.5*(1 - numpy.power((1 + m*x)/g,3) - numpy.power(x,3)*(2 - 3*t + numpy.power(t,3) ) - 3*m*numpy.power(x,2)*(1-t))
        #print('value',value)
        return value
    
    def func3(self,g):
        r""" It derive the pressure vaporisation :math:`P` at probability :math:`p = 0.5` by finding the root of the function
    
        .. math::
            \frac{16\pi\gamma^{3}_{gl}}{3P_{0.5}^{2}} f(m,x) - k_{B}T \ln\left(\frac{\Pi_{0}A\tau}{\ln(2)}\right) = 0
    
    
        :param P: Acoustic pressure :math:`P`, in Pa.
        :type command: float    
    
        :param T: Temperature :math:`T`, in K.
        :type command: float

        :param R:  Radius :math:`R` in meter of a droplet, in m.
        :type command: float

        :param f: Acoustic frequency :math:`f`, in Hz.
        :type command: float    

        :param A: Surface area :math:`A` of a droplet , in m\ :sup:`2`.
        :type command: float    

        :param n: number of droplets :math:`n`
        :type command: float    
       
        :return: Root value
        :rtype: float     
       
        """
        R =  numpy.array([0.1,2.5,5,20,30])*1e-6      
        P = numpy.array([2.49,2.39,2.35,2.08,1.34])*1e6
        n = numpy.array([2.4e10,1.5e6,1.9e5,2.9e3,8.8e2])
        g_gw = g[0]
        g_lw = g[1]
        g_gl = g[2]
        N = len(P)
        value = numpy.zeros(N)
        for i in range(0,N):
            value[i] = self.Wdiff(P[i],R[i],n[i],g_lw,g_gw,g_gl)    
        return value
###############################################################################
