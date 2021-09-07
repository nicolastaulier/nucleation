
###############################################################################
class PFC():
    r""" This class provides the properties of several perfluorocarbon liquids. 
    Namely, the density :math:`\rho` (in liquid phase), 
    the surface tension :math:`\gamma` (between PCF and air),
    the molecular weight :math:`M_{w}`,
    and the boiling temperature :math:`T_{b}`.

    :param abbr: perfluorocarbon to be considered, could be PFH
    :type command: string

    """
    
    def  __init__(self,abbr):
        r""" The values for PFH are 
        :math:`\rho = 1680` kg/m\ :sup:`3` , 
        :math:`\gamma = 0.056` N/m,  
        :math:`M_{w}=0.338` kg/mol, 
        and :math:`T_{b} = 330.15` K.
        """
        if abbr == 'PFH':
            self.rho = 1680 # in kg/m:sup:`3` at 25°C
            self.gamma_0 = 0.056 # in N/m
            self.Mw = 0.338 # in kg/mol
            self.Tb = 59 + 271.15 # in K
###############################################################################

