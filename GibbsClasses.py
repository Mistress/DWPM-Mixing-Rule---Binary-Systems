#!/usr/bin/env python

from scipy import *
import scipy.optimize
import GibbsFunctions



class Model:
    R = 8.314
    name = 'Generic'
       
class DWPM(Model):
    def __init__(self, param):
        self.name = 'DWPM'
        self.c12 = param[0]
        self.c21 = param[1]
        self.s1 = param[2]
        self.s2 = param[3]
        
                         
    def deltaGmix(self, x, T, c, M):
        
        delGmix = GibbsFunctions.gibbsfunctions.dwpm(x, self.s1, self.s2, self.c12, self.c21, c)
                   
        return delGmix
            
    def FirstDerivative(self, x, T, c, M):
        
        dGmix_dx = GibbsFunctions.gibbsfunctions.d_dwpm_dx(x, self.s1, self.s2, self.c12, self.c21, c)
        
        return dGmix_dx
    
class NRTL(Model):
    def __init__(self, param):
        self.name = 'NRTL'
        self.A12 = self.R*param[0]
        self.A21 = self.R*param[1]
          
    
    def deltaGmix(self, x, T, c, M):

        delGmix = GibbsFunctions.gibbsfunctions.nrtl(x, self.A12, self.A21, T, self.R)
                    
        return delGmix

    def FirstDerivative(self, x, T, c, M):        
        
        dGmix_dx = GibbsFunctions.gibbsfunctions.d_nrtl_dx(x, self.A12, self.A21, T, self.R)

        return dGmix_dx
          
    
class UNIQUAC(Model):
    def __init__(self, param):
        self.name = 'UNIQUAC'        
        self.A12 =  self.R*param[0]
        self.A21 =  self.R*param[1]
    
    def deltaGmix(self, x, T, c, M):
        
        r1 = M['r1']
        r2 = M['r2']
        q1 = M['q1']
        q2 = M['q2']

        delGmix = GibbsFunctions.gibbsfunctions.uniquac(x, T, self.R, self.A12, self.A21, r1, r2, q1, q2)
                
        return delGmix
        
    def FirstDerivative(self, x, T, c, M):
        
        r1 = M['r1']
        r2 = M['r2']
        q1 = M['q1']
        q2 = M['q2']

        dGmix_dx = GibbsFunctions.gibbsfunctions.d_uniquac_dx(x, T, self.R, self.A12, self.A21, r1, r2, q1, q2)
       
        return dGmix_dx
