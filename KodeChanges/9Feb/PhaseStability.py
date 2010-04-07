#!/usr/bin/env python

from scipy import *
import scipy.optimize

def Tieline(x, ModelInstance, T, c, M):
            slope = (ModelInstance.deltaGmix(x[0], T, c, M) - ModelInstance.deltaGmix(x[1], T, c, M))/(x[0]-x[1])
            return slope
        
def System(x, ModelInstance, T, c, M):
            SlopeGmix = ModelInstance.FirstDerivative(x[0], T, c, M)-ModelInstance.FirstDerivative(x[1], T, c, M)
            Tangent = Tieline(x, ModelInstance, T, c, M)-ModelInstance.FirstDerivative(x[0], T, c, M)
            system = append(SlopeGmix, Tangent)

            return system

#===================================================================
       
def CalcPhaseStability(ModelInstance, T, c, M):
    xo = array([0.01, 0.99])
    Pred_x_eq = scipy.optimize.fsolve(System, xo,(ModelInstance, T, c, M), None, 0, 0, 10**-4, 0, None, 0.0, 100, None,True)
  
    return Pred_x_eq
        
    




