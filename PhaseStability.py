#!/usr/bin/env python

from scipy import *
import scipy.optimize

def Tieline(x, ModelInstance, T, c, M):
    slope = (ModelInstance.deltaGmix(x[0], T, c, M) - ModelInstance.deltaGmix(x[1], T, c, M))/(x[0]-x[1])
    Tangent = slope - ModelInstance.FirstDerivative(x[0], T, c, M)
    return array([Tangent])
        
def System(x, ModelInstance, T, c, M):
    SlopeGmix = ModelInstance.FirstDerivative(x[0], T, c, M)-ModelInstance.FirstDerivative(x[1], T, c, M)
    return SlopeGmix**2

def NonEqConstr(x, ModelInstance, T, c, M):
   
    TangentTest = array([((ModelInstance.FirstDerivative(x[0], T, c,M))*(comp - x[0])+ ModelInstance.deltaGmix(x[0], T, c, M))-ModelInstance.deltaGmix(comp, T, c, M) for comp in arange(0, 1.001, 0.001)])
    DistinctTest = -1*array([abs(x[0]-x[1])+0.1])
    return -append(TangentTest, DistinctTest)
       

#===================================================================
       
def CalcPhaseStability(ModelInstance, T, c, M):
    xo = array([0.01,0.99])
    bounds = [(0, 1), (0, 1)]

    [Pred_x_eq, fx, its, imode, smode] = scipy.optimize.fmin_slsqp(System, xo,[], Tieline, [], NonEqConstr, bounds, None, None, None, (ModelInstance, T, c, M), 1000, 10e-4, 1, 1, 10e-8)
    
    return (Pred_x_eq, imode, smode)
        
    





