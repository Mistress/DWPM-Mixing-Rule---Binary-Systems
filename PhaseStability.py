#!/usr/bin/env python

from scipy import *
from openopt import *

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
    return append(TangentTest, DistinctTest)
       

#===================================================================
       
def CalcPhaseStability(ModelInstance, T, c, M):
    xo = array([0.01,0.5])
    

    problem = NLSP(lambda x:System(x, ModelInstance, T, c, M), xo, h = lambda x: Tieline(x, ModelInstance, T, c, M), c = lambda x: NonEqConstr(x, ModelInstance, T, c, M), lb = zeros(2), ub = ones(2), ftol = 10^-4, xtol =10^-4)
    #problem.args = ( ModelInstance, T, c, M)
    #solution = problem.solve('nlp:ralg')
    solution = problem.solve('nssolve')
    
    return (solution.xf, solution.xf)
        
    





