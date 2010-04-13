#!/usr/bin/env python

from scipy import *
import scipy.optimize

def Tieline(x, ModelInstance, T, c, M):
    slope = (ModelInstance.deltaGmix(x[0], T, c, M) - ModelInstance.deltaGmix(x[1], T, c, M))/(x[0]-x[1])
    Tangent = slope - ModelInstance.FirstDerivative(x[0], T, c, M)
    return Tangent
        
def System(x, ModelInstance, T, c, M):
    SlopeGmix = ModelInstance.FirstDerivative(x[0], T, c, M)-ModelInstance.FirstDerivative(x[1], T, c, M)
    
    Tangent = Tieline(x, ModelInstance, T, c, M)-ModelInstance.FirstDerivative(x[0], T, c, M)
    system = append(SlopeGmix, Tangent)
    return system
    

def NonEqConstr(x, ModelInstance, T, c, M):
   
    TangentTest = array([((ModelInstance.FirstDerivative(x[0], T, c,M))*(comp - x[0])+ ModelInstance.deltaGmix(x[0], T, c, M))-ModelInstance.deltaGmix(comp, T, c, M) for comp in arange(0, 1.001, 0.001)])
    DistinctTest = -1*array([abs(x[0]-x[1])+0.1])
    return append(TangentTest, DistinctTest)
       

#===================================================================
       
def CalcPhaseStability(ModelInstance, T, c, M):
   
    
    dtype = [('x_Pred1', float), ('x_Pred2',float), ('f_value', float)]
    OutArray = 100*array([ones(3)])

    for x1 in arange(0, 1.05, 0.05):
        for x2 in arange(0, 1.05, 0.05):
            
            xo = array([x1, x2])
            (Pred_x_eq, infodict, ier, msg) = scipy.optimize.fsolve(System, xo,(ModelInstance, T, c, M), None, 1, 0, 10**-4, 0, None, 0.0, 100, None,True)
            
            Test = NonEqConstr(Pred_x_eq, ModelInstance, T, c, M)
            if all(Test<0) & ier==1 :
                OutArray = append(OutArray, [[Pred_x_eq[0], Pred_x_eq[1], sum(infodict['fvec']**2)]], 0)
        
    OutList = [tuple(item) for item in OutArray]
    StructOutArray = array(OutList, dtype)
    BestFit = sort(StructOutArray, order=['f_value'])

    return array([BestFit[0][0], BestFit[0][1]])
