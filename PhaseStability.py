#!/usr/bin/env python

from scipy import *
import scipy.optimize

def Tieline(x, ModelInstance, T, c, M):
    slope = (ModelInstance.deltaGmix(x[0], T, c, M) - ModelInstance.deltaGmix(x[1], T, c, M))/(x[0]-x[1])
    return slope
        
def System(x, ModelInstance, T, c, M):
    SlopeGmix = ModelInstance.FirstDerivative(x[0], T, c, M)-ModelInstance.FirstDerivative(x[1], T, c, M)
    Tangent = Tieline(x, ModelInstance, T, c, M)-ModelInstance.FirstDerivative(x[0], T, c, M)
    system = array([SlopeGmix, Tangent])
    return system
    

def NonEqConstr(x, ModelInstance, T, c, M):
   
    TangentTest = array([((ModelInstance.FirstDerivative(x[0], T, c,M))*(comp - x[0])+ ModelInstance.deltaGmix(x[0], T, c, M))-ModelInstance.deltaGmix(comp, T, c, M) for comp in arange(0, 1.01, 0.01)])
    DistinctTest = -1*array([abs(x[0]-x[1])-0.1])
    return append(TangentTest, DistinctTest)
       

#===================================================================
       
def CalcPhaseStability(ModelInstance, T, c, M):
   
    
    dtype = [('x_Pred1', float), ('x_Pred2',float), ('Equal_Slopes', float),('Tangent', float) ]
    OutArray = array([0, 1, 100, 100])

    for x1 in append(0.5*rand(20), 0.001):
        for x2 in append(1-0.5*rand(20), 0.999):
            
            xo = array([x1, x2])
            (Pred_x_eq, infodict, ier, msg) = scipy.optimize.fsolve(System, xo,(ModelInstance, T, c, M), None, 1, 0, 10**-6, 0, None, 0.0, 100, None,True)
            
            Test = NonEqConstr(Pred_x_eq, ModelInstance, T, c, M)
            if all(Test<0) and ier==1 and all(Pred_x_eq>=0) and all(Pred_x_eq<=1) and (abs(infodict['fvec'][0])<OutArray[2] and abs(infodict['fvec'][1])<OutArray[3]):
                OutArray = array([Pred_x_eq[0], Pred_x_eq[1], abs(infodict['fvec'][0]),abs(infodict['fvec'][1])])
        
    
    BestFit = array([tuple(OutArray)], dtype)
    return array([BestFit[0][0], BestFit[0][1]])
