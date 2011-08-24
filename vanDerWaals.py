#!/usr/bin/env python

from scipy import *
import scipy.optimize
import pickle
import glob

def vdWaals(v, bvdW,a, T, R):
    PvdW = (R*T/(v - bvdW)) - a/(v**2)
    return PvdW
    
def vdWaalsIntegral(vl,vg, bvdW,a, T, R, PvdW):
    Integral = abs(R*T*log((vg-bvdW)/(vl -bvdW)) + a*((1/vg)- (1/vl))- PvdW*(vg-vl))
    return Integral

def vdWaalsPredict(m, Tc, Pc, T, P_act, R):
    b_vdW   = (R*Tc)/(8*Pc*100)
    ac_vdW  = (27*R**2*Tc**2)/(64*Pc*100)
    Tr  = T/Tc
    a = ac_vdW*exp(m*(1-Tr))

    k1 = 1
    k2 = -1*((P_act*b_vdW)+ (R*T))/P_act 
    k3 = a/P_act
    k4 = -1*(a*b_vdW)/P_act

    vdW_poly = array([k1, k2, k3, k4])
    vdW_roots = sort(roots(vdW_poly))

    vl = vdW_roots[0]
    vi = vdW_roots[1]
    vg = vdW_roots[2]
    
    if all(isreal(array([vl, vi, vg]))):
        P_prd = vdWaals(vl,b_vdW,a, T, R)
        pointerror = (P_prd -P_act)**2
    else:
       P_prd = NaN
       pointerror = NaN
             
    return array([P_prd, pointerror])

def mError(m, Tc, Pc, T, P, R):
    b_vdW   = (R*Tc)/(8*Pc*100)
    ac_vdW  = (27*R**2*Tc**2)/(64*Pc*100)
    Tr  = T/Tc
    a = ac_vdW*exp(m*(1-Tr))

    k1 = 1
    k2 = -1*((P*b_vdW)+ (R*T))/P 
    k3 = a/P
    k4 = -1*(a*b_vdW)/P

    vdW_poly = ([k1, k2, k3, k4])
    vdW_roots = sort(roots(vdW_poly))

    vl = vdW_roots[0]
    vi = vdW_roots[1]
    vg = vdW_roots[2]


    if all(isreal(array([vl, vi, vg]))): 
        P_vdW = vdWaals(vl, b_vdW,a, T, R)
        merror = abs(vdWaalsIntegral(vl,vg, b_vdW,a, T, R, P_vdW))

    else:
        merror = NaN

    return merror

def sqrError(m, Tc, Pc, Texp, Pexp, R,Tint):
    ErrorArea = 0.
    for count in arange(1,Tint):
        T = Texp[count]
        Predictions = vdWaalsPredict(m, Tc, Pc, T, interp(T,Texp,Pexp), R) 
        Pprd = Predictions[0]
        PointError = Predictions[1]
        if isreal(Pprd):
           ErrorArea = ErrorArea+ PointError*(Texp[1]-Texp[0])
   
    return ErrorArea
    

class VPData:
    
    def __init__(self, C, Tint):
        self.C = C
        self.TRange = linspace(C['Tmin'], C['Tmax'], Tint)

class Method1(VPData):
    
    def Generate(self):
        xRange =[1 - (Tpoint/self.C['Tc']*1.0) for Tpoint in self.TRange]
        logPvpRange =[(1/(1-x))*(self.C['VPa']*x + self.C['VPb']*x**1.5 + self.C['VPc']*x**3 + self.C['VPd']*x**6) for x in xRange]
        PexpRange = [100*self.C['Pc']*exp(logPvp) for logPvp in logPvpRange]
        
        return array([self.TRange, PexpRange])

class Method3(VPData):
    
    def Generate(self):
        logPvpRange =[self.C['VPa'] - self.C['VPb']/(Tpoint + self.C['VPc']) for Tpoint in self.TRange]
        PexpRange = [100*exp(logPvp) for logPvp in logPvpRange]
        
        return array([self.TRange, PexpRange])

class Slope:
      
    def __init__(self, Data, Compounds):
        self.Compounds = Compounds
        self.Data = Data
        

        
    def BestFit(self):
        Tint = 500
        R = 8.314
        self.VPData ={}
        self.Slope = dict((Compound, {}) for Compound in self.Compounds)
        for Compound in self.Compounds:
            C = self.Data[Compound]
            self.VPData[Compound] = eval('Method'+ str(C['Method'])+'(C, Tint)'+'.Generate()')
            self.Slope[Compound]['Range']= [scipy.optimize.fminbound(mError,0.,5.,(C['Tc'],C['Pc'], T, interp(T,self.VPData[Compound][0],self.VPData[Compound][1]),R),1e-4,500,0,0) for T in self.VPData[Compound][0]]
            self.Slope[Compound]['BestFit']= scipy.optimize.fminbound(sqrError,min(self.Slope[Compound]['Range']),max(self.Slope[Compound]['Range']),(C['Tc'],C['Pc'], self.VPData[Compound][0],self.VPData[Compound][1],R,Tint),1e-4,500,0,0) 
        
        return [self.Slope, self.VPData]

    def CompC(self, T):
        R = 8.314
        CompC = {}
        for Compound in self.Compounds:
            ac_vdW  = (27*R**2*self.Data[Compound]['Tc']**2)/(64*self.Data[Compound]['Pc']*100)
            b_vdW   = (R*self.Data[Compound]['Tc'])/(8*self.Data[Compound]['Pc']*100)
            Tr  = T/self.Data[Compound]['Tc']
            CompC[Compound] = -1*(ac_vdW/(R*T*b_vdW))*exp(self.Slope[Compound]['BestFit']*(1-(T/self.Data[Compound]['Tc'])))

        return CompC
