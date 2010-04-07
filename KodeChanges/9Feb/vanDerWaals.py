#!/usr/bin/env python

from scipy import *
import scipy.optimize
import pickle
import glob

#=========================#Functions#==============================#

def GenVPDataMeth1(Tint, C):
    TRange = linspace(C['Tmin'],C['Tmax'],Tint)
    xRange =[1 - (Tpoint/C['Tc']*1.0) for Tpoint in TRange]
    logPvpRange =[(1/(1-x))*(C['VPa']*x + C['VPb']*x**1.5 + C['VPc']*x**3 + C['VPd']*x**6) for x in xRange]
    PexpRange = [100*C['Pc']*exp(logPvp) for logPvp in logPvpRange]
    return array([TRange, PexpRange])

def GenVPDataMeth3(Tint, C):
    TRange = linspace(C['Tmin'],C['Tmax'],Tint)
    logPvpRange =[C['VPa'] - C['VPb']/(Tpoint + C['VPc']) for Tpoint in TRange]
    PexpRange = [100*exp(logPvp) for logPvp in logPvpRange]
    return array([TRange, PexpRange])

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
    
#=========================#MainCode#===============================#

def MainSlope(Data, Compounds):
    Tint = 500
    R = 8.314
    VPData = {}
    Slope ={}
    for Compound in Compounds:
        C = Data[Compound]
        if C['Method'] == 1:
            VPData[Compound] = GenVPDataMeth1(Tint,C)
        elif C['Method'] == 3:
            VPData[Compound] = GenVPDataMeth3(Tint,C)
        Slope[Compound]={}
        Slope[Compound]['Range']= [scipy.optimize.fminbound(mError,0.,5.,(C['Tc'],C['Pc'], T, interp(T,VPData[Compound][0],VPData[Compound][1]),R),1e-4,500,0,0) for T in VPData[Compound][0]]
        Slope[Compound]['BestFit']= scipy.optimize.fminbound(sqrError,min(Slope[Compound]['Range']),max(Slope[Compound]['Range']),(C['Tc'],C['Pc'], VPData[Compound][0],VPData[Compound][1],R,Tint),1e-4,500,0,0) 
        
    return [Slope, VPData]



def CompC(Compounds, Data, Slopes, T, R):
    CompC = {}
    for Compound in Compounds:
        ac_vdW  = (27*R**2*Data[Compound]['Tc']**2)/(64*Data[Compound]['Pc']*100)
        b_vdW   = (R*Data[Compound]['Tc'])/(8*Data[Compound]['Pc']*100)
        Tr  = T/Data[Compound]['Tc']
        a = ac_vdW*exp(Slopes[Compound]*(1-Tr))
        CompC[Compound] = -1*(ac_vdW/(R*T*b_vdW))*exp(Slopes[Compound]*(1-(T/Data[Compound]['Tc'])))

    return CompC
