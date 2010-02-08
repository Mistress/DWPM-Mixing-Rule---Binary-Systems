#!/usr/bin/env python

from scipy import *

def DWPM(x1, param, R, T, c):
    
    c12 = param[0][0]
    c21 = param[0][1]
    s = param[0][2]
    A12 = c[0]/c12
    A21 = c[1]/c21
    x2 =1-x1

    if x1 == 0:
       delGmix = 0
    elif x1 == 1:
       delGmix = 0
    else:
       delGmix = real(x1*log(x1) + x2*log(x2) - (1/s)*(x1*log(x1 + x2*A12**s) + x2*log(x1*A21**s + x2)))
    
    dGmix_dx = real(-(-log(x1*A21**s-x1+1)+(1-x1)*(A21**s-1)/(x1*A21**s-x1+1) + log((1-x1)*A12**s+x1)+x1*(1-A12**s)/((1-x1)*A12**s+x1))/s +log(x1)-log(1-x1))
    
    return (delGmix, dGmix_dx)


def NRTL(x1,param,R,T):

    A12 = R*param[1][0]
    A21 = R*param[1][1]
    Tau12 = A12/(R*T)
    Tau21 = A21/(R*T)
    alpha = 0.2
    G12 = exp(-alpha*Tau12)
    G21 = exp(-alpha*Tau21)
    x2 = 1 - x1
          
    if x1 == 0:
        delGmix = 0
    elif x1 == 1:
        delGmix == 0
    else:
        delGmix = real(x1*x2*(((Tau21*G21)/(x1+(x2*G21))) + ((Tau12*G12)/(x2+(x1*G12)))) + x1*log(x1) + x2*log(x2))

    dGmix_dx = real((x1*x2)*(-1*Tau21*G21*((x1 + x2*G21)**-2)*(1-G21) - Tau12*G12*((x2 + x1*G12)**-2)*(G12-1)) + (1-2*x1)*((Tau21*G21)/(x1+(x2*G21)) + (Tau12*G12)/(x2+(x1*G12))) +log(x1)-log(x2))
          
    return (delGmix, dGmix_dx)

def UNIQUAC(x1,param,R,T,M):
    
    A12 = R*param[2][0]
    A21 = R*param[2][1]
    x2 = 1-x1
    r1 = M['r1']
    r2 = M['r2']
    q1 = M['q1']
    q2 = M['q2']

    theta1 = (x1*q1)/(x1*q1+x2*q2)
    theta2 = (x2*q2)/(x1*q1+x2*q2)
    phi1 = (x1*r1)/(x1*r1 +x2*r2)
    phi2 = (x2*r2)/(x1*r1 +x2*r2)
    Tau_12 = exp(-A12/(R*T))
    Tau_21 = exp(-A21/(R*T))
    z = 10;

    Gr = -q1*x1*log(theta1 + theta2*Tau_21) - q2*x2*log(theta2 + theta1*Tau_12)
    Gc = x1*log(phi1/x1)+ x2*log(phi2/x2) + (z/2)*(q1*x1*log(theta1/phi1) + q2*x2*log(theta2/phi2))

    if x1 == 0:
      delGmix = 0
    elif x1 == 1:
      delGmix = 0
    else:
      delGmix = real(Gr + Gc + x1*log(x1) + x2*log(x2))

    A = -q1*log((q1*x1+q2*Tau_21*x2)/(q1*x1+q2*x2))- q1*x1*((q1-q2*Tau_21)/(q1*x1+q2*x2)-((q1*(q1-q2)*x1)+((q1-q2)*q2*Tau_21*x2))/(q1*x1+q2*x2)**2)/((q1*x1+ q2*Tau_21*x2)/(q1*x1+q2*x2))
    B = q2*log((q1*Tau_12*x1+q2*x2)/(q1*x1+q2*x2))-(q2*x2*(q1*x1+q2*x2)*((q1*Tau_12-q2)/(q1*x1+q2*x2)-((q1-q2)*(q1*Tau_12*x1+q2*x2))/(q1*x1+q2*x2)**2))/(q1*Tau_12*x1+q2*x2)
    C = -log(r2/(r1*x1+r2*x2))+log(r1/(r1*x1+r2*x2))-((r1-r2)*x1)/(r1*x1+r2*x2)-((r1-r2)*x2)/(r1*x1+r2*x2)

    D = q1*log(q1*(r1*x1+r2*x2)/(r1*(q1*x1+q2*x2)))+ r1*x1*(q1*x1+q2*x2)*(q1*(r1-r2)/(r1*(q1*x1+q2*x2))-(q1*(q1-q2)*(r1*x1+r2*x2))/(r1*(q1*x1+q2*x2)**2))/(r1*x1+r2*x2)
    E = (r2*x2*(q1*x1+q2*x2)*((q2*(r1-r2))/(r2*(q1*x1+q2*x2))-((q1-q2)*q2*(r1*x1+r2*x2))/(r2*(q1*x1+q2*x2)**2)))/(r1*x1+r2*x2)-q2*log((q2*(r1*x1+r2*x2))/(r2*(q1*x1+q2*x2)))

    dGmix_dx = real(A + B + C +(z/2)*(D+E) + log(x1)- log(x2));

    return (delGmix, dGmix_dx)


param1 = ([-90.4, -16.4, 0.5], [700,173],[5.4,600])
param2 = ([-130, -44.1, 0.5], [475,480],[46.4,427])
T = 288
R = 8.314
M = {'r1':2.0086, 'r2':4.0464, 'q1':1.8680, 'q2':3.2400}
c = (-13.7086,  -11.8103)
x1 = 0.5

DWPMPython_Param1 = DWPM(x1, param1, R, T, c)
DWPMPython_Param2 = DWPM(x1, param2, R, T, c)
NRTLPython_Param1 = NRTL(x1, param1,R,T)
NRTLPython_Param2 = NRTL(x1, param2,R,T)
UNIQUACPython_Param1 = UNIQUAC(x1, param1,R,T,M)
UNIQUACPython_Param2 = UNIQUAC(x1, param2,R,T,M)

DWPMMatlab_Param1 = (-0.2502, -0.1440)
DWPMMatlab_Param2 = (-0.0051, -0.1119)
NRTLMatlab_Param1 = (-0.0892, 0.2037)
NRTLMatlab_Param2 = (-3.5786e-4, -0.0021)
UNIQUACMatlab_Param1 = (0.0474, -0.1931)
UNIQUACMatlab_Param2 = (-4.5623e-4, -0.0068)

print('%15s | %15s | %15s | %15s |%15s | %15s | %15s |%16s' %('','ParamSet', 'GibbsDWPM','dGibbs/dxDWPM','GibbsNRTL','dGibbs/dxNRTL','GibbsUNIQUAC','dGibbs/dxUNIQUAC'))
print('%15s | %15i | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f | %16.4f' %(('Matlab', 1)+DWPMMatlab_Param1+ NRTLMatlab_Param1+ UNIQUACMatlab_Param1))
print('%15s | %15i | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f | %16.4f' %(('Python', 1)+DWPMPython_Param1+ NRTLPython_Param1+ UNIQUACPython_Param1))
print('%15s | %15i | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f | %16.4f' %(('Matlab', 2)+DWPMMatlab_Param2+ NRTLMatlab_Param2+ UNIQUACMatlab_Param2))
print('%15s | %15i | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f | %16.4f' %(('Python', 2)+DWPMPython_Param2+ NRTLPython_Param2+ UNIQUACPython_Param2))
