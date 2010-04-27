#!/usr/bin/env python

from scipy import *
import scipy.optimize



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
        x1 = x
        x2 =1-x1
        A12 = c[0]/self.c12
        A21 = c[1]/self.c21

        
        if x1 == 0:
            delGmix = 0
        elif x1 == 1:
            delGmix = 0
        else:
            delGmix = real(x1*log(x1) + x2*log(x2) - (x1/self.s1)*log(x1 + x2*A12**self.s1) - (x2/self.s2)*log(x1*A21**self.s2 + x2))
        return delGmix
            
    def FirstDerivative(self, x, T, c, M):
        x1 = x
        x2 =1-x1
        A12 = c[0]/self.c12
        A21 = c[1]/self.c21

        dGmix_dx = real((1/self.s2)*log(x1*A21**self.s2+x2)-(x2/self.s2)*(A21**self.s2-1)/(x1*A21**self.s2 +x2) -(1/self.s1)*log(x2*A12**self.s1+x1) - (x1/self.s1)*(1-A12**self.s1)/(x2*A12**self.s1+x1) +log(x1)-log(x2))
        return dGmix_dx
    
class NRTL(Model):
    def __init__(self, param):
        self.name = 'NRTL'
        self.A12 = self.R*param[0]
        self.A21 = self.R*param[1]
          
    
    def deltaGmix(self, x, T, c, M):        
        x1 = x
        Tau12 = self.A12/(self.R*T)
        Tau21 = self.A21/(self.R*T)
        alpha = 0.2
        G12 = exp(-alpha*Tau12)
        G21 = exp(-alpha*Tau21)
        x2 = 1 - x1 
        if x1 == 0:
            delGmix = 0
        elif x1 == 1:
            delGmix = 0
        else:
            delGmix = real(x1*x2*(((Tau21*G21)/(x1+(x2*G21))) + ((Tau12*G12)/(x2+(x1*G12)))) + x1*log(x1) + x2*log(x2))
        return delGmix

    def FirstDerivative(self, x, T, c, M):        
        x1 = x
        Tau12 = self.A12/(self.R*T)
        Tau21 = self.A21/(self.R*T)
        alpha = 0.2
        G12 = exp(-alpha*Tau12)
        G21 = exp(-alpha*Tau21)
        x2 = 1 - x1 
        
        dGmix_dx = real((x1*x2)*(-1* Tau21*G21*((x1 + x2*G21)**-2)*(1- G21) - Tau12*G12*((x2 +x1*G12)**-2)*(G12-1)) + (1-2*x1)*((Tau21*G21)/(x1+(x2*G21)) + (Tau12* G12)/(x2+(x1*G12))) +log(x1)-log(x2))

        return dGmix_dx
          
    
class UNIQUAC(Model):
    def __init__(self, param):
        self.name = 'UNIQUAC'        
        self.A12 =  self.R*param[0]
        self.A21 =  self.R*param[1]
    
    def deltaGmix(self, x, T, c, M):
        x1 = x
        r1 = M['r1']
        r2 = M['r2']
        q1 = M['q1']
        q2 = M['q2']
        x2 = 1- x1
        theta1 = (x1*q1)/(x1*q1+x2*q2)
        theta2 = (x2*q2)/(x1*q1+x2*q2)
        phi1 = (x1*r1)/(x1*r1 +x2*r2)
        phi2 = (x2*r2)/(x1*r1 +x2*r2)
        Tau_12 = exp(-self.A12/(self.R*T))
        Tau_21 = exp(-self.A21/(self.R*T))
        z = 10

        Gr = -q1*x1*log(theta1 + theta2*Tau_21) - q2*x2*log(theta2 + theta1*Tau_12)
        Gc = x1*log(phi1/x1)+ x2*log(phi2/x2) + (z/2)*(q1*x1*log(theta1/phi1) + q2*x2*log(theta2/phi2))

        if x1 == 0:
            delGmix = 0
        elif x1 == 1:
            delGmix = 0
        else:
            delGmix = real(Gr + Gc + x1*log(x1) + x2*log(x2))
        return delGmix
        
    def FirstDerivative(self, x, T, c, M):
        x1 = x
        r1 = M['r1']
        r2 = M['r2']
        q1 = M['q1']
        q2 = M['q2']
        x2 = 1- x1
        theta1 = (x1*q1)/(x1*q1+x2*q2)
        theta2 = (x2*q2)/(x1*q1+x2*q2)
        phi1 = (x1*r1)/(x1*r1 +x2*r2)
        phi2 = (x2*r2)/(x1*r1 +x2*r2)
        Tau_12 = exp(-self.A12/(self.R*T))
        Tau_21 = exp(-self.A21/(self.R*T))
        z = 10

        A = -q1*log((q1*x1+q2*Tau_21*x2)/(q1*x1+q2*x2))- q1*x1*((q1-q2*Tau_21)/(q1*x1+q2*x2)-((q1*(q1-q2)*x1)+((q1-q2)*q2*Tau_21*x2))/(q1*x1+q2*x2)**2)/((q1*x1+ q2*Tau_21*x2)/(q1*x1+q2*x2))
        B = q2*log((q1*Tau_12*x1+q2*x2)/(q1*x1+q2*x2))-(q2*x2*(q1*x1+q2*x2)*((q1*Tau_12-q2)/(q1*x1+q2*x2)-((q1-q2)*(q1*Tau_12*x1+q2*x2))/(q1*x1+q2*x2)**2))/(q1*Tau_12*x1+q2*x2)
        C = -log(r2/(r1*x1+r2*x2))+log(r1/(r1*x1+r2*x2))-((r1-r2)*x1)/(r1*x1+r2*x2)-((r1-r2)*x2)/(r1*x1+r2*x2)

        D = q1*log(q1*(r1*x1+r2*x2)/(r1*(q1*x1+q2*x2)))+ r1*x1*(q1*x1+q2*x2)*(q1*(r1-r2)/(r1*(q1*x1+q2*x2))-(q1*(q1-q2)*(r1*x1+r2*x2))/(r1*(q1*x1+q2*x2)**2))/(r1*x1+r2*x2)
        E = (r2*x2*(q1*x1+q2*x2)*((q2*(r1-r2))/(r2*(q1*x1+q2*x2))-((q1-q2)*q2*(r1*x1+r2*x2))/(r2*(q1*x1+q2*x2)**2)))/(r1*x1+r2*x2)-q2*log((q2*(r1*x1+r2*x2))/(r2*(q1*x1+q2*x2)))
    
        dGmix_dx = real(A +B + C +(z/2)*(D+E) + log(x1)- log(x2))
        return dGmix_dx


#x = 0.5
#c = (-13.7086, -11.8103)
#M = {'r1':2.0086, 'r2':4.0464, 'q1':1.8680, 'q2':3.2400}
#T = 288
#param = (-90.4, -16.4, 0.5)
#a = DWPM(x, T, param, c, M) #Nie almal van die classes het M en c nodig nie maar as ek die "inputs" nie identies maak nie dan is daar probleme in die lyn van "init takes 4 inputs, 3 given" etc Ek kan seker net 'n if statement inbou wat die regte inputs stuur wanneer dit 'n spesifieke model gebruik... which might be better, aangesien daar nie onnodige "copies" van veranderlikes gemaak word of onnodig aangestuur word nie (memory effieciency gewys)
#gibbs = a.deltaGmix()
#ddxGmix = a.FirstDerivative()

#a = NRTL(x, T, param, c, M)
