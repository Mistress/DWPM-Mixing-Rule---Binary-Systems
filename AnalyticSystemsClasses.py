#!/usr/bin/env python

from scipy import *


class SystemDWPM:
    def __init__(self, Cell_s):
        self.s_1 = Cell_s[0]
        self.s_2 = Cell_s[1]

    def SystemEquations(self, Params, Actual, R, T):
        
        lambda_12 = Params[0]
        lambda_21 = Params[1]
        m = Params[2]
        b = Params[3]
       
        x = Actual[0]
        
        Eq1 = -(1-x)*log(lambda_21**self.s_2*x-x+1)/self.s_2-x*log(x+lambda_12**self.s_1*(1-x))/self.s_1+x*log(x)-m*x+log(1-x)*(1-x)-b
        Eq3 = log(lambda_21**self.s_2*x-x+1)/self.s_2-log(x+lambda_12**self.s_1*(1-x))/self.s_1+log(x)-(lambda_21**self.s_2-1)*(1-x)/(self.s_2*(lambda_21**self.s_2*x-x+1))\
            -(1-lambda_12**self.s_1)*x/(self.s_1*(x+lambda_12**self.s_1*(1-x)))-log(1-x)-m

        x = Actual[1]
        
        Eq2 =  -(1-x)*log(lambda_21**self.s_2*x-x+1)/self.s_2-x*log(x+lambda_12**self.s_1*(1-x))/self.s_1+x*log(x)-m*x+log(1-x)*(1-x)-b
        Eq4 = log(lambda_21**self.s_2*x-x+1)/self.s_2-log(x+lambda_12**self.s_1*(1-x))/self.s_1+log(x)-(lambda_21**self.s_2-1)*(1-x)/(self.s_2*(lambda_21**self.s_2*x-x+1))\
            -(1-lambda_12**self.s_1)*x/(self.s_1*(x+lambda_12**self.s_1*(1-x)))-log(1-x)-m

        return array([Eq1, Eq2, Eq3, Eq4])

    def SystemEquationsJac(self, Params, Actual, R, T):
        
        lambda_12 = Params[0]
        lambda_21 = Params[1]
        m = Params[2]
        b = Params[3]
        
        x = Actual[0]
                
        Row1 = [-lambda_12**(self.s_1-1)*(1-x)*x/(x+lambda_12**self.s_1*(1-x)), -lambda_21**(self.s_2-1)*(1-x)*x/(lambda_21**self.s_2*x-x+1), -x, -1.]
        Row3 = [lambda_12**(self.s_1-1)*x/(x+lambda_12**self.s_1*(1-x))-lambda_12**(self.s_1-1)*(1-x)/(x+lambda_12**self.s_1*(1-x))+\
        lambda_12**(self.s_1-1)*(1-lambda_12**self.s_1)*(1-x)*x/(x+lambda_12**self.s_1*(1-x))**2,\
        lambda_21**(self.s_2-1)*x/(lambda_21**self.s_2*x-x+1)-lambda_21**(self.s_2-1)*(1-x)/(lambda_21**self.s_2*x-x+1)+\
        lambda_21**(self.s_2-1)*(lambda_21**self.s_2-1)*(1-x)*x/(lambda_21**self.s_2*x-x+1)**2,-1., 0.]

        x = Actual[1]
        
        Row2 = [-lambda_12**(self.s_1-1)*(1-x)*x/(x+lambda_12**self.s_1*(1-x)), -lambda_21**(self.s_2-1)*(1-x)*x/(lambda_21**self.s_2*x-x+1), -x, -1.]
        Row4 = [lambda_12**(self.s_1-1)*x/(x+lambda_12**self.s_1*(1-x))-lambda_12**(self.s_1-1)*(1-x)/(x+lambda_12**self.s_1*(1-x))+\
        lambda_12**(self.s_1-1)*(1-lambda_12**self.s_1)*(1-x)*x/(x+lambda_12**self.s_1*(1-x))**2,\
        lambda_21**(self.s_2-1)*x/(lambda_21**self.s_2*x-x+1)-lambda_21**(self.s_2-1)*(1-x)/(lambda_21**self.s_2*x-x+1)+\
        lambda_21**(self.s_2-1)*(lambda_21**self.s_2-1)*(1-x)*x/(lambda_21**self.s_2*x-x+1)**2,-1., 0.]

        return array([Row1, Row2, Row3, Row4])

class SystemNRTL:
    def __init__(self, alpha):
        self.alpha = alpha
    
    def SystemEquations(self, Params, Actual, R, T):

        A12 = Params[0]*R
        A21 = Params[1]*R
        m = Params[2]
        b = Params[3]

        x = Actual[0]

        Eq1 = 1000*((1-x)*x*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x)*exp(-self.alpha*A21/(R*T\
        ))+x))+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x*exp(-self.alpha*A12/(R*T))-\
        x+1)))+x*log(x)-m*x+log(1-x)*(1-x)-b)
        
        Eq3 = 1000*(-x*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x)*exp(-self.alpha*A21/(R*T))+x)\
        )+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x*exp(-self.alpha*A12/(R*T))-x+1))\
        )+(1-x)*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x)*exp(-self.alpha*A21/(\
        R*T))+x))+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x*exp(-self.alpha*A12/(R*T\
        ))-x+1)))+(1-x)*x*(-A21*(1-exp(-self.alpha*A21/(R*T)))*exp(-self.alpha*A21\
        /(R*T))/(R*T*((1-x)*exp(-self.alpha*A21/(R*T))+x)**2)-A12*(exp(-self.alpha\
        *A12/(R*T))-1)*exp(-self.alpha*A12/(R*T))/(R*T*(x*exp(-self.alpha*A12/(\
        R*T))-x+1)**2))+log(x)-log(1-x)-m)

        x = Actual[1]

        Eq2 = 1000*((1-x)*x*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x)*exp(-self.alpha*A21/(R*T\
        ))+x))+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x*exp(-self.alpha*A12/(R*T))-\
        x+1)))+x*log(x)-m*x+log(1-x)*(1-x)-b)
        
        Eq4 = 1000*(-x*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x)*exp(-self.alpha*A21/(R*T))+x)\
        )+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x*exp(-self.alpha*A12/(R*T))-x+1))\
        )+(1-x)*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x)*exp(-self.alpha*A21/(\
        R*T))+x))+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x*exp(-self.alpha*A12/(R*T\
        ))-x+1)))+(1-x)*x*(-A21*(1-exp(-self.alpha*A21/(R*T)))*exp(-self.alpha*A21\
        /(R*T))/(R*T*((1-x)*exp(-self.alpha*A21/(R*T))+x)**2)-A12*(exp(-self.alpha\
        *A12/(R*T))-1)*exp(-self.alpha*A12/(R*T))/(R*T*(x*exp(-self.alpha*A12/(\
        R*T))-x+1)**2))+log(x)-log(1-x)-m)
           
        return array([Eq1, Eq2, Eq3, Eq4])

    def SystemEquationsJac(self, Params, Actual, R, T):

        A12 = Params[0]*R
        A21 = Params[1]*R
        m = Params[2]
        b = Params[3]

        x = Actual[0]

        Row1 = [1000*((1-x)*x*(exp(-self.alpha*A12/(R*T))/(R*T*(x*exp(-self.alpha*A12/(R*T))-x+1)\
        )-self.alpha*A12*exp(-self.alpha*A12/(R*T))/(R**2*T**2*(x*exp(-self.alpha*A12/\
        (R*T))-x+1))+self.alpha*x*A12*exp(-2*self.alpha*A12/(R*T))/(R**2*T**2*(x*\
        exp(-self.alpha*A12/(R*T))-x+1)**2))),1000*((1-x)*x*(exp(-self.alpha*A21/(R*T))/\
        (R*T*((1-x)*exp(-self.alpha*A21/(R*T))+x))-self.alpha*A21*exp(-self.alpha*A21/\
        (R*T))/(R**2*T**2*((1-x)*exp(-self.alpha*A21/(R*T))+x))+self.alpha*(1-x)*\
        A21*exp(-2*self.alpha*A21/(R*T))/(R**2*T**2*((1-x)*exp(-self.alpha*A21/(R\
        *T))+x)**2))),-1000*x,-1000]
 
        Row3 = [1000*(-x*(exp(-self.alpha*A12/(R*T))/(R*T*(x*exp(-self.alpha*A12/(R*T))-x+1))-self.alpha\
        *A12*exp(-self.alpha*A12/(R*T))/(R**2*T**2*(x*exp(-self.alpha*A12/(R*T)\
        )-x+1))+self.alpha*x*A12*exp(-2*self.alpha*A12/(R*T))/(R**2*T**2*(x*exp(-\
        self.alpha*A12/(R*T))-x+1)**2))+(1-x)*(exp(-self.alpha*A12/(R*T))/(R*T*(x\
        *exp(-self.alpha*A12/(R*T))-x+1))-self.alpha*A12*exp(-self.alpha*A12/(R*T))/(\
        R**2*T**2*(x*exp(-self.alpha*A12/(R*T))-x+1))+self.alpha*x*A12*exp(-2*self.alpha\
        *A12/(R*T))/(R**2*T**2*(x*exp(-self.alpha*A12/(R*T))-x+1)**2))+(1-x\
        )*x*(-(exp(-self.alpha*A12/(R*T))-1)*exp(-self.alpha*A12/(R*T))/(R*T*(x*\
        exp(-self.alpha*A12/(R*T))-x+1)**2)+self.alpha*A12*(exp(-self.alpha*A12/(R*T))-\
        1)*exp(-self.alpha*A12/(R*T))/(R**2*T**2*(x*exp(-self.alpha*A12/(R*T))-x+\
        1)**2)+self.alpha*A12*exp(-2*self.alpha*A12/(R*T))/(R**2*T**2*(x*exp(-self.alpha\
        *A12/(R*T))-x+1)**2)-2*self.alpha*x*A12*exp(-2*self.alpha*A12/(R*T))*(exp\
        (-self.alpha*A12/(R*T))-1)/(R**2*T**2*(x*exp(-self.alpha*A12/(R*T))-x+1\
        )**3))),1000*(-x*(exp(-self.alpha*A21/(R*T))/(R*T*((1-x)*exp(-self.alpha*A21/(R*\
        T))+x))-self.alpha*A21*exp(-self.alpha*A21/(R*T))/(R**2*T**2*((1-x)*exp(-\
        self.alpha*A21/(R*T))+x))+self.alpha*(1-x)*A21*exp(-2*self.alpha*A21/(R*T))/(\
        R**2*T**2*((1-x)*exp(-self.alpha*A21/(R*T))+x)**2))+(1-x)*(exp(-self.alpha\
        *A21/(R*T))/(R*T*((1-x)*exp(-self.alpha*A21/(R*T))+x))-self.alpha*A21*exp\
        (-self.alpha*A21/(R*T))/(R**2*T**2*((1-x)*exp(-self.alpha*A21/(R*T))+x))+\
        self.alpha*(1-x)*A21*exp(-2*self.alpha*A21/(R*T))/(R**2*T**2*((1-x)*exp(-\
        self.alpha*A21/(R*T))+x)**2))+(1-x)*x*(-(1-exp(-self.alpha*A21/(R*T)))*exp\
        (-self.alpha*A21/(R*T))/(R*T*((1-x)*exp(-self.alpha*A21/(R*T))+x)**2)+self.alpha\
        *A21*(1-exp(-self.alpha*A21/(R*T)))*exp(-self.alpha*A21/(R*T))/(R**2*\
        T**2*((1-x)*exp(-self.alpha*A21/(R*T))+x)**2)-self.alpha*A21*exp(-2*self.alpha*\
        A21/(R*T))/(R**2*T**2*((1-x)*exp(-self.alpha*A21/(R*T))+x)**2)-2*self.alpha\
        *(1-x)*A21*exp(-2*self.alpha*A21/(R*T))*(1-exp(-self.alpha*A21/(R*T)))/\
        (R**2*T**2*((1-x)*exp(-self.alpha*A21/(R*T))+x)**3))),-1000,0]

        x = Actual[1]

        Row2 = [1000*((1-x)*x*(exp(-self.alpha*A12/(R*T))/(R*T*(x*exp(-self.alpha*A12/(R*T))-x+1)\
        )-self.alpha*A12*exp(-self.alpha*A12/(R*T))/(R**2*T**2*(x*exp(-self.alpha*A12/\
        (R*T))-x+1))+self.alpha*x*A12*exp(-2*self.alpha*A12/(R*T))/(R**2*T**2*(x*\
        exp(-self.alpha*A12/(R*T))-x+1)**2))),1000*((1-x)*x*(exp(-self.alpha*A21/(R*T))/\
        (R*T*((1-x)*exp(-self.alpha*A21/(R*T))+x))-self.alpha*A21*exp(-self.alpha*A21/\
        (R*T))/(R**2*T**2*((1-x)*exp(-self.alpha*A21/(R*T))+x))+self.alpha*(1-x)*\
        A21*exp(-2*self.alpha*A21/(R*T))/(R**2*T**2*((1-x)*exp(-self.alpha*A21/(R\
        *T))+x)**2))),-1000*x,-1000]
 
        Row4 = [1000*(-x*(exp(-self.alpha*A12/(R*T))/(R*T*(x*exp(-self.alpha*A12/(R*T))-x+1))-self.alpha\
        *A12*exp(-self.alpha*A12/(R*T))/(R**2*T**2*(x*exp(-self.alpha*A12/(R*T)\
        )-x+1))+self.alpha*x*A12*exp(-2*self.alpha*A12/(R*T))/(R**2*T**2*(x*exp(-\
        self.alpha*A12/(R*T))-x+1)**2))+(1-x)*(exp(-self.alpha*A12/(R*T))/(R*T*(x\
        *exp(-self.alpha*A12/(R*T))-x+1))-self.alpha*A12*exp(-self.alpha*A12/(R*T))/(\
        R**2*T**2*(x*exp(-self.alpha*A12/(R*T))-x+1))+self.alpha*x*A12*exp(-2*self.alpha\
        *A12/(R*T))/(R**2*T**2*(x*exp(-self.alpha*A12/(R*T))-x+1)**2))+(1-x\
        )*x*(-(exp(-self.alpha*A12/(R*T))-1)*exp(-self.alpha*A12/(R*T))/(R*T*(x*\
        exp(-self.alpha*A12/(R*T))-x+1)**2)+self.alpha*A12*(exp(-self.alpha*A12/(R*T))-\
        1)*exp(-self.alpha*A12/(R*T))/(R**2*T**2*(x*exp(-self.alpha*A12/(R*T))-x+\
        1)**2)+self.alpha*A12*exp(-2*self.alpha*A12/(R*T))/(R**2*T**2*(x*exp(-self.alpha\
        *A12/(R*T))-x+1)**2)-2*self.alpha*x*A12*exp(-2*self.alpha*A12/(R*T))*(exp\
        (-self.alpha*A12/(R*T))-1)/(R**2*T**2*(x*exp(-self.alpha*A12/(R*T))-x+1\
        )**3))),1000*(-x*(exp(-self.alpha*A21/(R*T))/(R*T*((1-x)*exp(-self.alpha*A21/(R*\
        T))+x))-self.alpha*A21*exp(-self.alpha*A21/(R*T))/(R**2*T**2*((1-x)*exp(-\
        self.alpha*A21/(R*T))+x))+self.alpha*(1-x)*A21*exp(-2*self.alpha*A21/(R*T))/(\
        R**2*T**2*((1-x)*exp(-self.alpha*A21/(R*T))+x)**2))+(1-x)*(exp(-self.alpha\
        *A21/(R*T))/(R*T*((1-x)*exp(-self.alpha*A21/(R*T))+x))-self.alpha*A21*exp\
        (-self.alpha*A21/(R*T))/(R**2*T**2*((1-x)*exp(-self.alpha*A21/(R*T))+x))+\
        self.alpha*(1-x)*A21*exp(-2*self.alpha*A21/(R*T))/(R**2*T**2*((1-x)*exp(-\
        self.alpha*A21/(R*T))+x)**2))+(1-x)*x*(-(1-exp(-self.alpha*A21/(R*T)))*exp\
        (-self.alpha*A21/(R*T))/(R*T*((1-x)*exp(-self.alpha*A21/(R*T))+x)**2)+self.alpha\
        *A21*(1-exp(-self.alpha*A21/(R*T)))*exp(-self.alpha*A21/(R*T))/(R**2*\
        T**2*((1-x)*exp(-self.alpha*A21/(R*T))+x)**2)-self.alpha*A21*exp(-2*self.alpha*\
        A21/(R*T))/(R**2*T**2*((1-x)*exp(-self.alpha*A21/(R*T))+x)**2)-2*self.alpha\
        *(1-x)*A21*exp(-2*self.alpha*A21/(R*T))*(1-exp(-self.alpha*A21/(R*T)))/\
        (R**2*T**2*((1-x)*exp(-self.alpha*A21/(R*T))+x)**3))),-1000,0]
   
        return array([Row1, Row2, Row3, Row4])

class SystemUNIQUAC:
    def __init__(self, M, z):
        self.r1 = M['r1']
        self.r2 = M['r2']
        self.q1 = M['q1']
        self.q2 = M['q2']
        self.z = z

    def SystemEquations(self, Params, Actual, R, T):

        A12 = Params[0]*R
        A21 = Params[1]*R
        m = Params[2]
        b = Params[3]
               
        x = Actual[0]

        Eq1 = -self.q1*x*log(self.q2*(1-x)*exp(-A21/(R*T))/(self.q1*x+self.q2*(1-x))+self.q1*x/(self.q1*x+self.q2*(\
        1-x)))-self.q2*(1-x)*log(self.q1*x*exp(-A12/(R*T))/(self.q1*x+self.q2*(1-x))+self.q2*(1-\
        x)/(self.q1*x+self.q2*(1-x)))+(self.q2*(1-x)*log(self.q2*(self.r1*x+self.r2*(1-x))/(self.r2*(self.q1*x+\
        self.q2*(1-x))))+self.q1*x*log(self.q1*(self.r1*x+self.r2*(1-x))/(self.r1*(self.q1*x+self.q2*(1-x)))))*\
        self.z/2.0E+0+(1-x)*log(self.r2/(self.r1*x+self.r2*(1-x)))+x*log(self.r1/(self.r1*x+self.r2*(1-x))\
        )+x*log(x)-m*x+log(1-x)*(1-x)-b

        Eq3 = -self.q1*log(self.q2*(1-x)*exp(-A21/(R*T))/(self.q1*x+self.q2*(1-x))+self.q1*x/(self.q1*x+self.q2*(1-\
        x)))+self.q2*log(self.q1*x*exp(-A12/(R*T))/(self.q1*x+self.q2*(1-x))+self.q2*(1-x)/(self.q1*x\
        +self.q2*(1-x)))-self.q1*x*(-self.q2*exp(-A21/(R*T))/(self.q1*x+self.q2*(1-x))-(self.q1-self.q2)*self.q2\
        *(1-x)*exp(-A21/(R*T))/(self.q1*x+self.q2*(1-x))**2+self.q1/(self.q1*x+self.q2*(1-x))-self.q1\
        *(self.q1-self.q2)*x/(self.q1*x+self.q2*(1-x))**2)/(self.q2*(1-x)*exp(-A21/(R*T))/(self.q1*x\
        +self.q2*(1-x))+self.q1*x/(self.q1*x+self.q2*(1-x)))-self.q2*(1-x)*(self.q1*exp(-A12/(R*T))/(\
        self.q1*x+self.q2*(1-x))-self.q1*(self.q1-self.q2)*x*exp(-A12/(R*T))/(self.q1*x+self.q2*(1-x))**2-\
        self.q2/(self.q1*x+self.q2*(1-x))-(self.q1-self.q2)*self.q2*(1-x)/(self.q1*x+self.q2*(1-x))**2)/(self.q1*x*exp\
        (-A12/(R*T))/(self.q1*x+self.q2*(1-x))+self.q2*(1-x)/(self.q1*x+self.q2*(1-x)))+(-self.q2*log\
        (self.q2*(self.r1*x+self.r2*(1-x))/(self.r2*(self.q1*x+self.q2*(1-x))))+self.q1*log(self.q1*(self.r1*x+self.r2*\
        (1-x))/(self.r1*(self.q1*x+self.q2*(1-x))))+self.r1*x*(self.q1*x+self.q2*(1-x))*(self.q1*(self.r1-self.r2)/(\
        self.r1*(self.q1*x+self.q2*(1-x)))-self.q1*(self.q1-self.q2)*(self.r1*x+self.r2*(1-x))/(self.r1*(self.q1*x+self.q2*(1-\
        x))**2))/(self.r1*x+self.r2*(1-x))+self.r2*(1-x)*(self.q1*x+self.q2*(1-x))*(self.q2*(self.r1-self.r2)/(\
        self.r2*(self.q1*x+self.q2*(1-x)))-(self.q1-self.q2)*self.q2*(self.r1*x+self.r2*(1-x))/(self.r2*(self.q1*x+self.q2*(1-\
        x))**2))/(self.r1*x+self.r2*(1-x)))*self.z/2.0E+0-log(self.r2/(self.r1*x+self.r2*(1-x)))+log(\
        self.r1/(self.r1*x+self.r2*(1-x)))+log(x)-(self.r1-self.r2)*x/(self.r1*x+self.r2*(1-x))-(self.r1-self.r2)*(1\
        -x)/(self.r1*x+self.r2*(1-x))-log(1-x)-m
        
        x = Actual[1]

        Eq2 = -self.q1*x*log(self.q2*(1-x)*exp(-A21/(R*T))/(self.q1*x+self.q2*(1-x))+self.q1*x/(self.q1*x+self.q2*(\
        1-x)))-self.q2*(1-x)*log(self.q1*x*exp(-A12/(R*T))/(self.q1*x+self.q2*(1-x))+self.q2*(1-\
        x)/(self.q1*x+self.q2*(1-x)))+(self.q2*(1-x)*log(self.q2*(self.r1*x+self.r2*(1-x))/(self.r2*(self.q1*x+\
        self.q2*(1-x))))+self.q1*x*log(self.q1*(self.r1*x+self.r2*(1-x))/(self.r1*(self.q1*x+self.q2*(1-x)))))*\
        self.z/2.0E+0+(1-x)*log(self.r2/(self.r1*x+self.r2*(1-x)))+x*log(self.r1/(self.r1*x+self.r2*(1-x))\
        )+x*log(x)-m*x+log(1-x)*(1-x)-b

        Eq4 = -self.q1*log(self.q2*(1-x)*exp(-A21/(R*T))/(self.q1*x+self.q2*(1-x))+self.q1*x/(self.q1*x+self.q2*(1-\
        x)))+self.q2*log(self.q1*x*exp(-A12/(R*T))/(self.q1*x+self.q2*(1-x))+self.q2*(1-x)/(self.q1*x\
        +self.q2*(1-x)))-self.q1*x*(-self.q2*exp(-A21/(R*T))/(self.q1*x+self.q2*(1-x))-(self.q1-self.q2)*self.q2\
        *(1-x)*exp(-A21/(R*T))/(self.q1*x+self.q2*(1-x))**2+self.q1/(self.q1*x+self.q2*(1-x))-self.q1\
        *(self.q1-self.q2)*x/(self.q1*x+self.q2*(1-x))**2)/(self.q2*(1-x)*exp(-A21/(R*T))/(self.q1*x\
        +self.q2*(1-x))+self.q1*x/(self.q1*x+self.q2*(1-x)))-self.q2*(1-x)*(self.q1*exp(-A12/(R*T))/(\
        self.q1*x+self.q2*(1-x))-self.q1*(self.q1-self.q2)*x*exp(-A12/(R*T))/(self.q1*x+self.q2*(1-x))**2-\
        self.q2/(self.q1*x+self.q2*(1-x))-(self.q1-self.q2)*self.q2*(1-x)/(self.q1*x+self.q2*(1-x))**2)/(self.q1*x*exp\
        (-A12/(R*T))/(self.q1*x+self.q2*(1-x))+self.q2*(1-x)/(self.q1*x+self.q2*(1-x)))+(-self.q2*log\
        (self.q2*(self.r1*x+self.r2*(1-x))/(self.r2*(self.q1*x+self.q2*(1-x))))+self.q1*log(self.q1*(self.r1*x+self.r2*\
        (1-x))/(self.r1*(self.q1*x+self.q2*(1-x))))+self.r1*x*(self.q1*x+self.q2*(1-x))*(self.q1*(self.r1-self.r2)/(\
        self.r1*(self.q1*x+self.q2*(1-x)))-self.q1*(self.q1-self.q2)*(self.r1*x+self.r2*(1-x))/(self.r1*(self.q1*x+self.q2*(1-\
        x))**2))/(self.r1*x+self.r2*(1-x))+self.r2*(1-x)*(self.q1*x+self.q2*(1-x))*(self.q2*(self.r1-self.r2)/(\
        self.r2*(self.q1*x+self.q2*(1-x)))-(self.q1-self.q2)*self.q2*(self.r1*x+self.r2*(1-x))/(self.r2*(self.q1*x+self.q2*(1-\
        x))**2))/(self.r1*x+self.r2*(1-x)))*self.z/2.0E+0-log(self.r2/(self.r1*x+self.r2*(1-x)))+log(\
        self.r1/(self.r1*x+self.r2*(1-x)))+log(x)-(self.r1-self.r2)*x/(self.r1*x+self.r2*(1-x))-(self.r1-self.r2)*(1\
        -x)/(self.r1*x+self.r2*(1-x))-log(1-x)-m

        return array([Eq1, Eq2, Eq3, Eq4])

    def SystemEquationsJac(self, Params, Actual, R, T):

        A12 = Params[0]*R
        A21 = Params[1]*R
        m = Params[2]
        b = Params[3]
               
        x = Actual[0]

        Row1 = [self.q1*self.q2*(1-x)*x*exp(-A12/(R*T))/((self.q1*x+self.q2*(1-x))*R*T*(self.q1*x*exp(-A12\
        /(R*T))/(self.q1*x+self.q2*(1-x))+self.q2*(1-x)/(self.q1*x+self.q2*(1-x)))),self.q1*self.q2*(1-x)*\
        x*exp(-A21/(R*T))/((self.q1*x+self.q2*(1-x))*R*T*(self.q2*(1-x)*exp(-A21/(R*T)\
        )/(self.q1*x+self.q2*(1-x))+self.q1*x/(self.q1*x+self.q2*(1-x)))),-x,-1]

        Row3 = [-self.q1*self.q2*x*exp(-A12/(R*T))/((self.q1*x+self.q2*(1-x))*R*T*(self.q1*x*exp(-A12/(R*T\
        ))/(self.q1*x+self.q2*(1-x))+self.q2*(1-x)/(self.q1*x+self.q2*(1-x))))-self.q1*self.q2*(1-x)*x*exp\
        (-A12/(R*T))*(self.q1*exp(-A12/(R*T))/(self.q1*x+self.q2*(1-x))-self.q1*(self.q1-self.q2)*x*exp\
        (-A12/(R*T))/(self.q1*x+self.q2*(1-x))**2-self.q2/(self.q1*x+self.q2*(1-x))-(self.q1-self.q2)*self.q2\
        *(1-x)/(self.q1*x+self.q2*(1-x))**2)/((self.q1*x+self.q2*(1-x))*R*T*(self.q1*x*exp(-A12/\
        (R*T))/(self.q1*x+self.q2*(1-x))+self.q2*(1-x)/(self.q1*x+self.q2*(1-x)))**2)-self.q2*(1-x)*(\
        self.q1*(self.q1-self.q2)*x*exp(-A12/(R*T))/((self.q1*x+self.q2*(1-x))**2*R*T)-self.q1*exp(-A12\
        /(R*T))/((self.q1*x+self.q2*(1-x))*R*T))/(self.q1*x*exp(-A12/(R*T))/(self.q1*x+self.q2\
        *(1-x))+self.q2*(1-x)/(self.q1*x+self.q2*(1-x))),self.q1*self.q2*(1-x)*exp(-A21/(R*T))/(\
        (self.q1*x+self.q2*(1-x))*R*T*(self.q2*(1-x)*exp(-A21/(R*T))/(self.q1*x+self.q2*(1-x))+self.q1\
        *x/(self.q1*x+self.q2*(1-x))))-self.q1*self.q2*(1-x)*x*exp(-A21/(R*T))*(-self.q2*exp(-A21\
        /(R*T))/(self.q1*x+self.q2*(1-x))-(self.q1-self.q2)*self.q2*(1-x)*exp(-A21/(R*T))/(self.q1*\
        x+self.q2*(1-x))**2+self.q1/(self.q1*x+self.q2*(1-x))-self.q1*(self.q1-self.q2)*x/(self.q1*x+self.q2*(1-x))**2\
        )/((self.q1*x+self.q2*(1-x))*R*T*(self.q2*(1-x)*exp(-A21/(R*T))/(self.q1*x+self.q2*(1-\
        x))+self.q1*x/(self.q1*x+self.q2*(1-x)))**2)-self.q1*x*(self.q2*exp(-A21/(R*T))/((self.q1*x+self.q2\
        *(1-x))*R*T)+(self.q1-self.q2)*self.q2*(1-x)*exp(-A21/(R*T))/((self.q1*x+self.q2*(1-x))\
        **2*R*T))/(self.q2*(1-x)*exp(-A21/(R*T))/(self.q1*x+self.q2*(1-x))+self.q1*x/(self.q1*x+\
        self.q2*(1-x))),-1,0]

        x = Actual[1]

        Row2 = [self.q1*self.q2*(1-x)*x*exp(-A12/(R*T))/((self.q1*x+self.q2*(1-x))*R*T*(self.q1*x*exp(-A12\
        /(R*T))/(self.q1*x+self.q2*(1-x))+self.q2*(1-x)/(self.q1*x+self.q2*(1-x)))),self.q1*self.q2*(1-x)*\
        x*exp(-A21/(R*T))/((self.q1*x+self.q2*(1-x))*R*T*(self.q2*(1-x)*exp(-A21/(R*T)\
        )/(self.q1*x+self.q2*(1-x))+self.q1*x/(self.q1*x+self.q2*(1-x)))),-x,-1]

        Row4 = [-self.q1*self.q2*x*exp(-A12/(R*T))/((self.q1*x+self.q2*(1-x))*R*T*(self.q1*x*exp(-A12/(R*T\
        ))/(self.q1*x+self.q2*(1-x))+self.q2*(1-x)/(self.q1*x+self.q2*(1-x))))-self.q1*self.q2*(1-x)*x*exp\
        (-A12/(R*T))*(self.q1*exp(-A12/(R*T))/(self.q1*x+self.q2*(1-x))-self.q1*(self.q1-self.q2)*x*exp\
        (-A12/(R*T))/(self.q1*x+self.q2*(1-x))**2-self.q2/(self.q1*x+self.q2*(1-x))-(self.q1-self.q2)*self.q2\
        *(1-x)/(self.q1*x+self.q2*(1-x))**2)/((self.q1*x+self.q2*(1-x))*R*T*(self.q1*x*exp(-A12/\
        (R*T))/(self.q1*x+self.q2*(1-x))+self.q2*(1-x)/(self.q1*x+self.q2*(1-x)))**2)-self.q2*(1-x)*(\
        self.q1*(self.q1-self.q2)*x*exp(-A12/(R*T))/((self.q1*x+self.q2*(1-x))**2*R*T)-self.q1*exp(-A12\
        /(R*T))/((self.q1*x+self.q2*(1-x))*R*T))/(self.q1*x*exp(-A12/(R*T))/(self.q1*x+self.q2\
        *(1-x))+self.q2*(1-x)/(self.q1*x+self.q2*(1-x))),self.q1*self.q2*(1-x)*exp(-A21/(R*T))/(\
        (self.q1*x+self.q2*(1-x))*R*T*(self.q2*(1-x)*exp(-A21/(R*T))/(self.q1*x+self.q2*(1-x))+self.q1\
        *x/(self.q1*x+self.q2*(1-x))))-self.q1*self.q2*(1-x)*x*exp(-A21/(R*T))*(-self.q2*exp(-A21\
        /(R*T))/(self.q1*x+self.q2*(1-x))-(self.q1-self.q2)*self.q2*(1-x)*exp(-A21/(R*T))/(self.q1*\
        x+self.q2*(1-x))**2+self.q1/(self.q1*x+self.q2*(1-x))-self.q1*(self.q1-self.q2)*x/(self.q1*x+self.q2*(1-x))**2\
        )/((self.q1*x+self.q2*(1-x))*R*T*(self.q2*(1-x)*exp(-A21/(R*T))/(self.q1*x+self.q2*(1-\
        x))+self.q1*x/(self.q1*x+self.q2*(1-x)))**2)-self.q1*x*(self.q2*exp(-A21/(R*T))/((self.q1*x+self.q2\
        *(1-x))*R*T)+(self.q1-self.q2)*self.q2*(1-x)*exp(-A21/(R*T))/((self.q1*x+self.q2*(1-x))\
        **2*R*T))/(self.q2*(1-x)*exp(-A21/(R*T))/(self.q1*x+self.q2*(1-x))+self.q1*x/(self.q1*x+\
        self.q2*(1-x))),-1,0]

        return array([Row1, Row2, Row3, Row4])
