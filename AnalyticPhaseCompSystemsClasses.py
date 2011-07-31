 #!/usr/bin/env python

from scipy import *


class SystemDWPM:
    def __init__(self, Cell_s):
        self.s_1 = Cell_s[0]
        self.s_2 = Cell_s[1]

    def SystemEquations(self, Predicted, Params, R, T):
        
        lambda_12 = Params[0]
        lambda_21 = Params[1]
        
        x_1 = Predicted[0]
        x_2 = Predicted[1]
        m = Predicted[2]
        b = Predicted[3]


        Eq1 = -(1-x_1)*log(lambda_21**self.s_2*x_1-x_1+1)/self.s_2-x_1*log(x_1+lambda_12**self.s_1*(1-x_1))/self.s_1+x_1*log(x_1)-m*x_1+log(1-x_1)*(1-x_1)-b
        
        Eq2 = log(lambda_21**self.s_2*x_1-x_1+1)/self.s_2-log(x_1+lambda_12**self.s_1*(1-x_1))/self.s_1+log(x_1)-(lambda_21**self.s_2-1)*(1-x_1)/(self.s_2*(lambda_21**self.s_2*x_1\
        -x_1+1))-(1-lambda_12**self.s_1)*x_1/(self.s_1*(x_1+lambda_12**self.s_1*(1-x_1)))-log(1-x_1)-m
        
        Eq3 = -(1-x_2)*log(lambda_21**self.s_2*x_2-x_2+1)/self.s_2-x_2*log(x_2+lambda_12**self.s_1*(1-x_2))/self.s_1+x_2*log(x_2)-m*x_2+log(1-x_2)*(1-x_2)-b
        
        Eq4 = log(lambda_21**self.s_2*x_2-x_2+1)/self.s_2-log(x_2+lambda_12**self.s_1*(1-x_2))/self.s_1+log(x_2)-(lambda_21**self.s_2-1)*(1-x_2)/(self.s_2*(lambda_21**self.s_2*x_2\
        -x_2+1))-(1-lambda_12**self.s_1)*x_2/(self.s_1*(x_2+lambda_12**self.s_1*(1-x_2)))-log(1-x_2)-m
       
        return real(array([Eq1, Eq2, Eq3, Eq4]))

    def SystemEquationsJac(self, Predicted, Params, R, T):
        
        lambda_12 = Params[0]
        lambda_21 = Params[1]
        
        x_1 = Predicted[0]
        x_2 = Predicted[1]
        m = Predicted[2]
        b = Predicted[3]

        Row1 = [log(lambda_21**self.s_2*x_1-x_1+1)/self.s_2-log(x_1+lambda_12**self.s_1*(1-x_1))/self.s_1+log(x_1)-(lambda_21**self.s_2-1)*(1-x_1)/(self.s_2*(lambda_21**self.s_2*\
        x_1-x_1+1))-(1-lambda_12**self.s_1)*x_1/(self.s_1*(x_1+lambda_12**self.s_1*(1-x_1)))-log(1-x_1)-m,0,-x_1,-1]       
        
        Row2 = [2*(lambda_21**self.s_2-1)/(self.s_2*(lambda_21**self.s_2*x_1-x_1+1))+(lambda_21**self.s_2-1)**2*(1-x_1)/(self.s_2*(lambda_21**self.s_2*x_1-x_1+1)**2)-2*(1-lambda_12\
        **self.s_1)/(self.s_1*(x_1+lambda_12**self.s_1*(1-x_1)))+(1-lambda_12**self.s_1)**2*x_1/(self.s_1*(x_1+lambda_12**self.s_1*(1-x_1))**2)+1/x_1+1/(1-x_1),0,-1,0]
        
        Row3 = [0,log(lambda_21**self.s_2*x_2-x_2+1)/self.s_2-log(x_2+lambda_12**self.s_1*(1-x_2))/self.s_1+log(x_2)-(lambda_21**self.s_2-1)*(1-x_2)/(self.s_2*(lambda_21**self.s_2*\
        x_2-x_2+1))-(1-lambda_12**self.s_1)*x_2/(self.s_1*(x_2+lambda_12**self.s_1*(1-x_2)))-log(1-x_2)-m,-x_2,-1]

        Row4 = [0,2*(lambda_21**self.s_2-1)/(self.s_2*(lambda_21**self.s_2*x_2-x_2+1))+(lambda_21**self.s_2-1)**2*(1-x_2)/(self.s_2*(lambda_21**self.s_2*x_2-x_2+1)**2)-2*(1-lambda_12\
        **self.s_1)/(self.s_1*(x_2+lambda_12**self.s_1*(1-x_2)))+(1-lambda_12**self.s_1)**2*x_2/(self.s_1*(x_2+lambda_12**self.s_1*(1-x_2))**2)+1/x_2+1/(1-x_2),-1,0]
        
        return real(array([Row1, Row2, Row3, Row4]))

class SystemNRTL:
    def __init__(self, alpha):
        self.alpha = alpha
    
    def SystemEquations(self, Predicted, Params, R, T):

        A12 = Params[0]*R
        A21 = Params[1]*R
        
        x_1 = Predicted[0]
        x_2 = Predicted[1]
        m = Predicted[2]
        b = Predicted[3]

        Eq1 = (1-x_1)*x_1*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_1)*exp(-self.alpha*A21/(R*T))+x_1))+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x_1*\
        exp(-self.alpha*A12/(R*T))-x_1+1)))+x_1*log(x_1)-m*x_1+log(1-x_1)*(1-x_1)-b
        
        Eq2 = -x_1*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_1)*exp(-self.alpha*A21/(R*T))+x_1))+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x_1*exp(-self.alpha*A12/\
        (R*T))-x_1+1)))+(1-x_1)*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_1)*exp(-self.alpha*A21/(R*T))+x_1))+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x_1*\
        exp(-self.alpha*A12/(R*T))-x_1+1)))+(1-x_1)*x_1*(-A21*(1-exp(-self.alpha*A21/(R*T)))*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_1)*exp(-self.alpha*A21/(\
        R*T))+x_1)**2)-A12*(exp(-self.alpha*A12/(R*T))-1)*exp(-self.alpha*A12/(R*T))/(R*T*(x_1*exp(-self.alpha*A12/(R*T))-x_1+1)**2))+log(x_1)-log(1-x_1)-m

        Eq3 = (1-x_2)*x_2*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_2)*exp(-self.alpha*A21/(R*T))+x_2))+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x_2*\
        exp(-self.alpha*A12/(R*T))-x_2+1)))+x_2*log(x_2)-m*x_2+log(1-x_2)*(1-x_2)-b

        Eq4 = -x_2*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_2)*exp(-self.alpha*A21/(R*T))+x_2))+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x_2*exp(-self.alpha*A12/\
        (R*T))-x_2+1)))+(1-x_2)*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_2)*exp(-self.alpha*A21/(R*T))+x_2))+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x_2*\
        exp(-self.alpha*A12/(R*T))-x_2+1)))+(1-x_2)*x_2*(-A21*(1-exp(-self.alpha*A21/(R*T)))*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_2)*exp(-self.alpha*A21/(\
        R*T))+x_2)**2)-A12*(exp(-self.alpha*A12/(R*T))-1)*exp(-self.alpha*A12/(R*T))/(R*T*(x_2*exp(-self.alpha*A12/(R*T))-x_2+1)**2))+log(x_2)-log(1-x_2)-m
        
        return real(array([Eq1, Eq2, Eq3, Eq4]))

    def SystemEquationsJac(self, Predicted,  Params, R, T):

        A12 = Params[0]*R
        A21 = Params[1]*R

        x_1 = Predicted[0]
        x_2 = Predicted[1]
        m = Predicted[2]
        b = Predicted[3]

        Row1 = [-x_1*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_1)*exp(-self.alpha*A21/(R*T))+x_1))+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x_1*exp(-self.alpha*A12/\
        (R*T))-x_1+1)))+(1-x_1)*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_1)*exp(-self.alpha*A21/(R*T))+x_1))+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x_1*\
        exp(-self.alpha*A12/(R*T))-x_1+1)))+(1-x_1)*x_1*(-A21*(1-exp(-self.alpha*A21/(R*T)))*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_1)*exp(-self.alpha*A21/\
        (R*T))+x_1)**2)-A12*(exp(-self.alpha*A12/(R*T))-1)*exp(-self.alpha*A12/(R*T))/(R*T*(x_1*exp(-self.alpha*A12/(R*T))-x_1+1)**2))+log(x_1)-log(1-x_1)-m,0,-x_1,-1]

        Row2 = [-2*A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_1)*exp(-self.alpha*A21/(R*T))+ x_1))-2*x_1*(-A21*(1-exp(-self.alpha*A21/(R*T)))*exp(-self.alpha*A21/\
        (R*T))/(R*T*((1-x_1)*exp(-self.alpha*A21/(R*T))+x_1)**2)-A12*(exp(-self.alpha*A12/(R*T))-1)*exp(-self.alpha*A12/(R*T))/(R*T*(x_1*exp(-self.alpha*A12/(\
        R*T))-x_1+1)**2))+2*(1-x_1)*(-A21*(1-exp(-self.alpha*A21/(R*T)))*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_1)*exp(-self.alpha*A21/(R*T))+x_1)**2)\
        -A12*(exp(-self.alpha*A12/(R*T))-1)*exp(-self.alpha*A12/(R*T))/(R*T*(x_1*exp(-self.alpha*A12/(R*T))-x_1+1)**2))+(1-x_1)*x_1*(2*A21*(1-exp(-self.alpha*\
        A21/(R*T)))**2*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_1)*exp(-self.alpha*A21/(R*T))+x_1)**3)+2*A12*(exp(-self.alpha*A12/(R*T))-1)**2*exp(\
        -self.alpha*A12/(R*T))/(R*T*(x_1*exp(-self.alpha*A12/(R*T))-x_1+1)**3))-2*A12*exp(-self.alpha*A12/(R*T))/(R*T*(x_1*exp(-self.alpha*A12/(R*T))-x_1+1))+1/\
        x_1+1/(1-x_1),0,-1,0]

        Row3 = [0,-x_2*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_2)*exp(-self.alpha*A21/(R*T))+x_2))+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x_2*exp(-self.alpha*A12/\
        (R*T))-x_2+1)))+(1-x_2)*(A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_2)*exp(-self.alpha*A21/(R*T))+x_2))+A12*exp(-self.alpha*A12/(R*T))/(R*T*(x_2*\
        exp(-self.alpha*A12/(R*T))-x_2+1)))+(1-x_2)*x_2*(-A21*(1-exp(-self.alpha*A21/(R*T)))*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_2)*exp(-self.alpha*A21\
        /(R*T))+x_2)**2)-A12*(exp(-self.alpha*A12/(R*T))-1)*exp(-self.alpha*A12/(R*T))/(R*T*(x_2*exp(-self.alpha*A12/(R*T))-x_2+1)**2))+log(x_2)-log(1-x_2)-m,-x_2,-1]

        Row4 = [0,-2*A21*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_2)*exp(-self.alpha*A21/(R*T))+x_2))-2*x_2*(-A21*(1-exp(-self.alpha*A21/(R*T)))*exp(-self.alpha*A21/\
        (R*T))/(R*T*((1-x_2)*exp(-self.alpha*A21/(R*T))+x_2)**2)-A12*(exp(-self.alpha*A12/(R*T))-1)*exp(-self.alpha*A12/(R*T))/(R*T*(x_2*exp(-self.alpha*A12\
        /(R*T))-x_2+1)**2))+2*(1-x_2)*(-A21*(1-exp(-self.alpha*A21/(R*T)))*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_2)*exp(-self.alpha*A21/(R*T))+x_2)**2\
        )-A12*(exp(-self.alpha*A12/(R*T))-1)*exp(-self.alpha*A12/(R*T))/(R*T*(x_2*exp(-self.alpha*A12/(R*T))-x_2+1)**2))+(1-x_2)*x_2*(2*A21*(1-exp(\
        -self.alpha*A21/(R*T)))**2*exp(-self.alpha*A21/(R*T))/(R*T*((1-x_2)*exp(-self.alpha*A21/(R*T))+x_2)**3)+2*A12*(exp(-self.alpha*A12/(R*T))-1)**2*exp(\
        -self.alpha*A12/(R*T))/(R*T*(x_2*exp(-self.alpha*A12/(R*T))-x_2+1)**3))-2*A12*exp(-self.alpha*A12/(R*T))/(R*T*(x_2*exp(-self.alpha*A12/(R*T))-x_2+1))+1/\
        x_2+1/(1-x_2),-1,0]

        return real(array([Row1, Row2, Row3, Row4]))

class SystemUNIQUAC:
    def __init__(self, M, z):
        self.r1 = M['r1']
        self.r2 = M['r2']
        self.q1 = M['q1']
        self.q2 = M['q2']
        self.z = z

    def SystemEquations(self, Predicted, Params, R, T):

        A12 = Params[0]*R
        A21 = Params[1]*R

        x_1 = Predicted[0]
        x_2 = Predicted[1]
        m = Predicted[2]
        b = Predicted[3]

        Eq1 = -self.q1*x_1*log(self.q2*(1-x_1)*exp(-A21/(R*T))/(self.q1*x_1+self.q2*(1-x_1))+self.q1*x_1/\
        (self.q1*x_1+self.q2*(1-x_1)))-self.q2*(1-x_1)*log(self.q1*x_1*exp(-A12/(R*T))/(self.q1*\
        x_1+self.q2*(1-x_1))+self.q2*(1-x_1)/(self.q1*x_1+self.q2*(1-x_1)))+(self.q2*(1-x_1)*\
        log(self.q2*(self.r1*x_1+self.r2*(1-x_1))/(self.r2*(self.q1*x_1+self.q2*(1-x_1))))+self.q1*x_1*log(self.q1\
        *(self.r1*x_1+self.r2*(1-x_1))/(self.r1*(self.q1*x_1+self.q2*(1-x_1)))))*self.z/2.0E+0+(1-x_1\
        )*log(self.r2/(self.r1*x_1+self.r2*(1-x_1)))+x_1*log(self.r1/(self.r1*x_1+self.r2*(1-x_1)))+x_1\
        *log(x_1)-m*x_1+log(1-x_1)*(1-x_1)-b

        Eq2 = -self.q1*log(self.q2*(1-x_1)*exp(-A21/(R*T))/(self.q1*x_1+self.q2*(1-x_1))+self.q1*x_1/(self.q1*\
        x_1+self.q2*(1-x_1)))+self.q2*log(self.q1*x_1*exp(-A12/(R*T))/(self.q1*x_1+self.q2*(1-x_1\
        ))+self.q2*(1-x_1)/(self.q1*x_1+self.q2*(1-x_1)))-self.q1*x_1*(-self.q2*exp(-A21/(R*T))\
        /(self.q1*x_1+self.q2*(1-x_1))-(self.q1-self.q2)*self.q2*(1-x_1)*exp(-A21/(R*T))/(self.q1*x_1\
        +self.q2*(1-x_1))**2+self.q1/(self.q1*x_1+self.q2*(1-x_1))-self.q1*(self.q1-self.q2)*x_1/(self.q1*x_1+self.q2\
        *(1-x_1))**2)/(self.q2*(1-x_1)*exp(-A21/(R*T))/(self.q1*x_1+self.q2*(1-x_1))+\
        self.q1*x_1/(self.q1*x_1+self.q2*(1-x_1)))-self.q2*(1-x_1)*(self.q1*exp(-A12/(R*T))/(self.q1*\
        x_1+self.q2*(1-x_1))-self.q1*(self.q1-self.q2)*x_1*exp(-A12/(R*T))/(self.q1*x_1+self.q2*(1-x_1\
        ))**2-self.q2/(self.q1*x_1+self.q2*(1-x_1))-(self.q1-self.q2)*self.q2*(1-x_1)/(self.q1*x_1+self.q2*(1-\
        x_1))**2)/(self.q1*x_1*exp(-A12/(R*T))/(self.q1*x_1+self.q2*(1-x_1))+self.q2*(1-x_1\
        )/(self.q1*x_1+self.q2*(1-x_1)))+(-self.q2*log(self.q2*(self.r1*x_1+self.r2*(1-x_1))/(self.r2*(self.q1*\
        x_1+self.q2*(1-x_1))))+self.q1*log(self.q1*(self.r1*x_1+self.r2*(1-x_1))/(self.r1*(self.q1*x_1+self.q2*\
        (1-x_1))))+self.r1*x_1*(self.q1*x_1+self.q2*(1-x_1))*(self.q1*(self.r1-self.r2)/(self.r1*(self.q1*x_1+self.q2\
        *(1-x_1)))-self.q1*(self.q1-self.q2)*(self.r1*x_1+self.r2*(1-x_1))/(self.r1*(self.q1*x_1+self.q2*(1-x_1\
        ))**2))/(self.r1*x_1+self.r2*(1-x_1))+self.r2*(1-x_1)*(self.q1*x_1+self.q2*(1-x_1))*(self.q2\
        *(self.r1-self.r2)/(self.r2*(self.q1*x_1+self.q2*(1-x_1)))-(self.q1-self.q2)*self.q2*(self.r1*x_1+self.r2*(1-x_1)\
        )/(self.r2*(self.q1*x_1+self.q2*(1-x_1))**2))/(self.r1*x_1+self.r2*(1-x_1)))*self.z/2.0E+0-\
        log(self.r2/(self.r1*x_1+self.r2*(1-x_1)))+log(self.r1/(self.r1*x_1+self.r2*(1-x_1)))+log(x_1)-\
        (self.r1-self.r2)*x_1/(self.r1*x_1+self.r2*(1-x_1))-(self.r1-self.r2)*(1-x_1)/(self.r1*x_1+self.r2*(1-x_1\
        ))-log(1-x_1)-m

        Eq3 = -self.q1*x_2*log(self.q2*(1-x_2)*exp(-A21/(R*T))/(self.q1*x_2+self.q2*(1-x_2))+self.q1*x_2/\
        (self.q1*x_2+self.q2*(1-x_2)))-self.q2*(1-x_2)*log(self.q1*x_2*exp(-A12/(R*T))/(self.q1*\
        x_2+self.q2*(1-x_2))+self.q2*(1-x_2)/(self.q1*x_2+self.q2*(1-x_2)))+(self.q2*(1-x_2)*\
        log(self.q2*(self.r1*x_2+self.r2*(1-x_2))/(self.r2*(self.q1*x_2+self.q2*(1-x_2))))+self.q1*x_2*log(self.q1\
        *(self.r1*x_2+self.r2*(1-x_2))/(self.r1*(self.q1*x_2+self.q2*(1-x_2)))))*self.z/2.0E+0+(1-x_2\
        )*log(self.r2/(self.r1*x_2+self.r2*(1-x_2)))+x_2*log(self.r1/(self.r1*x_2+self.r2*(1-x_2)))+x_2\
        *log(x_2)-m*x_2+log(1-x_2)*(1-x_2)-b

        Eq4 = -self.q1*log(self.q2*(1-x_2)*exp(-A21/(R*T))/(self.q1*x_2+self.q2*(1-x_2))+self.q1*x_2/(self.q1*\
        x_2+self.q2*(1-x_2)))+self.q2*log(self.q1*x_2*exp(-A12/(R*T))/(self.q1*x_2+self.q2*(1-x_2\
        ))+self.q2*(1-x_2)/(self.q1*x_2+self.q2*(1-x_2)))-self.q1*x_2*(-self.q2*exp(-A21/(R*T))\
        /(self.q1*x_2+self.q2*(1-x_2))-(self.q1-self.q2)*self.q2*(1-x_2)*exp(-A21/(R*T))/(self.q1*x_2\
        +self.q2*(1-x_2))**2+self.q1/(self.q1*x_2+self.q2*(1-x_2))-self.q1*(self.q1-self.q2)*x_2/(self.q1*x_2+self.q2\
        *(1-x_2))**2)/(self.q2*(1-x_2)*exp(-A21/(R*T))/(self.q1*x_2+self.q2*(1-x_2))+\
        self.q1*x_2/(self.q1*x_2+self.q2*(1-x_2)))-self.q2*(1-x_2)*(self.q1*exp(-A12/(R*T))/(self.q1*\
        x_2+self.q2*(1-x_2))-self.q1*(self.q1-self.q2)*x_2*exp(-A12/(R*T))/(self.q1*x_2+self.q2*(1-x_2\
        ))**2-self.q2/(self.q1*x_2+self.q2*(1-x_2))-(self.q1-self.q2)*self.q2*(1-x_2)/(self.q1*x_2+self.q2*(1-\
        x_2))**2)/(self.q1*x_2*exp(-A12/(R*T))/(self.q1*x_2+self.q2*(1-x_2))+self.q2*(1-x_2\
        )/(self.q1*x_2+self.q2*(1-x_2)))+(-self.q2*log(self.q2*(self.r1*x_2+self.r2*(1-x_2))/(self.r2*(self.q1*\
        x_2+self.q2*(1-x_2))))+self.q1*log(self.q1*(self.r1*x_2+self.r2*(1-x_2))/(self.r1*(self.q1*x_2+self.q2*\
        (1-x_2))))+self.r1*x_2*(self.q1*x_2+self.q2*(1-x_2))*(self.q1*(self.r1-self.r2)/(self.r1*(self.q1*x_2+self.q2\
        *(1-x_2)))-self.q1*(self.q1-self.q2)*(self.r1*x_2+self.r2*(1-x_2))/(self.r1*(self.q1*x_2+self.q2*(1-x_2\
        ))**2))/(self.r1*x_2+self.r2*(1-x_2))+self.r2*(1-x_2)*(self.q1*x_2+self.q2*(1-x_2))*(self.q2\
        *(self.r1-self.r2)/(self.r2*(self.q1*x_2+self.q2*(1-x_2)))-(self.q1-self.q2)*self.q2*(self.r1*x_2+self.r2*(1-x_2)\
        )/(self.r2*(self.q1*x_2+self.q2*(1-x_2))**2))/(self.r1*x_2+self.r2*(1-x_2)))*self.z/2.0E+0-\
        log(self.r2/(self.r1*x_2+self.r2*(1-x_2)))+log(self.r1/(self.r1*x_2+self.r2*(1-x_2)))+log(x_2)-\
        (self.r1-self.r2)*x_2/(self.r1*x_2+self.r2*(1-x_2))-(self.r1-self.r2)*(1-x_2)/(self.r1*x_2+self.r2*(1-x_2\
        ))-log(1-x_2)-m
               
        return real(array([Eq1, Eq2, Eq3, Eq4]))

    def SystemEquationsJac(self, Predicted, Params, R, T):

        A12 = Params[0]*R
        A21 = Params[1]*R

        x_1 = Predicted[0]
        x_2 = Predicted[1]
        m = Predicted[2]
        b = Predicted[3]

        Row1 = [-self.q1*log(self.q2*(1-x_1)*exp(-A21/(R*T))/(self.q1*x_1+self.q2*(1-x_1))+self.q1*x_1/(self.q1\
        *x_1+self.q2*(1-x_1)))+self.q2*log(self.q1*x_1*exp(-A12/(R*T))/(self.q1*x_1+self.q2*(1-x_1\
        ))+self.q2*(1-x_1)/(self.q1*x_1+self.q2*(1-x_1)))-self.q1*x_1*(-self.q2*exp(-A21/(R*T)\
        )/(self.q1*x_1+self.q2*(1-x_1))-(self.q1-self.q2)*self.q2*(1-x_1)*exp(-A21/(R*T))/(self.q1*x_1\
        +self.q2*(1-x_1))**2+self.q1/(self.q1*x_1+self.q2*(1-x_1))-self.q1*(self.q1-self.q2)*x_1/(self.q1*x_1+\
        self.q2*(1-x_1))**2)/(self.q2*(1-x_1)*exp(-A21/(R*T))/(self.q1*x_1+self.q2*(1-x_1))\
        +self.q1*x_1/(self.q1*x_1+self.q2*(1-x_1)))-self.q2*(1-x_1)*(self.q1*exp(-A12/(R*T))/(self.q1\
        *x_1+self.q2*(1-x_1))-self.q1*(self.q1-self.q2)*x_1*exp(-A12/(R*T))/(self.q1*x_1+self.q2*(1-x_1\
        ))**2-self.q2/(self.q1*x_1+self.q2*(1-x_1))-(self.q1-self.q2)*self.q2*(1-x_1)/(self.q1*x_1+self.q2*(1\
        -x_1))**2)/(self.q1*x_1*exp(-A12/(R*T))/(self.q1*x_1+self.q2*(1-x_1))+self.q2*(1-x_1\
        )/(self.q1*x_1+self.q2*(1-x_1)))+(-self.q2*log(self.q2*(self.r1*x_1+self.r2*(1-x_1))/(self.r2*(self.q1\
        *x_1+self.q2*(1-x_1))))+self.q1*log(self.q1*(self.r1*x_1+self.r2*(1-x_1))/(self.r1*(self.q1*x_1+self.q2\
        *(1-x_1))))+self.r1*x_1*(self.q1*x_1+self.q2*(1-x_1))*(self.q1*(self.r1-self.r2)/(self.r1*(self.q1*x_1+\
        self.q2*(1-x_1)))-self.q1*(self.q1-self.q2)*(self.r1*x_1+self.r2*(1-x_1))/(self.r1*(self.q1*x_1+self.q2*(1-x_1\
        ))**2))/(self.r1*x_1+self.r2*(1-x_1))+self.r2*(1-x_1)*(self.q1*x_1+self.q2*(1-x_1))*(self.q2\
        *(self.r1-self.r2)/(self.r2*(self.q1*x_1+self.q2*(1-x_1)))-(self.q1-self.q2)*self.q2*(self.r1*x_1+self.r2*(1-x_1\
        ))/(self.r2*(self.q1*x_1+self.q2*(1-x_1))**2))/(self.r1*x_1+self.r2*(1-x_1)))*self.z/2.0E+0-\
        log(self.r2/(self.r1*x_1+self.r2*(1-x_1)))+log(self.r1/(self.r1*x_1+self.r2*(1-x_1)))+log(x_1)\
        -(self.r1-self.r2)*x_1/(self.r1*x_1+self.r2*(1-x_1))-(self.r1-self.r2)*(1-x_1)/(self.r1*x_1+self.r2*(1-\
        x_1))-log(1-x_1)-m,0,-x_1,-1]

        Row2 = [-2*self.q1*(-self.q2*exp(-A21/(R*T))/(self.q1*x_1+self.q2*(1-x_1))-(self.q1-self.q2)*self.q2*(1-x_1)\
        *exp(-A21/(R*T))/(self.q1*x_1+self.q2*(1-x_1))**2+self.q1/(self.q1*x_1+self.q2*(1-x_1))-\
        self.q1*(self.q1-self.q2)*x_1/(self.q1*x_1+self.q2*(1-x_1))**2)/(self.q2*(1-x_1)*exp(-A21/(R*\
        T))/(self.q1*x_1+self.q2*(1-x_1))+self.q1*x_1/(self.q1*x_1+self.q2*(1-x_1)))-self.q1*x_1*(2*(\
        self.q1-self.q2)*self.q2*exp(-A21/(R*T))/(self.q1*x_1+self.q2*(1-x_1))**2+2*(self.q1-self.q2)**2*self.q2\
        *(1-x_1)*exp(-A21/(R*T))/(self.q1*x_1+self.q2*(1-x_1))**3-2*self.q1*(self.q1-self.q2)/(\
        self.q1*x_1+self.q2*(1-x_1))**2+2*self.q1*(self.q1-self.q2)**2*x_1/(self.q1*x_1+self.q2*(1-x_1))**3\
        )/(self.q2*(1-x_1)*exp(-A21/(R*T))/(self.q1*x_1+self.q2*(1-x_1))+self.q1*x_1/(self.q1*x_1\
        +self.q2*(1-x_1)))+self.q1*x_1*(-self.q2*exp(-A21/(R*T))/(self.q1*x_1+self.q2*(1-x_1))\
        -(self.q1-self.q2)*self.q2*(1-x_1)*exp(-A21/(R*T))/(self.q1*x_1+self.q2*(1-x_1))**2+self.q1/(\
        self.q1*x_1+self.q2*(1-x_1))-self.q1*(self.q1-self.q2)*x_1/(self.q1*x_1+self.q2*(1-x_1))**2)**2/(self.q2\
        *(1-x_1)*exp(-A21/(R*T))/(self.q1*x_1+self.q2*(1-x_1))+self.q1*x_1/(self.q1*x_1+self.q2\
        *(1-x_1)))**2+2*self.q2*(self.q1*exp(-A12/(R*T))/(self.q1*x_1+self.q2*(1-x_1))-self.q1*(\
        self.q1-self.q2)*x_1*exp(-A12/(R*T))/(self.q1*x_1+self.q2*(1-x_1))**2-self.q2/(self.q1*x_1+self.q2\
        *(1-x_1))-(self.q1-self.q2)*self.q2*(1-x_1)/(self.q1*x_1+self.q2*(1-x_1))**2)/(self.q1*x_1*\
        exp(-A12/(R*T))/(self.q1*x_1+self.q2*(1-x_1))+self.q2*(1-x_1)/(self.q1*x_1+self.q2*(1-x_1)\
        ))-self.q2*(1-x_1)*(-2*self.q1*(self.q1-self.q2)*exp(-A12/(R*T))/(self.q1*x_1+self.q2*(1-x_1)\
        )**2+2*self.q1*(self.q1-self.q2)**2*x_1*exp(-A12/(R*T))/(self.q1*x_1+self.q2*(1-x_1))**3\
        +2*(self.q1-self.q2)*self.q2/(self.q1*x_1+self.q2*(1-x_1))**2+2*(self.q1-self.q2)**2*self.q2*(1-x_1)/(self.q1\
        *x_1+self.q2*(1-x_1))**3)/(self.q1*x_1*exp(-A12/(R*T))/(self.q1*x_1+self.q2*(1-x_1\
        ))+self.q2*(1-x_1)/(self.q1*x_1+self.q2*(1-x_1)))+self.q2*(1-x_1)*(self.q1*exp(-A12/(R*T\
        ))/(self.q1*x_1+self.q2*(1-x_1))-self.q1*(self.q1-self.q2)*x_1*exp(-A12/(R*T))/(self.q1*x_1+self.q2\
        *(1-x_1))**2-self.q2/(self.q1*x_1+self.q2*(1-x_1))-(self.q1-self.q2)*self.q2*(1-x_1)/(self.q1*x_1\
        +self.q2*(1-x_1))**2)**2/(self.q1*x_1*exp(-A12/(R*T))/(self.q1*x_1+self.q2*(1-x_1))\
        +self.q2*(1-x_1)/(self.q1*x_1+self.q2*(1-x_1)))**2+(2*self.r1*(self.q1*x_1+self.q2*(1-x_1))*(\
        self.q1*(self.r1-self.r2)/(self.r1*(self.q1*x_1+self.q2*(1-x_1)))-self.q1*(self.q1-self.q2)*(self.r1*x_1+self.r2*(1-x_1\
        ))/(self.r1*(self.q1*x_1+self.q2*(1-x_1))**2))/(self.r1*x_1+self.r2*(1-x_1))+(self.q1-self.q2)*self.r1\
        *x_1*(self.q1*(self.r1-self.r2)/(self.r1*(self.q1*x_1+self.q2*(1-x_1)))-self.q1*(self.q1-self.q2)*(self.r1*x_1+self.r2\
        *(1-x_1))/(self.r1*(self.q1*x_1+self.q2*(1-x_1))**2))/(self.r1*x_1+self.r2*(1-x_1))-self.r1*(\
        self.r1-self.r2)*x_1*(self.q1*x_1+self.q2*(1-x_1))*(self.q1*(self.r1-self.r2)/(self.r1*(self.q1*x_1+self.q2*(1-x_1\
        )))-self.q1*(self.q1-self.q2)*(self.r1*x_1+self.r2*(1-x_1))/(self.r1*(self.q1*x_1+self.q2*(1-x_1))**2)\
        )/(self.r1*x_1+self.r2*(1-x_1))**2-2*self.r2*(self.q1*x_1+self.q2*(1-x_1))*(self.q2*(self.r1-self.r2)/(\
        self.r2*(self.q1*x_1+self.q2*(1-x_1)))-(self.q1-self.q2)*self.q2*(self.r1*x_1+self.r2*(1-x_1))/(self.r2*(self.q1*\
        x_1+self.q2*(1-x_1))**2))/(self.r1*x_1+self.r2*(1-x_1))+(self.q1-self.q2)*self.r2*(1-x_1)*(self.q2\
        *(self.r1-self.r2)/(self.r2*(self.q1*x_1+self.q2*(1-x_1)))-(self.q1-self.q2)*self.q2*(self.r1*x_1+self.r2*(1-x_1)\
        )/(self.r2*(self.q1*x_1+self.q2*(1-x_1))**2))/(self.r1*x_1+self.r2*(1-x_1))-(self.r1-self.r2)*self.r2*(\
        1-x_1)*(self.q1*x_1+self.q2*(1-x_1))*(self.q2*(self.r1-self.r2)/(self.r2*(self.q1*x_1+self.q2*(1-x_1)))\
        -(self.q1-self.q2)*self.q2*(self.r1*x_1+self.r2*(1-x_1))/(self.r2*(self.q1*x_1+self.q2*(1-x_1))**2))/(self.r1\
        *x_1+self.r2*(1-x_1))**2+self.r1*x_1*(self.q1*x_1+self.q2*(1-x_1))*(2*self.q1*(self.q1-self.q2)**2\
        *(self.r1*x_1+self.r2*(1-x_1))/(self.r1*(self.q1*x_1+self.q2*(1-x_1))**3)-2*self.q1*(self.q1-self.q2)*\
        (self.r1-self.r2)/(self.r1*(self.q1*x_1+self.q2*(1-x_1))**2))/(self.r1*x_1+self.r2*(1-x_1))+self.r2*(1-\
        x_1)*(self.q1*x_1+self.q2*(1-x_1))*(2*(self.q1-self.q2)**2*self.q2*(self.r1*x_1+self.r2*(1-x_1))/(\
        self.r2*(self.q1*x_1+self.q2*(1-x_1))**3)-2*(self.q1-self.q2)*self.q2*(self.r1-self.r2)/(self.r2*(self.q1*x_1+self.q2*\
        (1-x_1))**2))/(self.r1*x_1+self.r2*(1-x_1)))*self.z/2.0E+0+(self.r1-self.r2)**2*x_1/(self.r1*\
        x_1+self.r2*(1-x_1))**2+(self.r1-self.r2)**2*(1-x_1)/(self.r1*x_1+self.r2*(1-x_1))**2+1/\
        x_1+1/(1-x_1),0,-1,0]

        Row3 = [0,-self.q1*log(self.q2*(1-x_2)*exp(-A21/(R*T))/(self.q1*x_2+self.q2*(1-x_2))+self.q1*x_2/(\
        self.q1*x_2+self.q2*(1-x_2)))+self.q2*log(self.q1*x_2*exp(-A12/(R*T))/(self.q1*x_2+self.q2*(1\
        -x_2))+self.q2*(1-x_2)/(self.q1*x_2+self.q2*(1-x_2)))-self.q1*x_2*(-self.q2*exp(-A21/(R*\
        T))/(self.q1*x_2+self.q2*(1-x_2))-(self.q1-self.q2)*self.q2*(1-x_2)*exp(-A21/(R*T))/(self.q1*\
        x_2+self.q2*(1-x_2))**2+self.q1/(self.q1*x_2+self.q2*(1-x_2))-self.q1*(self.q1-self.q2)*x_2/(self.q1*x_2\
        +self.q2*(1-x_2))**2)/(self.q2*(1-x_2)*exp(-A21/(R*T))/(self.q1*x_2+self.q2*(1-x_2\
        ))+self.q1*x_2/(self.q1*x_2+self.q2*(1-x_2)))-self.q2*(1-x_2)*(self.q1*exp(-A12/(R*T))/(\
        self.q1*x_2+self.q2*(1-x_2))-self.q1*(self.q1-self.q2)*x_2*exp(-A12/(R*T))/(self.q1*x_2+self.q2*(1\
        -x_2))**2-self.q2/(self.q1*x_2+self.q2*(1-x_2))-(self.q1-self.q2)*self.q2*(1-x_2)/(self.q1*x_2+self.q2*\
        (1-x_2))**2)/(self.q1*x_2*exp(-A12/(R*T))/(self.q1*x_2+self.q2*(1-x_2))+self.q2*(1-\
        x_2)/(self.q1*x_2+self.q2*(1-x_2)))+(-self.q2*log(self.q2*(self.r1*x_2+self.r2*(1-x_2))/(self.r2*(\
        self.q1*x_2+self.q2*(1-x_2))))+self.q1*log(self.q1*(self.r1*x_2+self.r2*(1-x_2))/(self.r1*(self.q1*x_2+\
        self.q2*(1-x_2))))+self.r1*x_2*(self.q1*x_2+self.q2*(1-x_2))*(self.q1*(self.r1-self.r2)/(self.r1*(self.q1*x_2\
        +self.q2*(1-x_2)))-self.q1*(self.q1-self.q2)*(self.r1*x_2+self.r2*(1-x_2))/(self.r1*(self.q1*x_2+self.q2*(1\
        -x_2))**2))/(self.r1*x_2+self.r2*(1-x_2))+self.r2*(1-x_2)*(self.q1*x_2+self.q2*(1-x_2))*\
        (self.q2*(self.r1-self.r2)/(self.r2*(self.q1*x_2+self.q2*(1-x_2)))-(self.q1-self.q2)*self.q2*(self.r1*x_2+self.r2*(1-x_2\
        ))/(self.r2*(self.q1*x_2+self.q2*(1-x_2))**2))/(self.r1*x_2+self.r2*(1-x_2)))*self.z/2.0E+0\
        -log(self.r2/(self.r1*x_2+self.r2*(1-x_2)))+log(self.r1/(self.r1*x_2+self.r2*(1-x_2)))+log(x_2\
        )-(self.r1-self.r2)*x_2/(self.r1*x_2+self.r2*(1-x_2))-(self.r1-self.r2)*(1-x_2)/(self.r1*x_2+self.r2*(\
        1-x_2))-log(1-x_2)-m,-x_2,-1]

        Row4 = [0,-2*self.q1*(-self.q2*exp(-A21/(R*T))/(self.q1*x_2+self.q2*(1-x_2))-(self.q1-self.q2)*self.q2*(1-x_2\
        )*exp(-A21/(R*T))/(self.q1*x_2+self.q2*(1-x_2))**2+self.q1/(self.q1*x_2+self.q2*(1-x_2)\
        )-self.q1*(self.q1-self.q2)*x_2/(self.q1*x_2+self.q2*(1-x_2))**2)/(self.q2*(1-x_2)*exp(-A21/(\
        R*T))/(self.q1*x_2+self.q2*(1-x_2))+self.q1*x_2/(self.q1*x_2+self.q2*(1-x_2)))-self.q1*x_2*(2\
        *(self.q1-self.q2)*self.q2*exp(-A21/(R*T))/(self.q1*x_2+self.q2*(1-x_2))**2+2*(self.q1-self.q2)**2\
        *self.q2*(1-x_2)*exp(-A21/(R*T))/(self.q1*x_2+self.q2*(1-x_2))**3-2*self.q1*(self.q1-self.q2)\
        /(self.q1*x_2+self.q2*(1-x_2))**2+2*self.q1*(self.q1-self.q2)**2*x_2/(self.q1*x_2+self.q2*(1-x_2))\
        **3)/(self.q2*(1-x_2)*exp(-A21/(R*T))/(self.q1*x_2+self.q2*(1-x_2))+self.q1*x_2/(self.q1\
        *x_2+self.q2*(1-x_2)))+self.q1*x_2*(-self.q2*exp(-A21/(R*T))/(self.q1*x_2+self.q2*(1-x_2\
        ))-(self.q1-self.q2)*self.q2*(1-x_2)*exp(-A21/(R*T))/(self.q1*x_2+self.q2*(1-x_2))**2+self.q1\
        /(self.q1*x_2+self.q2*(1-x_2))-self.q1*(self.q1-self.q2)*x_2/(self.q1*x_2+self.q2*(1-x_2))**2)**2/\
        (self.q2*(1-x_2)*exp(-A21/(R*T))/(self.q1*x_2+self.q2*(1-x_2))+self.q1*x_2/(self.q1*x_2+\
        self.q2*(1-x_2)))**2+2*self.q2*(self.q1*exp(-A12/(R*T))/(self.q1*x_2+self.q2*(1-x_2))-self.q1\
        *(self.q1-self.q2)*x_2*exp(-A12/(R*T))/(self.q1*x_2+self.q2*(1-x_2))**2-self.q2/(self.q1*x_2+\
        self.q2*(1-x_2))-(self.q1-self.q2)*self.q2*(1-x_2)/(self.q1*x_2+self.q2*(1-x_2))**2)/(self.q1*x_2*\
        exp(-A12/(R*T))/(self.q1*x_2+self.q2*(1-x_2))+self.q2*(1-x_2)/(self.q1*x_2+self.q2*(1-x_2\
        )))-self.q2*(1-x_2)*(-2*self.q1*(self.q1-self.q2)*exp(-A12/(R*T))/(self.q1*x_2+self.q2*(1-x_2\
        ))**2+2*self.q1*(self.q1-self.q2)**2*x_2*exp(-A12/(R*T))/(self.q1*x_2+self.q2*(1-x_2))**3\
        +2*(self.q1-self.q2)*self.q2/(self.q1*x_2+self.q2*(1-x_2))**2+2*(self.q1-self.q2)**2*self.q2*(1-x_2)/\
        (self.q1*x_2+self.q2*(1-x_2))**3)/(self.q1*x_2*exp(-A12/(R*T))/(self.q1*x_2+self.q2*(1-x_2\
        ))+self.q2*(1-x_2)/(self.q1*x_2+self.q2*(1-x_2)))+self.q2*(1-x_2)*(self.q1*exp(-A12/(R\
        *T))/(self.q1*x_2+self.q2*(1-x_2))-self.q1*(self.q1-self.q2)*x_2*exp(-A12/(R*T))/(self.q1*x_2\
        +self.q2*(1-x_2))**2-self.q2/(self.q1*x_2+self.q2*(1-x_2))-(self.q1-self.q2)*self.q2*(1-x_2)/(self.q1*x_2\
        +self.q2*(1-x_2))**2)**2/(self.q1*x_2*exp(-A12/(R*T))/(self.q1*x_2+self.q2*(1-x_2\
        ))+self.q2*(1-x_2)/(self.q1*x_2+self.q2*(1-x_2)))**2+(2*self.r1*(self.q1*x_2+self.q2*(1-x_2))\
        *(self.q1*(self.r1-self.r2)/(self.r1*(self.q1*x_2+self.q2*(1-x_2)))-self.q1*(self.q1-self.q2)*(self.r1*x_2+self.r2*(1-\
        x_2))/(self.r1*(self.q1*x_2+self.q2*(1-x_2))**2))/(self.r1*x_2+self.r2*(1-x_2))+(self.q1-self.q2)*\
        self.r1*x_2*(self.q1*(self.r1-self.r2)/(self.r1*(self.q1*x_2+self.q2*(1-x_2)))-self.q1*(self.q1-self.q2)*(self.r1*x_2+\
        self.r2*(1-x_2))/(self.r1*(self.q1*x_2+self.q2*(1-x_2))**2))/(self.r1*x_2+self.r2*(1-x_2))-self.r1\
        *(self.r1-self.r2)*x_2*(self.q1*x_2+self.q2*(1-x_2))*(self.q1*(self.r1-self.r2)/(self.r1*(self.q1*x_2+self.q2*(1-\
        x_2)))-self.q1*(self.q1-self.q2)*(self.r1*x_2+self.r2*(1-x_2))/(self.r1*(self.q1*x_2+self.q2*(1-x_2))**2\
        ))/(self.r1*x_2+self.r2*(1-x_2))**2-2*self.r2*(self.q1*x_2+self.q2*(1-x_2))*(self.q2*(self.r1-self.r2)\
        /(self.r2*(self.q1*x_2+self.q2*(1-x_2)))-(self.q1-self.q2)*self.q2*(self.r1*x_2+self.r2*(1-x_2))/(self.r2*(self.q1\
        *x_2+self.q2*(1-x_2))**2))/(self.r1*x_2+self.r2*(1-x_2))+(self.q1-self.q2)*self.r2*(1-x_2)*(\
        self.q2*(self.r1-self.r2)/(self.r2*(self.q1*x_2+self.q2*(1-x_2)))-(self.q1-self.q2)*self.q2*(self.r1*x_2+self.r2*(1-x_2\
        ))/(self.r2*(self.q1*x_2+self.q2*(1-x_2))**2))/(self.r1*x_2+self.r2*(1-x_2))-(self.r1-self.r2)*self.r2\
        *(1-x_2)*(self.q1*x_2+self.q2*(1-x_2))*(self.q2*(self.r1-self.r2)/(self.r2*(self.q1*x_2+self.q2*(1-x_2)\
        ))-(self.q1-self.q2)*self.q2*(self.r1*x_2+self.r2*(1-x_2))/(self.r2*(self.q1*x_2+self.q2*(1-x_2))**2))/\
        (self.r1*x_2+self.r2*(1-x_2))**2+self.r1*x_2*(self.q1*x_2+self.q2*(1-x_2))*(2*self.q1*(self.q1-self.q2)\
        **2*(self.r1*x_2+self.r2*(1-x_2))/(self.r1*(self.q1*x_2+self.q2*(1-x_2))**3)-2*self.q1*(self.q1-self.q2\
        )*(self.r1-self.r2)/(self.r1*(self.q1*x_2+self.q2*(1-x_2))**2))/(self.r1*x_2+self.r2*(1-x_2))+self.r2*(\
        1-x_2)*(self.q1*x_2+self.q2*(1-x_2))*(2*(self.q1-self.q2)**2*self.q2*(self.r1*x_2+self.r2*(1-x_2))\
        /(self.r2*(self.q1*x_2+self.q2*(1-x_2))**3)-2*(self.q1-self.q2)*self.q2*(self.r1-self.r2)/(self.r2*(self.q1*x_2+self.q2\
        *(1-x_2))**2))/(self.r1*x_2+self.r2*(1-x_2)))*self.z/2.0E+0+(self.r1-self.r2)**2*x_2/(self.r1\
        *x_2+self.r2*(1-x_2))**2+(self.r1-self.r2)**2*(1-x_2)/(self.r1*x_2+self.r2*(1-x_2))**2+\
        1/x_2+1/(1-x_2),-1,0]
             
        return real(array([Row1, Row2, Row3, Row4]))
