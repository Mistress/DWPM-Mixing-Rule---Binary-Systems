#!/usr/bin/env python

from scipy import *


class SystemDWPM:
    def __init__(self, Cell_s):
        self.s_1 = Cell_s[0]
        self.s_2 = Cell_s[1]

    def SystemEquations(self, Params, Actual):
        
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

    def SystemEquationsJac(self, Params, Actual):
        
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

#class SystemNRTL:
#    def __init__(self, alpha):
#        self.s_1 = Cell_s[0]
#        self.s_2 = Cell_s[1]

#    def SystemEquations(self, Params, Actual):

#    def SystemEquations(self, Params, Actual):
