#!/usr/bin/env python

from scipy import *

class ErrorMeasure:
    def __init__ (self, Predicted, Actual):       
        self.Predicted = Predicted
        self.Actual = Actual

        
class SumSquare(ErrorMeasure):
    def Error(self):
        SqrError = sum((self.Predicted - self.Actual)**2)
        return SqrError

class SumAbs(ErrorMeasure):
    def Error(self):
        AbsoluteError = sum(abs(self.Predicted - self.Actual))
        return AbsoluteError

class SumRelative(ErrorMeasure):
    def Error(self):
        RelativeError = sum((self.Predicted - self.Actual)/self.Actual)
        return RelativeError
            



    
