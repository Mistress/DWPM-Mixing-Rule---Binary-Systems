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

class SquareError(ErrorMeasure):
    def Error(self):
        SquareError = (self.Predicted - self.Actual)**2
        return SquareError

class SumAbs(ErrorMeasure):
    def Error(self):
        AbsoluteError = sum(abs(self.Predicted - self.Actual))
        return AbsoluteError

class AbsError(ErrorMeasure):
    def Error(self):
        AbsError = abs(self.Predicted - self.Actual)

class SumRelative(ErrorMeasure):
    def Error(self):
        RelativeError = sum(abs(self.Predicted - self.Actual)/self.Actual)
        return RelativeError

class RelError(ErrorMeasure):
    def Error(self):
        RelError = abs(self.Predicted - self.Actual)/self.Actual
            



    
