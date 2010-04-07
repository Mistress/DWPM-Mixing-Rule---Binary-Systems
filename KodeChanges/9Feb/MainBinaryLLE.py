#!/usr/bin/env python

from scipy import *
import scipy.optimize
import pickle
import glob
import vanDerWaals
import ErrorClasses
import GibbsClasses
import PhaseStability


class Mixture:
        
    def __init__(self, Compounds, MixtureDataDir, PureDataDir):

        Data ={}
        MixtureName = Compounds[0]
        self.AdachiLuParam ={}
        self.Compounds = Compounds
        for Compound in Compounds:
            PureCompFile = glob.glob(PureDataDir+'/'+Compound+'.dat')
            Data[Compound] = pickle.load(file(PureCompFile[0]))
            if not(Compound == Compounds[0] or Compound == Compounds[-1]):
                MixtureName = MixtureName+'-'+Compound
            if Compound == Compounds[-1]:
                MixtureName = MixtureName+'-'+Compound
        self.Name = MixtureName
        self.Data = Data
        MixtureFileName = MixtureDataDir+'/'+MixtureName + '.dat'

        MixtureFile = glob.glob(MixtureFileName)
        self.M = pickle.load(file(MixtureFile[0]))
        [Slope, VPData] = vanDerWaals.MainSlope(Data, Compounds)
        for Compound in Compounds:
            self.AdachiLuParam[Compound] = Slope[Compound]['BestFit']
    
    def OptFunction(self, params, Model, Actual, T, c):
        
        ModelInstance = getattr(GibbsClasses, Model)(params)        
        Predicted = PhaseStability.CalcPhaseStability(ModelInstance, T, c, self.M)
        Error = ErrorClasses.SumSquare(Predicted ,Actual).Error()

        return Error
    def NonEqConstr(self, params, Model, Actual, T, c):
        
        ModelInstance = getattr(GibbsClasses, Model)(params) 
        deltaGibbsMixTest = array([ModelInstance.deltaGmix(x, T, c, self.M) for x in arange(0.001, 1, 0.001)])
        TangentTestComps = PhaseStability.CalcPhaseStability(ModelInstance, T, c, self.M)
        TangentTest = array([((ModelInstance.FirstDerivative(TangentTestComps[0], T, c, self.M))*(x - TangentTestComps[0])+ ModelInstance.deltaGmix(TangentTestComps[0], T, c, self.M))-ModelInstance.deltaGmix(x, T, c, self.M) for x in linspace(TangentTestComps[0], TangentTestComps[1],1000, True, False)])
        DistinctTest = -1*array([abs(TangentTestComps[0]-TangentTestComps[1])+0.01])
       
        return -1*append(append(deltaGibbsMixTest, TangentTest), DistinctTest)
        
    def BestFitParams(self, Model, InitParams, Bounds):

       # for T in self.M['T']:
        R = 8.314
       # T = self.M['T']
        T = 303.0
        Actual =  array([interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][0])),interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][0]))])
        CompC = vanDerWaals.CompC(self.Compounds, self.Data, self.AdachiLuParam, T, R)
        c = [CompC[Compound] for Compound in Compounds]
            
        params  = scipy.optimize.fmin_slsqp(self.OptFunction, InitParams,[], None, [], self.NonEqConstr, Bounds, None, None, None, (Model, Actual, T, c), 500, 10e-8, 1, 0, 10e-8)
        ModelInstance = getattr(GibbsClasses, Model)(params) 
        print PhaseStability.CalcPhaseStability(ModelInstance, T, c, self.M)
        
        return params
            
            
        


##=============================================================##
#Models = ('DWPM', 'NRTL', 'UNIQUAC')
Model = 'DWPM'
MixtureDataDir = 'Data/Mixtures'
PureDataDir = 'Data/PureComps'
Compounds = ('methanenitro', 'nonanol')
InitParams = (-29.28946, -200.35796, 0.5)
Bounds = [(-1000, 0), (-1000, 0), (0.5, 0.5)]

a = Mixture(Compounds, MixtureDataDir, PureDataDir)
a.BestFitParams(Model, InitParams, Bounds)



#    for Method in Methods:
#        print('Fitting data using ', Method)
#        [Best_Params{m}, x_eq{m}] = CalcBestFitParam(R, T, M, Best_Params, c, delGmixfun, ddelGmix_dx, m )
        


#delGmix_plot(m,:) = arrayfun(@(x) delGmixfun{m}(R, T, M, Best_Params{m}, c, x), xplot);
#Gmix_eq(m,:) = arrayfun(@(x) delGmixfun{m}(R, T, M, Best_Params{m}, c, x), x_eq{m});

