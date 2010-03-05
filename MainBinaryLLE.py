#!/usr/bin/env python

from scipy import *
from pylab import *
from os import mkdir, path, remove, listdir
import scipy.optimize
import glob
import tables
import vanDerWaals
import ErrorClasses
import GibbsClasses
import PhaseStability


class Mixture:
        
    def __init__(self, Compounds, MixtureDataDir, PureDataDir):
        
        MixtureName = Compounds[0]
        self.Compounds = Compounds
        self.Data =dict((Compound, {}) for Compound in Compounds)
        for Compound in Compounds:
            h5file = tables.openFile(PureDataDir+'/'+Compound+'.h5', 'r')
            properties = h5file.root.Properties
            self.Data[Compound] = dict(((field, row[field]) for row in properties.iterrows() for field in properties.colnames))        
        self.Name = '-'.join(Compounds)
        h5file = tables.openFile(MixtureDataDir+'/'+ self.Name +'.h5', 'r')
        self.M = dict(Compounds = h5file.root.Compounds.read(), ExpComp = h5file.root.ExperimentalData.ExpComp.read(), T = h5file.root.ExperimentalData.T.read())
        for row in h5file.root.UNIQUACParams.iterrows():
            for field in h5file.root.UNIQUACParams.colnames:
                self.M[field] = row[field]
        self.vdWaalsInstance = vanDerWaals.Slope(self.Data, Compounds)
        [Slope, VPData] = self.vdWaalsInstance.BestFit()
        self.AdachiLuParam = dict((Compound,  Slope[Compound]['BestFit']) for Compound in Compounds)
    
    def Plotter(self, BestParams, Model, Fit, ModelInstance, Actual, c, T):
    
        deltaGibbsmix = array([ModelInstance.deltaGmix(x, T, c, self.M) for x in arange(0.001, 1, 0.001)])
        TangentComps = PhaseStability.CalcPhaseStability(ModelInstance, T, c, self.M)
        TangentGibbs = [ModelInstance.deltaGmix(x, T, c, self.M)for x in TangentComps]
        AbsError = ErrorClasses.SumAbs(TangentComps, Actual).Error()
        matplotlib.rc('text', usetex = True)
        fig = matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(arange(0.001, 1, 0.001), deltaGibbsmix, 'b-')
        matplotlib.pyplot.plot(TangentComps, TangentGibbs, 'r-')
        matplotlib.pyplot.xlabel(r'Mole Fraction of '+self.Compounds[0].capitalize(), fontsize = 14)
        matplotlib.pyplot.ylabel(r'$\displaystyle\frac{\Delta G_{mix}}{RT}$', fontsize = 14)
        matplotlib.pyplot.title(r'\textbf{Predicted Phase Equilibrium at %3.2f}'%T, fontsize = 14)
        ax = fig.add_axes([0,0,1,1])
        ax.text(0,0, 'Model: %s, Params: %s, Abslote Error(sum): %4.4e'%(Model, str(around(BestParams, 2)),AbsError), fontsize=12, transform=ax.transAxes)
        ax.set_axis_off()
        #show()
        savefig('Results/'+self.Name+'/'+Model+'/'+ Fit+'/T_'+str(T) +'.pdf')
        matplotlib.pyplot.close()
    
    def OptFunctionIndvT(self, params, Model, Actual, T, c):
        
        ModelInstance = getattr(GibbsClasses, Model)(params)        
        Predicted = PhaseStability.CalcPhaseStability(ModelInstance, T, c, self.M)
        Error = ErrorClasses.SumSquare(Predicted ,Actual).Error()

        return Error
  
    def NonEqConstrIndvT(self, params, Model, Actual, T, c):
        
        ModelInstance = getattr(GibbsClasses, Model)(params) 
        deltaGibbsMixTest = array([ModelInstance.deltaGmix(x, T, c, self.M) for x in arange(0.001, 1, 0.001)])
        TangentTestComps = PhaseStability.CalcPhaseStability(ModelInstance, T, c, self.M)
        TangentTest = array([((ModelInstance.FirstDerivative(TangentTestComps[0], T, c, self.M))*(x - TangentTestComps[0])+ ModelInstance.deltaGmix(TangentTestComps[0], T, c, self.M))-ModelInstance.deltaGmix(x, T, c, self.M) for x in arange(0.001, 1, 0.001)])
        DistinctTest = -1*array([abs(TangentTestComps[0]-TangentTestComps[1])+0.01])
       
        return -1*append(append(deltaGibbsMixTest, TangentTest), DistinctTest)

    def BestFitParamsIndvT(self, Model, InitParams, Bounds):
        Fit = 'IndividualT'
        if not(path.exists('Results/'+self.Name+'/'+Model)):
            mkdir('Results/'+self.Name+'/'+Model)        
        if path.exists('Results/'+self.Name+'/'+Model+'/'+Fit):
            fileList = listdir('Results/'+self.Name+'/'+Model+'/'+Fit)
            for fn in fileList: 
                remove(path.join('Results/'+self.Name+'/'+Model+'/'+Fit, fn))
        else:
            mkdir('Results/'+self.Name+'/'+Model+'/'+Fit)
        
        for T in self.M['T']:
                    
            Actual =  array([interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][0])),interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][1]))])
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in self.Compounds]
            
            [params, fx, its, imode, smode] = scipy.optimize.fmin_slsqp(self.OptFunctionIndvT, InitParams,[], None, [], self.NonEqConstrIndvT, Bounds, None, None, None, (Model, Actual, T, c), 500, 10e-8, 1, 1, 10e-8)
            InitParams = params
            ModelInstance = getattr(GibbsClasses, Model)(params)      
            self.Plotter(params, Model, Fit, ModelInstance, Actual, c, T)
            
        return params
    
    def OptFunctionOvrlT(self, params, Model,R):

        ModelInstance = getattr(GibbsClasses, Model)(params) 
        Errors = zeros(size(self.M['T']))
        
        for T in self.M['T']:
            Actual =  array([interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][0])),interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][1]))])
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in Compounds]
            Predicted = PhaseStability.CalcPhaseStability(ModelInstance, T, c, self.M)
            Errors[where(self.M['T']== T)] = ErrorClasses.SumAbs(Predicted ,Actual).Error()
        
        OvrlError = sum(Errors**2)

        return OvrlError

    def NonEqConstrOvrlT(self, params, Model,R):
        
        ModelInstance = getattr(GibbsClasses, Model)(params) 
        Test = array([])
        
        for T in self.M['T']:
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in Compounds]
            deltaGibbsMixTest = array([ModelInstance.deltaGmix(x, T, c, self.M) for x in arange(0.001, 1, 0.001)])
            TangentTestComps = PhaseStability.CalcPhaseStability(ModelInstance, T, c, self.M)
            TangentTest = array([((ModelInstance.FirstDerivative(TangentTestComps[0], T, c, self.M))*(x - TangentTestComps[0])+ ModelInstance.deltaGmix(TangentTestComps[0], T, c, self.M))-ModelInstance.deltaGmix(x, T, c, self.M) for x in arange(0.001, 1, 0.001)])
            DistinctTest = -1*array([abs(TangentTestComps[0]-TangentTestComps[1])+0.01])
            Test = append(Test, append(append(deltaGibbsMixTest, TangentTest), DistinctTest))
       
        return -1*Test
        
    def BestFitParamsOvrlT(self, Model, InitParams, Bounds):
        
        Fit = 'OverallT'
        if not(path.exists('Results/'+self.Name+'/'+Model)):
            mkdir('Results/'+self.Name+'/'+Model)        
        if path.exists('Results/'+self.Name+'/'+Model+'/'+Fit):
            fileList = listdir('Results/'+self.Name+'/'+Model+'/'+Fit)
            for fn in fileList: 
                remove(path.join('Results/'+self.Name+'/'+Model+'/'+Fit, fn))
        else:
            mkdir('Results/'+self.Name+'/'+Model+'/'+Fit)

        [params, fx, its, imode, smode] = scipy.optimize.fmin_slsqp(self.OptFunctionOvrlT, InitParams,[], None, [], self.NonEqConstrOvrlT, Bounds, None, None, None, (Model, R), 100, 10e-4, 1, 1, 10e-2)
       
        ModelInstance = getattr(GibbsClasses, Model)(params)      
        for T in self.M['T']:
            Actual =  array([interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][0])),interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][1]))])
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in self.Compounds]
            self.Plotter(params, Model, Fit, ModelInstance, Actual, c, T)

        return params

    
       

##=============================================================##
Models = ('DWPM', 'NRTL', 'UNIQUAC')
MixtureDataDir = 'Data/Mixtures'
PureDataDir = 'Data/PureComps'
Compounds = ('methanenitro', 'nonanol')
Bounds = [((-1000, 0), (-1000, 0), (0.5, 0.5)), ((-800, 800), (-800, 800)), ((-800, 800), (-800, 800))]
InitParams = [(-29.0,-200.0, 0.5), (700.0, 173.3), (35.00, 370.00)]
Model = Models[0]
R = 8.314

if not(path.exists('Results/')):
    mkdir('Results/')        

a = Mixture(Compounds, MixtureDataDir, PureDataDir)
if not(path.exists('Results/'+a.Name)):
    mkdir('Results/'+a.Name)
a.BestFitParamsIndvT(Model, InitParams[0], Bounds[0])
a.BestFitParamsOvrlT(Model, InitParams[0], Bounds[0])



#    for Method in Methods:
#        print('Fitting data using ', Method)
#        [Best_Params{m}, x_eq{m}] = CalcBestFitParam(R, T, M, Best_Params, c, delGmixfun, ddelGmix_dx, m )
        


#delGmix_plot(m,:) = arrayfun(@(x) delGmixfun{m}(R, T, M, Best_Params{m}, c, x), xplot);
#Gmix_eq(m,:) = arrayfun(@(x) delGmixfun{m}(R, T, M, Best_Params{m}, c, x), x_eq{m});

