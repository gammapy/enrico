import os
try:
    from packaging.version import Version
except ImportError:
    from distutils.version import LooseVersion as Version
import numpy as np
try:
    import astropy.io.fits as fits
except ImportError:
    import pyfits as fits
import pyLikelihood
import matplotlib
matplotlib.use('Agg')
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 15})
matplotlib.rc('text', usetex=False)
import matplotlib.pyplot as plt
from enrico.constants import MEV_TO_ERG, ERG_TO_MEV
from enrico.config import get_config
from enrico import utils
from enrico import Loggin
from enrico.extern.astropy_bayesian_blocks import bayesian_blocks

class Params:
    """Collection of Plotting parameters like Energy bounds,
    colors, file name, etc...."""
    def __init__(self, srcname, Emin=100, Emax=3e5,
                 PlotName="LAT_SED", LineColor=2,
                 PointColor = 1, N = 2000, SaveResData=False):
        self.Emin = Emin #Energy bounds
        self.Emax = Emax
        self.N = N #Number of points for the TGraph
        self.srcname = srcname # Source of interest
        self.PlotName = PlotName #file name
        #color options
        self.LineColor = LineColor
        self.PointColor = PointColor
        self.SaveResData = SaveResData


class Result(Loggin.Message):
    """Helper class to get the results from a (Un)BinnedAnalysis object
    and compute the SED and errors"""
    def __init__(self, Fit, pars):
        super(Result,self).__init__()
        Loggin.Message.__init__(self)
        
        self.Fit = Fit
        self.Model = Fit[pars.srcname].funcs['Spectrum'].genericName()
        self.ptsrc = pyLikelihood.PointSource_cast(Fit[pars.srcname].src)
        covar_matrix,covar_pars = utils.GetCovar(pars.srcname, self.Fit, verbose = False, with_par_map = True)
        self.covar = np.array(covar_matrix)
        self.covar_pars = np.array(covar_pars)
        self.srcpars = pyLikelihood.StringVector()
        Fit[pars.srcname].src.spectrum().getFreeParamNames(self.srcpars)

    def GetDecorrelationEnergy(self,par):
        self.E, self.SED = self.MakeSED(par)
        self.Err         = self.MakeSEDError(par)
        i=np.argmin(self.Err/self.SED)
        self.decE       = self.E[i]
        self.decFlux    = self.SED[i]/self.E[i]**2*ERG_TO_MEV
        self.decFluxerr = self.Err[i]/self.E[i]**2*ERG_TO_MEV
        self.decSED     = self.SED[i]
        self.decSEDerr  = self.Err[i]

    def _WriteFitPars(self,par):
        header  = '#### Parameter matrix. ###\n#Par Name, Value, Error, Scale:\n'
        spectrum = self.Fit[par.srcname].funcs['Spectrum']
        ParName = spectrum.paramNames
        OutFile = par.PlotName+'.fitpars.dat'
        data = []
        for p in ParName:
            ParValue = '{:.3e}'.format(float(spectrum.getParam(p).value()))
            ParError = '{:.3e}'.format(float(spectrum.getParam(p).error()))
            ParScale = '{:.3e}'.format(float(spectrum.getParam(p).getScale()))
            data.append([p,ParValue,ParError,ParScale])
        
        try:
            np.savetxt(OutFile, data, header=header, fmt='%s', comments='', delimiter=', ')    
        except FileNotFoundError:
            self.warning("Cannot write fitpars file: {}".format(OutFile))
        
        

    def _WriteCovMatrix(self,par):
        header  = '#### Covariance matrix. ###\n#Parameters:\n'
        header += ''.join(['#'+str(s)+'\n' for s in self.covar_pars])
        try:
            np.savetxt(par.PlotName+'.cov.dat', self.covar, header=header, fmt='%.3e', comments='', delimiter=', ')    
        except FileNotFoundError:
            self.warning("Cannot write cov file: {}".format(par.PlotName+'.cov.dat'))

    def _DumpSED(self,par):
        """Save the energy, E2.dN/dE, and corresponding  error in an ascii file
        The count and residuals plot vs E is also made"""

        try:
            self.decE
        except NameError:
            self.GetDecorrelationEnergy(par)

        self.info("Decorrelation energy : %4.2e MeV"% self.decE)
        self.info("Diffential flux  at the Decorrelation energy : %2.2e +/-  %2.2e ph/cm2/s/MeV" \
                %(self.decFlux, self.decFluxerr))
        self.info("SED value at the Decorrelation energy : %2.2e +/-  %2.2e erg/cm2/s" \
                %(self.decSED, self.decSEDerr))

        try:
            self.CountsPlot(par)
        except Exception as e:
            print((type(e)))    # the exception instance
            print((e.args))     # arguments stored in .args
            print(e)          # __str__ allows args to be printed directly,
            #raise
        
        # Save all in ascii file
        # log(E)  log (E**2*dN/dE)   log(E**2*dN/dE_err)  is_dot (0,1) is_upper (0,1)
        save_file = open(par.PlotName + '.dat', 'w')
        save_file.write("#E [MeV]  E**2*dN/dE [erg/cm2/s]  Error on (E**2*dN/dE) [erg/cm2/s]   \n")
        for i in range(par.N):
            save_file.write("%12.4e  %12.4e  %12.4e \n" % (self.E[i], self.SED[i], self.Err[i]))
        save_file.close()

    def MakeFlux(self, params):
        """Compute differential Flux distribution and
        corresponding energy and return a numpy array"""
        E = np.logspace(np.log10(params.Emin), np.log10(params.Emax), params.N)
        Flux = np.zeros(params.N)
        for i in range(params.N):
            Flux[i] = self.dNde(E[i])
        return E, Flux

    def MakeSED(self, pars):
        """Compute Spectral energy distribution and corresponding energy
        and return a numpy array"""
        E = np.logspace(np.log10(pars.Emin), np.log10(pars.Emax), pars.N)
        nuFnu = np.zeros(pars.N)
        for i in range(pars.N):
            nuFnu[i] = MEV_TO_ERG  * E[i] ** 2 * self.dNde(E[i]) #Mev to Ergs
        return E, nuFnu

    def MakeSEDError(self, pars):
        """@todo: document me"""
        estep = np.log(pars.Emax / pars.Emin) / (pars.N - 1)
        energies = pars.Emin * np.exp(estep * np.arange(np.float(pars.N)))
        err = np.zeros(pars.N)
        j = 0
        for ene in energies:
            arg = pyLikelihood.dArg(ene)
            partials = np.zeros(len(self.srcpars))
            for i in range(len(self.srcpars)):
                x = self.srcpars[i]
                partials[i] = self.ptsrc.spectrum().derivByParam(arg, x)
            err[j] = np.sqrt(np.dot(partials, np.dot(self.covar, partials)))
            j += 1

        return MEV_TO_ERG  * energies ** 2 * err #Mev to Ergs

    def dNde(self, energy):
        arg = pyLikelihood.dArg(energy)
        return self.ptsrc.spectrum()(arg)

    def CountsPlot(self, Parameter):
        """@todo: document me"""
        imName = "tmp.fits"
        filebase = Parameter.PlotName

        total   = np.array([])
        obs     = np.array([])
        obs_err = np.array([])
        emax    = np.array([])
        emin    = np.array([])
        src     = np.array([])

        # Summed Likelihood has no writeCountsSpectra
        # but we can do it component by component
        for comp in self.Fit.components:
            #self.Fit.writeCountsSpectra(imName)
            try:
                comp.writeCountsSpectra(imName)
                image = fits.open(imName)

                #loop on the source names to find the good one
                j = 0
                for ID in image[1].data.names:
                    if ID == Parameter.srcname:
                        indice = j
                    j += 1

                for jn in range(len(image[3].data.field(0))):
                    energymin = image[3].data.field(1)[jn]
                    energymax = image[3].data.field(0)[jn]
                    if energymax in emax and energymin in emin:
                        k = np.where(energymax==emax)
                        obs[k]     = obs[k] + image[1].data.field(0)[jn]
                        obs_err[k] = np.sqrt(obs[k])
                        src[k]     = src[k] + image[1].data.field(indice)[jn]
                        for i in range(len(image[1].data.names) - 1):
                            total[k] = total[k] + image[1].data.field(i + 1)[jn]
                    else:
                        emax    = np.append(emax, max(energymax,energymin))
                        emin    = np.append(emin, min(energymax,energymin))
                        obs     = np.append(obs,image[1].data.field(0)[jn])
                        obs_err = np.append(obs_err,\
                                            np.sqrt(image[1].data.field(0)[jn]))
                        src     = np.append(src, image[1].data.field(indice)[jn])
                        total   = np.append(total,0)
                        for i in range(len(image[1].data.names) - 1):
                            total[-1] = total[-1] + image[1].data.field(i + 1)[jn]
            except RuntimeError as e:
                print("Exception RuntimeError ocurred: ")
                print((type(e)))
                print((e.args))
                print(e)
                break
            except IndexError:
                print("Exception IndexError ocurred (component unavailable): ")
                print((type(e)))
                print((e.args))
                print(e)
                continue

        # Sort by energy 
        energy_order = np.argsort(emin)
        src     = src[energy_order]
        obs     = obs[energy_order]
        obs_err = obs_err[energy_order]
        total   = total[energy_order]
        emin    = emin[energy_order]
        emax    = emax[energy_order]

        other = np.array(total - src)
        Nbin  = len(src)
        E = 10**((np.log10(emin)+np.log10(emax))/2.)
        #E = np.asarray((emax + emin) / 2.)
        err_E = np.asarray([E-emin,emax-E])
        total = np.asarray(total)
        residual = np.zeros(Nbin)
        Dres = np.zeros(Nbin)

        srcname = Parameter.srcname.replace("_"," ")

        print('Generating counts plot')
        plt.figure()
        plt.loglog()
        plt.title('Counts plot')
        plt.xlabel("E (MeV) ")
        plt.ylabel("Counts / bin")
        plt.errorbar(E,obs,xerr=err_E,yerr=obs_err,fmt='o',color="red",ls='None',label="Data")
        plt.plot(E,src,ls='dashed',color="blue",label=srcname)
        plt.plot(E,other,ls='solid',color="green",label="Other Sources")
        plt.plot(E,total,lw=1.5,ls='solid',label="All Sources")
        plt.legend()
        plt.tight_layout()
        plt.savefig(filebase + "_CountsPlot.png", dpi=150, facecolor='w', edgecolor='w',
            orientation='portrait', format=None,
            transparent=False, pad_inches=0.1)

        print('Generating residuals plot')
        plt.figure()
        plt.title('Residuals plot')
        plt.semilogx()
        for i in range(Nbin):
            try:
                residual[i] = (obs[i] - total[i]) / total[i]
                Dres[i] = (obs_err[i] / total[i])
            except:
                residual[i] = 0.
                Dres[i] = 0.
            if residual[i] == -1.:
               residual[i] = 0.

        # Write residuals to csv file
        if Parameter.SaveResData:
            print('Writing also residuals data')
            residual_array = np.asarray([E,err_E[0],err_E[1],obs,obs_err,residual,Dres]).transpose()
            np.savetxt(filebase + '.ResData.dat', residual_array, 
                       header='E, E_err-, E_err+, Counts, Counts_err, Residuals,Residuals_err', fmt='%.3e', delimiter=',')
        
        ymin = min(residual) - max(Dres)
        ymax = max(residual) + max(Dres)
        plt.ylim(ymax = ymax, ymin = ymin)
        plt.xlim(xmin = min(E)*0.3, xmax = max(E)*2)
        plt.xlabel("E (MeV) ")
        plt.ylabel("(counts-model)/model")
        plt.errorbar(E,residual,xerr=err_E,yerr=Dres,fmt='o',color="red",ls='None',label="Data")
        zero = np.zeros(2)
        Ezero = np.array([1e-5, 1e10])
        plt.plot(Ezero,zero,lw=1.5,ls='solid',color='black')
        plt.tight_layout()
        plt.savefig(filebase + "ResPlot.png", dpi=150, facecolor='w', 
            edgecolor='w', orientation='portrait', format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1)
        os.system("rm " + imName)
        image.close()

def GetDataPoints(config,pars,ignore_missing_bins=False):
    """Collect the data points/UL and generate a TGraph for the points
    and a list of TArrow for the UL. All is SED format"""

    #Preparation + declaration of arrays
    arrows = []
    NEbin = int(config['Ebin']['NumEnergyBins'])
    lEmax = np.log10(float(config['energy']['emax']))
    lEmin = np.log10(float(config['energy']['emin']))
    Epoint = np.zeros(NEbin)
    EpointErrp = np.zeros(NEbin)
    EpointErrm = np.zeros(NEbin)
    Fluxpoint = np.zeros(NEbin)
    FluxpointErrp = np.zeros(NEbin)
    FluxpointErrm = np.zeros(NEbin)
    uplim = np.zeros(NEbin,dtype=int)
    ener = np.logspace(lEmin, lEmax, NEbin + 1)

    mes = Loggin.Message()
    mes.info("Save Ebin results in ",pars.PlotName+".Ebin.dat")
    dumpfile = open(pars.PlotName+".Ebin.dat",'w')
    dumpfile.write("# Energy (MeV)\tEmin (MeV)\tEmax (MeV)\tE**2. dN/dE (erg.cm-2s-1)\tGaussianError\tMinosNegativeError\tMinosPositiveError\n")

    from enrico.constants import EbinPath
    for i in range(NEbin):#Loop over the energy bins
        #E = int(pow(10, (np.log10(ener[i + 1]) + np.log10(ener[i])) / 2))
        filename = (config['out'] + '/'+EbinPath+str(NEbin)+'/' + config['target']['name'] +
                    "_" + str(i) + ".conf")

        try:#read the config file of each data points
            CurConf = get_config(filename)
            mes.info("Reading "+filename)
            results = utils.ReadResult(CurConf)
        except:
            if not ignore_missing_bins:
                mes.warning("cannot read the Results of energy bin "+ str(i))
            continue
        #fill the energy arrays
        #Epoint[i] = results.get("Scale") 
        #if Epoint[i] in [results.get("Emin"),results.get("Emax")]: 
            #### <---- is this a mistake?? does not make much sense to me
            Epoint[i] = 10**((np.log10(results.get("Emin"))+np.log10(results.get("Emax")))/2.)
            #Epoint[i] = int(pow(10, (np.log10(ener[i + 1]) + np.log10(ener[i])) / 2))
        
        Epoint[i] = 10**((np.log10(results.get("Emin"))+np.log10(results.get("Emax")))/2.)

        EpointErrm[i] = Epoint[i] - results.get("Emin")
        EpointErrp[i] = results.get("Emax") - Epoint[i]
        dprefactor = 0

        #Compute the flux or the UL (in SED format)
        if 'Ulvalue' in results:
            PrefUl = utils.Prefactor(results.get("Ulvalue"),results.get("Index"),
                                    results.get("Emin"),results.get("Emax"),Epoint[i])
            Fluxpoint[i] = MEV_TO_ERG  * PrefUl * Epoint[i] ** 2
            uplim[i] = 1
        else : #Not an UL : compute points + errors
            Fluxpoint[i] = MEV_TO_ERG  * results.get("Prefactor") * Epoint[i] ** 2

        dprefactor = results.get("dPrefactor")
        try:
            down = abs(results.get("dPrefactor-"))
            up = results.get("dPrefactor+")
            if down==0 or  up ==0 :
              mes.error("cannot get Error value")
            FluxpointErrp[i] = MEV_TO_ERG  * up * Epoint[i] ** 2
            FluxpointErrm[i] = MEV_TO_ERG  * down * Epoint[i] ** 2
        except:
            try:
                err = MEV_TO_ERG  * dprefactor * Epoint[i] ** 2
                FluxpointErrp[i] = err
                FluxpointErrm[i] = err
            except:
                pass

        mes.info("Energy bins results")
        print(("Energy = ",Epoint[i]))
        #Save the data point in a ascii file
        if 'Ulvalue' in results:
            dumpfile.write(str(Epoint[i])+"\t"+\
                           str(results.get("Emin"))+"\t"+\
                           str(results.get("Emax"))+"\t"+\
                           str(Fluxpoint[i])+"\t0\t0\t0\n")
            print(("E**2. dN/dE = ",Fluxpoint[i]))
        else:
            dumpfile.write(str(Epoint[i])+"\t"+\
                           str(results.get("Emin"))+"\t"+\
                           str(results.get("Emax"))+"\t"+\
                           str(Fluxpoint[i])+"\t"+\
                           str(MEV_TO_ERG  * dprefactor * Epoint[i] ** 2)+"\t"+\
                           str(FluxpointErrm[i])+"\t"+\
                           str(FluxpointErrp[i])+"\n")
            print(("E**2. dN/dE = ",Fluxpoint[i]," + ",FluxpointErrp[i]," - ",FluxpointErrm[i]))
    dumpfile.close()
    return Epoint, Fluxpoint, EpointErrm, EpointErrp, FluxpointErrm, FluxpointErrp, uplim

def plot_errorbar_withuls(x,xerrm,xerrp,y,yerrm,yerrp,uplim,bblocks=False):
    """ plot an errorbar plot with upper limits. Optionally compute and draw bayesian blocks (bblocks) """
    # plt.errorbar(Epoint, Fluxpoint, xerr=[EpointErrm, EpointErrp], yerr=[FluxpointErrm, FluxpointErrp],fmt='o',color='black',ls='None',uplims=uplim)
    uplim = np.asarray(uplim,dtype=bool) # It is an array of 1 and 0s, needs to be a bool array.
    # make sure that the arrays are numpy arrays and not lists.
    x = np.asarray(x)
    xerrm = np.asarray(xerrm)
    xerrp = np.asarray(xerrp)
    y = np.asarray(y)
    yerrm = np.asarray(yerrm)
    yerrp = np.asarray(yerrp)
    # Get the strict upper limit (best fit value + error, then set the error to 0 and the lower error to 20% of the value)
    y[uplim] += yerrp[uplim]
    yerrm[uplim] = 0
    yerrp[uplim] = 0

    optimal_markersize = (0.5+4./(1.+np.log10(len(y))))
    optimal_errorlinewidth = (0.2+2./(1.+4.*np.log10(len(y))))

    # Plot the significant points
    plt.errorbar(x[~uplim], y[~uplim],
        xerr=[xerrm[~uplim], xerrp[~uplim]],
        yerr=[yerrm[~uplim], yerrp[~uplim]],
        lw=optimal_errorlinewidth,
        fmt='o',ms=optimal_markersize,capsize=0,zorder=10,
        color='black',ls='None',uplims=False,label='LAT data')

    # Plot the upper limits. For some reason, matplotlib draws the arrows inverted for uplim and lolim [?]
    # This is a known issue fixed in matplotlib 1.4: https://github.com/matplotlib/matplotlib/pull/2452
    if Version(matplotlib.__version__) < Version("1.4.0"):
        plt.errorbar(x[uplim], y[uplim],
            xerr=[xerrm[uplim], xerrp[uplim]],
            yerr=[yerrm[uplim], yerrp[uplim]],
            fmt='o',markersize=0,capsize=0,zorder=-1,
            lw=optimal_errorlinewidth,
            color='0.50',ls='None',lolims=False)
        plt.errorbar(x[uplim], 0.8*y[uplim],
            yerr=[0.2*y[uplim], 0.2*y[uplim]],
            fmt='o',markersize=0,capsize=optimal_markersize/1.5,zorder=-1,
            lw=optimal_errorlinewidth,
            color='0.50',ls='None',lolims=True)
    else:
        plt.errorbar(x[uplim], y[uplim],
            xerr=[xerrm[uplim], xerrp[uplim]],
            yerr=[yerrm[uplim], yerrp[uplim]],
            lw=optimal_errorlinewidth,
            fmt='o',markersize=0,capsize=0,zorder=-1,
            color='0.50',ls='None',uplims=False)
        plt.errorbar(x[uplim], y[uplim],
            yerr=[0.2*np.abs(y[uplim])+1e-40, 0.2*np.abs(y[uplim])+1e-40],
            lw=optimal_errorlinewidth,
            fmt='o',markersize=0,capsize=optimal_markersize/1.5,zorder=-1,
            color='0.50',ls='None',uplims=True)

    if bblocks and len(x[~uplim])>2:
        yerr = 0.5*(yerrm+yerrp)
        # Set the value and error for the uls.
        yerr[uplim] = y[uplim] #min(y[yerr>0]+yerr[yerr>0])
        y[uplim] = 0
        edges = bayesian_blocks(x,y,yerr,fitness='measures',p0=0.5)
        #edges = bayesian_blocks(x[yerr>0],y[yerr>0],yerr[yerr>0],fitness='measures',p0=0.1)
        xvalues = 0.5*(edges[:-1]+edges[1:])
        xerrors = 0.5*(edges[1:]-edges[:-1])
        yvalues = []
        yerrors = []
        for k in range(len(edges)-1):
            xmin,xmax = edges[k],edges[k+1]
            filt = (x>=xmin)*(x<=xmax)*(yerr>0)
            sum_inv_square = np.sum(1./yerr[filt]**2)
            yvalues.append(np.sum(y[filt]/yerr[filt]**2)/sum_inv_square)
            yerrors.append(1./np.sqrt(sum_inv_square))

        yvalues = np.asarray(yvalues)
        yerrors = np.asarray(yerrors)

        # Plot the significant points
        ystep = []
        ystepmin = []
        ystepmax = []
        xstep = []
        for k in range(len(xvalues)):
            for _ in range(2):
                ystep.append(yvalues[k]) # 3 values, to mark the minimum and center
                ystepmin.append(yvalues[k]-yerrors[k]) # 3 values, to mark the minimum and center
                ystepmax.append(yvalues[k]+yerrors[k]) # 3 values, to mark the minimum and center
            xstep.append(xvalues[k]-xerrors[k])
            xstep.append(xvalues[k]+xerrors[k])

        plt.step(xstep, ystep,
            color='#d62728',zorder=-10,
            ls='solid')
        plt.fill_between(xstep, ystepmin, ystepmax,
            color='#d62728',zorder=-10, alpha=0.5)
        plt.errorbar(xvalues, yvalues,
            xerr=xerrors,yerr=yerrors,
            marker=None,ms=0,capsize=0,color='#d62728',zorder=-10,
            ls='None',label='bayesian blocks')

        plt.legend(loc=0,fontsize='small',numpoints=1)

def plot_bayesianblocks(xmin, xmax, y, yerrm, yerrp, uplim):
    # Set the value and error for the uls.
    yerrm[uplim] = y[uplim]
    yerrp[uplim] = y[uplim]
    y[uplim] = 0.

    xvalues = 0.5*(xmax+xmin)
    xerrors = 0.5*(xmax-xmin)

    # Plot the significant points
    ystep = []
    ystepmin = []
    ystepmax = []
    xstep = []
    for k in range(len(xvalues)):
        for _ in range(2):
            ystep.append(y[k]) # 3 values, to mark the minimum and center
            ystepmin.append(y[k]-yerrm[k]) # 3 values, to mark the minimum and center
            ystepmax.append(y[k]+yerrp[k]) # 3 values, to mark the minimum and center
        xstep.append(xmin[k])
        xstep.append(xmax[k])

    plt.step(xstep, ystep,
             color='#d62728',zorder=-10,
             ls='solid')
    plt.fill_between(xstep, ystepmin, ystepmax,
                     color='#d62728',zorder=-10, alpha=0.5)
    plt.errorbar(xvalues, y,
                 xerr=xerrors,yerr=[yerrm, yerrp],
                 marker=None,ms=0,capsize=0,color='#d62728',zorder=-10,
                 ls='None')




def PlotSED(config,pars,ignore_missing_bins=False):
    """plot a nice SED with a butterfly and points"""

    # Read the ascii file where the butterfly is stored
    filebase = pars.PlotName#utils._SpecFileName(config)
    print("\033[33m PrepareEbin",filebase," \033[0m")
    lines = open(filebase + '.dat', 'r').readlines()
    SED = []
    E = []
    Err = []
    for i in range(len(lines) - 1):
        words = lines[i + 1].split()
        if float(words[0])<pars.Emax :
            E.append(float(words[0]))
            SED.append(float(words[1]))
            Err.append(float(words[2]))
    ilen = len(SED)

    #From dN/dE to SED
    Fluxp = np.array(SED)*np.exp(np.array(Err)/np.array(SED))
    Fluxm =  np.array(SED)*np.exp(-np.array(Err)/np.array(SED))
    ErrorFlux = np.zeros(2 * ilen + 1)
    ErrorE = np.zeros(2 * ilen + 1)

    #Compute the butterfly and close it
    for i in range(ilen):
        ErrorFlux[i] = Fluxp[i]
        ErrorE[i] = E[i]
    for i in range(ilen):
        ErrorFlux[ilen + i] = Fluxm[ilen - i - 1]
        ErrorE[ilen + i] = E[ilen - i - 1]
    ErrorFlux[-1] = Fluxp[0]
    ErrorE[-1] = E[0]

    #Actually make the plot
    plt.figure()
    plt.title(pars.PlotName.split("/")[-1])
    name = pars.PlotName.split("/")[-1]
    plt.loglog()

    plt.xlabel(r"Energy (MeV)")
    plt.ylabel(r"$\mathrm{E^2\ dN/dE}\ \mathrm{(erg\ cm^{-2} s^{-1})}$")
    plt.plot(E,SED,"-r",label='LAT model')
    plt.plot(ErrorE,ErrorFlux,"-r")

    #Plot points
    NEbin = int(config['Ebin']['NumEnergyBins'])
    if NEbin > 0:
        Epoint, Fluxpoint, EpointErrm, EpointErrp, FluxpointErrm, FluxpointErrp, uplim = GetDataPoints(config,pars,ignore_missing_bins) #collect data points
        plot_errorbar_withuls(Epoint,EpointErrm,EpointErrp,Fluxpoint,FluxpointErrm,FluxpointErrp,uplim)

    #print uplim
    #print FluxpointErrm
    #print FluxpointErrp

    #Set meaningful axes limits
    xlim = plt.xlim()
    ylim = plt.ylim()
    xlim = (max([20,xlim[0]]),min([2e6,xlim[1]]))
    ylim = (max([1e-14,ylim[0]]),min([1e-8,ylim[1]]))
    plt.xlim(xlim)
    plt.ylim(ylim)
    # turn them into log10 scale
    #xticks = plt.xticks()[0]
    #xticklabels = np.array(np.log10(xticks),dtype=int)
    #plt.xticks(xticks,xticklabels)
    #plt.xlabel('$\mathrm{\log_{10}\mathbf{(Energy)} \\ \\ [MeV]}$')

    plt.legend(fontsize='small',ncol=1,\
               loc=3,numpoints=1)#,framealpha=0.75)


    #Upper horizontal secondary axis with frequency
    #Plt2 = plt.twiny()
    #Plt2.set_xscale('log')
    #Plt2.set_xlim(2.417990504024163e+20 *np.array(xlim))
    #Plt2.set_xticklabels(np.array(np.log10(Plt2.get_xticks()),dtype=int))
    #Plt2.set_xlabel('$\mathrm{\log_{10}\mathbf{(Frequency)} \\ \\ [Hz]}$')

    #save the canvas
    #plt.grid()
    plt.tight_layout()
    plt.savefig("%s.png" %filebase, dpi=150, facecolor='w', edgecolor='w',
            orientation='portrait', format=None,
            transparent=False, pad_inches=0.1)

def PlotUL(pars,config,ULFlux,Index):

    #Compute the SED
    E = np.logspace(np.log10(pars.Emin), np.log10(pars.Emax), pars.N)
    SED = MEV_TO_ERG  * E ** 2 * (-Index+1)*ULFlux* np.power(E,-Index)/(np.power(pars.Emax,-Index+1)-np.power(pars.Emin,-Index+1))

    #Actually make the plot
    plt.xlabel(r"E [MeV]")
    plt.ylabel(r"$\mathrm{E^2\ dN/dE}\ \mathrm{(erg\ cm^{-2} s^{-1})}$")
    plt.loglog()
    plt.plot(E,SED,"-",color='black')

    # Plot the upper limits. For some reason, matplotlib draws the arrows inverted for uplim and lolim [?]
    # This is a known issue fixed in matplotlib 1.4: https://github.com/matplotlib/matplotlib/pull/2452
    if Version(matplotlib.__version__) < Version("1.4.0"):
        plt.errorbar([E[0],E[-1]], [SED[0],SED[-1]],  yerr=[SED[0]*0.8,SED[-1]*0.8],fmt='.',color='black',ls='None',lolims=[1,1])
    else:
        plt.errorbar([E[0],E[-1]], [SED[0],SED[-1]],  yerr=[SED[0]*0.8,SED[-1]*0.8],fmt='.',color='black',ls='None',uplims=[1,1])

    #save the plot
    filebase = utils._SpecFileName(config)
    plt.tight_layout()
    plt.savefig(filebase + '.png', dpi=150, facecolor='w', edgecolor='w',
            orientation='portrait', format=None,
            transparent=False, pad_inches=0.1)


def plot_sed_fromconfig(config,ignore_missing_bins=False):
    config = get_config(config)
    utils.mkdir_p(config["out"]+"/Spectrum")
    srcname = config['target']['name']
    Emin = config['energy']['emin']
    Emax = config['energy']['emax']
    filename = utils._SpecFileName(config)
    Param = Params(srcname, Emin=Emin, Emax=Emax, PlotName=filename)
    Result = utils.ReadResult(config)

    # if the TS > ts limit plot the butterfly, if not draw UL
    if Result["TS"]> config['UpperLimit']['TSlimit'] :
        # try :
        PlotSED(config,Param,ignore_missing_bins)
        # except:
            # pass
    else :
        try :
            PlotUL(Param,config,Result['Ulvalue'],config['UpperLimit']['SpectralIndex'])
        except :
            print("Not able to plot an upper limit in a SED diagram. UL computed?")
