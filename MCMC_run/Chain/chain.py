import numpy as np
from aux import K, Windows, Cl, clxc
from multiprocessing import Pool
import sys
sys.path.insert(1, '../../')
import emcee



filename = '../../Data/chains/No_Om_prior/chain_5.h5'
n_cores = 10
max_n = 50000

# load files
Wpars = np.load('../data/window.npy')
fiducial1, fiducial2, fiducial3 = np.load('../data/fiducial.npy')
larray = np.load('../data/larray.npy')


def normal(x, c, s):
    return -((x-c)**2.)/2./(s**2.)

def lnlike(params):
    b, K0_bsqr, z_pivot, Om = params

    clxcarray = clxc(Om)

    a = np.array(clxcarray[0])
    z = 1./a - 1.
    KarrayCl, Karrayn = K(z, K0_bsqr, z_pivot, Om, b)
    W1array, W2array, W3array = Windows(z, Wpars)


    Bl1 = np.abs(np.trapz(Karrayn*W1array,z)) # noise bias
    Bl2 = np.abs(np.trapz(Karrayn*W2array,z)) # noise bias
    Bl3 = np.abs(np.trapz(Karrayn*W3array,z)) # noise bias
    Cl1 = Cl(larray, KarrayCl, W1array, clxcarray, b, True) # Cl
    Cl2 = Cl(larray, KarrayCl, W2array, clxcarray, b, True) # Cl
    Cl3 = Cl(larray, KarrayCl, W3array, clxcarray, b, True) # Cl

    return (-0.5 * ( ((Cl1 + Bl1) - fiducial1)**2./(Cl1 + Bl1)**2.  \
                   + ((Cl2 + Bl2) - fiducial2)**2./(Cl2 + Bl2)**2.  \
                   + ((Cl3 + Bl3) - fiducial3)**2./(Cl3 + Bl3)**2.)* \
                   (2.*larray +1.)/2. ).sum()


def lnpost(params):
    b, K0_bsqr, z_pivot, Om = params

    if( (1e-2 < b < 100.) and (1e-2 < K0_bsqr < 100.) and (0.5 < z_pivot < 1.5) and (0.1 < Om < 0.7)):
        return lnlike(params)
    return -np.inf


print "initialize chain"
p0 = np.array([0.5, 0.5, 0.8, 0.5])
ndim, nwalkers = len(p0), n_cores-1
pos = [p0 + p0*np.random.randn(ndim)*1e-1 for i in range(nwalkers)]
backend = emcee.backends.HDFBackend(filename)
first_time = False


try:
    backend.get_chain()
except AttributeError:
    first_time = True
    print "first time"

print "start sampling"
pool = Pool(n_cores)
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnpost, backend=backend, pool=pool)
if(first_time):
    sampler.run_mcmc(pos, nsteps=max_n)
else:
    sampler.run_mcmc(pos0=None, nsteps=max_n)
