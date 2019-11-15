import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import norm
import CrossCorr_Cells

#THESE ARE FILES FOR CLXC interpolation
clxcfolder = '../Data/grid/'
Om_baryons = 0.0493
num_Om_M = 60
Omarray = np.round(np.linspace(0.1, 0.7, num_Om_M), 4)
lOm = len(Omarray)
As = 2.4736e-9


# This is the kernel function
def K(z, K0_bsqr, z_pivot, Om, bGW):
    # Returns K(r) and K(r)*r^2 * nbar(r) accoding to our model.
    result1 = np.tanh( (-z + z_pivot)*10 )*K0_bsqr/2. + K0_bsqr/2.
    denom = FlatLambdaCDM(Om0=Om, H0=67.32).comoving_distance(z).to('Mpc').value
    #denom = denom + 1e-5#Add a small number so that the inverse is not INF
    denom[-1] = denom[-2]
    result2 = result1/denom**2./0.1/bGW/bGW

    return result1, result2

# Two window functions, very simple Normal dist.
def Windows(z, Wpars):
    std1, avg1, std3, avg3 = Wpars
    return norm.pdf(z, scale=std1, loc=avg1), norm.pdf(z, scale=std3, loc=avg3)

# This is the Cl function
def Cl(larray, Karray, Warray, clxcarray, bGW, spectrum = 2, Fortran = True):

    if(Fortran):
        a_unique_lst, k_unique_lst, eta_unique_lst, bess_int_keta_reshaped, confH, sigma, \
        deltaCDM, k_length, eta_length = clxcarray

        # We force bGW = bgal = 1
        z = 1/np.array(a_unique_lst) - 1.

        bias = np.ones(eta_length)*bGW

        jacobian = confH[0]/a_unique_lst

        if spectrum == 1: #autocorrelation
            W_kernel = Warray
            K_kernel = Karray
        else: #crosscorrelation
            W_kernel = Warray*jacobian
            K_kernel = Karray

        Clarray = np.ones(len(larray))

        CrossCorr_Cells.c_ells(Clarray, larray, k_unique_lst, eta_unique_lst, bess_int_keta_reshaped,\
                               confH, sigma, deltaCDM, W_kernel, K_kernel, bias, len(larray), k_length, eta_length)


    return As*4.*np.pi*Clarray

# This is the nearest neighbour in the Om grid, next best thing
def clxc(Om, Fortran=True):
    idx = np.abs(Om-Omarray).argmin()
    if(Fortran):
        return np.load(clxcfolder+'F/'+str(idx)+'.npy')
