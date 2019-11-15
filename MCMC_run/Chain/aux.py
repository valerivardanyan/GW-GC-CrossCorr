import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import norm
import CrossCorr_Cells

#THESE ARE FILES FOR CLXC interpolation
clxcfolder = '../data/grid/'
Omarray = np.load('../data/Om_grid.npy')
lOm = len(Omarray)
As = 2.4736e-9


# This is the kernel function
def K(z, K0_bsqr, z_pivot, Om, b):
    # Returns K(r) and K(r)*r^2 * nbar(r) accoding to our model.
    result1 = np.tanh( (-z+z_pivot)*10 )*K0_bsqr/2. + K0_bsqr/2.
    denom = FlatLambdaCDM(Om0=Om, H0=67.32).comoving_distance(z).to('Mpc').value
    #denom = denom + 1e-5#Add a small number so that the inverse is not INF
    denom[-1] = denom[-2]
    result2 = result1/denom**2./0.1/b/b

    return result1, result2

# Two window functions, very simple Normal dist.
def Windows(z, Wpars):
    std1, avg1, std2, avg2, std3, avg3 = Wpars
    return norm.pdf(z, scale=std1, loc=avg1), norm.pdf(z, scale=std2, loc=avg2), norm.pdf(z, scale=std3, loc=avg3)

# This is the Cl function
def Cl(larray, Karray, Warray, clxcarray, b, Fortran=True):

    if(Fortran):
        a_unique_lst, k_unique_lst, eta_unique_lst, bess_int_keta_reshaped, confH, sigma, \
        deltaCDM, k_length, eta_length = clxcarray

        # We force bGW = bgal = 1
        z = 1/np.array(a_unique_lst) - 1.

        bias = np.ones(eta_length)

        jacobian = confH[0]/a_unique_lst
        W_kernel = Warray*jacobian
        K_kernel = Karray

        Clarray = b*np.ones(len(larray))

        CrossCorr_Cells.c_ells(Clarray, larray, k_unique_lst, eta_unique_lst, bess_int_keta_reshaped,\
                               confH, sigma, deltaCDM, W_kernel, K_kernel, bias, len(larray), k_length, eta_length)

    else:
        a_unique, clxc, j_int, H, sigma, k, N_k, N_l, N_tau = clxcarray


        k = np.tile(k.reshape(N_k, N_tau), (N_l, 1, 1))
        H = np.tile(H.reshape(N_k, N_tau), (N_l, 1, 1))
        clxc = np.tile(clxc.reshape(N_k, N_tau), (N_l, 1, 1))
        sigma = np.tile(sigma.reshape(N_k, N_tau), (N_l, 1, 1))
        k_unique = k[0, :, 0]
        H_unique = H[0, 0, :]

        W_kernel = Warray*H_unique/a_unique
        K_kernel = Karray
        W = np.tile(np.tile(W_kernel, (N_k, 1)), (N_l, 1,1))
        K = np.tile(np.tile(K_kernel, (N_k, 1)), (N_l, 1,1))

        Omega1 = W * (bGW*clxc - 3. * H * sigma / k )* j_int
        Omega2 = K * (bGW*clxc - 3. * H * sigma / k )* j_int

        Clarray = np.trapz(np.sum(Omega1, axis=2)*np.sum(Omega2, axis=2), np.log(k_unique), axis=1)


    return As*4.*np.pi*Clarray

# This is the nearest neighbour in the Om grid, next best thing
def clxc(Om, Fortran=True):
    idx = np.abs(Om-Omarray).argmin()
    if(Fortran):
        return np.load(clxcfolder+'F/'+str(idx)+'.npy')
