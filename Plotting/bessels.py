import numpy as np
from scipy import special
import scipy.integrate as integrate
from scipy import interpolate

def bessel_calculation(bessel_tabulate, ell_lst, x, x_Recompute, ell_Recompute, ktau_lst, ktau_Recompute_index):

    bess_int_lst = np.zeros((len(ell_lst), len(x)))
    bess_int_ktau_lst = np.zeros((len(ell_lst), len(ktau_lst)))

    bessPrime_int_lst = np.zeros((len(ell_lst), len(x)))
    bessPrime_int_ktau_lst = np.zeros((len(ell_lst), len(ktau_lst)))

    dx = x[2] - x[1]
    dx_Recompute = x_Recompute[2] - x_Recompute[1]
    #print len(ell_lst)


    for ell_index in xrange(len(ell_lst)):

        #print "Performing Bessel Computation for l= ", ell_lst[ell_index]
        bess_int_lst[ell_index, :] = np.cumsum(special.spherical_jn(ell_lst[ell_index], x, derivative = False))*dx
        bessPrime_int_lst[ell_index, :] = np.cumsum(special.spherical_jn(ell_lst[ell_index], x, derivative = True))*dx

        interpolation = interpolate.interp1d(x, bess_int_lst[ell_index, :])
        bess_int_ktau_lst[ell_index, :] = interpolation(ktau_lst)

        interpolationDer = interpolate.interp1d(x, bessPrime_int_lst[ell_index, :])
        bessPrime_int_ktau_lst[ell_index, :] = interpolationDer(ktau_lst)


    if bessel_tabulate == 1:

        for ell_index in xrange(ell_Recompute):

            #print "Performing Bessel ReComputation for l= ", ell_lst[ell_index]

            bess_int_lst[ell_index, :] = np.cumsum(special.spherical_jn(ell_lst[ell_index],\
                                                   x_Recompute, derivative = False))*dx_Recompute
            bessPrime_int_lst[ell_index, :] = np.cumsum(special.spherical_jn(ell_lst[ell_index],\
                                                        x_Recompute, derivative = True))*dx_Recompute


            for i in xrange(ktau_Recompute_index):
                interpolation = interpolate.interp1d(x_Recompute, bess_int_lst[ell_index, :])
                bess_int_ktau_lst[ell_index, i] = interpolation(ktau_lst[i])

                interpolationDer = interpolate.interp1d(x_Recompute, bessPrime_int_lst[ell_index, :])
                bessPrime_int_ktau_lst[ell_index, i] = interpolationDer(ktau_lst[i])

    #np.save('Files/bessel_integrals.npy', bess_int_ktau_lst)
    #np.save('Files/besselPrime_integrals.npy', bessPrime_int_ktau_lst)

    return bess_int_ktau_lst, bessPrime_int_ktau_lst
