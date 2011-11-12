import os

import numpy as np
import tables as tb

from nuc_sets import serp_nucs, model_nucs 

nuclides = set(['U234',  'U235',  'U236',  'U238',  'NP237', 'PU238', 'PU239', 
                'PU240', 'PU241', 'PU242', 'AM241', 'AM242', 'AM243', 'CM242', 
                'CM243', 'CM244', 'CM245', 'CM246', 'SE79',  'KR85',  'SR90',  
                'ZR93',  'TC99',  'I129',  'PD107', 'CS134', 'CS135', 'CS137', 
                'SM151', 'EU155', 'H1',    'H3',    'O16',   'HE4',   'NA23',])

# Ensure only serpent nuclides
nuclides = (nuclides & serp_nucs)
nuclides = serp_nucs


def load_sigma_s(nucs, lib, p):
    """Loads scattering xs from a Char lib for a set of nucs for the pth row."""
    sig_s_g = {}
    sig_s_gh = {}

    f = tb.openFile(lib, 'r')
    for nuc in nucs:
        sig_s_g[nuc] = np.array(f.getNode('/sigma_s/' + nuc)[p])
        sig_s_gh[nuc] = np.array(f.getNode('/sigma_s_gh/' + nuc)[p])

    f.close()
    return sig_s_g, sig_s_gh

def group_transfer_prob(sigma_s_g, sigma_s_gh):
    """Computes the group transfer probabilities from scattering dictionaries."""
    gtp = {}
    for nuc in sigma_s_g:
        gtp[nuc] = sigma_s_gh[nuc] / sigma_s_g[nuc]
    return gtp 


def interp_data(r, r0, r1, data0, data1):
    idata = {}
    ir_dr = (r - r0)/(r1 - r0)
    for nuc in data0:
        ddata = data1[nuc] - data0[nuc]
        idata[nuc] = ddata*ir_dr + data0[nuc]
    return idata


def interp_xs(r, r0, r1, sig_g0, sig_g1, sig_gh0, sig_gh1):
    isig_g = interp_data(r, r0, r1, sig_g0, sig_g1)
    isig_gh = interp_data(r, r0, r1, sig_gh0, sig_gh1)
    igtp = group_transfer_prob(isig_g, isig_gh)
    return isig_g, isig_gh, igtp


def frac_diff(interp, meas, norm):
    fdiff = {}
    for nuc in interp:
        fdiff[nuc] = (interp[nuc] - meas[nuc])/norm[nuc]
    return fdiff

def sort_frac_diff(fdiff):
    key_func = lambda x: np.max(np.abs(x[1]))
    s = sorted(fdiff.items(), key=key_func, reverse=True)
    s = [(k, np.argmax(np.abs(v))) for k,v in s]
    s = [(k, i, fdiff[k].flat[i]) for k,i in s]
    return s



def main():
    r0 = 0.369
    r025 = 0.3895
    r05 = 0.41
    r075 = 0.4305
    r1 = 0.451

    p_offset = 0

    # r_fuel = 0.9 base
    s_g0, s_gh0 = load_sigma_s(nuclides, '/home/scopatz/lwr_physor2012.h5', 0 + p_offset)
    gtp0 = group_transfer_prob(s_g0, s_gh0)

    # r_fuel = 1.1 base
    s_g1, s_gh1 = load_sigma_s(nuclides, '/home/scopatz/lwr_physor2012.h5', 5 + p_offset)
    gtp1 = group_transfer_prob(s_g1, s_gh1)

    # r_fuel = 0.95 base
    s_g025, s_gh025 = load_sigma_s(nuclides, '/home/scopatz/lwr_physor2012.h5', 10 + p_offset)
    gtp025 = group_transfer_prob(s_g025, s_gh025)

    # Interpolate to 0.95 base
    is_g025, is_gh025, igtp025 = interp_xs(r025, r0, r1, s_g0, s_g1, s_gh0, s_gh1)
    fdiff025 = frac_diff(is_gh025, s_gh025, s_g025)
    sfdiff025 = sort_frac_diff(fdiff025)
    #for nuc, ind, val in sfdiff025:
    #    if (s_gh0[nuc].flat[ind] != 0.0)   and (s_gh1[nuc].flat[ind] != 0.0) and \
    #       (s_gh025[nuc].flat[ind] != 0.0) and (is_gh025[nuc].flat[ind] != 0.0):
    #        print nuc, val, (ind/19, ind%19), s_gh025[nuc].flat[ind]


    # r_fuel = 1.05 base
    s_g075, s_gh075 = load_sigma_s(nuclides, '/home/scopatz/lwr_physor2012.h5', 15 + p_offset)
    gtp075 = group_transfer_prob(s_g075, s_gh075)

    # Interpolate to 1.05 base
    is_g075, is_gh075, igtp075 = interp_xs(r075, r0, r1, s_g0, s_g1, s_gh0, s_gh1)
    fdiff075 = frac_diff(is_gh075, s_gh075, s_g075)
    sfdiff075 = sort_frac_diff(fdiff075)
    for nuc, ind, val in sfdiff075:
        if (s_gh0[nuc].flat[ind] != 0.0)   and (s_gh1[nuc].flat[ind] != 0.0) and \
           (s_gh075[nuc].flat[ind] != 0.0) and (is_gh075[nuc].flat[ind] != 0.0):
            print nuc, val, (ind/19, ind%19), s_gh075[nuc].flat[ind]


if __name__ == '__main__':
    main()

