import os

import numpy as np
import tables as tb

from nuc_sets import serp_nucs, model_nucs 

nuclides = set(['U234',  'U235',  'U236',  'U238',  'NP237', 'PU238', 'PU239', 
                'PU240', 'PU241', 'PU242', 'AM241', 'AM242', 'AM243', 'CM242', 
                'CM243', 'CM244', 'CM245', 'CM246', 'SE79',  'KR85',  'SR90',  
                'ZR93',  'TC99',  'I129',  'PD107', 'CS134', 'CS135', 'CS137', 
                'SM151', 'EU155', 'H1'])

# Ensure only serpent nuclides
nuclides = (nuclides & serp_nucs)


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

def main():
    r0 = 0.369
    r025 = 0.3895
    r075 = 0.4305
    r1 = 0.451

    # r_fuel = 0.9 base
    s_g0, s_gh0 = load_sigma_s(nuclides, '/home/scopatz/lwr_physor2012.h5', 0)
    gtp0 = group_transfer_prob(s_g0, s_gh0)

    # r_fuel = 1.1 base
    s_g1, s_gh1 = load_sigma_s(nuclides, '/home/scopatz/lwr_physor2012.h5', 5)
    gtp1 = group_transfer_prob(s_g1, s_gh1)

    # r_fuel = 0.95 base
    s_g025, s_gh025 = load_sigma_s(nuclides, '/home/scopatz/lwr_physor2012.h5', 10)
    gtp025 = group_transfer_prob(s_g025, s_gh025)

    is_g025 = interp_data(r025, r0, r1, s_g0, s_g1)
    is_gh025 = interp_data(r025, r0, r1, s_gh0, s_gh1)
    igtp025 = group_transfer_prob(is_g025, is_gh025)
    #igtp025 = interp_data(r025, r0, r1, gtp0, gtp1)
    print is_g025['H1']
    #print igtp025['H1']
    #print igtp025['H1'] - gtp025['H1']
    print np.max(np.abs((is_gh025['H1'] - s_gh025['H1'])/s_g025['H1']))
    print np.max(np.abs(igtp025['H1'] - gtp025['H1']))


if __name__ == '__main__':
    main()
