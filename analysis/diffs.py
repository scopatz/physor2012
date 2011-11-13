import os

import numpy as np
import tables as tb
from scipy.stats import kendalltau

from nuc_sets import serp_nucs, model_nucs 

NUCLIDES = set(['U234',  'U235',  'U236',  'U238',  'NP237', 'PU238', 'PU239', 
                'PU240', 'PU241', 'PU242', 'AM241', 'AM242', 'AM243', 'CM242', 
                'CM243', 'CM244', 'CM245', 'CM246', 'SE79',  'KR85',  'SR90',  
                'ZR93',  'TC99',  'I129',  'PD107', 'CS134', 'CS135', 'CS137', 
                'SM151', 'EU155', 'H1',    'H3',    'O16',   'HE4',   'NA23',])

# Ensure only serpent nuclides
NUCLIDES = (NUCLIDES & serp_nucs)
NUCLIDES = serp_nucs

BASEPATH = '../data/lwr_base.h5'
PHYSPATH = '../data/lwr_physor2012.h5'


def load_sigma(nucs, lib, p):
    """Loads xs from a Char lib for a set of nucs for the pth row."""
    sig_a_g = {}
    sig_s_g = {}
    sig_s_gh = {}

    f = tb.openFile(lib, 'r')
    for nuc in nucs:
        sig_a_g[nuc] = np.array(f.getNode('/sigma_a/' + nuc)[p])
        sig_s_g[nuc] = np.array(f.getNode('/sigma_s/' + nuc)[p])
        sig_s_gh[nuc] = np.array(f.getNode('/sigma_s_gh/' + nuc)[p])

    f.close()
    return sig_a_g, sig_s_g, sig_s_gh

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


def dtype2template(dt):
    conv = []
    for i, name in enumerate(dt.names):
        t = dt.fields[name][0].type
        if issubclass(t, np.float_):
            conv.append("{" + str(i) + ":.3F}")
        else:
            conv.append("{" + str(i) + "}")
    template = " & ".join(conv) + "\\\\\n"
    return template


def array2tabular(a, path=""):
    end = "\\\\\n"

    header = "\\begin{tabular}{|l|" + "c"*(len(a.dtype)-1) + "|}\n"
    header += "\\hline\n"
    header += " & ".join(a.dtype.names) + end
    header += "\\hline\n"

    body = ""
    dtem = dtype2template(a.dtype)
    for row in a[:44]:
        body += dtem.format(*row) 

    footer = "\\hline\n\\end{tabular}"
    
    s = header + body + footer
    with open(path, 'w') as f:
        f.write(s)


avl_dtype = np.dtype([
    ('nuclide', 'S6'),
    ('$\\epsilon$', float),
    ('$\\tau$', float),
    ('$g$', int),
    ('$h$', int),
    ('$\\sigma_{s,g\\to h,i}$', float),
    ('$P_{g\\to h,i}$', float),
    ('$R_{a/s,g}$', float),
    ])


def analyze_vs_lib(r, libpath, p_start, p_offset, label):
    r0 = 0.369
    r1 = 0.451

    # r_fuel = 0.9 base
    s_a0, s_g0, s_gh0 = load_sigma(NUCLIDES, PHYSPATH, 0 + p_offset)
    gtp0 = group_transfer_prob(s_g0, s_gh0)

    # r_fuel = 1.1 base
    s_a1, s_g1, s_gh1 = load_sigma(NUCLIDES, PHYSPATH, 5 + p_offset)
    gtp1 = group_transfer_prob(s_g1, s_gh1)

    # r_fuel from lib
    s_a, s_g, s_gh = load_sigma(NUCLIDES, libpath, p_start + p_offset)
    gtp = group_transfer_prob(s_g, s_gh)

    # Interpolated values
    is_g, is_gh, igtp = interp_xs(r, r0, r1, s_g0, s_g1, s_gh0, s_gh1)
    fdiff = frac_diff(is_gh, s_gh, s_g)
    sfdiff = sort_frac_diff(fdiff)

    res = []
    for nuc, ind, val in sfdiff:
        g = ind / 19
        h = ind % 19
        if (s_gh[nuc][g,h] != 0.0)  and (is_gh[nuc][g,h] != 0.0) and \
           (s_gh0[nuc][g,h] != 0.0) and (s_gh1[nuc][g,h] != 0.0):
            row = (nuc, val, 
                   kendalltau(is_gh[nuc], s_gh[nuc])[0],
                   g, h, 
                   s_gh[nuc][g,h], 
                   gtp[nuc][g,h], 
                   s_a[nuc][g] / s_g[nuc][g], 
                   )
            res.append(row)
        #else:
        #    print nuc, kendalltau(is_gh[nuc], s_gh[nuc])[0]
    res = np.array(res, dtype=avl_dtype)
    print res
    array2tabular(res, label + '.tex')


def main():
    r25 = 0.3895
    r50 = 0.41
    r75 = 0.4305

    p_offset = 0

    analyze_vs_lib(r25, PHYSPATH, 10, p_offset, '../paper/r25')
    analyze_vs_lib(r50, BASEPATH, 0,  p_offset, '../paper/r50')
    analyze_vs_lib(r75, PHYSPATH, 15, p_offset, '../paper/r75')
    

if __name__ == '__main__':
    main()

