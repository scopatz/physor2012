import numpy as np
import tables as tb

from nuc_sets import serp_nucs, model_nucs 

nuclides = set(['U234',  'U235',  'U236',  'U238',  'NP237', 'PU238', 'PU239', 'PU240', 'PU241', 'PU242', 
                'AM241', 'AM242', 'AM243', 'CM242', 'CM243', 'CM244', 'CM245', 'CM246', 'SE79',  'KR85',  
                'SR90',  'ZR93',  'TC99',  'I129',  'PD107', 'CS134', 'CS135', 'CS137', 'SM151', 'EU155'])

nuclides = (nuclides & serp_nucs)

