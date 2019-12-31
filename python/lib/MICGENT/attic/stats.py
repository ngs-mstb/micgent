from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *
from past.utils import old_div
import numpy as np
import scipy as sp

def betabinom(k,n,p,ro):
    """Beta-binomial distribution parameterized by mean probability and cluster correlation.
    ro - cluster correlation has a range (0,1); the larger values mean higher overdispersion
    relative to binomial distribution, which is achieved when ro -> 0.
    The implementation uses vectorized SciPy and Numpy functions, so inputs can be a
    combination of scalars and Numpy 1D arrays.
    """
    theta = old_div(1,ro**2) - 1
    alpha = p*theta
    beta = (1-p)*theta
    part_1 = sp.special.comb(n,k)
    part_2 = sp.special.betaln(k+alpha,n-k+beta)
    part_3 = sp.special.betaln(alpha,beta)
    result = (np.log(part_1) + part_2)- part_3
    return np.exp(result)

