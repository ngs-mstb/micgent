### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""UUID generation functions"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import *
import random
import string
import uuid

__all__ = ["maxIdLen","idDtype","genId","zeroId","random_string"]

# maximum length of UUID
maxIdLen = 32

# Numpy dtype for UUID
idDtype = 'S%i' % maxIdLen

def genId(n=None):
    """Generate standard conforming UUIDs as a 32 byte printable hex string.
    @param n if not None, must be a number and a Numpy 'S32' array of length n will be returned - otherwise, a single scalar.
    @todo This uses Python hex.uuid4() method that generates random numbers.
    Takes 1 min to generate 10^6 UUIDs on a 3GHz machine.
    In my installation, the alternative uuid.uuid1() apparently uses Linux libuuid through
    ctypes and segfaults every now and then.
    Some implementations read /proc/sys/kernel/random/uuid where available,
    which might be safer than using ctypes"""
    if n is None:
        return uuid.uuid4().hex
    else:
        import numpy
        a = numpy.empty(n,dtype=idDtype)
        for i in range(n):
            a[i] = uuid.uuid4().hex
        return a

def zeroId(n,val=None):
    import numpy
    a = numpy.zeros(n,dtype=idDtype)
    if val is not None:
        a[:] = val
    return a

_random_symbols_choice = string.ascii_lowercase + string.digits

def random_string(length=6, symbols=_random_symbols_choice):
    return ''.join(random.choice(symbols) for x in range(length))

