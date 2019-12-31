### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MDDV package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Configuration and resources of this package"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *
import pkg_resources

_pkg_data = "data" ## if paths are used, they should be unix style - passed to Requirement
_pkg_test_data = "test_data"

## see https://pythonhosted.org/setuptools/pkg_resources.html#basic-resource-access
## for a description of what can be passed as `pkg` spec.
## The extraction of default package name will probably break
## if the top level package is a `namespace` package as per docs above
def get_pkg_data_dir(pkg,pkg_data=_pkg_data):
    return pkg_resources.resource_filename(pkg,pkg_data)

def get_pkg_test_data_dir(pkg,pkg_test_data=_pkg_test_data):
    return pkg_resources.resource_filename(pkg,pkg_test_data)

def get_pkg_data_file(pkg,name,pkg_data=_pkg_data):
    return pkg_resources.resource_filename(pkg,pkg_data+"/"+name)

def get_pkg_data_string(pkg,name,pkg_data=_pkg_data,encoding="utf-8"):
    raw_content = pkg_resources.resource_string(pkg,pkg_data+"/"+name)
    return raw_content.decode(encoding=encoding)

