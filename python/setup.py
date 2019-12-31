distr_name = "micgent"
script_pref = "mcg"
version = "1.0.1"
pkg_subdir = "lib"

import os
from fnmatch import fnmatch

from setuptools import setup, find_packages
import logging

def logging_config(detail="high",level=logging.DEBUG):
    """Some common logging configuration scenarious"""
    if detail == "high":
        logging.basicConfig(level=level, datefmt="%y-%m-%d %H:%M:%S",
            filemode="a",
            format=("%(asctime)s [%(levelname)5.5s] pid:%(process)-5s\t"
                + "thread:%(threadName)10.10s\t"
                + "%(module)20.20s.:%(funcName)-12.12s:%(lineno)-5s:\t"
                + "%(message)s"))
    else:
        raise ValueError("Unknown value detail={}".format(detail))

logging_config()
log = logging.getLogger("MAIN")

def _entry_point(script_name,pkg_path):
    return script_pref+script_name+' = '+pkg_path

def _entry_point_argh(pkg_path,use_pref=True):
    script_name = pkg_path.rsplit(".",1)[-1].replace("_","-")
    s = script_name+' = '+pkg_path+":_main"
    if use_pref:
        s = script_pref + s
    print (s)
    return s

def iter_files_tree(dir,to_base=True,patt=None):
    dir_base = os.path.dirname(dir)
    if patt is not None:
        if isinstance(patt,str):
            patt = [ patt ]
    for root,dirs,files in os.walk(dir):
        for f in files:
            if patt:
                match = False
                for p in patt:
                    if fnmatch(f,p):
                        match = True
                        break
                if not match:
                    continue
            path = os.path.join(root,f)
            if to_base:
                path = os.path.relpath(path,dir_base)
            yield path

def _list_pkg_files(pkg_name,dir,patt=None):
    return list(iter_files_tree(os.path.join(pkg_subdir,pkg_name,dir),patt=patt))


packages = find_packages(pkg_subdir)

log.info("packages={}".format(packages))

package_data = dict((
    (
    pkg_name,
    _list_pkg_files(pkg_name,"data") +
    _list_pkg_files(pkg_name,"test",patt="*.py") +
    _list_pkg_files(pkg_name,"test_data")
    ) \
    for pkg_name in packages
))

log.info("package_data={}".format(package_data))

setup(
    name = distr_name,
    version = version,
    packages = packages,
    package_dir = {'':pkg_subdir},
    #argh is used to auto-generate command line argument 
    #processing in entry points
    install_requires = [
        'future',
        'six',
        'argh',
        'argcomplete',
        'pytest',
        'biopython',
        'pybedtools',
        'intervals',
        'cutadapt',
        'intervaltree',
        'petl',
        'petlx',
        'pandas',
        'numpy',
        'ruamel.yaml',
        'jinja2',
        'beautifulsoup4',
        'lxml',
        'bioblend'
        ],
    extras_require = dict(denovo_mic_assembly=['jsonnet']),
	#this will install pytest module
    tests_require=['pytest'],
    package_data = package_data,
    entry_points = {
        #'console_scripts' is a fixed group name - it will cause
        #creation of scripts
        'console_scripts': [
            _entry_point_argh('MICGENT.plasmid_finder'),
            ##module that is just script_pref[:-1] must exist;
            ##it will be exposed as a script w/o prefix
            #_entry_point_argh(script_pref[:-1],use_pref=False)
            ]
        }
)

