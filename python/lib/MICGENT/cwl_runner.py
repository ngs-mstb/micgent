#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from future import standard_library

standard_library.install_aliases()
from builtins import *
from . import util
from . import resources
from . import workflow_util
from . import yaml_util
from . import jinja_util
from . import pysteps

import os
import shutil
import subprocess
import glob
import shlex

from MGT.Logging import *


def add_dep_resolver(dep_resolver_fn,
                     dep_type="galaxy_packages",
                     mode="a",
                     pkg_name=None,
                     pkg_version="1.0",
                     conda_root=None,
                     conda_env=None,
                     conda_activate="conda activate",
                     keep_all_pkg_versions=False):
    """
    Create or update Cwltool (or Galaxy) dependency-resolvers-configuration file and/or associated env.sh file.

    The main objective of this facility is to use entire pre-existing Conda environments in CWL requirements section.
    Example of CWL:
    ```
    hints:
      SoftwareRequirement:
        packages:
        - package: 'ngs-mstb'
          version:
          - '1.0'
    ```
    Example of corresponding call here (assuming you have already created the Conda environments and are running this from under
    your target Conda installation (in any environment):
    ```
    add_dep_resolver("my_dep_resolver.yaml",conda_env="ngs-mstb",mode="w")
    add_dep_resolver("my_dep_resolver.yaml",conda_env="ngs-mstb-py2")
    ```
    or, from with the CLI (imagining that you need two Conda environments in your workflow):
    ```
    python -m MICGENT.cwl_runner add-dep-resolver --conda-env ngs-mstb --mode w my_dep_resolver.yaml # `--mode w` resets dep_resolver config
    python -m MICGENT.cwl_runner add-dep-resolver --conda-env ngs-mstb-py2 my_dep_resolver.yaml # default `--mode` is `a`, to append
    ```
    Then, run `cwl-runner --beta-dependency-resolvers-configuration my_dep_resolver.yaml ...`.


    :param dep_resolver_fn: Path to dep-resolver yaml file
    :param dep_type: type of resolver record to add (currently only galaxy_packages
    :param mode: ["w","a"] - edit mode for dep-resolver file. "a" means appending to existing content
    :param pkg_name: name of dependency package to add; if None, use conda_env if that one is defined; if neither is defined,
    only edit dep-resolver file and, possibly, create root directrory for package dependency files
    :param pkg_version: version of dep-resolver package to use as label
    :param conda_root: Path to Conda install (see conda_activate_cmd)
    :param conda_env: Name of Conda environment to activate when pkg_name is requested (see conda_activate_cmd)
    :param conda_activate: Command to activate Conda environment with (see conda_activate_cmd)
    :param keep_all_pkg_versions: if False, erase any prior existing versions of pkg_name if such exist
    """
    assert dep_type == "galaxy_packages", "Unimplemented dependecy type: {}".format(dep_type)
    dep_resolver_fn = os.path.abspath(dep_resolver_fn)
    dr = []
    if mode == "a" and os.path.exists(dep_resolver_fn):
        dr = yaml_util.load_yaml(dep_resolver_fn)
    deps_root = os.path.abspath(os.path.join(os.path.dirname(dep_resolver_fn), "deps"))
    dr_el_ind = util.list_find(dr, lambda x: x["type"] == dep_type)
    dr_updated = False
    if dr_el_ind < 0:
        dr_el = {"type": dep_type,
                 "base_path": deps_root}
        dr.append(dr_el)
        dr_updated = True
    else:
        dr_el = dr[dr_el_ind]
        if not dr_el["base_path"] == deps_root:
            log.info("Record for dependency type {} already exists in {}, but with a path {} different from requested path {}. Replacing.". \
                format(dep_type, dep_resolver_fn, dr_el["base_path"], deps_root))
            dr_el["base_path"] = deps_root
            dr_updated = True
    if os.path.exists(deps_root):
        if not mode == "a":
            shutil.rmtree(deps_root, ignore_errors=True)
    os.makedirs(deps_root, exist_ok=True)
    if pkg_name is None:
        if conda_env is not None:
            pkg_name = conda_env
    if pkg_name:
        pkg_dir = os.path.join(deps_root, str(pkg_name))
        if not os.path.exists(pkg_dir):
            os.makedirs(pkg_dir, exist_ok=True)
        else:
            if not keep_all_pkg_versions:
                shutil.rmtree(pkg_dir, ignore_errors=True)
                os.makedirs(pkg_dir, exist_ok=True)
        env_dir = os.path.join(pkg_dir, str(pkg_version))
        os.makedirs(env_dir, exist_ok=True)
        cmd_env_activate = pysteps.conda_activate_cmd(conda_root=conda_root, conda_env=conda_env,
                                                      conda_activate=conda_activate)
        with open(os.path.join(env_dir, "env.sh"), "w") as env_sh:
            env_sh.write(cmd_env_activate + "\n")
    if dr_updated:
        yaml_util.dump_yaml(dr, dep_resolver_fn)


def run_toil(wf,
             wf_inputs,
             dep_resolver=None,
             outdir=None,
             workdir=None,
             restart=False,
             leave_workdir=False,
             clean_jobstore=False,
             runner_use_conda=False,
             runner_conda_root=None,
             runner_conda_env="toil",
             runner_conda_activate="conda activate",
             dry_run=False,
             batchSystem=None,
             disableCaching=False,
             logLevel="INFO",
             maxLocalJobs=4,
             maxCores=None,
             ## try to plow through select() exceptions in large workflows
             ## that are probably related to file descriptor limits:
             retryCount=5,
             max_wf_tries=3,
             max_job_time=60 * 60 * 12,  # seconds
             use_container=False,
             toil_args=None,
             runner_log_dir=None
             ):
    no_container = not use_container
    curdir = os.getcwd()

    assert wf and os.path.isfile(wf), "Workflow file '{}' must exist".format(wf)
    wf = os.path.abspath(wf)

    assert wf_inputs and os.path.isfile(wf_inputs), "Workflow inputs file '{}' must exist".format(wf_inputs)
    wf_inputs = os.path.abspath(wf_inputs)

    if dep_resolver:
        dep_resolver = os.path.abspath(dep_resolver)
        assert os.path.isfile(dep_resolver), "Dependency resolver file '{}' must exist".format(dep_resolver)

    ## cwltool/draft2tool.py needs to be patched to avoid undefined dockerimg var access when cachdir is defined
    # subprocess.check_call(["cwltool","--tmpdir-prefix",curdir,"--no-container","--cachedir",cachedir, os.path.join(cwl_dir,"gene_extractor.cwl"), wf_inputs])
    jobStore = os.path.join(curdir, "jobStore")

    # keep it short, or we get "path too long" errors in tasks
    if not workdir:
        workdir = os.path.join(curdir, "wd")
    else:
        workdir = os.path.abspath(workdir)
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    if not outdir:
        outdir = os.path.join(curdir, "out")
    else:
        outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not runner_log_dir:
      runner_log_dir = os.path.join(curdir,"runner_logs")
    else:
      runner_log_dir = os.path.abspath(runner_log_dir)
    if not os.path.exists(runner_log_dir):
        os.makedirs(runner_log_dir)

    env = {}
    env["TMPDIR"] = workdir  # trying to prevent torque-wrapper.sh in /tmp
    env["TOIL_WORKDIR"] = workdir
    env["TMP"] = workdir
    env["TEMP"] = workdir

    run_args = None
    try:
      for wf_try in range(1, max_wf_tries + 1):
          try:
              # Tried to use "--maxCores", "8", to keep local singleMachine used by cwltoil under control,
              # but that also throttles Torque backend used by actual tasks
              cmd = ["cwltoil"]
              if batchSystem:
                  cmd += ["--batchSystem", batchSystem]
              if no_container:
                  cmd += ["--no-container"]
              ## cleaning fo jobStore is very slow - uses shutil.rmtree
              ## we always skip cleaning here, and do it with our own
              ## method.
              ## cleaning of the job work dirs is probably distributed,
              ## and has to be done anyway to keep disk usage in check
              ## during run
              cmd += ["--jobStore", jobStore,
                     "--workDir", workdir,
                     "--outdir", outdir,
                     "--cleanWorkDir", "always" if not leave_workdir else "never",
                     "--clean", "never"]
              if dep_resolver:
                  cmd += ["--beta-dependency-resolvers-configuration", dep_resolver]
              cmd += ["--tmpdir-prefix", workdir,
                      "--tmp-outdir-prefix", workdir,
                      "--retryCount", str(retryCount),
                      "--maxJobDuration", str(max_job_time),
                      "--stats", # too many open files error
                      "--clusterStats",os.path.join(runner_log_dir,"toil.clusterstats.try_{}.json".format(wf_try)),
                      "--writeLogs", runner_log_dir, #worker logs (only failed ones unless logDebug)
                      "--maxLogFileSize", "1M",
                      "--logLevel", logLevel,
                      "--logFile", os.path.join(runner_log_dir,"toil.try_{}.log".format(wf_try))]
              if maxCores:
                  cmd += ["--maxCores", str(maxCores)]
              if maxLocalJobs:
                  cmd += ["--maxLocalJobs", str(maxLocalJobs)]                  
              if restart:
                  cmd.append("--restart")
              if disableCaching:
                  cmd.append("--disableCaching")
              if toil_args:
                  cmd += shlex.split(toil_args)
              cmd += [wf, wf_inputs]
              log.info("Toil command is: " + " ".join(cmd))
              ## changing cwd seems to be the only way to get all cwltoil temp files
              ## isolated, but it assumes that all paths in input yaml are absolute
              run_args = pysteps.run_step(cmd=cmd,
                               step_dir=workdir,
                               descr="Toil run attempt {}".format(wf_try),
                               env=env,
                               use_conda=runner_use_conda,
                               conda_root=runner_conda_root,
                               conda_env=runner_conda_env,
                               conda_activate=runner_conda_activate,
                               dry_run=dry_run)
          except subprocess.CalledProcessError:
              if wf_try < max_wf_tries:
                  log.exception(
                      "Toil has exited with an error, will try running it {} more time(s).".format(max_wf_tries - wf_try))
                  restart = True
              else:
                  raise
          else:          
              break
    finally:
      try:
        pysteps.run_step(cmd=["toil","stats","--raw",
           "--outputFile",os.path.join(runner_log_dir,"toil.stats.json"),
           jobStore],
           step_dir=workdir,
           descr="Toil stats collection from jobStore",
           env=env,
           use_conda=runner_use_conda,
           conda_root=runner_conda_root,
           conda_env=runner_conda_env,
           conda_activate=runner_conda_activate,
           dry_run=dry_run)
      except:
        log.exception("Toil stats collection failed")
      if not leave_workdir:
          util.rmtree(workdir, ignore_errors=True, threads=8, thread_targets="*")
          ## cwltoil still puts some junk into outdir. Apparently, it sets runtime.tmpdir to outdir, hence
          ## cleanup extra here:
          util.rmtree([os.path.join(outdir, patt) for patt in ["out_tmpdir*", "tmp*"]], 
              ignore_errors=True, threads=8, thread_targets="", file_ok=True)
      if clean_jobstore:
          util.rmtree(jobStore, ignore_errors=True, threads=8, thread_targets="*/*")
    return run_args

## import package module and add argh entry points

def _main():
    from . import arg_parsing
    parser = arg_parsing.ArghParserChainedConfig()
    parser.add_commands([
        add_dep_resolver,
        run_toil
    ])
    parser.dispatch()

if __name__ == "__main__":
    _main()
