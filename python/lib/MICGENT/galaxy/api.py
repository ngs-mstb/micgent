"""Helper function for interacting with Galaxy through its API
"""
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from future import standard_library

standard_library.install_aliases()
from builtins import str
from builtins import *

from .. import util
from .. import yaml_util
from .. import sig

import bioblend
from bioblend.galaxy import GalaxyInstance
from bioblend.galaxy.jobs import JobsClient
#from bioblend.galaxy.objects import GalaxyInstance
import argh
import sys
import os
from collections import OrderedDict, deque
import warnings

def keep_keys(o,keys):
    return OrderedDict([(key,o[key]) for key in keys if key in o ])

def drop_keys(o,keys):
    keys = set(keys)
    return OrderedDict([(key,val) for key,val in o.items() if key not in keys ])

def add_keyvals(o,kv):
    o = o.copy()
    o.update(kv)
    return o

def select_case(o):
    conditional = o
    value = conditional["test_param"]["value"]
    for case in conditional["cases"]:
        if case["value"] == value:
            return case["inputs"]

def select_dataset_values(o):
    data = o
    options = data["options"]
    ret = []
    value = data["value"]
    if value:
        for value in value["values"]:
            for x in options[value["src"]]:
                if value["id"] == x["id"]:
                    x = x.copy()
                    x["type"] = "data_value"
                    ret.append(x)
    return ret

def api_get_data_values_jobs(x,gi,dep_datasets,verify_cert=False):
    """Obtain job ID that created each value of a data parameter.
    Updates the data_value objects by setting the job_id key.
    The value will be None if no job ID is available.
    Also checks and updates cache of already processed objects to
    minimize the number of API calls."""
    with warnings.catch_warnings():
        if not verify_cert:
            import urllib3
            warnings.filterwarnings("ignore",category=urllib3.exceptions.InsecureRequestWarning) 
        for data_value in x["value"]:
            data_id = data_value["id"]
            if data_id in dep_datasets:
                data_value["job_id"] = dep_datasets[data_id]["job_id"]
            else:
                ## TODO: need special code here for getting job out of ldda - see how this is handled in
                ## the GUI 'i' icon code
                query_res = gi.datasets.show_dataset(data_id,hda_ldda=data_value["src"])
                job_id = None
                if util.is_mapping(query_res):
                    job_id = query_res.get("creating_job",None)
                else:
                    raise ValueError("Galaxy API show_dataset call returned a string instead of a dict. The string was: {}. The input dataset dict was: {}".\
                        format(query_res,data_value))
                data_value["job_id"] = job_id
                dep_datasets[data_id] = data_value.copy()

def api_get_output_properties(x,gi,dep_datasets,verify_cert=False):
    """Obtain properties for the output dataset.
    Updates the input object by setting the keys.
    Also checks cache of already processed objects to
    minimize the number of API calls. Cache is not updated because
    the calling code uses the cache to collect input dataset objects
    that need to be further processed.
    """
    with warnings.catch_warnings():
        if not verify_cert:
            import urllib3
            warnings.filterwarnings("ignore",category=urllib3.exceptions.InsecureRequestWarning) 
        data_id = x["id"]
        if data_id in dep_datasets:
            src = dep_datasets[data_id]
        else:
            src = gi.datasets.show_dataset(data_id,hda_ldda=x["src"])
        x.update(keep_keys(src,("name",)))


def drop_identical_key(o,key1,key2):
    drop = False
    if key1 in o:
        if key2 in o:
            if o[key2] == o[key1]:
                return drop_keys(o,(key2,))
    return o

def convert_value_type(o):
    if "value" in o and "type" in o:
        o = o.copy()
        value = o["value"]
        type_ = o["type"]
        try:
            if type_ == "integer":
                if value == "":
                    o["value"] = None
                else:
                    o["value"] = int(value)
            elif type_ == "float":
                if value == "":
                    o["value"] = None
                else:
                    o["value"] = float(value)
            elif type_ == "boolean":
                if value == "":
                    o["value"] = None
                else:
                    if value == "true":
                        o["value"] = True
                    elif value == "false":
                        o["value"] = False
        except ValueError:
            pass
    return o


def simplify_rerun_json_inputs(x,recursive=False,gi=None,verify_cert=False,dep_datasets=None):
    keys_keep_always = ("name","type","title","label","test_param","value","text_value","inputs")
    keys_keep_dataset_values = keys_keep_always + ("id","src","hid")
    if util.is_sequence(x):
        ret = []
        for y in x: 
            ret.append(simplify_rerun_json_inputs(y,recursive=recursive,gi=gi,verify_cert=verify_cert,dep_datasets=dep_datasets))
    elif util.is_mapping(x):
        x = x.copy()
        type_ = x.get("type","")
        if type_ == "conditional":
            x["inputs"] = select_case(x)
            x = keep_keys(x,keys_keep_always)
        elif type_ == "select":
            x = keep_keys(x,keys_keep_always)
        elif type_ == "data":
            x["value"] = [ keep_keys(_,keys_keep_dataset_values) for _ in select_dataset_values(x) ]
            x = drop_keys(x,("text_value",))
            x = keep_keys(x,keys_keep_always)
            if recursive:
                api_get_data_values_jobs(x,dep_datasets=dep_datasets,gi=gi,verify_cert=verify_cert)
        elif type_ == "data_value":
            ## keys are already finalized
            pass
        else:
            x = keep_keys(x,keys_keep_always)
        x = drop_identical_key(x,"value","text_value")
        x = convert_value_type(x)
        ret = {}
        for key,val in x.items():
            ret[key] = simplify_rerun_json_inputs(val,recursive=recursive,gi=gi,verify_cert=verify_cert,dep_datasets=dep_datasets)
    else:
        ret = x
    return ret


def simplify_rerun_json(job,recursive=False,gi=None,verify_cert=False,dep_datasets=None):
    keys_keep_always = ("panel_section_name","panel_section_id","name","id","version","job_id","history_id")
    job_orig = job
    job = job.copy()    
    inputs = job.get("inputs",None)
    if "inputs" in job:
        del job["inputs"]
    job = keep_keys(job,keys_keep_always)
    job.update(keep_keys(
        job_orig.get("job_status",OrderedDict()),
        ("user_email","create_time","update_time","external_id","job_metrics")))
    if inputs is not None:
        job["inputs"] = simplify_rerun_json_inputs(inputs,recursive=recursive,gi=gi,verify_cert=verify_cert,dep_datasets=dep_datasets)
    job.update(keep_keys(
        job_orig.get("job_status",OrderedDict()),
        ("outputs",)))
    outputs = job.get("outputs",None)
    if outputs is not None:
        for data_val in outputs.values():
            api_get_output_properties(data_val,gi=gi,dep_datasets=dep_datasets,verify_cert=verify_cert)
    return job

def simplify_rerun_json_top(job_id,recursive=False,gi=None,verify_cert=False):
    dep_datasets = {}
    jobs_to_do = deque([job_id])
    jobs_done = {job_id:False}
    job_top = None
    job_simp_top = None
    while True:
        if not len(jobs_to_do):
            break
        job_id = jobs_to_do.popleft()
        job = api_get_job_info(gi=gi,job_id=job_id,verify_cert=verify_cert)
        dep_datasets_len_last = len(dep_datasets)
        job_simp = simplify_rerun_json(job,recursive=recursive,gi=gi,verify_cert=verify_cert,dep_datasets=dep_datasets)
        if not job_top:
            job_top = job
            job_simp_top = job_simp
        else:
            job_simp_top.setdefault("parent_jobs",OrderedDict())[job_id] = job_simp
        jobs_done[job_id] = True
        ## only iterate through dep_datasets if anything was added to them
        if len(dep_datasets) > dep_datasets_len_last:
            for data_id, data_value in dep_datasets.items():
                job_id_dep = data_value.get("job_id",None)
                if job_id_dep:
                    if not job_id_dep in jobs_done:
                        jobs_to_do.append(job_id_dep)
                        jobs_done[job_id_dep] = False

    return job_top,job_simp_top

_job_user_comment_header = """\
Information for re-populating the Galaxy tool form that is needed to re-run a job.
The information is saved in a structured YAML format that is both human-readable
and machine-parsable. The user can consult the information in this file in order
to re-run a particular job through the Galaxy user interface (UI) and reproduce
the previous results.

The YAML file is indented to mirror the hierarchy of data and visual layout of the tool in the UI.
The top-level keys provide information about the tool and the job that executed the tool,
including the user account identified by the email, start and end time of the job, as well
as some job metrics collected from the batch job runner.

The input values that were used by the job are contained in the `inputs` key.

The `type` key describes the type of the UI elements.
The sections in the UI layout correspond to elements with type `section`.
Choices made in the conditional elements are defined in `test_param` objects.

Keys such as `label` provide names of the elements that are presented by the UI.
Their corresponding `name` keys provide the names of the internal parameters. The
exception to that rule is the elements with a type `data_value`. Those describe
datasets, with the `name` key containing the dataset descriptive label from the UI,
and the `hid` key containing the dataset ID within a given Galaxy history as shown
through the UI.

Various other `id` fields are internal Galaxy IDs provided here for provenance information.

This data structure contains all input parameters of the tool, even those that were set to
their default values and then hidden from the UI. To reproduce the run, the UI user needs
to follow the data parameter hierarchy in the `inputs` key from top to bottom and set only the 
parameters that stay visible in the UI after selection of the previous parameters has been made.
"""

_job_user_comment_recursive_header = """\
The key `parent_jobs` contains information about all jobs that generated inputs
for the top-level job in this file, then about all jobs that generated the inputs
of those parent jobs, and so on recursively till the jobs that imported the initial datasets (such
as jobs that ran the `upload1` Get Data tool).

The parent jobs are keyed by their Galaxy job IDs. These job IDs are also referenced within
records describing the input datasets of the children jobs. This allows determining
which parent job has generated which input dataset. The dataset IDs are also listed within
input and output dataset records thus additionally linking parent and child jobs to each other.

To fully reproduce the final outputs from scratch, the user would have to start 
rerunning jobs from the bottom of this file, and work their way up to the job
at the top.
"""

def api_gal_instance(url,key,verify_cert=True):
    ## catch_warnings is not thread-safe
    with warnings.catch_warnings():
        if not verify_cert:
            import urllib3
            warnings.filterwarnings("ignore",category=urllib3.exceptions.InsecureRequestWarning) 
        return GalaxyInstance(url=url, key=key, verify=verify_cert)

def api_gal_instance_base_url(gi):
    """Extract the base URL of the current Galaxy instance object."""
    ##This uses internal field, and so is factored out into a separate method to only
    ##make a single change if the Bioblend implementation changes
    return gi.base_url

def api_gal_instance_verify(gi):
    """Extract the `verify` parameter of the current Galaxy instance object."""
    ##This uses internal field, and so is factored out into a separate method to only
    ##make a single change if the Bioblend implementation changes
    return gi.verify

def api_gal_instance_change_key(gi,key):
    return api_gal_instance(url=api_gal_instance_base_url(gi),key=key,verify_cert=api_gal_instance_verify(gi))

def api_gal_instance_works(gi):
    try:
        version = gi.config.get_version()
        return util.is_mapping(version) and 'version_major' in version
    except bioblend.ConnectionError:
        return False

def api_get_user_info(gi,user_email):
    users = gi.users.get_users(f_email=user_email)
    if not (util.is_sequence(users) and len(users) == 1):
        raise ValueError("Galaxy API get_users query by user email should have returned a list with a single record.")
    return users[0]

def api_get_user_key(gi,user_email):
    user_id = api_get_user_info(gi,user_email=user_email)["id"]
    key = gi.users.get_user_apikey(user_id)
    ## The problem with the API call above is that it returns a hex string on success, and a string "Not Available" if the API
    ## key does not yet exist. To avoid relying on the return error string that might change, we test if the returned value
    ## works as a key by creating new Galaxy instance and making a test call
    user_gi = api_gal_instance_change_key(gi,key)
    if not api_gal_instance_works(user_gi):
        key = gi.users.create_user_apikey(user_id)
    return key

def simplify_rerun_json_yaml(job_id, recursive=False, gi=None, verify_cert=False, out_job_yaml="job.yaml", out_job_user_yaml="job_user.yaml"):
    util.make_pardir(out_job_yaml)
    ## Extract user info from the job and create another Galaxy instance with the user API key because
    ## the admin cannot access information about user datasets. 
    ## TODO: We would still need the original admin-keys Galaxy instance to get those metrics from user jobs
    ## that might be made unavailable to regular users due to the Galaxy config settings (resource consumptions,
    ## for example). Right now we simply use user_gi instead of the original gi. Long term, we need
    ## to pass both down the chain of functions if the job metrics data is requested, and use the original,
    ## presumably admin-level gi when the job API is queried, and use user_gi for dataset queries.
    job_info = api_get_job_info(gi=gi,job_id=job_id,verify_cert=verify_cert)
    user_email=job_info["job_status"]["user_email"]    
    user_key = api_get_user_key(gi,user_email=user_email)
    user_gi = api_gal_instance_change_key(gi,user_key)

    job, job_user = simplify_rerun_json_top(job_id, recursive=recursive, gi=user_gi, verify_cert=verify_cert)
    job_user["user_email"] = user_email

    job_user = yaml_util.to_yaml(job_user)
    yaml_util.dump_yaml(job,out_job_yaml)    
    comment = _job_user_comment_header
    if recursive:
        comment += "\n"+_job_user_comment_recursive_header
    job_user.yaml_set_start_comment(comment)        
    util.make_pardir(out_job_user_yaml)
    yaml_util.dump_yaml(job_user,out_job_user_yaml)

def api_get_job_info(gi, job_id, verify_cert=False):
    """Use Galaxy API calls to collect information about a given job"""
    ## catch_warnings is not thread-safe
    with warnings.catch_warnings():
        if not verify_cert:
            import urllib3
            warnings.filterwarnings("ignore",category=urllib3.exceptions.InsecureRequestWarning) 
        jc = JobsClient(gi)
        job_info_main = jc.show_job(job_id, full_details=True)        
        url = jc.gi._make_url(jc)
        try:
            job_info_rerun = jc._get('/'.join([job_id,"build_for_rerun"]))
        except bioblend.ConnectionError:
            ## we will get this for the upload tool job with a traceback string payload of ConfigDoesNotAllowException server-side exception
            job_info_rerun = OrderedDict([("job_id",job_id),("id",job_info_main["tool_id"])])
        job_info_rerun["job_status"] = job_info_main        
        return job_info_rerun

def job_info_export(job_id=None, job_yaml=None, url=None, api_key=None, verify_cert=False, 
    recursive=False,
    out_job_yaml="job.yaml", out_job_user_yaml="job_user.yaml", 
    check_sig=False,make_sig=False,key=None,msg_sig="Job YAML signature mismatch"):
    """Extract information about Galaxy job sufficient to re-run it, and save as a YAML file.

    Instead of a job ID, a previously saved YAML file can be provided, in order to use the job ID
    recorded in this file, and create a new file with the refreshed information.

    :param check_sig: If reading job_yaml, verify signature from file job_yaml+".sig"
    :param make_sig: Create signature in a file out_job_yaml+".sig"
    """
    assert bool(job_id) != bool(job_yaml), "Either one of job_id or the previously created job YAML file is required"

    if job_yaml:
        if check_sig:
            sig.file_sig_cmp_msg(job_yaml+".sig",job_yaml,key=key,msg=msg_sig)
            ## nothing is copied from job_user.yaml, so no need to verify that signature
        job = yaml_util.load_yaml(job_yaml)
        job_id = job["job_id"]

    gi = api_gal_instance(url=url,key=api_key,verify_cert=verify_cert)
    simplify_rerun_json_yaml(job_id,recursive=recursive,gi=gi,verify_cert=verify_cert,
        out_job_yaml=out_job_yaml,out_job_user_yaml=out_job_user_yaml)
    if make_sig:
        sig.file_sig(out_job_yaml,key=key,out_file=out_job_yaml+".sig")
        sig.file_sig(out_job_user_yaml,key=key,out_file=out_job_user_yaml+".sig")

## import package module and add argh entry points

def _main():
    from .. import arg_parsing
    parser = arg_parsing.ArghParserChainedConfig()
    parser.add_commands([
        job_info_export
    ])
    parser.dispatch()


if __name__ == "__main__":
    _main()
