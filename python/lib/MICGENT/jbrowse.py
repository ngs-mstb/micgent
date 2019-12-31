from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division


from future import standard_library
standard_library.install_aliases()
from builtins import *
from builtins import object
from subprocess import check_call
import logging
from . import config
from . import util
import os, glob, shutil, tempfile

log = logging.getLogger(__name__)

class jbrowse(object):


    def __init__(self,
            jbrowse_url,
            jbrowse_bin_dir,
            jbrowse_index_html_tpl=None,
            jbrowse_data_subdir="data"):
        #the line below should be the first in
        #order to get all parameters into a dict
        self.opt = locals().copy()
        self.opt.pop("self")
        self.opt["jbrowse_bin_dir"] = util.abspath(self.opt["jbrowse_bin_dir"])

    def gff_to_jbrowse(self,
            gff_file,
            index_html,
            data_dir_out,
            jbrowse_url=None
            ):
        opt = self.opt
        if jbrowse_url is None:
            jbrowse_url = opt["jbrowse_url"]
        env = None
        if self.opt.get("jbrowse_bin_dir",None):
            env = os.environ.copy()
            util.add_to_path(opt["jbrowse_bin_dir"],
                    prepend=True,
                    env=env)
        if not os.path.exists(data_dir_out):
            os.makedirs(data_dir_out)
        gff_file = util.abspath(gff_file)
        #fasta_file = util.abspath(fasta_file)
        jbrowse_out_dir = os.path.join(data_dir_out,opt["jbrowse_data_subdir"])
        #can use os.devnull to discard all output
        jbrowse_conv_log_base = os.path.join(os.getcwd(),"jbrowse_conv_log")
        with open(jbrowse_conv_log_base+".out","w") as stdout,\
            open(jbrowse_conv_log_base+".err","w") as stderr:

            check_call(["prepare-refseqs.pl","--gff",gff_file,"--out",jbrowse_out_dir],
                    env=env,
                    stdout=stdout,
                    stderr=stderr)
            #@todo use biodb-to-json instead with flat file input, and accept config
            #file as a parameter (provide a default one too). See volvox.json config
            #in the distribution. Also add dropped_features param to load everything
            #unique in field 3 of GFF and check that only dropped_features are missing
            #from the config
            check_call(["flatfile-to-json.pl","--gff",gff_file,"--out",jbrowse_out_dir,
                "--trackLabel","Genes",
                "--cssClass","feature5",
                "--type","gene",
                "--autocomplete","all"
                "--getLabel",
                "--getType"
                ],
                env=env,
                stdout=stdout,
                stderr=stderr)

            check_call(["flatfile-to-json.pl","--gff",gff_file,"--out",jbrowse_out_dir,
                "--trackLabel","CDS",
                "--cssClass","generic_parent",
                "--subfeatureClasses",'{ "exon" : "exon" }',
                "--type","CDS",
                "--type","exon",
                "--autocomplete","all"
                "--getLabel",
                "--getType",
                "--getSubs",
                "--getPhase"
                ],
                env=env,
                stdout=stdout,
                stderr=stderr)

            check_call(["flatfile-to-json.pl","--gff",gff_file,"--out",jbrowse_out_dir,
                "--trackLabel","Peptides",
                "--cssClass","est",
                "--subfeatureClasses",'{ "mat_peptide" : "transcript-CDS" }',
                "--type","mat_peptide",
                "--autocomplete","all"
                "--getLabel",
                "--getType",
                "--getSubs",
                "--getPhase"
                ],
                env=env,
                stdout=stdout,
                stderr=stderr)

            check_call(["flatfile-to-json.pl","--gff",gff_file,"--out",jbrowse_out_dir,
                "--trackLabel","Misc",
                "--cssClass","feature3",
                "--type","misc_feature",
                "--autocomplete","all"
                "--getLabel",
                "--getType",
                "--getSubs",
                "--getPhase"
                ],
                env=env,
                stdout=stdout,
                stderr=stderr)

            check_call(["generate-names.pl","--out",jbrowse_out_dir],
                env=env,
                stdout=stdout,
                stderr=stderr)

        tracks_conf_file = os.path.join(jbrowse_out_dir,"trackList.json")
        tracks_conf = util.load_config_json(tracks_conf_file)
        tracks_conf["refSeqDropdown"] = True #show pull-down even for very many sequences
        util.save_config_json(tracks_conf,tracks_conf_file)
        #create index.html that redirects to JBrowse index.html with correct data param etc
        _jbrowse_dataset_index_html = \
                config.get_data_string(self.opt["jbrowse_galaxy_index_html_tpl"],
                    "galaxy.index.html")
        jbrowse_url_params = util.to_url_params(dict(
            tracks=",".join(("DNA","Genes","CDS","Peptides","Misc")),
            tracklist=0
            ))
        with open(index_html,"w") as f:
            f.write(_jbrowse_dataset_index_html.\
                format(jbrowse_url=jbrowse_url.rstrip("/"),
                    jbrowse_data_subdir=opt["jbrowse_data_subdir"],
                    jbrowse_url_params=jbrowse_url_params))
