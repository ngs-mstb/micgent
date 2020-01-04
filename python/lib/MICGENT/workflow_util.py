"""Utility functions around generating and transforming input and output files in workflows"""
from MICGENT.py23 import *
from . import util
import argh
import os, re, glob
import warnings
import functools

def make_samp_id_extractor(samp_id_extractor):
    samp_id_extractor = util.none_from_str(samp_id_extractor)
    if samp_id_extractor:
        if not '(' in samp_id_extractor:
            samp_id_extractor = "({})".format(samp_id_extractor)
        def _extract_samp_id(patt, s):
            samp_id = re.findall(patt, s)
            if len(samp_id) < 1:
                return None
            if len(samp_id) > 1 or not samp_id[0].strip():
                raise ValueError("Pattern '{}' for samp_id_extractor should extract one " \
                                 "or zero non blank matches from '{}', but it extracted this: '{}'". \
                                 format(patt, s, samp_id))
            samp_id = samp_id[0]
            if len(samp_id.split()) > 1:
                raise ValueError("Pattern '{}' for samp_id_extractor extracted from '{}' a string with spaces: '{}'". \
                                 format(patt, s, samp_id))                
            return samp_id
        return functools.partial(_extract_samp_id,samp_id_extractor)
    else:
        return None

def iterate_sequence_run(input_dir,
                         input_seqfile_ext=".fastq.gz",
                         forward="_R1_",
                         reverse="_R2_",
                         samp_id_extractor=r'^[^_]+',
                         dir_name_is_sample_id=False,
                         sample_id_prefix="",
                         to_abspath=False,
                         out_file=None,
                         out_file_format="json",
                         out_file_sep="tab",
                         samp_id_nomatch="error",
                         from_dir=None,
                         allowed_roots=None,
                         out_file_rejected=None,
                         require_friendly_paths=True,
                         allow_repeated_sample_id=False,
                         files_per_group=2):
    """
    Iterate through files in a sequencing batch.

    Either one of the two layouts should be present:

    - each sample is in its own subdirectory under input_dir. Subdirectory name 
      will be used as sample ID (optionally processed with pattern extraction with
      samp_id_extractor). Files in the subdirectory will be selected by a glob
      *.input_seqfile_ext

    - input_dir matches many files. Sample ID will be extracted from each file name
      with samp_id_extractor, files will be grouped by sample ID

    Once the files are grouped by sample ID in either layout, they will be ordered
    by forward and reverse patterns and returned

    @param input_dir Glob as in file_globs parameter to @see util.glob_files,
    listing either subdirectories for each sample or sample files

    @param input_seqfile_ext extension of sequence files, used only in sudirectory 
    layout

    @param forward pattern for forward reads file; can match anywhere inside file name

    @param reverse pattern for reverse reads file; can match anywhere inside file name

    @param to_abspath return all paths as abspath

    @param samp_id_extractor regular expression to extract sample ID from file name,
    optional in subdirectory layout but mandatory otherwise. Samples where this does 
    not matched will either cause exception and be skipped, depending on the value
    of samp_id_nomatch.

    @param dir_name_is_sample_id base dir name can be used as sample ID

    @param sample_id_prefix add this extra prefix to sample ID after the latter is extracted

    @param out_file write results to this file

    @param out_file_format format for the output file ('json' as list of objects,
    or 'csv' as a sample sheet)

    @param out_file_sep ['tab','comma'] If output file format is 'csv', use tab or comma
    as a field separator.

    @param samp_id_nomatch ["error","skip","warn"] What to do when samp_id_extractor
    does not match anything in the file name. "error" is the default and causes
    the exception to be raised; "skip" ignores non-matching files - use it when
    you want to use the samp_id_extractor to filter samples; "warn" generates
    a warning and skips the file.

    @param require_friendly_paths Generate an exception if any of the final output paths
    are non Linux-friendly (contain spaces spaces or other characters that can be interpreted
    by the shell unless quoted. This to ensure that these paths behave as expected when
    passed to various analysis tools that might be implemented as shell scripts that do not
    take care of quoting their arguments.

    @param files_per_group How many files to expect per group, where the group is defined
    as one unique sample ID per containing directory. The default of 2 corresponds to
    paired-end sequencing. Allowed values are [1,2].

    @return iterator of {samp_id,samp_dir,seqfiles=[forward_file,reverse_file])} if
    out_file is None, else list of such records dumped into the out_file
    """

    assert files_per_group in (1,2), "Allowed target number of files per group of (directory, sample ID) should be 1 or 2."

    assert out_file_format in ("json","csv"), "The only supported formats for the output file are json and csv"
    assert out_file_sep in ("tab", "comma"), "The only supported separator specifications are tab and comma"

    sample_id_prefix = str(sample_id_prefix).strip() if sample_id_prefix else ''
    if sample_id_prefix:
        util.assert_id_friendly(sample_id_prefix)

    if not from_dir is None:
        assert util.is_rel_subdir(input_dir), "Input path {} must conform to a format of simple relative subdirectory".\
            format(input_dir)
    
    if not allowed_roots is None:
        allowed_roots = [ util.absrealpath(_) for _ in util.glob_files(files_globs=allowed_roots) ]

    samp_id_extractor_str = samp_id_extractor
    samp_id_extractor = make_samp_id_extractor(samp_id_extractor)

    with util.chdir(from_dir):

        seq_samps = {}
        samp_match_rejected = []
        samp_match_error_msg = None
        for samp_match in util.glob_files(files_globs=input_dir):

            if os.path.isdir(samp_match):

                rec = {}

                rec["samp_dir"] = samp_match

                samp_seqfiles = glob.glob(os.path.join(samp_match, "*" + input_seqfile_ext))
                assert len(samp_seqfiles), \
                    "No sequence files found for subdirectory {}".format(samp_match)

                rec["samp_seqfiles"] = samp_seqfiles

                assert dir_name_is_sample_id
                samp_id = os.path.basename(samp_match)
                if samp_id_extractor:
                    samp_id = samp_id_extractor(samp_id)
                if samp_id:
                    dir_recs = seq_samps.setdefault(samp_id,dict())
                    dir_recs[rec["samp_dir"]] = rec

            else:
                samp_id = os.path.basename(samp_match)
                assert samp_id_extractor
                samp_id = samp_id_extractor(samp_id)
                if samp_id:
                    samp_dir = os.path.dirname(samp_match)
                    dir_recs = seq_samps.setdefault(samp_id,dict())
                    rec = dir_recs.setdefault(samp_dir,dict(samp_dir=samp_dir))
                    samp_seqfiles = rec.setdefault("samp_seqfiles",list())
                    samp_seqfiles.append(samp_match)
                else:
                    msg = "Unable to extract sample ID with pattern {} from the file name {}".\
                                                        format(samp_id_extractor_str,samp_match)
                    if samp_id_nomatch == "error":
                        samp_match_error_msg = msg
                        if not out_file_rejected:
                            raise ValueError(samp_match_error_msg)
                        else:
                            samp_match_error_msg += " See file '{}' for all unmatched file names".format(out_file_rejected)
                    elif samp_id_nomatch == "warn":
                        warnings.warn(msg)
                    if out_file_rejected:
                        samp_match_rejected.append(samp_match)

        if out_file_rejected:
            if out_file_format == "json":
                import json
                with util.open_text_py23(out_file_rejected, "w") as _:
                    json.dump(samp_match_rejected, _)
            elif out_file_format == "csv":
                with util.open_text_py23(out_file_rejected, "w") as _:
                    _.write("file\n")
                    _.write("\n".join(samp_match_rejected))

        if samp_id_nomatch == "error":
            if samp_match_error_msg:
                raise ValueError(samp_match_error_msg)

        out_res = []
        for samp_id,dir_recs in seq_samps.items():

            if not allow_repeated_sample_id:
                assert len(dir_recs) <= 1, "Sample ID {} is not unique".format(samp_id)

            samp_id = "{}".format(sample_id_prefix)+samp_id

            for samp_dir, samp_rec in dir_recs.items():
                samp_seqfiles = samp_rec["samp_seqfiles"]
                assert len(samp_seqfiles) == files_per_group, \
                    "Expected number of {} files are not found for sample ID {} in dir {}, found instead these files: {}".\
                        format(files_per_group, samp_id, samp_dir, samp_seqfiles)
                if len(samp_seqfiles) == 2:
                    ## reorder files in forward, reverse pattern order
                    for ff in (samp_seqfiles, tuple(reversed(samp_seqfiles))):
                        if re.search(forward, ff[0]):
                            assert re.search(reverse, ff[1]), "Forward read file name pattern matched, but reverse pattern did not: {}".format(ff)
                            samp_seqfiles = ff
                            break
                    else:
                        raise ValueError("Did not find forward and reverse file name patterns: {}".format(samp_seqfiles))
                if not allowed_roots is None:
                    for f in (samp_seqfiles + [samp_rec["samp_dir"]]):
                        assert util.is_real_subdir_any(f,allowed_roots,_upper_prepared=True),\
                            "Real path of {} is outside of allowed roots {}".format(f,allowed_roots)
                if to_abspath:
                    samp_seqfiles = [os.path.abspath(_) for _ in samp_seqfiles]
                    samp_rec["samp_dir"] = os.path.abspath(samp_rec["samp_dir"])
                if require_friendly_paths:
                    for f in samp_seqfiles:
                        util.assert_path_unix_friendly(f)

                rec = dict(samp_dir=samp_rec["samp_dir"], seqfiles=samp_seqfiles, samp_id = samp_id)
                out_res.append(rec)

    if out_file is None:
        for rec in out_res:
            yield rec
    else:
        if out_file_format == "json":
            import json
            with util.open_text_py23(out_file,"w") as _:
                json.dump(out_res,_)
        elif out_file_format == "csv":
            import csv
            fields = ("SampleID","file1","file2")
            with util.open_text_py23(out_file,"w") as _:
                out = csv.DictWriter(_,fieldnames=fields,
                        dialect="excel-tab" if out_file_sep=="tab" else "excel",
                                     lineterminator="\n")
                out.writeheader()
                for rec in out_res:
                    row = dict(
                            file1=rec["seqfiles"][0],
                            file2=rec["seqfiles"][1],
                            SampleID=rec["samp_id"]
                            )
                    out.writerow(row)


def join_files_with_ids(id_file,reads_file,out_file):
    """Join sample ID file with file names file.
    This is a helper method for building inputs to WDL tasks.
    @param id_file Headless tab-delimited table - sample ID per line
    @param reads_file Headless tab-delimited table - forward and reverse file names per line
    @param out_file Tab-delimited table with header SampleID,file1,file2. Records
    from the input table are joined line-by-line.
    
    Input files are assumed to be created by WDL write_tsv. Output should be read
    by WDL read_objects. This services a WDL coding pattern when multiple outputs
    are returned from inside a `scatter` block, they are passed to another
    task, and there they are joined by this function."""

    with open(id_file,'r') as id_inp,\
        open(reads_file,'r') as reads_inp,\
        util.open_text_py23(out_file,'w') as out:
            out.write("\t".join(["SampleID","file1","file2"])+"\n")
            for  samp_id in id_inp:
                samp_id = samp_id.strip()
                reads = reads_inp.readline()
                out.write("{}\t{}".format(samp_id,reads))

def load_manifest(manifest,sep="\t",dtype={"SampleID":str,"SeqID":str, "SeqID_RADAR":str},*l,**kw):
    """Load a sample manifest file as pandas.DataFrame.

    This is a thin wrapper around Pandas.DataFrame.
    The default list of fixed string dtypes can be expanded - field names that are not
    present in the data frame are ignored by pandas.read_table. The objective of using
    this argument here is to preserve various ID fields that should be always treated 
    as strings.
    """
    import pandas as pd
    dtype = dtype.copy()
    return pd.read_table(manifest,dtype=dtype,sep=sep,*l,**kw)

def manifest_to_cwl(manifest,cwl,datadir=None,annot_sep="\t",require_friendly_paths=True,
    allowed_roots=None):
    from ruamel import yaml
    import pandas as pd
    
    if not allowed_roots is None:
        allowed_roots = [ util.absrealpath(_) for _ in util.glob_files(files_globs=allowed_roots) ]

    df = pd.read_table(manifest,sep=annot_sep)
    
    def create_cwl_rec(row,datadir=None,allowed_roots=None):
        file1 = row.file1
        file2 = row.file2
        if datadir:
            assert util.is_rel_subdir(file1) and \
                   util.is_rel_subdir(file2), "Received unsafe file paths: {}".format(row)
            file1 = os.path.join(datadir,file1)
            file2 = os.path.join(datadir,file2)
        if require_friendly_paths:
            for f in (file1,file2):
                util.assert_path_unix_friendly(f)

        if not allowed_roots is None:
            for f in (file1,file2):
                assert util.is_real_subdir_any(f,allowed_roots,_upper_prepared=True),\
                    "Real path of {} is outside of allowed roots {}".format(f,allowed_roots)

        out = {}
        out['file1'] = {}
        out['file1']['class'] = 'File'
        out['file1']['path'] = file1

        out['file2'] = {}
        out['file2']['class'] = 'File'
        out['file2']['path'] = file2

        out['SampleID'] = str(row.SampleID)

        return out

    res = [create_cwl_rec(row,datadir=datadir,allowed_roots=allowed_roots) for row in df.itertuples()]
    with util.open_text_py23(cwl,'w') as out:
        out.write(yaml.dump(res))

def manifest_multi_read_to_cwl(manifest,cwl,datadir=None,annot_sep="\t",require_friendly_paths=True):
    from ruamel import yaml
    import pandas as pd
    df = pd.read_table(manifest,sep=annot_sep)
    def create_cwl_rec(SampleID,rows,datadir=None):
        files1 = rows.file1
        files2 = rows.file2
        out = {}
        out['files1'] = []
        out['files2'] = []

        for file1, file2 in zip(files1, files2):
            if datadir:
                assert util.is_rel_subdir(file1) and \
                       util.is_rel_subdir(file2), "Received unsafe file paths: {}".format(row)
                file1 = os.path.join(datadir,file1)
                file2 = os.path.join(datadir,file2)
            if require_friendly_paths:
                for f in (file1,file2):
                    util.assert_path_unix_friendly(f)
            out['files1'].append({'class':'File','path':file1})
            out['files2'].append({'class':'File','path':file2})

        out['SampleID'] = str(SampleID)

        return out

    res = [create_cwl_rec(SampleID,rows,datadir=datadir) for SampleID,rows in df.groupby("SampleID")]
    with util.open_text_py23(cwl,'w') as out:
        out.write(yaml.dump(res))



def _prep_manifest_output_file(fname,datadir=None,require_friendly_paths=True):
    if datadir:
        assert util.is_rel_subdir(fname), "Received unsafe file path: {}".format(fname)
        fname = os.path.join(datadir, fname)
    if require_friendly_paths:
        util.assert_path_unix_friendly(fname)
    return fname

def _path_to_cwl(fname):
    return {"path":fname,"class":"File"}

def manifest_to_atacseq_cwl(manifest, cwl, samp_id_exp_extractor, datadir=None, annot_sep="\t", require_friendly_paths=True):
    """Convert flat read manifest into CWL input for ATACSeq experiment with replicates"""
    from ruamel import yaml
    import pandas as pd

    samp_id_exp_extractor_str = samp_id_exp_extractor
    samp_id_exp_extractor = make_samp_id_extractor(samp_id_exp_extractor)

    df = pd.read_table(manifest, sep=annot_sep)

    def create_cwl_rec(SampleIDExp,rows, datadir=None):
        out = dict(SampleID=SampleIDExp)
        for i_rep, (SampleID, rep_rows) in enumerate(rows.groupby("SampleID")):
            #import pdb; pdb.set_trace()
            file1 = rep_rows.file1.apply(lambda x: _prep_manifest_output_file(x,datadir,require_friendly_paths))
            file2 = rep_rows.file2.apply(lambda x: _prep_manifest_output_file(x,datadir,require_friendly_paths))
            for i_read, f_read in enumerate((file1,file2)):
                out['files{}_{}'.format(i_rep+1,i_read+1)] = [ _path_to_cwl(fname) for fname in f_read ]
        return out

    df["SampleIDExp"] = df.SampleID.apply(samp_id_exp_extractor)
    res = [create_cwl_rec(SampleIDExp,rows, datadir=datadir) \
           for SampleIDExp, rows \
           in df.groupby("SampleIDExp")]
    with util.open_text_py23(cwl, 'w') as out:
        out.write(yaml.dump(res,Dumper=yaml.RoundTripDumper,default_flow_style=False))

def manifest_to_atacseq_wdl(manifest, inp_json, out_json, samp_id_exp_extractor,
                            samples_key="atac_many.sample_fastqs",
                            datadir=None, annot_sep="\t",
                            require_friendly_paths=True):
    """Convert flat read manifest into JSON WDL input for ATACSeq experiment with replicates.
    Samples are represented as:
      Array[Pair[String,Array[Array[Array[File]]]]] sample_fastqs # samp_id => [rep_id][merge_id][read_end_id]
    """
    import pandas as pd

    wd_inp = util.load_json(inp_json)
    samp_id_exp_extractor = make_samp_id_extractor(samp_id_exp_extractor)

    df = pd.read_table(manifest, sep=annot_sep)

    def create_wdl_rec(SampleIDExp, rows, datadir=None):
        ## dict(Left:X,Right:Y) is converted into Pair[X,Y] in WDL
        out = dict(Left=SampleIDExp)
        files = []
        for i_rep, (SampleID, rep_rows) in enumerate(rows.groupby("SampleID")):
            file1 = rep_rows.file1.apply(lambda x: _prep_manifest_output_file(x, datadir, require_friendly_paths))
            file2 = rep_rows.file2.apply(lambda x: _prep_manifest_output_file(x, datadir, require_friendly_paths))
            files.append([ fr for fr in zip(file1,file2) ])
        out["Right"] = files
        return out

    df["SampleIDExp"] = df.SampleID.apply(samp_id_exp_extractor)
    res = [create_wdl_rec(SampleIDExp, rows, datadir=datadir) \
           for SampleIDExp, rows \
           in df.groupby("SampleIDExp")]
    wd_inp[samples_key] = res
    util.save_json(wd_inp,out_json,pretty=True)


## import package module and add argh entry points

def _main():
    from . import arg_parsing
    parser = arg_parsing.ArghParserChainedConfig()
    parser.add_commands([
        iterate_sequence_run,
        join_files_with_ids,
        manifest_to_cwl,
        manifest_to_atacseq_cwl,
        manifest_to_atacseq_wdl,
        manifest_multi_read_to_cwl
    ])
    parser.dispatch()

if __name__ == "__main__":
    _main()




