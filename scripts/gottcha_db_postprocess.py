#!/usr/bin/env python
import argparse as ap, textwrap as tw
import pandas as pd
import numpy as np
import subprocess
from tqdm import tqdm
import re
import sys
import os

def parse_params():
    p = ap.ArgumentParser( prog='gottcha_db_build.py', description="""Building GOTTCHA2 database.""" )

    #p.add_argument(
    #    '-dp', '--dbPath', metavar='[PATH]', type=str, default=None,
    #    help="""Path of taxonomy databases [default: <BIN_DIR>/database/]""")

    p.add_argument(
        '-tar', '--targetAssemblyReport', metavar='[PATH]', type=str, default=None, required=True,
        help="""Get target genomes from assembly report file [default: None]""")

    p.add_argument(
        '-l', '--list', metavar='[PATH]', type=str, default=None, required=True,
        help="""Get all input genomes from gottcha_db list [default: None]""")

    p.add_argument(
        '-gdb', '--gottchaDatabaseDir', metavar='[PATH]', type=str,
        help="""local directory for FASTA files [default: <PREFIX>_refseq_genomes]""")

    p.add_argument(
        '-r', '--rank', metavar='[STR]', type=str, default='species',
        help="""Rank of signatures [default: species]""")

    p.add_argument(
        '-p', '--prefix', metavar='[STR]', type=str, default="gottcha",
        help="""Output file [default: gottcha]""")


    args_parsed = p.parse_args()

    #if not args_parsed.dbPath:
    #    bin_dir = os.path.dirname(os.path.realpath(__file__))
    #    args_parsed.dbPath = bin_dir + "/database"

    return args_parsed

def parse_asr(t_asr, glist, gdb_dir, rank):
    #split each line in assembly_summary_refseq.txt:
    #  0- 4  assembly_accession    bioproject      biosample         wgs_master           refseq_category
    #  5- 9  taxid                 species_taxid   organism_name     infraspecific_name   isolate 
    # 10-14  version_status        assembly_level  release_type      genome_rep           seq_rel_date
    # 15-19  asm_name              submitter       gbrs_paired_asm   paired_asm_comp      ftp_path        
    # 20-21  excluded_from_refseq  relation_to_type_material

    def commonprefix(m):
        """Given a list of pathnames, returns the longest common leading component"""
        if not m: return ''
        s1 = min(m)
        s2 = max(m)
        for i, c in enumerate(s1):
            if c != s2[i]:
                return s1[:i]
        return s1

    def conv_basename(path):
        return path.split('/')[-1]+"_genomic.fna.gz" if path.startswith('ftp') else path.split('/')[-1]

    def get_real_signature_file_path(path):
        filename = path.split('/')[-1]
        local_filepath = "%s/%s/%s-%s"%(gdb_dir, path.lstrip(cprefix).rstrip(filename), rank, filename)
        return local_filepath

    sys.stderr.write( "[INFO] Loading %s..."%t_asr )
    t_asr_df = pd.read_csv(t_asr, sep='\t', 
                           low_memory=False,
                           skip_blank_lines=True,
                           comment='#',
                           header=None,
                           names=['assembly_accession', 'bioproject', 'biosample', 'wgs_master', 'refseq_category', 'taxid', 'species_taxid', 'organism_name', 'infraspecific_name', 'isolate', 'version_status', 'assembly_level', 'release_type', 'genome_rep', 'seq_rel_date', 'asm_name', 'submitter', 'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path', 'excluded_from_refseq', 'relation_to_type_material'],
                           )
    t_asr_df = t_asr_df.dropna(subset=['ftp_path'])
    sys.stderr.write( "%s target genome(s) loaded.\n"%len(t_asr_df.index) )

    sys.stderr.write( "[INFO] Loading %s..."%glist )
    l_df = pd.read_csv(glist, sep='\t', header=None, names=['name','path','topo','taxid'])
    cprefix = commonprefix( l_df['path'].tolist() )
    cprefix = re.sub(r'/\w*$', '/', cprefix)
    l_df = l_df[l_df['name']!='SQUASH']
    sys.stderr.write( "%s gottcha_db genome(s) loaded.\n"%len(l_df.index) )

    sys.stderr.write( "[INFO] Converting path to basename..." )
    # convert ftp path to file basename for target ASR file
    t_asr_df['fn'] = t_asr_df['ftp_path'].apply( conv_basename )
    # get file basename for gottcha_db list
    l_df['fn'] = l_df['path'].apply(lambda x: x.split('/')[-1])
    sys.stderr.write( "Done\n" )

    # get genomes within the target list
    sys.stderr.write( "[INFO] Getting genomes within the target list..." )
    l_df = l_df.loc[l_df['fn'].isin(t_asr_df['fn'])]
    sys.stderr.write( "%s genome(s) found.\n"%len(l_df.index) )

    sys.stderr.write( "[INFO] Converting to real signature path..." )
    l_df['path'] = l_df['path'].apply( get_real_signature_file_path )
    l_df = l_df.dropna(subset=['path'])
    sys.stderr.write( "%s file(s) verified.\n"%len(l_df.index) )

    return l_df

if __name__ == '__main__':
    argvs = parse_params()
    df = parse_asr(argvs.targetAssemblyReport, argvs.list, argvs.gottchaDatabaseDir, argvs.rank)

    outfile = "%s.%s.fna.gz"%(argvs.prefix, argvs.rank)
    if os.path.isfile(outfile):
        os.remove(outfile)

    cnt=0
    all_cnt=len(df.index)
    filelist = df['path'].tolist()

    sys.stderr.write( "[INFO] Combining signatures...\n" )
    
    for path in filelist:
        cnt+=1
        subprocess.run("cat %s >> %s"%(path, outfile), shell=True, check=True)
        if not cnt%1000:
            sys.stderr.write( "[INFO] %.1f%% processed.\r"%(cnt/all_cnt*100) )
    sys.stderr.write( "[INFO] %.1f%% processed.\n"%(cnt/all_cnt*100) )

    sys.stderr.write( "\nAll done!\n" )
