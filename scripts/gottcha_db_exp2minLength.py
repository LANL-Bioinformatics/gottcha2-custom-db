#!/usr/bin/env python
from pyfaidx import Fasta
from tqdm import tqdm
import multiprocessing as mp
import argparse as ap, textwrap as tw
import pandas as pd
import numpy as np
import subprocess
import re
import sys
import os

def parse_params():
    p = ap.ArgumentParser( prog=sys.argv[0], description="""Expanding signatures to a minimal length.""" )

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
        '-e', '--expand', metavar='[INT]', type=int, default=120,
        help="""Expand signature to at a minimal length (bp) [default: 120]""")

    p.add_argument(
        '-c', '--cpu', metavar='[CPU_NUMBER]', type=int, default=2,
        help="""Number of CPU [default: 2]""")

    p.add_argument(
        '-sp', '--secondpass', action='store_true',
        help="""2nd pass gottcha signature [default: false]""")

    p.add_argument(
        '-p', '--prefix', metavar='[STR]', type=str, default="gottcha",
        help="""Output file [default: gottcha]""")


    args_parsed = p.parse_args()

    return args_parsed

def parse_asr(t_asr, glist, gdb_dir, rank, second_pass):
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
        if second_pass:
            local_filepath = "%s/%s/%s-%s-%s"%(gdb_dir, path.lstrip(cprefix).rstrip(filename), rank, rank, filename)
        else:
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
    sys.stderr.write( "%s gottcha_db genome(s) loaded.\n"%len(l_df.index) )

    sys.stderr.write( "[INFO] Converting path to basename..." )
    # convert ftp path to file basename for target ASR file
    t_asr_df['fn'] = t_asr_df['ftp_path'].apply( conv_basename )
    # get file basename for gottcha_db list
    l_df['fn'] = l_df['path'].apply(lambda x: x.split('/')[-1])
    # for 2nd pass
    l_df['fn_orig'] = l_df['fn'].str.replace('%s-'%rank,'')

    sys.stderr.write( "Done\n" )

    # get genomes within the target list
    sys.stderr.write( "[INFO] Getting genomes within the target list..." )
    l_df = l_df.loc[l_df['fn_orig'].isin(t_asr_df['fn'])]
    sys.stderr.write( "%s genome(s) found.\n"%len(l_df.index) )

    sys.stderr.write( "[INFO] Converting to real signature path..." )
    l_df['genome_path'] = l_df['path']
    l_df['path'] = l_df['path'].apply( get_real_signature_file_path )
    l_df = l_df.dropna(subset=['path'])
    sys.stderr.write( "%s file(s) verified.\n"%len(l_df.index) )

    return l_df

def depressSeq(filename_gz):
    filename_fa=filename_gz.replace(".gz","")
    if os.path.isfile(filename_gz):
        if not os.path.isfile(filename_fa):
            subprocess.run("gzip -dc %s > %s"%(filename_gz,filename_fa), shell=True, check=True)
        return filename_fa
    else:
        sys.stderr.write( "[WARNNING] %s not found.\n"%filename_gz )
        return None

def output_expanded_sig(genome_gz, sig_gz, expand_len):
    genome_fa = depressSeq(genome_gz)
    sig_fa    = depressSeq(sig_gz)
    fi_genome = Fasta(genome_fa)
    fi_sig    = None
    outbuff   = ""

    # if genomes or signatures not found
    if not (genome_fa and sig_fa):
        return sig_fa,outbuff

    if os.path.getsize(sig_fa) > 0:
        fi_sig = Fasta(sig_fa)
    else:
        return sig_fa,outbuff

    def expandPos(start, end, genome_end):
        """
        Recalculate signature positions if the lenght < expand_len
        """
        slen = end-start+1
        if slen < expand_len:
            offset = int((expand_len-slen)/2)
            if start-offset < 1:
                return [1, expand_len]
            elif end+expand_len-slen-offset > genome_end:
                return [genome_end-expand_len+1, genome_end]
            else:
                return [start-offset, end+expand_len-slen-offset]
        else:
            return [start, end]

    def getSeq(pf_genome, gid, es, ee):
        return pf_genome[gid][es:ee].seq

    for record in fi_sig:
        (sid, start, end, taxid) = (0,0,0,0)
        tmp = record.name.split('|')

        # for 2nd pass signatures
        if len(tmp)==8:
            (sid, start, end, taxid, ingrp_start, ingrp_end, taxid2, null) = record.name.split('|')
            start = int(start)+int(ingrp_start)-1
            end   = int(start)+int(ingrp_end)-1
        # for original signatures
        else:
            (sid, start, end, taxid, null) = record.name.split('|')

        (exp_start, exp_end) = expandPos(int(start), int(end), len(fi_genome[sid]))
        seq = fi_genome[sid][(exp_start-1):exp_end].seq
        outbuff += ">%s|%s|%s|%s|\n%s\n"%(sid, exp_start, exp_end, taxid, seq)

    f = open("%s.exp120"%sig_fa, "w")
    f.write(outbuff)
    f.close()

    return sig_fa, outbuff

if __name__ == '__main__':
    argvs = parse_params()
    df = parse_asr(argvs.targetAssemblyReport, argvs.list, argvs.gottchaDatabaseDir, argvs.rank, argvs.secondpass)

    outfile = "%s.%s.fna"%(argvs.prefix, argvs.rank)
    if os.path.isfile(outfile):
        os.remove(outfile)

    sys.stderr.write( "[INFO] Processing signatures...\n" )

    filelists = df[['genome_path','path']].values.tolist()

    ## debug code
    #cnt=0
    #tol_jobs=1
    #for (genome_path, sig_path) in filelists:
    #    (sig_fa, outbuff) = output_expanded_sig(genome_path, sig_path, argvs.expand)
    #    cnt+=1
    #    sys.stderr.write( "[INFO] Progress: %s/%s (%.1f%%) %s done.\n"%(cnt, tol_jobs, cnt/tol_jobs*100, sig_fa))

    #pool
    pool = mp.Pool(processes=argvs.cpu)
    jobs = [pool.apply_async(output_expanded_sig, (genome_path, sig_path, argvs.expand)) for (genome_path, sig_path) in filelists ]
    tol_jobs = len(jobs)
    sys.stderr.write( "[INFO] %s jobs pushed, processing in %s CPU(s)...\n"%(tol_jobs, argvs.cpu) )

    #wait for all jobs to finish
    cnt = 0
    f = open(outfile, "w")
    for job in jobs:
        (sig_fa, outbuff) = job.get()
        f.write(outbuff) #flushing
        cnt+=1
        sys.stderr.write( "[INFO] Progress: %s/%s (%.1f%%) %s done.\n"%(cnt, tol_jobs, cnt/tol_jobs*100, sig_fa))

    #wait
    pool.close()
    f.close()

    sys.stderr.write( "\nAll done!\n" )
