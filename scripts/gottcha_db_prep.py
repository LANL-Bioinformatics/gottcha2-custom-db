#!/usr/bin/env python
import argparse as ap, textwrap as tw
import multiprocessing as mp
import taxonomy as t
import time
import copy
import gzip
import sqlite3
import pandas as pd
import sys
import os
import subprocess
from re import search
from io import StringIO

def parse_params():
	p = ap.ArgumentParser( prog=sys.argv[0], description="""Preparing genome sequences for GOTTCHA2 database and creating lists for gottcha_db.""" )

	p.add_argument(
		'-ar', '--assemblyReport', metavar='[PATH]', type=str, default=None,
		help="""Use assembly report to define strain [default: None]""")

	p.add_argument(
		'-dp', '--dbPath', metavar='[PATH]', type=str, default=None,
		help="""Path of taxonomy databases [default: <BIN_DIR>/taxanomy_db/]""")

	p.add_argument(
		'-sdb','--sqlitedb', metavar='[FILE]', nargs='?',
		help="Custom SQLite3 taxonomy db [default: <PREFIX>.custom.db]")

	p.add_argument(
		'-scf', '--skipCheckingFiles', action="store_true",
		help="Skip checking files [default: None]")

	p.add_argument(
		'-asd', '--assemblyRefseqDir', metavar='[PATH]', type=str,
		help="""local directory for FASTA files [default: <PREFIX>_refseq_genomes]""")

	p.add_argument(
		'-rc', '--refseqCategoryOnly', nargs='*', metavar='[STR]', type=str,
		help="""Only parse specified refseq category(ies), like 'reference' or 'representative'. [default: None]""")

	p.add_argument(
		'-al', '--assemblyLevelOnly', nargs='*', metavar='[STR]', type=str,
		help="""Only parse specified assembly level(s), like 'complete'. [default: None]""")

	p.add_argument(
		'-ssf', '--strainSeparateFasta', action="store_true",
		help="""Each FASTA file represents a strain. [default: FALSE]""")

	p.add_argument(
		'-sr', '--splitDatabaseByRank', metavar='[RANK]', type=str, default="superkingdom",
		help="""Split database by rank [default: superkingdom]""")

	p.add_argument(
		'-tid', '--taxid', metavar='[TAXID]', nargs='*', type=str,
		help="""Only process genomes belong to specified taxid(s). The option "--excludeTaxid" is prioritied over "--taxid". [default: None]""")

	p.add_argument(
		'-xt', '--excludeTaxid', metavar='[TAXID]', nargs='*', type=str, default=None,
		help="""Exclude genomes belong to specified taxid(s). The option "--excludeTaxid" is prioritied over "--taxid". [default: None]""")

	p.add_argument(
		'-ngl','--noGottchaDbList', action="store_true",
		help="""Do not generate gottcha_db input list. [default: False]""")

	p.add_argument(
		'-p', '--prefix', metavar='[STR]', type=str, default="gottcha_db",
		help="""Output prefix [default: gottcha_db]""")

	p.add_argument(
		'-t', '--cpu', metavar='[CPU_NUMBER]', type=int, default=2,
		help="""Number of CPU [default: 2]""")

	p.add_argument(
		'-rs', '--redundantStrain', action="store_true",
		help="""Allow redundant strains [default: FALSE]""")

	p.add_argument(
		'-v', '--verbose', action="store_true",
		help="""Verbose message [default: FALSE]""")

	p.add_argument(
		'-iso', '--includeIsolate', action="store_true",
		help="""Different isolates will be treated as differnet strains [default: FALSE]""")

	args_parsed = p.parse_args()

	if not args_parsed.dbPath:
		bin_dir = os.path.dirname(os.path.realpath(__file__))
		args_parsed.dbPath = bin_dir + "/taxanomy_db"

	if not args_parsed.assemblyRefseqDir:
		args_parsed.assemblyRefseqDir = args_parsed.prefix + "ncbi_genomes"

	if not args_parsed.sqlitedb:
		args_parsed.sqlitedb = args_parsed.prefix + ".custom.db"

	return args_parsed

def getCustomStrain(taxid, name, desc, c, str_name, iso_name, include_iso):
	# get full custom strain name
	cus_name = name
	if str_name and not str_name in cus_name:
		cus_name += " strain " + str_name
	if iso_name and include_iso and not iso_name in cus_name:
		cus_name += " isolate " + iso_name

	#check custom name's existence
	c.execute( 'SELECT * FROM taxonomy_custom WHERE NAME="%s"' % cus_name )
	tax_rec = c.fetchone()

	if tax_rec:
		cus_taxid = tax_rec[0]
		new_strain = False
	else:
		cus_taxid = createNewCusTaxRecord( c, taxid, cus_name )
		new_strain = True

	return cus_taxid, cus_name, str_name, iso_name, new_strain

def createNewCusTaxRecord( c, taxid, cus_name, retry=0 ):
	#insert new custom name
	cus_taxid = ""
	str_id = ""
	depth = ""
	c.execute('SELECT * FROM taxonomy_custom WHERE P_TAXID="%s" ORDER BY STR_ID DESC' % taxid )
	tax_rec = c.fetchone()
	# if parent taxid (usually a species) already in the database
	if tax_rec:
		str_id = tax_rec[2]+1
		depth  = tax_rec[4]
	else:
		str_id = 1
		depth  = int(t.taxid2depth(taxid)) + 1

	cus_taxid = "%s.%s" % (taxid, str_id)
	c.execute('INSERT INTO taxonomy_custom VALUES ("%s", "%s", %s, "%s", %s, "%s")' % ( cus_taxid, taxid, str_id, cus_name, depth, "no rank") )

	return cus_taxid

def checkStrainNameExist(c, name):
	c.execute('SELECT * FROM strain_taxonomy WHERE NAME like "%s"' % tid )
	tax_rec = c.fetchone()
	# if another assembly is named the same strain, skip this file (assembly)
	if tax_rec:
		print_message( "[INFO] Strain name exists (Taxid: %s, Name: %s).\n" % (tid, tax_rec[2]) , argvs.verbose, begin_t, logfile )
		return False
	else:
		return True

def checkStrainTaxidExist(c, tid):
	c.execute('SELECT * FROM strain_taxonomy WHERE TAXID="%s"' % tid )
	tax_rec = c.fetchone()
	# if another assembly is named the same strain, skip this file (assembly)
	if tax_rec:
		print_message( "[INFO] Strain exists (Taxid: %s, Name: %s)." % (tid, tax_rec[2]) , argvs.verbose, begin_t, logfile )
		return False
	else:
		return True

def addStrain( c, tid, ptid, name, str_name="", iso_name="", aacc="", filename="" ):
	#    TAXID      CHAR(30)    NOT NULL, /* taxid */
	#    P_TAXID    CHAR(30)    NOT NULL, /* parent taxid */
	#    NAME       CHAR(200)   NOT NULL, /* taxonomy name */
	#    STR_NAME   CHAR(200)   NULL,     /* strain name */
	#    ISO_NAME   CHAR(200)   NULL,     /* isolate name */
	#    ASSEM_ACC  CHAR(50)    NULL,     /* assembly acc# */
	#    FILE       CHAR(50)    NULL,     /* file name */
	#    PRIMARY KEY (TAXID)
	c.execute('INSERT INTO strain_taxonomy VALUES ("%s", "%s", "%s", "%s", "%s", "%s", "%s")' % (
		tid,
		ptid,
		name,
		str_name,
		iso_name,
		aacc,
		filename
		)
	)
	return True

def initTaxaDB(c):
	c.execute("""
	CREATE TABLE IF NOT EXISTS taxonomy_custom(
	   TAXID      CHAR(30)    NOT NULL,
	   P_TAXID    CHAR(20)    NOT NULL,
	   STR_ID     INT         NOT NULL,
	   NAME       CHAR(200)   NOT NULL,
	   DEPTH      INT         NOT NULL,
	   RANK       CHAR(20)    NOT NULL,
	   PRIMARY KEY (TAXID)
	);
	""")

	c.execute("""
	CREATE TABLE IF NOT EXISTS strain_taxonomy(
	   TAXID      CHAR(30)    NOT NULL, /* taxid */
	   P_TAXID    CHAR(30)    NOT NULL, /* parent taxid */
	   NAME       CHAR(200)   NOT NULL, /* taxonomy name */
	   STR_NAME   CHAR(200)   NULL,     /* strain name */
	   ISO_NAME   CHAR(200)   NULL,     /* isolate name */
	   ASSEM_ACC  CHAR(50)    NULL,     /* assembly acc# */
	   FILE       CHAR(50)    NULL,     /* file name */
	   PRIMARY KEY (TAXID)
	);
	""")

	conn.commit()

def processASR(c, tid2lineage, asr, rc, al, asr_dir, target_tid, exclude_tid, skip_checking_files, flag_red_str, flag_ind_iso=0):
	"""
	ARGVS:
	   c             OBJ   sqlite3 connection obj
	   asr           STR   asr file
	   rc            LIST  refseq category
	   al            LIST  assembly level
	   asr_dir       STR   asr local directory
	   target_tid    LIST  target taxid
	   flag_red_str  BOOL  flag for allowing redundant strains
	   flag_ind_iso  BOOL  flag for including isolate in strain name
	RETURN:
	   asr_info      DICT
	"""

	# init vars
	asr_info = t._autoVivification()
	cnt_q    = 0 #qualified assembly
	cnt_tol  = 0 #total genomes

	#df = pd.read_csv(asr, sep='\t')
	#df = df.replace('NA', np.nan)

	#q_df = df.loc[
	#    (df.LEVEL           == rank.value) &
	#    (df.READ_COUNT      >= max_r_raw.value) &
	#    (df.READ_COUNT_RSNB >= max_r_rsnb.value) &
	#    (df.LINEAR_COV      >= min_cov.value) &
	#    (df.SCORE           >= min_score.value) &
	#    (df.DEPTH_COV       >= min_dc.value/1000) &
	#    (df.RS_DEPTH_COV_NR >= min_rsdc.value/1000)
	#]

	# parsing assembly_summary_refseq.txt file
	with open(asr) as f:
		for line in f.readlines():
			if line.startswith('#'):
				continue
			elif line.startswith('assembly_accession'):
				continue
			else:
				cnt_tol += 1

			line = line.strip('\n')

			#split each line in assembly_summary_refseq.txt:
			#  0- 4  assembly_accession    bioproject      biosample         wgs_master           refseq_category
			#  5- 9  taxid                 species_taxid   organism_name     infraspecific_name   isolate 
			# 10-14  version_status        assembly_level  release_type      genome_rep           seq_rel_date
			# 15-19  asm_name              submitter       gbrs_paired_asm   paired_asm_comp      ftp_path        
			# 20-21  excluded_from_refseq  relation_to_type_material
			tmp = line.split('\t')
			filename = ""
			local_path = ""
			lineage = {}

			print_message( "[INFO] Processing: %s %s..."%(cnt_tol,tmp[0]) , argvs.verbose, begin_t, logfile )

			# try to get taxonomy lineage
			try:
				lineage = t.taxid2lineageDICT(tmp[5], 1, 1)
			except:
				print_message( "skipped. Removed TaxID (%s) found for %s." % (tmp[5], tmp[0]) , argvs.verbose, begin_t, logfile )
				continue
			# SKIPPING following records
			# 1) not specified refseq_category, example: reference 
			if rc:
				flag=0
				for cate in rc:
					if cate.lower() in tmp[4].lower():
						flag=1
						break
				if not flag:
					print_message( "skipped. Not belongs to specific refseq_category %s."%rc , argvs.verbose, begin_t, logfile )
					continue
			# 2) not specified assembly_level, example: complete
			if al:
				flag=0
				for l in al:
					if l.lower() in tmp[11].lower():
						flag=1
						break
				if not flag:
					print_message( "skipped. Not belongs to specific assembly_level %s."%al , argvs.verbose, begin_t, logfile )
					continue
			# 3) assemblies that marked "excluded_from_refseq"
			if tmp[20]:
				print_message( "skipped. Marked as excluded_from_refseq." , argvs.verbose, begin_t, logfile )
				continue
			# 4) "only include" or "exclude" specified tax id
			if tmp[5] and target_tid:
				flag_tid_in_lineage=False
				for t_tid in target_tid:
					rank = t.taxid2rank(t_tid)
					if lineage[rank]['taxid'] == t_tid:
						flag_tid_in_lineage=True
						break
				if not flag_tid_in_lineage:
					print_message( "skipped. Not belongs to specified taxid(s)." , argvs.verbose, begin_t, logfile )
					continue

			if tmp[5] and exclude_tid:
				flag_tid_in_lineage=False
				for ex_tid in exclude_tid:
					rank = t.taxid2rank(ex_tid)
					if lineage[rank]['taxid'] == ex_tid:
						flag_tid_in_lineage=True
						break
				if flag_tid_in_lineage:
					print_message( "skipped. Excluded taxtid(s) found." , argvs.verbose, begin_t, logfile )
					continue

			# 5) sequence location is not available
			if tmp[19] != "na":
				if tmp[19].startswith('ftp'):
					# real filename of assembly
					filename = tmp[19].split('/')[-1] + "_genomic.fna.gz"
					# directory structure is retained to keep local filesystem healthy
					local_path = asr_dir + tmp[19].split('genomes')[1]
				else:
					filename = tmp[19].split('/')[-1]
					local_path = "/".join(tmp[19].split('/')[:-1])
			else:
				print_message( "skipped. No URL provided." , argvs.verbose, begin_t, logfile )
				continue

			# 6) not an unique strain
			tid = tmp[5]
			rank = t.taxid2rank(tid)
			name = t.taxid2name(tid)
			if rank != "strain" and rank != "unknown":
				# generate custom strain if not exists and add strains to database
				(cus_str_taxid, cus_str_name, str_name, iso_name, flag_new_strain) = getCustomStrain( tid, name, '', c, tmp[8].replace("strain=",""), tmp[9], flag_ind_iso )
				#add custom strain to lineage
				lineage["strain"]['taxid'] = cus_str_taxid
				lineage["strain"]['name'] = cus_str_name

			flag_new_strain = checkStrainTaxidExist(c, lineage["strain"]['taxid'])

			if flag_new_strain:
				addStrain(c, lineage["strain"]['taxid'], lineage["species"]['taxid'], lineage["strain"]['name'], tmp[8].replace("strain=",""), tmp[9], tmp[0], filename )
			else:
				print_message( "skipped. Not an unique strain." , argvs.verbose, begin_t, logfile )
				continue

			print_message( "qualified." , argvs.verbose, begin_t, logfile )

			# use assembly_accession as key
			cnt_q += 1
			tid = lineage["strain"]['taxid']
			asr_info[tmp[0]]['taxid']         = tid
			asr_info[tmp[0]]['full_str_name'] = lineage["strain"]['name']
			asr_info[tmp[0]]['ftp_path']      = tmp[19]
			asr_info[tmp[0]]['local_path']    = local_path
			asr_info[tmp[0]]['filename']      = filename
			asr_info[tmp[0]]['type_material'] = True if tmp[21] else False
			tid2lineage[tid] = copy.deepcopy(lineage)
	f.close()
	return asr_info, tid2lineage, cnt_q, cnt_tol

def prepareGenomeWorker(ftpurl,local_path,filename):
	local_fn = local_path+"/"+filename
	url = ftpurl+"/"+filename

	if not os.path.isfile( local_fn ) or os.stat(local_fn).st_size == 0:
		wget( url, local_path )
		if not os.path.isfile( local_fn ) or os.stat(local_fn).st_size == 0:
			print_message( "[ERROR] Failed to download genome:\nFrom: %s\nTo: %s."%(url, local_path), argvs.verbose, begin_t, logfile, True )
	else:
		print_message( "[INFO] Found local file: %s."%local_fn , argvs.verbose, begin_t, logfile )

def prepareGenomes(asr_info, cpunum):
	pool = mp.Pool(processes=cpunum)
	for asm in asr_info:
		asm_asr_info = asr_info[asm]
		pool.apply_async(prepareGenomeWorker, (asm_asr_info['ftp_path'],asm_asr_info['local_path'],asm_asr_info['filename']))
	#clean up
	pool.close()
	pool.join()

def generateGottchaDbList(asr_info, output_gdbl_buffer, tid2lineage, prefix):
	ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
	for rank in ranks:
		output_gdbl_buffer=""
		for asm in asr_info:
			asm_asr_info = asr_info[asm]
			tid = asm_asr_info["taxid"]
			lineage = tid2lineage[tid]
			fn = "%s/%s"%(asm_asr_info['local_path'], asm_asr_info['filename'])
			text = "%s\t%s\tlinear\t%s\n"%( lineage[rank]['name'], fn, tid )
			output_gdbl_buffer += text

		with open( "%s.%s.list"%(prefix, rank), "w") as f:
			f.write( output_gdbl_buffer )
			f.close()

def wget(url, lpath="./"):
	subprocess.run("mkdir -p %s"%lpath, shell=True, check=True)
	subprocess.run("wget '%s' -P %s"%(url, lpath), shell=True, check=True)

def init_sqlite_db(db_file):
	# Read database to tempfile
	con = sqlite3.connect(db_file)
	tempfile = StringIO()
	for line in con.iterdump():
		tempfile.write('%s\n' % line)
	con.close()
	tempfile.seek(0)

	# Create a database in memory and import from tempfile
	mcon = sqlite3.connect(":memory:")
	mcon.cursor().executescript(tempfile.read())
	mcon.commit()
	mcon.row_factory = sqlite3.Row
	return mcon

def save_sqlite_db(conn, db_file):
	# Read database to tempfile
	tempfile = StringIO()
	for line in conn.iterdump():
		tempfile.write('%s\n' % line)
	tempfile.seek(0)
	
	# Create a database in memory and import from tempfile
	ocon = sqlite3.connect(db_file)
	ocon.cursor().executescript(tempfile.read())
	ocon.commit()
	ocon.row_factory = sqlite3.Row

	conn.close()
	ocon.close()

def timeSpend( start ):
	done = time.time()
	elapsed = done - start
	return time.strftime( "%H:%M:%S", time.gmtime(elapsed) )

def print_message(msg, verbose, start, logfile, errorout=False):
	message = "[%s] %s\n" % (timeSpend(start), msg)
	#loging
	logfile.write( message )
	
	if errorout:
		sys.exit( message )
	elif verbose:
		sys.stderr.write(message)

if __name__ == '__main__':
	argvs    = parse_params()

	cnt         = 0
	cnt_fa      = 0 # accepted genomes
	cnt_fa_tol  = 0 # total genomes
	cnt_seq     = 0
	cnt_seq_tol = 0
	tid2lineage = t._autoVivification()
	output      = sys.stdout
	output_gdbl_buffer = {}
	begin_t     = time.time()
	logfile_fn  = argvs.prefix+".gottcha_db_prep.log"
	logfile     = open(logfile_fn, "w")

	# loading taxonomy
	print_message( "Loading taxonomy..." , True, begin_t, logfile )
	t.loadTaxonomy( argvs.dbPath )
	print_message( "completed." , True, begin_t, logfile )

	# init the taxonomy sqlite3 db file
	print_message( "Initiating sqlite DB..." , True, begin_t, logfile )
	conn = init_sqlite_db(argvs.sqlitedb)
	print_message( "completed." , True, begin_t, logfile )
	c = conn.cursor()
	initTaxaDB(c)

	# input from assembly summary refseq report
	if not argvs.assemblyReport:
		print_message( "Downloading Refseq assembly report...", True, begin_t, logfile )
		url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"
		wget( url )
		argvs.assemblyReport = "assembly_summary_refseq.txt"

	# input from assembly summary refseq report
	if argvs.assemblyReport and os.path.isfile(argvs.assemblyReport):
		print_message( "Parsing assembly report...", True, begin_t, logfile )
		(asr_info, tid2lineage, cnt_fa, cnt_fa_tol) = processASR( c, 
																tid2lineage, 
																argvs.assemblyReport, 
																argvs.refseqCategoryOnly, 
																argvs.assemblyLevelOnly, 
																argvs.assemblyRefseqDir, 
																argvs.taxid, 
																argvs.excludeTaxid, 
																argvs.skipCheckingFiles, 
																argvs.redundantStrain, 
																argvs.includeIsolate)

		if not argvs.skipCheckingFiles:
			print_message( "Checking genome sequences...", True, begin_t, logfile )
			prepareGenomes(asr_info, argvs.cpu)
			print_message( "Done. %s/%s genome(s) checked."%(cnt_fa, cnt_fa_tol), True, begin_t, logfile )

		if not argvs.noGottchaDbList:
			print_message( "Generating gottcha_db input lists..." , True, begin_t, logfile )
			output_gdbl_buffer = generateGottchaDbList(asr_info, output_gdbl_buffer, tid2lineage, argvs.prefix)
			print_message( "Done.", True, begin_t, logfile )

	print_message( "Generating custom taxonomies tsv file..." , True, begin_t, logfile )
	# output custom database from sqlite3 to tsv
	with open(argvs.sqlitedb+".tsv", "w") as tsv_output:
		try:
			c.execute( 'SELECT * FROM taxonomy_custom ORDER BY NAME' )
			for row in c.fetchall():
				tsv_output.write( "%s\t%s\t%s\t%s\t%s\n"%(row[0], row[4], row[1], row[5], row[3]) )
		except sqlite3.Error as e:
			print_message( "[ERROR] An error occurred: %s"%e.args[0] , argvs.verbose, begin_t, logfile )
	tsv_output.close()
	print_message( "Done." , True, begin_t, logfile )

	print_message( "Generating strains list..." , True, begin_t, logfile )
	# output custom database from sqlite3 to tsv
	with open(argvs.sqlitedb+".all_strains.tsv", "w") as tsv_output:
		try:
			c.execute( 'SELECT * FROM strain_taxonomy ORDER BY NAME' )
			for row in c.fetchall():
				tsv_output.write( "%s\t%s\t%s\t%s\t%s\n"%(row[0], row[1], row[2], row[4], row[5]) )
		except sqlite3.Error as e:
			print_message( "[ERROR] An error occurred: %s"%e.args[0], argvs.verbose, begin_t, logfile, True )
	tsv_output.close()
	print_message( "Done." , True, begin_t, logfile )

	print_message( "Saving sqlite DB..." , True, begin_t, logfile )
	save_sqlite_db(conn, argvs.sqlitedb)
	print_message( "Done." , True, begin_t, logfile )

	print_message( "All completed!" , True, begin_t, logfile )

	logfile.close()
