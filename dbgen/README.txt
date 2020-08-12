*** How to run the GOTTCHA database tool (gottcha_db) **

July 24, 2018
J. D. Gans (jgans@lanl.gov)
Bioscience Division
Los Alamos National Laboratory

The GOTTCHA database program (gottcha_db) is software that reads in DNA sequences from any
number of oragnisms, splits those DNA sequences into short fragments (k-mers), *removes*
the fragments that are shared between different taxonomic groups, and then writes the now
group-specific (unique, "signature") sequences back to disk.

For information on how to build and install the gottcha_db program, please see the INSTALL.txt
file.

COMMAND LINE OPTIONS:

The gottcha_db program is a command line program that accepts the following arguments:

	--strain <strain-level genome mapping file>
	--species <species-level genome mapping file>
	--genus <genus-level genome mapping file>
	--family <family-level genome mapping file>
	--order <order-level genome mapping file>
	--class <class-level genome mapping file>
	--phylum <phylum-level genome mapping file>
	--kingdom <kingdom-level genome mapping file>
		Specify the (optionally compressed) ASCII text, genome mapping file for the given
		taxonomic level. The genome mapping files provide the locations and taxonomic labels 
		for all of the input sequences to gottcha_db at a given taxonomic level. 
		The format for these files is described in more detail below. 
		If specifying a compressed file, the compression format must be "gzip".
		
	--strain.prefix <strain-level output file prefix>
	--species.prefix <species-level output file prefix>
	--genus.prefix <genus-level output file prefix>
	--family.prefix <family-level output file prefix>
	--order.prefix <order-level output file prefix>
	--class.prefix <class-level output file prefix>
	--phylum.prefix <phylum-level output file prefix>
	--kingdom.prefix <kingdom-level output file prefix>
		The taxonomic level-prefix that will be appended to the output fasta files. 
		These user defined filename prefixes are needed to make the output fasta files 
		distinct from both the input fasta files and the output fasta files at other
		taxonomic levels.
		
	[--squash <taxa label>] (prevent writing; can appear multiple times)
	
		This optional argument *excludes* the specified taxonomic label (taxonomic
		labels are strings that are defined in the genome mapping files) from being
		included in the *output*. These genomes *are* still included in the uniqueness
		calculation (that is, k-mers that are shared with squashed genomes are still
		removed). This is useful if you want to identify target genomic signatures that are 
		not found in a large set of background genomes, but don't care about signatures
		that are unique to the background genomes themselves.
		
		**Please note that squashed genomes are excluded at the specified taxonomic level and
		all levels below**.
		
	[-w <digestion word size in bp>] (default is 24)
		
		All input genomes are split in DNA "words" (i.e. k-mers) of the specified size.
		Unique words are then identified by performing set different operations on the
		set of words that represent each taxonomic group.
		
	[-f <min output fragment length in bp>] (default is 30)
		
		Only write fragments that are greater than or equal to the specified length
		to the output file fasta file (this avoids having fasta files that are cluttered
		with many very short DNA sequences).
		
	[--compress (compress the output fasta files)]
		
		Save space by compressing the output fasta files during writing. Files are compressed
		using the gzip algorithm and can be uncompressed using the popular "gunzip" 
		file utility.
		
	[--linear (force all genomes to have linear topology)]
		
		Over-ride the topology information specified in the genome mapping file and force
		all genomes to have linear topology.
		
	[--circular (force all genomes to have circular topology)]
		
		Over-ride the topology information specified in the genome mapping file and force
		all genomes to have circular topology.
		
	[--silent|--informative|--verbose] (default is slient)
		
		Control the amount of information written to STDERR.
		
	[--log <log file>]
		
		Write progress information to the specified file. This is useful for tracking the
		(often) lengthly gottcha_db calculations as they progress on a cluster computer.
		
	[--root <root directory for output files>]
		
		Over-ride the path information provided in the genome mapping file. Output fasta files
		will be written to the specified directory, using the specified file prefix (i.e. see 
		the --prefix option above) and file name provided in the genome mapping file.
		
		To avoid writing thousands of files to a single directory, and to avoid the chance of
		accidentaly creating multiple files with the same file name, only the largest common
		path prefix (shared by all entries in the genome mapping files) are replaced by the
		specified root directory. In other words, the root directory simply redirects the output
		file system path to point to a different root directory (and preserves the subdirectories
		needed to maintain file uniqueness).
		
	[--RAM <max RAM per MPI rank>] (default is 80% of available)
	
		Running the gottcha_db program is very RAM intensive. To avoid running out of RAM, this argument
		limits the amount of RAM that will be used by each MPI rank. RAM amounts are in bytes by default, 
		but can be specified in MB (i.e. --RAM 100MB) or GB (i.e. --RAM 250 GB).

THE GENOME MAPPING FILE:

The genome mapping file is a four column, tab-delimited file that specifies one or more sequence file records. 
The genome mapping file does *not* use a header and may be (optionally) compressed in gzip-format. 
The fields in the genome mapping file are defined as follows:

	Column 1: Name of taxonomic group (string)
		
		The name of the taxonomic group. Input sequence files are grouped into taxonomic groups 
		by this name. Spaces are allowed.
		
	Column 2: Path to fasta or GBK (GenBank) file
		
		The fasta or GBK (GenBank) file that contains the sequences of interest. 
		This file may be (optionally) compressed with gzip-based compression.
		
	Column 3: Topology ("linear" or "circular")
		
		DNA molecules have either "linear" (i.e. human chromosomes) or "circular" (i.e. bacterial plasmids)
		topology. For gottcha_db, only difference between the two topologies is that the k-mers will be
		tiled across the junction between the beginning and end of circular sequences.
		
	Column 4: Taxonomic Id (integer [.integer])
		
		The taxonomic Id to be included in fasta file output of gottcha_db. ** This id is NOT used to define
		groups for gottcha_db **. (Groups are defined by the string-based name in the first column). The
		downstreams steps in the GOTTCHA algorithm require fasta file deflines that have a taxonomic Id.
		
		The taxonomic Id of the specified accession (or accessions). Note that gottcha_db defines a
		composite taxonomic Id. The first part is an unsigned integer (just like the NCBI taxonomy database) and the
		optional second part is a period, ".", followed by a second unsigned integer. This optional second integer is
		used to specify additional taxonomic granularity for strain-level signature calculations (since the 
		commonly used NCBI taxonomy does not provide a number taxonomic Id to separate different strains).

Here is a toy example of a genome mapping file

Yersinia similis	genomes/Yersinia_similis_582515/NZ_CP007231.1.gbk.gz	circular	367190
Yersinia sp. n.51	genomes/Yersinia_sp._n.51_493635/wgs.CBLF.1.gbff.gz	linear	1341639
Yersinia wautersii	genomes/Yersinia_wautersii_1319825/wgs.CVMG.1.gbff.gz	linear	1341643.1
Yersinia wautersii	genomes/Yersinia_wautersii_493435/wgs.CBLJ.1.gbff.gz	linear	1341643.2

In this example, three different taxonomic group names are provided (i.e. "Yersinia similis", "Yersinia sp. n.51" and
"Yersinia wautersii"), so gottcha_db will compute the signatures that are unique to each of these three groups. Since
the "Yersinia sp. n.51" and the "Yersinia wautersii" genomes are unfinished (and may be in multiple contigs), we specify
a "linear" topology (to avoid artificially circularizing a linear contig). Note that the taxonomic Ids (in the last 
column) for the "Yersinia wautersii" genomes have strain-specific suffixes (i.e. ".1" and ".2"). This taxonomic Id 
information will be written to the output fasta files (containing the unique signatures), but does not other-wise 
effect the operation of the gottcha_db program (which defines the taxonomic groups using the string-based name
provided in the first column).

	
RUNNING ON A CLUSTER:

While gottcha_db can run on a single processor workstation or laptop (for toy programs), it likes a *lot* of memory for
processing the tens of thousands of genomes currently stored in GenBank. Additional memory can be provided by 
running gottcha_db on more cluster nodes (since gottcha_db is essentially a custom-build distributed memory machine -- 
additional cluster nodes are needed for their RAM and only a little bit of their CPUs). The best performance is
usually obtained by running one MPI rank per CPU core ("real" cores only, hyperthreaded cores don't count!).

Running gottcha_db on a cluster requires the "mpirun" helper program (included with every version of MPI) that manages 
the interactions between the parallel program (i.e. gottcha_db) and all of the backend system stuff needed to run
in parallel on a cluster (i.e. network communications, start up, shut down, ...).

For example: "mpirun -np 15 ./gottcha_db <gottcha_db command line argument here>" would run gottcha_db on 15 nodes of your
cluster. Typically this command will be run as part of a batch script that is submitted to a cluster schedualing system
(i.e. MOAB, Torque, SLURM). The gottcha_db program is schedualer agnostic.
