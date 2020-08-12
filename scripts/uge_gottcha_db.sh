#!/bin/bash
#$ -cwd
#$ -l h_vmem=5000G
#$ -pe openmpi 24
#$ -m abe
#$ -M po-e@lanl.gov
#$ -j y
#$ -V
#$ -q 517mem.q,1Tmem.q,503Gmem.q,all.q@n1005,all.q@n1006,all.q@n1007,all.q@n1008,all.q@n1009,all.q@n1010
#$ -S /bin/bash
#$ -N gottcha_db

set -xe

module load openmpi/gcc/64/2.1.3
export OMP_NUM_THREADS=8

NCPU=24
INDIR=$1
OUTDIR=$2
PREFIX=$3

mpirun \
    -np $NCPU \
    --oversubscribe \
    -x OMP_NUM_THREADS \
    --mca btl self,sm,tcp \
  gottcha_db \
    --strain  $INDIR/$PREFIX.strain.list \
    --species $INDIR/$PREFIX.species.list \
    --genus   $INDIR/$PREFIX.genus.list \
    --family  $INDIR/$PREFIX.family.list \
    --order   $INDIR/$PREFIX.order.list \
    --class   $INDIR/$PREFIX.class.list \
    --phylum  $INDIR/$PREFIX.phylum.list \
    --kingdom $INDIR/$PREFIX.superkingdom.list \
    --strain.prefix  "strain-" \
    --species.prefix "species-" \
    --genus.prefix   "genus-" \
    --family.prefix  "family-" \
    --order.prefix   "order-" \
    --class.prefix   "class-" \
    --phylum.prefix  "phylum-" \
    --kingdom.prefix "superkingdom-" \
    --log $PREFIX.gottcha_db.log \
    --root $OUTDIR \
    --squash SQUASH \
    --compress \
    -w 24 \
    -f 30 \
    --RAM 150GB \
    --verbose

set +xe
#  --strain <strain-level genome mapping file>
#  --species <species-level genome mapping file>
#  --genus <genus-level genome mapping file>
#  --family <family-level genome mapping file>
#  --order <order-level genome mapping file>
#  --class <class-level genome mapping file>
#  --phylum <phylum-level genome mapping file>
#  --kingdom <kingdom-level genome mapping file>
#
#  --strain.prefix <strain-level output file prefix>
#  --species.prefix <species-level output file prefix>
#  --genus.prefix <genus-level output file prefix>
#  --family.prefix <family-level output file prefix>
#  --order.prefix <order-level output file prefix>
#  --class.prefix <class-level output file prefix>
#  --phylum.prefix <phylum-level output file prefix>
#  --kingdom.prefix <kingdom-level output file prefix>
#          The taxonomic level-prefix that will be appended to the output fasta files.
#          These user defined filename prefixes are needed to make the output fasta files
#          distinct from both the input fasta files and the output fasta files at other
#          taxonomic levels.
#
#  [--squash <taxa label>] (prevent writing; can appear multiple times)
#  [-w <digestion word size in bp>] (default is 24)
#  [-f <min output fragment length in bp>] (default is 30)
#  [--compress (compress the output fasta files)]
