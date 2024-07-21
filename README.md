# SELECTINGS
A pipeline to detect positive Darwinian selection in large datasets

The detection of positive Darwinian selection enables the identification of important changes in protein-coding regions of a genome that helped shape the evolutionary history of a particular lineage. Despite the availability of coding sequence data for many thousands of species, it remains a challenge to conduct large-scale selection analyses on broad collections of sequence data. To address this gap, we have developed a pipeline that allows for the analysis of an unlimited number of genes from an unlimited number of assembled transcriptomes and/or gene models gathered from genomic datasets. The SELECTINGS (Screening Evolutionary Lineages for Exceptional Coding Transcripts In Next Generation Sequence) Pipeline incorporates several existing bioinformatic packages for translating transcriptomes (Transdecoder), identifying orthologous sequences (OrthoFinder), aligning sequences (MAFFT), building gene trees (FastTree or IQTree), pruning trees (PhyloPyPruner), converting protein alignments to nucleotide alignments (PAL2NAL), detecting gene and site-based evidence of positive selection (comparable branch-site tests in PAML and HyPhy). This repo includes a set of scripts, documentation, sample data, installation instructions , and a vignette. Compared to other similar tools, the SELECTINGS Pipeline allows for greater flexibility in terms of methodology, number of genes that can be analyzed, and the phylogenetic breadth that can be sampled.

## Getting started with SELECTINGS
You need to install the following prerequesites. We will create 2 conda images (~2.3Gb disk space required). This is necessary as hyphy does not work when some of these other packages are installed.

1. Create hyphy environment

```bash
conda create --name hyphy
conda activate hyphy
conda install -c bioconda hyphy
conda deactivate
```

2. Create general environment

```bash
conda create --name selectings
conda activate selectings

# note: the following commands will take a while
conda install -y -c conda-forge perl perl-uri r-ape scipy
conda install -y -c bioconda orthofinder transdecoder pal2nal paml perl-db-file perl-math-cdf perl-json-parse perl-set-intervaltree perl-uri
pip install phylopypruner   

3. Clone SELECTINGS (~1.5Gb disk space required) and install scripts/modules

```bash
git clone https://github.com/josephryan/SELECTINGS
cd SELECTINGS/scripts
perl Makefile.PL
make
make install
```

4. Uncompress sample data

```bash
cd ..
gzip -d sample_data/*.gz 
```

5. Install SwissProt Database

```bash
lwp-download ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gzip -d uniprot_sprot.fasta.gz
diamond makedb --in uniprot_sprot.fasta -d swissprot
```

## Vignette 1

#### To run your own data make a directory at the same level as sample_data and make sure your deflines are simple like the ones in our sample_data.

1. Create vignette1 directory and change to this directory

```bash
mkdir vignette1
cd vignette1
```

##### NOTE: Many of the commands below can run much faster if the number of threads/CPUs is increased from 1.  For commands where threads/CPUs can be adjusted there is a note in brackets that starts with "THREADS: ". 
For example, [THREADS: adjust --threads]

2. Identify longest orfs (~2 minutes per FASTA)

```bash
TransDecoder.LongOrfs -t ../sample_data/Species_A.fasta > tdA.out 2> tdA.err
TransDecoder.LongOrfs -t ../sample_data/Species_B.fasta > tdB.out 2> tdB.err
TransDecoder.LongOrfs -t ../sample_data/Species_C.fasta > tdC.out 2> tdC.err
TransDecoder.LongOrfs -t ../sample_data/Species_D.fasta > tdD.out 2> tdD.err
TransDecoder.LongOrfs -t ../sample_data/Species_E.fasta > tdE.out 2> tdE.err
TransDecoder.LongOrfs -t ../sample_data/Species_F.fasta > tdF.out 2> tdF.err
TransDecoder.LongOrfs -t ../sample_data/Species_G.fasta > tdG.out 2> tdG.err
```

3. Run diamond blastp vs. swissprot on each longest_orfs file [THREADS: adjust --threads]

```bash
diamond blastp --threads 1 -e 1e-5 -d ../swissprot -q Species_A.fasta.transdecoder_dir/longest_orfs.pep > Sp_A.dmd.out 2> Sp_A.dmd.err 
diamond blastp --threads 1 -e 1e-5 -d ../swissprot -q Species_B.fasta.transdecoder_dir/longest_orfs.pep > Sp_B.dmd.out 2> Sp_B.dmd.err 
diamond blastp --threads 1 -e 1e-5 -d ../swissprot -q Species_C.fasta.transdecoder_dir/longest_orfs.pep > Sp_C.dmd.out 2> Sp_C.dmd.err 
diamond blastp --threads 1 -e 1e-5 -d ../swissprot -q Species_D.fasta.transdecoder_dir/longest_orfs.pep > Sp_D.dmd.out 2> Sp_D.dmd.err 
diamond blastp --threads 1 -e 1e-5 -d ../swissprot -q Species_E.fasta.transdecoder_dir/longest_orfs.pep > Sp_E.dmd.out 2> Sp_E.dmd.err 
diamond blastp --threads 1 -e 1e-5 -d ../swissprot -q Species_F.fasta.transdecoder_dir/longest_orfs.pep > Sp_F.dmd.out 2> Sp_F.dmd.err 
diamond blastp --threads 1 -e 1e-5 -d ../swissprot -q Species_G.fasta.transdecoder_dir/longest_orfs.pep > Sp_G.dmd.out 2> Sp_G.dmd.err 
```

4. Predict likely coding regions [THREADS: adjust --cpu]

```bash
TransDecoder.Predict -t ../sample_data/Species_A.fasta --retain_blastp_hits Sp_A.dmd.out --cpu 1 > Sp_A.td.p.out 2> Sp_A.td.p.err
TransDecoder.Predict -t ../sample_data/Species_B.fasta --retain_blastp_hits Sp_B.dmd.out --cpu 1 > Sp_B.td.p.out 2> Sp_B.td.p.err
TransDecoder.Predict -t ../sample_data/Species_C.fasta --retain_blastp_hits Sp_C.dmd.out --cpu 1 > Sp_C.td.p.out 2> Sp_C.td.p.err
TransDecoder.Predict -t ../sample_data/Species_D.fasta --retain_blastp_hits Sp_D.dmd.out --cpu 1 > Sp_D.td.p.out 2> Sp_D.td.p.err
TransDecoder.Predict -t ../sample_data/Species_E.fasta --retain_blastp_hits Sp_E.dmd.out --cpu 1 > Sp_E.td.p.out 2> Sp_E.td.p.err
TransDecoder.Predict -t ../sample_data/Species_F.fasta --retain_blastp_hits Sp_F.dmd.out --cpu 1 > Sp_F.td.p.out 2> Sp_F.td.p.err
TransDecoder.Predict -t ../sample_data/Species_G.fasta --retain_blastp_hits Sp_G.dmd.out --cpu 1 > Sp_G.td.p.out 2> Sp_G.td.p.err
```

5. Create directories of peps and cds files

```bash
mkdir peps cds
mv *.pep peps
mv *.cds cds
```

6. Run orthofinder (long-running process) [THREADS: adjust -t]

```bash
orthofinder -X -z -t 1 -f peps -M msa > of.out 2> of.err 
```

7. Get data associated with orthogroups containing all seven species [NOTE:adjust Results_MonDD in command]. 

```bash
get_fasta_and_tree_w_min_number.pl --out_dir=GFWMN.out --min_taxa=7 --fa_dir=peps/OrthoFinder/Results_MonDD/MultipleSequenceAlignments --tree_dir=peps/OrthoFinder/Results_MonDD/Gene_Trees
```

8. Run PhyloPyPruner (identify 1-to-1 orthologs) [THREADS: adjust --threads] [NOTE: pip installs phylopypruner as phylopypruner.1.1.6]

```bash
phylopypruner --threads 1 --output=ppp.out --dir GFWMN.out --mask longest --min-support 0.5 --min-taxa 7 --prune MI > phylopypruner.out 2> phylopypruner.err
```

(TMPNOTE: DO WE NEED TO RUN remove_blank_seqs_and_fewer_than_n.pl ???)

9. Get CDS files that with each of the pruned AA files

```bash
get_corresponding_cds.pl --cds_dir=cds --aa_dir=ppp.out/phylopypruner_output/output_alignments --outdir=cds.subset
```

10. Adjust names of sequences to only include species names

```bash
perl -pi -e 's/^>([^|]+)\|.*/>$1/' cds.subset/* ppp.out/phylopypruner_output/output_alignments/*
```

11. Run pal2nal on the sequences in the cds and aa directories:

```bash
run_pal2nal_on_cds_and_aa_dirs.pl --aa_dir=ppp.out/phylopypruner_output/output_alignments --cds_dir=cds.subset --outdir=p2n.out > runp2n.out 2> runp2n.err
```

12. Create an unrooted version of the species tree generated by orthonfinder [NOTE:adjust Results_MonDD in command]

```bash
Rscript ../scripts/unroot.R peps/OrthoFinder/Results_MonDD/Species_Tree/SpeciesTree_rooted.txt > unrooted.tree
```

13. Edit species name in unrooted.tree

```bash
perl -pi -e 's/.fasta.transdecoder//g' unrooted.tree
```

##### PAML

14. Annotate PAML tree ('unrooted.tree') by adding #1 after species name on foreground branches; branches A,B,C,D.

```bash
# Before
cat unrooted.tree

perl -pi -e 's/(Species_[ABCD])/$1#1/g' unrooted.tree

# After
cat unrooted.tree
```

15. Run CODEML (program within PAML that tests for selection; estimated time = 2+ hours?)

```bash
cd p2n.out
run_codeml.pl --tree=../unrooted.tree --null --alt --aln_suf=phy > rc.out 2> rc.err
cd ..
```

##### HYPHY

16. deactivate the selectings conda environment and activate the hyphy environment

```bash
conda deactivate
conda activate hyphy
```

17. Annotate HYPHY tree ('unrooted.tree') by adding {Foreground} after species name on foreground branches for tree file; branches A,B,C,D.

```bash
mkdir hyphy
cd hyphy
cp unrooted.tre .
perl -pi -e 's/(Species_[ABCD])#1/$1\{Foreground\}/g' unrooted.tree

# create symbolic links to FASTA alignments (only those with non-zero file sizes)
perl -e 'opendir(D, "../p2n.out"); foreach $ff (readdir(D)) { $gg = "../p2n.out/$ff"; symlink $gg, $ff if (-f $gg && -s $gg && $gg =~ /fa$/); }'
```

18. Run a single BUSTED

```bash
hyphy busted --alignment OGXXXXXXXXX_pruned.cds.fa_align.fa --tree unrooted.tree --branches Foreground --output OGXXXXXX.busted.json
```

19. Run a single aBSREL (adjust OGXXXXXXXXX to correspond with real file)

```bash
hyphy aBSREL --alignment OGXXXXXXXXX_pruned.cds.fa_align.fa --tree unrooted.tree --branches Foreground --output OGXXXXXXXXX.absrel.json
```

20. Run a single MEME (adjust OGXXXXXXXXX to correspond with real file)

```bash
hyphy meme --alignment OGXXXXXXXXX_pruned.cds.fa_align.fa --tree unrooted.tree --branches Foreground --output OGXXXXXXXXX.meme.json
```

21. Run a single RELAX (note: RELAX use --test instead of --branches to specify branches)  (adjust OGXXXXXXXXX to correspond with real file)

```bash
hyphy relax --alignment OGXXXXXXXXX_pruned.cds.fa_align.fa --tree unrooted.tree --test Foreground --output OGXXXXXXXXX.relax.json
```

22. Run Busted, Absrel, and Meme (meme=7 hrs; busted=13 hrs; absrel=9 hrs):
```bash

perl ../scripts/run_hyphy.pl --absrel --busted --meme --aln_dir=ALNDIR --out_dir=02-OUT --tree=unrooted.tree --pre=hyphy --require_num_seqs=7 &

```

# Parsing results

PAML ALT versus NULL models

Once complete you will have two CODEML MCL Results files for each of your CDS gene alignments (ALT versus NULL). Using the codeml_chisquare.pl you can generate p-values. The script calculates the cumulative probability of the chi-square distribution, given the degrees of freedom (DF = number of sequences) and the chi-square test statistic (`X`) which is: 2*(lnL1(ALT)-lnL0(NULL)). The p-value is then computed as 1 - chisqrprob(DF, X).

NOTE: --df (degrees of freedom) = the number of parameters in the alternative model minus the number of parameters in the null model).

NOTE: --max_pval limits the output to those results with p-value less than the value supplied

```bash
conda deactivate hyphy
conda activate selectings
codeml_chisquare.pl --codeml_dir=CODEMLDIR --alt_suf=alt.codeml --null_suf=null.codeml --df=1 --max_pval=0.05
