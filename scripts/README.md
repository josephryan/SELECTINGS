# Utility scripts used in Polar Genomics Workshop Pipeline

### Makefile.PL

script to generate Makefile

### codeml_chisquare.pl

run ChiSquare analyses on codeml outputs

### get_corresponding_cds.pl

Get CDS files that correspond with pruned AA files

### get_fasta_and_tree_w_min_number.pl

copy alignments and corresponding gene trees into a new folder

### remove_blank_seqs_and_fewer_than_n.pl

remove blank sequences from PhyloPyPruner output

### run_codeml.pl

run codeml (from PAML) on many sequence alignments

### run_hyphy.pl

run hyphy on many sequence alignments

### run_pal2nal_on_cds_and_aa_dirs.pl

run pal2nal on many files

### selectoprune.pl (requires https://github.com/josephryan/JFR-PerlModules)

prune species tree to match subset of represented species in CDS alignment

### version

used by Makefile.PL during install
