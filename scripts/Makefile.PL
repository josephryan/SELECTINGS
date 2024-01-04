use 5.008004;
use ExtUtils::MakeMaker;

WriteMakefile(
  NAME             => 'tmppolar',
  VERSION_FROM     => 'version',
  AUTHOR           => 'Joseph Ryan',
  LICENSE          => 'lgpl_3_0',
  MIN_PERL_VERSION => '5.008004',

  EXE_FILES        => ['get_corresponding_cds.pl', 'get_fasta_and_tree_w_min_number.pl', 'remove_blank_seqs_and_fewer_than_n.pl', 'run_codeml.pl', 'run_hyphy.pl', 'run_pal2nal_on_cds_and_aa_dirs.pl', 'codeml_chisquare.pl', 'parse_absrel.pl', 'parse_busted.pl', 'selectoprune.pl'],
  PREREQ_PM        => {
    'strict'   => 0,
    'warnings' => 0,
  },
  (eval { ExtUtils::MakeMaker->VERSION(6.46) } ? (META_MERGE => {
      'meta-spec' => { version => 2 },
      resources => {
          repository => {
              type => 'git',
              url  => 'https://github.com/josephryan/tmppolar',
              web  => 'https://github.com/josephryan/tmppolar',
          },
      }})
   : ()
  ),
);
