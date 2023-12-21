#!/usr/bin/perl

$|++;

use strict;
use warnings; 
use File::Copy;
use File::Spec;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

our $VERSION  = '0.02';
our $AUTHOR  = 'Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>';

our $HYPHY    = 'hyphy';
our @PROGS    = qw(absrel busted meme relax);

MAIN: {
    my $rh_o     = get_options();
    my $aln_dir  = $rh_o->{'aln_dir'};
    my $out_dir  = $rh_o->{'out_dir'};
    my $ra_alns  = get_files($aln_dir);
    my $pre      = $rh_o->{'pre'};
    my $tree     = $rh_o->{'tree'};

    check_outdir($out_dir);

    my $rh_skip = get_seqs_to_skip($rh_o) if ($rh_o->{'require_num_seqs'});
    foreach my $prog (@PROGS) {
        mkdir "$out_dir/$pre.$prog" or die "cannot make $out_dir/$pre.$prog:$!";
        run_hyphy($ra_alns,$aln_dir,$out_dir,$prog,$pre,$tree);
    }
}

sub check_outdir {
    my $outdir = shift;
    $outdir =~ s/\/\s*$//;
    if (-d $outdir) {
        my $timestamp = time();
        my $newdir = $outdir . ".$timestamp";
        File::Copy::move($outdir,$newdir);
        warn "warning: $outdir exists and has been moved to $newdir\n";
    }
    mkdir $outdir or die "cannot open $outdir";
}

sub get_seqs_to_skip {
    my $rh_o     = shift;
    my %to_skip  = ();
    my $dir      = $rh_o->{'aln_dir'};
    my $num      = $rh_o->{'require_num_seqs'};
    my $aln_dir  = $rh_o->{'aln_dir'};
    my $ra_files = get_files($dir);

    foreach my $file (@{$ra_files}) {
        my $fa_num = 0;
        open IN, "$aln_dir/$file" or die "cannot open $aln_dir/$file:$!";
        while (my $line = <IN>) {
            $fa_num++ if ($line =~ m/^>/);
        }
        $to_skip{$file} = 1 unless ($fa_num == $num);
    }
    return \%to_skip;
}

sub run_hyphy {
    my $ra_a = shift;
    my $adir = shift;
    my $odir = shift;
    my $prog = shift;
    my $pre  = shift;
    my $tree = shift;

    my $branch_option = "--branches Foreground";
    $branch_option = "--test Foreground" if ($prog eq 'relax');
    foreach my $aln (@{$ra_a}) {
        my $cmd = "$HYPHY $prog --alignment $adir/$aln --tree $tree $branch_option --output $odir/$pre.$prog/$aln.$prog.jason > $odir/$pre.$prog/$aln.out 2> $odir/$pre.$prog/$aln.err";
        print "running: $cmd\n";
        system $cmd;
        print "\n";
    }
}

sub usage {
    print "usage: run_hyphy.pl --aln_dir=ALIGNMENT_DIR --out_dir=OUT_DIR --tree=TREEFILE --pre=PREFIX_FOR_OUTPUT {--absrel --busted --meme --relax} [--require_num_seqs=NUM_SEQS_IN_ALIGNMENT] [--help] [--version]\n";
    exit;
}

sub get_files {
    my $dir = shift;
    opendir DIR, $dir or die "cannot opendir $dir:$!";
    my @files = grep { !/^\./ } readdir DIR;
    return \@files;
}

sub get_options {
    my $rh_opts = {};
    my $opt_results = Getopt::Long::GetOptions(
                                   'version' => \$rh_opts->{'version'},
                                         'v' => \$rh_opts->{'version'},
                                      'help' => \$rh_opts->{'help'},
                                         'h' => \$rh_opts->{'help'},
                                 'aln_dir=s' => \$rh_opts->{'aln_dir'},
                                 'out_dir=s' => \$rh_opts->{'out_dir'},
                                    'absrel' => \$rh_opts->{'absrel'},
                                    'busted' => \$rh_opts->{'busted'},
                                      'meme' => \$rh_opts->{'meme'},
                                     'relax' => \$rh_opts->{'relax'},
                        'require_num_seqs=i' => \$rh_opts->{'require_num_seqs'},
                                    'tree=s' => \$rh_opts->{'tree'},
                                    'pre=s'  => \$rh_opts->{'pre'});

    pod2usage({-exitval => 0, -verbose => 2}) if($rh_opts->{'help'});
    if ($rh_opts->{'version'}) {
        print "run_hyphy.pl version $VERSION\n";
        exit;
    }
    print "missing --aln_dir\n" unless ($rh_opts->{'aln_dir'});
    print "missing --out_dir\n" unless ($rh_opts->{'out_dir'});
    print "missing --tree\n" unless ($rh_opts->{'tree'});
    print "missing --pre\n" unless ($rh_opts->{'pre'});
    usage() unless ($rh_opts->{'aln_dir'});
    usage() unless ($rh_opts->{'out_dir'});
    usage() unless ($rh_opts->{'tree'});
    usage() unless ($rh_opts->{'pre'});
    unless ($rh_opts->{'absrel'} || $rh_opts->{'busted'} || $rh_opts->{'meme'} || $rh_opts->{'relax'}) {
        warn "$0 requires at least 1 of --absrel, --busted, --meme, --relax\n"; 
        usage();
    }
    return $rh_opts;
}

__END__

=head1 NAME

B<run_hyphy.pl> - utility script to run many instances of hyphy

=head1 AUTHOR

Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>

=head1 SYNOPSIS

run_hyphy.pl --aln_dir=ALIGNMENT_DIR --out_dir=OUT_DIR --tree=TREEFILE --pre=PREFIX_FOR_OUTPUT {--absrel --busted --meme --relax} [--require_num_seqs=NUM_SEQS_IN_ALIGNMENT] [--help] [--version]

=head1 OPTIONS

=over

=item B<--aln_dir>

directory with FASTA formatted alignment files

=item B<--out_dir>

directory where output will be written. If this directory exists, it will be renamed with a timestamp at the end and a new directory will be created.

=item B<--tree>

Newick formatted species tree

=item B<--pre>

Prefix for outfiles

=item B<--absrel>

Run aBSREL

=item B<--busted>

Run BUSTED

=item B<--meme>

Run MEME

=item B<--relax>

Run RELAX

=item B<--require_num_seqs>

minimum number of taxa represented in a sequence alignment in order for that alignment to be processed. optional.

=item B<--help>

Print this manual.

=item B<--version>

Print the version and exit.

=back

=head1 DESCRIPTION

This program is a helper script that we developed as part of the Polar Genomics Workshop. It automates running multiple hyphy instances on a number of sequences

=head1 BUGS

Please report them to the author.

=head1 ACKNOWLEDGEMENT

This material is based upon work supported by the National Science Foundation under Grant Numbers (1935672 and 1935635) awarded to Joseph Ryan and Scott Santagata.

=head1 COPYRIGHT

Copyright (C) 2023 Joseph F. Ryan, Scott Santagata

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

