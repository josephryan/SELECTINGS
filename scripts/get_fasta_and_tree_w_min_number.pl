#!/usr/bin/perl

use strict;
use warnings;
use File::Copy;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

our $VERSION = '0.02';
our $AUTHOR  = 'Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>';

MAIN: {
    my $rh_o    = process_options();
    my $fa_dir  = $rh_o->{'fa_dir'};
    my $treedir = $rh_o->{'tree_dir'};
    my $outdir  = $rh_o->{'out_dir'};
    my $min     = $rh_o->{'min_taxa'};
    my $count   = 0;

    check_outdir($outdir);

    opendir DIR, $fa_dir or die "cannot open $fa_dir:$!";
    my @files = readdir DIR;
    foreach my $f (@files) {
        next if ($f =~ m/SpeciesTreeAlignment.fa/);
        open IN, "$fa_dir/$f" or die "cannot open $fa_dir/$f:$!";
        my $count = 0;
        my $seqs = '';
        my %species = ();
        while (my $line = <IN>) {
            $seqs .= $line;
            next unless ($line =~ m/^>([^|]+)/);
            my $sp = $1;
            $count++ unless ($species{$sp});
            $species{$sp}++;
        }
        if ($count >= $min) {
            write_seqs($outdir,$f,$seqs);
            write_trees($treedir,$outdir,$f);
        }
    }
}

sub process_options {
    my $rh_opts = {};
    my $res = Getopt::Long::GetOptions(
                               "version"    => \$rh_opts->{'version'},
                               "v"          => \$rh_opts->{'version'},
                               "h"          => \$rh_opts->{'help'},
                               "help"       => \$rh_opts->{'help'},
                               "fa_dir=s"   => \$rh_opts->{'fa_dir'},
                               "fadir=s"    => \$rh_opts->{'fa_dir'},
                               "tree_dir=s" => \$rh_opts->{'tree_dir'},
                               "treedir=s"  => \$rh_opts->{'tree_dir'},
                               "out_dir=s"  => \$rh_opts->{'out_dir'},
                               "outdir=s"   => \$rh_opts->{'out_dir'},
                               "min_taxa=s" => \$rh_opts->{'min_taxa'},
                               "mintaxa=s"  => \$rh_opts->{'min_taxa'},
                                     );
    pod2usage({-exitval => 0, -verbose => 2}) if($rh_opts->{'help'});
    die "get_fasta_and_tree_w_min_number.pl version $VERSION\n" if ($rh_opts->{'version'});
    unless ($rh_opts->{'fa_dir'} && $rh_opts->{'tree_dir'} 
        && $rh_opts->{'out_dir'} && $rh_opts->{'min_taxa'}) {
        warn "--fa_dir is required\n" unless ($rh_opts->{'fa_dir'});
        warn "--tree_dir is required\n" unless ($rh_opts->{'tree_dir'});
        warn "--out_dir is required\n" unless ($rh_opts->{'out_dir'});
        warn "--min_taxa is required\n" unless ($rh_opts->{'min_taxa'});
        usage();
    } 
    return $rh_opts;
}

sub write_trees {
    my $tdir   = shift;
    my $outdir = shift;
    my $fasta  = shift;

    $fasta =~ m/^([^\/]+).fa$/ or die "unexpected format of fasta: $fasta";
    my $id = $1;
    File::Copy::copy("$tdir/${id}_tree.txt","$outdir/$id.tree");
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

sub write_seqs {
    my $dir = shift;
    my $file = shift;
    my $seqs = shift;
    open OUT, ">$dir/$file" or die "cannot open >$dir/$file:$!";
    print OUT $seqs;
    close OUT; 
}

sub usage {
    print "usage: get_fasta_and_tree_w_min_number.pl --fa_dir=FASTA_DIR --tree_dir=TREE_DIR --out_dir=OUT_DIR --min_taxa=MINIMUM_SEQS\n";    
    exit;
}

__END__

=head1 NAME

B<get_fasta_and_tree_w_min_number.pl> - get fasta and tree w min!

=head1 AUTHOR

Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>

=head1 SYNOPSIS

get_fasta_and_tree_w_min_number.pl --fa_dir=FASTA_DIR --tree_dir=TREE_DIR --out_dir=OUT_DIR --min_taxa=MINIMUM_SEQS [--help] [--version]

=head1 OPTIONS

=over

=item B<--fa_dir>

MultipleSequenceAlignments directory from an OrthoFinder run. Usually, OrthoFinder/Results_MonDD/MultipleSequenceAlignments where MonDD is 3 letter code for month followed by the date (e.g. Jul24).

=item B<--tree_dir>

Gene_Trees directory from an OrthoFinder run. Usually, OrthoFinder/Results_MonDD/Gene_Trees where MonDD is 3 letter code for month followed by the date (e.g. Jul24).

=item B<--out_dir>

directory where output will be written. If this directory exists, it will be renamed with a timestamp at the end and a new directory will be created.

=item B<--min_taxa>

minimum number of taxa represented in a sequence alignment in order for that alignment to be processed

=item B<--version>

Print the version and exit.

=item B<--help>

Print this manual.

=back

=head1 DESCRIPTION

This program is a helper script that we developed as part of the Polar Genomics Workshop. It automates a set of tasks that are part of utilizing the output of OrthoFinder for a large-scale selection analysis.

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

