#!/usr/bin/perl

use strict;
use warnings;
use File::Copy;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

our $VERSION = '0.03';
our $AUTHOR  = 'Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>';

MAIN: {
    my $rh_o    = process_options();
    my $outdir  = $rh_o->{'out_dir'};
    my $min     = $rh_o->{'min_taxa'};
    my $ofout   = $rh_o->{'of_out'} || '';;
    my $ofrdir  = $rh_o->{'of_resdir'} || '';;
    my $count   = 0;
    my $tdir = '';
    my $fdir  = '';

    ($tdir,$fdir) = get_tree_and_fa_dir($ofout,$ofrdir);
    check_outdir($outdir);

    opendir DIR, $fdir or die "cannot open $fdir:$!";
    my @files = readdir DIR;
    foreach my $f (@files) {
        next if ($f =~ m/SpeciesTreeAlignment.fa/);
        open IN, "$fdir/$f" or die "cannot open $fdir/$f:$!";
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
            write_trees($tdir,$outdir,$f);
        }
    }
}

sub get_tree_and_fa_dir {
    my $of     = shift;
    my $ofrdir = shift;
    my $tdir   = '';
    my $fdir   = '';
    if ($of) {
        open IN, $of or die "cannot open $of:$!";
        my $flag = 0;
        while (my $line = <IN>) {
            chomp $line;
            if ($flag) {
                chomp $line;
                $line =~ s/^\s+//;
                $ofrdir = $line;
                $flag = 0;
            } elsif ($line =~ m/^Results:/) {
                $flag = 1;
            }
        }
        die "$of (--of_out option) should be standard output of OrthoFinder run\n" unless ($ofrdir);
    }
    $tdir = "$ofrdir/Gene_Trees";
    $fdir = "$ofrdir/MultipleSequenceAlignments";
    return ($tdir,$fdir);
}

sub process_options {
    my $rh_opts = {};
    my $res = Getopt::Long::GetOptions(
                               "version"    => \$rh_opts->{'version'},
                               "v"          => \$rh_opts->{'version'},
                               "h"          => \$rh_opts->{'help'},
                               "help"       => \$rh_opts->{'help'},
                               "of_out=s"   => \$rh_opts->{'of_out'},
                              "of_resdir=s" => \$rh_opts->{'of_resdir'},
                               "out_dir=s"  => \$rh_opts->{'out_dir'},
                               "outdir=s"   => \$rh_opts->{'out_dir'},
                               "min_taxa=s" => \$rh_opts->{'min_taxa'},
                               "mintaxa=s"  => \$rh_opts->{'min_taxa'},
                                     );
    pod2usage({-exitval => 0, -verbose => 2}) if($rh_opts->{'help'});
    die "get_fasta_and_tree_w_min_number.pl version $VERSION\n" if ($rh_opts->{'version'});
    unless ( ($rh_opts->{'of_out'} || $rh_opts->{'of_resdir'}) 
             && $rh_opts->{'out_dir'} && $rh_opts->{'min_taxa'}) {
        unless ($rh_opts->{'of_out'} || $rh_opts->{'of_resdir'}) {
            warn "Either --of_resdir or --of_out are required\n";
        }
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
    print "usage: get_fasta_and_tree_w_min_number.pl --out_dir=OUT_DIR --min_taxa=MINIMUM_SEQS {--of_out|--of_resdir} [--help] [--version]\n";    
    exit;
}

__END__

=head1 NAME

B<get_fasta_and_tree_w_min_number.pl> - get fasta and tree w min!

=head1 AUTHOR

Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>

=head1 SYNOPSIS

get_fasta_and_tree_w_min_number.pl --out_dir=OUT_DIR --min_taxa=MINIMUM_SEQS {--of_out|--of_resdir} [--help] [--version]

=head1 OPTIONS

=over

=item B<--of_resdir>

Results directory from an OrthoFinder run. Usually something like, OrthoFinder/Results_MonDD where MonDD is 3 letter code for month followed by the date (e.g. Jul24). --of_out can be used instead.

=item B<--of_out>

If you save the standard output of an OrthoFinder run, you can supply this as an option to this script. The script will then parse the Results directory from this file. --of_resdir can be used instead.

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

