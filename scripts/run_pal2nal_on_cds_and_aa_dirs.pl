#!/usr/bin/perl

use strict;
use warnings; 
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

our $VERSION = '0.02';
our $AUTHOR  = 'Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>';

MAIN: {
    my $rh_o = get_options();
    my $aa   = $rh_o->{'aa_dir'};
    my $cds  = $rh_o->{'cds_dir'};
    my $ra_c = get_files($cds);
    my $ra_a = get_files($aa);
    my $out  = $rh_o->{'outdir'};

    unless (-d $out) {
        mkdir $out or die "cannot make $out:$!";
    }

    for (my $i = 0; $i < @{$ra_c}; $i++) {
        print "pal2nal.pl $aa/$ra_a->[$i] $cds/$ra_c->[$i] -output paml -nomismatch -nogap -codontable 1 > $out/$ra_c->[$i]_align.phy 2> $out/$ra_c->[$i].paml.stderr\n";
        system "pal2nal.pl $aa/$ra_a->[$i] $cds/$ra_c->[$i] -output paml -nomismatch -nogap -codontable 1 > $out/$ra_c->[$i]_align.phy 2> $out/$ra_c->[$i].paml.stderr\n";
        print "pal2nal.pl $aa/$ra_a->[$i] $cds/$ra_c->[$i] -output fasta -nomismatch -nogap -codontable 1 > $out/$ra_c->[$i]_align.fa 2> $out/$ra_c->[$i].fasta.stderr\n";
        system "pal2nal.pl $aa/$ra_a->[$i] $cds/$ra_c->[$i] -output fasta -nomismatch -nogap -codontable 1 > $out/$ra_c->[$i]_align.fa 2> $out/$ra_c->[$i].fasta.stderr\n";
    }
}

sub usage {
    print "usage: run_pal2nal_on_cds_and_aa_dirs.pl --cds_dir=CDS_DIR --aa_dir=AA_DIR --outdir=OUTDIR [--help] [--version]\n";
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
                                 'version'   => \$rh_opts->{'version'},
                                 'v'         => \$rh_opts->{'version'},
                                 'cds_dir=s' => \$rh_opts->{'cds_dir'},
                                 'aa_dir=s'  => \$rh_opts->{'aa_dir'},
                                 'outdir=s'  => \$rh_opts->{'outdir'},
                                 "h"         => \$rh_opts->{'help'},
                                 "help"      => \$rh_opts->{'help'});

    pod2usage({-exitval => 0, -verbose => 2}) if($rh_opts->{'help'});
    die "run_pal2nal_on_cds_and_aa_dirs.pl version $VERSION\n" if ($rh_opts->{'version'});
    print "missing --cds_dir\n" unless ($rh_opts->{'cds_dir'});
    print "missing --aa_dir\n" unless ($rh_opts->{'aa_dir'});
    print "missing --outdir\n" unless ($rh_opts->{'outdir'});
    usage() unless ($rh_opts->{'cds_dir'});
    usage() unless ($rh_opts->{'aa_dir'});
    usage() unless ($rh_opts->{'outdir'});
    return $rh_opts;
}

__END__

=head1 NAME

B<run_pal2nal_on_cds_and_aa_dirs.pl> - utility script for running pal2nal

=head1 AUTHOR

Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>

=head1 SYNOPSIS

run_pal2nal_on_cds_and_aa_dirs.pl --cds_dir=CDS_DIR --aa_dir=AA_DIR --outdir=OUTDIR [--help] [--version]

=head1 OPTIONS

=over

=item B<--cds_dir>

directory with CDS sequences (e.g. from TransDecoder)

=item B<--aa_dir>

directory with pruned (e.g. by PhyloPyPruner) amino acid sequences

=item B<--out_dir>

directory where output will be written. If this directory exists, it will be renamed with a timestamp at the end and a new directory will be created.

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


