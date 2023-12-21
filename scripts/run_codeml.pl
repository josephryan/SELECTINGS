#!/usr/bin/perl

$|++;

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

our $VERSION = '0.02';
our $AUTHOR  = 'Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>';

MAIN: {
    my $rh_o    = get_options();
    my $tree    = $rh_o->{'tree'};
    my $aln_suf = $rh_o->{'aln_suf'};

    opendir DIR, "." or die "cannot opendir .:$!";
    my @files = grep { /\.$aln_suf$/ } readdir DIR;

    foreach my $file (@files) {
        if ($rh_o->{'alt'}) {
            write_ctl_file(0,'alt',$file,$tree);
            print "running: codeml ${file}.alt.ctl > ${file}.alt.out 2> ${file}.alt.err...";
            system "codeml ${file}.alt.ctl > ${file}.alt.out 2> ${file}.alt.err";
            print "\n";
        } 
        if ($rh_o->{'null'}) {
            write_ctl_file(1,'null',$file,$tree);
            system "codeml ${file}.null.ctl > ${file}.null.out 2> ${file}.null.err";
        } 
    }
}

sub write_ctl_file {
    my $fix_omega = shift;
    my $type      = shift;
    my $file      = shift;
    my $tree      = shift;

    my $ctl_file  = "${file}.$type.ctl";
    open OUT, ">$ctl_file" or die "cannot open $ctl_file:$!";
    print OUT qq~seqfile = $file
treefile = $tree
outfile = ${file}.${type}.codeml
noisy = 4
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
estFreq = 0
ndata = 1
clock = 0
aaDist = 0
model = 2
NSsites = 2
icode = 0
Mgene = 0
fix_kappa = 0
kappa = 2
fix_omega = $fix_omega
omega = 1
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 5
getSE = 0
RateAncestor = 0
Small_Diff = 5e-7
cleandata = 0
fix_blength = 1
method = 0
~;

}

sub usage {
    print "usage: run_codeml.pl --tree=TREEFILE --aln_suf=SUFFIX_OF_ALIGNMENT_FILES {--null|--alt} [--version]\n";
    print "use --null to run null models and/or --alt to run alternative models\n";
    print "alignment files should be in the directory from which you are running the script.\n";
    exit;
}

sub get_options {
    my $rh_opts = {};
    my $opt_results = Getopt::Long::GetOptions(
                                 'version'   => \$rh_opts->{'version'},
                                 'v'         => \$rh_opts->{'version'},
                                 'tree=s'    => \$rh_opts->{'tree'},
                                 'aln_suf=s' => \$rh_opts->{'aln_suf'},
                                 'alt'       => \$rh_opts->{'alt'},
                                 'null'      => \$rh_opts->{'null'},
                                 "h"         => \$rh_opts->{'help'},
                                 "help"      => \$rh_opts->{'help'});
    pod2usage({-exitval => 0, -verbose => 2}) if($rh_opts->{'help'});
    die "run_codeml.pl version $VERSION\n" if ($rh_opts->{'version'});
    print "missing --tree\n" unless ($rh_opts->{'tree'});
    print "missing --aln_suf\n" unless ($rh_opts->{'aln_suf'});
    print "--alt and/or --null is required\n" unless ($rh_opts->{'alt'});
    usage() unless ($rh_opts->{'tree'});
    usage() unless ($rh_opts->{'aln_suf'});
    usage() unless ($rh_opts->{'alt'} || $rh_opts->{'null'});
    return $rh_opts;
}

__END__

=head1 NAME

B<run_codeml.pl> - utility script for running codeml

=head1 AUTHOR

Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>

=head1 SYNOPSIS

run_codeml.pl --tree=TREEFILE --aln_suf=SUFFIX_OF_ALIGNMENT_FILES {--null|--alt} [--version]

=head1 OPTIONS

=over

=item B<--tree>

Newick formatted species tree

=item B<--aln_suf>

Suffix of alignment files, which should be in the directory from which you are running the script.

=item B<--null>

use --null to run null models

=item B<--alt>

use --alt to run alternative models\n

tem B<--version>

Print the version and exit.

=item B<--help>

Print this manual.

=back

=head1 DESCRIPTION

This program is a helper script that we developed as part of the Polar Genomics Workshop. It automates a set of tasks necessary to run codeml on a number of sequences.

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

