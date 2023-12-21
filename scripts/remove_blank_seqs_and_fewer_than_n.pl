#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Copy;

our $VERSION = '0.02';
our $AUTHOR  = 'Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>';

MAIN: {
    my $rh_o     = get_options();
    my $alndir  = $rh_o->{'aln_dir'};
    my $outdir  = $rh_o->{'out_dir'};
    my $minseq  = $rh_o->{'min_seq'};
    my $ra_alns  = get_files($alndir);

    check_outdir($outdir);
    
    foreach my $aln (@{$ra_alns}) {
        if (check_aln("$alndir/$aln",$minseq)) {
            File::Copy::copy("$alndir/$aln","$outdir/$aln")
                or die "copy failed: $alndir/$aln $outdir/$aln:$!";
        }
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

sub check_aln {
    my $file = shift;
    my $min  = shift;
    open IN, $file or die "cannot open $file:$!";
    my $count = 0;
    my $seq = '';
    while (my $line = <IN>) {
        if ($line =~ m/^>/) {
            return 0 if ($count && !$seq);
            $count++;
            $seq = '';
        } else {
            chomp $line;
            $line =~ s/\s+//g;
            $seq .= $line;
        }
    }
    return 0 unless ($count >= $min && $seq);
    return 1;
}

sub get_files {
    my $dir = shift;
    opendir DIR, $dir or die "cannot opendir $dir:$!";
    my @files = grep { !/^\./ } readdir DIR;
    return \@files;
}

sub usage {
    print "usage: remove_blank_seqs_and_fewer_than_n.pl --aln_dir=ALIGNMENT_DIR --out_dir=OUT_DIR --min_seq=MINIMUM_NUMBER_OF_SEQS [--help] [--version]\n";
    exit;
}

sub get_options {
    my $rh_opts = {};
    my $opt_results = Getopt::Long::GetOptions(
                                 'version'   => \$rh_opts->{'version'},
                                 'v'         => \$rh_opts->{'version'},
                                 'help'      => \$rh_opts->{'help'},
                                 'h'         => \$rh_opts->{'help'},
                                 'aln_dir=s' => \$rh_opts->{'aln_dir'},
                                 'out_dir=s' => \$rh_opts->{'out_dir'},
                                 'min_seq=s' => \$rh_opts->{'min_seq'});

    pod2usage({-exitval => 0, -verbose => 2}) if($rh_opts->{'help'});
    die "remove_blank_seqs_and_fewer_than_n.pl version $VERSION\n" if ($rh_opts->{'version'});
    print "missing --aln_dir\n" unless ($rh_opts->{'aln_dir'});
    print "missing --out_dir\n" unless ($rh_opts->{'out_dir'});
    print "missing --min_seq\n" unless ($rh_opts->{'min_seq'});
    usage() unless ($rh_opts->{'aln_dir'});
    usage() unless ($rh_opts->{'out_dir'});
    usage() unless ($rh_opts->{'min_seq'});
    return $rh_opts;
}

__END__

=head1 NAME

B<remove_blank_seqs_and_fewer_than_n.pl> - utility script for phylopypruner

=head1 AUTHOR

Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>

=head1 SYNOPSIS

remove_blank_seqs_and_fewer_than_n.pl --aln_dir=ALIGNMENT_DIR --out_dir=OUT_DIR --min_seq=MINIMUM_NUMBER_OF_SEQS [--help] [--version]

=head1 OPTIONS

=over

=item B<--aln_dir>

directory of alignments in FASTA format. 

=item B<--out_dir>

directory where output will be written. If this directory exists, it will be r
enamed with a timestamp at the end and a new directory will be created.

=item B<--min_taxa>

minimum number of taxa represented in a sequence alignment in order for that alignment to be processed.

=item B<--version>

Print the version and exit.

=item B<--help>

Print this manual.

=back

=head1 DESCRIPTION

This program is a helper script that we developed as part of the Polar Genomics Workshop. PhyloPyPruner alignments sometimes have no residues and fewer than the expected number of sequences (although this may have been fixed in more recent versions). This program checks that and only copies alignments that pass the checks to the output directory.

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

