#!/usr/bin/perl

package main;

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Copy;
#use JFR::Fasta;
use Data::Dumper;

our $VERSION = '0.02';
our $AUTHOR  = 'Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>';

MAIN: {
    my $rh_o     = get_options();
    my $cds_dir  = $rh_o->{'cds_dir'};
    my $aa_dir   = $rh_o->{'aa_dir'};
    my $outdir   = $rh_o->{'out_dir'};
#    my $cds_dir = $ARGV[0] or die "usage: $0 CDS_DIR PRUNED_AA_DIR OUTDIR\n";
#    my $aa_dir  = $ARGV[1] or die "usage: $0 CDS_DIR PRUNED_AA_DIR OUTDIR\n";
#    my $outdir = $ARGV[2] or die "usage: $0 CDS_DIR PRUNED_AA_DIR OUTDIR\n";

    check_outdir($outdir);

    my $ra_cds = get_files($cds_dir,'cds');
    my $ra_aa  = get_files($aa_dir,'fa');
    my $rh_cds = get_cds_hash($cds_dir,$ra_cds);
    make_pruned_cds_files($rh_cds,$ra_aa,$aa_dir,$outdir);
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

sub make_pruned_cds_files {
    my $rh_cds = shift;
    my $ra_aa  = shift;
    my $aa_dir = shift;
    my $dir    = shift;
    foreach my $file (@{$ra_aa}) {
        $file =~ m/(.*)\.fa$/;
        my $outcds = "$dir/$1.cds.fa";
        open OUT, ">$outcds" or die "cannot open $outcds:$!";
        my $fp = JFR::Fasta->new("$aa_dir/$file");
        while (my $rec = $fp->get_record()) {
            if ($rh_cds->{$rec->{'def'}}) {
                print OUT "$rec->{'def'}\n$rh_cds->{$rec->{'def'}}\n";
            } else {
                die "cannot find $rec->{'def'}\nWARNING: $outcds is corrupt";
            }
        }
    }
}

sub get_cds_hash {
    my $dir    = shift;
    my $ra_cds = shift;
    my %cds    = ();
    foreach my $cds (@{$ra_cds}) {
        my $fp = JFR::Fasta->new("$dir/$cds");
        while (my $rec = $fp->get_record()) {
            $rec->{'def'} =~ s/:/_/g; 
            $rec->{'def'} =~ s/ .*//;
            $cds{$rec->{'def'}} = $rec->{'seq'};
        }
    }
    return \%cds;
}

sub get_files {
    my $dir = shift;
    my $suf = shift;
    opendir DIR, $dir or die "cannot opendir $dir:$!";
    my @files = grep { /\.$suf/ } readdir DIR;
    return \@files;
}

sub get_options {
    my $rh_opts = {};
    my $opt_results = Getopt::Long::GetOptions(
                                 'version'   => \$rh_opts->{'version'},
                                 'v'         => \$rh_opts->{'version'},
                                 'help'      => \$rh_opts->{'help'},
                                 'h'         => \$rh_opts->{'help'},
                                 'cds_dir=s' => \$rh_opts->{'cds_dir'},
                                 'cdsdir=s'  => \$rh_opts->{'cds_dir'},
                                 'out_dir=s' => \$rh_opts->{'out_dir'},
                                 'outdir=s'  => \$rh_opts->{'out_dir'},
                                 'aadir=s'   => \$rh_opts->{'aa_dir'},
                                 'aa_dir=s'  => \$rh_opts->{'aa_dir'});

    pod2usage({-exitval => 0, -verbose => 2}) if($rh_opts->{'help'});
    die "get_corresponding_cds.pl version $VERSION\n" if ($rh_opts->{'version'});
    print "missing --cds_dir\n" unless ($rh_opts->{'cds_dir'});
    print "missing --aa_dir\n" unless ($rh_opts->{'aa_dir'});
    print "missing --out_dir\n" unless ($rh_opts->{'out_dir'});
    usage() unless ($rh_opts->{'cds_dir'});
    usage() unless ($rh_opts->{'aa_dir'});
    usage() unless ($rh_opts->{'out_dir'});
    return $rh_opts;
}

sub usage {
    print "get_corresponding_cds.pl --cds_dir=CDS_DIR --aa_dir=AA_DIR --out_dir=OUT_DIR [--help] [--version]\n";
    exit;
}

package JFR::Fasta;

use strict;
use FileHandle;
use Carp;
use IO::Uncompress::Gunzip qw($GunzipError);
use IO::Uncompress::Bunzip2 qw($Bunzip2Error);
use IO::Uncompress::Unzip qw($UnzipError);

# based off code by Ken Trout named NHGRI::FastaParser
$JFR::Fasta::AUTHOR  = 'Joseph Ryan';
$JFR::Fasta::VERSION = '1.4'; # match Git release version

sub get_def_w_o_gt {
    my $self = shift;
    my $def = shift;
    $def =~ m/^>(.*)/;
    return $1;
}

sub get_record{
    my $self = shift;
    my @record_lines = ();

    return '' if (! $self->{'filehandle'} );

    my %h;
    my $seq  = '';
    my $line = '';

    my $fh = $ { $self->{'filehandle'} };

    while( $line = $fh->getline() ) {
	chomp $line;
	if ($line =~ m/^(>.+)/) {
            if ($seq) {
		$h{'def'}  = $self->{'stored_def'};
                $h{'seq'} = $seq;
                $self->{'stored_def'} = $1;   
                return \%h;
	    }
	    $self->{'stored_def'} = $1;   
	}   else {
            $seq .= $line;
        }
    }
    $h{'def'}  = $self->{'stored_def'};
    $h{'seq'} = $seq;
    undef $self->{'filehandle'};
    return \%h;
}

sub new {
    my $invocant = shift;
    my $file     = shift;
    my $fh       = {};

    if ($file =~ m/\.gz$/i) {
        $fh = IO::Uncompress::Gunzip->new($file)
            or croak "IO::Uncompress::Gunzip of $file failed: $GunzipError\n";
    } elsif ($file =~ m/\.bz2$/i) {
        $fh = IO::Uncompress::Bunzip2->new($file)
            or croak "IO::Uncompress::Bunzip2 of $file failed: $Bunzip2Error\n";
    } elsif ($file =~ m/\.zip$/i) {
        $fh = IO::Uncompress::Unzip->new($file)
            or croak "IO::Uncompress::Unzip of $file failed: $UnzipError\n";
    } else {
        $fh = FileHandle->new($file,'r');
        croak "cannot open $file: $!\n" unless(defined $fh);
    }
    my $class = ref($invocant) || $invocant;
    my $self = { 'filehandle'  => \$fh, 'stored_def' => '' };
    return bless $self, $class;
}

1;

__END__

=head1 NAME

B<get_corresponding_cds.pl> - utility script for selection analyses

=head1 AUTHOR

Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>

=head1 SYNOPSIS

get_corresponding_cds.pl --cds_dir=CDS_DIR --aa_dir=AA_DIR --out_dir=
OUT_DIR [--help] [--version]

=head1 OPTIONS

=over

=item B<--cds_dir>

directory of CDS sequences

=item B<--aa_dir>

directory of amino acid sequences

=item B<--out_dir>

directory where output will be written. If this directory exists, it will be r
enamed with a timestamp at the end and a new directory will be created.

=item B<--version>

Print the version and exit.

=item B<--help>

Print this manual.

=back

=head1 DESCRIPTION

This program is a helper script that we developed as part of the Polar Genomics Workshop. This usually involves getting CDS sequences from TransDecoder and removing sequences that were removed during PhyloPyPruner run.

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
