#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Math::CDF;
use Data::Dumper;

our $VERSION = '0.01';
our $AUTHOR  = 'Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>';

MAIN: {
    my $rh_o = process_options();    
    my $ra_alt = get_files($rh_o->{'codeml_dir'},$rh_o->{'alt_suf'});
    my $ra_nul = get_files($rh_o->{'codeml_dir'},$rh_o->{'null_suf'});
    check_files($ra_alt,$ra_nul,$rh_o->{'alt_suf'},$rh_o->{'null_suf'});
    my $rh_lnls = get_lnls($rh_o->{'codeml_dir'},$ra_alt,$ra_nul,$rh_o->{'null_suf'});
    calculate_chisquare($rh_lnls,$rh_o->{'df'},$rh_o);
}

sub calculate_chisquare {
    my $rh_l = shift;
    my $df   = shift;
    my $rh_o = shift;
    print "Alignment,Chi-square test statistic,P-value\n";
    foreach my $key (keys %{$rh_l}) {
        my $lnl_alt  = $rh_l->{$key}->[0];
        my $lnl_null = $rh_l->{$key}->[1];
        # Calculate the chi-square test statistic
        my $chi = 2 * ($lnl_alt - $lnl_null);

        my $p_value = 1 - Math::CDF::pchisq($chi, $df);
        next if ($rh_o->{'max_pval'} && $p_value > $rh_o->{'max_pval'});
        print "$key,$chi,$p_value\n";
    }
}

sub get_lnls {
    my $dir  = shift;
    my $ra_a = shift;
    my $ra_n = shift;
    my $nsuf = shift;
    my %lnl  = ();
   
    for (my $i = 0; $i < @{$ra_a}; $i++) {
        $ra_n->[$i] =~ m/^(.*).$nsuf/ or die "unexpected suffix: $ra_n->[$i]";
        my $root = $1;
        my $alnl = get_lnl_from_file("$dir/$ra_a->[$i]");
        my $nlnl = get_lnl_from_file("$dir/$ra_n->[$i]");
        $lnl{$root} = [$alnl,$nlnl] if ($alnl && $nlnl);
    }
    return \%lnl;
}

sub get_lnl_from_file {
    my $file = shift;
    open IN, $file or die "cannot open $file:$!";
    while (my $line = <IN>) {
        next unless ($line =~ m/lnL/);
        $line =~ s/\s*$//;
        my @ff = split /\s+/, $line;
        return $ff[-2];
    }
    return '';
}

sub check_files {
    my $ra_a = shift;
    my $ra_n = shift;
    my $asuf = shift;
    my $nsuf = shift;
    my $na = scalar(@{$ra_a});
    my $nn = scalar(@{$ra_n});
    die "diff number of alt ($na) and null ($nn) files" unless ($na == $nn);
    for (my $i = 0; $i < @{$ra_a}; $i++) {
        $ra_a->[$i] =~ m/^(.*).$asuf/ or die "unexpected suffix: $ra_a->[$i]";
        my $aroot = $1;
        $ra_n->[$i] =~ m/^(.*).$nsuf/ or die "unexpected suffix: $ra_n->[$i]";
        my $nroot = $1;
        unless ($aroot eq $nroot) {
            die "$ra_a->[$i] $ra_n->[$i] file name roots don't match\n";
        }
    }
}

sub get_files {
    my $dir = shift;
    my $suf = shift;
    opendir DIR, $dir or die "cannot opendir $dir:$!";
    my @files = grep { /$suf/ } readdir(DIR);
    return \@files;
}

sub process_options {
    my $rh_opts = {};
    my $res = Getopt::Long::GetOptions(
                               "version"      => \$rh_opts->{'version'},
                               "v"            => \$rh_opts->{'version'},
                               "h"            => \$rh_opts->{'help'},
                               "help"         => \$rh_opts->{'help'},
                               "codeml_dir=s" => \$rh_opts->{'codeml_dir'},
                               "df=s"         => \$rh_opts->{'df'},
                             "max_pval=f"     => \$rh_opts->{'max_pval'},
                               "alt_suf=s"    => \$rh_opts->{'alt_suf'},
                               "null_suf=s"   => \$rh_opts->{'null_suf'},
                                     );
    pod2usage({-exitval => 0, -verbose => 2}) if($rh_opts->{'help'});
    die "codeml_chisquare.pl version $VERSION\n" if ($rh_opts->{'version'});
    unless ($rh_opts->{'codeml_dir'} && $rh_opts->{'alt_suf'} 
        && $rh_opts->{'null_suf'} ) {
        warn "--codeml_dir is required\n" unless ($rh_opts->{'codeml_dir'});
        warn "--alt_suf is required\n" unless ($rh_opts->{'alt_suf'});
        warn "--null_suf is required\n" unless ($rh_opts->{'null_suf'});
        warn "--df is required\n" unless ($rh_opts->{'df'});
        usage();
    } 
    return $rh_opts;
}
sub usage {
    print "codeml_chisquare.pl --codeml_dir=CODEML_DIR --alt_suf=ALT_SUFFIX --null_suf=NULL_SUFFIX --df=DEGREES_OF_FREEDOM [--max_pval=MAX_PVAL] [--help] [--version]\n";
    exit;
}

__END__

=head1 NAME

B<codeml_chisquare.pl> - generate chisquare values from codeml outputs

=head1 AUTHOR

Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>

=head1 SYNOPSIS

codeml_chisquare.pl --codeml_dir=CODEML_DIR --alt_suf=ALT_SUFFIX --null_suf=NULL_SUFFIX --df=DEGREES_OF_FREEDOM [--max_pval=MAX_PVAL] [--help] [--version]

=head1 OPTIONS

=over

=item B<--codeml_dir>

directory with codeml outputs

=item B<--alt_suf>

suffix of files generated by codeml that include lnL vals for alternative models

=item B<--null_suf>

suffix of files generated by codeml that include lnL vals for null models

=item B<--df>

degrees of freedom (the difference in the number of free parameters between the two models). For example, M1a has 2 free parameters and M2a has 4 so df=2. The M7-M8 comparison also has 2 degrees of freedom. (see https://doi.org/10.1093/molbev/msad041)

=item B<--max_pval=MAX_PVAL>

only print records with pval <= MAX_PVAL (e.g. --max_pval=0.05)

=item B<--version>

Print the version and exit.

=item B<--help>

Print this manual.

=back

=head1 DESCRIPTION

This program is a helper script that we developed as part of the Polar Genomics Workshop. It automates the calculation of chisquare values from codeml runs

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
