#!/usr/bin/perl

# requires https://github.com/josephryan/JFR-PerlModules

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use JFR::Fasta;
use Data::Dumper;

our $PROGRAM = 'selectoprune.pl';
our $VERSION = '0.01';

our $PHYUTILITY = 'phyutility';  # adjust to path to phyutility if not in path
our $VERBOSE    = 0;  # set to 1 to print phyutility command line

MAIN: {
    my $rh_o  = get_opts();
    my $rh_ss = get_subset($rh_o->{'subset_fa'});
    my $ra_pr = get_pruners($rh_ss,$rh_o->{'full_fa'},$rh_o->{'subset_fa'});
    prune($ra_pr,$rh_o->{'tree'},$rh_o->{'pre'});
}

sub prune {
    my $ra_p = shift;
    my $tree = shift;
    my $pre  = shift;
    my $outfile = "${pre}.pruned.tre";
    my $names = '';
    foreach my $name (@{$ra_p}) {
        $names .= "-names $name ";
    }
    my $cmd = "phyutility -pr $names -in $tree -out $outfile";
    print "$cmd\n" if ($VERBOSE);
    system $cmd;
}

sub get_pruners {
    my $rh_ss   = shift;
    my $file    = shift;
    my $subset  = shift;
    my @pruners = ();
    my $overlap = 0;
    my $fp   = JFR::Fasta->new($file);
    while (my $rec = $fp->get_record()) {
        my $id = JFR::Fasta->get_def_w_o_gt($rec->{'def'});
        if ($rh_ss->{$id}) {
            $overlap++;
        } else {
            push @pruners, $id; 
        }
    }
    if ($overlap == 0) {
        die "there is no overlap between $file and $subset. Check that deflines are exactly the same.";
    } elsif (scalar(@pruners) == 0) {
        warn "WARNING: all the taxa between $file and $subset are the same. There is nothing to prune. Writing file anyway in case this is part of a pipeline.";
    }
    return \@pruners;
}

sub get_subset {
    my $file = shift;
    my %ss   = ();
    my $fp   = JFR::Fasta->new($file);
    while (my $rec = $fp->get_record()) {
        my $id = JFR::Fasta->get_def_w_o_gt($rec->{'def'});
        die "$id occurs more than once in $file" if ($ss{$id});
        $ss{$id} = 1;
    } 
    return \%ss;
}

sub usage {
    print "$PROGRAM --full_fa=FULL_MATRIX_FASTA --tree=NEWICK_TREE --subset_fa=SUBSET_ALIGNMENT_FASTA --pre=PREFIX_FOR_OUTFILE\n";
    exit;
}

sub get_opts {
    my $rh_o = {};
    my $opt_results = Getopt::Long::GetOptions(
                                       'full_fa=s'    => \$rh_o->{'full_fa'},
                                       'tree=s'       => \$rh_o->{'tree'},
                                       'subset_fa=s'  => \$rh_o->{'subset_fa'},
                                       'pre=s'        => \$rh_o->{'pre'},
                                       'help'         => \$rh_o->{'help'},
                                       'version'      => \$rh_o->{'version'});
    die "$PROGRAM version $VERSION\n" if ($rh_o->{'version'});
    pod2usage({-exitval => 0, -verbose => 2}) if($rh_o->{'help'});
    print "missing --full_fa\n" unless ($rh_o->{'full_fa'});
    print "missing --tree\n" unless ($rh_o->{'tree'});
    print "missing --subset_fa\n" unless ($rh_o->{'subset_fa'});
    print "missing --pre\n" unless ($rh_o->{'pre'});
    usage() unless ($rh_o->{'full_fa'} && $rh_o->{'tree'} && $rh_o->{'subset_fa'} && $rh_o->{'pre'});

    return $rh_o;
}

__END__

=head1 NAME

B<selectoprune.pl> build a pruned treefile to represent a subset of sequences


