#!/usr/bin/env perl

################################################################################
# Global Variables and Modules
################################################################################

use lib "../lib";

use strict;
use File::Find;
use File::Path;
use File::Spec::Functions;
use File::Basename;
use File::Copy;
use Getopt::Long;
use Data::Dumper;

use Config::Adhesions qw(ParseConfig);

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l")
  or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

################################################################################
# Main Program
################################################################################

my @folders = sort <$cfg{individual_results_folder}/*>;

die "Unable to find image results folders" if scalar(@folders) == 0;

our @file_list;
find(\&collect_all, ($folders[0]));

my @first_file_set = @file_list;
my $first_file_count = scalar(@file_list);
foreach (@folders) {
    @file_list = ();
    find(\&collect_all, ($_));
	

	if (scalar(@file_list) != $first_file_count) {
		my $die_str = "\nProblem with Directory: $_\nFound " . scalar(@file_list) . " files in set, while $first_file_count in the first directory.\n\n";
		die $die_str;
	};
}

################################################################################
# Functions
################################################################################

sub collect_all {
    push @file_list, $File::Find::name;
}

sub diff_matrices {
	my @mat_1 = @{$_[0]};
	my @mat_2 = @{$_[1]};
	
	my %file_marks;

	if (scalar(@mat_1) > scalar(@mat_2)) {
		
	} else {
	
	}
}

################################################################################
#Documentation
################################################################################

=head1 NAME

check_file_complement.pl - check the count of files deposited in each individual
results directory

=head1 SYNOPSIS

check_file_complement.pl -cfg FA_config

=head1 Description

Running MATLAB jobs on emerald is hit or miss, as such, I need a short program to search for and determine if a full set of files have been deposited by the latest program run. This program simply counts the number of files in each image directory and fails with an error when the counts differ for any of the folders. Future versions may pass to image number back, this version simply fails.

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=back

Optional parameter(s):

=head1 EXAMPLES

check_file_complement.pl -cfg FA_config

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 7/24/2009 

=cut
