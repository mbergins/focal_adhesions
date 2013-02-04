#!/usr/bin/perl -w

use File::Basename;
use Config::General qw(ParseConfig);
use Getopt::Long;

my %opt;
GetOptions(\%opt, "debug|d") or die;

###############################################################################
# Main
###############################################################################

my $find_results = `find /home/mbergins/Documents/Projects/focal_adhesions/trunk/data/FAAS_*/*.cfg`;
my @cfgs = split("\n",$find_results);

my @target_zips = map {
	my $cfg_file = basename($_);
	$cfg_file =~ s/cfg/zip/;
	$cfg_file
} @cfgs;

my $search_location = "/var/log/apache2/access*";

my %hit_counts;
foreach (@target_zips) {
	my $search_results = `zgrep $_ $search_location;`;
	my @search_results = split("\n",$search_results);
	
	# if (scalar(@search_results) > 1) {
	# 		print "################### $_ ################\n";
	# 		print $search_results, "\n";
	# }
	
	if (scalar(@search_results) == 0) {
		$hit_counts{$_} = scalar(@search_results);
	}
}

my @hit_sort = sort {$hit_counts{$a} <=> $hit_counts{$b}} keys %hit_counts;

for (@hit_sort) {
	print "$_ => $hit_counts{$_}\n";
}

###############################################################################
# Functions
###############################################################################
