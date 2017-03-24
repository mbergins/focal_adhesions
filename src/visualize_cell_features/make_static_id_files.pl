#!/usr/bin/env perl

################################################################################
# Global Variables and Modules
################################################################################

use lib "../lib";

use strict;
use File::Path;
use File::Spec::Functions;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use POSIX;
use Text::CSV::Simple;

use Config::Adhesions;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

################################################################################
# Main Program
################################################################################

#required files and reading in data

my $parser = Text::CSV::Simple->new;

my @individual_folders = <$cfg{individual_results_folder}/*>;

for my $image_folder (@individual_folders) {
	my $centroid_file = "$image_folder/raw_data/Centroid.csv";
	warn "Couldn't find $centroid_file" if not -e $centroid_file;
	next if not -e $centroid_file;
	my @centroid = $parser->read_file($centroid_file);
	
	my $highlights_file = "$image_folder/highlights.png" ;
	
	my $target_file = "$image_folder/cell_id.svg";
	
	my $source_file = basename($highlights_file);

	open OUTPUT, ">$target_file" or die $!;
	select OUTPUT;
	&output_header;

	print "<g style=\"opacity:0.75\"><image xlink:href=\"$source_file\" /></g>\n";
	
	for (0..$#centroid) {
		next if $centroid[$_] eq "NaN";

		my $ID_start_zero = $_ + 1;

		print "<text y=\"$centroid[$_][1]\" x=\"$centroid[$_][0]\" style=\"font-size:1px;\">$ID_start_zero</text>\n";
	}

	print "</svg>";

	close OUTPUT;

}

################################################################################
# Functions
################################################################################

sub output_header {
	print "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
	print "\n";
	print "<svg xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";
}
