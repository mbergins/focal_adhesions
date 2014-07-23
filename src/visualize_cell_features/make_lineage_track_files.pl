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

my $centroid_x_file = "$cfg{exp_results_folder}/adhesion_props/lin_time_series/Centroid_x.csv";
die "Couldn't find $centroid_x_file" if not -e $centroid_x_file;
my @centroid_x = $parser->read_file($centroid_x_file);

my $centroid_y_file = "$cfg{exp_results_folder}/adhesion_props/lin_time_series/Centroid_y.csv";
die "Couldn't find $centroid_y_file" if not -e $centroid_y_file;
my @centroid_y = $parser->read_file($centroid_y_file);

my @tracking_files = <$cfg{exp_results_folder}/visual*/lineage_id/*.png>;
my @lineage_colors_file = <$cfg{exp_results_folder}/visual*/lineage_id/*.csv>;
my @lineage_colors = $parser->read_file($lineage_colors_file[0]);

#

my @lineage_colors_hex = &convert_colors_to_hex(@lineage_colors);

for my $file_num (0..$#tracking_files) {
	my $target_file = $tracking_files[$file_num];
	$target_file =~ s/\.png$/\.svg/;
	
	my $source_file = basename($tracking_files[$file_num]);

	open OUTPUT, ">$target_file" or die $!;
	select OUTPUT;
	&output_header;

	print "<g style=\"opacity:0.75\"><image xlink:href=\"$source_file\" /></g>\n";
	
	for (0..$#centroid_x) {
		next if $centroid_x[$_][$file_num] eq "NaN";

		my $ID_start_zero = $_ + 1;

		print "<text y=\"$centroid_x[$_][$file_num]\" x=\"$centroid_y[$_][$file_num]\" style=\"font-size:1px;\">$ID_start_zero</text>\n";
	}

	print "</svg>";

	close OUTPUT;
}

################################################################################
# Functions
################################################################################

sub convert_colors_to_hex {
	my @hex;
	for (0..$#lineage_colors) {
		push @hex, sprintf("%02x%02x%02x", $lineage_colors[$_][0]*255, 
						  $lineage_colors[$_][1]*255, 
						  $lineage_colors[$_][2]*255);
	}
	return(@hex);
}

sub output_header {
	print "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
	print "\n";
	print "<svg xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";
}


