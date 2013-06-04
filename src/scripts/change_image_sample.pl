#!/usr/bin/env perl

################################################################################
# Global Variables and Modules
################################################################################

use lib "../lib";

use strict;
use POSIX;
use File::Path;
use File::Spec::Functions;
use File::Basename;
use File::Copy;
use Image::ExifTool;
use Math::Matlab::Local;
use Getopt::Long;
use Data::Dumper;

use Config::Adhesions qw(ParseConfig);
use Image::Stack;
use Math::Matlab::Extra;
use Emerald;
use FA_job;

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

#check for the presence of a "/" at the end of the data folder config variable,
#if present remove it, we will need a "/"-free ending for building the new
#folder names
#
#Also note that I'm using a m^^ form of the regular expression matching, see
#page 115 in Learning Perl, 3rd ed for details
if ($cfg{data_folder} =~ m^(.*)/$^) {
	$cfg{data_folder} = $1;
}

my @image_sets = ([qw(raw_mask_folder raw_mask_file)], [qw(adhesion_image_folder adhesion_image_file)]);
for (@image_sets) {
	&move_target_image_set(@{$_});
}

################################################################################
#Functions
################################################################################

sub move_target_image_set {
	my $target_dir = shift @_;
	my $source_file = shift @_;
		
	if (not exists $cfg{$target_dir} || not exists $cfg{$source_file}) {
		die;
	}

	my @image_files = sort <$cfg{individual_results_folder}/*/$cfg{$source_file}>;
	my @image_dirs = sort <$cfg{individual_results_folder}/*>;
	my @image_nums = map basename($_), @image_dirs;
	
	if (scalar(@image_files) == 0 && scalar(@image_dirs) != 0) {
		next;
	}

	if (scalar(@image_files) != scalar(@image_dirs) || 
		scalar(@image_nums) != scalar(@image_files)) {
		die "Expected same number of image files as image dirs, instead found " .
			 scalar(@image_files) . " and " . scalar(@image_dirs);
	}

	my $max_digit_count = length($image_nums[floor($#image_nums/2)]);

	for (1 .. floor($#image_files/2)) {
		my $image_name = sprintf("%0" . $max_digit_count . "d", $_);
		my $target_file = catfile($cfg{data_folder} . "_reduced", $cfg{exp_name}, $cfg{$target_dir}, $image_name . ".png");
		
		mkpath(dirname($target_file));
		copy($image_files[2*$_], $target_file);
	}

	#also copy over the config file, which will need to be updated by hand
	copy($opt{cfg}, catfile($cfg{data_folder} . "_reduced", $cfg{exp_name}, basename($opt{cfg}))); 
}

################################################################################
#Documentation
################################################################################

=head1 NAME

change_image_sample.pl: sample the images from an experiment to create a new
experiment with half as many images, used for testing robustness of image
sampling rate

=head1 SYNOPSIS

change_image_sample.pl -cfg config_file

=head1 Description

In order to test the effect of reducing the image sampling rate on the outcome
of experiment, some automated method to downsample an experiment's raw images
was needed. This program looks for cell mask and FA image files in a given
results folder and copies those files back into a new experimental data folder
in the data directory. The name of the new experiment will be $cfg{data_folder}
. "_reduced". All the appropriate config files will also be copied, but those
will need to be modified to reflect the new data file locations.

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=back

Optional parameter(s):

=over 

=item * debug or d: print debuging information during program execution

=back

=head1 EXAMPLES

./change_image_sample.pl -cfg ../../data/focal_adhesions/time_series_01/analysis.cfg

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 8/18/2009 

=cut
