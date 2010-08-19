#!/usr/bin/env perl

################################################################################
# Global Variables and Modules
################################################################################

use lib "../lib";
use lib "../lib/perl";

use strict;
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
GetOptions(\%opt, "cfg|c=s", "debug|d")
  or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

die "The split number must be specified on config file with the option \"split_num\"" if not exists $cfg{split_num};

################################################################################
# Main Program
################################################################################

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

	my $max_digit_count = length($image_nums[-1]);

	my @split_index = grep $image_nums[$_] == $cfg{split_num}, (0 .. $#image_nums);
	die if (scalar(@split_index) > 1);

	for (0 .. $split_index[0]) {
		my $image_name = sprintf("%0" . $max_digit_count . "d", $_ + 1);
		my $target_file = catdir($cfg{data_folder}, $cfg{exp_name} . "_pre", 
								 $cfg{$target_dir}, $image_name . ".png");
		mkpath(dirname($target_file));
		copy($image_files[$_], $target_file);
	}
	my $new_exp_name = $cfg{exp_name} . "_pre";
	my $new_config_file = catfile($cfg{data_folder},$new_exp_name , basename($opt{cfg}));
	copy($opt{cfg}, $new_config_file); 
	system("sed -i -r 's/$cfg{exp_name}/$new_exp_name/' $new_config_file");

	for (($split_index[0] + 1) .. $#image_files) {
		my $image_name = sprintf("%0" . $max_digit_count . "d", $_ + 1 - $cfg{split_num});
		my $target_file = catdir($cfg{data_folder}, $cfg{exp_name} . "_post", 
								 $cfg{$target_dir}, $image_name . ".png");
		mkpath(dirname($target_file));
		copy($image_files[$_], $target_file);
	}
	$new_exp_name = $cfg{exp_name} . "_post";
	$new_config_file = catfile($cfg{data_folder},$new_exp_name , basename($opt{cfg}));
	copy($opt{cfg}, $new_config_file); 
	system("sed -i -r 's/$cfg{exp_name}/$new_exp_name/' $new_config_file");
}

################################################################################
#Documentation
################################################################################

=head1 NAME

split_exp.pl - Move all the raw data files into the proper locations
in the results folder

=head1 SYNOPSIS

split_exp.pl -cfg FA_config -split_num ###

=head1 Description

Occasionally it is helpful to analyze a single experiment as though it were two
independent experiments, but splitting an experimental data set by hand is time
consuming. This script analyzes the directory structure and files output by
setup_results_folder.pl and puts the raw data images back into the data
directory with the proper numbering and directory structure. The script only
looks for adhesion files and cell mask data files, but that can be extended. The
split folders are named simply $cfg{exp_name} . "_pre" or "_post", with the pre
data being before the specified split image number.

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file
=item * split_num: the focal adhesion analysis config file

=back

Optional parameter(s):

=over 

=item * debug or d: print debuging information during program execution

=back

=head1 EXAMPLES

split_exp.pl -cfg FA_config.cfg -split_num 50

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 7/15/2009 

=cut
