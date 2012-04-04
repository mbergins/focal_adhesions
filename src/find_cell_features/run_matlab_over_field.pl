#!/usr/bin/env perl

################################################################################
# Global Variables and Modules
################################################################################

use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Path;
use File::Spec::Functions;
use File::Spec;
use File::Basename;
use Image::ExifTool;
use Getopt::Long;
use Data::Dumper;

use Config::Adhesions;
use Math::Matlab::Extra;
use Emerald;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "script=s", "debug|d", "lsf|l", "queue=s",
	"extra|e=s") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};
die "Can't find script specified on the command line" if not exists $opt{script};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

if ($opt{script} =~ /(.*)\.m/) {
	$opt{script} = $1;
}

$opt{script_dir} = dirname($opt{script});
chdir $opt{script_dir};
if ($opt{script} =~ /$opt{script_dir}\/(.*)/) {
	$opt{script} = $1;
}
$opt{root_mwd} = File::Spec->rel2abs($opt{script_dir});

################################################################################
# Main Program
################################################################################

###########################################################
# Run Decisions
###########################################################
if ($opt{script} =~ /apply_bleaching_correction/) {
	#using this to check decide whether or not to run the photo bleach
	#correction, if undefined or not true, exit from the program without
	#running any matlab jobs
	if (! (defined $cfg{photo_bleach_correction} && $cfg{photo_bleach_correction})) {
		exit;
	} 	
}

if ($opt{script} =~ /find_background_intensity/) {
	if (! (defined $cfg{background_correction} && $cfg{background_correction})) {
		exit;
	} 	
}

###########################################################
# Build and Execute
###########################################################

my $extra = &build_extra_command_line_opts;

my @matlab_code = ("$opt{script}('$cfg{exp_results_folder}'$extra)\n");

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, $opt{script});
$opt{error_file} = catfile($opt{error_folder}, 'error.txt');
$opt{runtime} = "1";
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

################################################################################
# Functions
################################################################################

sub build_extra_command_line_opts {
	my $extra = '';
	
	if ($opt{script} eq "make_filtered_vis") {
		if ($opt{extra} eq "lifetime") {
			$extra = ",'type','lifetime','min_value',10";
		} elsif ($opt{extra} eq "ratio") {
			$extra = ",'type','ratio'";
		} elsif ($opt{extra} eq "FA_angle") {
			$extra = ",'type','FA_angle'";
		} else {
			$extra = ",'type','FA_dist'";
		}
	}
	
	if ($opt{script} eq "recenter_FA_angles") {
		if (defined $cfg{by_hand_direction}) {
			$extra = ",'by_hand_direction',$cfg{by_hand_direction}";
		}
	}

	return $extra;
}
