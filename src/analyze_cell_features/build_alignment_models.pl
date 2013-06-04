#!/usr/bin/env perl

################################################################################
# Global Variables and Modules
################################################################################

use lib "../lib";

use strict;
use File::Path;
use File::Spec::Functions;
use File::Basename;
use File::Copy;
use Getopt::Long;
use Data::Dumper;

use Config::General;
use Config::Adhesions qw(ParseConfig);
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
if ($opt{lsf}) {
    my @commands;
    push @commands, "$0 -cfg $opt{cfg}";
    
    $opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'alignment_models');
    if (defined $cfg{job_group}) {
        $opt{job_group} = $cfg{job_group};
    }
    
    &FA_job::send_general_lsf_program(\@commands,\%opt);

    exit(0);
}

#######################################
# Output file setup
#######################################

my $data_dir = catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder},'lin_time_series');
if (! -e $data_dir) {
	&copy_over_single_image_props;
	if (! -e $data_dir) {
		warn "Unable to find lin time series dir ($data_dir), exiting.";
		exit;
	}
}

$opt{error_folder} = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'alignment_models');
if (! -e $opt{error_folder}) {
	mkpath($opt{error_folder});
}

my $output_file = catfile($opt{error_folder}, 'R_out.txt');

my $best_angle = '';
if (defined $cfg{fixed_best_angle}) {
	$best_angle=" fixed_best_angle=$cfg{fixed_best_angle}";
}

my $min_ratio = '';
if (defined $cfg{FAAI_min_ratio}) {
	$min_ratio=" min.ratio=$cfg{FAAI_min_ratio}";
}

#######################################
# Output file setup
#######################################
my @R_cmds;
push @R_cmds, "R CMD BATCH --vanilla \"--args time_series_dir=$data_dir $best_angle $min_ratio\" FA_alignment_search.R $output_file";

$opt{error_file} = catfile($opt{error_folder}, 'error.txt');

for (@R_cmds) {
	if ($opt{debug}) {
		print "$_\n";
	} else {
		system($_);
	}
}

################################################################################
#Functions
################################################################################

sub copy_over_single_image_props {
	my @image_folders = <$cfg{individual_results_folder}/*/raw_data>;
	if (scalar(@image_folders) == 1) {
		my $lin_time_series_dir = catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder},'lin_time_series');
		mkdir $lin_time_series_dir;

		copy(catfile($image_folders[0],'Orientation.csv'),$lin_time_series_dir);
		copy(catfile($image_folders[0],'Area.csv'),$lin_time_series_dir);
		copy(catfile($image_folders[0],'MajorAxisLength.csv'),$lin_time_series_dir);
		copy(catfile($image_folders[0],'MinorAxisLength.csv'),$lin_time_series_dir);
	}
}

################################################################################
#Documentation
################################################################################
