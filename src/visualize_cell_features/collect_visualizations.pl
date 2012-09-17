#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################
use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Path;
use File::Basename;
use File::Find;
use Image::ExifTool;
use Math::Matlab::Local;
use Getopt::Long;
use Data::Dumper;
use File::Spec::Functions;
use Benchmark;
use POSIX;

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
GetOptions(\%opt, "cfg=s", "debug|d", "movie_debug", "config_only|only_config", 
                  "lsf|l", "single_ad_folders") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my %cfg = ParseConfig(\%opt);

###############################################################################
#Main Program
###############################################################################
my @matlab_code;
push @matlab_code, "make_tracking_visualization('$cfg{exp_results_folder}')";
push @matlab_code, "make_single_ad_frames('$cfg{exp_results_folder}','min_longevity',10)";

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'visualization');
$opt{error_file} = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'visualization', 'error.txt');
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

###############################################################################
#Functions
###############################################################################
