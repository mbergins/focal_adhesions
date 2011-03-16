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
    
	$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'signif_R_models');
    $opt{resource} = 'blade';
    if (defined $cfg{job_group}) {
        $opt{job_group} = $cfg{job_group};
    }
    
    &FA_job::send_general_lsf_program(\@commands,\%opt);

    exit(0);
}
my $dir_glob;
if ($cfg{exp_results_folder} =~ /$cfg{results_folder}(.*)/) {
	my @split_dirs = &File::Spec::Functions::splitdir($1);
	$dir_glob = "/*" x scalar(@split_dirs);
} else {
	die "Couldn't find results folder in exp_results_folder variable";
}

my $output_base = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'signif_R_models');
if (! -e $output_base) {
	mkpath($output_base);
}
$opt{error_file} = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'signif_R_models', 'error.txt');

my @R_cmds;
push @R_cmds, "R CMD BATCH --vanilla \"--args exp_dir=$cfg{results_folder}$dir_glob/adhesion*/models/\" find_significant_rate_models.R " . catfile($output_base,'signif_models.Rout');


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
