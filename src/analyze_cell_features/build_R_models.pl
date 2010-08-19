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
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l", "model_type=s")
  or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

################################################################################
# Main Program
################################################################################
my @model_types = qw(average cell_background local_background area box_intensity);
#push @model_types, "background_correlation_model";
if ($opt{lsf}) {
    my @commands;
    foreach (@model_types) {
        #$0 - the name of the program currently running, used to protect against
        #future file name changes
        push @commands, "$0 -cfg $opt{cfg} -model_type $_";
    }
    
    $opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'R_models');
    if (defined $cfg{job_group}) {
        $opt{job_group} = $cfg{job_group};
    }
    
    &FA_job::send_general_lsf_program(\@commands,\%opt);

    exit(0);
}

my $data_dir = catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_ts_folder});

my $output_base = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'R_models');
if (! -e $output_base) {
	mkpath($output_base);
}

my @R_cmds;
if (defined($opt{model_type})) {
	my $output_file = catfile($output_base, 'R_out_' . $opt{model_type} . '.txt');
	push @R_cmds, "R CMD BATCH --vanilla \"--args data_dir=$data_dir " .
	  "model_type=$opt{model_type} min_length=$cfg{min_linear_model_length}\" " . 
	  "FA_analysis_lib.R $output_file";
} else {
	for (@model_types) {
		my $output_file = catfile($output_base, 'R_out_' . $_ . '.txt');
		push @R_cmds, "R CMD BATCH --vanilla \"--args data_dir=$data_dir " .
		  "model_type=$_ min_length=$cfg{min_linear_model_length}\" " . 
		  "FA_analysis_lib.R $output_file";
	}
}

$opt{error_file} = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'R_models', 'error.txt');

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


################################################################################
#Documentation
################################################################################

=head1 NAME

build_R_models.pl - Run the R commands needed to produce the assembly and
disasembly rate models

=head1 SYNOPSIS

build_R_models.pl -cfg FA_config

=head1 Description

The code needed to build and output the R model files is contained in FA_analysis_lib.R. The R program takes three command line options:

=over

=item * data_dir - the location of the lineage time series files
=item * model_type - the type of model that will be built
=item * debug - turns on debuging mode (optional)

=back

All of these parameters can be set by changing the perl program

=head1 Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=back

Optional parameter(s):

=over 

=item * debug or d: print debuging information during program execution

=item * lsf: submit jobs through the emerald queuing system

=back

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 7/7/2009

=cut
