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
use Math::Matlab::Local;
use Getopt::Long;
use Data::Dumper;

use Config::Adhesions qw(ParseConfig);
use Math::Matlab::Extra;
use Emerald;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l", "data_file=s")
  or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

################################################################################
# Main Program
################################################################################
my @data_files = qw(Average_adhesion_signal.csv Average_adhesion_signal.csv.gz);
if ($opt{lsf}) {
    my @commands;
    foreach (@data_files) {
        #$0 - the name of the program currently running, used to protect against
        #future file name changes
        push @commands, "$0 -cfg $opt{cfg} -data_file $_";
    }
    
    $opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'build_rate_models');
    if (defined $cfg{job_group}) {
        $opt{job_group} = $cfg{job_group};
    }
    
    &FA_job::send_general_lsf_program(\@commands,\%opt);

    exit(0);
}

my $data_dir = catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder});

my $output_base = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'R_models');
if (! -e $output_base) {
	mkpath($output_base);
}

my @R_cmds;
if (defined($opt{data_file})) {
	my $this_data_file = catfile($data_dir,'lin_time_series',$cfg{data_file});
	if (not -e $this_data_file) {
		warn "Could not find file for $cfg{data_file} ($this_data_file), skipping";
		next;
	}
	my $output_file = catfile($output_base, 'R_out_' . $opt{data_file} . '.txt');
	push @R_cmds, "R CMD BATCH --vanilla \"--args data_dir=$data_dir " .
	  "data_file=$opt{data_file} time_spacing=$cfg{time_spacing}  " . 
	  "min_length=$cfg{min_linear_model_length}\" " . 
	  "bilinear_modeling.R $output_file";
} else {
	for (@data_files) {
		my $this_data_file = catfile($data_dir,'lin_time_series',$_);

		if (not -e $this_data_file) {
			# warn "Could not find file for $_ ($this_data_file), skipping";
			next;
		}

		my $output_file = catfile($output_base, 'R_out_' . $_ . '.txt');
		push @R_cmds, "R CMD BATCH --vanilla \"--args data_dir=$data_dir " .
		  "data_file=$_ time_spacing=$cfg{time_spacing} " .  
		  "min_length=$cfg{min_linear_model_length}\" " . 
		  "bilinear_modeling.R $output_file";
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

The code needed to build and output the R model files is contained in
bilinear_modeling.R. The R program takes three command line options:

=over

=item * data_dir - the location of the lineage time series files
=item * data_file - the file with the data that will be used to build the model
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
