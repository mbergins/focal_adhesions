#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################

use lib "../lib";

use strict;
use File::Path;
use File::Spec::Functions;
use File::Basename;
use Math::Matlab::Local;
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
$opt{debug}   = 0;
$opt{emerald} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

###############################################################################
# Main Program
###############################################################################

my @image_folders = <$cfg{individual_results_folder}/*>;
my @image_files   = <$cfg{individual_results_folder}/*/$cfg{raw_mask_file}>;
if (scalar(@image_files) == 0) {
    print "Didn't find any cell mask image files, aborting.\n";
    exit 0;
}
warn "Expected to find the same number of image files as folders in the results directory ($cfg{individual_results_folder})." if (scalar(@image_files) != scalar(@image_folders));

if ($opt{debug}) {
    if (scalar(@image_files) > 1) {
        print "Cell mask files found: $image_files[0] - $image_files[$#image_files]\n";
    } elsif (scalar(@image_files) == 0) {
        warn "Couldn't find any cell mask files in $cfg{individual_results_folder} subfolders\n\n";
    } else {
        print "Cell mask file found: $image_folders[0]\n";
    }
}

my @matlab_code = &create_single_matlab_command;

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'mask_set');
$opt{error_file} = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'mask_set', 'error.txt');
$opt{runtime} = "0:5";
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

###############################################################################
#Functions
###############################################################################

sub create_all_matlab_commands {
    my @matlab_code;

    my @image_files   = <$cfg{individual_results_folder}/*/$cfg{raw_mask_file}>;
    foreach my $file_name (@image_files) {
        $matlab_code[0] .= "find_cell_mask('$file_name')\n";
    }

    return @matlab_code;
}

sub create_single_matlab_command {
    my @matlab_code;
	
	my $extra_opts = '';
	if (defined $cfg{mask_threshold}) {
		$extra_opts .= ",'mask_threshold',$cfg{mask_threshold}";
	}
	if (defined $cfg{median_cell_mask_filter}) {
		$extra_opts .= ",'median_filter',$cfg{median_cell_mask_filter}";
	}
	if (defined $cfg{single_mask_threshold}) {
		$extra_opts .= ",'single_threshold',$cfg{single_mask_threshold}";
	}

	$matlab_code[0] = "find_cell_mask_full_exp('$cfg{exp_results_folder}'$extra_opts)\n";

    return @matlab_code;
}

###############################################################################
#Documentation
###############################################################################

=head1 NAME

collect_mask_set.pl - Executes the MATLAB programs designed collect the binary
masks which define the intracellular region of a cell

=head1 SYNOPSIS

collect_mask_set.pl -cfg FA_config

=head1 DESCRIPTION

This program is used to create all the binary mask files which define the
location of the intracellular space, which will be used in subsequent steps. The
primary logic of the program is in a set of MATLAB scripts which do all the
threshold selection and image creation/writing. This perl script is more or less
an intelligent wrapper that sets up all the MATLAB commands.  

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=back

Optional parameter(s):

=over 

=item * debug or d: print debuging information during program execution

=item * emerald: submit jobs through the emerald queuing system

=back

=head1 EXAMPLES

collect_mask_set.pl -cfg FA_config

OR

collect_mask_set.pl -cfg FA_config -d

=head1 SEE ALSO

collect_fa_image_set.pl: similar program designed to collect the binary mask
that locates focal adhesions

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 4/10/2008

=cut
