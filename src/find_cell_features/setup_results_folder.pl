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
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l","convert_to_png")
  or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

################################################################################
# Main Program
################################################################################

if (not $opt{debug}) {
	mkpath($cfg{individual_results_folder});
}
my @image_sets = (
    [qw(raw_mask_folder raw_mask_file)],
    [qw(raw_kinase_folder raw_kinase_file)],
    [qw(adhesion_image_folder_generic adhesion_image_file)],
    [qw(adhesion_image_folder_secondary adhesion_image_secondary_file)],
    [qw(adhesion_image_folder_pax adhesion_image_file)],
    [qw(adhesion_image_folder_vin adhesion_image_file)],
    [qw(adhesion_image_folder_fak adhesion_image_file)],
    );
my @matlab_code;
my $all_images_empty = 1;

foreach (@image_sets) {
    my $folder   = $cfg{ $_->[0] };
    my $out_file = $cfg{ $_->[1] };
 	
	next if (not(defined($folder)));

    #Remove an ' marks used in config file to keep the folder name together
    $folder =~ s/\'//g;
	
	if ($opt{debug}) {
		print "Searching: $cfg{exp_data_folder}/$folder/\n";
	}
    my $search_dir = catdir($cfg{exp_data_folder},$folder);

    if (-e $search_dir) {
		$all_images_empty = 0;
        print "For Config Variable: ", $_->[0], "\nFound $search_dir\n\n" if $opt{debug};
    } else {
		print "\n" if $opt{debug};
        next;
    }
	
	if ($opt{convert_to_png}) {
		&convert_data_file_to_png($search_dir);
	}

    push @matlab_code, "setup_results_folder('$search_dir','$cfg{individual_results_folder}', '$out_file');";
}
die "Unable to find any images to include in the new experiment" if $all_images_empty;

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'setup');
$opt{error_file} = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'setup', 'error.txt');
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

################################################################################
#Functions
################################################################################

sub convert_data_file_to_png {
	my $data_dir = $_[0];
	my @data_files = <$data_dir/*>;
	
	#Only run this if there is a single file in the data directory
	if (scalar(@data_files) == 1) {
		my $command = "convert $data_files[0] " . catdir($data_dir,"data_%05d.png") . ";";
		if ($opt{debug}) {
			print "$command\n";
		} else {
			system($command);
			unlink $data_files[0];
		}
	}
}

################################################################################
#Documentation
################################################################################

=head1 NAME

setup_results_folder.pl - Move all the raw data files into the proper locations
in the results folder

=head1 SYNOPSIS

setup_results_folder.pl -cfg FA_config

=head1 Description

Since the the data sets being used come in multiple forms, mostly split and
stacked image sequences, it became easier to write a single function to move all
the data files into a standard results directory. Then all the downstream
programs would have an easier time trying to find specific images. This program
scans through a the data directory for an experiment and moves all the files
into the correct locations.

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

setup_results_folder.pl -cfg FA_config

OR

setup_results_folder.pl -cfg FA_config -d

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 4/10/2008 

=cut
