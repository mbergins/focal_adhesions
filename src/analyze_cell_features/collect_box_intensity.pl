#!/usr/bin/env perl

################################################################################
# Global Variables and Modules
################################################################################

use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Path;
use File::Find;
use File::Spec::Functions;
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
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

################################################################################
# Main Program
################################################################################
my $tracking_folder = $cfg{tracking_output_file};
$tracking_folder =~ s/(.*)\.csv/$1/;

our @files;
find(\&include_in_vis, catdir($cfg{exp_results_folder}, $cfg{movie_output_folder}, $tracking_folder));
my @matlab_code = "find_box_intensity('$files[0]')";

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'Box_intensity');
$opt{error_file} = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'Box_intensity', 'error.txt');
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

################################################################################
#Functions
################################################################################

sub include_in_vis {
    my $tracking_file = $cfg{tracking_output_file};
    $tracking_file =~ s/(.*)\.csv/$1/;

    if (   $File::Find::name =~ /\.m/
        && $File::Find::name =~ /$tracking_file/ ) {
        push @files, $File::Find::name;
    }
}

################################################################################
#Documentation
################################################################################

=head1 NAME

collect_box_intensity.pl - Executes the MATLAB programs designed collect the
boxed intensity values for each focal adhesion

=head1 SYNOPSIS

collect_box_intensity.pl -cfg FA_config

=head1 Description

This program runs the MATLAB code that collects the intensity over a box which
encloses the entire focal adhesion's life cycle. This process mimics the methods
used by other researchers. The output is placed in the same folder as the other
time series properties collections.

Required parameter(s):

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

Documentation last updated: 7/3/2008 

=cut
