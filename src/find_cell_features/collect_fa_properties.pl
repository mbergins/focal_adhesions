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

my @image_folders = <$cfg{individual_results_folder}/*>;
my @image_files   = <$cfg{individual_results_folder}/*/$cfg{adhesion_image_file}>;
die "Expected to find the same number of image files as folders in the results directory ($cfg{individual_results_folder})."
  if (scalar(@image_files) != scalar(@image_folders));

if ($opt{debug}) {
    if (scalar(@image_files) > 1) {
        print "Focal image files found: $image_files[0] - $image_files[$#image_files]\n";
    } elsif (scalar(@image_files) == 0) {
        warn "Couldn't find any focal image files in $cfg{individual_results_folder} subfolders\n\n";
    } else {
        print "Focal image file found: $image_folders[0]\n";
    }
}

my @matlab_code = &create_all_matlab_commands;

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'FA_props');
$opt{error_file} = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'FA_props', 'error.txt');
$opt{runtime} = "1";
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

################################################################################
#Functions
################################################################################

sub create_all_matlab_commands {
    my @matlab_code;

    my @raw_image_files = grep -e $_, <$cfg{individual_results_folder}/*/$cfg{adhesion_image_file}>;
    my @adhesion_image_files = grep -e $_, <$cfg{individual_results_folder}/*/adhesions.png>;
    if (scalar(@raw_image_files) != scalar(@adhesion_image_files)) {
        die "Expected to find equal number of raw data adhesion image and adhesion mask images." . 
            "Instead found " . scalar(@raw_image_files) . " and " . scalar(@adhesion_image_files);
    }
    
    foreach (0..$#raw_image_files) {
        my $raw_image_file = $raw_image_files[$_];
        my $adhesion_image_file = $adhesion_image_files[$_];
        my $cell_mask = catfile(dirname($raw_image_file), $cfg{cell_mask_file});
        
	    #protrusion_data_file is not defined in any config file
        my $protrusion_file = catfile($cfg{exp_results_folder}, $cfg{protrusion_data_file});
        
        my $extra_opt = "";
        if (-e $cell_mask) {
            $extra_opt .= ",'cell_mask','$cell_mask'";
        }
        if (-f $protrusion_file) {
            $extra_opt .= ",'protrusion_file','$protrusion_file'";
        }
        
        my $this_command = "find_adhesion_properties('$raw_image_file','$adhesion_image_file'";
        $this_command .= "$extra_opt);\n";

        $matlab_code[0] .= $this_command;
    }

    return @matlab_code;
}
