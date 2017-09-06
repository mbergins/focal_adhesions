#!/usr/bin/env perl

################################################################################
# Global Variables and Modules
################################################################################

use lib "../lib";

use strict;
use File::Path;
use File::Spec::Functions;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use POSIX;

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
die "Expected to find the same number of image files as folders in the results " . 
	"directory ($cfg{individual_results_folder})."
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

my @matlab_code;
if (defined($cfg{min_independent_size}) && ! $cfg{no_ad_splitting}) {
	@matlab_code = &create_bundled_matlab_commands(10);
} else {
	@matlab_code = &create_bundled_matlab_commands;
}

$opt{queue} = "day";
$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'FA');
$opt{error_file} = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'FA', 'error.txt');
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

################################################################################
#Functions
################################################################################

sub create_all_matlab_commands {
    my @matlab_code;

    my @image_files = <$cfg{individual_results_folder}/*/$cfg{adhesion_image_file}>;
	my $extra_opt = &build_extra_opts;
    foreach my $file_name (@image_files) {
		$matlab_code[0] .= "find_focal_adhesions('$file_name'$extra_opt)\n";
    }

    return @matlab_code;
}

sub create_bundled_matlab_commands {
	my $bundled_count = 150;
	if (scalar(@_) > 0) {
		$bundled_count = $_[0];
	}
	
    my @matlab_code;

    my @image_files = <$cfg{individual_results_folder}/*/$cfg{adhesion_image_file}>;
	my $extra_opt = &build_extra_opts;
	
	my @image_nums = 1..scalar(@image_files);
	my $data_sets_runs = floor(scalar(@image_files)/$bundled_count);
	$data_sets_runs = 1 if ($data_sets_runs < 1);
	my $min_image_count = ceil(scalar(@image_files)/$data_sets_runs);
	while (@image_nums) {
		my @this_set;
		for (1..$min_image_count) {
			if (@image_nums) {
				push @this_set, shift @image_nums;
			}
		}
		my $num_str = '[' . join(",",@this_set) . ']';
			
		my $these_opt = "$extra_opt,'i_nums',$num_str";

		push @matlab_code, "find_focal_adhesions_full_exp('$cfg{exp_results_folder}'$these_opt)\n";
	}

    return @matlab_code;
}

sub create_single_matlab_command {
    my @matlab_code;
	
	my $extra_opt = &build_extra_opts;

	$matlab_code[0] .= "find_focal_adhesions_full_exp('$cfg{exp_results_folder}'$extra_opt)\n";

    return @matlab_code;
}

sub build_extra_opts {
	my $extra_opt = "";
	if (defined $cfg{stdev_thresh}) {
		my @split_stdev_vals = split(/\s+/,$cfg{stdev_thresh});
		$extra_opt .= ",'stdev_thresh',[" . join(",",@split_stdev_vals) . "]";
	}
	if (defined $cfg{no_ad_splitting}) {
		$extra_opt .= ",'no_ad_splitting',$cfg{no_ad_splitting}";
	}
	if (defined $cfg{min_adhesion_size}) {
		$extra_opt .= ",'min_adhesion_size',$cfg{min_adhesion_size}";
	}
	if (defined $cfg{max_adhesion_size}) {
		$extra_opt .= ",'max_adhesion_size',$cfg{max_adhesion_size}";
	}
	if (defined $cfg{min_independent_size}) {
		$extra_opt .= ",'min_independent_size',$cfg{min_independent_size}";
	}
	if (defined $cfg{max_adhesion_count}) {
		$extra_opt .= ",'max_adhesion_count',$cfg{max_adhesion_count}";
	}
	if (defined $cfg{proximity_filter}) {
		$extra_opt .= ",'proximity_filter',$cfg{proximity_filter}";
	}
	if (defined $cfg{min_seed_size}) {
		$extra_opt .= ",'min_seed_size',$cfg{min_seed_size}";
	}
	if (defined $cfg{per_image_thresh}) {
		$extra_opt .= ",'per_image_thresh',$cfg{per_image_thresh}";
	}
	if (defined $cfg{confocal_mode}) {
		$extra_opt .= ",'confocal_mode',$cfg{confocal_mode}";
	}
	if (defined $cfg{atrous_segmentation}) {
		$extra_opt .= ",'atrous_segmentation',$cfg{atrous_segmentation}";
	}
	if (defined $cfg{structure_element_size}) {
		$extra_opt .= ",'structure_element_size',$cfg{structure_element_size}";
	}
	if (defined $cfg{atrous_export_level}) {
		$extra_opt .= ",'atrous_export_level',$cfg{atrous_export_level}";
	}
	
	return $extra_opt;
}

################################################################################
#Documentation
################################################################################

=head1 NAME

collect_fa_image_set.pl - Executes the MATLAB programs designed collect the
masks which define the focal adhesions

=head1 SYNOPSIS

collect_fa_image_set.pl -cfg FA_config

=head1 Description

This program is used to create all the binary mask files which define the
location of the focal adhesions, which will be used in subsequent steps. The
primary logic of the program is in a set of MATLAB scripts which do all the
image analysis/writing and also collects properties of the focal adhesions and
writes those to file.

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

collect_fa_image_set.pl -cfg FA_config

OR

collect_fa_image_set.pl -cfg FA_config -d

=head1 SEE ALSO

collect_mask_set.pl: similar program designed to collect the binary mask that
locates the intracellular area

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 7/3/2008 

=cut
