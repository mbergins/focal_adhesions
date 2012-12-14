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
use Math::Matlab::Extra;
use Emerald;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
$opt{first} = -1;
$opt{max} = -1;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l", 'first|f=i', 'max|m=i') or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

################################################################################
# Main Program
################################################################################

my @matlab_code;

#my $folder   = $cfg{adhesion_image_folder};
my $first_image = catfile($cfg{individual_results_folder}, '001', $cfg{raw_mask_file});
my $out_file = $cfg{adhesion_image_file};

next if (not(defined($first_image)));
#next if (not(defined($folder)));

#Remove an ' marks used in config file to keep the folder name together
$first_image =~ s/\'//g;
#$folder =~ s/\'//g;

my @image_files = sort <$cfg{individual_results_folder}/*>;
#my @image_files = sort <$cfg{exp_data_folder}/$folder/*>;
my @image_files = map { $_ =~ s/\'//g; $_; } @image_files;

if ($opt{debug}) {
    if (scalar(@image_files) > 1) {
        print "Image files found: $image_files[0] - $image_files[$#image_files]\n";
    } elsif (scalar(@image_files) == 0) {
        #print "No image files found matching $cfg{exp_data_folder}/$folder, moving onto next image set.\n";
        print "No image files found matching $cfg{individual_results_folder}, moving onto next image set.\n";
        next;
    } else {
        print "Image file found: $image_files[0]\n";
    }
} else {
    next if (not @image_files);
}

push @matlab_code, &create_matlab_code(\@image_files, $out_file);

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'edge_velo');
$opt{error_file} = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'edge_velo', 'error.txt');

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

################################################################################
#Functions
################################################################################

sub create_matlab_code {
    my @image_files = @{ $_[0] };
    my $out_file    = $_[1];

    #my @image_stack_count = map { Image::Stack::get_image_stack_number($_) } @image_files;

    my @matlab_code;
    #if (grep { $_ > 1 } @image_stack_count) {
    #    if (scalar(@image_files) > 1) {
    #        die "Found more than one image stack in: ", join(", ", @image_files), "\n",
    #          "Expected single image stack or multiple non-stacked files\n";
    #    }
    #    die "Found image stack, need to implement way to use stack files in edge velocity code";
    #    @matlab_code = &create_matlab_code_stack(\@image_files, $out_file);
    #} else {
        @matlab_code = &create_matlab_code_single(\@image_files, $out_file);
    #}
    return @matlab_code;
}

sub create_matlab_code_single {
    my @image_files = @{ $_[0] };
    
    my $first_dir = $image_files[0];
    if ($opt{first} > 0) {
        foreach my $dir (@image_files) {
            if (basename($dir) >= $opt{first}) {
                $first_dir = $dir;
                last;
            }
        }
    }
    my $first_file = catfile($first_dir, $cfg{raw_mask_file});
    
    my $max_img = scalar(@image_files);
    if ($opt{max} > 0) {
        $max_img = $opt{max};
    }
    
    my $results_folder = catdir($cfg{exp_results_folder}, 'pax_edge_velocity');

    my $matlab_code = "edge_velocity_wrapper('contr',0,'protrusion',1,'t_step',1,'file','$first_file','results','$results_folder','max_img',$max_img";

    my @exclude_img_nums = @{ $cfg{exclude_image_nums} };

    if (scalar(@exclude_img_nums) > 0) {
        my $exclude_imgs = join(',', @exclude_img_nums);
        $matlab_code .= ",'exclude_imgs','$exclude_imgs'";
    }
    
    $matlab_code .= ")";

    return $matlab_code;
}

################################################################################
#Documentation
################################################################################
