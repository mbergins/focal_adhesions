#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Path;
use File::Find;
use File::Spec::Functions;
use Getopt::Long;
use IO::File;
use Benchmark;

use lib "../lib";
use Config::Adhesions;
use Image::Data::Collection;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|config=s", "debug|d") or die;
die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Overall Configuration\n\n" if $opt{debug};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg = $ad_conf->get_cfg_hash;

###############################################################################
# Main Program
###############################################################################
my $debug_string = "";
$debug_string = "-d" if $opt{debug};

my ($t1, $t2);

open STATUS, ">" . catfile($cfg{exp_data_folder},"status.txt");
print STATUS "RUNNING";
close STATUS;

if (-e $cfg{exp_results_folder}) {
    File::Path::rmtree($cfg{exp_results_folder});
}


#Find Features
chdir "../find_cell_features";
print "\n\nSetting Up Results Directory\n\n" if $opt{debug};
$t1 = new Benchmark;
system "./setup_results_folder.pl -cfg $opt{cfg} $debug_string";
$t2 = new Benchmark;
print "Runtime: ",timestr(timediff($t2,$t1)), "\n" if $opt{debug};

print "\n\nCollecting Cell Mask Set\n\n" if $opt{debug};
$t1 = new Benchmark;
system "./collect_mask_set.pl -cfg $opt{cfg} $debug_string";
$t2 = new Benchmark;
print "Runtime: ",timestr(timediff($t2,$t1)), "\n" if $opt{debug};

print "\n\nFinding Focal Adhesions\n\n" if $opt{debug};
$t1 = new Benchmark;
system "./collect_fa_image_set.pl -cfg $opt{cfg} $debug_string";
$t2 = new Benchmark;
print "Runtime: ",timestr(timediff($t2,$t1)), "\n" if $opt{debug};


#Building the Tracking Matrix
chdir "../analyze_cell_features";
print "\n\nTracking Focal Adhesions\n\n" if $opt{debug};
$t1 = new Benchmark;
system "./track_adhesions.pl -cfg $opt{cfg} -o data.stor -i data.stor $debug_string";
$t2 = new Benchmark;
print "Runtime: ",timestr(timediff($t2,$t1)), "\n" if $opt{debug};

#Analyze the Adhesions
print "\n\nAnalyzing Focal Adhesions\n\n" if $opt{debug};
$t1 = new Benchmark;
system "./gather_tracking_results.pl -cfg $opt{cfg} $debug_string";
$t2 = new Benchmark;
print "Runtime: ",timestr(timediff($t2,$t1)), "\n" if $opt{debug};

#Filtering the Tracking Matrix
print "\n\nFiltering the Tracking Matrix\n\n" if $opt{debug};
$t1 = new Benchmark;
system "./filter_tracking_matrix.pl -cfg $opt{cfg} $debug_string";
$t2 = new Benchmark;
print "Runtime: ",timestr(timediff($t2,$t1)), "\n" if $opt{debug};

#Collecting Visualizations
chdir "../visualize_cell_features";
print "\n\nBuilding Tracking Visualization\n\n" if $opt{debug};
$t1 = new Benchmark;
system "./collect_visualizations.pl -cfg $opt{cfg} $debug_string";
$t2 = new Benchmark;
print "Runtime: ",timestr(timediff($t2,$t1)), "\n" if $opt{debug};

#Build Movies 
my @image_numbers = &Image::Data::Collection::gather_sorted_image_numbers(\%cfg);
my $image_num_length = length(scalar(@image_numbers));

print "\n\nBuild Movies\n\n" if $opt{debug};
$t1 = new Benchmark;

our @movie_dirs; 
find(\&add_to_movie_dir, (catdir($cfg{exp_results_folder},'movies')));

foreach my $f1 (@movie_dirs) {
    foreach my $f2 (@{ $cfg{movie_output_prefix} }) {
        my $input_folder = catdir($f1,$f2);
        system "ffmpeg -v 0 -y -r $cfg{movie_frame_rate} -b $cfg{movie_bit_rate} -i $input_folder/%0" . $image_num_length . "d.png $input_folder.mov 2>&1";
    }
}
$t2 = new Benchmark;
print "Runtime: ",timestr(timediff($t2,$t1)), "\n" if $opt{debug};

open STATUS, ">" . catfile($cfg{exp_data_folder},"status.txt");
print STATUS "DONE";
close STATUS;

###############################################################################
# Functions
###############################################################################

sub add_to_movie_dir {
    if ($File::Find::name =~ /$cfg{vis_config_file}/) {
        push @movie_dirs, $File::Find::dir;
    }
}
