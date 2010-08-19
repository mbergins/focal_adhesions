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
my ($t1, $t2);

#Collecting Visualizations
chdir "../visualize_cell_features";

#Build Movies 
my @image_numbers = &Image::Data::Collection::gather_sorted_image_numbers(\%cfg);
my $image_num_length = length(scalar(@image_numbers));

print "\n\nBuild Movies\n\n" if $opt{debug};
$t1 = new Benchmark;

our @movie_dirs; 
find(\&add_to_movie_dir, (catdir($cfg{exp_results_folder},$cfg{movie_output_folder})));

foreach my $f1 (@movie_dirs) {
    foreach my $f2 (@{ $cfg{movie_output_prefix} }) {
        my $input_folder = catdir($f1,$f2);
        system "ffmpeg -v 0 -y -r $cfg{movie_frame_rate} -i $input_folder/%0" . $image_num_length . "d.png -sameq $input_folder.mov 2>&1";
    }
}
$t2 = new Benchmark;
print "Runtime: ",timestr(timediff($t2,$t1)), "\n" if $opt{debug};

###############################################################################
# Functions
###############################################################################
sub add_to_movie_dir {
    if ($File::Find::name =~ /$cfg{vis_config_file}/) {
        push @movie_dirs, $File::Find::dir;
    }
}
