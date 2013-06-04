#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################
use lib "../lib";

use strict;
use File::Path;
use File::Basename;
use File::Find;
use Image::ExifTool;
use Math::Matlab::Local;
use Getopt::Long;
use Data::Dumper;
use File::Spec::Functions;
use Benchmark;
use POSIX;
use Text::CSV::Simple;

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
GetOptions(\%opt, "cfg=s", "debug|d") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my %cfg = ParseConfig(\%opt);

###############################################################################
#Main Program
###############################################################################

my $base_dir = catdir($cfg{exp_results_folder}, $cfg{movie_output_folder}, 'single_ad_files');

my $color_file = catfile($base_dir, 'c_map.csv');
my $parser = Text::CSV::Simple->new;
my @color_map = $parser->read_file($color_file);
my @hex_color_map = &convert_rbg_to_hex(@color_map);

my $max_image = &find_max_image_num;

my @full_svg_file;
my @svg_header;
for my $i_num (1 .. $max_image) {
    my ($files_ref, $montage_num_ref) = &get_image_data($i_num);

    my @bmp_files = &build_bitmaps(@{$files_ref});
    my @svg_files = &convert_using_potrace(\@bmp_files,$montage_num_ref,\@hex_color_map);
    
    if (scalar(@svg_header) == 0) {
        @svg_header = &get_svg_header($svg_files[0]);
    }

    push @full_svg_file, &build_full_svg_file(@svg_files);
    
    my $output_file = catfile($cfg{exp_results_folder}, $cfg{movie_output_folder}, "ghost_unique.svg");
    open SVG_OUT, ">$output_file";
    print SVG_OUT @svg_header;
    print SVG_OUT @full_svg_file;
    print SVG_OUT "</g>\n";
    print SVG_OUT '</svg>';
    close SVG_OUT;

    unlink(@bmp_files);
    unlink(@svg_files);

    print $i_num, " " if $opt{debug};
}

my $output_file = catfile($cfg{exp_results_folder}, $cfg{movie_output_folder}, "ghost_unique.svg");
open SVG_OUT, ">$output_file";
print SVG_OUT @svg_header;
print SVG_OUT @full_svg_file;
print SVG_OUT "</g>\n";
print SVG_OUT '</svg>';
close SVG_OUT;

#my $output_png = catfile($cfg{exp_results_folder}, $cfg{movie_output_folder}, "ghost.png");
#my $output_png_small = catfile($cfg{exp_results_folder}, $cfg{movie_output_folder}, "ghost_small.png");
#
#system "convert -density 100x100 $output_file $output_png_small";
#system "convert -density 300x300 $output_file $output_png";
#

###############################################################################
#Functions
###############################################################################
sub find_max_image_num {
    my $max_num = 0;
    for (<$base_dir/*/*.png>) {
        $_ =~ m#$base_dir/\d+/(\d+)\.png#;
        if ($1 > $max_num) {
            $max_num = $1;
        }
    }
    return $max_num;
}

sub get_image_data {
    my $i_num = $_[0];
    my @files;
    my @montage_nums;
    for (<$base_dir/*/*.png>) {
        $_ =~ m#$base_dir/(\d+)/(\d+)\.png#;
        if ($i_num == $2) {
            push @files, $_;
            push @montage_nums, $1;
        }
    }
    return (\@files, \@montage_nums);
}

sub build_bitmaps {
    my @files = @_;
    my @bmp_files;
    for (@files) {
        my $bmp_file = $_;
        $bmp_file =~ s/\.png/\.bmp/;
        push @bmp_files, $bmp_file;
        system "convert $_ $bmp_file";
    }
    return @bmp_files;
}

sub convert_using_potrace {
    my @files = @{$_[0]};
    my @montage_num = @{$_[1]};
    my @colors = @{$_[2]};
    my @svg_files;
    for (0 .. $#files) {
        my $bmp_file = $files[$_];
        my $color = $colors[$montage_num[$_] - 1];
        my $svg_file = $bmp_file;
        
        $svg_file =~ s/\.bmp/\.svg/;
        push @svg_files, $svg_file;
        system "potrace -t 0 -s --fillcolor=#$color $bmp_file\n";
    }
    return @svg_files;
}

sub get_svg_header {
    my $file = $_[0];

    my @header;
    my $path_count = 0;

    open SVG_FILE, $file;
    while (<SVG_FILE>) {
        if (/\<path/) {
            $path_count++;
        }
        if ($path_count < 2) {
            push @header, $_;
        }
    }
    close SVG_FILE;

    return @header;
}

sub build_full_svg_file {
    my @svg_files = @_;
    
    my @svg_path_data;

    for (@svg_files) {
        my $path_count = 0;
        my $out_of_path = 0;
        open SVG_FILE, $_;
        for my $line (<SVG_FILE>) {
            if ($line =~ /\<path/) {
                $path_count++;
            }
            if ($line =~ /\<\/g\>/) {
                $out_of_path = 1;
            }
            if ($path_count > 1 && not($out_of_path)) {
                push @svg_path_data, $line;
            }
        }
        close SVG_FILE;
    }

    return @svg_path_data;
}

sub convert_rbg_to_hex {
    my @colors = @_;
    
    my @hex_colors;
    
    for (@colors) {
        my @this_color = @{$_};
        my $in_hex;
        for (@this_color) {
            my $hex_num = sprintf('%x',round($_ * 255));
            if (length($hex_num) == 1) {
                $hex_num = "0" . $hex_num;
            }
            $in_hex .= $hex_num;
        }
        push @hex_colors, $in_hex;
    }
    return @hex_colors;
}

sub round {
    my($number) = shift;
    return int($number + .5);
}

###############################################################################
#Documentation
###############################################################################

=head1 NAME

collect_visualizations.pl - build the visualizations of the focal adhesion
movies

=head1 SYNOPSIS

collect_mask_set.pl -cfg FA_config

=head1 DESCRIPTION

This program builds a series of movies based on files available in the tracking
matrices folder. Each file in the tracking matrices folder which ends with
'.csv' and does not contain 'no_movie' is used to build a visualization of the
tracked focal adhesions.

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=back

Optional parameter(s):

=over 

=item * debug or d: print debuging information during program execution

=item * movie_debug: pass along the debug flag to the MATLAB visualization
program, causing only a small subset of the tracked adhesions to be visualized
in a single frame

=item * config_only: only write the MATLAB config files out, do not execute the
MATLAB program

=item * emerald: build and execute long commands throught the LSF job system

=back

=head1 EXAMPLES

collect_visualizations.pl -cfg FA_config

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 6/30/2008

=cut
