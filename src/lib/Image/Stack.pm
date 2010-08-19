#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################

use strict;
use warnings;
use Image::ExifTool;
use Data::Dumper;

package Image::Stack;

###############################################################################
# Functions
###############################################################################

sub get_image_stack_number {
    my $image_file = shift;
    my $image_info = new Image::ExifTool;
    $image_info->ExtractInfo($image_file);
    my @tag_list    = $image_info->GetFoundTags($image_file);
    my $image_count = 0;
    foreach (@tag_list) {
        if (/\((\d+)\)/) {
            if ($1 > $image_count) {
                $image_count = $1;
            }
        }
    }

    $image_count++;
    return $image_count;
}

1;
