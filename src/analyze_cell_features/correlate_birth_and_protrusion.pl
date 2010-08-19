#!/usr/bin/perl -w

###############################################################################
# This script collets the protrusion history prior to the birth of each focal
# adhesion in a time series.
#
# Input:
#   Birth matrix: The focal adhesions born in each frame of the time series.
#     The row index is the frame index and each row contains the indices into
#     the tracking matrix for the focal adhesions that were born in the
#     corresponding frame.
#   Tracking matrix: Local IDs of focal adhesions in each frame of the time
#     series. The row index is the global focal adhesion ID and the column
#     index is the frame index. The columns are the local IDs of that focal 
#     adhesion in that frame.
#   Focal adhesion ID/nearest edge pixel map: There is one of these for each
#     frame in the time series. It maps the local focal adhesion IDs in a
#     frame with the coordinates of the nearest edge pixel.
#   Protrusion marker pixels: There is one of set of these for each frame in
#     the series. These are the edge pixels used by the protrusion tracking 
#     system. 
#   Protrusion vectors: There is one set of these for each consecutive pair
#     of frames in the series. It maps each marker pixel at time T to the X
#     and Y components of the magnitude of the vector that translates that
#     pixel to its location at time T+1.
#
# Algorithm:
#   1. For each frame in the time series, get the tracking matrix indices of
#   the focal adhesions born in that frame.
#   2. The column and row in the birth matrix gives us the row and column in
#   the tracking matrix to find the local ID of a focal adhesion in its birth
#   frame.
#   3. Get the nearest edge pixel to the focal adhesion in the birth frame of
#   that focal adhesion by indexing into the map using the local ID obtained in
#   step 2.
#   4. Find the marker in the birth frame nearest to the pixel found in step 3.
#   5. Working backward, find the set of vectors that evolve to the marker 
#   found in step 4.
#
# Output: This script writes out two files:
#   1. Sequential protrusion history: Each output row is a focal adhesion (the 
#   indices are the same as those in the tracking matrix). The first number is 
#   the distance of the focal adhesion from its nearest edge pixel in the birth
#   frame. Each pair of numbers are the normal and parallel components of the 
#   vector to the direction of protrusion. If the protrusions could not be traced
#   all the way back to the first frame, pairs of zeros are used to pad the list
#   out to the first frame back to which the protrusion could be traced. For 
#   example, if a focal adhesion was born in frame 3, the output would look like:
#
#   37,2.1,8.7,-1.8,6.5,1.1,4.8,0.7,-2.1
#
#   This is interpreted as: The protrusion vector from frame 1->frame 2 has a normal 
#   component of 2.1 and a parallel component of 8.7. The protrusion vector from 
#   frame 2->frame 3 has a normal component of -1.8 and a parallel component of 6.5.
#   2. Relative protrusion history: Similar to the first file, but each row is
#   the relative protrusion growth in the current frame compared to the overall 
#   protrusion size, going backward from the birth frame. For example, given the data
#   above, the output would be:
#
#   37,.895,1.316,1.747
###############################################################################

###############################################################################
# Global Variables and Modules
###############################################################################
use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Path;
use File::Basename;
use File::Spec::Functions;
use Text::CSV;
use IO::File;
use Math::Trig;
use Image::Data::Collection;
use Getopt::Long;
use Text::CSV::Simple::Extra;

use Config::Adhesions qw(ParseConfig);

# Perl built-in variable that controls buffering print output, 1 turns off buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|config=s", "debug|d", "i|first=i", "n|images=i") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration\n" if $opt{debug};

my %cfg = &ParseConfig(\%opt);

my $csv=Text::CSV->new();

###############################################################################
#Main Program
###############################################################################

sub find_nearest_pixel_idx {
    my @test_pixel = @{ $_[0] };
    my @pixel_mat =  @{ $_[1] };
    my $col = $_[2]*2;
    my $nearest_idx;
    my $nearest_dist = 2 << 31;

    for (my $idx = 1; $idx < $#pixel_mat; $idx++) {
	my @row = @{ $pixel_mat[$idx] };
        my @pixel = @row[$col..$col+1];
        my $dist = sqrt( ($test_pixel[0]-$pixel[0])**2 + ($test_pixel[1]-$pixel[1])**2 );
        if ($dist < $nearest_dist) {
            $nearest_idx = $idx;
            $nearest_dist = $dist;
        }
    }

    return $nearest_idx;
}

sub normalize_vectors {
    my @vectors = @{ $_[0] };
    my $n_ang = $_[1];
    my $norm_factor = $_[2];
    my @normalized;

    for (my $i = 2; $i < scalar(@vectors); $i+=2) {
	next if $vectors[$i] < 0;
        my $vx = $vectors[$i]-$vectors[$i-2];
        my $vy = $vectors[$i+1]-$vectors[$i-1];
        my $v_mag = sqrt($vx**2 + $vy**2);

        my $v_nrm = 0;
        my $v_par = 0;

        if ($v_mag > 0) {
            my $v_ang = asin($vy/$v_mag);

            if ($vx < 0) {
                $v_ang = pi - $v_ang;
            }
            elsif ($vy < 0) {
                $v_ang = (2 * pi) + $v_ang;
            }
            
            my $a = $v_ang - $n_ang;
            $v_nrm = $v_mag * sin($a);
            $v_par = $v_mag * cos($a) / $norm_factor;
        }
        
#        push @normalized, ($v_nrm, $v_par);
        push @normalized, $v_par;
    }

    return @normalized;
}

# tracking matrix (each row is the ID of a single focal adhesion in each frame)
my $tracking_mat_file = catfile($cfg{exp_results_folder}, $cfg{tracking_folder}, $cfg{tracking_output_file});
my @tracking_mat = input_mat_csv($tracking_mat_file);

# protrusion vector history (each row is the vectors between each consecutive pairs of
# frames for a single marker point on the cell edge)
my $pr_vectors_file = catfile($cfg{exp_results_folder}, $cfg{protrusion_folder}, $cfg{protrusion_vectors_file});
my @pr_vectors = input_mat_csv($pr_vectors_file);

# birth matrix (each row contains the tracking matrix indices of the focal adhesions born 
# during a single frame)
my $birth_map_file = catfile($cfg{exp_results_folder}, $cfg{tracking_folder}, $cfg{birth_map_output_file});
my $birth_handle = new IO::File $birth_map_file || die "Can't open file $birth_map_file";

# the sorted image numbers actually used in this analysis
my @data_keys = Image::Data::Collection::gather_sorted_image_numbers(\%cfg);

my @sequential_history;
my @relative_history;

my $first_img = $opt{i} || 1;
my $num_frames = ($opt{n} || 0) - 1;
my $processed = 0;

for (my $frame = 0; <$birth_handle>; $frame++) {
    next if $frame < $first_img;
    last if $num_frames > 0 && $processed >= $num_frames;
    print("frame $frame\r");

    my $row = $_;

    # parse the current row of the birth matrix
    die "could not parse row $frame of $birth_map_file" if (!$csv->parse($row) || scalar($csv->fields) == 0);
    
    my @tracking_mat_indices = $csv->fields;

    # continue if there were no births in this frame
    next if $tracking_mat_indices[0] < 0;

    # get the image number for this frame index
    my $img_num = $data_keys[$frame];
    my $img_data_folder = catfile($cfg{individual_results_folder}, $img_num, $cfg{raw_data_folder});

    # get the nearest edge pixels for the focal adhesions in this frame
    my $centroid_file = catfile($img_data_folder, 'Centroid.csv');
    my $centroid_handle = new IO::File $centroid_file || die "Can't open file $centroid_file";
    die "Could not parse $centroid_file" if (!$csv->parse(<$centroid_handle>) || scalar($csv->fields) == 0);
    my @centroid_pixels = $csv->fields;
    $centroid_handle->close;

    # get the distance-to-edge for the focal adhesions in this frame
    my $dist_from_edge_file = catfile($img_data_folder, 'Centroid_dist_from_edge.csv');
    my $dfe_handle = new IO::File $dist_from_edge_file || die "Can't open file $dist_from_edge_file";
    die "Could not parse $dist_from_edge_file" if (!$csv->parse(<$dfe_handle>) || scalar($csv->fields) == 0);
    my @dists_from_edge = $csv->fields;
    $dfe_handle->close;

    # get the angle of position of focal adhesions relative to the 
    # center of the cell in this frame
    my $angle_to_center_file = catfile($img_data_folder, 'Angle_to_center.csv');
    my $atc_handle = new IO::File $angle_to_center_file || die "Can't open file $angle_to_center_file";
    die "Could not parse $angle_to_center_file" if (!$csv->parse(<$atc_handle>) || scalar($csv->fields) == 0);
    my @angles_to_center = $csv->fields;
    $atc_handle->close;

    # iterate through the tracking matrix indices of the focal 
    # adhesions born in this frame
    foreach my $tracking_mat_index (@tracking_mat_indices) {
        my $local_fa_id = $tracking_mat[$tracking_mat_index][$frame];
	my $centroid_idx = $local_fa_id * 2;
	my @centroid_pixel = @centroid_pixels[$centroid_idx..$centroid_idx+1];
        my $dist_from_edge = $dists_from_edge[$local_fa_id];
        my $angle_to_center = $angles_to_center[$local_fa_id];
        # find the protrusion vector origin closest to the edge
        # pixel closest to the FA. 
	my $pr_vectors_idx = $frame - $first_img + 1;
        my $pr_pixel_idx = find_nearest_pixel_idx(\@centroid_pixel, \@pr_vectors, $pr_vectors_idx);
        # get the sequence of protrusion vectors from the first
        # frame to the current frame
	my @pr_vector_row = @{ $pr_vectors[$pr_pixel_idx] };
        my @pr_vector_history = @pr_vector_row[0..($pr_vectors_idx*2)+1];

        # get only the components of the vectors parallel to this
        # FA's angle from the center of the cell
        my @normalized_pr_vectors = normalize_vectors(\@pr_vector_history, $angle_to_center, $dist_from_edge);
#        my @seq = ($dist_from_edge);
#        push(@seq, @normalized_pr_vectors);
        push(@sequential_history, [@normalized_pr_vectors]);

        # calculate the growth/shrinkage of each vector relative
        # to the total size of the protrusion
	my $norm_vec_size = scalar(@normalized_pr_vectors);
	if ($norm_vec_size > 0) {
	        my @rel;
        	my $total_mag = 0;
	        for (my $i = 0; $i < scalar @normalized_pr_vectors; $i++) {
        	    my $vec_mag = $normalized_pr_vectors[$i];
            	    if ($i == 0) {
                	$total_mag = $vec_mag;
            	    }
	            else {
        	        my $growth = 0;
			if ($total_mag != 0) {
				$growth = (($total_mag + $vec_mag) / $total_mag) - 1;
			}
        	        $total_mag += $vec_mag;
                	push(@rel, $growth);
	            }
        	}
#        my @temp = ($dist_from_edge);
		my @rev = reverse(@rel);
#        push(@temp, @rev);
	        push(@relative_history, [@rev]);   
   	 }
    }

    $processed++;
}

$birth_handle->close;

my $sequential_history_output_file = catfile($cfg{exp_results_folder}, $cfg{protrusion_folder}, $cfg{sequential_history_output_file});
output_mat_csv(\@sequential_history, $sequential_history_output_file);

my $relative_history_output_file = catfile($cfg{exp_results_folder}, $cfg{protrusion_folder}, $cfg{relative_history_output_file});
output_mat_csv(\@relative_history, $relative_history_output_file);
