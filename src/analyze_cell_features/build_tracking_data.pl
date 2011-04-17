#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Path;
use File::Basename;
use File::Spec::Functions;
use Getopt::Long;
use Data::Dumper;
use Storable;
use Text::CSV;
use IO::File;
use Math::Matrix;

use Config::Adhesions;
use Image::Data::Collection;
use Text::CSV::Simple::Extra;
use Emerald;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
$opt{output} = "data.stor";
GetOptions(\%opt, "cfg|config=s", "debug|d", "output|o=s", "image_num=s", 
                  "lsf|l") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration\n" if $opt{debug};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg = $ad_conf->get_cfg_hash;

###############################################################################
#Main Program
###############################################################################
my %data_sets;
if ($opt{lsf}) {
    my @image_nums = &Image::Data::Collection::gather_sorted_image_numbers(\%cfg);

    my @commands;
    foreach (@image_nums) {
        #$0 - the name of the program currently running, used to protect against
        #future file name changes
        push @commands, "$0 -cfg $opt{cfg} -o $opt{output} -image_num $_";
    }
    
    $opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'tracking_data');
    $opt{resource} = "blade RH5";
    $opt{queue} = "week";
	#mem specified in kb
    $opt{mem} = (1024**2)*4;
    if (defined $cfg{job_group}) {
        $opt{job_group} = $cfg{job_group};
    }
    
    &FA_job::send_general_lsf_program(\@commands,\%opt);

    exit(0);
} else {
    print "\n\nGathering Data Files\n" if $opt{debug};

    my @data_files;
    push @data_files, @{$cfg{general_data_files}};
    push @data_files, @{$cfg{tracking_files}};

    %data_sets = &Image::Data::Collection::gather_data_sets(\%cfg, \%opt, \@data_files);
    
	print "\n\nMaking Comparison Matrices\n" if $opt{debug};
    &make_comp_matices;
}

###############################################################################
#Functions
###############################################################################

#######################################
# Process Data Sets
#######################################
sub make_comp_matices {
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for (0 .. $#data_keys) {
        if (exists $opt{image_num} && $data_keys[$_] != $opt{image_num}) {
            next;
        }

        #The last image can not be compared to a future image, so we skip
        #calculations on it, but still save the image data if the output
        #option is specified
        if ($_ == $#data_keys) {
            if (defined $opt{output}) {
                my $key_1 = $data_keys[$_];
                store \%{ $data_sets{$key_1} }, catfile($cfg{individual_results_folder}, $key_1, $opt{output});
                delete $data_sets{$key_1};
            }
            next;
        }

        #These are the keys we will use for all the subsequent matrix creation
        my ($key_1, $key_2) = @data_keys[ $_, $_ + 1 ];
        
        print "Working on image Number: $key_1 - " if $opt{debug};

        #Gather the Centroid distance matrix
        my @x1 = @{ $data_sets{$key_1}{Centroid_x} };
        my @y1 = @{ $data_sets{$key_1}{Centroid_y} };
        my @x2 = @{ $data_sets{$key_2}{Centroid_x} };
        my @y2 = @{ $data_sets{$key_2}{Centroid_y} };
        @{ $data_sets{$key_1}{Cent_dist} } = &make_dist_mat(\@x1, \@y1, \@x2, \@y2);
        print "Cent_dist Collected - " if $opt{debug};
		
        #Gather the Pixel Similarity matrix
        my @pix_id1 = @{ $data_sets{$key_1}{PixelIdxList} };
        my @pix_id2 = @{ $data_sets{$key_2}{PixelIdxList} };
        @{ $data_sets{$key_1}{Pix_sim} } = &calc_pix_sim(\@pix_id1, \@pix_id2, $data_sets{$key_1}{Cent_dist});
        @{ $data_sets{$key_1}{Recip_pix_sim} } = &calc_pix_sim(\@pix_id2, \@pix_id1);
       
        if ($opt{debug}) {
            @{ $data_sets{$key_1}{Pix_sim_f} } = &calc_pix_sim(\@pix_id1, \@pix_id2);
            my $mat_1 = new Math::Matrix @{$data_sets{$key_1}{Pix_sim}};
            my $mat_2 = new Math::Matrix @{$data_sets{$key_1}{Pix_sim_f}};
            die "Problem with new pixel sim calc method\n", join(" ", $mat_1->size), 
                "  ",join(" ", $mat_2->size) if not $mat_1->equal($mat_2);
        }


        delete $data_sets{$key_1}{PixelIdxList};
        print "Pix_sim Collected" if $opt{debug};
        print "\r"                if $opt{debug};

        if (defined $opt{output}) {
            store \%{ $data_sets{$key_1} }, catfile($cfg{individual_results_folder}, $key_1, $opt{output});
            delete $data_sets{$key_1};
        }
    }
}

sub make_dist_mat {
    my ($ref_1, $ref_2, $ref_3, $ref_4) = @_;
    my @x1 = @$ref_1;
    my @y1 = @$ref_2;
    my @x2 = @$ref_3;
    my @y2 = @$ref_4;

    my @dist_mat;

    for my $i (0 .. $#x1) {
        for my $j (0 .. $#x2) {
            $dist_mat[$i][$j] = sqrt(($x1[$i] - $x2[$j])**2 + ($y1[$i] - $y2[$j])**2);
        }
    }
    return @dist_mat;
}

sub calc_pix_sim {
    my @pix_id1 = @{ $_[0] };
    my @pix_id2 = @{ $_[1] };
  	my @cent_dists;
    if (scalar(@_) > 2) {
    	@cent_dists = @{ $_[2] };
    }

    my @sim_percents;
    for my $i (0 .. $#pix_id1) {
        our @current_pix_list = @{ $pix_id1[$i] };
        my $current_pix_list_length = scalar(@current_pix_list);
        
        die "\n\nProblem with pixel ID list ($i), the length is zero" if not $current_pix_list_length;
        
        my @search_order;
        if (@cent_dists) {
            my @dist_to_next_ads = @{$cent_dists[$i]};
            @search_order = sort {$dist_to_next_ads[$a] <=> $dist_to_next_ads[$b]} (0 .. $#dist_to_next_ads);
        } else {
            @search_order = (0 .. $#pix_id2);
        }

        for my $j (@search_order) {
            die "\n\nProblem with the number of entries in the second pixel " . 
				"ID list index ($j), the whole list is:", Dumper(\@pix_id2) 
				if (not($pix_id2[$j]));
            my @next_pix_list = @{ $pix_id2[$j] };
            my $match_count = 0;
            for (0 .. $#current_pix_list) {
                my $poss_match = pop @current_pix_list;
                my $temp_match_count = grep $poss_match == $_, @next_pix_list;
                if ($temp_match_count == 0) {
                    unshift @current_pix_list, $poss_match;
                } else {
                    $match_count += $temp_match_count;
                }
            }
            $sim_percents[$i][$j] = $match_count / $current_pix_list_length;
        }
    }
    return @sim_percents;
}
