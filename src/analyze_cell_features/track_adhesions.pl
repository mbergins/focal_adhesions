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

use Config::Adhesions qw(ParseConfig);
use Image::Data::Collection;
use Text::CSV::Simple::Extra;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
$opt{input} = "data.stor";
$opt{keep_data_files} = 0;
GetOptions(\%opt, "cfg|config=s", "debug|d", "input|i|input_data_files=s", 
                  "lsf|l", "keep_data_files") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration\n" if $opt{debug};

my %cfg = &ParseConfig(\%opt);

###############################################################################
#Main Program
###############################################################################

if ($opt{lsf}) {
    my $extra = "";
    $extra .= " -keep_data_files" if $opt{keep_data_files};
    my @command = "$0 -cfg $opt{cfg} -input $opt{input} $extra";
    $opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'tracking');    
    if (defined $cfg{job_group}) {
        $opt{job_group} = $cfg{job_group};
    }
    
    &FA_job::send_general_lsf_program(\@command,\%opt);
    
    exit;
}

my %data_sets;

print "\n\nDetermining Tracking Matrix\n" if $opt{debug};
my @tracking_mat;
my %tracking_probs;
my %tracking_facts;
&make_tracking_mat;

#&output_mat_csv(\@{$tracking_facts{"no_pix_sim_dist"}}, "no_sim_dist.csv");
#&output_mat_csv(\@{$tracking_facts{"no_pix_sim_ratio"}}, "no_sim_ratio.csv");

print "\n\nTracking Results\n" if $opt{debug};
#&output_tracking_facts         if $opt{debug};
print "\n"                     if $opt{debug};

print "\n\nOutputing Tracking Problem Data\n" if $opt{debug};
#&output_tracking_probs;

print "\n\nOutputing Tracking Matrix\n" if $opt{debug};
my $tracking_output_file = catfile($cfg{exp_results_folder}, $cfg{tracking_folder}, $cfg{tracking_output_file});
mkpath(dirname($tracking_output_file));
&output_mat_csv(\@tracking_mat, $tracking_output_file);

###############################################################################
#Functions
###############################################################################

#######################################
# Tracking Matrix Production
#######################################
sub make_tracking_mat {
    my @data_keys;
    if (scalar(keys %data_sets) == 0) {
        @data_keys = Image::Data::Collection::gather_sorted_image_numbers(\%cfg);
    } else {
        @data_keys = sort { $a <=> $b } keys %data_sets;
    }
	
	#Stop operation when there is only one image set, can't do any tracking in
	#that case
	if (scalar(@data_keys) == 1) {
        if (not($opt{keep_data_files})) {
            unlink(catfile($cfg{individual_results_folder}, $data_keys[0], $opt{input}));
        }
		print "Found only one image set for tracking, quitting." if $opt{debug};
		exit;
	}
	
    #This loop contains all the commands used to create the tracking matrix.
    #During each loop another time point is added to all lineages. There are
    #four discrete steps in this process, each one addressing a specific
    #operation in adhesion model, either birth, merging/death or survival.
    #
    #	Tracking Matrix Building Operations:
    #
    #	1. Initialize the tracking matrix with all the adhesion numbers present
    #	   in the first image number (initial births)
    #	2. Try to find a match for all the live adhesions in the next frame of
    #	   the image set (survival)
    #	3. Resolve the adhesion lineages that appear to merge, either by
    #	   picking a lineage as the best match or declaring the other lineages
    #	   dead (death/merging)
    #	4. Find any adhesions in the next frame not assigned to a lineage and
    #	   form a new lineage for each adhesion (birth)
    #	5. Repeat 2-5 until the data from each image has been used to build the
    #	   tracking matrix
    #
    #This loop also builds a birth matrix. Each row in the matrix is a frame in
    #the movie (data set) and contains the numbers of the adhesions (i.e. 
    #indices of the tracking matrix) that were born in that frame.

    for (0 .. $#data_keys - 1) {
        my $i_num      = $data_keys[$_];
        my $next_i_num = $data_keys[ $_ + 1 ];
        
        #Collect data seta and comparison matrices if not present
        if (not(defined $data_sets{$i_num})) {
            if (defined $opt{input}) {
                %{ $data_sets{$i_num} } = %{ retrieve(catfile($cfg{individual_results_folder}, $i_num, $opt{input})) };
            } else {
                die "Unable to find the data sets and comparison matrices for image number \"$i_num\".";
            }
        }
        if (not(defined $data_sets{$next_i_num})) {
            if (defined $opt{input}) {
                %{ $data_sets{$next_i_num} } =
                  %{ retrieve catfile($cfg{individual_results_folder}, $next_i_num, $opt{input}) };
            } else {
                die "Unable to find the data sets and comparison matrices for image number \"$next_i_num\".";
            }
        }
        
        #STEP 1
        #Start the tracking matrix, if this is the first time throught the loop
        &initialize_tracking_mat($data_keys[0]) if $_ == 0;
        &check_tracking_mat_integrity if $_ == 0;

        #Begin tracking
        print "\r", " "x80, "\rImage #: $i_num - " if $opt{debug};

        #STEP 2
        &track_live_adhesions($i_num);
        print "# Tracked: $tracking_facts{$i_num}{live_adhesions} - " if $opt{debug};

        #STEP 3
        &detect_merged_adhesions($i_num, $next_i_num);
        print "# Merged: $tracking_facts{$i_num}{merged_count}/$tracking_facts{$i_num}{merged_prob_count} - "
          if $opt{debug};

        #STEP 4
        &detect_new_adhesions($i_num, $next_i_num);
        print "# New: $tracking_facts{$i_num}{new_count}" if $opt{debug};

        &check_tracking_mat_integrity;

        delete $data_sets{$i_num};
        if (not($opt{keep_data_files})) {
            unlink(catfile($cfg{individual_results_folder}, $i_num, $opt{input}));
        }

        #Remove the last image's data file, if keep_data_files is not true
        if (not($opt{keep_data_files}) && ($_ == ($#data_keys - 1))) { 
            unlink(catfile($cfg{individual_results_folder}, $next_i_num, $opt{input}));
        }
    }
}

sub initialize_tracking_mat {
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    my $first_key = $data_keys[0];
    my $num_adhesions = $#{ $data_sets{$first_key}{Area} };

    for (0 .. $num_adhesions) {
        push @{ $tracking_mat[$_] }, $_;
    }
    print "Intialized the tracking matrix on image #: $first_key\n" if $opt{debug};
}

sub track_live_adhesions {
    my ($i_num) = @_;
    my $cur_step = $#{ $tracking_mat[0] };

    my @dists        = @{ $data_sets{$i_num}{Cent_dist} };
    my @p_sims       = @{ $data_sets{$i_num}{Pix_sim} };
    my @recip_p_sims = @{ $data_sets{$i_num}{Recip_pix_sim} };
    my @areas        = @{ $data_sets{$i_num}{Area} };
    
    my $prop_indeter_percent = 0.8;
    $prop_indeter_percent = $cfg{prop_indeter_percent} if defined $cfg{prop_indeter_percent};

    #This function makes a tracking guess for all the adhesion lineages living
    #in the current image. The tracking guess is based on two factors: the
    #distance to the adhesions in the next image and the percentage of pixels
    #that overlap with the each of the adhesions in the next image (pixel
    #similarity). There are several cases we need to deal with in making the
    #tracking guess:
    #
    #    Live Tracking Cases:
    #
    #    1. There is a single pixel similarity value that is much higher than
    #    the other values, select the adhesion with the best similarity,
    #    assuming that the reciprocal pixel similarity is close as well
    #
    #    2. There are multiple adhesions with close pixel similarity value to
    #    the current adhesion, select the adhesion with the smallest distance
    #    between the centriods
    #
    #    3. There is no pixel similarity to the next adhesion set, select the
    #    adehsion with the smallest distance between the centriods
    #
    #    Notes:
    #
    #    -how close the similarity measures have to be to trigger using centroid
    #    distance between adhesions is defined in $prop_indeter_percent, which
    #    defaults to 0.8

    for my $i (0 .. $#tracking_mat) {
        
        #The tracking matrix code for dead lineage is any number less than or
        #equal to -1, add another -1 to those lineages to make sure they stay
        #the proper length, then skip to the next lineage
        if ($tracking_mat[$i][$cur_step] <= -1) {
            push @{ $tracking_mat[$i] }, -1;
            next;
        }
        
        $tracking_facts{$i_num}{live_adhesions}++;
        my $adhesion_num      = ${ $tracking_mat[$i] }[$cur_step];
        my @dist_to_next_ads  = @{ $dists[$adhesion_num] };
        my @p_sim_to_next_ads = @{ $p_sims[$adhesion_num] };
        my @recip_p_sims      = map $recip_p_sims[$_][$adhesion_num], (0 .. $#dist_to_next_ads);

        my @sorted_dist_indexes =
          sort { $dist_to_next_ads[$a] <=> $dist_to_next_ads[$b] } (0 .. $#dist_to_next_ads);

        my @sorted_p_sim_indexes =
          sort { $p_sim_to_next_ads[$b] <=> $p_sim_to_next_ads[$a] } (0 .. $#p_sim_to_next_ads);

        my $high_p_sim = $p_sim_to_next_ads[ $sorted_p_sim_indexes[0] ];

        if ($sorted_p_sim_indexes[0] != $sorted_dist_indexes[0] && $high_p_sim > 0) {
            $tracking_facts{$i_num}{dist_p_sim_guess_diff}++;
        }

        my $tracking_guess;

        #Cases 1 and 2
        #Check if the highest pixel sim (p_sim) is above zero, indicating that
        #this adhesion overlaps with at least one adhesion in the next image
        if ($high_p_sim > 0) {

            #Find how many adhesions in the next image overlap with this
            #adhesion
            my @p_sim_close_ad_nums = grep {
                if ($p_sim_to_next_ads[$_] >= $high_p_sim * $prop_indeter_percent) {
                    1;
                } else {
                    0;
                }
            } (0 .. $#p_sim_to_next_ads);
            
            #Of the adhesions with close pixel similarities, sort them by the
            #centroid distances, most of the time, this matrix will only contain
            #one value, either way, the first entry in this matrix will always
            #be our winner
            my @close_p_sim_by_dist = sort { $dist_to_next_ads[$a] <=> $dist_to_next_ads[$b] } @p_sim_close_ad_nums;
            
            $tracking_guess = $close_p_sim_by_dist[0];

            #Now for some additional data collection: detecting split birth
            #events and cases where the adhesion with the best pixel similarity
            #was not selected
            
            #Find adhesions that are close in pixel similiarity to the winning
            #adhesion, but not close enough to trigger a tracking decision based
            #on distance. As the variable name implies it will be used to find
            #split birth events.
            my @split_birth_p_sim_close = grep {
                if ($p_sim_to_next_ads[$_] > 0) {
                    1;
                } else {
                    0;
                }
            } (0 .. $#p_sim_to_next_ads);
            my $split_length_before = scalar(@split_birth_p_sim_close);

            #Remove the winning adhesion from the split birth distance list, as
            #it should not be included in split birth decisions
            @split_birth_p_sim_close = map {
                if ($close_p_sim_by_dist[0] != $_) {
                    $_;
                }
            } @split_birth_p_sim_close;
            my $split_length_after = scalar(@split_birth_p_sim_close);
            
            die "Problem with removing winning adhesion from split birth list." 
              if (($split_length_after + 1) == $split_length_before);
            
            foreach my $ad_num (@split_birth_p_sim_close) {
                if (exists $tracking_facts{$i_num}{split_birth_overlap}[$ad_num]) {
                    if ($tracking_facts{$i_num}{split_birth_overlap}[$ad_num] < $p_sim_to_next_ads[$ad_num]) {
                        $tracking_facts{$i_num}{split_birth}[$ad_num] = $tracking_guess;
                        $tracking_facts{$i_num}{split_birth_overlap}[$ad_num] = $p_sim_to_next_ads[$ad_num];
                        $tracking_facts{$i_num}{split_birth_ness}[$ad_num] = $p_sim_to_next_ads[$ad_num]/$high_p_sim;
                    }
                } else {
                    $tracking_facts{$i_num}{split_birth}[$ad_num] = $tracking_guess;
                    $tracking_facts{$i_num}{split_birth_overlap}[$ad_num] = $p_sim_to_next_ads[$ad_num];
                    $tracking_facts{$i_num}{split_birth_ness}[$ad_num] = $p_sim_to_next_ads[$ad_num]/$high_p_sim;
                }
            }

            if ($close_p_sim_by_dist[0] != $sorted_p_sim_indexes[0]) {
                $tracking_facts{$i_num}{best_pix_sim_not_selected}++;
            }

            if (scalar(@close_p_sim_by_dist) > 0) {
                $tracking_facts{$i_num}{multiple_good_p_sims}++;
            }
        } else {

            #Case 3
            $tracking_guess = $sorted_dist_indexes[0];

            if ($dist_to_next_ads[$tracking_guess]/$areas[$adhesion_num] >= 5) {
                $tracking_guess = -1;
            } else {
                push @{$tracking_facts{"no_pix_sim_dist"}}, $dist_to_next_ads[$tracking_guess];
                push @{$tracking_facts{"no_pix_sim_ratio"}}, $dist_to_next_ads[$tracking_guess]/$areas[$adhesion_num];
            }
        }

        push @{ $tracking_mat[$i] }, $tracking_guess;
    }
}

sub detect_merged_adhesions {
    my ($i_num, $next_i_num) = @_;

    my $num_ad_lineages = $#tracking_mat;
    my $cur_step        = $#{ $tracking_mat[0] };
    my @areas           = @{ $data_sets{$i_num}{Area} };
    my @next_areas      = @{ $data_sets{$next_i_num}{Area} };
    my @dists           = @{ $data_sets{$i_num}{Cent_dist} };
    my @pix_sims        = @{ $data_sets{$i_num}{Pix_sim} };

    # %dest_adhesions will act as the lookup table for adhesions predicted to
    # merge
    my %dest_adhesions;
    for my $i (0 .. $num_ad_lineages) {
        next if $tracking_mat[$i][$cur_step] <= -1;

        my $this_ad = $tracking_mat[$i][$cur_step];
        push @{ $dest_adhesions{$this_ad}{lineage_nums} }, $i;
        push @{ $dest_adhesions{$this_ad}{starting_ad} },  $tracking_mat[$i][ $cur_step - 1 ];
        push @{ $dest_adhesions{$this_ad}{ending_ad} },    $this_ad;
        push @{ $dest_adhesions{$this_ad}{starting_lifetime} }, &determine_ahesion_lifetime(@{$tracking_mat[$i]});
    }

    $tracking_facts{$i_num}{merged_count}      = 0;
    $tracking_facts{$i_num}{merged_prob_count} = 0;
    for my $i (keys %dest_adhesions) {
        next if scalar(@{ $dest_adhesions{$i}{lineage_nums} }) == 1;

        $tracking_facts{$i_num}{merged_count}++;

        my @lineage_nums = @{ $dest_adhesions{$i}{lineage_nums} };
        my @starting_ads = @{ $dest_adhesions{$i}{starting_ad} };
        my @ending_ad    = @{ $dest_adhesions{$i}{ending_ad} };
        for my $j (1 .. $#ending_ad) {
            die "Problem with merge detection, attempted to merge with different predicted end adhesions." . 
                join(" ",@ending_ad) if ($ending_ad[0] != $ending_ad[$j]);
        }

        my @merged_areas = @areas[@starting_ads];
        my @dist_shifts  = map { $dists[ $starting_ads[$_] ][ $ending_ad[$_] ] } (0 .. $#lineage_nums);
        my @pix_sims = map { $pix_sims[ $starting_ads[$_] ][ $ending_ad[$_] ] } (0 .. $#lineage_nums);
        my $ending_area = @next_areas[$ending_ad[0]];
        my @lifetimes    = @{ $dest_adhesions{$i}{starting_lifetime} };

        my @merge_decisions = &select_best_merge_decision(\@merged_areas, \@dist_shifts, \@pix_sims, \$ending_area, \@lifetimes);
        
        foreach (@merge_decisions) {
            if (${$_}{case} != 0) {
                my $case_str = "merge_case_" . ${$_}{case};
                $tracking_facts{$i_num}{$case_str}++;
            }
        }
        
        my $all_dead = 1;
        foreach (0 .. $#lineage_nums) {
            $all_dead = 0 if (not($merge_decisions[$_]{dead}));
        }
        foreach (0 .. $#lineage_nums) {
            if ($merge_decisions[$_]{dead}) {
                $tracking_mat[ $lineage_nums[$_] ][$cur_step] = -1;
            } elsif (not $merge_decisions[$_]{winner}) {
                $tracking_mat[ $lineage_nums[$_] ][$cur_step] =
                  -1 * ($tracking_mat[ $lineage_nums[$_] ][$cur_step] + 2);
            }
            
            if ($merge_decisions[$_]{winner} && $merge_decisions[$_]{dead} && not($all_dead)) {
                $tracking_facts{$i_num}{dead_winners}++;
            }
        }
    }
}

sub select_best_merge_decision {
    my @merged_areas = @{ $_[0] };
    my @dist_shifts  = @{ $_[1] };
    my @pix_sims = @{ $_[2] };
    my $ending_area = ${$_[3]};
    my @lifetimes = @{$_[4]};

    #There are several cases to deal with in picking the adhesion which will
    #continue. Also note that "close" will be defined as within the percentage
    #specified by "$prop_indeter_percent"
    #    
    # 1. Find the adhesion whose pixels overlap the greatest percentage of the
    # pixels in the merged adhesion, if there are multiple adhesions with close
    # overlap percentages, try case 2, otherwise, select the adhesion with the
    # greatest overlap
    #
    # 2. Find the adhesion from the high overlap percentage with the greatest
    # area, if the areas are close, try case 3, otherwise, select the adhesion
    # with the greatest area from the high overlap percentages
    #
    # 3. Find the adhesion from the adhesions with close areas and high overlap
    # percentages and pick the adhesion whose centroid is the shortest distance
    # from the merged adhesion's centroid, if more than one adhesion is close,
    # try case 4
    #
    # 4. Find the adhesion with the longest life from the adhesions with close
    # areas, high overlap percentages and similar centroid distances, pick one
    # adhesion regardless of how close the lifetimes are
    #
    # 5. The lifetimes used in step 5 are equal and I don't know what else to
    # use for deciding the merge decision, select the adhesion that happens to
    # come up first 
    
    #Merge decision data structure, holds which adhesion wins the merge event,
    #which adhesions die and which case from the above comment made the
    #decision
    my @merge_decisions = map { { winner => 0, dead => 0, case => 0 } } (0 .. $#merged_areas);
    
    ###################################
    #Property Calculations   
    ###################################
    my @percent_end_overlap = map {($pix_sims[$_] * $merged_areas[$_])/$ending_area } (0 .. $#merged_areas);
    my $sum = 0;
    for (0 .. $#percent_end_overlap) {
        $sum += $percent_end_overlap[$_];
    }
    die "\nProblem with determining percent of each merging adhesion that " . 
	  "overlaps with the merged adhesion." if ($sum > 1.01);
    my $high_overlap = (sort {$b <=> $a} (@percent_end_overlap))[0];

    my $prop_indeter_percent = 0.8;
    $prop_indeter_percent = $cfg{prop_indeter_percent} if defined $cfg{prop_indeter_percent};
    
    my @overlap_close = grep {
        if ($percent_end_overlap[$_] >= $high_overlap * $prop_indeter_percent) {
            1;
        } else {
            0;
        }
    } (0 .. $#merged_areas);
    die "Merge problem: no overlap close results." if scalar(@overlap_close) == 0;
    
    my $high_area = (sort {$b <=> $a} (@merged_areas[@overlap_close]))[0];
    my @areas_close = grep {
        if ($merged_areas[$_] >= $high_area * $prop_indeter_percent) {
            1;
        } else {
            0;
        }
    } (@overlap_close);
    die "Merge problem: no area close results." if scalar(@areas_close) == 0;
    
    my $min_dist = (sort {$a <=> $b} @dist_shifts[@areas_close])[0];
    my @dists_close = grep {
        if ($dist_shifts[$_] <= $min_dist * (2 - $prop_indeter_percent)) {
            1;
        } else {
            0;
        }
    } (@areas_close);
    die "No dist close results.$min_dist\n", join(" ", @areas_close) if scalar(@dists_close) == 0;
    
    my @sorted_lifetime_indexes = sort {$lifetimes[$b] <=> $lifetimes[$a]} (@dists_close);
    my $lifetimes_close_count = grep {
        if ($lifetimes[$_] == $lifetimes[$sorted_lifetime_indexes[0]] * $prop_indeter_percent) {
            1;
        } else {
            0;
        }
    } (@dists_close);
    foreach (@sorted_lifetime_indexes) {
        die "Problem with sorted lifetime indexes.\n Values: ", join(" ",@sorted_lifetime_indexes) if $_ > $#merged_areas;
    }

    ###################################
    # Make Decisions  
    ###################################
    if (scalar(@overlap_close) == 1) {
        $merge_decisions[$overlap_close[0]]{winner} = 1;
        $merge_decisions[$overlap_close[0]]{case}   = 1;
    } else {
        if (scalar(@areas_close) == 1) {
            $merge_decisions[$areas_close[0]]{winner} = 1;
            $merge_decisions[$areas_close[0]]{case}   = 2;
        } else {
            if (scalar(@dists_close) == 1) {
                $merge_decisions[$dists_close[0]]{winner} = 1;
                $merge_decisions[$dists_close[0]]{case}   = 3;
             } else {
                $merge_decisions[$sorted_lifetime_indexes[0]]{winner} = 1;
                $merge_decisions[$sorted_lifetime_indexes[0]]{case}   = 4;
                $merge_decisions[$sorted_lifetime_indexes[0]]{case}   = 5 if $lifetimes_close_count > 1;
            }
        }
    } 
    
    #To detect dead adhesions we will examine the pixel similarities counts, if
    #the similarity to the next adhesion is zero, then we know there wasn't any
    #overlap with the merged adhesion and it is unlikely that the event
    #constitutes a merge. Instead, it is more likely that we are observing a
    #death event. Also, make sure there is only one and no more winners
    my $winner_count = 0;
    for (0 .. $#merge_decisions) {
        if ($pix_sims[$_] == 0) {
            $merge_decisions[$_]{dead} = 1;
        }
        $winner_count++ if ($merge_decisions[$_]{winner});
    }
    die "Too many winners found in merge event." if ($winner_count > 1);
    die "No winners found in merge event." if ($winner_count <= 0);

    return @merge_decisions;
}

sub detect_new_adhesions {
    my ($i_num, $next_i_num) = @_;

    my $expected_ad_count = $#{ $data_sets{$next_i_num}{Area} };
    my %expected_ad_nums  = map { $_ => 0 } (0 .. $expected_ad_count);
    my $lineage_length    = $#{ $tracking_mat[0] };
    my @births;

    for my $i (0 .. $#tracking_mat) {
        next if ($tracking_mat[$i][$lineage_length] <= -1);
        $expected_ad_nums{ $tracking_mat[$i][$lineage_length] } = 1;
    }

    for my $i (sort { $a <=> $b } keys %expected_ad_nums) {
        if (not($expected_ad_nums{$i})) {
            $tracking_facts{$i_num}{new_count}++;
            my @temp;
            for (0 .. $lineage_length - 1) {
                push @temp, -1;
            }
            #Check to see if this unassigned adhesion was part of a split birth,
            #if so replace the last entry in the tracking matrix with the
            #adhesion number that it merged from + 1 times -1
            if (defined $tracking_facts{$i_num}{split_birth}[$i] && $tracking_facts{$i_num}{split_birth}[$i] >= 0) {
                $temp[$#temp] = -1*($tracking_facts{$i_num}{split_birth}[$i] + 2);
            }
            push @temp,         $i;
            push @tracking_mat, \@temp;
            push @births, $#tracking_mat;
        }
    }

    # if no FAs were born in this frame, write a -1 to the birth matrix
    if (scalar @births == 0) {
        push @births, -1;
    }
}

sub check_tracking_mat_integrity {
    if (grep $#{$tracking_mat[$_]} != $#{$tracking_mat[0]}, (0 .. $#tracking_mat)) {
        die "Tracking matrix entries not at same length on cycle ", $#{$tracking_mat[0]};
    }

    my @assigned_ad_nums = map {
        if ($tracking_mat[$_][-1] >= 0) {
            $tracking_mat[$_][-1];
        } else {
            ();
        }
    } (0 .. $#tracking_mat);
    @assigned_ad_nums = sort {$a <=> $b} @assigned_ad_nums;
    
    for (0 .. $#assigned_ad_nums) {
        if ($_ != $assigned_ad_nums[$_]) {
            print Dumper(@assigned_ad_nums);
            print "$_ => $assigned_ad_nums[$_]\n";
            die "Not all ad numbers filled in on cycle ",  $#{$tracking_mat[0]} + 1;
        }
    }
}

#######################################
# Output Tracking Facts
#######################################
sub output_tracking_facts {
    print "# of Adhesion Lineages: ", scalar(@tracking_mat), "\n";
    print "# of Live Adhesions Tracked: ",      &get_all_i_num_count("live_adhesions"),            "\n";
    print "# of Best Pixel Sim Not Selected: ", &get_all_i_num_count("best_pix_sim_not_selected"), "\n";
    print "# of Multiple Good P Sims: ",        &get_all_i_num_count("multiple_good_p_sims"),      "\n";
    print "# of Dist and P Sim Guesses Diff: ", &get_all_i_num_count("dist_p_sim_guess_diff"),     "\n";
    print "# of Merge Operations/# Dead Merge Winners: ",
      &get_all_i_num_count("merged_count"), "/", &get_all_i_num_count("dead_winners"), "\n";
    print "# of Merge Operations - Case 1: ",
      &get_all_i_num_count("merge_case_1")/&get_all_i_num_count("merged_count"),"\n" ;
    print "# of Merge Operations - Case 2: ",
      &get_all_i_num_count("merge_case_2")/&get_all_i_num_count("merged_count"),"\n" ;
    print "# of Merge Operations - Case 3: ",
      &get_all_i_num_count("merge_case_3")/&get_all_i_num_count("merged_count"),"\n" ;
    print "# of Merge Operations - Case 4: ",
      &get_all_i_num_count("merge_case_4")/&get_all_i_num_count("merged_count"),"\n" ;
    print "# of Merge Operations - Case 5: ",
      &get_all_i_num_count("merge_case_5")/&get_all_i_num_count("merged_count"),"\n" ;
}

sub get_all_i_num_count {
    my $prop  = $_[0];
    my $total = 0;
    for (keys %tracking_facts) {
        next if not defined $tracking_facts{$_}{$prop};
        $total += $tracking_facts{$_}{$prop};
    }
    return $total;
}

#######################################
# Output Tracking Issues
#######################################
sub output_tracking_probs {
    if (not defined $cfg{tracking_probs_folder}) {
        print "Could not find variable tracking_probs_folder in config, skipping output of tracking problems\n"
          if $cfg{debug};
    } else {
        &output_merge_problems;
    }
}

sub output_merge_problems {
    my $full_probs_folder = catdir($cfg{exp_results_folder}, $cfg{tracking_probs_folder});

    mkpath(catdir($full_probs_folder, "merge"));

    for my $i (keys %{ $tracking_probs{merge} }) {
        for my $j (0 .. $#{ $tracking_probs{merge}{$i} }) {
            mkpath(catdir($full_probs_folder, "merge", $i, $j));
            my %merge_prob = %{ $tracking_probs{merge}{$i}->[$j] };
            for my $k (keys %merge_prob) {
                open OUTPUT, ">" . catfile($full_probs_folder, "merge", $i, $j, "$k.csv");
                print OUTPUT join(",", @{ $merge_prob{$k} });
                close OUTPUT;
            }
        }
    }
}

#######################################
# Misc Helper Functions
#######################################

sub determine_ahesion_lifetime {
    my @tracking_row = @_;

    return scalar( grep $tracking_row[$_] >= 0, (0 .. $#tracking_row) );
}

################################################################################
#Documentation
################################################################################

=head1 NAME

track_adhesions.pl - Track the focal adhesions identified in prior steps 

=head1 SYNOPSIS

track_adhesions.pl -cfg FA_config

=head1 Description

This program handles tracking the FA from birth until death. The process is
split into four parts:

=over 

=item 1. Initialize the tracking sequence on the FA present in the first image

=item 2. Find the best match for each FA in the next image

=item 3. Resolve any conflicts where two FA are predicted to become the same FA,
either by merging or death

=item 4. Add any unassigned FA into the list of current FA and execute steps 2-4
until all the image are cycled through

=back

The results of this algorithm is a matrix with assigns each adhesion number in
each image to part of a sequence, which identifies the identified adhesions in
each image that make up the entire adhesion life cycle. Functionally, the matrix
is output as a csv file that is used by other programs to determine the
properties of the adhesions and visualize them.

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=back

Optional parameter(s):

=over 

=item * input or i: Look for files in each inidividual image directory with this
file name, which should be a Storable object created by an earlier program
execution 

=item * output or o: Create files using the spectified file name in each
individual directory with Storable containing FA information; can be loaded
using the -i flag

=item * debug or d: print debuging information during program execution

=item * emerald: setups and runs a job tailored for the LSF job system on emerald

=back

=head1 EXAMPLES

track_adhesions.pl -cfg FA_config

OR

track_adhesions.pl -cfg FA_config -i data.stor -o data.stor

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 6/6/2008 

=cut
