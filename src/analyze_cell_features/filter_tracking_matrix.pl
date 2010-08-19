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
use Text::CSV::Simple;

use Config::Adhesions qw(ParseConfig);
use Image::Data::Collection;
use Text::CSV::Simple::Extra;
use Emerald;
use FA_job;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg=s", "debug|d", "lsf|l") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration\n" if $opt{debug};

my %cfg = ParseConfig(\%opt);

###############################################################################
# Main Program
###############################################################################
if ($opt{lsf}) {
    my @command = "$0 -cfg $opt{cfg}";
    $opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'track_filter');
    &FA_job::send_general_lsf_program(\@command,\%opt);
    
    exit;
}

print "\n\nCollecting Tracking Matrix\n" if $opt{debug};
my @tracking_mat = &read_in_tracking_mat(\%cfg, \%opt);

print "\n\nGathering Adhesion Lineage Properties\n", if $opt{debug};
my %lin_props = &collect_ad_props;

print "\n\nFiltering Tracking Matrix\n" if $opt{debug};
my %filtered_matrix_set = &filter_tracking_matrix;

print "\n\nOutputing Tracking Matrix\n" if $opt{debug};
&output_filtered_matrices(\%filtered_matrix_set, "");

###############################################################################
# Functions
###############################################################################

sub collect_ad_props {
    my $data_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_summary_props_file});

    my $parser = Text::CSV::Simple->new();
    my @data   = $parser->read_file($data_file);

    my @header_line = @{ shift @data };

    my %lin_props;
    foreach (@data) {
        my @data_line = @{$_};
        foreach my $i (0 .. $#header_line) {
            push @{ $lin_props{ $header_line[$i] } }, $data_line[$i];
        }
    }
    return %lin_props;
}

sub filter_tracking_matrix {
    my %matrix_set;

    for my $required_longevity (split(/\s+/, $cfg{required_longevity_filter})) {
        for my $i (0 .. $#{ $lin_props{'longevity'} }) {
            my $this_longev = $lin_props{'longevity'}[$i];
            next if ($this_longev eq "NA");
            if ($this_longev >= $required_longevity) {
                push @{ $matrix_set{'longevity'}{$required_longevity} }, $tracking_mat[$i];
            }
            if ($this_longev >= $required_longevity && $lin_props{'death_status'}[$i]) {
                push @{ $matrix_set{'dead'}{$required_longevity} }, $tracking_mat[$i];
            }
            if ($this_longev >= $required_longevity && $lin_props{'split_birth_status'}[$i]) {
                push @{ $matrix_set{'split_birth'}{$required_longevity} }, $tracking_mat[$i];
            }
        }
        @{ $matrix_set{'split_birth'}{$required_longevity} } =
          &add_split_birth_parents(@{ $matrix_set{'split_birth'}{$required_longevity} });
    }

    #ad-hoc line to pick out specific lineages
    #@{$matrix_set{'special'}{'high_speed'}} = map $tracking_mat[$_], (146,262,516);

    #Filter the tracking matrix if a 'for_vis' folder is available
    my $R_sq_folder = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, 'for_vis');
    if (-d $R_sq_folder) {
        for my $file_name (<$R_sq_folder/*.csv>) {
            $file_name =~ /$R_sq_folder\/(.*).csv/;

            my $parser = Text::CSV::Simple->new;
            eval {
                my @data = $parser->read_file($file_name);
                @{ $matrix_set{$1} } = map $tracking_mat[ $_->[0] - 1 ], (@data);
            };
        }
    }

    return %matrix_set;
}

sub add_split_birth_parents {
    my @mat_set = @_;

    my @set_with_parents;
    for my $i (0 .. $#mat_set) {
        my $pre_birth_index =
          (grep $mat_set[$i][$_] >= 0 && $mat_set[$i][ $_ - 1 ] <= -2, (1 .. $#{ $mat_set[$i] }))[0] - 1;
        die "Error identifying matrix index of split birth event" if ($mat_set[$i][$pre_birth_index] >= -1);

        my $ad_parent_num = -1 * ($mat_set[$i][$pre_birth_index] + 2);

        my @parent_lin_num = grep $tracking_mat[$_][ $pre_birth_index + 1 ] == $ad_parent_num, (0 .. $#tracking_mat);
        die "Expected to only find one parent lineage" if (scalar(@parent_lin_num) > 1);
        die "Unable to find any parent lineages\n", join(" ", @{ $mat_set[$i] }), "\n"
          if (scalar(@parent_lin_num) == 0);

        push @set_with_parents, $tracking_mat[ $parent_lin_num[0] ];
        push @set_with_parents, $mat_set[$i];
    }
    die if 2 * scalar(@mat_set) != scalar(@set_with_parents);

    return @set_with_parents;
}

sub output_filtered_matrices {
    my %mat_set = %{ $_[0] };
    my $prefix  = $_[1];

    my $base_folder = catdir($cfg{exp_results_folder}, $cfg{tracking_folder}, 'filtered');
    mkpath(catdir($base_folder, $prefix));

    foreach my $i (keys %mat_set) {
        if (ref($mat_set{$i}) eq "ARRAY") {
            if (scalar(@{ $mat_set{$i} }) > 0) {
                output_mat_csv(\@{ $mat_set{$i} }, catfile($base_folder, $prefix, $i . '.csv'));
            }
        } elsif (ref($mat_set{$i}) eq "HASH") {
            output_filtered_matrices(\%{ $mat_set{$i} }, catdir($prefix, $i));
        } else {
            die "unexpected data type $prefix";
        }
    }
}

################################################################################
#Documentation
################################################################################

=head1 NAME

filter_tracking_matrix.pl - Filter and output the tracking matrix based on collected FA properties 

=head1 SYNOPSIS

filter_tracking_matrix.pl -cfg FA_config

=head1 Description

After collecting the properties of the FA, this program uses those collected properties to filter the tracking matrix to only those adhesions which meet certain criteria. The criteria can be any property collected in earlier stages.

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=back

Optional parameter(s):

=over 

=item * debug or d: print debuging information during program execution

=item * emerald: setups and runs a job tailored for the LSF job system on emerald

=back

=head1 EXAMPLES

filter_tracking_matrix.pl -cfg FA_config

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 6/6/2008 

=cut
