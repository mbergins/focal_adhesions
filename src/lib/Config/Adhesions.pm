#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
{
package Config::Adhesions;
use strict;
use warnings;
use File::Spec;
use File::Basename;

use base qw(Config::General);

our @EXPORT = qw(ParseConfig);
use Exporter;
our @ISA = qw(Exporter);

my %derived_vars = (
    individual_results_folder => [qw(results_folder exp_name single_image_folder)],
    exp_results_folder        => [qw(results_folder exp_name)],
    exp_data_folder           => [qw(data_folder exp_name)],
);

my @vars_to_split = qw(exclude_image_nums general_data_files tracking_files
    lineage_analysis_data_files movie_output_prefix);

###############################################################################
# Module Definition
###############################################################################

sub new {
    my $class        = $_[0];
    my %opt          = %{ $_[1] };

    my %cfg = Config::General::ParseConfig(
        -ConfigFile            => $opt{cfg},
        -MergeDuplicateOptions => 1,
        -IncludeRelative       => 1,
    );
    $cfg{opt} = \%opt;

    my $cfg_ref = \%cfg;

    bless $cfg_ref, $class;
	
    $cfg_ref->build_derived_parameters;
    $cfg_ref->split_config_variables;

    return $cfg_ref;
}

sub ParseConfig {
    my %opt = %{ $_[0] };
    
    my %cfg = Config::General::ParseConfig(
        -ConfigFile            => $opt{cfg},
        -MergeDuplicateOptions => 1,
        -IncludeRelative       => 1,
    );
    $cfg{opt} = \%opt;
    
    my $cfg_ref = \%cfg;

    bless $cfg_ref, "Config::Adhesions";
	
    $cfg_ref->build_derived_parameters;
    $cfg_ref->split_config_variables;

    return $cfg_ref->get_cfg_hash;
}

sub split_config_variables {
    my $cfg = shift;
    
    foreach (@vars_to_split) {
        if (defined $cfg->{$_}) {
            $cfg->{$_} = [split(/\s+/, $cfg->{$_})];
        }
    }

    if (not defined $cfg->{exclude_image_nums}) {
        $cfg->{exclude_image_nums} = [];
    }
}

sub build_derived_parameters {
    my $cfg = $_[0];
	
	my %hash_cfg = %{$cfg};
	
	if (not defined $cfg->{exp_name}) {
		if ($cfg->{opt}{cfg} =~ /$cfg->{data_folder}(.*)/) {
			$cfg->{exp_name} = File::Basename::dirname($1);
		} else {
			die "Unable to find $cfg->{data_folder} in $cfg->{opt}{cfg}.";
		}
	}

    foreach my $this_key (keys %derived_vars) {
        my $all_present = 1;
        foreach my $var_name (@{ $derived_vars{$this_key} }) {
            $all_present = 0 if not defined $cfg->{$var_name};
        }

        if ($all_present) {
            $cfg->{$this_key} = File::Spec->catdir(map $cfg->{$_}, @{ $derived_vars{$this_key} });
        }
    }
}

sub get_cfg_hash {
    my $self = shift;
    return map { $_ => ${$self}{$_} } keys %{$self};
}

}
1;
