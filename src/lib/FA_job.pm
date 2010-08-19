#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
{
package FA_job;
use strict;
use warnings;
use File::Spec;
use File::Path;
use Data::Dumper;
use Emerald;
use Math::Matlab::Extra;

###############################################################################
# Module Definition
###############################################################################
sub run_matlab_progam {
    my @matlab_code = @{$_[0]};
    my %opt = %{$_[1]};

    if (not defined $opt{error_folder}) {
        $opt{error_folder} = '.';
    }

    mkpath($opt{error_folder});
    
    if ($opt{lsf}) {
        my @commands = &Emerald::create_LSF_Matlab_commands(\@matlab_code, \%opt);
        if ($opt{debug}) {
            print "\n", join("\n", @commands), "\n";
        } else {
            &Emerald::send_LSF_commands(\@commands);
        }
    } else {
        if ($opt{debug}) {
            print "\n", join("\n", @matlab_code), "\n";
        } else {
            &Math::Matlab::Extra::execute_commands(\@matlab_code, $opt{error_file});
        }
    }
}

sub send_general_lsf_program {
    my @commands = @{$_[0]};
    my %opt = %{$_[1]};

    if (not defined $opt{error_folder}) {
        $opt{error_folder} = '.';
    }

    mkpath($opt{error_folder});
    
    @commands = &Emerald::create_general_LSF_commands(\@commands,\%opt);
    if ($opt{debug}) {
        print join("\n", @commands);
    } else {
        &Emerald::send_LSF_commands(\@commands);
    }
}

1;
}
