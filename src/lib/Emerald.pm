#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
{
package Emerald;
use strict;
use warnings;
use File::Spec;
use Data::Dumper;

###############################################################################
# Module Definition
###############################################################################
my %opt = ("queue" => "day", "output_filename" => "out.txt", 
           "error_filename" => "error.txt", "error_folder" => "./",
           "runtime" => "24",);

sub send_LSF_commands {
    my @commands = @{$_[0]};
    
    foreach (@commands) {
        system("$_ > /dev/null 2>/dev/null");
    }
}

sub create_LSF_Matlab_commands {
    my @commands = @{$_[0]};
    if (scalar(@_) > 1) {
        my %temp = %{$_[1]};
        $opt{$_} = $temp{$_} foreach (keys %temp);
    }
    $opt{output_file} = File::Spec->catfile($opt{error_folder},$opt{output_filename});
    $opt{error_file}  = File::Spec->catfile($opt{error_folder},$opt{error_filename});
    unlink($opt{output_file}, $opt{error_file});
    
    my $bsub_command   = "bsub -mig 0 -r -q $opt{queue} -o $opt{output_file} -e $opt{error_file} -We $opt{runtime}";
    if (defined $opt{job_group}) {
        $bsub_command .= " -g $opt{job_group}";
    }
    if (defined $opt{resource}) {
        foreach (split(/\s/, $opt{resource})) {
            $bsub_command .= " -R $_";
        }
    }
    if (defined $opt{memsize}) {
		$bsub_command .= " -memsize $opt{memsize}";
    }
    my $matlab_command = "matlab -singleCompThread -nodisplay -nojvm -nosplash -r";

    @commands = map { split(/\n/, $_) } @commands;
    @commands = map { "$bsub_command $matlab_command \"$_\""} @commands;
    return @commands;   
}

sub create_general_LSF_commands {
    my @commands = @{$_[0]};
    if (scalar(@_) > 1) {
        my %temp = %{$_[1]};
        $opt{$_} = $temp{$_} foreach (keys %temp);
    }
    $opt{output_file} = File::Spec->catfile($opt{error_folder},$opt{output_filename});
    $opt{error_file}  = File::Spec->catfile($opt{error_folder},$opt{error_filename});
    unlink($opt{output_file}, $opt{error_file});
    
    my $bsub_command   = "bsub -mig 0 -r -q $opt{queue} -o $opt{output_file} -e $opt{error_file} -We $opt{runtime}";
    if (defined $opt{job_group}) {
        $bsub_command .= " -g $opt{job_group}";
    }
    if (defined $opt{resource}) {
        foreach (split(/\s/, $opt{resource})) {
            $bsub_command .= " -R $_";
        }
    }
    if (defined $opt{memsize}) {
		$bsub_command .= " -memsize $opt{memsize}";
    }
    @commands = map { split(/\n/, $_) } @commands;
    @commands = map { "$bsub_command \"$_\""} @commands;
    return @commands;   
}

1;
}
