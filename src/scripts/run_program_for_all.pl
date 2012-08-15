#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use lib "../lib";
use lib "../lib/perl";

use File::Basename;
use File::Spec::Functions;
use Cwd 'abs_path';
use Getopt::Long;
use File::Find::Rule;

use Config::Adhesions qw(ParseConfig);

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
$opt{extra} = "";
GetOptions(\%opt, "cfg|config=s", "debug|d", "program|p=s", "extra|e=s", 
                  "run_all_debug", "exp_filter=s") or die;

die "Can't find cfg file specified on the command line" if not(exists $opt{cfg});
die "Can't find program to execute on the command line" if not(exists $opt{program});

print "Collecting Configuration\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

###############################################################################
#Main Program
###############################################################################
chdir(dirname($opt{program}));

my $program_base = basename($opt{program});

my $debug_string = ($opt{debug}) ? "-d" : "";

my $cfg_suffix = basename($opt{cfg});
$cfg_suffix =~ s/.*\.(.*)/$1/;

my @config_files = File::Find::Rule->file()->name( "*.$cfg_suffix" )->in( ($cfg{data_folder}) );
@config_files = sort @config_files;
if (exists($opt{exp_filter})) {
   @config_files = grep $_ =~ /$opt{exp_filter}/, @config_files;
}

my $parallel_return = system("which parallel");

my @command_set;
foreach (@config_files) {
    next if /config\/default/;

    my $command = "./$program_base -cfg $_ $debug_string $opt{extra}";
    $command =~ s/"/\\"/g;
    push @command_set, "\"$command; echo $_;\"";
	if ($parallel_return != 0) {
		if ($opt{debug}) {
			print "$command\n";
		} else {
			system("ionice -c3 nice -n 20 $command");
		}
	}
}

if ($parallel_return == 0) {
	if ($opt{debug}) {
		my $parallel_cmd = "time parallel --ungroup --nice 20 ::: \n\t" . join("\n\t", @command_set) . "\t\n";
		print $parallel_cmd;
	} else {
		my $parallel_cmd = "time parallel --ungroup --nice 20 ::: " . join(" ", @command_set);
		system($parallel_cmd);
	}
}

###############################################################################
#Documentation
###############################################################################

=head1 NAME

run_perl_program_for_all.pl

=head1 SYNOPSIS

run_program_for_all.pl -cfg config.cfg -p program 

=head1 Description

This program provides the means to run a program for each of the experiments in a set. This is helpful when a single step of a processing pipeline needs to be executed again. Config files (those with the same suffix as the config file provided on the command line) are searched for in subdirectories of the $cfg{data_folder} directory.

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=item * program or p: the path to a program that should be run using the standard "-cfg exp_config.cfg" parameter

=back

Optional parameter(s):

=over 

=item * debug or d: if present do not run the program specified in the program parameter, instead print the commands that would be run

=back

=head1 EXAMPLES

run_program_for_all.pl -cfg FA_config -p ../analyze_cell_features/gather_tracking_results.pl

OR

run_program_for_all.pl -cfg FA_config -p ../analyze_cell_features/gather_tracking_results.pl -e '-skip_pix_props'

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 6/9/2008 

=cut
