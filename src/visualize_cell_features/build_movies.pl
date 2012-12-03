#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Path;
use File::Find;
use File::Spec::Functions;
use Getopt::Long;
use IO::File;
use Benchmark;

use lib "../lib";
use Config::Adhesions;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|config=s", "debug|d") or die;
die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Overall Configuration\n\n" if $opt{debug};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg = $ad_conf->get_cfg_hash;

###############################################################################
# Main Program
###############################################################################
my ($t1, $t2);

print "\n\nBuild Movies\n\n" if $opt{debug};
$t1 = new Benchmark;

my $tracking_folder = catdir($cfg{exp_results_folder},'visualizations','tracking');

my $unfoldered_exp_name = $cfg{exp_name};
$unfoldered_exp_name =~ s#/#_#g;
my $output_file = catdir($cfg{exp_results_folder},'..',"$unfoldered_exp_name.mov");

my $command = "avconv -v 0 -y -r 15 -i $tracking_folder/%05d.png -qscale 0 $output_file > /dev/null 2>&1";
if ($opt{debug}) {
	print "$command\n";
} else {
	system($command);
}

$t2 = new Benchmark;
print "Runtime: ",timestr(timediff($t2,$t1)), "\n" if $opt{debug};

###############################################################################
# Functions
###############################################################################
