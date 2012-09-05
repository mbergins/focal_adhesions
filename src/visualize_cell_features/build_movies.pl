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

#Build Movies 
my @image_folders = <$cfg{individual_results_folder}/*>;
my $image_num_length = length(scalar(@image_folders));

print "\n\nBuild Movies\n\n" if $opt{debug};
$t1 = new Benchmark;

our @movie_dirs; 
find(\&add_to_movie_dir, (catdir($cfg{exp_results_folder},$cfg{movie_output_folder})));

foreach my $f1 (@movie_dirs) {
    foreach my $f2 (@{ $cfg{movie_output_prefix} }) {
		my $unfoldered_exp_name = $cfg{exp_name};
		$unfoldered_exp_name =~ s#/#_#g;
        my $input_folder = catdir($f1,$f2);
		my $output_file = catdir($cfg{exp_results_folder},'..',"$unfoldered_exp_name.mov");
		
		my $command = "avconv -v 0 -y -r $cfg{movie_frame_rate} -i $input_folder/%0" . 
			$image_num_length . "d.png -qscale 5 $output_file > /dev/null 2>&1";
		if ($opt{debug}) {
        	print "$command\n";
		} else {
        	system($command);
		}
    }
}
$t2 = new Benchmark;
print "Runtime: ",timestr(timediff($t2,$t1)), "\n" if $opt{debug};

###############################################################################
# Functions
###############################################################################
sub add_to_movie_dir {
    if ($File::Find::name =~ /$cfg{vis_config_file}/) {
        push @movie_dirs, $File::Find::dir;
    }
}
