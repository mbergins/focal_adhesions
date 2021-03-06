#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Basename;
use File::Spec::Functions;
use Getopt::Long;

use lib "../lib";
use Config::Adhesions;
use Image::Data::Collection;

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
my $vis_dir = catdir($cfg{exp_results_folder},$cfg{movie_output_folder});

my @filtered_movie_folders = qw(axis_ratio leading_vs_trailing live_at_15);

foreach (@filtered_movie_folders) {
	my @match_folders = <$vis_dir/$_*>;
	
	foreach my $this_vis_folder (@match_folders) {
		my $match_folder_name = basename($this_vis_folder);

		my $unfoldered_exp_name = $cfg{exp_name};
		$unfoldered_exp_name =~ s#/#_#g;

		my $output_file = catfile($cfg{exp_results_folder}, "..",
			$unfoldered_exp_name . "_" . $match_folder_name . ".mov");
		my $command = "ffmpeg -v 0 -y -r 10 -i  " . 
			"$this_vis_folder/%04d.png -sameq $output_file > /dev/null 2>&1";

		if ($opt{debug}) {
        	print "$command\n";
		} else {
        	system($command);
		}
	}
}
