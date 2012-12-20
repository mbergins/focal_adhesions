#!/usr/bin/perl -wT

use strict;
use File::Path;
use File::Basename;
use File::Spec::Functions;
use File::Copy;
use File::Temp qw(tempfile);
use POSIX;
use CGI::Pretty;
use CGI::Carp;
use IO::Handle;
use Config::General;
use HTML::Template;

my $start = time;

$| = 1;
###############################################################################
# Configuration
###############################################################################

my $upload_dir = 'upload';
my $running_results_dir = '/home/mbergins/Documents/Projects/focal_adhesions/trunk/results/';
my $final_results_dir = '/var/www/FA_webapp/results/';

###############################################################################
# Main Program
###############################################################################

my $q = CGI->new();

print $q->header();

my $exp_id = $q->param('exp_id');

my $template = HTML::Template->new(filename => 'template/exp_status.tmpl');

if (not defined $q->param('exp_id')) {
	$template->param('no_id' => 1);
	print $template->output;
	exit;
}

#There are three places an experiment might be: in the queue, being processed or
#in the final results directory. I'll check them in that order.
my @final_results_files = <$final_results_dir/*.zip>;
if (grep $_ =~ /$exp_id\.zip/, @final_results_files) {
	$template->param(download_link => "/results/$exp_id.zip",
					 exp_name => $exp_id,
				 	 exp_finished => 1);
	print $template->output;
} else {
	my @running_results_files = <$running_results_dir/*>;
	if (grep $_ =~ /$exp_id$/, @running_results_files) {
		$template->param(exp_name => $exp_id,exp_running => 1);
		print $template->output;
	} else {
		my @upload_zips = <$upload_dir/*.zip>;
		if (grep $_ =~ /$exp_id\.zip/, @upload_zips) {
			my $queue_pos = &find_exp_position($exp_id,\@upload_zips);
			
			$template->param(exp_name => $exp_id, 
							 queue_count => scalar(@upload_zips), 
							 queue_position => $queue_pos,
						 	 exp_in_queue => 1);
			print $template->output;
		} else {
			$template->param(exp_name => $exp_id,exp_wrong_id => 1);
			print $template->output;
		}
	}
}

###############################################################################
# Functions
###############################################################################

sub find_exp_position {
	my $exp_id = $_[0];
	my @upload_zips = @{$_[1]};
	
	my %upload_data;
	for my $file (@upload_zips) {
		my $name = basename($file);
		$name =~ s/\.zip//;
		$upload_data{$name}{age} = -C $file;
		$upload_data{$name}{source} = $file;
	}
	my @sorted_zips = sort {$upload_data{$b}{age} <=> $upload_data{$a}{age}} keys %upload_data;

	# for (@sorted_zips) {
	# 	print "<BR>$_ - $upload_data{$_}{age}\n";
	# }

	my $position;
	for (1..scalar(@sorted_zips)) {
		if ($sorted_zips[$_ - 1] =~ /$exp_id/) {
			$position = $_;
		}
	}
	return $position;
}
