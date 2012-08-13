#!/usr/bin/perl -wT

use strict;
use CGI;
use CGI::Carp;
use HTML::Template;

my $start = time;

$| = 1;
###############################################################################
# Configuration
###############################################################################

my $upload_dir = 'upload';
my $run_exp_dir = '/home/mbergins/Documents/Projects/focal_adhesions/trunk/src/webapp/';

###############################################################################
# Main Program
###############################################################################

my $q = CGI->new();

my $template = HTML::Template->new(filename => 'template/server_status.tmpl');

my @upload_zips = <$upload_dir/*.zip>;
$template->param(queue_count => scalar(@upload_zips));

my @run_files = <$run_exp_dir/fa_webapp.*.run>;
$template->param(run_count => scalar(@run_files),
				 worker_count => &count_upload_workers);

my %uptime_props = &process_uptime_reading;

$template->param(server_uptime => $uptime_props{runtime},
				 server_load => $uptime_props{load});

print $q->header();
print $template->output;

###############################################################################
# Functions
###############################################################################

sub process_uptime_reading {
	$ENV{PATH} = '/usr/bin/';
	my $uptime = `/usr/bin/uptime`;
	
	my %uptime_props;
		
	if ($uptime =~ /up (.*?),/) {
		$uptime_props{runtime} = $1;
	}

	if ($uptime =~ /load average: (.*)/) {
		$uptime_props{load} = $1;	
	}

	return %uptime_props;
}

sub count_upload_workers {
	open INPUT, "current_cron";
	my @cron = <INPUT>;
	close INPUT;

	@cron = grep !($_ =~ /^#/), @cron;
	@cron = grep $_ =~ /run_uploaded_exp/, @cron;

	return(scalar(@cron));
}
