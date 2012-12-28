#!/usr/bin/perl -w

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

my $uptime_file = "/home/mbergins/Documents/uptime_readings.txt";
###############################################################################
# Main Program
###############################################################################

&build_server_load_day_plot;
&build_server_load_week_plot;

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

sub build_server_load_day_plot {
	system("tail -n 1440 $uptime_file > last_day_load.txt;");
	system("sed -i \"s/.*average: //\" last_day_load.txt");

	system("R --vanilla --silent -f plot_webapp_load.R \"--args data_file=last_day_load.txt plot_type=day target_dir=/var/www/FA_webapp/images/\" > /dev/null");

	unlink("last_day_load.txt");
}

sub build_server_load_week_plot {
	system("tail -n 10080 $uptime_file > last_week_load.txt;");
	system("sed -i \"s/.*average: //\" last_week_load.txt");

	system("R --vanilla --silent -f plot_webapp_load.R \"--args data_file=last_week_load.txt plot_type=week target_dir=/var/www/FA_webapp/images/\" > /dev/null");

	unlink("last_week_load.txt");
}
