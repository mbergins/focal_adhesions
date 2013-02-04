package server_status;
use Dancer ':syntax';
use strict;
use warnings;
use Cwd;
use Sys::Hostname;
use File::Temp qw/ tempfile /;
use File::Basename qw/ dirname /;
use File::Spec::Functions;
use File::Path qw/ make_path /;
use File::Copy qw/ move /;
use File::Find;
use File::Basename;
use Config::General;

my $upload_dir = catdir('..','uploaded_experiments');

my $run_exp_dir = '/home/mbergins/Documents/Projects/focal_adhesions/trunk/src/webapp_dancer/utilities/';

my $uptime_file = "/home/mbergins/Documents/uptime_readings.txt";
my $cron_file = "/home/mbergins/Documents/current_cron";
###############################################################################
# Main
###############################################################################

get '/server_status' => sub {
	&build_server_load_day_plot;
	&build_server_load_week_plot;
	
	my %template;

	my @upload_zips = <$upload_dir/*>;
	$template{queue_count} = scalar(@upload_zips);

	my @run_files = <$run_exp_dir/fa_webapp.*.run>;
	$template{run_count} = scalar(@run_files);
	$template{worker_count} = &count_upload_workers($cron_file);
	
	my %uptime_props = &process_uptime_reading;

	$template{server_uptime} = $uptime_props{runtime};
	$template{server_load} = $uptime_props{load};

	template 'server_status', \%template;
};

###############################################################################
# Functions
###############################################################################

sub process_uptime_reading {
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
	open INPUT, $_[0] or die $!;
	my @cron = <INPUT>;
	close INPUT;
	
	@cron = grep !($_ =~ /^#/), @cron;
	@cron = grep $_ =~ /run_uploaded_exp/, @cron;

	return(scalar(@cron));
}

sub build_server_load_day_plot {
	system("tail -n 1440 $uptime_file > last_day_load.txt;");
	system("sed -i \"s/.*average: //\" last_day_load.txt");
	
	my $image_dir = "/home/mbergins/Documents/Projects/focal_adhesions/trunk/src/webapp_dancer/public/images";

	system("R --vanilla --silent -f plot_webapp_load.R \"--args data_file=last_day_load.txt plot_type=day target_dir=$image_dir\" > /dev/null");

	unlink("last_day_load.txt");
}

sub build_server_load_week_plot {
	system("tail -n 10080 $uptime_file > last_week_load.txt;");
	system("sed -i \"s/.*average: //\" last_week_load.txt");

	my $image_dir = "/home/mbergins/Documents/Projects/focal_adhesions/trunk/src/webapp_dancer/public/images";
	
	system("R --vanilla --silent -f plot_webapp_load.R \"--args data_file=last_week_load.txt plot_type=week target_dir=$image_dir\" > /dev/null");

	unlink("last_week_load.txt");
}

true;
