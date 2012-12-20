#!/usr/bin/perl -w

use File::Basename;
use Getopt::Long;

my %opt;
GetOptions(\%opt, "debug|d", "type=s") or die;

###############################################################################
# Main
###############################################################################

if ($opt{type} eq "day") {
	system("tail -n 1440 ~/Documents/uptime_readings.txt > ~/Documents/Projects/focal_adhesions/trunk/src/webapp/last_day_load.txt;");
	system("sed -i \"s/.*average: //\" last_day_load.txt");

	system("R --vanilla --silent -f plot_webapp_load.R \"--args data_file=last_day_load.txt plot_type=day target_dir=/var/www/FA_webapp/images/\" > /dev/null");

	unlink("last_day_load.txt");
}

if ($opt{type} eq "week") {
	system("tail -n 10080 ~/Documents/uptime_readings.txt > ~/Documents/Projects/focal_adhesions/trunk/src/webapp/last_week_load.txt;");
	system("sed -i \"s/.*average: //\" last_week_load.txt");

	system("R --vanilla --silent -f plot_webapp_load.R \"--args data_file=last_week_load.txt plot_type=week target_dir=/var/www/FA_webapp/images/\" > /dev/null");

	unlink("last_week_load.txt");
}

###############################################################################
# Functions
###############################################################################
