#!/usr/bin/perl -wT

use strict;
use File::Path;
use File::Basename;
use File::Spec::Functions;
use File::Copy;
use File::Temp qw(tempfile);
use POSIX;
use CGI;
use CGI::Carp;
use IO::Handle;
use Config::General;

my $start = time;

$| = 1;
###############################################################################
# Configuration
###############################################################################

my $upload_dir = 'upload';
my $run_exp_dir = '/home/mbergins/Documents/Projects/fa_webapp_run/src/webapp/';

###############################################################################
# Main Program
###############################################################################

my $q = CGI->new();

print $q->header,
	  $q->start_html(-title=>'Focal Adhesion Analysis Server',-style=>'/FA_webapp/css/screen.css');
	  
print "<div class=\"container\">\n";
print $q->h1('Focal Adhesion Analysis Server - Server Status');

my @upload_zips = <$upload_dir/*.zip>;
print $q->p, "<b>Experiments in queue:</b> ", scalar(@upload_zips);

my @run_files = <$run_exp_dir/fa_webapp.*.run>;
print $q->p, "<b>Experiments being processed:</b> ", scalar(@run_files);

my %uptime_props = &process_uptime_reading;

print $q->p, "<b>The server has been running for:</b> ", $uptime_props{runtime};
print $q->p, "<b>The server load levels are:</b> ", $uptime_props{load};

print $q->br,$q->br,$q->hr, $q->p, 
	"Return to the ", $q->a({href=>"/FA_webapp"},'home page');

print "</div>\n";
print $q->end_html;

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
