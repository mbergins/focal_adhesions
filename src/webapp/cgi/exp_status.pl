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

use lib './';
use webserver_funcs qw(print_html_end);

my $start = time;

$| = 1;
###############################################################################
# Configuration
###############################################################################

my $upload_dir = 'upload';
my $results_dir = '/var/www/FA_webapp/results/';

###############################################################################
# Main Program
###############################################################################

my $q = CGI->new();

print $q->header,
	  $q->start_html(-title=>'Focal Adhesion Analysis Server - Experiment Status',
		  -style=>'/FA_webapp/css/screen.css');

print "<div class=\"container\">\n";
print $q->h1('Focal Adhesion Analysis Server - Experiment Status');

if (not defined $q->param('exp_id')) {
	print $q->p, "I'm sorry, no experiment ID was specified. If you followed a
	link please make sure it appears correct in the URL bar.";
	&print_html_end($q);
	exit;
}

my $exp_id = $q->param('exp_id');

print $q->p, "<b>Experiment ID:</b> ", $exp_id;

my @upload_zips = <$upload_dir/*.zip>;
my @results_files = <$results_dir/*>;
if (not grep $_ =~ /$exp_id/, @upload_zips) {
	if (not grep $_ =~ /$exp_id/, @results_files) {
		print $q->p, "I'm sorry, but the experiment ID you specified ($exp_id)
		wasn't found. If you followed a link please make sure it appears correct
		in the URL bar.";
		&print_html_end($q);
		exit;
	}
}

my $queue_pos = &find_exp_position($exp_id,\@upload_zips);
if (defined $queue_pos) {
	print $q->p, "<b>Experiments in queue:</b> ", scalar(@upload_zips);
	print $q->p, "<b>Position in queue:</b> ", $queue_pos;
	&print_html_end($q);
	exit;
}

my @results_zips = <$results_dir/*.zip>;
if (grep $_ =~ /$exp_id/, @results_zips) {
	print $q->p, "Your experiment has finished processing. You can download your
	results ", $q->a({href=>"/FA_webapp/results/$exp_id.zip"}, "here") , ".";
} else {
	print $q->p, "Your experiment is currently processing.";
}
&print_html_end($q);

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
