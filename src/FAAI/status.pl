#!/usr/bin/perl -wT

###############################################################################
# Setup
###############################################################################

use strict;
use File::Path;
use File::Basename;
use File::Spec::Functions;
use File::Copy;
use File::Temp;
use POSIX;
use CGI;
use CGI::Carp;
use IO::Handle;
use Config::General;
use Cwd;

$| = 1;

my $upload_dir = catdir(getcwd,'..','..','upload');
$upload_dir =~ /(.*)/;
$upload_dir = $1;
if (not -e $upload_dir) {
    mkdir($upload_dir) or die $!;
}

my $results_dir = catdir(getcwd,'..','..','results');

###############################################################################
# Main
###############################################################################

my $q = CGI->new();

print $q->header,
	  $q->start_html(-title=>'Focal Adhesion Alignment Index Server - Status Page'), 
	  $q->h1('Focal Adhesion Alignment Index Server - Status Page');

my $vars = $q->Vars;
my $exp_id = $vars->{exp_id};

#Deal with the case that the experiment ID was not specified in the URL bar
if (not defined $exp_id) {
    print "The experimental id was not specified, please go back and check the link provided through email.";
    print $q->end_html;
    exit();
}

#deal with the case where an exp id was provided, but isn't present in the
#results directory
if (! -e catdir($results_dir,$exp_id)) {
    my @upload_files = <$upload_dir/*>;
    
    my %files_ages;
    for (@upload_files) {
        my $file_name = basename($_);
        if ($file_name =~ /.cfg/) {
            $files_ages{$file_name} = (stat($_))[9];
        }
    }

    my @cfg_sort = sort { $files_ages{$a} <=> $files_ages{$b} } keys %files_ages;
    
    #if the cfg file is not present this list will be empty, so first we check
    #for that and then exit with a notice to the user
    my @cfg_position = grep {$cfg_sort[$_] =~ /$exp_id/} (0..$#cfg_sort);
    if (not defined $cfg_position[0]) {
        print "Your experiment ($exp_id) was not found in the system, please check the email sent to you after submitting your job.";
        print $q->end_html;
        exit();
    }
    
    my $cfg_position = $cfg_position[0]+1;
    print "Your experiment ($exp_id) is in the queue in position $cfg_position/" . scalar(@cfg_sort) . ", you will be notified via email when it begins and when it ends.", $q->p();
    print $q->end_html;
    exit();
}

#now we know that the experiment has been found, the next step is to determine
#if is is done running, finished experiments will have a zip file present in the
#results directory with the same config name
if (-e catdir($results_dir,"$exp_id.zip")) {
    print "Your experiment ($exp_id) is done, you can download your results from this link:", $q->p();
    print $q->a({href=>"http://balder.bme.unc.edu/~mbergins/FAAI_results/$exp_id.zip"},
        "Results for exp $exp_id");
} else {
    #no zip file was found, provide a message to the user that their experiment
    #is running
    print "Your experiment ($exp_id) is currently running, you will receive an email when it is done processing.", $q->p();
}

print $q->end_html;
