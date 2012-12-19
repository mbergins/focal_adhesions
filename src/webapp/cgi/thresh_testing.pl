#!/usr/bin/perl -w

use strict;
use File::Path;
use File::Basename;
use File::Spec::Functions;
use File::Copy;
use File::Temp qw(tempfile tempdir);
use POSIX;
use CGI::Pretty;
use CGI::Carp qw(fatalsToBrowser);
use IO::Handle;
use Config::General;
use Data::Dumper;
use HTML::Template;

my $upload_start = time;

$| = 1;

###############################################################################
# Configuration
###############################################################################

my $output_dir = "$ENV{DOCUMENT_ROOT}/FA_webapp/threshold_testing";
$output_dir =~ /(.*)/;
$output_dir = $1;
if (! -e $output_dir) {
	mkpath($1);
}

###############################################################################
# Main Program
###############################################################################

my $q = CGI->new();

my $lightweight_fh = $q->upload('uploaded_file');
if (defined $lightweight_fh) {
	my $start_time = time;
    my $io_handle = $lightweight_fh->handle();
    $q->param('uploaded_file') =~ /(.*)/;
    
	my $dir = tempdir( DIR => $output_dir );
	chmod 0777, $dir;
	my $output_file = catfile($dir,'focal_image');
	my $fh = FileHandle->new(">$output_file");

    my $buffer;
    while (my $bytesread = $io_handle->read($buffer,1024)) {
        print $fh $buffer or die;
    }
    close $fh;
    print "Done loading file.";
    chmod 0777, "$output_file" or die "$!";
	my $upload_end = time;

	my $start_process = time;
	system("octave --eval \"build_thresholded_image_sets('$output_file')\" > /dev/null 2> /dev/null");
	my $end_process = time;
	
	my $template = HTML::Template->new(filename => 'template/threshold_results.tmpl');
	$template->param(run_time => $end_process - $start_process,upload_time => $upload_end - $upload_start);
	
	open OUTPUT, ">$dir/index.html" or die $!;
	print OUTPUT $template->output;
	close OUTPUT;
	
	my $url_out_dir = $dir;
	$url_out_dir =~ s/$ENV{'DOCUMENT_ROOT'}//g;
	my $hostname = $ENV{'HTTP_ORIGIN'};
	my $final_url = "$ENV{'HTTP_ORIGIN'}/$url_out_dir";

	print $q->redirect($final_url);
} else {
    #no valid file handle was available, instead present the user with the
    #regular upload page
	my $template = HTML::Template->new(filename => 'template/threshold_start.tmpl');
	
	print $q->header();
	print $template->output;
}


###############################################################################
# Functions
###############################################################################
