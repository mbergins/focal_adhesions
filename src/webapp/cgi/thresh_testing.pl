#!/usr/bin/perl -wT

use strict;
use File::Path;
use File::Basename;
use File::Spec::Functions;
use File::Copy;
use File::Temp qw(tempfile tempdir);
use POSIX;
use CGI;
use CGI::Carp;
use IO::Handle;
use Config::General;
use Data::Dumper;

use lib './';
use webserver_funcs qw(print_html_end);

my $start = time;

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
	
	#place a file in the upload dir to let the worker script know that the
	#upload is complete and then wait for the one of the last images to be
	#written to disk
	open OUTPUT, ">$dir/upload_done";
	print OUTPUT "done";
	close OUTPUT;
	my $sleep_count = 0;
	while (! -e "$dir/thumbs/focal_image_4.jpg") {
		sleep 1;
		$sleep_count++;
		if ($sleep_count > (10*60)) {
			die "Threshold processing on $dir took over ten minutes, exiting.";
		}
	}

	#start a new CGI object to build the browsing page for the thresholded image
	#set
	my $output_cgi = CGI->new();
	open OUTPUT, ">$dir/index.html";
	select OUTPUT;

	print OUTPUT $output_cgi->start_html(
		-title=>'Focal Adhesion Analysis Server - Threshold Testing',
		-style=>'/FA_webapp/css/screen.css');
	print "<div class=\"container\">\n";
	print $output_cgi->h1('Focal Adhesion Analysis Server - Threshold Testing');
	print $output_cgi->p, "It took ", time - $start_time, 
		" seconds to process your image.", $output_cgi->br;
	
	print $output_cgi->p(),"\n";
	print "<a href='focal_image_norm.png'><img src='thumbs/focal_image_norm.jpg'></a>\n";
	print "<a href='focal_image_1.png'><img src='thumbs/focal_image_1.jpg'></a>\n";
	print "<a href='focal_image_1.5.png'><img src='thumbs/focal_image_1.5.jpg'></a>\n";
	print "<a href='focal_image_2.png'><img src='thumbs/focal_image_2.jpg'></a>\n";
	print "<a href='focal_image_2.5.png'><img src='thumbs/focal_image_2.5.jpg'></a>\n";
	print "<a href='focal_image_3.png'><img src='thumbs/focal_image_3.jpg'></a>\n";
	print "<a href='focal_image_3.5.png'><img src='thumbs/focal_image_3.5.jpg'></a>\n";
	print "<a href='focal_image_4.png'><img src='thumbs/focal_image_4.jpg'></a>\n";
	print_html_end($output_cgi);
	close OUTPUT;
	
	select STDOUT;

	my $url_out_dir = $dir;
	$url_out_dir =~ s/$ENV{'DOCUMENT_ROOT'}//g;
	my $hostname = $ENV{'HTTP_ORIGIN'};
	my $final_url = "$ENV{'HTTP_ORIGIN'}/$url_out_dir";
	
	print $q->redirect($final_url);
} else {
    #no valid file handle was available, instead present the user with the
    #regular upload page
	print $q->header,
		  $q->start_html(
			  -title=>'Focal Adhesion Analysis Server - Threshold Testing',
			  -style=>'/FA_webapp/css/screen.css');

	print "<div class=\"container\">\n";
	print $q->h1('Focal Adhesion Analysis Server - Threshold Testing');

    print $q->start_form(-method=>"POST",
                     -action=>"thresh_testing.pl",
                     -enctype=>"multipart/form-data");
    print $q->h3('Adhesion Image File'), 
          $q->filefield('uploaded_file','',50,80);
    
    print $q->br,$q->br,
          $q->submit(-name=>"Submit Data");

    print $q->end_form;

    print $q->br,$q->hr;

	print $q->p, "Most image formats should work in the thresholding algorithm,
	but I've only tested it with tiff (unstacked please) and png files. Also
	keep in mind that in a full run of the processing pipeline, the entire image
	set is used to determine the actual threshold value, here only the image you
	provide is used. Thus, there will be some slight differences in the exact
	threshold selected, but in general I've found them to be very similar.";
	
	print $q->br,$q->br, "Don't have an image to upload as a test? Try this ", 
		$q->a({href=>"/FA_webapp/sample/sample_image.png"},'one'), ".";
	
	print $q->br,$q->br, "See a sample result page ", 
		$q->a({href=>"/FA_webapp/threshold_testing_sample"},'here'), ".";
}

print_html_end($q);

###############################################################################
# Functions
###############################################################################
