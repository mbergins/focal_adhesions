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

my $start = time;

$| = 1;

###############################################################################
# Configuration
###############################################################################

my $prefix = '/var/www/';
my $output_dir = '/var/www/FA_webapp/threshold_testing';
if (! -e $output_dir) {
	mkpath($output_dir);
}

###############################################################################
# Main Program
###############################################################################

my $q = CGI->new();

#First check for a valid file handle, if so, present user with acknowledgement
#of the upload
my $lightweight_fh = $q->upload('uploaded_file');
if (defined $lightweight_fh) {
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
	
	my $output_cgi = CGI->new();

	open OUTPUT, ">$dir/index.html";

	print OUTPUT $output_cgi->start_html(-title=>'Focal Adhesion Analysis Server - Threshold Testing'), 
		  $output_cgi->h1('Focal Adhesion Analysis Server - Threshold Testing');
	
	#print OUTPUT "matlab -nojvm -nodisplay -nosplash -r \"build_thresholded_set('$output_file');exit;\"";
	print OUTPUT "Can't see any images below? Please wait about 10 seconds and refresh the webpage.";
	
	print OUTPUT $output_cgi->p(),"\n";
	print OUTPUT "<a href='focal_image_norm.png'><img src='thumbs/focal_image_norm.jpg'></a>\n";
	print OUTPUT "<a href='focal_image_1.png'><img src='thumbs/focal_image_1.jpg'></a>\n";
	print OUTPUT "<a href='focal_image_1.5.png'><img src='thumbs/focal_image_1.5.jpg'></a>\n";
	print OUTPUT "<a href='focal_image_2.png'><img src='thumbs/focal_image_2.jpg'></a>\n";
	print OUTPUT "<a href='focal_image_2.5.png'><img src='thumbs/focal_image_2.5.jpg'></a>\n";
	print OUTPUT "<a href='focal_image_3.png'><img src='thumbs/focal_image_3.jpg'></a>\n";
	print OUTPUT "<a href='focal_image_3.5.png'><img src='thumbs/focal_image_3.5.jpg'></a>\n";
	print OUTPUT "<a href='focal_image_4.png'><img src='thumbs/focal_image_4.jpg'></a>\n";
	print OUTPUT $output_cgi->end_html;
	close OUTPUT;
	
	my $url_out_dir = $dir;
	$url_out_dir =~ s/$prefix//g;
	my $final_url = "http://snotra.bme.unc.edu/$url_out_dir";
	
	sleep 10;

	print $q->redirect($final_url);
} else {
    #no valid file handle was available, instead present the user with the
    #regular upload page
	print $q->header,
		  $q->start_html(-title=>'Focal Adhesion Analysis Server - Threshold Testing',-style=>'../../FA_webapp/css/screen.css');

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

	print "Most image formats should work in the thresholding algorithm, but
	I've only tested it with tiff (unstacked please) and png files. Also keep
	in mind that in a full run of the processing pipeline, the entire image set
	is used to determine the actual threshold value, here only the image you
	provide is used. Thus, there will be some slight differences in the exact
	threshold selected, but in general I've found them to be very similar.";

	print $q->br,$q->br,$q->hr, $q->p, 
		"Return to the ", $q->a({href=>"/FA_webapp"},'home page');
	print "</div>\n";
}

print $q->end_html;

###############################################################################
# Functions
###############################################################################
