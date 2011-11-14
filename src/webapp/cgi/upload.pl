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


my $start = time;

$| = 1;

###############################################################################
# Main Program
###############################################################################

my $q = CGI->new();

print $q->header,
	  $q->start_html(-title=>'Focal Adhesion Analysis Server', 
		  -style=>'../../FA_webapp/css/screen.css',
	  );
	  	  # -script=>{-type=>'text/javascript', -src=>'http://ajax.googleapis.com/ajax/libs/jquery/1.7.0/jquery.min.js'},
# 	  	  -script=>{-type=>'text/javascript', -src=>'http://dev.jquery.com/view/trunk/plugins/validate/jquery.validate.js'});
# print "<script>\n",
#   "\$(document).ready(function(){\n",
#   "\$(\"#fa_upload\").validate();\n",
#   "});\n",
#   "</script>\n";

# print "<script type=\"text/javascript\"> </script>";

print "<div class=\"container\">\n";
print $q->h1( 'Focal Adhesion Analysis Server');
#First check for a valid file handle, if so, present user with acknowledgement
#of the upload
my $lightweight_fh = $q->upload('uploaded_file');
if (defined $lightweight_fh) {
    # Upgrade the handle to one compatible with IO::Handle:

    my $io_handle = $lightweight_fh->handle();
    $q->param('uploaded_file') =~ /(.*)/;
    
    my $temp_output = catdir('upload');
    if (! -e $temp_output) {
        mkpath($temp_output);
    }

    my ($output_handle, $output_file) = tempfile('FAAS_XXXXXX',DIR=>catdir('upload')) or die "$!";
	
    my $buffer;
    while (my $bytesread = $io_handle->read($buffer,1024)) {
        print $output_handle $buffer or die;
    }
    close $output_handle;
    print "Done loading file.";

	
	my $end = time;
	
	print $q->p;
	printf "Total Upload Time: %d seconds.", $end - $start;

    #config creation
    my %new_cfg;
    $new_cfg{email} = $q->param('email_address');
    $new_cfg{self_note} = $q->param('self_note');
    $new_cfg{stdev_thresh} = $q->param('filter_thresh');
    $new_cfg{upload_time} = $end-$start;
    $new_cfg{exp_ID} = basename($output_file);
	&output_config(\%new_cfg,"$output_file.cfg");
	
	my $file_type = &determine_file_type($output_file);
	if ($file_type eq 'TIFF') {
		&move_tiff_file($output_file);
	} elsif ($file_type eq 'zip') {
		move($output_file, $output_file . ".zip");
		chmod 0777, "$output_file.zip" or die "$!";
	} else {
		my $file_out = `/usr/bin/file $output_file`;
		print $q->p, "Unexpected file type, got: $file_out\n";
		print $q->p, "Please upload either a zip file or a stacked TIFF, see directions on upload page.";
		unlink $output_file, "$output_file.cfg";
	}
	
	if ($file_type ne 'unknown') {
		print $q->p, 'Thanks for the file, expect an email when the image processing begins.';
	}
} else {
    #no valid file handle was available, instead present the user with the
    #regular upload page
    print $q->start_form(-method=>"POST",
                     -action=>"upload.pl",
                     -enctype=>"multipart/form-data",
				 	 -name=>'fa_upload',
				 	 -onSubmit=>'displaymessage()');
    print $q->h3({class=>'thin'},'Adhesion Image Zip File'), 
          $q->filefield('uploaded_file','',50,80);
    print $q->h3({class=>'thin'},'Your Email Address'),
          $q->textfield('email_address','',50,80);
    print $q->h3({class=>'thin'},'Note to Self About Experiment'),
          $q->textfield('self_note','',50,80);

    print $q->h2('Advanced Settings'), 
	   	  $q->h3({class=>'thin'},'Adhesion Detection Threshold'),
          $q->textfield('filter_thresh','2',50,80);
    
    print $q->br,$q->br,
          $q->submit(-name=>"Submit Data");

    print $q->end_form;

    print $q->br, $q->hr;
    
    &print_analysis_instrucitons
}

print $q->br,$q->br,$q->hr, $q->p, 
	"Return to the ", $q->a({href=>"/FA_webapp"},'home page') ;

print "</div>\n";
print $q->end_html;

###############################################################################
# Functions
###############################################################################

sub print_analysis_instrucitons {
    print $q->h1('Instructions');
    
    print "Thank you for helping to test the focal adhesion analysis webserver.
    Currently the service is only running on one CPU of my lab's server, so it
    might take some time for many jobs to make it through the system. If you
    encounter any problems, feel free to email <a href=mailto:mbergins\@unc.edu>me</a>.";

    print $q->h2('Adhesion Image Zip File');

	print "The program expects that you will submit a single zip file or a
	stacked TIFF. If you upload a zip file, it should contain one folder with
	all the images the images from your experiment.  The images can be in a
	single stack or in separate files.";
    
    print $q->h2('Email Address');

    print "This one is self expanitory, but please put a valid email address
    here. After you submit your files, notification of where to download the
    results will be sent through email.";

    print $q->h2('Note to Self About Experiment');

    print "Whatever you put in this box will be send back to you in any email
    the system sends concerning your experiment. It is limited to 80
    characters.";
        
	print $q->h2('Adhesion Detection Threshold');

	print "This number is used by the adhesion detection script to determine
	when a pixel is or isn't part of an adhesion. After appling a high pass
	filter to the images, pixels above this level are considered part of an
	adhesion, while the pixels below are classified as background. The lower
	this number, the more pixels will be classified as part of an adhesion.  The
	default value of 2 works well in most cases, but values down to around 1 may
	be appropriate. Also be aware that lower values will lengthen the runtime.
	If you want to see what one of your images looks like when processed with a
	specific threshold try out the threshold <a
	href=thresh_testing.pl>tester</a>.";
}

sub output_config {
	my %cfg = %{$_[0]};
	my $target_file = $_[1];
	
	open OUTPUT, ">$target_file" or die $!;
	print OUTPUT "<<include ../../config/webapp_default.cfg>>\n\n";
	foreach (keys %cfg) {
		print OUTPUT "$_ = $cfg{$_}\n";
	}
	close OUTPUT;

    chmod 0777, "$target_file" or die "$!";
}

sub determine_file_type {
	my $file_name = shift;
	$ENV{PATH} = '/usr/bin/';
	my $file_out = `/usr/bin/file $file_name`;

	my $file_type;
	if ($file_out =~ /TIFF/) {
		$file_type = 'TIFF';
	} elsif ($file_out =~ /Zip/) {
		$file_type = 'zip';
	} else {
		$file_type = 'unknown';
	}
	return $file_type;
}

sub move_tiff_file {
	my $file = shift;

	$file = basename($file);

	chdir 'upload';
	mkdir($file);

	move($file,$file);
	$ENV{PATH} = '/usr/bin/';
	system("/usr/bin/zip -q -0 $file.zip $file");
	
	chmod 0777, "$file.zip" or die "$!";
	unlink $file;
	chdir '..';
}
