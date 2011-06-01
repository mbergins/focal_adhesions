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

###############################################################################
# Main
###############################################################################

my $q = CGI->new();

print $q->header,
	  $q->start_html(-title=>'Focal Adhesion Alignment Index Server'), 
	  $q->h1('Focal Adhesion Alignment Index Server');

my $lightweight_fh = $q->upload('uploaded_file');
# undef may be returned if it's not a valid file handle
if (defined $lightweight_fh) {
    # Upgrade the handle to one compatible with IO::Handle:

    my $io_handle = $lightweight_fh->handle();
    $q->param('uploaded_file') =~ /(.*)/;
    
    my ($output_handle, $output_file) = File::Temp::tempfile(DIR=>$upload_dir) or die "$!";
    
    my %new_cfg;
    $new_cfg{email} = $q->param('email_address');
    $new_cfg{self_note} = $q->param('self_note');
    $new_cfg{exp_id} = basename($output_file);
    if ($q->param('color_blind')) {
        $new_cfg{color_blind} = 1;
    } else {
        $new_cfg{color_blind} = 0;
    }
    $new_cfg{stdev_thresh} = $q->param('stdev_thresh');
    $new_cfg{min_adhesion_size} = $q->param('min_adhesion_size');
    $new_cfg{min_axial_ratio} = $q->param('min_axial_ratio');
    my $conf = new Config::General(\%new_cfg);
    $conf->save_file("$output_file.cfg");
    chmod 0666, "$output_file.cfg" or die "$!";
    
    #print $q->h2('Data Loaded So Far');
    my $data_read = 0;
    my $buffer;
    while (my $bytesread = $io_handle->read($buffer,1024)) {
        print $output_handle $buffer or die;
        $data_read++;
    }
    close $output_handle;
    chmod 0666, "$output_file" or die "$!";
    
    print $q->p, 'Thanks for the file, you can expect two emails: one when your ' .
    'images start processing and another when they finish.';
} else {
    print $q->start_form(-method=>"POST",
                     -action=>basename($0),
                     -enctype=>"multipart/form-data");
    print $q->h2('Required Options');
    print $q->h3('Adhesion Image Zip File'), 
          $q->filefield('uploaded_file','',50,80);
    print $q->h3('Your Email Address'),
          $q->textfield('email_address','',50,80);
    print $q->h3('Note to Self About Experiment'),
          $q->textfield('self_note','',50,80);
    
    print $q->h2('Adhesion Segmentation Options');
    print $q->h3('Identifcation Threshold (generally between 2-4)'),
          $q->textfield('stdev_thresh',2,10,10);
    print $q->h3('Minimum Adhesion Size (in pixels, generally either 2 or 3)'),
          $q->textfield('min_adhesion_size',3,10,10);
    print $q->h3('Minimum Major/Minor Axis Ratio (generally 3 or greater)'),
          $q->textfield('min_axial_ratio',3,10,10);

    print $q->h2('Miscelleous Options');

    print $q->checkbox(-name=>'color_blind',
                       -checked=>0,
                       -label=>'Are You Color Blind?');;
    
    print $q->p,
          $q->submit(-name=>"Submit Data");

    print $q->end_form;

    print $q->hr;

    print $q->h1('Instructions');
    
    print $q->h2('Adhesion Image Zip File');

    print "The program expects that you will submit a single zip file. Inside
    the zip file will be one folder containing all the images in your
    experiment.  The images can be in a single stack or in separate files.";
    
    print $q->p;

    print $q->h2('Email Address');

    print "After you submit your files, notification of where to download the
    results will be sent through email.";

    print $q->h2('Note to Self About Experiment');

    print "Whatever you put in this box will be send back to you in any email
    the system sends concerning your experiment. It is limited to 80
    characters.";
    
    if ($q->param('advanced') == 1) {
    }
}

print $q->end_html;
