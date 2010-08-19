#!/usr/bin/perl -wT

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

my $q = CGI->new();

print $q->header,                    # create the HTTP header
	  $q->start_html(-title=>'Focal Adhesion Analysis Server (Version 0.0.1 Alpha)'), # start the HTML
	  $q->h1('Focal Adhesion Analysis Server (Version 0.0.1 Alpha)');         # level 1 header

my $lightweight_fh = $q->upload('uploaded_file');
# undef may be returned if it's not a valid file handle
if (defined $lightweight_fh) {
    # Upgrade the handle to one compatible with IO::Handle:

    my $io_handle = $lightweight_fh->handle();
    $q->param('uploaded_file') =~ /(.*)/;

    my ($output_handle, $output_file) = File::Temp::tempfile(DIR=>catdir('/tmp','FA_webapp'));
    
    open OUTPUT, ">$output_file.working";
    print OUTPUT "wait for it";
    close OUTPUT;
    chmod 0666, "$output_file.working" or die "$!";
    
    my %new_cfg;
    $new_cfg{email} = $q->param('email_address');
    $new_cfg{self_note} = $q->param('self_note');
    if (defined $q->param('filter_thresh')) {
        $new_cfg{filter_thresh} = $q->param('filter_thresh');
    }
    my $conf = new Config::General(\%new_cfg);
    $conf->save_file("$output_file.cfg");
    chmod 0666, "$output_file.cfg" or die "$!";
    
    #print $q->h2('Data Loaded So Far');
    print $q->h2('Working...');
    my $data_read = 0;
    my $buffer;
    while (my $bytesread = $io_handle->read($buffer,1024)) {
        print $output_handle $buffer or die;
        $data_read++;
        if ($data_read % (1024*5) == 0) {
            #print $data_read/1024, " megs read in so far.";
            #print $q->br;
        }
    }
    close $output_handle;
    print "Done loading file.";
    chmod 0666, "$output_file" or die "$!";
    
    unlink "$output_file.working";
    
    print $q->p, 'Thanks for the file, expect an email soon with further information about tracking your job.';
} else {
    print $q->start_form(-method=>"POST",
                     -action=>"fa_webapp_upload.pl",
                     -enctype=>"multipart/form-data");
    print $q->h2('Adhesion Image Zip File'), 
          $q->filefield('uploaded_file','',50,80);
    print $q->h2('Your Email Address'),
          $q->textfield('email_address','',50,80);
    print $q->h2('Note to Self About Experiment'),
          $q->textfield('self_note','',50,80);

    if ($q->param('advanced') == 1) {
        print $q->h2('Advanced Settings - See Notes Below for Description'), $q->h2('Adhesion Detection Threshold'),
              $q->textfield('filter_thresh','0.1',50,80);
    }
    print $q->p,
          $q->submit(-name=>"Submit Data");

    print $q->end_form;

    print $q->hr;

    print $q->h1('Instructions');
    
    print "Thank you for helping to test the focal adhesion analysis webserver.
    Currently the service is only running on one CPU of my lab's server, so it
    might take some time for many jobs to make it through the system. If you
    encounter any problems, feel free to email <a href=mailto:mbergins\@unc.edu>me</a>.";

    print $q->h2('Adhesion Image Zip File');

    print "The program expects that you will submit a single zip file. Inside
    the zip file will be one folder containing all the images in your
    experiment.  The images can be in a single stack or in separate files.";
    
    print $q->p;

    print "Zip files can be made in Windows XP by right clicking on the folder
    that directly contains your experimental images and selecting the 'send to'
    menu.  Inside that menu, is the option 'Compressed (zipped) Folder'.";
    
    print $q->h2('Email Address');

    print "This one is self expanitory, but please put a valid email address
    here. After you submit your files, notification of where to download the
    results will be sent through email.";

    print $q->h2('Note to Self About Experiment');

    print "Whatever you put in this box will be send back to you in any email
    the system sends concerning your experiment. It is limited to 80
    characters.";
    
    if ($q->param('advanced') == 1) {
        print $q->h2('Adhesion Detection Threshold');

        print "This number is used by the adhesion detection script to determine
        when a pixel is or isn't part of an adhesion. After appling a high pass
        filter to the images, pixels above this level are considered part of an
        adhesion, while the pixels below are classified as background. The lower
        this number, the more pixels will be classified as part of an adhesion.
        The default value of 0.1 works well in most cases, but values down to
        around 0.05 may be appropriate for image sets with more subtle
        differences in the fluoresence intensities between adhesions and
        background. Also be aware that lower values will lengthen the runtime.";
    }
}

print $q->end_html;                  # end the HTML
