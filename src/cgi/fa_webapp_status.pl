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

#these two varaibles specify what needs to be removed to allow the proper file
#links to be made
my $file_location_strip = "/Library/WebServer/Documents/";
my $file_location_prefix = "../../";

my $base_data_dir = '/Volumes/Data/projects/fa_webapp/data/';
my $base_results_dir = '/Volumes/Data/projects/fa_webapp/results/';

################################################################################
# Main Program
################################################################################

my $q = CGI->new();

print $q->header,                    # create the HTTP header
	  $q->start_html(-title=>'Focal Adhesion Analysis Server Status Page (Version 0.0.1 Alpha)'),
      $q->h1('Focal Adhesion Analysis Server Status Page');

#Has an experiment number been specified in the address bar?
if (defined $q->param('exp_num')) {
    
    #collect all the files that are available in the webserver file results
    #directory
    my @available_files = </Library/WebServer/Documents/FA_webapp_results/*>;

    #print join("\n",@available_files), $q->p;
    
    #parse out only the numbers from the files found
    my %available_nums;
    for (@available_files) {
        my $file_name = basename($_);
        if (/(\d+)\.zip/) {
            $available_nums{$1} = $_;
        }
    }

    #print join("\n",keys %available_nums), $q->p;
    
    #See if there are any matches with the experiment number specified in the
    #browser
    my @final_matches = grep {
        if ($_ == $q->param('exp_num')) {
            1;
        }
    } keys %available_nums;

    #print join("\n",@matches), $q->p;
    
    #if there were matches in the final results directory, return a page with a
    #link to them
    if (@final_matches) {
        #print $available_nums{$matches[0]};
        
        my $link_address = &convert_local_to_web_filename($available_nums{$final_matches[0]});

        print $q->h2('Your Results are Ready');
        
        print "<A HREF=$link_address>Download Your Results</A>";
        
        my $vis_file = $available_nums{$final_matches[0]};
        $vis_file =~ /(.*)\.zip/;
        $vis_file = $1 . ".png";

        if (-e "$vis_file") {
            print $q->h2("Visualization of Your Experiment");
            print $q->p, "<img src=" . &convert_local_to_web_filename($vis_file) . ">";
        }
    } else {
        #The query experiment number was not found in the final matches, now
        #we see if we can determine how far back in the queue the experiment is
        my @data_dirs = <$base_data_dir*>;
        my @data_dirs = grep !($_ =~ /config/), @data_dirs;
        my @data_basenames = map basename($_), @data_dirs;
        
        my @results_dirs = <$base_results_dir*>;
        my @results_basenames = map basename($_), @results_dirs;

        my @data_matches = grep $_ == $q->param('exp_num'), @data_basenames;
        my @results_matches = grep $_ == $q->param('exp_num'), @results_basenames;
        
        #there weren't any matches in the data directory, this shouldn't happen,
        #but if it does, throw an error
        if (scalar(@data_matches) == 0) {
            print $q->h1("Warning Error");

            print "This page shows the status of individual jobs in the Focal
            Adhesion Analysis Server. Please go back to the email you recieved
            and make sure the entire link is present in the address bar of your
            browser.";
        } elsif (@results_matches) {
            print $q->h1("Getting Close");

            print "Your job is currently being processed in the queue, when it
            is done you should expect an email.";
        } else {
            #we know the experiment is in the data directory and isn't in the
            #results directory yet, so lets figure out how far down the queue it
            #is
            my %current_dir_ages;
            foreach my $this_dir (@data_dirs) {
                if (-d $this_dir) {
                    if (-e catfile($this_dir, "working")) {
                        next;
                    }
                    if (grep $this_dir =~ /$_/, @results_basenames) {
                        next;
                    }
                    my($dev, $ino, $mode, $nlink, $uid, $gid, $rdev ,
                       $size, $atime, $mtime, $ctime, $blksize, $blocks) = stat($this_dir);
                    $current_dir_ages{$this_dir} = $mtime;
                }
            }
            my @sorted_dir_names = sort { $current_dir_ages{$a} <=> $current_dir_ages{$b} } keys %current_dir_ages;
            
            #print join("\n", @sorted_dir_names), $q->p;
            #print join("\n", @data_basenames), $q->p;
            #print join("\n", @results_basenames), $q->p;
            
            my @match_positions = grep basename($sorted_dir_names[$_]) == $q->param('exp_num'), 
                (0..$#sorted_dir_names);
            $match_positions[0]++;

            print $q->h2("Your job is in queue position number $match_positions[0]");
        }
    }
} else {
    #no experiment number was specified in the address bar, return an error
    print $q->h1("Warning Error");

    print "This page shows the status of individual jobs in the Focal Adhesion
    Analysis Server. Please go back to the email you recieved and make sure the
    entire link is present in the address bar of your browser.";
}

print $q->end_html;                  # end the HTML

################################################################################
# Functions
################################################################################

sub convert_local_to_web_filename {
    my $link_address = $_[0];
    if ($link_address =~ /$file_location_strip(.*)/) {
        $link_address = $1;
    } else {
        print $q->h1('File name creation error');
    }

    #print $q->h1($link_address);
    $link_address = $file_location_prefix . $link_address;
    #print $q->h1($link_address);
    return $link_address;
}
