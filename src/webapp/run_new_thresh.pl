#!/usr/bin/perl -w

use Getopt::Long;
my %opt;
GetOptions(\%opt,"debug|d") or die;

my $run_file = 'thresh_testing.run';
&process_run_file($run_file);

my $threshold_dir = '/var/www/FA_webapp/threshold_testing/';
while (1) {
	my @thresh_dirs = <$threshold_dir/*>;
	@thresh_dirs = grep -d $_, @thresh_dirs;

	for (@thresh_dirs) {
		print "Checking: $_\n" if $opt{debug};
		my @files = <$_/*>;
		#check for the index.html file, if present, then the file has been fully uploaded.
		if (! grep $_ =~ /upload_done/, @files) {
			print "Didn't find upload_done, exiting\n" if $opt{debug};
			next;
		}
		#now remove the index.html file, if only one file is left it must be the
		#uploaded file, so run the matlab script on it
		@files = grep !($_ =~ /upload_done/), @files;
		if (scalar(@files) > 1) {
			print "Found more than one other file.\n" if $opt{debug};
			next;
		} else {
			my $command = "matlab -nojvm -nodisplay -nosplash -r \"build_thresholded_set('$files[0]');exit;\" > /dev/null";
			if ($opt{debug}) {
				print "$command\n";
			} else {
				system($command);
			}
		}
	}
	if ($opt{debug}) {
		exit;
	}
	sleep 1;
}

###############################################################################
# Functions
###############################################################################

sub process_run_file {
    my $run_file = shift @_; 
    if (-e $run_file) {
        open INPUT, $run_file;
        my $process_ID = <INPUT>;
        chomp($process_ID);
        close INPUT;

        my $exists = kill 0, $process_ID;
        if ($exists) {
            if ($opt{debug}) {
                print "Found running process\n";
            }   
            exit;
        } else {
            unlink $run_file;
        }   
    }   

    open OUTPUT, ">$run_file" or die "$!";
    print OUTPUT $$; 
    close OUTPUT;
}

sub delete_run_file {
    my $run_file = shift;
    unlink $run_file;
}

