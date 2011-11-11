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
		my @files = <$_/*>;
		@files = grep !($_ =~ /index.html/), @files;
		if (scalar(@files) > 1) {
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

	sleep 1;
}

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

