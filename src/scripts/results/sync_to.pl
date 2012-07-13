#!/usr/bin/perl -w

use Getopt::Long;

my %opt;
$opt{debug} = 0;
$opt{server} = "NOSERVER";
$opt{no_time} = 0;
$opt{repeat} = 1;
$opt{delay} = 0;
$opt{exclude} = "";
GetOptions(\%opt,"server=s", "debug|d", "progress", "no_time|nt", "repeat=s", 
	"delay=s", "exclude=s") or die;

my $progress_str = "";
if ($opt{progress}) {
    $progress_str = '--progress ';
}

my $time_str = "time ";
if ($opt{no_time}) {
    $time_str = '';
}

my $rsync_command = "rsync";
if (-e "/nas02/home/m/b/mbergins/bin/rsync") {
	$rsync_command = "/nas02/home/m/b/mbergins/bin/rsync";
}

my $exclude_str = '';
if ($opt{exclude}) {
	for my $i (split(",",$opt{exclude})) {
		$exclude_str .= " --exclude=**$i**";
	}
}

my $command = "$time_str $rsync_command $progress_str $exclude_str " . 
	"--exclude=**data.stor** -a ../../results/* " .
	"$opt{server}:~/Documents/Projects/focal_adhesions/trunk/results/";
if ($opt{server} eq "NOSERVER" || $opt{debug}) {
    print "$command\n";
} else {
    while($opt{repeat} != 0) {
		my $return = system "$command";
	  	print "Done with sync $opt{repeat} Return code: $return\n\n";
	    if ($return == 65280) { 
			print "Caught exit code.";
			last; 
		}
        $opt{repeat}--;
        sleep $opt{delay}*60;
    }
}

