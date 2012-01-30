#!/usr/bin/perl -w

use Getopt::Long;

my %opt;
$opt{debug} = 0;
$opt{server} = "NOSERVER";
$opt{no_time} = 0;
$opt{repeat} = 1;
$opt{delay} = 0;
GetOptions(\%opt,"server=s", "debug|d", "progress", "no_time|nt", "repeat=s", "delay=s") or die;

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

my $command = "$time_str" . "$rsync_command $progress_str" . "-a --exclude data.stor ../../results/* $opt{server}:~/Documents/Projects/focal_adhesions/trunk/results/";
if ($opt{server} eq "NOSERVER" || $opt{debug}) {
    print "$command\n";
} else {
    while($opt{repeat} != 0) {
        system "$command";
        print "Done with sync $opt{repeat}\n\n";
        $opt{repeat}--;
        sleep $opt{delay}*60;
    }
}

