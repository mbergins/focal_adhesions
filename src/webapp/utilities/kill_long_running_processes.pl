#!/usr/bin/perl -w

use Proc::ProcessTable;

$t = new Proc::ProcessTable( 'cache_ttys' => 1 ); 

foreach $p ( @{$t->table} ){
	if ($p->cmndline =~ /perl.*dispatch.fcgi/) {
		my $run_time_min = $p->time/(1000000*60);
		if ($run_time_min >= 15) {
			# print "Found this job to kill: ". $p->pid . $p->cmndline."\n". $run_time_min . "\n";
			kill 'KILL', $p->pid;
		}
	}
}
