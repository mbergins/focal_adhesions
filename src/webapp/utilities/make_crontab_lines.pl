#!/usr/bin/perl -w

use Cwd;

if (scalar(@ARGV) != 1) {
	die "Specify number of workers as only paramter.";
}

my $number_of_workers = $ARGV[0];

my @lines;
my $worker_num = 0;
for (1..59) {
	push @{$lines[$worker_num]}, $_;
	$worker_num++;
	if ($worker_num > ($number_of_workers - 1)) {
		$worker_num = 0;
	}
}

my $cwd = getcwd();

my $id = 1;

foreach (@lines) {
	print join(",",@$_), " * * * * cd $cwd; ./run_uploaded_exp.pl -ID $id -fullnice\n";
	$id++;
}
