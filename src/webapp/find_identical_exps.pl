#!/usr/bin/perl -w

# my $md5_results = `md5sum /var/www/FA_webapp/results/*Z*`;
my $md5_results = `md5sum /home/mbergins/Documents/Projects/focal_adhesions/trunk/data/*/*.zip`;

my @md5_lines = split(/\n/,$md5_results);

# die join("\n",@md5_lines);

my %results_md5;
for (@md5_lines) {
	my @results_line = split(/\s+/,$_);
	
	push @{$results_md5{$results_line[0]}}, $results_line[1];
}

for (keys %results_md5) {
	if (scalar(@{$results_md5{$_}} > 1)) {
		print "\n";
		print join("\n",@{$results_md5{$_}}), "\n";
		print "\n";
	}
}
