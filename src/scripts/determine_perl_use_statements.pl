#!/usr/bin/perl

use File::Find;

our @files;
find(\&wanted, ('../'));

my %use_statements;
foreach (@files) {
	open INPUT, "$_";
	while (<INPUT>) {
		chomp;
		if (/^use/) {
			my @split_line = split("", $_);
			if ($split_line[$#split_line] eq ";") {
				my $temp = pop @split_line;
			}
			my $joined_line = join("",@split_line);
			$use_statements{$joined_line}++;
		}
	}
	close INPUT;
}


open OUTPUT, ">perl_use_statements.txt";
print OUTPUT join("\n", sort keys %use_statements);
close OUTPUT;


sub wanted {
	if (($_ =~ /.pl$/ || $_ =~ /.plx$/) &&
	 	!($File::Find::name =~ /lib\/perl/)) {
		push @files, $File::Find::name;
	}
}

sub get_use_statements {
}
