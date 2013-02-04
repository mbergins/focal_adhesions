#!/usr/bin/perl -w

my @df_lines = split(/\n/,`df -h 2>/dev/null`);

for (@df_lines) {
	if ($_ =~ /sda1.*(\d+)%/) {
		if ($1 > 90) {
			system("echo \"$_\" | mail -s \"Disk Usage Hit\" matthew.berginski\@gmail.com");
		}
	}
}

