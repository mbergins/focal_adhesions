#!/usr/bin/perl -w

if (scalar(@ARGV) == 0) {
	die "Expected one parameter: an email address";
}

my $email_address = $ARGV[0];

my @df_lines = split(/\n/,`df -h 2>/dev/null`);

for (@df_lines) {
	if ($_ =~ /vg-root.*?(\d+)%/) {
		if ($1 > 90) {
			system("echo \"$_\" | mail -s \"Disk Usage Hit ($1%)\" $email_address");
		}
	}
}

