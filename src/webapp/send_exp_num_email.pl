#!/usr/bin/perl -w

use Config::General qw(ParseConfig);

###############################################################################
# Main
###############################################################################

my $find_results = `find /home/mbergins/Documents/Projects/focal_adhesions/trunk/data/FAAS_*/*.cfg -ctime -7`;
my @last_7_day_cfgs = split("\n",$find_results);
my $cfg_count = scalar(@last_7_day_cfgs);

%emails = &get_email_addresses_and_counts(@last_7_day_cfgs);
my @count_sort = sort {$emails{$b} <=> $emails{$a}} keys %emails;

my $email_text = "";

for (@count_sort) {
	$email_text .= "$_ => $emails{$_}\n";
}

my $subject = "The FAAS server processed $cfg_count experiments in the past week";

system("echo \"$email_text\" | mail -s \"$subject\" matthew.berginski\@gmail.com");

###############################################################################
# Functions
###############################################################################

sub get_email_addresses_and_counts {
	my @cfgs = @_;
	my %emails;
	foreach (@cfgs) {
		my %config = ParseConfig(-ConfigFile => $_, -IncludeRelative => 1);
		if (defined $config{email}) {
			$emails{$config{email}}++;
		} else {
			$emails{"No Email"}++;
		}
	}
	

	return %emails;
}
