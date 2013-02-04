#!/usr/bin/perl -w

use Config::General qw(ParseConfig);
use Statistics::Descriptive;

use Getopt::Long;
use Geo::IP;

my %opt;
$opt{days} = 7;
GetOptions(\%opt, "days=s") or die;

###############################################################################
# Main
###############################################################################

my $find_results = `find /home/mbergins/Documents/Projects/focal_adhesions/trunk/data/FAAS_*/*.cfg -ctime -$opt{days}`;
my @last_7_day_cfgs = split("\n",$find_results);
my $cfg_count = scalar(@last_7_day_cfgs);

%emails = &get_email_addresses_and_counts(@last_7_day_cfgs);
my @count_sort = sort {$emails{$b}{count} <=> $emails{$a}{count}} keys %emails;

my $email_text = "";

for (@count_sort) {
	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@{$emails{$_}{runtime}});

	$email_text .= sprintf("%s => %d mean (%0.1f)\n", $_, $emails{$_}{count}, $stat->mean());
}

my $subject = "The FAAS server processed $cfg_count experiments in the past $opt{days} days";

system("echo \"$email_text\" | mail -s \"$subject\" matthew.berginski\@gmail.com");

###############################################################################
# Functions
###############################################################################

sub get_email_addresses_and_counts {
	my @cfgs = @_;
	my %emails;
	foreach (@cfgs) {
		my %config = ParseConfig(-ConfigFile => $_, -IncludeRelative => 1);
		my $email_str; 
		if (defined $config{email}) {
			$email_str = $config{email};
		} else {
			$email_str = "No Email";
			# if (defined $config{submitter_ip}) {
			# 	$email_str = "No Email ($config{submitter_ip})";
			# } else {
			# 	$email_str = "No Email (No IP)";
			# }
		}
		
		if (defined $config{submitter_ip}) {
			$email_str .= " ($config{submitter_ip})";
		}

		$emails{$email_str}{count}++;

		if (defined $config{runtime}) {
			push @{$emails{$email_str}{runtime}}, $config{runtime};
		}	
	}

	return %emails;
}
