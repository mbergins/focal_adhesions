#!/usr/bin/perl -w

use Config::General qw(ParseConfig);
use Statistics::Descriptive;

use Getopt::Long;
use Geo::IP;

my %opt;
$opt{days} = 7;
$opt{email} = "matthew.berginski\@gmail.com";
GetOptions(\%opt, "days=s", "email=s") or die;

###############################################################################
# Main
###############################################################################

my $find_results = `find ../../../data/FAAS_*/*.cfg`;
my @all_cfgs = split("\n",$find_results);

#I have periodically needed to clear old experimental data from the server and
#we didn't transfer any data from the original hardware, but wanted to maintain
#the submission/image count. I've maintained all the config files from 2016
#onward, but need to add in the number of experiments from the original
#hardware:
#
#	Original Hardware: 19319 experiments
my $all_cfg_count = scalar(@all_cfgs) + 19319;

$find_results = `find ../../../data/FAAS_*/*.cfg -ctime -$opt{days}`;
my @recent_cfgs = split("\n",$find_results);
my $cfg_count = scalar(@recent_cfgs);

my $image_count_line = `find ../../../data/ -iregex .*png.* | wc`;
my @image_count = split(/\s+/,$image_count_line);

#All of the image data has been cleaned periodically as well, here are the counts:
#
# Original Hardware: 614705
# 2016: 318392
# 2017: 264893
# 2018: 361919

my $total_images = $image_count[1] + 614705 + 318392 + 264893 + 361919;

my $IP_count_line  = `grep -h ip ../../../data/*/*.cfg | sort | uniq | wc`;
my @IP_count = split(/\s+/,$IP_count_line);
my $total_IPs = $IP_count[1];

%emails = &get_email_addresses_and_counts(@recent_cfgs);
my @count_sort = sort {$emails{$b}{count} <=> $emails{$a}{count}} keys %emails;

my $email_text = "In total $total_images images have been processed in $all_cfg_count experiments from $total_IPs IP adresses.\n\n";

for (@count_sort) {
	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@{$emails{$_}{runtime}});

	$email_text .= sprintf("%s => %d mean (%0.1f)\n", $_, $emails{$_}{count}, $stat->mean());
}

my $disk_usage = `df -h`;

$email_text .= "\n\nDisk Usage\n\n";
$email_text .= $disk_usage;

my $subject = "The FAAS processed $cfg_count experiments in the past $opt{days} days and $all_cfg_count overall";

system("echo \"$email_text\" | mail -s \"$subject\" $opt{email}");

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
