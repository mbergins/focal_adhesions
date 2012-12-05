#!/usr/bin/perl -w

use File::Basename;
use Config::General qw(ParseConfig);

my $target_email = "mbergins\@unc.edu";

###############################################################################
# Main
###############################################################################

my $find_results = `find /home/mbergins/Documents/Projects/focal_adhesions/trunk/data/FAAS_*/*.cfg`;
my @cfgs = split("\n",$find_results);

%emails = &get_email_addresses_and_counts(@cfgs);
@cfg_files = @{$emails{$target_email}};

die "No Exps found for $target_email" if (scalar(@cfg_files) == 0);

my @exp_ids;
for (@cfg_files) {
	if (basename($_) =~ /(.*)\./) {
		push @exp_ids, $1;
	}
}

my $data_dir = "/home/mbergins/Documents/Projects/focal_adhesions/trunk/data/";
my $processing_results_dir = "/home/mbergins/Documents/Projects/focal_adhesions/trunk/results/";
my $results_dir = "/var/www/FA_webapp/results/";

my $file_list;

for (@exp_ids) {
	$file_list .= "$data_dir" . "$_* ";
	$file_list .= "$processing_results_dir" . "$_* ";
	$file_list .= "$results_dir" . "$_* ";
}

system("rm -rf $file_list");

###############################################################################
# Functions
###############################################################################

sub get_email_addresses_and_counts {
	my @cfgs = @_;
	my %emails;
	foreach (@cfgs) {
		my %config = ParseConfig(-ConfigFile => $_, -IncludeRelative => 1);
		if (defined $config{email}) {
			push @{$emails{$config{email}}}, $_;
		}
	}

	return %emails;
}
