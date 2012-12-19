#!/usr/bin/perl -w

use File::Basename;
use Config::General qw(ParseConfig);
use Getopt::Long;

my %opt;
GetOptions(\%opt, "debug|d", "configs=s") or die;

my @cfgs_to_remove;
if (defined $opt{configs}) {
	@cfgs_to_remove = split(",",$opt{configs});
}

my %search_targets = (
	email => "matthew.berginski\@gmail.com",
	submitter_ip => "152.19.101.236",
);

###############################################################################
# Main
###############################################################################

my $find_results = `find /home/mbergins/Documents/Projects/focal_adhesions/trunk/data/FAAS_*/*.cfg`;
my @cfgs = split("\n",$find_results);

my @cfg_hits;
if (@cfgs_to_remove) {
	foreach my $this_cfg (@cfgs) {
		if (grep $this_cfg =~ /$_/, @cfgs_to_remove) {
			push @cfg_hits, $this_cfg;
		}
	} 
} else {
	@cfg_hits = &search_targets(@cfgs);
}

die "No Exps found." if (scalar(@cfg_hits) == 0);

my @exp_ids;
for (@cfg_hits) {
	if (basename($_) =~ /(.*)\./) {
		push @exp_ids, $1;
	}
}

if ($opt{debug}) {
	print "Found these configs:\n" . join("\n",@exp_ids) ;
	die;
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

sub get_email_addresses_and_ids {
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

sub search_targets {
	my @cfgs = @_;

	my %hits;
	foreach (@cfgs) {
		my %config = ParseConfig(-ConfigFile => $_, -IncludeRelative => 1);
		foreach my $id (keys %search_targets) {
			if (defined $config{$id} && 
				$config{$id} =~ /$search_targets{$id}/) {
				$hits{$_}++;
			}
		}
	}
	
	my @hits = keys %hits;
	return @hits;
}
