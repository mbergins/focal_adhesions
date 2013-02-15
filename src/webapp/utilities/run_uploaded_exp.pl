#!/usr/bin/perl -w

use strict;
use lib "../lib";
use Data::Dumper;
use File::Basename;
use File::Copy;
use File::Path qw(make_path);
use File::Spec::Functions qw(catdir catfile rel2abs);
use Getopt::Long;
use Cwd;
use Config::General qw(ParseConfig);

my %opt;
$opt{ID} = 1;
GetOptions(\%opt, "fullnice", "ID=s","debug|d") or die;

$| = 1;

if ($opt{fullnice}) {
	system("renice -n 20 -p $$ > /dev/null");
	system("ionice -c 3 -p $$");
}

###############################################################################
# Configuration
###############################################################################
my $start_time = time;

my $hostname = "faas.bme.unc.edu";

my $webapp_dir = "../";

my %dir_locations = (
	upload => catdir($webapp_dir, "uploaded_experiments"),
	data_proc => "../../../data/",
	results => "../../../results/",
	public_output => catdir($webapp_dir,"public","results"),
);

foreach (keys %dir_locations) {
	$dir_locations{$_} = rel2abs($dir_locations{$_});
}

if (! -e $dir_locations{public_output}) {
	make_path($dir_locations{public_output});
}

my $run_file = "fa_webapp.$opt{ID}.run";
&process_run_file($run_file);

###############################################################################
# Main
###############################################################################

###########################################################
# Preliminary Setup
###########################################################
my @uploaded_folders = <$dir_locations{upload}/*>;
if (scalar(@uploaded_folders) == 0) {
	&delete_run_file($run_file);
	if ($opt{debug}) {
		print "No new experiments found\n";
	}
	exit;
}

my %oldest_data;
for my $folder (@uploaded_folders) {
	#folders marked with temp haven't been completely setup yet, skip those
	if ($folder =~ /temp$/) {
		next;
	}
	
	if (not defined $oldest_data{upload_folder}) {
		$oldest_data{upload_folder} = $folder;
	} else {
		if (-C $folder > -C $oldest_data{upload_folder}) {
			$oldest_data{upload_folder} = $folder;
		}
	}
}

if ($opt{debug}) {
	print "Found oldest upload file: $oldest_data{upload_folder}\n";
}

if (basename($oldest_data{upload_folder}) =~ /FAAS_(.*)/) {
	$oldest_data{ID} = $1;
}

###########################################################
# Processing
###########################################################

$oldest_data{data_folder} = catdir($dir_locations{data_proc},basename($oldest_data{upload_folder}));
move($oldest_data{upload_folder}, $oldest_data{data_folder});
$oldest_data{cfg_file} = rel2abs(catfile($oldest_data{data_folder},"analysis.cfg"));

my %temp = ParseConfig(
	-ConfigFile => $oldest_data{cfg_file},
	-MergeDuplicateOptions => 1,
	-IncludeRelative       => 1,
);
$oldest_data{cfg} = \%temp;

$oldest_data{results_folder} = rel2abs(catdir($dir_locations{results},basename($oldest_data{upload_folder})));

# &send_start_email(%oldest_file);

&setup_exp(%oldest_data);
&run_processing_pipeline(%oldest_data);
&build_vector_vis(%oldest_data);
$oldest_data{public_zip} = &zip_results_folder(%oldest_data);

###########################################################
# Notifications, Cleanup
###########################################################
if (defined $oldest_data{cfg}{email}) {
	&send_done_email(%oldest_data);
}

&delete_run_file($run_file);
&add_runtime_to_config(\%oldest_data,$start_time);

###############################################################################
# Functions
###############################################################################

sub add_runtime_to_config {
	my %oldest_data = %{$_[0]};
	my $start_time = $_[1];
	
	my $end_time = time;
	my $total_time = $end_time - $start_time;

	open OUTPUT, ">>$oldest_data{cfg_file}" or die $!;
	print OUTPUT "runtime = $total_time\n";
	close OUTPUT;
}

###########################################################
# Run File Processing
###########################################################

sub process_run_file {
	my $run_file = shift @_;
	if (-e $run_file) {
		open INPUT, $run_file;
		my $process_ID = <INPUT>;
		chomp($process_ID);
		close INPUT;

		my $exists = kill 0, $process_ID;
		if ($exists) {
			if ($opt{debug}) {
				print "Found running process\n";
			}
			exit;
		} else {
			unlink $run_file;
		}
	}

	open OUTPUT, ">$run_file" or die "$!";
	print OUTPUT $$;
	close OUTPUT;
}

sub delete_run_file {
	my $run_file = shift;
	unlink $run_file;
}

###########################################################
# Email
###########################################################

sub send_email {
	my %email_data = @_;
	
	if (defined $email_data{self_note}) {
		$email_data{body} = "$email_data{body}\n" . 
			"Your note to yourself about this experiment:\n\n$email_data{self_note}";
	}
	
	my $from_str = "\"From: noreply\@mimir.bme.unc.edu (FAAS Notification)\"";

	my $command = "echo \"$email_data{body}\" | mail -a $from_str -s \"$email_data{subject}\" $email_data{address}";
	# print $command;
	system $command;
}

sub send_start_email {
	my %config = @_;

	my %start_email = (
		'address' => "$config{email}",
		'body' => "Your experiment ($config{name}) has started processing. You will receive an email later with a link to download your results.",
		'subject' => "Your experiment has started processing ($config{name})",
		'self_note' => "$config{self_note}",
	);

	if (defined $config{self_note}) {
		$start_email{self_note} = $config{self_note};
	}

	&send_email(%start_email);
}

sub send_done_email {
	my %config = @_;
	
	my $full_id = basename($oldest_data{upload_folder});
	
	my $body = "Your experiment ($full_id) has finished processing. " . 
		"You can download your results here:\n\n" .
		"http://$hostname/results/$oldest_data{public_zip}\n\n" . 
		"You can find help with understanding the results here:\n\n" .
		"http://$hostname/results_understanding/\n\n";

	my %done_email = (
		'address' => "$oldest_data{cfg}{email}",
		'body' => $body,
		'subject' => "Your experiment has finished processing ($full_id)",
	);

	if (defined $oldest_data{cfg}{note}) {
		$done_email{self_note} = $oldest_data{cfg}{note};
	}

	&send_email(%done_email);
}

sub send_done_text {
	my %config = @_;
	
	my $provider_email;
	if (defined $config{provider}) {
		if ($config{provider} eq "AT&T") {
			$provider_email = 'txt.att.net';
		} elsif ($config{provider} eq "Verizon") {
			$provider_email = 'vtext.com';
		} elsif ($config{provider} eq "Sprint") {
			$provider_email = 'messaging.nextel.com';
		} else {
			print "Unrecognized provider code: $config{provider}\n" if $opt{debug};
			return;
		}
	}

	if (defined $config{phone} && $config{phone} =~ /\d+/) {
		my %text_email = (
			'address' => $config{phone} . "\@$provider_email",
			'body' => "Your exp ($config{name}) has finished. Self note: $config{self_note}",
			'subject' => "",
			'self_note' => undef,
		);
		&send_email(%text_email);
	}
}

###########################################################
# Experiment Processing
###########################################################
sub setup_exp {
	my %oldest_data = @_;
		
	my $starting_dir = getcwd;
	chdir "../../find_cell_features";
	my $command = "./setup_results_folder.pl -cfg $oldest_data{cfg_file}";
	
	if ($opt{debug}) {
		print "Running: $command\n";
	}	
	system($command);
	chdir $starting_dir;
}

sub run_processing_pipeline {
	my %oldest_data = @_;

	my $starting_dir = getcwd;
	chdir "../../scripts";
	
	my $output_status = catfile($oldest_data{results_folder},'run_status.txt');
	my $output_error = catfile($oldest_data{results_folder},'run_error.txt');

	my $command = "./build_all_results.pl -cfg $oldest_data{cfg_file} -exp_filter $oldest_data{ID} > $output_status 2> $output_error";
	if ($opt{debug}) {
		print "Running: $command\n";
	}
	system($command);
	chdir $starting_dir;
}

sub zip_results_folder {
	my %oldest_data = @_;
	
	my $zip_filename = basename($oldest_data{upload_folder}).".zip";
	my $output_zip = catfile(rel2abs($dir_locations{public_output}),$zip_filename);
	
	my $starting_dir = getcwd;
	chdir $dir_locations{results};
	my $command = "zip -q -r $output_zip " . basename($oldest_data{upload_folder});
	if ($opt{debug}) {
		print "Running: $command\n";
	}
	system($command);
	chdir $starting_dir;

	return $zip_filename; 
}

sub build_vector_vis {
	my %oldest_data = @_;

	my $starting_dir = getcwd;
	chdir "../../visualize_cell_features";
	my $command = "./build_vector_vis.pl -cfg $oldest_data{cfg_file} -white_background";
	
	if ($opt{debug}) {
		print "Running: $command\n";
	}
	system($command);
	chdir $starting_dir;
}
