#!/usr/bin/perl -w

use strict;
use lib "../lib";
use Data::Dumper;
use File::Basename;
use File::Copy;
use File::Spec::Functions qw(catdir catfile);
use Getopt::Long;
use Cwd;

use Config::Adhesions qw(ParseConfig);

my %opt;
$opt{ID} = 1;
GetOptions(\%opt, "ID=s","debug|d") or die;

if ($opt{fullnice}) {
	system("renice -n 20 p $$");
	system("ionice -c 3 p $$");
}

###############################################################################
# Configuration
###############################################################################
my $start_time = time;

my $upload_dir = "/usr/lib/cgi-bin/FA_webapp/upload/";
my $data_proc_dir = "../../data/fa_webapp/";

my $public_output_folder = "/var/www/FA_webapp/results/";
my $run_file = "fa_webapp.$opt{ID}.run";
&process_run_file($run_file);

###############################################################################
# Main
###############################################################################

###########################################################
# Preliminary Setup
###########################################################
my %upload_data;
my @uploaded_files = <$upload_dir/*.zip>;
if (scalar(@uploaded_files) == 0) {
	&delete_run_file($run_file);
	if ($opt{debug}) {
		print "No new experiments found\n";
	}
	exit;
}

for my $file (@uploaded_files) {
	my $name = basename($file);
	$name =~ s/\.zip//;
	$upload_data{$name}{age} = -C $file;
	$upload_data{$name}{source} = $file;

	my $cfg = $file;
	$cfg =~ s/zip/cfg/;
	$upload_data{$name}{orig_cfg} = $cfg;
}

# my @data_names = <$data_proc_dir/*>;
# @data_names = grep -d $_, @data_names;
# @data_names = map basename($_), @data_names;

my %youngest_file = &find_oldest_new_upload_file(\%upload_data);
$youngest_file{target_dir} = catdir($data_proc_dir,$youngest_file{name});

%youngest_file = &setup_and_unzip_file(%youngest_file);

&add_image_dir_to_config(%youngest_file);

my %config = ParseConfig(\%youngest_file);
$youngest_file{email} = $config{email};
$youngest_file{self_note} = $config{self_note};

###########################################################
# Processing
###########################################################

&send_start_email(%youngest_file);
&setup_exp(%youngest_file);
&run_processing_pipeline(\%youngest_file,\%config);
&build_vector_vis(%youngest_file);
$youngest_file{public_zip} = &zip_results_folder(\%youngest_file,\%config);
&send_done_email(%youngest_file);
&delete_run_file($run_file);

&add_runtime_to_config(\%youngest_file,$start_time);
###############################################################################
# Functions
###############################################################################

sub find_oldest_new_upload_file {
	my %upload_data = %{$_[0]};

	my @sorted_files = sort {
		if ($upload_data{$a}{age} < $upload_data{$b}{age}) {
			1;
		} else {
			0;
		}
	} keys %upload_data;
	
	my %youngest_file = %{$upload_data{$sorted_files[0]}};
	$youngest_file{name} = $sorted_files[0];
	return %youngest_file;
}

sub setup_and_unzip_file {
	my %file_data = @_;
	
	mkdir $file_data{target_dir};
	
	$file_data{output_zip} = catfile($file_data{target_dir},"$file_data{name}.zip");
	$file_data{output_cfg} = catfile($file_data{target_dir},"$file_data{name}.cfg");
	$file_data{cfg} = catfile($file_data{target_dir},"$file_data{name}.cfg");

	move($file_data{source},$file_data{output_zip}) or die "$!";
	move($file_data{orig_cfg},$file_data{output_cfg}) or die "$!";
	
	system("unzip -q -d $file_data{target_dir} $file_data{output_zip}");

	return %file_data;
}

sub determine_image_folder {
	my $target_dir = $_[0];
	
	my @files = <$target_dir/*>;
	my @dirs = grep -d $_, @files;
	@dirs = map {
		if ($_ =~ /$target_dir(.*)/) {
			$1;
		}
	} @dirs;

	if (scalar(@dirs) > 1) {
		print "Found more than one folder after unzipping:";
		print join ("\n", @dirs);
		exit;
	}
	if (scalar(@dirs) == 0) {
		print "Didn't find any folders after unzipping";
		exit;
	}

	return $dirs[0];
}

sub add_image_dir_to_config {
	my %file_data = @_;

	my $image_folder = &determine_image_folder($file_data{target_dir});

	open OUTPUT, ">>$file_data{output_cfg}" or die $!;
	print OUTPUT "adhesion_image_folder = $image_folder\n";
	close OUTPUT;
}

sub add_runtime_to_config {
	my %file_data = %{$_[0]};
	my $start_time = $_[1];
	
	my $end_time = time;
	my $total_time = $end_time - $start_time;

	open OUTPUT, ">>$file_data{output_cfg}" or die $!;
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
	
	$email_data{body} = "$email_data{body}\n\n" . 
		"Your note to yourself about this experiment:\n$email_data{self_note}";

	my $command = "echo \"$email_data{body}\" | mail -s \"$email_data{subject}\" $email_data{address}";

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

	&send_email(%start_email);
}

sub send_done_email {
	my %config = @_;
	
	my $striped_output = $public_output_folder;
	$striped_output =~ s#/var/www/##;
	
	my $body = "Your experiment ($config{name}) has finished processing. " . 
		"You can download your results here:\n\n" .
		"http://snotra.bme.unc.edu/$striped_output/$config{public_zip}\n\n";

	my %done_email = (
		'address' => "$config{email}",
		'body' => $body,
		'subject' => "Your experiment has finished processing ($config{name})",
		'self_note' => "$config{self_note}",
	);

	&send_email(%done_email);
}

###########################################################
# Experiment Processing
###########################################################

sub setup_exp {
	my %youngest_file = @_;
		
	my $starting_dir = getcwd;
	chdir "../find_cell_features";
	my $command = "./setup_results_folder.pl -cfg $youngest_file{output_cfg}";
	
	if ($opt{debug}) {
		print "$command\n";
	} else {
		system($command);
	}
		system($command);
	chdir $starting_dir;
}

sub run_processing_pipeline {
	my %youngest_file = %{$_[0]};
	my %config = %{$_[1]};

	my $starting_dir = getcwd;
	chdir "../scripts";
	
	my $output_status = catfile($config{exp_results_folder},'run_status.txt');
	my $output_error = catfile($config{exp_results_folder},'run_error.txt');

	my $command = "./build_all_results.pl -cfg $youngest_file{output_cfg} -exp_filter $youngest_file{name} > $output_status 2> $output_error";
	if ($opt{debug}) {
		print "$command\n";
	} else {
		system($command);
	}
	chdir $starting_dir;
}

sub zip_results_folder {
	my %youngest_file = %{$_[0]};
	my %config = %{$_[1]};
	
	my $output_zip = catfile($public_output_folder,"$youngest_file{name}.zip");
	
	my $starting_dir = getcwd;
	chdir $config{results_folder};
	my $command = "zip -q -r $output_zip $youngest_file{name}";
	if ($opt{debug}) {
		print "$command\n";
	} else {
		system($command);
	}
	chdir $starting_dir;
	return "$youngest_file{name}.zip"; 
}

sub build_vector_vis {
	my %youngest_file = @_;
		
	my $starting_dir = getcwd;
	chdir "../visualize_cell_features";
	my $command = "./build_vector_vis.pl -cfg $youngest_file{output_cfg}";
	
	if ($opt{debug}) {
		print "$command\n";
	} else {
		system($command);
	}
		system($command);
	chdir $starting_dir;
}
