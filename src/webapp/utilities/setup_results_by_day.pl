#!/usr/bin/perl -w

use POSIX;
use File::Path;
use File::Basename;
use File::Spec::Functions;
use Config::General;
use Data::Dumper;

my $target_folder;

###############################################################################
# Main
###############################################################################

###########################################################
# Results Files Organization
###########################################################

$target_folder = "results_by_day";

File::Path::remove_tree($target_folder);

my @results_folders = <../public/results/*>;

foreach (@results_folders) {
	my $date = &POSIX::strftime("%Y-%m-%d", localtime( ( stat $_ )[9]));
	mkpath("$target_folder/$date");
	my $abs_target = &File::Spec::Functions::rel2abs($_);
	system("ln -s $abs_target $target_folder/$date");

	$date = &POSIX::strftime("%Y-%m", localtime( ( stat $_ )[9]));
	mkpath("$target_folder/$date");
	system("ln -s $abs_target $target_folder/$date");

	$date = &POSIX::strftime("%Y", localtime( ( stat $_ )[9]));
	mkpath("$target_folder/$date");
	system("ln -s $abs_target $target_folder/$date");
}

##############################
# Session ID 
##############################
&make_results_links_by_config_var("session_user_id");

###########################################################
# Data Folder Organization
###########################################################
$target_folder = "data_by_day";

File::Path::remove_tree($target_folder);

@data_folders = <../../../data/*>;

foreach (@data_folders) {
	next if ($_ =~ /config/);
	my $date = &POSIX::strftime("%Y-%m-%d", localtime( ( stat $_ )[9]));
	mkpath("$target_folder/$date");
	my $abs_target = &File::Spec::Functions::rel2abs($_);
	system("ln -s $abs_target $target_folder/$date");

	$date = &POSIX::strftime("%Y-%m", localtime( ( stat $_ )[9]));
	mkpath("$target_folder/$date");
	system("ln -s $abs_target $target_folder/$date");

	$date = &POSIX::strftime("%Y", localtime( ( stat $_ )[9]));
	mkpath("$target_folder/$date");
	system("ln -s $abs_target $target_folder/$date");
}

##############################
# IP Address Folder Organization
##############################
&make_data_links_by_config_var("submitter_ip");

##############################
# Email Address Folder Organization
##############################
&make_data_links_by_config_var("email");

##############################
# Session ID Folder Organization
##############################
&make_data_links_by_config_var("session_user_id");

###############################################################################
# Functions
###############################################################################

sub make_data_links_by_config_var {
	my $cfg_var = $_[0];
	
	my $target_folder = "data_by_$cfg_var";

	File::Path::remove_tree($target_folder);

	@data_folders = <../../../data/*>;

	foreach (@data_folders) {
		next if ($_ =~ /config/);
		my @cfg = <$_/*.cfg>;
		$conf = new Config::General((-ConfigFile => $cfg[0],
				-IncludeRelative => 1));
		my %config = $conf->getall;
		
		my $abs_target = &File::Spec::Functions::rel2abs($_);
		if (defined $config{$cfg_var}) {
			$config{$cfg_var} =~ s/^\s+//;
			mkpath("$target_folder/$config{$cfg_var}");
			system("ln -s \"$abs_target\" \"$target_folder/$config{$cfg_var}\"");
		} else {
			mkpath("$target_folder/no_$cfg_var");
			system("ln -s \"$abs_target\" \"$target_folder/no_$cfg_var\"");
		}
	}
}

sub make_results_links_by_config_var {
	my $cfg_var = $_[0];
	
	$target_folder = "results_by_$cfg_var";

	File::Path::remove_tree($target_folder);

	my @results_files = <../public/results/*>;

	@data_folders = <../../../data/*>;

	foreach (@data_folders) {
		next if ($_ =~ /config/);
		my @cfg = <$_/*.cfg>;
		$conf = new Config::General((-ConfigFile => $cfg[0],
				-IncludeRelative => 1));
		my %config = $conf->getall;
		
		my $exp_name = basename($_);

		my @results_files = grep /$exp_name/, @results_files;
		if (scalar(@results_files) > 1) {
			die "Problem matching $exp_name";
		}
		if (scalar(@results_files) == 0) {
			next;
		}

		my $abs_target = &File::Spec::Functions::rel2abs($results_files[0]);
		if (defined $config{$cfg_var}) {
			$config{$cfg_var} =~ s/^\s+//;
			mkpath("$target_folder/$config{$cfg_var}");
			system("ln -s \"$abs_target\" \"$target_folder/$config{$cfg_var}\"");
		} else {
			mkpath("$target_folder/no_$cfg_var");
			system("ln -s \"$abs_target\" \"$target_folder/no_$cfg_var\"");
		}
	}
}
