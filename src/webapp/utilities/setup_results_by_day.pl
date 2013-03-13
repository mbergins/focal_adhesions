#!/usr/bin/perl -w

use POSIX;
use File::Path;
use File::Spec::Functions;
use Config::General;

my $target_folder;

$target_folder = "results_by_day";

File::Path::remove_tree($target_folder);

my @results_folders = <../../../results/*>;

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

###########################################################
# IP Address Folder Organization
###########################################################
&make_data_links_by_config_var("submitter_ip");

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
			mkpath("$target_folder/$config{$cfg_var}");
			system("ln -s $abs_target $target_folder/$config{$cfg_var}");
		} else {
			mkpath("$target_folder/no_$cfg_var");
			system("ln -s $abs_target $target_folder/no_$cfg_var");
		}
	}
}
