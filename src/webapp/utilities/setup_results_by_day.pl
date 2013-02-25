#!/usr/bin/perl -w

use POSIX;
use File::Path;
use File::Spec::Functions;

my $target_folder = "results_by_day";

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

@results_folders = <../../../data/*>;

foreach (@results_folders) {
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

