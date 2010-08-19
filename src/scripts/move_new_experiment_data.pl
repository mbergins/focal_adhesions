#!/usr/bin/perl

################################################################################
# Global Variables and Modules
################################################################################
use lib "../lib";
use lib "../lib/perl";

use threads;
use threads::shared;
use File::Spec::Functions;
use File::Basename;
use File::Find;
use File::Copy;
use Benchmark;
use Getopt::Long;
use Cwd;
use Data::Dumper;

use Config::Adhesions qw(ParseConfig);

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d") or die;

my $base_search_dir = '/tmp/FA_webapp/';
my $base_data_dir = '/Volumes/Data/projects/fa_webapp/data/';

################################################################################
# Main
################################################################################
my @files = <$base_search_dir*>;
if (scalar(@files) == 0) {
    die "No Files to work on" if $opt{debug};
    exit;
}

my %current_file_ages;
foreach (@files) {
    if (/\.cfg/ || /\.working/) {
        next;
    } else {
        if (-e "$_.working") {
            next;
        }
        my($dev, $ino, $mode, $nlink, $uid, $gid, $rdev ,
           $size, $atime, $mtime, $ctime, $blksize, $blocks) = stat($_);
        $current_file_ages{$_} = $mtime;
    }
}
if (scalar keys %current_file_ages == 0) {
    die "No finished uploads to work on" if $opt{debug};
    exit;
}

my @sorted_file_names = sort { $current_file_ages{$a} <=> $current_file_ages{$b} } keys %current_file_ages;

my $lowest_dir_num = &find_lowest_available_directory_number;
print "Lowest dir number available: $lowest_dir_num\n";

my $exp_data_dir = catdir($base_data_dir, $lowest_dir_num);
print "Target Experimental Dir: $exp_data_dir\n";

mkdir $exp_data_dir or die "$!";
open OUTPUT, ">$exp_data_dir/working";
print OUTPUT "wait for it";
close OUTPUT;

my $image_file_directory = &move_image_data($sorted_file_names[0], $exp_data_dir);

my %cfg = &print_exp_config_file($sorted_file_names[0],$exp_data_dir, $image_file_directory);

unlink("$exp_data_dir/working") or die;

&send_files_moved_email(\%cfg);

exit(0);

################################################################################
# Functions
################################################################################

sub find_lowest_available_directory_number {
    my @all_data_directories = <$base_data_dir*>;
    @all_data_directories = grep {
        if (/config/) {
            0;
        } else {
            1;
        }
    } @all_data_directories;
    @all_data_directories = reverse sort @all_data_directories;

    my $lowest_dir_num;
    if (@all_data_directories) {
        $lowest_dir_num = sprintf('%09d', basename($all_data_directories[0]) + 1);
    } else {
        $lowest_dir_num = sprintf('%09d', 1);
    }

    if (-e $lowest_dir_num) {
        die "Lowest dir ($lowest_dir_num) is already taken";
    }
    return $lowest_dir_num;
}

sub move_image_data {
    my $data_file = $_[0];
    my $exp_data_dir = $_[1];
    
    my $target_name = catdir($exp_data_dir, basename($data_file));
    move $data_file, $target_name or die "$!";
    
    system("unzip -q -d $exp_data_dir $target_name");
    unlink $target_name;

    my @data_files = <$exp_data_dir/*>;
    @data_files = grep -d $_, @data_files;
    die "More than one directory found" if (scalar(@data_files) != 1);
    die "File found is not a directory" if (! -d $data_files[0]);

    return basename($data_files[0]);
}

sub print_exp_config_file {
    my $conf_file = new Config::General("$_[0].cfg");
    unlink("$_[0].cfg");

    my %cfg = $conf_file->getall;
 
    $cfg{exp_name} = $lowest_dir_num;
    $cfg{adhesion_image_folder} = $_[2];
    
    my $output_str = $conf_file->save_string(\%cfg);
    
    $output_str = "<<include ../config/webapp_default.cfg>>\n" . $output_str;
    
    my $cfg_file = catfile($_[1], 'analysis.cfg');
    open OUTPUT, ">$cfg_file" or die "$!";
    print OUTPUT $output_str;
    close OUTPUT;
    
    return %cfg;
}

sub send_files_moved_email {
    my %cfg = %{$_[0]};
    
    my $short_exp_name = sprintf('%d', $cfg{exp_name});
    
    my $URL_address = "http://balder.bme.unc.edu/cgi-bin/mbergins/fa_webapp_status.pl?exp_num=$short_exp_name";

    my $body_text = "Your job (#$short_exp_name) has been moved into the processing queue. You can see the status of your job at:\n\n$URL_address\n";
    
    if (defined $cfg{self_note} && $cfg{self_note} ne "") {
        $body_text .= "\nThe note you submitted with this experiment was:\n\n$cfg{self_note}\n";
    }

    $body_text .= "\nThank you for using the focal adhesion analysis server.";
   
    my $system_command = "echo '$body_text' | mail -s 'Focal Adhesion Analysis Job #$short_exp_name' '$cfg{email}'";

    my $return_code = system $system_command;
    print $system_command;
    print "EMAIL RETURN CODE: $return_code";
}
