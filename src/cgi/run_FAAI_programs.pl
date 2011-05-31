#!/usr/bin/perl -w

################################################################################
# Global Variables and Modules
################################################################################

use lib "../lib";

use File::Copy;
use File::Basename;
use Data::Dumper;
use Cwd;
use Getopt::Long;

use Config::General;

my %opt;
$opt{debug} = 0;
$opt{worker_ID} = 1;
GetOptions(\%opt, "debug|d", "worker_ID=s") or die;

my $PID_file = "/var/run/FAAI_webapp/run_$opt{worker_ID}.pid";
my $data_upload_dir = "/tmp/FAAI_webapp/";

################################################################################
# Main
################################################################################

############################################################
# Initial Checks
############################################################
if (&processing_occuring) {
    # print "Processing ongoing\n";
    exit();
}

my @upload_files = </tmp/FAAI_webapp/*>;
if (scalar(@upload_files) == 0) {
    # print "Didn't find any files to process, quitting.\n";
    exit();
}

############################################################
# Find the youngest config file
############################################################
&make_processing_PID_file;  

my %files_ages;
for (@upload_files) {
    if ($_ =~ /.cfg/) {
        $files_ages{$_} = (stat($_))[9];
    }
}

my $youngest_cfg = (sort { $files_ages{$a} <=> $files_ages{$b} } keys %files_ages)[0];
print "Found youngest config $youngest_cfg\n";

my $conf = new Config::General($youngest_cfg);
my %cfg = $conf->getall;

############################################################
# Setup the data folder and process the data
############################################################
&send_processing_started_email;

my $config_file = &move_uploaded_data($youngest_cfg);

my @command_seq = (
    ['../find_cell_features','setup_results_folder.pl'],
    ['../find_cell_features','run_matlab_over_field.pl -script find_fa_angles.m'],
    ['../analyze_cell_features','build_alignment_models.pl'],
    );

print "\nStarting on processing pipeline:\n\n";

my $starting_dir = getcwd();
for my $command_ref (@command_seq) {
    my ($working_dir,$command) = @{$command_ref};
    
    chdir($working_dir);
    
    if ($opt{debug}) {
       print("$command -cfg $config_file\n");
    } else {
       print "Working on $command\n";
       system("./$command -cfg $config_file\n");
       print "Done with $command\n\n";
    }

    chdir($starting_dir);
}

open OUTPUT, ">../../results/$cfg{exp_id}/Exp Note.txt" or die "$!";
print OUTPUT $cfg{self_note};
close OUTPUT;

if (not $opt{debug}) {
    chdir('../../results');
    system("zip -r $cfg{exp_id}.zip $cfg{exp_id}");
} else {
    print("zip -r $cfg{exp_id}.zip $cfg{exp_id}");
}

&send_processing_finished_email;

&remove_PID_file;

################################################################################
# Functions
################################################################################

###########################################################
# Process Control
###########################################################

sub processing_occuring {
    if (-e $PID_file) {
        open PID_DATA, $PID_file;
        my @PID_num = <PID_DATA>;
        close PID_DATA;

        my @ps_output = `ps $PID_num[0]`;
        
        if (scalar(@ps_output) > 1) {
            return 1;
        } else {
            print "Found a run.pid file, but no corresponding program, removing PID file\n";
            unlink $PID_file;
            return 0;
        }
    } else {
        return 0;
    }
}

sub make_processing_PID_file {
    open PID_DATA, ">$PID_file" or die $!;
    print PID_DATA $$;
    close PID_DATA;
}

sub remove_PID_file {
    unlink $PID_file;
}

###########################################################
# Exp Setup
###########################################################

sub move_uploaded_data {
    my $cfg_file = shift;
    
    my $data_file;
    if ($cfg_file =~ /(.*).cfg/) {
        $data_file = $1;
    } else {
        die "Couldn't find data file to correspond to: $cfg_file";
    }

    my $exp_name = basename($data_file);
    my $exp_data_dir = "../../data/$exp_name";
    
    print "Making dir $exp_data_dir\n";
    if (not -e $exp_data_dir) {
        mkdir $exp_data_dir or die $!;
    }

    print "Moving $data_file to $exp_data_dir\n";
    move($data_file,$exp_data_dir);
    move($cfg_file,$exp_data_dir);
    $cfg_file = "$exp_data_dir/$exp_name.cfg";
    system("echo '<<include ../config/FAAI_default.cfg>>' >> $cfg_file");

    $data_file = "$exp_data_dir/$exp_name";
    print "Changing data file to: $data_file\n";
    my $unzip_return = system("unzip -o -q -d $exp_data_dir $data_file");
    if ($unzip_return) {
        die "unzip returned $unzip_return, failed";
    }

    my $image_folder;
    for (<$exp_data_dir/*>) {
        if (-d $_) {
            if (defined $image_folder) {
                die "Found multiple image folders on unzip.";
            }
            if ($_ =~ /$exp_data_dir(.*)/) {
                $image_folder = $1;
            } else {
                die "Found image folder, but couldn't find exp_data_dir: $_"
            }
        }
    }
    system("echo 'adhesion_image_folder = $image_folder' >> $cfg_file");

    return $cfg_file;
}

###########################################################
# Mail
###########################################################

sub send_processing_started_email {
    my $subject = "Started processing your job: $cfg{exp_id}";
    my $body = "Your focal adhesion analysis job has started processing.\n\nThe " . 
    "note you left yourself with this expeirment is: $cfg{self_note}.";
   
    &send_email($subject, $body);
}

sub send_processing_finished_email {
    my $subject = "Finished processing your job: $cfg{exp_id}";
    my $body = "Your focal adhesion analysis job has finished processing.\n\nThe " . 
    "note you left yourself with this expeirment is: $cfg{self_note}.\n\nYou can " .
    "download your results here: " .
    "\n\nhttp://balder.bme.unc.edu/~mbergins/FAAI_webapp/results/$cfg{exp_id}.zip";
   
    &send_email($subject, $body);
}

sub send_email {
    my $subject = $_[0];
    my $body = $_[1];
    
    if ($opt{debug}) {
        print("echo -e \"$body\" | mail -s \"$subject\" $cfg{email}");
    } else {
        system("echo -e \"$body\" | mail -s \"$subject\" $cfg{email}");
    }
}
