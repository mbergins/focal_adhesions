#!/usr/bin/perl

################################################################################
# Global Variables and Modules
################################################################################
use lib "../lib";
use lib "../lib/perl";

use strict;
use threads;
use threads::shared;
use File::Spec::Functions;
use File::Basename;
use File::Find::Rule;
use Benchmark;
use Getopt::Long;
use Cwd;
use Data::Dumper;
use POSIX;

use Config::Adhesions qw(ParseConfig);

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l", "exp_filter=s", "no_email") or die;

if (-e '/opt/lsf/bin/bjobs' && not $opt{lsf}) {
	die "LSF appears to be installed on this machine, don't you want to use it?" 
}	

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

my $t1 = new Benchmark;
$|  = 1;

################################################################################
# Main
################################################################################

#This data structure acts as the overall control mechanism for which programs
#are run to perform the analysis. The first layer of the array are another set
#of arrays with sets of commands that can be executed simultaneously. The next
#layer holds all of those commands with the appropriate directory to execute the
#commands in.
my @overall_command_seq = (
	[ [ "../find_cell_features",      "./collect_mask_image_set.pl" ], ],
	[ [ "../find_cell_features",      "./run_matlab_over_field.pl -script find_exp_thresholds"], ],
	[ [ "../find_cell_features",      "./collect_fa_image_set.pl" ], ],
	[ [ "../find_cell_features",      "./collect_fa_properties.pl" ], ],
	[ [ "../analyze_cell_features",   "./build_tracking_data.pl" ], ],
	[ [ "../analyze_cell_features",   "./track_adhesions.pl" ], ],
	[ [ "../analyze_cell_features",   "./gather_tracking_results.pl" ], ],
	[ [ "../analyze_cell_features",   "./build_R_models.pl" ], ],
	[ [ "../analyze_cell_features",   "./find_signif_models.pl" ], ],
	[ [ "../analyze_cell_features",   "./build_alignment_models.pl" ], ],
	[ [ "../visualize_cell_features", "./collect_visualizations.pl" ], ],
	[ [ "../find_cell_features",      "./run_matlab_over_field.pl -script ../visualize_cell_features/max_intent_project" ], ],
);

if (defined $cfg{photo_bleach_correction} && $cfg{photo_bleach_correction}) {
	@overall_command_seq = (@overall_command_seq[0..2], 
		[ [ "../find_cell_features", "./run_matlab_over_field.pl -script apply_adhesion_bleaching_correction.m" ], ],
		@overall_command_seq[3..$#overall_command_seq]);
}

#some of the scripts only need to be run once for each experiment, this will
#rely on being able to find an experiment with one of the following values in
#it's name
my @run_only_once = qw(find_signif_models);

my @skip_check = qw(find_exp_thresholds track_adhesions gather_tracking_results
	build_R_models build_alignment_models collect_visualizations
	max_intent_project make_eccen_filtered_vis);

my $cfg_suffix = basename($opt{cfg});
$cfg_suffix =~ s/.*\.(.*)/$1/;

#config file processing
my @config_files = File::Find::Rule->file()->name( "*.$cfg_suffix" )->in( ($cfg{data_folder}) );
@config_files = sort @config_files;
if (exists($opt{exp_filter})) {
   @config_files = grep $_ =~ /$opt{exp_filter}/, @config_files;
}
die "No config files left after filtering." if (scalar(@config_files) == 0);

my @run_once_configs = $config_files[0];

#######################################
# Program Running
#######################################
our @running_jobs;

if (not($opt{debug}) && $opt{lsf}) {
	my $job_queue_thread = threads->create('shift_idle_jobs','');
	$job_queue_thread->detach;
}

my $starting_dir = getcwd;
for (@overall_command_seq) {
	my $command_start_bench = new Benchmark;
	my @command_seq = @{$_};
	@command_seq = map { [ $_->[0], $_->[1] . " -lsf" ] } @command_seq if $opt{lsf};
	print "Starting on $command_seq[0][1]\n";
	if (grep $command_seq[0][1] =~ /$_/, @run_only_once) {
		&execute_command_seq(\@command_seq, $starting_dir,\@run_once_configs);
	} else {
		&execute_command_seq(\@command_seq, $starting_dir);
	}

	#If debugging is not on, we want to wait till the current jobs finish
	#and then check the file complements of the experiments for completeness
	&wait_till_LSF_jobs_finish if ($opt{lsf} && not($opt{debug}));
	if (not($opt{debug}) && not(grep $command_seq[0][1] =~ /$_/, @skip_check)) {
		print "Checking for all output files on command $command_seq[0][1]\n";
		my %exp_sets = &check_file_sets(\@config_files);

		for (1..2) {
			#if no experiments are left to retry, we break out and continue
			#to the next command set
			if (not(@{$exp_sets{retry}})) {
				last;
				next;
			}
			print "Retrying these experiments:\n".
				  join("\n", @{$exp_sets{retry}}) . "\n";
			&execute_command_seq(\@command_seq, $starting_dir, \@{$exp_sets{retry}});
			&wait_till_LSF_jobs_finish if ($opt{lsf});
			my %these_exp_sets = &check_file_sets(\@{$exp_sets{retry}});

			push @{$exp_sets{good}}, @{$these_exp_sets{good}};
			@{$exp_sets{retry}} = @{$these_exp_sets{retry}};
		}

		#check if there are any files left in the retry set, if so, clear
		#out the failed experiments from the next round
		if (@{$exp_sets{retry}}) {
			print "\nProblem with collecting full file complement on experiments after three tries:\n\t" .
				join("\n\t", @{$exp_sets{retry}}) . "\nRemoving them from the next run.\n\n";
			@config_files = @{$exp_sets{good}};  
		} else {
			print "Output file set complete, moving on.\n";
		}
	}

	if (not $opt{debug}) {
		my $command_end_bench = new Benchmark;
		my $td = timediff($command_end_bench, $command_start_bench);
		print "The command took:",timestr($td),"\n";
	}
	
	#we are done with the current command, so we print a few spacer lines to
	#indicate we have moved onto a new command
	print "\n\n";
	
	#we expect the running jobs listed to be emptied after every command
	@running_jobs = ();
}

if (not($opt{debug})) {
	my $t_bsub = new Benchmark;
	my $td = timediff($t_bsub, $t1);
	
	my $time_diff_str = "\"Took:" . timestr($td) . "\"";

	my $command = "bsub -J \"Job Finished: $opt{cfg}\" echo $time_diff_str";
	
	if (not $opt{no_email}) {
		system($command) if $opt{lsf};
	}
}

my $t2 = new Benchmark;
my $td = timediff($t2, $t1);
print "\nThe pipeline took:",timestr($td),"\n\n";

################################################################################
# Functions
################################################################################

sub execute_command_seq {
    my @command_seq  = @{ $_[0] };
    my $starting_dir = $_[1];
    my @these_config_files = @config_files;
    if (scalar(@_) > 2) {
        @these_config_files = @{$_[2]};
    }
	
    foreach my $set (@command_seq) {
        my $dir     = $set->[0];
        my $command = $set->[1];
        my $executed_scripts_count = 0;
        foreach my $cfg_file (@these_config_files) {
            my $config_command = "$command -cfg $cfg_file";
            chdir $dir;
            my $return_code = 0;
            if ($opt{debug}) {
				# print "Working in directory: $dir\n";
                print $config_command, "\n";
            } else {
                $executed_scripts_count++;
                print "RUNNING: $config_command $executed_scripts_count/" . scalar(@these_config_files) . "\n";
                $return_code = system($config_command);
				print "RETURN CODE: $return_code\n";
            }
            chdir $starting_dir;

            #if the return code was any number beside zero, indicating a problem
            #with the program exit, remove that config file from the run and
            #continue
            if ($return_code) {
				print "REMOVING: $cfg_file\n";
                @config_files = grep $cfg_file ne $_, @config_files;
            }
        }
    }
}

#######################################
# LSF
#######################################

sub wait_till_LSF_jobs_finish {
    #After each step of the pipeline, we want to wait till all the individual
    #jobs are completed, which will be checked three times
	my $total_checks = 3;
    for (1 .. $total_checks) {
        print "LSF finished check number $_/$total_checks\n";
        my $sleep_time = 10;
        do {
            sleep($sleep_time);
        } while (&running_LSF_jobs);
    }
    print "\n";
}

sub running_LSF_jobs {
    my $bjobs_command = "bjobs";
    if (defined $cfg{job_group}) {
        $bjobs_command .= " -g $cfg{job_group} 2>/dev/null";
    }

    my @lines = `$bjobs_command`;
    
    #find and add the job ids onto the running job list
    my @job_ids;
    for (@lines) {
        if (/^(\d+) /) {
            push @job_ids,$1;
        }
    }
    @job_ids = sort @job_ids;
    
    push @running_jobs, \@job_ids;
    &kill_long_running_jobs();

    if (scalar(@lines) <= 1) {
        return 0;
    } else {
        return 1;
    }
}

sub kill_long_running_jobs {
	#we only need to worry about killing jobs that have run for at least 6
	#hours, we check for running jobs every 10 seconds, so...
    if (scalar(@running_jobs) > (6*60*10)) {
        my @long_running_jobs = @{$running_jobs[$#running_jobs]};
        my @last_hour = @running_jobs[($#running_jobs - 59) .. $#running_jobs];
        foreach my $jobs_ref (@last_hour) {
            #check to make sure the length of the jobs ref is the same as the
            #long running jobs set, if it isn't, break out of the function
            if (scalar(@$jobs_ref) != scalar(@long_running_jobs)) {
                return 0;
            }
            
            #scan through this entry in the jobs list, if an entry is present
            #the current list of jobs, but not in the long_running_jobs list,
            #break away from the function
            for my $this_job (@$jobs_ref) {
                if (not grep $this_job == $_, @long_running_jobs) {
                    return 0;
                }
            }
        }

        #if we made is this far in the function, the entries in
        #long_running_jobs should have each been running for at least an hour
        #beyond the last other jobs in the queue, so we kill them

        for (@long_running_jobs) {
            print "Killing job $_\n";
            system("bkill $_");
        }
    }
}

sub shift_idle_jobs {
    my $bjobs_command = "bjobs";
    if (defined $cfg{job_group}) {
        $bjobs_command .= " -g $cfg{job_group}";
    }
    
    my @lines = `$bjobs_command 2>/dev/null`;
    my @running_lines = `$bjobs_command -r 2>/dev/null`;
    
    if (scalar(@lines) > 1) {
        #dealing with the strange case where all the idle jobs refuse to be
        #shifted into a running queue position, so we slowly shift all the
        #pending jobs into the week queue, where they will certainly get a
        #running jobs spot
        if (scalar(@lines) > scalar(@running_lines)) {
            &move_job_to_week_queue($lines[-1]);
        }
    }
    sleep(60);
    &shift_idle_jobs();
}

sub move_job_to_week_queue {
    my $line = pop @_;
    if ($line =~ /^(\d+)/) {
       system("bmod -q week $1 > /dev/null 2>/dev/null");
    }
}

sub check_file_sets {
    my @config_files = @{$_[0]};

    print "Out of " . scalar(@config_files) . ", # done:";
    my $number_configs_run = 0;
    
    my %exp_sets;
    #intialize these as the code that checks these matrices in the main program
    #assumes that a matrix will be present, even if it is empty
    @{$exp_sets{good}} = ();
    @{$exp_sets{retry}} = ();
    foreach my $this_config_file (@config_files) {
        my $return_code = system "./check_file_complement.pl -cfg $this_config_file";
        
        #if the return code is anything besides zero, add that config file back
        #to the exp_to_retry list
        if ($return_code) {
            push @{$exp_sets{retry}}, $this_config_file ;
        } else {
            push @{$exp_sets{good}}, $this_config_file ;
        }

        $number_configs_run++;
        if ($number_configs_run % ceil(scalar(@config_files)/20) == 0) {
            print " $number_configs_run";
        }
    }
    print "\n";

    return %exp_sets;
}
