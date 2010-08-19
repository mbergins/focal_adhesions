#!/usr/bin/perl

&wait_till_LSF_jobs_finish;

################################################################################
# Functions
################################################################################

sub wait_till_LSF_jobs_finish {
    for (1 .. 3) {
        my $sleep_time = 5;
        do {
            sleep($sleep_time);
            $sleep_time++;
        } while (&running_LSF_jobs);
    }
}

sub running_LSF_jobs {
    my @lines = `bjobs`;
    if (scalar(@lines) > 1) {
        return scalar(@lines) - 1;
    } else {
        0;
    }
}
