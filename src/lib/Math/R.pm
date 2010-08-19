#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################

use strict;
use warnings;
use File::Path;
use File::Basename;
use File::Temp ();

package Math::R;

###############################################################################
# Functions
###############################################################################

sub execute_commands {
    my @r_code = @{$_[0]};
    
    my $r_command = "R --slave";
    
    my $tmp = File::Temp->new(); 
    my $r_code_filename = $tmp->filename;
    open R_CODE, ">$r_code_filename";
    print R_CODE join("\n",@r_code);
    close R_CODE;
    
    my $r_out = `$r_command <$r_code_filename`;
}

1;
