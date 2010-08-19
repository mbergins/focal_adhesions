#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
{ 
package Text::CSV::Simple::Extra;
use strict;
use warnings;
use IO::File;
use Text::CSV;
use File::Basename;
use File::Path;

our @EXPORT = qw( input_mat_csv output_mat_csv output_hash_csv );
use Exporter;
our @ISA = qw(Exporter);
###############################################################################
# Functions
###############################################################################

sub input_mat_csv {
    my @mat;
    my $file = $_[0];

    my $in_hand = new IO::File $file or die "Unable to create csv file: $file";
    my $csv = Text::CSV->new();

    while (<$in_hand>) {
        if ($csv->parse($_)) {
            my @row = $csv->fields;
            push(@mat, [@row]);            
        }
    }

    $in_hand->close;

    return @mat;
}

sub output_mat_csv {
    my @mat = @{$_[0]};
    my $file = $_[1];

    if (ref($mat[0]) ne "ARRAY") {
        @mat = [@mat];
    }

    if (! -e dirname($file)) {
        mkpath(dirname($file))
    }
    
    my $out_hand = new IO::File ">" . $file or die "Unable to create csv file: $file";
    my $csv = Text::CSV->new();

    for my $i (0 .. $#mat) {
        $csv->print($out_hand, \@{ $mat[$i] });
        print $out_hand "\n";
    }

    $out_hand->close;
}

sub output_hash_csv {
    my %hash = %{$_[0]};
    my $file = $_[1];

    my $out_hand = new IO::File ">" . $file or die "Unable to create csv file: $file";
    my $csv = Text::CSV->new();

    while(my @row = each %hash) {
        $csv->print($out_hand, \@row);
        print $out_hand "\n";
    }

    $out_hand->close;
}

}
1;
