#!/usr/bin/perl -w

use Getopt::Long;
use File::Find;

my %opt;
$opt{debug} = 0;
$opt{resample} = 300;
$opt{dir} = '.';
$opt{convert_to} = "png";
GetOptions(\%opt, "debug|d", "resample=s", "dir=s", "convert_to=s");

find(\&find_vect_files, $opt{dir});

sub find_vect_files {
    if (   $_ =~ /.*\.svg/
        && not($File::Find::name =~ /no_conv/)
        && not($File::Find::name =~ /\.svn/)) {
       &convert_vect_images($_);
    }
}

sub convert_vect_images {
    my $image_name = $_[0];
    my $output_name = $image_name;
    $output_name =~ s/(.*)\..*/$1.png/;
    
    my $command = "time inkscape -z $image_name -d $opt{resample} --export-png=$output_name --export-background-opacity=1.0";
    
	if ($opt{debug}) {
        print $command, "\n";
		if (defined $opt{convert_to}) {
			my $converted_file = $output_name;
			$converted_file =~ s/\.png/\.$opt{convert_to}/;
			print "convert $output_name $converted_file\n";
			print "unlink $output_name\n";
		}
    } else {
        system($command);
		if (defined $opt{convert_to}) {
			my $converted_file = $output_name;
			$converted_file =~ s/\.png/\.$opt{convert_to}/;
			system "convert $output_name $converted_file";
			# system "unlink $output_name";
		}
    }
}
