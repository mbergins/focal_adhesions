#!/usr/bin/perl -wT

use strict;
use File::Path qw(make_path remove_tree);
use File::Basename;
use File::Spec::Functions;
use File::Copy;
use File::Temp qw(tempfile);
use POSIX;
use CGI::Pretty;
use CGI::Carp;
use IO::Handle;
use Config::General;
use HTML::Template;

my $start = time;

$| = 1;

###############################################################################
# Main Program
###############################################################################

my $q = CGI->new();

#First check for a valid file handle, if so, present user with acknowledgement
#of the upload
my $lightweight_fh = $q->upload('uploaded_file');
if (defined $lightweight_fh) {
    # Upgrade the handle to one compatible with IO::Handle:
    my $io_handle = $lightweight_fh->handle();
    $q->param('uploaded_file') =~ /(.*)/;
    
    my $temp_output = catdir('upload');
    if (! -e $temp_output) {
        make_path($temp_output);
		chmod 0777, $temp_output;
    }

    my ($output_handle, $output_file) = tempfile('FAAS_XXXXXX',DIR=>catdir('upload')) or die "$!";
	
    my $buffer;
    while (my $bytesread = $io_handle->read($buffer,1024)) {
        print $output_handle $buffer or die;
    }
    close $output_handle;
	my $end = time;
	
	##############################################################################
    # Config File Setup
	##############################################################################
    my %new_cfg;

	#standard required parameters
    $new_cfg{stdev_thresh} = $q->param('filter_thresh');
    $new_cfg{min_independent_size} = 14;
    $new_cfg{upload_time} = $end-$start;
    $new_cfg{exp_ID} = basename($output_file);
    $new_cfg{submitter_ip} = $q->remote_addr();
	if ($q->param('self_note') ne '') {
			$new_cfg{self_note} = $q->param('self_note');
	}
	if (defined $q->param('phone')) {
		my $phone = $q->param('phone');
		$phone =~ s/-//g;
		if ($phone ne 'XXXXXXXXXX') {
			$new_cfg{phone} = $phone;
			$new_cfg{provider} = $q->param('provider');
		}
	}
	if ($q->param('email') ne '') {
		$new_cfg{email} = $q->param('email');
	}
	if ($q->param('min_adhesion_size') ne '') {
		$new_cfg{min_adhesion_size} = $q->param('min_adhesion_size');
	}
	if ($q->param('max_adhesion_size') ne '') {
		$new_cfg{max_adhesion_size} = $q->param('max_adhesion_size');
	}
	if (defined $q->param('short_exp')) {
		$new_cfg{short_exp} = 1;
	}
	if (defined $q->param('no_ad_splitting')) {
		$new_cfg{no_ad_splitting} = 1;
	}
	&output_config(\%new_cfg,"$output_file.cfg");
	
	my $file_type = &determine_file_type($output_file);
	if ($file_type eq 'TIFF') {
		&move_tiff_file($output_file);
	} elsif ($file_type eq 'zip') {
		move($output_file, $output_file . ".zip");
		chmod 0777, "$output_file.zip" or die "$!";
	} else {
		my $file_out = `/usr/bin/file $output_file`;
		unlink $output_file, "$output_file.cfg";
		
		my $template = HTML::Template->new(filename => 'template/upload_problem.tmpl');
		$template->param(file_type => $file_out);
		print $q->header();
		print $template->output;
	}
	
	if ($file_type ne 'unknown') {
		my $template = HTML::Template->new(filename => 'template/upload_success.tmpl');
		$template->param(upload_time => $end - $start,
						 status_link => "exp_status.pl?exp_id=$new_cfg{exp_ID}");
		print $q->header();
		print $template->output;
	}
} else {
	my $template = HTML::Template->new(filename => 'template/upload_start.tmpl');
	if ($q->param('advanced_options')) {
		$template->param(advanced_options => 1);
	}
	print $q->header();
	print $template->output;
}

###############################################################################
# Functions
###############################################################################

sub output_config {
	my %cfg = %{$_[0]};
	my $target_file = $_[1];
	
	open OUTPUT, ">$target_file" or die $!;
	print OUTPUT "<<include ../config/webapp_default.cfg>>\n\n";
	foreach (keys %cfg) {
		print OUTPUT "$_ = $cfg{$_}\n";
	}
	close OUTPUT;

    chmod 0777, "$target_file" or die "$!";
}

sub determine_file_type {
	my $file_name = shift;
	$ENV{PATH} = '/usr/bin/';
	my $file_out = `/usr/bin/file $file_name`;

	my $file_type;
	if ($file_out =~ /TIFF/) {
		$file_type = 'TIFF';
	} elsif ($file_out =~ /Zip/) {
		$file_type = 'zip';
	} else {
		$file_type = 'unknown';
	}
	return $file_type;
}

sub move_tiff_file {
	my $file = shift;

	$file = basename($file);

	chdir 'upload';
	# move($file,"tmp_$file");
	make_path("Images/FA_marker");

	move($file,"Images/FA_marker/data.tif");
	$ENV{PATH} = '/usr/bin/';
	system("/usr/bin/zip -r -q -0 $file.zip Images");
	remove_tree("Images");
	
	chmod 0777, "$file.zip" or die "$!";
	chdir '..';
}
