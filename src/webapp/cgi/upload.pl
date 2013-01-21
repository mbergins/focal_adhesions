#!/usr/bin/perl -wT

use strict;
use File::Path qw(make_path remove_tree);
use File::Basename;
use File::Spec::Functions;
use File::Copy;
use File::Temp qw(tempfile);
use POSIX;
use CGI::Pretty;
use CGI::Carp "fatalsToBrowser";
use IO::Handle;
use Config::General;
use HTML::Template;

my $start = time;

$| = 1;

###############################################################################
# Main Program
###############################################################################

# my $q = CGI->new(\&hook);
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
	
	while ($output_file =~ /FAAS_.*_.*/) {
		unlink $output_file;
		($output_handle, $output_file) = tempfile('FAAS_XXXXXX',DIR=>catdir('upload')) or die "$!";
	}

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
    $new_cfg{time_spacing} = $q->param('time_spacing');
    $new_cfg{min_independent_size} = 14;
    $new_cfg{upload_time} = $end-$start;
    $new_cfg{exp_ID} = basename($output_file);
    $new_cfg{submitter_ip} = $q->remote_addr();

	$ENV{PATH} = '/bin/';
    $new_cfg{sub_date} = `/bin/date`;

	chomp $new_cfg{sub_date};
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
	
	##############################################################################
    # Zipped Image Set Setup
	##############################################################################
	my $file_type = &determine_file_type($output_file);
	if ($file_type eq 'TIFF') {
		&move_tiff_file($output_file);
	} else {
		unlink $output_file, "$output_file.cfg";
		
		my $template = HTML::Template->new(filename => 'template/upload_problem.tmpl');
		print $q->header();
		print $template->output;
	}
	
	#File type unknown was taken care of above, now we need to return the exp
	#status page
	if ($file_type ne 'unknown') {
		my $old_exps = $q->cookie('exp_ids');
		my $template = HTML::Template->new(filename => 'template/upload_success.tmpl');
		$template->param(upload_time => $end - $start,
						 status_link => "exp_status.pl?exp_id=$new_cfg{exp_ID}",
					 	 );
		
		if ($old_exps eq '') {
			$old_exps = "$new_cfg{exp_ID}";
		} else {
			$old_exps .= ",$new_cfg{exp_ID}";
		}

		my $cookie = $q->cookie(-name=>"exp_ids",-value=>$old_exps,expires=>"+99y");
		print $q->header(-cookie=>$cookie);
		print $template->output;
	}
} else {
	my $template = HTML::Template->new(filename => 'template/upload_start.tmpl');
	if ($q->param('debug_options')) {
		$template->param(debug_options => 1);
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
	} else {
		$file_type = 'unknown';
	}
	return $file_type;
}

sub move_tiff_file {
	my $file = shift;
	
	my $id = basename($file);

	make_path("upload/temp_$id/Images/FA_marker");
	move($file,"upload/temp_$id/Images/FA_marker/data.tif") or die $!;
	
	chdir "upload/temp_$id";
	$ENV{PATH} = '/usr/bin/';
	system("/usr/bin/zip -r -0 $id.zip Images");
	chmod 0777, "$id.zip" or die "$!";
	move("$id.zip","..") or die $!;
	chdir '../../';
	
	remove_tree("upload/temp_$id");
}

sub hook {
        my ($filename, $buffer, $bytes_read, $data) = @_;
		sleep 1;
}
