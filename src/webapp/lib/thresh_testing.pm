package upload;
use Dancer ':syntax';
use strict;
use warnings;
use Cwd;
use Sys::Hostname;
use File::Temp qw/ tempfile tempdir /;
use File::Basename qw/ dirname /;
use File::Spec::Functions;
use File::Path qw/ make_path /;
use File::Copy qw/ move /;
use File::Find;
use Config::General;

my $out_folder = catdir('threshold_testing');
my $start_upload = time;

###############################################################################
# Main
###############################################################################

get '/thresh_testing' => sub {
	template 'thresh_testing';
};

post '/thresh_testing' => sub {
	if (! -e $out_folder) {
		make_path $out_folder or die $!;
		chmod 0777, $out_folder;
	}

	my $input_file = upload('input_image') or die $!;
	my $out_dir = tempdir(DIR=>$out_folder);
	chmod 0777, $out_dir;
	my $out_file = catfile($out_dir,"input");
	$input_file->copy_to($out_file);
	my $end_upload = time;

	my $start_process = time;
	system("octave --eval \"cd ../misc_code; build_thresholded_image_sets('../public/$out_file')\" > /dev/null 2> /dev/null");
	my $end_process = time;
	
	my %cfg = (upload_time => $end_upload - $start_upload, 
		process_time => $end_process - $start_process);

	my $out_html = &generate_html(%cfg);
	open OUTPUT, ">" . catfile($out_dir,'index.html');
	print OUTPUT $out_html;
	close OUTPUT;
	
	redirect catfile($out_dir,'index.html');
};

###############################################################################
# Functions
###############################################################################

sub generate_html {
	my %cfg = @_;
	return template 'thresh_testing_results', \%cfg;
}

true;
