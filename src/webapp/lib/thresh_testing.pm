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

	my $start_process = time;
	system("octave --eval \"cd ../misc_code; build_thresholded_image_sets('../public/$out_file')\" > /dev/null 2> /dev/null");
	my $end_process = time;
	
	my $out_html = template 'thresh_testing_results', 
		{process_time => $end_process - $start_process};
	open OUTPUT, ">" . catfile($out_dir,'index.html');
	print OUTPUT $out_html;
	close OUTPUT;
	
	redirect catfile($out_dir,'index.html');
};

###############################################################################
# Functions
###############################################################################

true;
