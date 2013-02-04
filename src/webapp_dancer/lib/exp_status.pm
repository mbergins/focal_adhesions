package exp_status;
use Dancer ':syntax';
use strict;
use warnings;
use Cwd;
use Sys::Hostname;
use File::Spec::Functions;
use Data::Dumper;
 
my $upload_dir = catdir('..','uploaded_experiments');
my $results_dir = "results";
my $running_dir = catdir('..','..','..','data');

###############################################################################
# Main
###############################################################################
get '/exp_status/:exp_id' => sub {
	my %template_cfg;

	if (param('exp_id')) {
		$template_cfg{exp_status} = &get_exp_status(param('exp_id'));
		if ($template_cfg{exp_status} eq "queue") {
			$template_cfg{queue_position} = &get_queue_position(param('exp_id'));
		}
		if ($template_cfg{exp_status} eq "done") {
			$template_cfg{download_url} = "/results/".param('exp_id').".zip";
		}
	} else {
		$template_cfg{no_exp_id} = 1;
	}
	template 'exp_status', \%template_cfg;
};

###############################################################################
# Functions
###############################################################################

sub get_exp_status {
	my $exp_id = $_[0];
	
	my $exp_status;
	if (-e catdir($results_dir,"$exp_id.zip")) {
		$exp_status = "done";
	} else {
		if (-e catdir($upload_dir, $exp_id)) {
			$exp_status = "queue";
		} else {
			if (-e catdir($running_dir, $exp_id)) {
				$exp_status = "processing";
			} else {
				$exp_status = "missing";
			}
		}
	}
}

sub get_queue_position {
	my $exp_id = $_[0];

	my @uploaded_files = <$upload_dir/*>;
	my @sorted_upload = sort { -C $b <=> -C $a } @uploaded_files;
	
	my @position = grep $sorted_upload[$_] eq catdir($upload_dir,$exp_id), 0..$#sorted_upload;

	return $position[0] + 1;
}

true;
