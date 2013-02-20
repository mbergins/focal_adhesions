package all_exp_status;
use Dancer ':syntax';
use strict;
use warnings;
use Cwd;
use Sys::Hostname;
use File::Spec::Functions;
use Data::Dumper;
use Storable qw(lock_store lock_nstore lock_retrieve);
 
my $upload_dir = catdir('..','uploaded_experiments');
my $results_dir = "results";
my $running_dir = catdir('..','..','..','data');

my $user_exp_info_file = '../user_exp_info.stor';

###############################################################################
# Main
###############################################################################

get '/all_exp_status' => sub {
	my @exp_ids;
	if (defined session('user_id')) {
		if (-r $user_exp_info_file) {
			my %user_exp_data = %{lock_retrieve($user_exp_info_file)};
			@exp_ids = reverse @{$user_exp_data{session('user_id')}};
		}
	}
	
	my @cfg_status;
	
	for my $this_id (@exp_ids) {
		my %temp = (exp_status => &get_exp_status($this_id),
			exp_id => $this_id);
		if ($temp{exp_status} eq "processing") {
			$temp{exp_status} = "Processing";
			$temp{markup} = "warning";
		}
		if ($temp{exp_status} eq "queue") {
			$temp{exp_status} = "In Queue, position ".&get_queue_position($this_id);
			$temp{markup} = "info";
		}
		if ($temp{exp_status} eq "done") {
			my $download_url = "/results/".$this_id.".zip";
			$temp{exp_status} = "Done Processing, download <a href=$download_url>here</a>";
			$temp{markup} = "success";
		}
		if ($temp{exp_status} eq "missing") {
			$temp{exp_status} = "Experiment ID not found in data sets";
			$temp{markup} = "error";
		}

		push @cfg_status, \%temp;
	}
	
	my %template_cfg = (exp_list => \@cfg_status);

	template 'all_exp_status', \%template_cfg;
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
