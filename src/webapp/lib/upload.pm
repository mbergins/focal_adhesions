package upload;
use Dancer ':syntax';
use Dancer::Session::YAML;
use strict;
use warnings;
use Cwd;
use Sys::Hostname;
use File::Temp qw/ tempfile /;
use File::Basename qw/ dirname /;
use File::Spec::Functions;
use File::Path qw/ make_path /;
use File::Copy qw/ move /;
use File::Find;
use File::Basename;
use Config::General;
use Data::Dumper;
use Storable qw(lock_store lock_nstore lock_retrieve);
use Time::localtime;
use Archive::Extract;


my $out_folder = catdir('..','uploaded_experiments');
my $user_exp_info_file = '../user_exp_info.stor';

###############################################################################
# Main
###############################################################################

get '/upload' => sub {
	my %opt = (time_spacing => 1, stdev_thresh => 2, static => 0, 
		no_ad_splitting => 0, min_linear_model_length => 10, 
		min_adhesion_size => 2, max_adhesion_size => "", FAAI_min_ratio => 3,
		email => '');
	
	foreach (keys %opt) {
		$opt{$_} = param $_ if (defined param $_);
	}
	template 'upload', \%opt;
};

post '/upload' => sub {
	if (! -e $out_folder) {
		make_path $out_folder or die $!;
		chmod 0777, $out_folder;
	}

	my $fa_file = upload('adhesion_file') or die $!;

	if (! &is_file_TIFF($fa_file->tempname)) {
		template 'upload_problem';
	} else {
		my ($fa_fh, $fa_filename) = tempfile("FAAS_XXXXXX",DIR=>$out_folder);
		while ($fa_filename =~ /FAAS_.*_.*/) {
			unlink $fa_filename;
			($fa_fh, $fa_filename) = tempfile("FAAS_XXXXXX",DIR=>$out_folder);
		}
		$fa_file->copy_to($fa_filename);

		my $out_folder = &organize_uploaded_files($fa_filename);
		if (params->{cell_mask_file}) {
			&add_cell_mask_to_exp($out_folder);
		}

		#######################################################################
		# Config Processing
		#######################################################################
		my %cfg;
		$cfg{submitter_ip} = request->address();
		
		if (defined session('user_id')) {
			$cfg{session_user_id} = session('user_id');
		}

		my $date_str = `date`;
		chomp($date_str);
		$cfg{sub_date} = $date_str;

		my @copy_if_defined = qw(stdev_thresh no_ad_splitting min_adhesion_size
		max_adhesion_size email exp_note min_linear_model_length time_spacing
		static FAAI_min_ratio max_adhesion_count confocal_mode 
        atrous_segmentation structure_element_size atrous_export_level);
		foreach (@copy_if_defined) {
			my $val = param $_;
			if (defined $val && $val ne "") {
				$cfg{$_} = param $_;
				if ($cfg{$_} eq "on") {
					$cfg{$_} = 1;
				}
				if ($cfg{$_} eq "off") {
					$cfg{$_} = 0;
				}
			}
		}

		my $out_cfg = new Config::General(\%cfg);
		
		my $header = "<<include ../config/webapp_default.cfg>>";
		open CFG_OUT, ">" . catfile($out_folder, "analysis.cfg");
		print CFG_OUT "$header\n";
		print CFG_OUT $out_cfg->save_string;
		close CFG_OUT;
		chmod 0777, catfile($out_folder, "analysis.cfg");

		#######################################################################
		# Logged in user processing
		#######################################################################
		if (defined session('user_id')) {
			my %user_exp_data;
			if (-w $user_exp_info_file) {
				%user_exp_data = %{lock_retrieve($user_exp_info_file)};
			}
			push @{$user_exp_data{session('user_id')}}, basename($out_folder);
			lock_store \%user_exp_data, $user_exp_info_file;
		}
		
		#######################################################################
		# Return Page
		#######################################################################
		my $exp_status_url = "/exp_status/" . basename($out_folder);
		my $email = param 'email';
		if ($email =~ /gmail/) {
			template 'upload_success', { exp_status_link => $exp_status_url,
				gmail => 1};
		} else {
			template 'upload_success', { exp_status_link => $exp_status_url };
		}
	}
};

###############################################################################
# Functions
###############################################################################

sub is_file_TIFF {
	my $file = shift @_;

	my $type_output = `file $file`;
	if ($type_output =~ /TIFF/) {
		return 1;
	} else {
		return 0;
	}
}

sub organize_uploaded_files {
	my $fa_file = $_[0];
	
	$fa_file =~ /FAAS_(.*)/;
	my $out_folder = catdir(dirname($fa_file),"FAAS_$1_temp");
	
	my $fa_folder = catdir($out_folder,'Images','FA_marker');
	make_path($fa_folder);
	
	move($fa_file, catfile($fa_folder,'data.tif'));
	
	chmod 0777, $out_folder;
	find(\&change_file_perm, $out_folder);
	
	move($out_folder, catdir(dirname($fa_file),"FAAS_$1"));
	return catdir(dirname($fa_file),"FAAS_$1");
}

sub change_file_perm {
	chmod 0777, $_;
}

sub add_cell_mask_to_exp {
	my $out_folder = $_[0];

	my $cell_mask_folder = catdir($out_folder,'Images','cell_mask');
	make_path($cell_mask_folder);
	
	my $cm_filename = catfile($cell_mask_folder,'data.tif');

	my $cm_file = upload('cell_mask_file') or die $!;
	$cm_file->copy_to($cm_filename);
	
	chmod 0777, $out_folder;
	find(\&change_file_perm, $out_folder);
}

true;
