package login;
use Dancer ':syntax';
use Dancer::Session::YAML;
use strict;
use warnings;
use Cwd;
use Sys::Hostname;
use Crypt::SaltedHash;
use Storable qw(lock_store lock_nstore lock_retrieve);
use Data::Dumper;

my $start_time = time;
my $login_data_file = '../user_login.stor';

###############################################################################
# Main
###############################################################################

get '/login' => sub {
	template 'login';
};

post '/login' => sub {
	my $user_id = param('email');
	
	my %user_info;
	if (-w $login_data_file) {
		%user_info = %{lock_retrieve($login_data_file)};
	}	
	
	#deal with the case that the user_id exists, IE an account already exists
	#for this user
	if (defined $user_info{$user_id}) {
		my $csh = Crypt::SaltedHash->new(salt=>$user_info{$user_id}{salt});
		$csh->add(param('password'));
		my $salted = $csh->generate;
		
		if ($salted eq $user_info{$user_id}{password}) {
			session user_id => $user_id;
			template 'login', {login_success => 1, user_id => $user_id};
		} else {
			template 'login', {bad_login => 1, user_id => $user_id};
		}
	} else {
		my $csh = Crypt::SaltedHash->new();
		$csh->add(param('password'));
		my $salted = $csh->generate;

		$user_info{$user_id}{password} = $salted;
		$user_info{$user_id}{salt} = $csh->salt_bin();

		lock_store \%user_info, $login_data_file;
		session user_id => $user_id;
		template 'login', {login_success => 1, user_id => $user_id, new_user=>1};
	}
};

###############################################################################
# Functions
###############################################################################

true;
