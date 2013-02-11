package login;
use Dancer ':syntax';
use Dancer::Session::YAML;
use strict;
use warnings;
use Cwd;
use Sys::Hostname;

my $start_time = time;

###############################################################################
# Main
###############################################################################

get '/logout' => sub {
	session->destroy;
	template 'logout';
};

###############################################################################
# Functions
###############################################################################

true;
