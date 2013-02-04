package results_understanding; 
use Dancer ':syntax';
use strict;
use warnings;

###############################################################################
# Main
###############################################################################

get '/results_understanding' => sub {
	template 'results_understanding';
};

true;
