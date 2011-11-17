#!/usr/bin/perl -w

package webserver_funcs;

use strict;

use base 'Exporter';

our @EXPORT = qw/print_html_end/;

sub print_html_end {
	my $q = shift;
	print $q->br,$q->br,$q->hr, $q->p, 
		"Return to the ", $q->a({href=>"/FA_webapp"},'home page');
	print "</div>\n";
	print $q->end_html;
}

1;
