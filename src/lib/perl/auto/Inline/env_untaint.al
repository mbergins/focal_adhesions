# NOTE: Derived from blib/lib/Inline.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Inline;

#line 973 "blib/lib/Inline.pm (autosplit into blib/lib/auto/Inline/env_untaint.al)"
#==============================================================================
# Blindly untaint tainted fields in Inline object.
#==============================================================================
sub env_untaint {
    my $o = shift;

    for (keys %ENV) {
	($ENV{$_}) = $ENV{$_} =~ /(.*)/;
    }
    my $delim = $^O eq 'MSWin32' ? ';' : ':';
    $ENV{PATH} = join $delim, grep {not /^\./ and
				      not ((stat($_))[2] & 0022)
				  } split $delim, $ENV{PATH};
    map {($_) = /(.*)/} @INC;
}

# end of Inline::env_untaint
1;
