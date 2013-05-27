#!/usr/bin/perl -w

use Storable qw(lock_retrieve);
use Data::Dumper;

%user_info = %{lock_retrieve('../user_login.stor')};
print Dumper(\%user_info);

%user_info = %{lock_retrieve('../user_exp_info.stor')};
print Dumper(\%user_info);
