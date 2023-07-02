#!/usr/bin/env perl
use strict;
use warnings;
use utf8;

if ( $#ARGV < 2 ){
    print "Usage: ./bench.pl <numcores> <seq|tp|ws> <vec|novec>\n";
    exit;
}

my $numCores = $ARGV[0];
my $type =  0;
my $N = 100;
my $ver = '';

if ( $ARGV[1] eq 'seq' ){
    $type = 0;
} elsif (  $ARGV[1] eq 'tp'  ){
    $type = 1;
} elsif (  $ARGV[1] eq 'ws'  ){
    $type = 2;
}

if ( $ARGV[2] eq 'vec' ){
    $ver = 'vec';
} elsif (   $ARGV[2] eq 'novec' ){
    $ver = 'novec';
} else {
    print "Unvalid version. <vec|novec> are only valid ones.";
    exit;
}

print("# micro $numCores $ARGV[1]\n");

while ( $N < 8000000 ) {
    my $result;
    my $performance = '0.00';

    while ( $performance eq '0.00' ){
        $result =  `./micro-$ver $type $N`;
        $result =~ /([0-9.]+) ([0-9.]+)/;
        $performance = $2;
    }

    print $result;
    $N = int($N * 1.2);
}
