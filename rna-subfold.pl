#!/usr/bin/env perl
#===============================================================================
#
#         FILE:  rna-sub-fold.pl
#
#        USAGE:  ./rna-sub-fold.pl
#
#  DESCRIPTION:
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  YOUR NAME (),
#      COMPANY:
#      VERSION:  1.0
#      CREATED:  12/10/2014 06:08:48 PM
#     REVISION:  ---
#===============================================================================

use 5.10.0;
use utf8;
use strict;
use warnings;

use Getopt::Long;
use IPC::Open2;

our $COMMON = "--noLP";

our $RNA_FOLD   = "RNAfold -p -d2 $COMMON";
# our $RNA_SUBOPT = "RNAsubopt -s -C $COMMON";
our $RNA_SUBOPT = "RNAsubopt -s -C $COMMON";
our $kT         = 0.616321;

my ( $RADIUS, $VERBOSE, $CSV ) = ( 100, 0, 0 );

#-----------------------------------------------------------------------------
$SIG{CHLD} = 'IGNORE';

GetOptions(
	"v" => \$VERBOSE,
	"c" => \$CSV,
	"r=i" => \$RADIUS,
);

my $FOLD_CMD   = "$RNA_FOLD";
my $SUBOPT_CMD = "$RNA_SUBOPT -e $RADIUS";

{
	my ( $name, $sequence, $constraint );
	while ( my $line = <STDIN> ) {

		chomp $line;

		if ( $line =~ m/^>\s*(.*)$/ ) {
			$name       = $1;
			$sequence   = undef;
			$constraint = undef;
			print $line . "\n" unless $CSV;
		}

		elsif ( $line =~ m/^\s*[A-Za-z]+\s*$/ ) {
			$sequence = $line;
			print $line . "\n" unless $CSV;
		}

		elsif ( $line =~ m/^\s*[.<>x|()]+\s*$/ ) {
			$constraint = $line;
			print $line . "\n" unless $CSV;
			calculate_opt( $name, $sequence, $constraint );
		}
	}
}

sub calculate_opt
{
	my ( $name, $sequence, $constraint ) = @_;
	my ( $structure, $energy, $min_energy );

	my ( $fold_out, $fold_in );

	open2( $fold_out, $fold_in, $FOLD_CMD );

	print $fold_in $sequence . "\n";

	close $fold_in;

	while ( defined( my $line = <$fold_out> ) ) {

		if ( $line =~ m/\[([+-]?[0-9.]+)\]/ ) {
			$energy = 1 * $1;

			print $line if $VERBOSE;

			calculate_sub_opt( $name, $sequence, $constraint, $energy );

			last;
		}
	}

	close $fold_out;
}

sub calculate_sub_opt
{
	my ( $name, $sequence, $constraint, $energy ) = @_;
	my ( $subopt_out, $subopt_in );

	open2( $subopt_out, $subopt_in, $SUBOPT_CMD );

	print $subopt_in $sequence . "\n";
	print $subopt_in $constraint . "\n";

	close $subopt_in;

	my $sum = 0;
	while ( defined( my $line = <$subopt_out> ) ) {

		chomp $line;

		if ( $line =~ m/^([^ ]+) ([+-]?[0-9.]+)$/ ) {
			my ( $s_structure, $s_energy ) = ( $1, $2 );
			my $percent = 100 * exp( ( $energy - $s_energy ) / $kT );

			$sum += $percent;

			printf "%s % 6.3f % 4.3f%%\n", $s_structure, $s_energy, $percent
				if $VERBOSE;
		} else {
			#print "W: " . $line . "\n";
		}
	}

	if($CSV) {
		printf "%s;%.10f\n", $name, $sum;
	} else {
		printf "percent matching constraint: %6.10f\n", $sum;
	}

	close $subopt_out;
}
