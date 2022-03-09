#!/usr/bin/perl
#
use warnings;
use strict;
#
my $file = $ARGV[0] or die "Need to get a tab-delimited file on the command line\n";
open my $fh, '<', $file;
my $base_counter = 1;
my @line;
my @repaired;
my @temp;
while (<$fh>) {
 
	chomp($_);
	push (@line, $_);
}	
my @beginning = split ("\t", $line[1]);
$base_counter = $beginning[0];

while ($#line >= 0) {
	@temp = split ("\t", $line[0]);
	if ($temp[0] eq 'bp') {
		@repaired = join ("\t", @temp);
		shift @line;
		print "@repaired\n";
	}
	elsif ($temp[0] == $base_counter) {
		@repaired = join ("\t", @temp);
		$base_counter++;
		shift @line;
		print "@repaired\n";	
	} else {
		print "$base_counter\t0\t0\t0\t0\t0\t0\t0\t0\n";
		$base_counter++; 
	}

}

close $fh;
#
#
