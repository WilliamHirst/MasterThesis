#!/usr/bin/perl -w
#
# Configure script to configure DarkSUSY.
# NOTE. Run this script before compiling DarkSUSY. This is done automatically
# either at the configure stage or at the make stage.
#
# Author: Joakim Edsjo, edsjo@physto.se
# Date: April 12, 2001
# Modified: 
#   20080107 Paolo Gondolo split lines in dsdir.h
#   20080226 Joakim Edsjo, rewrite to use parameters instead


$deffile="include/dsdirver.h";

$dsroot = shift;
$dsroot =~ s#/$##;   # take away final / if any
$dsversion = $dsroot;
$dsversion =~ s#^.*/##;

$dsinstall=shift;
$dsinstall =~ s#/$##;   # take away final / if any
$dsinstall .= "/"; # put it back so that we know we have it.


# Add revision number to dsversion
$rev="";
if (open(IN,"svn info|")) {
    while(defined($line=<IN>)) {
	if ($line =~ /^Revision:\s+(\d+)$/) {
	    $rev=" (rev $1)";
	}
    }
    close(IN);
} else {
    print "Couldn't invoke svn. No revision number added to version tag.\n";
}
$dsversion .= $rev;

print "DarkSUSY version: $dsversion\n";


# OK, Now dsversion and dsintall contains the information needed. Now
# make the include file.

open(OUT,">$deffile") || die "Can't open $deffile for writing.\n";
print OUT <<END
*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                           dsdirver.h                             ***
***         this piece of code is needed as a separate file          ***
***            the rest of the code 'includes' dsdirver.h            ***
c----------------------------------------------------------------------c
END
;

$date=`date`;
print OUT "*** This file is created by config2.pl on $date\n";

$len = length ($dsversion);
print OUT "      character*${len} dsver\n";
$line=" " x 6 . "parameter(dsver='${dsversion}')\n";
$line=contline_string($line);
print OUT $line;

$len = length ($dsinstall);
print OUT "      character*${len} dsinstall\n";
$line=" " x 6 . "parameter(dsinstall='${dsinstall}')\n";
$line=contline_string($line);
print OUT $line;

$len2=int(length($dsversion)/10)*10+50;
print OUT " " x 6 . "character*${len2} dsversion\n";
print OUT " " x 6 . "common /dsv/dsversion\n";
print OUT " " x 6 . "save /dsv/\n";

print OUT "***" . " " x 66 . "***\n";
print OUT "*" x 26 . " end of dsdirver.h " . "*" x 27 . "\n";

### Split long lines

sub contline {
    my $line = $_[0];
    my $out;
    my $i;
 
    $out="";
    if (length($line) >= 71) {
        print "*** A line is too long. I will split it.\n";
        print "Line before: \n$line";
        while (length($line) >= 71) {
            for ($i=70; $i==20; $i--) {
                if (substr($line,$i,1) =~ m@(\*|\+|-|/)@ ) {
                    last;
                }
            }
            $out=$out . substr($line,0,$i-1) . "\n";
            $line = " " x 5 . "&  " . substr($line,$i-1)
        }
        $out = $out . $line;
        print "Line after: \n$out";
    } else {
        $out=$line;
    }
 
    return $out;
}


# The following only works for strings

sub contline_string {
    my $line = $_[0];
    my $out;
    my $max;
 
    $max=70;
    $out="";
    while (length($line) > $max) {
#        print "*** A line is too long. I will split it.\n";
#        print "Line before: \n$line";
	$max -=4 if (length($line) <$max+4); # to avoid cutting at wrong place
	$out .= substr($line,0,$max) . "'\n";
	$line= "     & //'" . substr($line,$max);
#        print "Line after: \n$out";
    }
    $out .= $line;
 
    return $out;
}
