#!/usr/bin/perl -w
#
# Configure script to configure DarkSUSY.
# NOTE. Run this script before compiling DarkSUSY.
#
# Author: Joakim Edsjo, edsjo@physto.se
# Date: April 12, 2001
# Modified: 
#   20080107 Paolo Gondolo split lines in dsdir.h

#$mufile="src/mu/dsmudir.h";
#$hafile="src/ha/dshadir.h";
$dsverfile="src/ini/dsversion.h";
#$ntsdfile="src/nt/dsntsddir.h";
$dsinstallfile="src/ini/dsdir.h";

$dsroot = shift;
$dsroot =~ s#/$##;   # take away final / if any
$dsversion = $dsroot;
$dsversion =~ s#^.*/##;

$dsinstall=shift;
$dsinstall =~ s#/$##;   # take away final / if any


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

# Take care of dsversion file.
if (open(IN,"$dsverfile")) {
    $line="";
    while ($in=<IN>) {
	$line .= $in;
    }
    close(IN);
} else {
	$line="";
}

$newline=" " x 6 . "data dsversion/'$dsversion'/\n";
$newline=contline($newline);

if ($newline ne $line) {
    open(OUT,">$dsverfile") || 
      die "Can't open $dsverfile for writing\n";
    print OUT $newline;
    close(OUT);
    print "$dsverfile updated.\n";
} else {
    print "$dsverfile is up-to-date.\n";
}

# Take care of DarkSUSY root directory file

if(open(IN,"$dsinstallfile")){
    $line="";
    while ($in=<IN>) {
	$line .= $in;
    }
    close(IN);
} else {
    $line="";
}

$newline=" " x 5 . "&'$dsinstall/'/\n";
$newline=contline($newline);

if ($newline ne $line) {
    open(OUT,">$dsinstallfile") || 
      die "Can't open $dsinstallfile for writing\n";
    print OUT " " x 6 . "data dsinstall/\n";
    print OUT $newline;
    close(OUT);
    print "$dsinstallfile updated.\n";
} else {
    print "$dsinstallfile is up-to-date.\n";
}


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
