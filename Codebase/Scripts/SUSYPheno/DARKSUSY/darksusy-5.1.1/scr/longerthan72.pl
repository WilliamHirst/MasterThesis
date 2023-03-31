#!/usr/bin/perl
#
# Script to search for lines longer than 72 columns in Fortran files

while (defined($file=shift)) {
    $lno=0;
    $printed=0;
    open(IN,"$file") || die "Can't open $file for reading.\n";
    while (defined($line=<IN>)) {
	$lno++;
	$line =~ s/\s+$//; #remove trailing spaces
	$exclpoint=index($line,"!");
        if (length($line) > 72 && substr($line,0,1) eq " " && ($exclpoint<0 || $exclpoint>72) ) {
	    print "File: $file\n" if $printed==0;
	    print "Line $lno: $line\n";
            $printed++;
	}
    }  
}
