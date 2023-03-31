#!/usr/bin/perl -w
#
# This script goes through the files given as arguments and changes
# the file names according to the supplied rule.


while (defined($file=shift)) {
    print "Taking care of file $file ...";
    $newfile = $file;
    $newfile =~ s/dshaIB/dsIB/;
    system ("svn move $file $newfile");
#    rename($file,$newfile);
    print " renamed to $newfile\n";
}
