#!/usr/bin/perl -w
#
# Script to go through all the subdirectories in src/ and create
# makefiles that includes all *.f files and properly defines all
# included files as dependencies.
#
# Author: Joakim Edsjo, edsjo@physto.se
# Date: August 31, 2000.
# This script is by now obsolete, use makemf.pl instead that compiles
# one file at a time. 

if (@ARGV >0) {
    @dirlist=@ARGV;
} else {
    @dirlist=qw(ac an an1l anstu 
      as bsg db dd ep ep2 ge ha hm hr ini mu nt
      pb rd rge rn su xcern xcmlib); 
}
$date=`date +"%b %d, %Y"`;
chomp($date);

print "Going through subdirectories and creating makefile.in's...\n";

chdir("src");
foreach $dir (@dirlist) {
    print "Now taking care of src/$dir...";
    chdir($dir) || die "Can't cd to $dir\n";
    @tmpfiles=<*.f *.F>;
    @files=();
    $suf=".f";
    foreach $file (@tmpfiles) {
        push(@files,$file) unless ($file eq "ds$dir.f" or $file eq "ds$dir.F");
	$suf=".F" if $file =~ /\.F$/;
    }

    @deps=getdeps(@files);
    rename("makefile.in","makefile.in.old");
    open(FILE,">makefile.in") || die "Can't open makefile.in in src/$dir.\n";
    print_header();
    $line="INC_DEP = " . join(" ",@deps);
    print_line($line);
    print FILE "vpath %.h \$(DINC)\n\n";
    $line="$dir = " . join(" ",@files);
    print_line($line);
    print FILE "all : ds$dir.o\n\n";
    print_compile();
    close(FILE);
    unlink("makefile.in.old");
    chdir("..") || die "Can't cd to ..\n";
    print " done.\n";
}


### getdeps ###
sub getdeps{
    %depstmp=();
    my @files_tmp = @_;
    foreach $file (@files_tmp) {
        open(DFILE,"$file");
        while(defined($line=<DFILE>)) {
	    if ($line =~ /include\s+'(.+)'/i) {
	        $depstmp{$1}++;
    	    }
        }
        close(DFILE);
    }
    return (keys %depstmp);
}


### print_header ###
sub print_header{
print FILE <<END;
# Makefile for $dir directory
# Author: Joakim Edsjo, edsjo\@physto.se
# This file is automatically created by makemfcat.pl on $date.

# Define fortran compiler and options (set when ./configure is run
# in the DarkSUSY root directory
FF=\@F77\@
FOPT=\@FOPT\@

# Dependencies and libraries
DINC=../../include

END
}

### print_line ###
sub print_line{
    my $line=$_[0];
    $cols=72;
    while(length($line) != 0) {
	if (length($line)>$cols) {
            $i=rindex($line," ",$cols);
	    print FILE substr($line,0,$i);
	    print FILE " \\\n";
            substr($line,0,$i+1)="";
	} else {
	    print FILE "$line\n";
	    $line="";
	}
    }
    print FILE "\n";
}

### print_compile ###
sub print_compile{
    print FILE "ds$dir.o : \$($dir) \$(INC_DEP) makefile\n";
    print FILE "\tcat \$($dir) > ds$dir$suf\n";
    print FILE "\t\$(FF) \$(FOPT) -c -I\$(DINC) -o ds$dir.o ds$dir$suf\n";
    print FILE "\trm ds$dir$suf\n";
}






