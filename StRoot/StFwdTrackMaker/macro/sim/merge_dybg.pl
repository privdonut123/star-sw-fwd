#!/usr/bin/perl -w

$dir="hist";
$name="pythia_dybg_vz0";

@hist = ("dilep","match");
#@hits =("dilep");
#@hits =("match");                                                                                                              

foreach $h (@hist) {
    $nfile=0;
    $out="$name.$h.root";
    $cmd="cd $dir; hadd -f $out ";
    if(-e "$out") {`/bin/rm $out`};
    
    print "Searching $dir for .$h.root\n";
    
    opendir(DIR,"$dir");
    my @files1= grep {/.$h.root/} readdir DIR;
    my @files2= grep {/$name/} @files1;
    $n=0;
    $nn=0;
    foreach $file (@files2) {
	if($n == 0){
	    $outn="$name.$h.$nn.root";
	    $cmdn="cd $dir; hadd -f $outn ";	
	}
	$cmdn = "$cmdn $file";
	$n++;
	$nfile++;
	print "nfile=$nfile grp=$nn n=$n for $name in $dir/$name\n";
	if($n==100){
	    print "cmdn=$cmdn\n";
	    $out=`$cmdn`;
	    print "$out\n";
	    $cmd = "$cmd $name.$h.$nn.root";
	    $n=0;
	    $nn++;
	}
    }
    print "cmd=$cmd\n";
    $out=`$cmd`;
    print "$out\n";
}
