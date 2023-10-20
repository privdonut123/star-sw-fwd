#!/usr/bin/perl -w

@pid = ("dy","JPsi");
#@pid = ("dy");
#@pid = ("JPsi");
#@pid = ("dybg");

@hist = ("dilep","match");
#@hits =("dilep");
#@hits =("match");

$dir="hist";

foreach $h (@hist) {
    foreach $p (@pid) {
	$nfile=0;
	$name="pythia_${p}_vz0";
	$out="$name.$h.root";
	$cmd="cd $dir; hadd -f $out ";
	if(-e "$out") {`/bin/rm $out`};
	
	print "Searching $dir for ${name}_*.$h.root\n";
	opendir(DIR,"$dir");
	my @files1= grep {/.$h.root/} readdir DIR;
	my @files2= grep {/${name}_run/} @files1;
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
#	    print "nfile=$nfile grp=$nn n=$n in $dir/$name file=$file\n";
	    if($n==100){
		print "cmdn=$cmdn\n";
		$out=`$cmdn`;
		print "$out\n";
		$cmd = "$cmd $name.$h.$nn.root";
		$n=0;
		$nn++;
	    }
	}
	print "cmdn=$cmdn\n";
	$out=`$cmdn`;
	print "$out\n";
	$cmd = "$cmd $name.$h.$nn.root";
	
	print "cmd=$cmd\n";
	$out=`$cmd`;
	print "$out\n";
    }
}
