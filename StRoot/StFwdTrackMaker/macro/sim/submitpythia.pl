#!/usr/bin/perl -w

#@pid = ("dy","JPsi","dybg"); 
#@pid = ("dy","JPsi"); 
@pid = ("dybg");
#@pid = ("dy"); 
#@pid = ("JPsi");
#@pid = ("jet");
@vz = ("0");
$indir="fcs2022";
$outdir="hist";
$nev = 1000;
#$nev = 10;
$run1= 1;
#$run2= 3;
#$run2= 100;
$run2= 500;
$nogeant=1;

my $opt="none";
if($#ARGV==0) {
    $opt=$ARGV[0];
    print "Option = $opt\n";
}else{
    print "Option = none, just creating job file. To submit job, add option submit\n";
}

my $exe="runpythia";
my $condordir="condor";
my $cwd=$ENV{'PWD'};
my $logdir="$outdir/log";

if (! -e "$outdir") {system("/bin/mkdir $outdir");}
if (! -e "$logdir") {system("/bin/mkdir $logdir");}
if (! -e "$condordir") {system("/bin/mkdir $condordir");}
my $condor="$condordir/submit_pythia.txt";
if (-e $condor) {system("/bin/rm $condor");}
print("Creating $condor\n");
open(OUT, "> $condor\n");
print(OUT "Executable   = $exe\n");
print(OUT "Universe     = vanilla\n");
print(OUT "notification = never\n");
print(OUT "getenv       = True\n");
print(OUT "Accounting_group = group_star.cas\n");
print(OUT "\n");

$njob=0;
foreach $p (@pid) {
    if( $p eq "dybg") {
	if( $run2 == 100 ) {$run2=500;}
    }
    foreach $vz (@vz) {	
	for($r=$run1; $r<=$run2; $r++){
	    $name="pythia_${p}_vz$vz";
	    $in="$indir/$name";
	    $out="$outdir";
	    $log="$cwd/$logdir/$name.run$r.log";
	    if (! -e "$cwd/$out")     {system("/bin/mkdir $cwd/$out");}
	    if (! -e "$cwd/$logdir")  {system("/bin/mkdir $cwd/$logdir");}
	    
	    print(OUT "Arguments = \"$nev $r $p $vz $in $out $nogeant\"\n");
	    print(OUT "Log    = $log\n");
	    print(OUT "Output = $log\n");
	    print(OUT "Error  = $log\n");
		print(OUT "Queue\n\n");    
	    $njob++;
	}
    }
}

print "$njob jobs created\n";

if($opt eq "submit" ){
    print("Submitting ${condor}\n");
    system("condor_submit ${condor}\n");
    system("running200.pl\n");
}
