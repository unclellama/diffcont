#!/usr/bin/perl -w
# Last-modified: 02 May 2011 08:27:30 PM
use File::Copy;

# @Author:
# Ying Zu

# @ChangeLog:
# Apr-29-2011 Major revision on options, and less dependence on the 
# file existence (Columbus, OH)
# Feb-08-2011 Added Redshift (Columbus, OH)
# Dec-29-2010 Fixed several bugs in msg printing (Columbus, OH)

# @Description:
# This script is designed to be an automate pipeline for analyzing 
# light curve variabilities with SPEAR. The sequence of the pipeline is
# as follows,

# 1. GRID (required)
# 1.1. Grid calculation of likelihood.
# 1.2. Identify peak on the grid.
 
# 2. PRED (optional)
# 2.1. Predict a simulated high resolution light curve for the data.

# 3. MCMC (optional)
# 3.1. Run an MCMC chain around the identified peak.
# 3.2. Find one sigma range of tau and sigmahat.


#------------------------------
#SUBROUTINES
#------------------------------
sub init()
{
    use Getopt::Std;
    my $opt_string = 'hvwaf:r:z:s:g:d:pm';
    getopts( "$opt_string", \%opt ) or usage();
    usage() if $opt{h};
}
#------------------------------
sub usage()
{
    printf STDERR << "EOF";
usage: $0 -f file [-h -v -w -a -r -p -m rootname -z redshift -s spear -g s1_s2_ds_t1_t2_dt -d grid.dat]
 -h         : this (help) message
 -v         : turn on verbose mode for script output, defaut: off
 -w         : turn on verbose mode for spear output, default: off
 -a         : show an ASCII plot of the likelihood grid on sigma-tau plane, default: off
 -f foo.dat           : light curve data (only the continuum)
 -r rootname          : root name for all the output files, default: null string
 -z redshift: redshift of the object for converting all time to rest-frame, default: 0
 -s pathtospear       : path to the spear executable, default: spear

 -g s1_s2_ds_t1_t2_dt : specify your own grid geometry in a '_' concatenated string, default: -3_1_41_0_4_41
 -d pathtogrid        : use existing grid file without running calculation, default: none
 -p         : predict hires light curves at the peak, default: off
 -m         : run a MCMC chain, default: off
example: $0 -f lc_loopdeloop_con.dat -v -p -m -a
EOF
    exit;
}
#------------------------------                                         
sub Assign123{                                                          
    my ($pnow) = shift;                                                 
    my ($pcut1)  = shift;                                               
    my ($pcut2)  = shift;                                               
    my ($pcut3)  = shift;                                               
    my ($v123);                                                         
    if($pnow > $pcut3){                                                 
        if($pnow > $pcut2){                                             
            if($pnow > $pcut1){                                         
                $v123 = 3;                                              
            }else{                                                      
                $v123 = 2;                                              
            }                                                           
        }else{                                                          
            $v123 = 1;                                                  
        }                                                               
    }else{                                                              
        $v123 = 0;                                                      
    }                                                                   
    return $v123;                                                       
}                                
#------------------------------                                         
sub PlotGrid{                                                           
    my ($ref) = shift;                                                  
    my ($nx)  = shift;                                                  
    my ($ny)  = shift;                                                  
    my ($xaxis)  = shift;                                               
    my ($yaxis)  = shift;                                               
    my ($i);                                                            
    my ($j);                                                            
    printf STDOUT "print the grid...\n";                                
    printf STDOUT "\n";                                                 
    for($j = $ny; $j > 0; $j--){                                        
        printf STDOUT "%-6.2f", $yaxis->[1][$j];                        
        for($i = 1; $i <= $nx; $i++){                                   
            printf STDOUT $ref->[$i][$j];                               
        }                                                               
        printf STDOUT "\n";                                             
    }                                                                   
    printf STDOUT "\n      ";                                           
    for($i = 1; $i <= $nx; $i++){                                       
        $imod5 = ($i-1) % 5;                                            
        if($imod5 == 0){                                                
            printf STDOUT '|    ';                                      
        }                                                               
    }                                                                   
    printf STDOUT "\n      ";                                           
    for($i = 1; $i <= $nx; $i++){                                       
        $imod5 = ($i-1) % 5;                                            
        if($imod5 == 0){                                                
            printf STDOUT "%-5.1f", $xaxis->[$i][1];                    
        }                                                               
    }                                                                   
    printf STDOUT "\n\n";                                               
}                                                                       
#------------------------------   


#------------------------------
#   Main Program
#------------------------------
init();

$format5     = " %13.6f %13.6f %13.6f %4d %4d\n";
$format6     = " %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n";

#------------------------------
# translate options, user-defined files/executables
if(! $opt{v}){ open STDOUT, '>', "/dev/null" }

$clarg = " ";
if(! $opt{w}){  $clarg = "> /dev/null" }

$iplotgrid = 0;
if($opt{a}){ 
    die "ERROR: require verbose mode (-v) to plot ASCII grid\n" if(! $opt{v});
    $iplotgrid = 1; 
}

if($opt{f}){
    $lc_dat = $opt{f};
    die "ERROR: no light curve data: $lc_dat\n" if(! -e $lc_dat);
    printf STDOUT "Using light curve data: $lc_dat\n";
}else{
    printf STDERR "ERROR: Please indicate a Light Curve Dat File!\n";
    usage();
}

if($opt{r}){
    $rootname = $opt{r};
    $gridoutput  = $rootname . '_grid.dat';
    $bestoutput  = $rootname . '_best3.txt';
    $predoutput  = $rootname . '_hires.dat';
    $mcmcoutput  = $rootname . '_mcmc.dat';
    $statoutput  = $rootname . '_stats.dat';
}else{
    $gridoutput  = 'grid.dat';
    $bestoutput  = 'best3.txt';
    $predoutput  = 'hires.dat';
    $mcmcoutput  = 'mcmc.dat';
    $statoutput  = 'stats.dat';
}

if($opt{z}){
    $redshift = $opt{z};
}else{
    $redshift = 0.0
}

if($opt{s}){
    $execute = $opt{s};
}else{
    $execute = "spear";
}

# check if spear is in $PATH using system "which"
system("which $execute") == 0 or die "ERROR: no $execute in current PATH:$?";

$igrid = 1;
if($opt{d}){
    die "ERROR: no -g option" if ($opt{g});
    $gridoutput = $opt{d};
    die "ERROR: cannot find $gridoutput file" if (! -e $gridoutput);
    # no grid calculation
    $igrid = 0;
}

printf STDOUT "setting up the sig-tau grid in log10 space...\n";

if($opt{g}){
    $geometry_str = $opt{g};
    @geometry = split '\_', $geometry_str;
    die "ERROR: cannot parse the -g option" if ($#geometry + 1 != 6);
    $sigmin = $geometry[0];
    $sigmax = $geometry[1];
    $sigstp = $geometry[2];
    $taumin = $geometry[3];
    $taumax = $geometry[4];
    $taustp = $geometry[5];
    printf STDOUT "$sigmin $sigmax $sigstp\n";
    printf STDOUT "$taumin $taumax $taustp\n";
    $ns = ($sigmax - $sigmin)/$sigstp + 1;
    $nt = ($taumax - $taumin)/$taustp + 1;
    if (not $ns =~ /^-?\d+$/) {die "num of sig stps not an integer\n"}
    if (not $nt =~ /^-?\d+$/) {die "num of tau stps not an integer\n"}
}else{
    # default grid
    $ns         = 41;
    $sigmin     =-3.0;
    $sigmax     = 1.0;
    $nt         = 41;
    $taumin     = 0.0;
    $taumax     = 4.0;
    $sigstp     = ($sigmax-$sigmin)/($ns - 1);
    $taustp     = ($taumax-$taumin)/($nt - 1);
    printf STDOUT "$sigmin $sigmax $sigstp\n";
    printf STDOUT "$taumin $taumax $taustp\n";
    printf STDOUT "$ns $nt\n";
}

printf STDOUT "Running log10(sig) from $sigmin to $sigmax with step size $sigstp\n";
printf STDOUT "Running log10(tau) from $taumin to $taumax with step size $taustp\n";

$ipred = 0;
if($opt{p}){ $ipred = 1; }

$imcmc = 0;
if($opt{m}){ $imcmc = 1; }

#------------------------------

if ($igrid == 1){
    $pmax = -1000000.0;
    open(RESULTS,">$gridoutput");
    for ($j=1; $j<=$nt; $j++) {
        $tau    = $taumin + ($j-1)*$taustp;
        # assuming the tau given from CL is in rest-frame, but we want
        # obs-frame here, tauuse
        $tauuse = (10.0**$tau)*(1.0+$redshift);
        for ($i=1; $i<=$ns; $i++) {
            $sigma    = $sigmin + ($i-1)*$sigstp;
            # again, assuming incoming sigma is in rest-frame
            $sigmause = (10.0**$sigma)/sqrt(1.0+$redshift);
            open(OTFILE,">stinput.single");
            printf OTFILE "$lc_dat\n";
            printf OTFILE "2 single mode\n";
            printf OTFILE "$sigmause $tauuse\n";
            printf OTFILE "0 0\n";
            printf OTFILE "1 do amoeba\n";
            printf OTFILE "0 predict light curves\n";
            printf OTFILE "0 do mcmc\n";
            close(OTFILE);
            system("$execute < stinput.single $clarg");
            open(INFILE,"<best3.new");
            $modeline  = <INFILE>; undef $modeline;
            $paramline = <INFILE>; undef $paramline;
            $likeline  = <INFILE>;  @entries = split ' ',$likeline ; $like = $entries[0];
            close(INFILE);
            # the current format is to comply with python meshgrid
            printf RESULTS $format5,$like,$sigma,$tau,$i-1,$j-1;
            # post-processing: peak finding
            if($like > $pmax){
                $pmax  =  $like;
                $speak =  $sigma;
                $tpeak =  $tau;
            }
            $plike2d[$i][$j] = $like;
            $sigma2d[$i][$j] = $sigma;
            $tau2d[$i][$j]   = $tau;
        }
    }
    close(RESULTS);
}else{
    $pmax = -1000000.0;
    open(OPENGRID,"<$gridoutput");
    while(<OPENGRID>){
        chomp;
        @entries  = split ' ';
        $pcurrent = $entries[0]; 
        $xid                 = $entries[3]+1;
        $yid                 = $entries[4]+1;
        $plike2d[$xid][$yid] = $pcurrent;
        $sigma2d[$xid][$yid] = $entries[1];
        $tau2d[$xid][$yid]   = $entries[2];
        if($pcurrent > $pmax){
            $pmax = $pcurrent;
            $speak =  $entries[1];
            $tpeak =  $entries[2];
        }
    }
    close(OPENGRID);
}


if ($iplotgrid == 1){                                                   
    $tolerance3sig = 11.8/2.0;                                          
    $tolerance2sig = 6.17/2.0;                                          
    $tolerance1sig = 2.30/2.0;                                          
    $pcut1 = $pmax - $tolerance1sig;                                    
    $pcut2 = $pmax - $tolerance2sig;                                    
    $pcut3 = $pmax - $tolerance3sig;                                    
    for($i = 1; $i <= $ns; $i++){                                       
        for($j = 1; $j <= $nt; $j++){                                   
            $idplike2d[$i][$j] = Assign123($plike2d[$i][$j],$pcut1,$pcut2,$pcut3);
        }                                                               
    }                                                                   
                                                                        
    PlotGrid(\@idplike2d,$ns,$nt,\@sigma2d,\@tau2d);                    
}


printf STDOUT "maximum log likelihood: $pmax\n";
# the peak values are in rest-frame, we want them to be in obs-frame
$speak = 10.**$speak/sqrt(1.0+$redshift);
$tpeak = 10.**$tpeak*(1.0+$redshift);
printf STDOUT "at obs-frame sig $speak\n";
printf STDOUT "at obs-frame tau $tpeak\n";

#---------------------------------------
if($ipred == 1){
    open   PRDINP, ">input.single.pred";
    printf PRDINP "$lc_dat\n";
    printf PRDINP "2 single mode\n";
    printf PRDINP "$speak $tpeak\n";
    printf PRDINP "1 1\n";
    printf PRDINP "1 do amoeba\n";
    printf PRDINP "1 predict light curves\n";
    printf PRDINP "0 do mcmc\n";
    close(PRDINP);
    system("$execute < input.single.pred $clarg");
    copy("fort.13","$predoutput") or die "Copy cont prediction failed: $!";
    copy("best3.new","$bestoutput") or die "Copy bestoutput failed: $!";
}
#---------------------------------------
if($imcmc == 1){
    open   QSOINP, ">input.single.mcmc";
    printf QSOINP "$lc_dat\n";
    printf QSOINP "2 single mode\n";
    printf QSOINP "$speak $tpeak\n";
    printf QSOINP "1 1\n";
    printf QSOINP "1 do amoeba\n";
    printf QSOINP "0 predict light curves\n";
    printf QSOINP "1 do mcmc\n";
    printf QSOINP "$mcmcoutput\n";
    close(QSOINP);
    system("$execute < input.single.mcmc $clarg");
}
#---------------------------------------
if(-e $mcmcoutput) {
    # open the mcmc output
    open INFILE,"<$mcmcoutput";
    $npt    = 0;
    #column number for tau
    $colid1 = 6;
    #column number for sigma
    $colid2 = 5;
    #comply with the 0-start array in perl
    $entid1 = $colid1 - 1;
    $entid2 = $colid2 - 1;
    while(<INFILE>) {
      $line         = $_;
      @entries      = split ' ',$line;
      $item1[$npt] = $entries[$entid1];
      $item2[$npt] = $entries[$entid2];
      $npt++;
    }
    close(INFILE);
    $nmed    = $npt/2;
    $nlow    = 0.158*$npt;
    $nhig    = (1-0.158)*$npt;
    
    @sortitem1= sort { $a <=> $b } @item1;
    $meditem1= $sortitem1[$nmed];
    $lowitem1= $sortitem1[$nlow];
    $higitem1= $sortitem1[$nhig];
    
    @sortitem2= sort { $a <=> $b } @item2;
    $meditem2= $sortitem2[$nmed];
    $lowitem2= $sortitem2[$nlow];
    $higitem2= $sortitem2[$nhig];
    
    open(OUTFILE,">$statoutput");
    printf OUTFILE $format6,$meditem1,$lowitem1,$higitem1,$meditem2,$lowitem2,$higitem2;
    printf $format6,$meditem1,$lowitem1,$higitem1,$meditem2,$lowitem2,$higitem2;
    close(OUTFILE);
    printf STDOUT "did $mcmcoutput with $npt points\n";
}
#------------------------------


