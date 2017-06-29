#!/usr/bin/perl -w
# Last-modified: 01 May 2011 06:44:09 PM
use File::Copy;

# @Author:
# Ying Zu

# @ChangeLog:
# May-01-2011 Use observed-frame time unit throughout, no redshift needed.
# Nov-11-2010 Added Connected Component Label sub to find peaks.

# @Description:
# This script is designed to be an automate pipeline for estimating
# quasar reverberation lags with SPEAR. The sequence of the pipeline is 
# as follows,

# 1. GRID (required)
# 1.1. Grid calculation of likelihood.
# 1.2. Identify peaks on the grid using a "Connected Component Labeling"
#      algorithm.

# 2. PRED (optional)                                                    
# 2.1. Predict a simulated high resolution light curve for the data at
#      each peak.
                                                                        

sub init()
{
    use Getopt::Std;
    my $opt_string = 'hvwaf:r:s:t:g:d:p';
    getopts( "$opt_string", \%opt ) or usage();
    usage() if $opt{h};
}
sub usage()
{
    print STDERR << "EOF";
usage: $0 -f file [-h -v -w -a -r rootname -s spear -t stats.dat -g t1_t2_dt_w1_w2_dw -p -d grid.dat]
 -h         : this (help) message
 -v         : verbose mode for script output
 -w         : verbose mode for spear output
 -a         : show an ASCII plot of the grid and peaks found

 -f foo.dat : light curve data (continuum + line)
 -r rootname          : root name for all the output files
 -s pathtospear       : path to the spear executable
 -t statsfile         : path to the prior file (if other than ./stats.dat 
                        or ../single/stats.dat)

 -g t1_t2_dt_w1_w2_dw : specify your own grid geometry in a '_' concatenated string
 -d pathtogrid        : use existing grid file without running calculation
 -p         : predict hires light curves at the peak

example: $0 -f lc_loopdeloop_con_hb.dat -v -a -g 0_50_1_0.01_10.01_10
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
    my ($lagaxis)  = shift;
    my ($widaxis)  = shift;
    my ($i);
    my ($j);
    printf STDOUT "print the grid...\n";
    printf STDOUT "\n";
    for($j = $ny; $j > 0; $j--){
        printf STDOUT "%-6.2f", $widaxis->[1][$j];
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
            printf STDOUT "%-5.1f", $lagaxis->[$i][1];
        }
    }
    printf STDOUT "\n\n";
}
#------------------------------
sub ContourToByte{
    my ($ref) = shift;
    my ($nx)  = shift;
    my ($ny)  = shift;
    my ($i);
    my ($j);
    for($i=1;$i<=$nx;$i++){
        for($j=1;$j<=$ny;$j++){
            if($ref->[$i][$j] > 0){
                $ref->[$i][$j]=1;
            }
        }
    }
}
#------------------------------
sub GetPeakNum{
    my ($ref) = shift;
    my ($nx)  = shift;
    my ($ny)  = shift;
    my ($i);
    my ($j);
    $npeaks = 0;
    for($i=1;$i<=$nx;$i++){
        for($j=1;$j<=$ny;$j++){
            if($ref->[$i][$j] > $npeaks){
                $npeaks = $ref->[$i][$j];
            }
        }
    }
    return $npeaks;
}
#------------------------------
sub BackgroundSubtraction{
    my ($foreimage) = shift;
    my ($backimage) = shift;
    my ($nx)  = shift;
    my ($ny)  = shift;
    my ($npeaks)    = shift;
    my ($peakpos)   = shift;
    my ($i);
    my ($j);
    my ($k);
    for($i=1;$i<=$npeaks;$i++){
        $blobmax[$i] = -1000000;
    }
    for($i=1;$i<=$nx;$i++){
        for($j=1;$j<=$ny;$j++){
            if($backimage->[$i][$j] > 0){
                for($k=1;$k<=$npeaks;$k++){
                    if($backimage->[$i][$j] == $k){
                        if($blobmax[$k] < $foreimage->[$i][$j]){
                            $blobmax[$k] = $foreimage->[$i][$j];
                            $peakpos->[$k][1] = $i;
                            $peakpos->[$k][2] = $j;
                            last;
                        }
                    }
                }
            }
        }
    }
}
#------------------------------
sub AddEquiv{
    my ($equiv) = shift;
    my ($i)     = shift;
    my ($j)     = shift;
    my ($k);
    if ($i==$j){
        return;
    }
    $k = $j;
    do {
      $k = $equiv->[$k];
    } while ( $k != $j && $k != $i );
    if ( $k == $j ) {
       $tmp = $equiv->[$i];
       $equiv->[$i] = $equiv->[$j];
       $equiv->[$j] = $tmp;
    }
}
#------------------------------
sub Connected8Component{
    my ($ref) = shift;
    my ($nx)  = shift;
    my ($ny)  = shift;
    # 1st pass counts max possible compts, 2nd records equivalences
    for ($pass = 1; $pass<=2; $pass++){
        if ($pass == 2){
            for($i=0;$i<=$newlabel;$i++){
                $equiv[$i]=$i;
            }
        }
        $newlabel = 1;
        for($i = 1 ; $i <= $nx ; $i++){
            for ($j = 1 ; $j <= $ny ; $j++ ){
                $nfound = 0; # number of neighbors
                $i1 = $i-1; $j1 = $j-1; $i2 = $i+1;
                if ($ref->[$i][$j] > 0){
                    if($i>1 && $ref->[$i1][$j]>0){
                        $neighbor[$nfound++] = $ref->[$i1][$j];
                    }
                    if($j>1 && $ref->[$i][$j1]>0){
                        $neighbor[$nfound++] = $ref->[$i][$j1];
                    }
                    if($j>1 && $i>1 && $ref->[$i1][$j1]>0){
                        $neighbor[$nfound++] = $ref->[$i1][$j1];
                    }
                    if($j>1 && $i<$nx && $ref->[$i2][$j1]>0){
                        $neighbor[$nfound++] = $ref->[$i2][$j1];
                    }
#if all four neighbors are 0, assign a new label to p, else
#if only one neighbor has V={1}, assign its label to p, else
#if more than one of the neighbors have V={1}, assign one of the labels to p and make a note of the equivalences.
                    if($nfound==0){
                        $ref->[$i][$j] = $newlabel++;
                    }else{
                        $ref->[$i][$j] = $neighbor[0];
                        if ($nfound>1 && $pass == 2){
                            for($k=1; $k<$nfound; $k++){
                                AddEquiv(\@equiv,$ref->[$i][$j],$neighbor[$k])
                            }
                        }
                    }
                }else{
                    $ref->[$i][$j]=0;
                }
            }
        }
    }
#    Replace each cycle by single label
    $count = 0;
    for ($i = 1; $i <= $newlabel; $i++){
        if($i <= $equiv[$i]) {
            $count++;
            $this=$i;
            while( $equiv[$this] != $i ){
                $next = $equiv[$this];
                $equiv[$this] = $count;
                $this = $next;
            }
            $equiv[$this] = $count;
        }
    }
#    Now remove equivalences 
    for($i=1; $i<$nx; $i++){
        for($j=1; $j<$ny; $j++) {
            $ref->[$i][$j]=$equiv[$ref->[$i][$j]] ;
        }
    }
}
#------------------------------


#------------------------------
# Main Program
#------------------------------
init();

$format8     = " %10.6f %10.6f %10.6f %4d %4d %10.6f %10.6f %10.6f\n";
$format6     = " %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f \n";
$statsfile   = "stats.dat";

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
    print STDERR "ERROR: Please indicate a Light Curve Dat File!\n";    
    usage();                                                            
}          

if($opt{r}){                                                            
    $rootname = $opt{r};                                                
    $gridoutput  = $rootname . '_grid.dat';
    $bestoutput  = $rootname . '_best3.txt';
    $predoutput  = $rootname . '_hires.dat';
    $peakoutput  = $rootname . '_peaks.grid';
    $repkoutput  = $rootname . '_peaks.pred';
}else{                                                                  
    $gridoutput  = 'grid.dat';
    $bestoutput  = 'best3.txt';
    $predoutput  = 'hires.dat';
    $peakoutput  = 'peaks.grid';
    $repkoutput  = 'peaks.pred';
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

printf STDOUT "setting up the lag-wid grid...\n";

if($opt{g}){
    $geometry_str = $opt{g};
    @geometry = split '\_', $geometry_str;
    die "ERROR: cannot parse the geometry" if ($#geometry + 1 != 6);
    $tmin = $geometry[0];
    $tmax = $geometry[1];
    $tstp = $geometry[2];
    $wmin = $geometry[3];
    $wmax = $geometry[4];
    $wstp = $geometry[5];
    printf STDOUT "$tmin $tmax $tstp\n";
    printf STDOUT "$wmin $wmax $wstp\n";
    $xlen = ($tmax - $tmin)/$tstp + 1;
    $ylen = ($wmax - $wmin)/$wstp + 1;
    if (not $xlen =~ /^-?\d+$/) {die "num of lag stps not an integer\n"}
    if (not $ylen =~ /^-?\d+$/) {die "num of wid stps not an integer\n"} 
}elsif($igrid == 1){
    # find the light curve span
    $timemin = 1000000;
    $timemax = 0;
    $ndata   = 0;
    open LCDAT, "<$lc_dat";
    while(<LCDAT>){
        chomp;
        @entries = split ' ';
        if($. == 1){
            $ncurve = $entries[0];
            die "ERROR: number of light curves different that 2" if($ncurve != 2);
        }
        $size = $#entries + 1;
        if($size == 3){
            $ndata   = $ndata + 1;
            $timenow = $entries[0];
            $timemin = $timenow if($timenow < $timemin); 
            $timemax = $timenow if($timenow > $timemax);
        }
        $#entries  = -1; # redef the array
    }
    close(LCDAT);
    
    $timerange = $timemax - $timemin;
    $timestep  = 2.0*$timerange/$ndata;
    printf STDOUT "light curves start at $timemin, and end at $timemax\n";
    printf STDOUT "range of light curves is $timerange\n";
    printf STDOUT "set the maximum lag to be half of the range of the light curves\n";
    printf STDOUT "set the step of lag to be the mean time gap between two points\n";
    $tmax = int($timerange/2.0);
    $tstp = int($timestep);
    $modmaxstp = $tmax % $tstp;
    $tmax = $tmax + $tstp - $modmaxstp if ($modmaxstp != 0);
    printf STDOUT "after pruning the modulus...\n";
    printf STDOUT "max lag time : $tmax\n";
    printf STDOUT "lag step size: $tstp\n";
    #------------------------------
    printf STDOUT "setting up the lag-width grid now...\n";
    $tmin = 0; # always be 0
    $xlen = ($tmax - $tmin)/$tstp + 1;
    $wmin = 0.001;  # slightly larger than 0 to avoid exotic results
    # simply go with 20 days, hard to gauge but should be within 20 days
    $wmax = $wmin + 19;
    $wstp = 1;
    $ylen = ($wmax - $wmin)/$wstp + 1;
}

if ($igrid == 1){
    printf STDOUT "Running lag from $tmin to $tmax with step size $tstp\n";
    printf STDOUT "Running wid from $wmin to $wmax with step size $wstp\n";
}

#------------------------------

if($opt{t}){
    $remotestatfile = $opt{s};
    copy("$remotestatfile","$statsfile") or die "Copy stats file failed: $!";
}elsif(-e $statsfile){
    printf STDOUT "using ./stats.dat as input for prior info\n";
}else{
    $sindir = '../single/';
    $statsfile = $sindir . $statsfile;
    if(! -e $sindir){
        printf STDERR "ERROR: cannot find stats.dat file, please set -t option\n";
        usage();
    }else{
        printf STDOUT "using ../single/stats.dat as input for prior info\n";
    }
}

open STAT,"<$statsfile";
$statline =  <STAT>;
@entries  =  split ' ',$statline;
$taucen   =  10.**$entries[0];
$sigcen   =  10.**$entries[3];
close(STAT);
printf STDOUT "initialize tau   to be: $taucen\n";
printf STDOUT "initialize sigma to be: $sigcen\n";

#------------------------------

$ipred = 0;                                                             
if($opt{p}){ $ipred = 1; }

#------------------------------

# scale quickly converges, so we don't care...
$scale= 1.0;

#------------------------------
if($igrid == 1){
    $pmax = -1000000.0;
    open(RESULTS,">$gridoutput");
    for($i = 1; $i <= $xlen; $i++){
        $lag = $tmin + ($i-1)*$tstp;
        for($j = 1; $j <= $ylen; $j++){
            $wid = $wmin + ($j-1)*$wstp;
            open(OTFILE,">lwinput.tophat");
            printf OTFILE "$lc_dat\n";
            printf OTFILE "1 tophat mode\n";
#            printf OTFILE "$sigcen $taucen $laguse $scale $widuse\n";
            printf OTFILE "$sigcen $taucen $lag $scale $wid\n";
            printf OTFILE "1 1 0 1 0\n";
            printf OTFILE "1 do amoeba\n";
            printf OTFILE "0 predict light curves\n";
            printf OTFILE "0 do mcmc\n";
            close(OTFILE);
            system("$execute < lwinput.tophat $clarg");
            open(INFILE,"<best3.new");
            $modeline  = <INFILE>; undef $modeline;
            $paramline = <INFILE>; undef $paramline;
            $likeline  = <INFILE>;  @entries = split ' ',$likeline ; $like = $entries[0];
            $chiline   = <INFILE>; undef $chiline;
            $ndofline  = <INFILE>; undef $ndofline;
            $ampline   = <INFILE>; undef $ampline;
            $mealine   = <INFILE>; undef $mealine;
            $linline   = <INFILE>; undef $linline;
            $params    = <INFILE>;  @entries = split ' ',$params   ; 
            # only scale is updated
            $sigma = $entries[0]; $tau = $entries[1]; $scale = $entries[3];
            close(INFILE);
            printf RESULTS $format8, $like,$lag,$wid,$i-1,$j-1,$sigma,$tau,$scale;
#            printf STDOUT $format8, $like,$lag,$wid,$i-1,$j-1,$sigma,$tau,$scale;
            # post-processing: peak finding
            if($like > $pmax){
                $pmax    = $like;
                $tpeak   = $lag;
                $wpeak   = $wid;
            }
            $plike2d[$i][$j] = $like;
            $lag2d[$i][$j]   = $lag;
            $wid2d[$i][$j]   = $wid;
            $sigma2d[$i][$j] = $sigma;
            $tau2d[$i][$j]   = $tau;
            $scale2d[$i][$j] = $scale;
        }
    }
    close(RESULTS)
}else{
    $pmax      = -1000000.0;
    $xlen      = 0;
    $ylen      = 0;
    open OPENGRID, "<$gridoutput";
    while(<OPENGRID>){
        chomp;
        @entries            = split ' ';
        $pcurrent           = $entries[0];
        $xid                = $entries[3]+1;
        $yid                = $entries[4]+1;
        $plike2d[$xid][$yid] = $pcurrent;
          $lag2d[$xid][$yid] = $entries[1];
          $wid2d[$xid][$yid] = $entries[2];
        $sigma2d[$xid][$yid] = $entries[5];
          $tau2d[$xid][$yid] = $entries[6];
        $scale2d[$xid][$yid] = $entries[7];
        $xlen = $xid if($xid > $xlen);
        $ylen = $yid if($yid > $ylen);
        if($pcurrent > $pmax){
            $pmax    = $pcurrent;
            $tpeak   = $entries[1];
            $wpeak   = $entries[2];
        }
    }
    close(OPENGRID);
}
#------------------------------

printf STDOUT "maximum log likelihood: $pmax\n";
printf STDOUT "at obs-frame lag $tpeak\n";                              
printf STDOUT "at obs-frame wid $wpeak\n"; 

$tolerance3sig = 11.8/2.0;
$tolerance2sig = 6.17/2.0;
$tolerance1sig = 2.30/2.0;
$pcut1 = $pmax - $tolerance1sig;
$pcut2 = $pmax - $tolerance2sig;
$pcut3 = $pmax - $tolerance3sig;
for($j = $ylen; $j > 0; $j--){
    for($i = 1; $i <= $xlen; $i++){
        $idplike2d[$i][$j] = Assign123($plike2d[$i][$j],$pcut1,$pcut2,$pcut3);
    }
}
if ($iplotgrid == 1){
    PlotGrid(\@idplike2d,$xlen,$ylen,\@lag2d,\@wid2d);
}
#------------------------------
# Peak finding 
$yedge = $ylen;
$xedge = $xlen;
while(1){
# two flag to test if there is a peak associated with aliasing problems
    $idangerx = 0;
    $idangery = 0;
    for($i = 1; $i <= $xlen; $i++){
        if($idplike2d[$i][$yedge] == 3){
            $idangery = 1;
            printf STDOUT "WARNING: detected danger of false peaks at wide width\n";
            last;
        }
    }
    for($j = 1; $j <= $ylen; $j++){
        if($idplike2d[$xedge][$j] == 3){
            $idangerx = 1;
            printf STDOUT "WARNING: detected danger of false peaks at long lag\n";
            last;
        }
    }
    # Change the ansic contour into a byte image with only values of 1 and 0
    ContourToByte(\@idplike2d,$xlen,$ylen);
    print "running connect component labeling...\n";
    Connected8Component(\@idplike2d,$xlen,$ylen);
    print "find number of peaks and their locations...\n";
    # initialize everything
    $npeak = GetPeakNum(\@idplike2d,$xlen,$ylen);
    print "number of peaks: $npeak\n";
    if(($npeak == 1)&&(($idangerx == 1)||($idangery == 1))){
        if($idangerx == 1){
            printf STDOUT "lag grid too long, eliminate the false blob and do it again\n";
            printf STDOUT "shrink the lag axis...";
            $xedgetemp = $xedge;
            for($j = 1; $j <= $yedge; $j++){
                for($i = 1; $i <= $xedge; $i++){
                    last if($idplike2d[$i][$j] == 1);
                }
                $xedgetemp = $i if($i < $xedgetemp);
            }
            $xedge = $xedgetemp - 1; # we donot expect anything exotic here, e.g., xedgetemp = 1
        }
        if($idangery == 1){
            printf STDOUT "wid grid too wide, eliminate the false blob and do it again\n";
            printf STDOUT "shrink the wid axis...";
            $yedgetemp = $yedge;
            for($i = 1; $i <= $xedge; $i++){
                for($j = 1; $j <= $yedge; $j++){
                    last if($idplike2d[$i][$j] == 1);
                }
                $yedgetemp = $j if($j < $yedgetemp);
            }
            $yedge = $yedgetemp - 1; # we donot expect anything exotic here, e.g., yedgetemp = 1
        }
        # refind the pmax
        $pmax = -1000000;
        for($i = 1; $i <= $xlen; $i++){
            for($j = 1; $j <= $ylen; $j++){
                if(($i > $xedge) || ($j > $yedge)){
                    $idplike2d[$i][$j] = 0;
                }else{
                    $pmax = $plike2d[$i][$j] if($plike2d[$i][$j] > $pmax);
                }
            }
        }
        $pcut1 = $pmax - $tolerance1sig;
        $pcut2 = $pmax - $tolerance2sig;
        $pcut3 = $pmax - $tolerance3sig;
        for($i = 1; $i <= $xedge; $i++){
            for($j = 1; $j <= $yedge; $j++){
                $idplike2d[$i][$j] = Assign123($plike2d[$i][$j],$pcut1,$pcut2,$pcut3);
            }
        }
        if ($iplotgrid == 1) {
            PlotGrid(\@idplike2d,$xlen,$ylen,\@lag2d,\@wid2d);
        }
    }else{
        # we are good
        for($i=1;$i<=$npeak;$i++){
            $peakpos[$i][1] = 0;$peakpos[$i][2] = 0;
        }
        BackgroundSubtraction(\@plike2d,\@idplike2d,$xlen,$ylen,$npeak,\@peakpos);
        if ($iplotgrid == 1) {
            PlotGrid(\@idplike2d,$xlen,$ylen,\@lag2d,\@wid2d);
        }
        last;
    }
}

open(PEAKS,">$peakoutput");
printf STDOUT "Peaks found:\b";
for($i=1;$i<=$npeak;$i++){
    printf STDOUT "Peak id: $i\n";
    $xpos = $peakpos[$i][1];
    $ypos = $peakpos[$i][2];
    printf STDOUT "Grid Location: $xpos,$xpos\n";
    printf STDOUT "Lag, Width: $lag2d[$xpos][$ypos],$wid2d[$xpos][$ypos]\n";
    printf PEAKS $format6,$plike2d[$xpos][$ypos],$lag2d[$xpos][$ypos],$wid2d[$xpos][$ypos],$sigma2d[$xpos][$ypos],$tau2d[$xpos][$ypos],$scale2d[$xpos][$ypos];
}
close(PEAKS);

#------------------------------
if(($ipred == 1)){
    open(REPKS,">$repkoutput");
    for($i  =  1; $i <= $npeak; $i++){
        $contpred  = $predoutput . ".cont." . "$i";
        $linepred  = $predoutput . ".line." . "$i";
        printf STDOUT "predicting light curves for peak $i\n";
        $xpos = $peakpos[$i][1];
        $ypos = $peakpos[$i][2];
        $peakinput = 'input.tophat.peak' . "$i";
        open(OTFILE,">$peakinput");
        printf OTFILE "$lc_dat\n";
        printf OTFILE "1 tophat mode\n";
        printf OTFILE "$sigma2d[$xpos][$ypos] $tau2d[$xpos][$ypos] $lag2d[$xpos][$ypos] $scale2d[$xpos][$ypos] $wid2d[$xpos][$ypos] \n";
        printf OTFILE "1 1 1 1 1\n";
        printf OTFILE "1 do amoeba\n";
        printf OTFILE "1 predict light curves\n";
        printf OTFILE "0 do mcmc\n";
        close(OTFILE);
        system("$execute < $peakinput $clarg");
        open(INFILE,"<best3.new");
        $modeline  = <INFILE>; undef $modeline;
        $paramline = <INFILE>; undef $paramline;
        $likeline  = <INFILE>;  @entries = split ' ',$likeline ; $like = $entries[0];
        $chiline   = <INFILE>; undef $chiline;
        $ndofline  = <INFILE>; undef $ndofline;
        $ampline   = <INFILE>; undef $ampline;
        $mealine   = <INFILE>; undef $mealine;
        $linline   = <INFILE>; undef $linline;
        $params    = <INFILE>;  @entries = split ' ',$params   ;
        $sigma = $entries[0]; $tau = $entries[1]; $lag = $entries[2]; $scale = $entries[3]; $wid = $entries[4]; 
        printf REPKS  $format6, $like,$sigma,$tau,$lag,$scale,$wid;
        printf STDOUT "Estimated Best Parameters Around Peak Position:\n";
        printf STDOUT $format6, $like,$sigma,$tau,$lag,$scale,$wid;
        close(INFILE);
        copy("fort.13","$contpred") or die "Copy cont prediction failed: $!";
        copy("fort.12","$linepred") or die "Copy line prediction failed: $!";
    }
    close(REPKS);
}
#------------------------------



