#!/usr/bin/perl -w
# Last-modified: 01 May 2011 07:33:30 PM

# @ChangeLog
# Feb-06-2011: Written by Ying Zu.

#------------------------------
# Subroutines
#------------------------------
sub init()
{
    use Getopt::Std;
    my $opt_string = 'hm:c:g:a:s:t:';
    getopts( "$opt_string", \%opt ) or usage();
    usage() if $opt{h};
}
#------------------------------
sub usage()
{
    print STDERR << "EOF";
usage: $0 [-h] [-m file] [-c num] [-g file] [-a x] [-s file] [-t stat] 
 -h            : this (help) message
 -m foo.dat    : mcmcoutput file
 -c num        : column number of the variable
 -g bar.dat    : gridoutput file
 -a x,y        : which variable (axis)
 -s foo.dat    : write a stats.dat for single mode mcmcoutput file
 -t stats.dat  : name of the sigma-tau prior file, default: stats.dat
example: $0 -m file -c num
example: $0 -g file -a x
example: $0 -s file -t stats.dat
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

#------------------------------
# begin main pro
#------------------------------
init();
#------------------------------
# pre-determined formats/files/executables
#------------------------------
if($opt{m}){
    $imcmc = 1;
    $igrid = 0;
    $istat = 0;
    $mcmcoutput = $opt{m};
    die "ERROR: no mcmcoutput data: $mcmcoutput\n" if(! -e $mcmcoutput);
    if($opt{c}){
        $colid = $opt{c};
        print STDOUT "Using column: $colid\n";
    }else{
        print STDERR "Please indicate a column id\n";
        usage();
    }
}elsif($opt{g}){
    $imcmc = 0;
    $igrid = 1;
    $istat = 0;
    $gridoutput = $opt{g};
    die "ERROR: no gridoutput data: $gridoutput\n" if(! -e $gridoutput);
    if($opt{a}){
        $axeid = $opt{a};
        print STDOUT "Using Axis: $axeid\n";
    }else{
        print STDERR "Please indicate an axis id\n";
        usage();
    }
}elsif($opt{s}){
    $imcmc = 0;
    $igrid = 0;
    $istat = 1;
    $mcmcoutput = $opt{s};
    die "ERROR: no mcmcoutput data: $mcmcoutput\n" if(! -e $mcmcoutput);
    if($opt{t}){
        $statsfile = $opt{t};
        print STDOUT "Writing to: $statsfile\n";
    }else{
        $statsfile = "stats.dat";
        print STDOUT "Writing to: $statsfile\n";
    }
}else{
    usage();
    die "ERROR: neither mcmc nor grid file is provided\n";
}
#------------------------------
# MCMC CL Determination
#------------------------------
if ($imcmc == 1){
    $entid  = $colid - 1;
    $npt    = 0;
    open MCFILE,"<$mcmcoutput";
    while(<MCFILE>) {
        $line         = $_;
        @entries      = split ' ',$line;
        $item[$npt] = $entries[$entid];
        $npt++;
    }
    close(MCFILE);
    $nmed    = $npt/2;
    $nlow    = 0.158*$npt;
    $nhig    = (1-0.158)*$npt;
    @sortitem= sort { $a <=> $b } @item;
    $meditem = $sortitem[$nmed];
    $lowitem = $sortitem[$nlow];
    $higitem = $sortitem[$nhig];
    printf STDOUT "CL: %13f %13f %13f\n",$meditem,$lowitem,$higitem;
}
#------------------------------
# GRID CL Determination
#------------------------------
if ($igrid == 1){
    if ($axeid eq 'x'){
        $iaxis = 1;
    }elsif($axeid eq 'y'){
        $iaxis = 2;
    }else{
        print STDERR "axis id should be either x or y";
        usage();
    }
    open OPENGRID, "<$gridoutput";
    # this should generally be fine, unless your grid has run-away
    # values, which most possibly is a bug.
#    $pmax      = -1000000.0;
    $pmax      = -1e+32;
    $xlen      = 0;
    $ylen      = 0;
    $xstart    = 5;
    $ystart    = 5;
    $npt    = 0;
# this is for 2-D grid only, 2 variables.
    $tolerance3sig = 11.8/2.0;
    $tolerance2sig = 6.17/2.0;
    $tolerance1sig = 2.30/2.0;
    while(<OPENGRID>){
        chomp;
        @entries         = split ' ';
        $npt++;
        $pcurrent        = $entries[0];
# all my grid files are indexed from zero, but due to my previous perl coding styles, perl array and sub only deal with indice starting from one. so...
        $xid             = $entries[3] + 1;
        $yid             = $entries[4] + 1;
        $x2d[$xid][$yid] = $entries[1];
        $y2d[$xid][$yid] = $entries[2];
        $plike2d[$xid][$yid] = $pcurrent;
        if($pcurrent > $pmax){
            $pmax    = $pcurrent;
            $maxxpos = $xid;
            $maxypos = $yid;
        }
        if($xid > $xlen){   $xlen   = $xid; }
        if($yid > $ylen){   $ylen   = $yid; }
        if($xid < $xstart){ $xstart = $xid; }
        if($yid < $ystart){ $ystart = $yid; }
    }
    close(OPENGRID);
# double check the grid dimension
    $xstart = $xstart - 1;
    $ystart = $ystart - 1;
    if (($xstart != 0) || ($ystart != 0)){
        printf STDERR "ERROR: Grid indices start from non-zero\n";
        die "ERROR: Not a Standard zygrid file\n";
    }
    if ($npt != $xlen*$ylen){
        die "ERROR: File Length does not match Grid Dimen\n";
    }
    $pcut1 = $pmax - $tolerance1sig;
    $pcut2 = $pmax - $tolerance2sig;
    $pcut3 = $pmax - $tolerance3sig;
    for($i = 1; $i <= $xlen; $i++){
        for($j = 1; $j <= $ylen; $j++){
            $idplike2d[$i][$j] = Assign123($plike2d[$i][$j],$pcut1,$pcut2,$pcut3);
        }
    }
    PlotGrid(\@idplike2d,$xlen,$ylen,\@x2d,\@y2d);
# find 3 sigma range (find the range of x ory enveloped by 1)
    $pout = 1;
    if($iaxis == 1){
        $xcen    = $x2d[$maxxpos][$maxypos];
        $uppxpos = $maxxpos;
        $lowxpos = $maxxpos;
        # find the upper bound
        for($i = $maxxpos; $i <= $xlen; $i++){
            for($j = $maxypos; $j <= $ylen; $j++){
                if($idplike2d[$i][$j] >= $pout){
                    if($i > $uppxpos){$uppxpos = $i;}
                }
            }
            for($j = $maxypos; $j >= 1; $j--){
                if($idplike2d[$i][$j] >= $pout){
                    if($i > $uppxpos){$uppxpos = $i;}
                }
            }
        }
        # find the lower bound
        for($i = $maxxpos; $i >= 1; $i--){
            for($j = $maxypos; $j <= $ylen; $j++){
                if($idplike2d[$i][$j] >= $pout){
                    if($i < $lowxpos){$lowxpos = $i;}
                }
            }
            for($j = $maxypos; $j >= 1; $j--){
                if($idplike2d[$i][$j] >= $pout){
                    if($i < $lowxpos){$lowxpos = $i;}
                }
            }
        }
        $xupp = $x2d[$uppxpos][1];
        $xlow = $x2d[$lowxpos][1];
        printf STDOUT "CL: %13f %13f %13f\n",$xcen,$xlow,$xupp;
    }
    if($iaxis == 2){
        $ycen    = $y2d[$maxxpos][$maxypos];
        $uppypos = $maxypos;
        $lowypos = $maxypos;
        # find the upper bound
        for($j = $maxypos; $j <= $ylen; $j++){
            for($i = $maxxpos; $i <= $xlen; $i++){
                if($idplike2d[$i][$j] >= $pout){
                    if($j > $uppypos){$uppypos = $j;}
                }
            }
            for($i = $maxxpos; $i >= 1; $i--){
                if($idplike2d[$i][$j] >= $pout){
                    if($j > $uppypos){$uppypos = $j;}
                }
            }
        }
        # find the lower bound
        for($j = $maxypos; $j >= 1; $j--){
            for($i = $maxxpos; $i <= $xlen; $i++){
                if($idplike2d[$i][$j] >= $pout){
                    if($j < $lowypos){$lowypos = $j;}
                }
            }
            for($i = $maxxpos; $i >= 1; $i--){
                if($idplike2d[$i][$j] >= $pout){
                    if($j < $lowypos){$lowypos = $j;}
                }
            }
        }
        $yupp = $y2d[1][$uppypos];
        $ylow = $y2d[1][$lowypos];
        printf STDOUT "CL: %13f %13f %13f\n",$ycen,$ylow,$yupp;
    }
}

#------------------------------
# generate a stats.dat file
#------------------------------

if($istat == 1) {                                                    
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

    $format6     = " %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n"; 
                                                                        
    open(OUTFILE,">$statsfile");                                       
    printf OUTFILE $format6,$meditem1,$lowitem1,$higitem1,$meditem2,$lowitem2,$higitem2;
    printf $format6,$meditem1,$lowitem1,$higitem1,$meditem2,$lowitem2,$higitem2;
    close(OUTFILE);                                                     
    printf STDOUT "did $mcmcoutput with $npt points\n";                 
}                                                       

#------------------------------




