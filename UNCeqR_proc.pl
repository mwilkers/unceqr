use Math::CDF;
use Math::Cephes;

use strict;

use vars qw($rnaf $tdnaf $ndnaf $preMinTumorCov $posf $startChr $endChr $maf $sampId $tumorLoose $filterTumor $testReg $regionRestrict $minTMapQ $minTBaseQ $maxNM $maxIH $trimEnd $flagFilter $st $flagStrike $regChrLess $minInNormFrac $filterScript $maxDepth $header $mainProc @hisArgv $preMinNormCov $snpFile $dense $fastaWchr $fastaWoChr $tdnafChr $ndnafChr $rnafChr $minNormCnt $maxNormPlural $minDnaTumorCnt $maxDnaBias $maxDnaStruckProp $minRnaTumorCnt $maxRnaBias $maxRnaStruckProp $trainNum $resultsDir $maxDnaTumorPluralProp $maxRnaTumorPluralProp $dnaOnly $rnaOnly $maxHomopolymer $normalLoose $indelShadow $verboseOut $normIndelShadowFrac $mafEnd $medStart $mode $unceqrver $R $Rscript $R_LIBS_USER $localPath $pvCut);

my @snpf;
my $chri;
my $pos=-1;

print STDERR "#UNCeqR_proc: running\n";
my %revCol;

####################
# functions and global variables.
####################
	my $minDnaCoSupport = 1; #minimum DNA maxVal for RNA variant inclusion


	my $regPadd = 10; #upstream/downstream padding for region intervals.
	my $posDist = 20; #grouping for indel position shifting between DNA and RNA.

	my $cacheLen=$indelShadow; #needs to be odd numbered, unless off.
	if ($cacheLen % 2 ==0 && $cacheLen != 0){
	  $cacheLen++;
	}
	$posDist = ($posDist > ($cacheLen-1)/2) ? ($cacheLen-1)/2 : $posDist;


	my $normIndelFrac=$normIndelShadowFrac; #1/10;
	my @bam_nt16_rev_table = ("A","C","G","T","V","del","ins");
	
	my $Inf =10**10**10;
	my $backInd = 3; #for good, strike, medDistToEnd counts.
        my $pad=5; #define... #first column of data to use from mpileup_unceqr.
		#21 cols
		#6 leader cols, starting arr_ind 5;
		#10 for alleles ATCGV 
		#2 for del, starting arr ind 15
		#2 for ins, starting arr ind 17
		#good #last minus 1
		#struck #last
	my $ry=0;
	my $dy=0;
	my $ny=0;
	my $dhqCnt=0;
	my $rhqCnt=0;
	my $dtn=0;
	my $rtn=0;
	my $mtn=0;
	my $fitted=0;

	sub max {
		my @v = @_;
		if ($v[0] < $v[1]){
			return $v[1];
		}
		return $v[0];
	}

	sub getMeanNoZero {
		my @a=@_;
		my $t=0;
		my $c=0;
		foreach my $i (@a){
			if($i==0){next;}
			$t += $i;
			$c++;
		}
		return $t/$c;
	}

	sub getMostCommonIndel{
		my @ins = shift;
		my %tmph;
		for(my $i=0;$i<scalar(@ins);$i++){
			$tmph{$ins[$i]}++;
		}
		my $maxV=0;
		my $maxK=-1;
		foreach my $k (keys(%tmph)){
			if($tmph{$k}>$maxV){
				$maxV=$tmph{$k};
				$maxK=$k;
			}
		}
		return($maxV,$maxK);
	}


	sub parseLine{
		#input: one line from mpileup
		my @y = @_;
		my @x=@y[$pad..($#y-$backInd)];  #take out leading columns that are not character counts

		###
		# process deletions
		###
		my @del0 = $x[10] =~ /(\w+),/g;
		my @del1 = $x[11] =~ /(\w+),/g;
		my @dcall = getMostCommonIndel((@del0,@del1));
		$x[10] = scalar(@del0); #separate indel counts from both strands, allows for strand bias calculation of all indels.
		$x[11] = scalar(@del1);
		my $maxDel = $dcall[1];		

		###
		# process insertions
		###
		my @ins0 = $x[12] =~ /(\w+),/g;
		my @ins1 = $x[13] =~ /(\w+),/g;
		my @icall = getMostCommonIndel((@ins0,@ins1));
		$x[12] = scalar(@ins0); #separate indel counts from both strands, allows for strand bias calculation of all indels.
		$x[13] = scalar(@ins1);
		my $maxIns = $icall[1];		

		my @bias;
		#order is char,char,...
		my $j=0;
		my @nx;
		for(my $i=0;$i<$#x;$i=$i+2){
			$nx[$j] = $x[$i]+$x[$i+1];
			$j++;
		}
		#struck prop
		my $struckProp=($y[$#y-1]+$y[$#y-2]); #numGood, numStruck
		$struckProp = ($struckProp>0) ? $y[$#y-1]/$struckProp : 0;
		$struckProp = sprintf("%.3f",$struckProp,3);

		#median dist to end
		my $distToEnd = $y[$#y];


		my @r = (  \@nx,
			   \@x, $struckProp,
			   $maxIns, $maxDel,$distToEnd  );			  
		return \@r;
	}


	sub getRef{
	   #get normal bases and reference calls

	   #input are fields from parse
	   my @xt = @_;

	   #universal reference
	   my $ur = $xt[2];
	   $ur =~ tr/[atcg]/[ATCG]/; #to upper;

	   my $ln = parseLine(@xt);

           @xt=@{$ln->[0]};
	   my @bias=@{$ln->[1]};

	   my $majInd=-1; 
	   my $majVal=-1;
	   my $minInd=-1;
	   my $minVal=-1; 
	   my $plural=0;
	   my %nchars;
	   my $sum=0;
	   my $str;
	   #for each character
	   for(my $i=0;$i<scalar(@xt);$i++){
	     if ($xt[$i]>0){ 
	       $sum += $xt[$i];
	       $nchars{$bam_nt16_rev_table[$i]} += $xt[$i];
	       $str .= $bam_nt16_rev_table[ $i]."[".$xt[$i]."]";
	       $plural++;
	     }
	   }

	   #determine if normal indel, This fraction is hard coded.
	   my $normIndel=0;
	   if( $sum>0 && (($nchars{"ins"} + $nchars{"del"})/$sum > $normIndelFrac) ){
		$normIndel=1;
	   }

	   #remove rarely occuring variants
	   my %origChars = %nchars;
	   my $tSum=$sum;
	   foreach my $k (keys(%nchars)){
             if ($nchars{$k} < $sum * $minInNormFrac){ #less than 1%.
		$tSum -= $nchars{$k};
		delete($nchars{$k});
	     }
	   }
	   
	   #add reference as a normal allele
	   $nchars{$ur}++;

	   #add snps
	   addSNP(\%nchars);

	   return ($sum,\%nchars,$plural,$str,$tSum,$ln->[3],$ln->[4],$normIndel);
	}

	sub addSNP{
	my $nchars = shift;
	###
	# SNP file
	###
                if($snpFile ne "-1"){
			while(!eof(SNP) && ($snpf[0] != $chri ||  $snpf[1] < $pos ) ){ #advance to next if less than SNP
			    my $l = <SNP>;
			    chomp $l;
                            @snpf = split /\,/, $l;
			  }
			  if($pos==$snpf[1] && $snpf[0] == $chri){ 
			    my @snpChars = split /\-/, $snpf[2];
			    foreach my $ki (@snpChars){
				$nchars->{$ki}++;
			    }		
          		  }
                }
	}

	sub mergeArr{
	  #for pooled 
	  my ($a,$b) = @_;

	  #keep leading columns
	  my @n = @$a; 
	  
	  #add columns of counts
	  for(my $i=$pad;$i<10+$pad;$i++){
	    $n[$i] = $$a[$i] + $$b[$i];
	  }

	  #concat indels
	  for(my $i=10+$pad;$i<=13+$pad;$i++){
	    $n[$i] = $$a[$i].$$b[$i];
	  }	  

	  #add struck columns
	  for(my $i=14+$pad;$i<=$#n-1;$i++){
	    $n[$i] = $$a[$i] + $$b[$i];
	  }	  
	  #medDist
	  $n[$#n] = ($$a[$#n] > $$b[$#n]) ? $$a[$#n] : $$b[$#n];
 	 
	  #output for countRef;
	  return \@n;	  
	}


sub countRef{
	  my ($a,$ref) = @_;
	  my $refCnt=0;
	  my $refCntNoIndel=0;
	  my $nonRefCntNoIndel=0;  
	  my $maxVal=0;
	  my $maxInd=-1;
	  my $str="";
	  
	  my $nonRefCnt=0;
	  my %nonRefChar;

	  my $maxNonRefChar;
	  my $maxNonRefCount=0;
	  my $qv = 0;

	   my $ln = parseLine(@$a);
           my @x=@{$ln->[0]};
	   my @bias=@{$ln->[1]};
	   my $struckProp=$ln->[2];
	   my $medEndDist = $ln->[5];

	   my $plural=0;
	   #calculate non-refs     
	   #take max non-ref   
	  my @refBias;	   

	   for(my $i=0;$i<scalar(@x);$i++){
	     if($x[$i]==0 || $i==15){ #no base or N character.
	       next;
	     }
	     if (exists($ref->[1]->{ $bam_nt16_rev_table[ $i]}  )){ #in reference or normal
	       
	       $refCnt += $x[$i];
	       if($i < 15){ #excludes N
	         $refCntNoIndel += $x[$i];
	       }
	     }else{
	       if($maxVal < $x[$i]){
	         $maxVal=$x[$i];
	         $maxInd=$i;
	       }
	       $nonRefCnt += $x[$i];
	       if($i< 15){
	         $nonRefCntNoIndel += $x[$i];
	       }
	       $plural++;
	       $str .= $bam_nt16_rev_table[ $i]."[".$x[$i]."]"; 
	     }
	   }#end allele loop
	   
	   my $majNonRef = ($maxInd != -1) ? ($bam_nt16_rev_table[$maxInd] eq "ins") ? "ins".$ln->[3] : ($bam_nt16_rev_table[$maxInd] eq "del") ? "del".$ln->[4] : $bam_nt16_rev_table[$maxInd] : "";

	   ###
	   #bias
 	   ###
	   #logic here for majorcount
	   #subtract from total.
	   my $refPlus = 0;
	   my $refNeg = 0;
           @x = @bias;
	   for(my $i=0;$i<$#x;$i=$i+2){
	     $refPlus += $x[$i];
	     $refNeg += $x[$i+1];
	   }
	   my $varPlus = $x[$maxInd*2];
	   my $varNeg = $x[($maxInd*2)+1];
	   $refPlus -= $varPlus;
	   $refNeg -= $varNeg;
	
	   return(	$refCnt,$refCntNoIndel,$nonRefCnt,$nonRefCntNoIndel,
	   		$maxVal,$majNonRef,
	   		$plural,
	   		$str,
			$refPlus,$refNeg,$varPlus,$varNeg,
			$struckProp,
			$ln->[3],$ln->[4],
			$medEndDist
	   		 );     
}


my $pln; #number of init lines seen.
my $ptab = {};

sub printOut{
  my $str = shift;
  my $ex = shift;
  print OF $str;
}

sub fitDist{
  #fit
  my $cmd = "$Rscript $localPath/UNCeqR_fit.r $resultsDir/unceqr_proc.all.csv $dnaOnly";
  print STDERR "$cmd\n";
  print STDERR `$cmd 2>&1`;
  #exec("$cmd");

  #read in table, just DNA for now.
  open(IFp,"$resultsDir/unceqr_proc.all.csv.fit") || die "cannot open $resultsDir/unceqr_proc.all.csv.fit\n";
  my $maxn;
  while(<IFp>){
    chomp;

    my @f = split /\,/;
    $ptab->{$f[0]}->{a} = $f[1];
    $ptab->{$f[0]}->{b} = $f[2];

  }
  close IFp;
  return;
}


my $qdBias = Math::CDF::qchisq(1-$maxDnaBias,1);
my $qrBias = Math::CDF::qchisq(1-$maxRnaBias,1);

sub chisq2x2 {
  #description: 2x2 chi-square with yates continuity correction.
  my ($a,$b,$c,$d) = @_;

  #pseudocounts.
  $a++;
  $b++;
  $c++;
  $d++;

  my $n = $a+$b+$c+$d;
  
  my $statNum = $n * ( (abs($a*$d - $b*$c) - $n/2)**2) ;
  my $statDenom =($a+$c)*($b+$d)*($a+$b)*($c+$d);
  my $stat = ($statDenom > 0) ? $statNum / $statDenom : 9999999;
  return $stat;
}


sub betabinCDFdec{
  my ($n,$k,$a,$b) = @_;
  my $v=0;
  for (my $i=$n;$i>=$k;$i--){
    my $x = betabinPMFG($n,$i,$a,$b);
    $v = $v + $x;
  }
  return $v;
}

sub betabinPMFG{
  #CDF by gamma.
  my ($n,$k,$a,$b) = @_;
  my $b1 = Math::Cephes::lgam($n+1);
  my $b2 = Math::Cephes::lgam($k+1) + Math::Cephes::lgam($n - $k + 1);

  my $b3 = Math::Cephes::lgam($k+$a) + Math::Cephes::lgam($n-$k+$b);
  my $b4 = Math::Cephes::lgam($n+$a+$b);

  my $b5 = Math::Cephes::lgam($a+$b);
  my $b6 = Math::Cephes::lgam($a) + Math::Cephes::lgam($b);

  my $v = $b1 - $b2 + $b3 - $b4 + $b5 - $b6;

  return(exp($v));
}

my $minVal = 1e-16;

sub pinch {
  my $x = shift;
  #pinch p-value range
  if($x<=$minVal){
    return($minVal)
  }
  if($x >= 1-$minVal){
    return( 1 - $minVal)
  }
  return($x);
}

sub combP{
  my ($dp,$dn,$rp,$rn) = @_;
  my $mp = "NA";
  if($dp eq "NA" || $dp eq "99"){
    $mp = $rp;
  }elsif($rp eq "NA" || $rp eq "99"){
    $mp = $dp;
  }else{
    
    my $dq = Math::CDF::qnorm(1-$dp);
    my $rq = Math::CDF::qnorm(1-$rp);

    my $stat = ($dn**0.5 * $dq + $rn**0.5 * $rq);
    $stat = $stat / (($dn+$rn)**0.5);

    $mp = pinch(1 - Math::CDF::pnorm($stat));
    #print "$dq $rq $stat $mp\n";  
  }
  $mtn++;
  return $mp;
}

sub addP{
  my ($s,$n,$k) = @_;
  if(!exists($ptab->{$s}->{pv}->{$n}->{$k})){
    $ptab->{$s}->{pv}->{$n}->{$k} = pinch(betabinCDFdec($n,$k,$ptab->{$s}->{a},$ptab->{$s}->{b})   );
  }
  return $ptab->{$s}->{pv}->{$n}->{$k};
}

sub evalDat{
  #description: one row at a time.
  my $cache = shift;
  my $dp = "NA";
  my $rp = "NA";

  my ($dnaDepth,$rnaDepth);

 ##DNA
 if(
	$cache->[1] == 0 && #indelFlag
	$cache->[2] == 1  #dna Usable

   ){
    $dtn++;
  
    if ( $cache->[0]->[$revCol{"DNA_maxVal"}] == 0 ){
      $dp=pinch(1); #return 1 if no variant alleles; 
      $dnaDepth = $cache->[0]->[$revCol{"DNA_refCnt"}] + $cache->[0]->[$revCol{"DNA_maxVal"}];
    }else{
	$dnaDepth = $cache->[0]->[$revCol{"DNA_refCnt"}] + $cache->[0]->[$revCol{"DNA_maxVal"}];
	my $dnaPlural = ($cache->[0]->[$revCol{"DNA_nonRefCnt"}] >0) ? $cache->[0]->[$revCol{"DNA_maxVal"}] / $cache->[0]->[$revCol{"DNA_nonRefCnt"}] : 0;
	my $dnaBias = chisq2x2($cache->[0]->[$revCol{"DNA_refPlus"}],$cache->[0]->[$revCol{"DNA_refNeg"}] , $cache->[0]->[$revCol{"DNA_varPlus"}] , $cache->[0]->[$revCol{"DNA_varNeg"}] );     

	if(       
		$dnaPlural >= $maxDnaTumorPluralProp && #not excessive plurality of variant alleles
		$cache->[0]->[$revCol{"DNA_medDist"}] >= $medStart && 
		$dnaBias < $qrBias &&
		($cache->[0]->[$revCol{"DNA_majNonRef"}] !~ /ins|del/ || ($cache->[0]->[$revCol{"DNA_varNeg"}] > 0 && $cache->[0]->[$revCol{"DNA_varPlus"}] > 0))
	){
	    	$dp = addP("dna",$dnaDepth,$cache->[0]->[$revCol{"DNA_maxVal"}]);
	}else{
		$dp=99; #marker for maximum p-value;
	}
    }
  }


 ##RNA
 if(
	$cache->[1] == 0 && #indelFlag
	$cache->[3] == 1  #rna Usable
   ){
    $rtn++;
  
    if ( $cache->[0]->[$revCol{"RNA_maxVal"}] == 0 ){
      $rp=pinch(1); #return 1 if no variant alleles; 
      $rnaDepth = $cache->[0]->[$revCol{"RNA_refCnt"}] + $cache->[0]->[$revCol{"RNA_maxVal"}];
    }else{
	$rnaDepth = $cache->[0]->[$revCol{"RNA_refCnt"}] + $cache->[0]->[$revCol{"RNA_maxVal"}];
	my $rnaPlural = ($cache->[0]->[$revCol{"RNA_nonRefCnt"}] >0) ? $cache->[0]->[$revCol{"RNA_maxVal"}] / $cache->[0]->[$revCol{"RNA_nonRefCnt"}] : 0;
	my $rnaBias = chisq2x2($cache->[0]->[$revCol{"RNA_refPlus"}],$cache->[0]->[$revCol{"RNA_refNeg"}] , $cache->[0]->[$revCol{"RNA_varPlus"}] , $cache->[0]->[$revCol{"RNA_varNeg"}] );     
	#my $rnaBias = 1;
	if(       
		$rnaPlural >= $maxRnaTumorPluralProp && #not excessive plurality of variant alleles
		$cache->[0]->[$revCol{"RNA_medDist"}] >= $medStart && 
		$rnaBias < $qrBias &&
		($cache->[0]->[$revCol{"RNA_majNonRef"}] !~ /ins|del/ || ($cache->[0]->[$revCol{"RNA_varNeg"}] > 0 && $cache->[0]->[$revCol{"RNA_varPlus"}] > 0))
	){
	    	$rp = addP("rna",$rnaDepth,$cache->[0]->[$revCol{"RNA_maxVal"}]);
	}else{
		$rp=99; #marker for maximum p-value;
	}
    }
  }

  ##Meta P-value;
  my $mp = (substr($cache->[0]->[$revCol{"DNA_majNonRef"}],0,3) eq substr($cache->[0]->[$revCol{"RNA_majNonRef"}],0,3)
		&& $cache->[0]->[$revCol{"DNA_maxVal"}]>=$minDnaCoSupport) ? combP($dp,$dnaDepth,$rp,$rnaDepth) : $dp;

  $cache->[4] = $dp;
  $cache->[5] = $rp;
  $cache->[6] = $mp;
  return($cache);
}


sub hqDataFlag{
  #description: evaluates high quality data rules for DNA and RNA reads, no variant allele consideration.
  my $cache = shift;
  my $df=0; 
  my $rf=0;

  my $normalDepth = $cache->[0]->[$revCol{"normCnt"}];

  #dna
  my $dnaDepth = $cache->[0]->[$revCol{"DNA_refCnt"}] + $cache->[0]->[$revCol{"DNA_nonRefCnt"}];
  if (
	($cache->[1] == 0 || $cacheLen==0)&& #not indel shadow
	$cache->[0]->[$revCol{"normCnt"}] >= $minNormCnt &&
 	$cache->[0]->[$revCol{"normPlural"}] <= $maxNormPlural &&
	$dnaDepth >= $minDnaTumorCnt && 
	$cache->[0]->[$revCol{"DNA_struckProp"}] < $maxDnaStruckProp
	){
	$df=1;
	$dhqCnt++;
  }
  #rna
  my $rnaDepth = $cache->[0]->[$revCol{"RNA_refCnt"}] + $cache->[0]->[$revCol{"RNA_nonRefCnt"}];
  if (
	($cache->[1] == 0 || $cacheLen==0)&& #not indel shadow
	$cache->[0]->[$revCol{"normCnt"}] >= $minNormCnt &&
 	$cache->[0]->[$revCol{"normPlural"}] <= $maxNormPlural &&
	$rnaDepth >= $minRnaTumorCnt && 
	$cache->[0]->[$revCol{"RNA_struckProp"}] < $maxRnaStruckProp
	){
	$rf=1;
	$rhqCnt++;
  }
  $cache->[2] = $df;
  $cache->[3] = $rf;
  return $cache;

}

my @cache;
my $fcache=0;
my $shadowPos = -1;
my $shadowChr = -1;

my @dn;
my @dk;
my @rn;
my @rk;

sub printOutCache{
  my $str = shift;
  my $normChar = shift;
  my $l;
  my $curPos;

  chomp $str;
  my @sa = split /\,/, $str;

  if ($cacheLen==0){
    @cache = [\@sa,$normChar];
    $curPos=0;
  }else{

   #position to be printed that is in the middle of the cache.
  $curPos = ($cacheLen-1)/2 - 1; #array index. 0-based. 
  $curPos = ($curPos > $#cache) ? $#cache : $curPos;

 #print STDERR "#$str $normChar $curPos $#cache $fcache $shadowPos $shadowChr $sa[0] $sa[1]";

  #add element
  my $tnormChar = $normChar;
  if($fcache>=1){
  	#update current parsed pos based on past.
    #check position that it within the cache distance, not just array index .
    if ( $sa[0] != $shadowChr || $sa[1] - $shadowPos > ($cacheLen-1)/2 ){
		$fcache=0;
		$shadowPos=$shadowChr=-1;
    }else{
      $tnormChar++;
      $fcache--; #update future from prior indel site.
    }
  }
  
  #add to cache.
  push(@cache,[\@sa,$tnormChar]);

  #update prior indel flag
  if($normChar==1){ #just added.
    $shadowPos=$cache[$#cache]->[0]->[1];
    $shadowChr = $cache[$#cache]->[0]->[0]; #chr~pos
    my $mm = (($cacheLen-1)/2 < $#cache) ? $#cache : ($cacheLen-1)/2; #floor at prior space available.
    REVCACHE: for(my $i=1;$i<=$mm && $i<$#cache;$i++){
	#check pos so that is not futher than distance of cache.
		if (!($cache[$#cache-$i]->[0]->[0] == $shadowChr && $shadowPos - $cache[$#cache-$i]->[0]->[1] <= ($cacheLen-1)/2 )){
			last REVCACHE;
		}
	    $cache[$#cache-$i]->[1] += 1; #increment shadow by one; could be larger than one if many indels nearby.
    }
    $fcache=($cacheLen)/2; #number of sites to look forward.
  }

  ##
  #indel align if needed at curPos.
  ##
  
    my $dmv = substr($cache[$curPos]->[0]->[13],0,3);#col i for DNA_majNonRef
    my $rmv = substr($cache[$curPos]->[0]->[27],0,3); #col i for RNA_majNonRef

    if(
	($dmv eq "ins" || $rmv eq "ins" || $dmv eq "del" || $rmv eq "del")
    ){

    my $dMax=0;
    my $rMax=0;
    my $dT="";
    my $rT="";
    my $rP=0;
    my $dP=0;

    #get max ins or del position within NEXT mm sites for RNA and DNA each
    #this aligns indel positions between DNA and RNA alignments; the current position is overwritten.

    my $mm = ($posDist + $curPos > $#cache) ? $#cache - $curPos : $posDist; #remaining space ceiling.

    my @dnaShoulder=();
    my @rnaShoulder=();
    my $dnaExtra=0;
    my $rnaExtra=0;

     for(my $ici=0;$ici<=$mm;$ici++){

        my $rc = substr($cache[$curPos+$ici]->[0]->[27],0,3); #RNA_str is column 27.
        if(
	    ($rc eq 'ins' || $rc eq 'del')){
             $rnaExtra += $cache[$curPos+$ici]->[0]->[26];
	     $rnaShoulder[++$#rnaShoulder] = $ici;
	 if(   $cache[$curPos+$ici]->[0]->[26] > $rMax
	   ){
			$rMax = $cache[$curPos+$ici]->[0]->[26];
			$rT = $rc;
			$rP = $ici;
	}}
        my $dc = substr($cache[$curPos+$ici]->[0]->[13],0,3); #DNA maxVal
        if(
	    ($dc eq 'ins' || $dc eq 'del')
	){
	   $dnaExtra += $cache[$curPos+$ici]->[0]->[12];
	   $dnaShoulder[++$#dnaShoulder] = $ici;
	if(    $cache[$curPos+$ici]->[0]->[12] > $dMax
	   ){
			$dMax = $cache[$curPos+$ici]->[0]->[12];
			$dT = $dc;
			$dP = $ici;
	}}
      }

      #move maxVal to maxSite from shoulder
      foreach my $icis (@dnaShoulder){
        if ($icis == $dP){ next;}

	if (substr($cache[$curPos+$icis]->[0]->[13],0,3) ne $dT){ next;} #must be ins or del same as max dna

	#subtract from germline count at max site
	$cache[$curPos+$dP]->[0]->[8] -= $cache[$curPos+$icis]->[0]->[12];

	#move maxVal count from shoulder to site.
        $cache[$curPos+$dP]->[0]->[12] += $cache[$curPos+$icis]->[0]->[12];

	#remove maxVal count from shoulder
	$cache[$curPos+$icis]->[0]->[12] = 0;
	print STDERR "#shoulder shift triggered at $cache[$curPos]->[0]->[0] $cache[$curPos]->[0]->[1], $dnaExtra \n";
      }
	$cache[$curPos+$dP]->[0]->[8] = max($cache[$curPos+$dP]->[0]->[8],0);

      #move maxVal to maxSite from shoulder
      foreach my $icis (@rnaShoulder){
        if ($icis == $rP){ next;}

	if (substr($cache[$curPos+$icis]->[0]->[27],0,3) ne $rT){ next;} #must be ins or del same as max rna

	#subtract from germline count at max site
	$cache[$curPos+$rP]->[0]->[22] -= $cache[$curPos+$icis]->[0]->[26];

	#move maxVal count from shoulder to site.
        $cache[$curPos+$rP]->[0]->[26] += $cache[$curPos+$icis]->[0]->[26];

	#remove maxVal count from shoulder
	$cache[$curPos+$icis]->[0]->[26] = 0;
	print STDERR "#shoulder shift triggered at $cache[$curPos]->[0]->[0] $cache[$curPos]->[0]->[1], $dnaExtra \n";
      }
	$cache[$curPos+$rP]->[0]->[22] = max($cache[$curPos+$rP]->[0]->[22],0);


      if($dP == 0 &&
	 $rP !=  0 &&
	 $rT eq $dT
      ){
	   my $ici = $rP;
	   print STDERR "#rna $rT $rP shift to current pos  $cache[$curPos]->[0]->[1] $cache[$curPos]->[0]->[2] DMV: $dmv $dMax $dT $dP\n";
	   #move RNA counts to this current position
	   for(my $reposi=22;$reposi<=35;$reposi++){
	     $cache[$curPos]->[0]->[$reposi] = $cache[$curPos+$ici]->[0]->[$reposi];
	     $cache[$curPos+$ici]->[0]->[$reposi] = "NA";
           } 
	   #POOL update
           for(my $reposi=36;$reposi<=49;$reposi++){
	     if($reposi == 41){
		#majNonRef is kept constant;
                $cache[$curPos+$ici]->[0]->[$reposi]  = "NA";
		next;
	     }
	     if($reposi == 44){ #update POOL_str
		my $npc =  $cache[$curPos]->[0]->[12] + $cache[$curPos]->[0]->[26] ;
		$cache[$curPos]->[0]->[$reposi] =~ s/$dmv\[\d+?\]/$dmv\[$npc\]/;
                $cache[$curPos+$ici]->[0]->[$reposi]  = "NA";
		next;
	     }
             $cache[$curPos]->[0]->[$reposi] += $cache[$curPos+$ici]->[0]->[$reposi];
             $cache[$curPos+$ici]->[0]->[$reposi] = "NA";
           }

	}elsif(
	 $dP != 0 &&
	 $rP ==  0 &&
	 $rT eq $dT
	){
           my $ici = $dP;
	   print STDERR "#dna $dT $dP shift to rna $cache[$curPos]->[0]->[1] $cache[$curPos]->[0]->[2] RMV: $rmv $rMax $rT $rP\n";
	   #move DNA counts to this current position
           for(my $reposi=8;$reposi<=21;$reposi++){
             $cache[$curPos]->[0]->[$reposi] = $cache[$curPos+$ici]->[0]->[$reposi];
             $cache[$curPos+$ici]->[0]->[$reposi] = "NA";
           }
	   for(my $reposi=34;$reposi<=49;$reposi++){ #clear pool.
             if($reposi == 41){ #majNonRef is kept constant;
               $cache[$curPos+$ici]->[0]->[$reposi]  = "NA";
	       next;
             }
             if($reposi == 44){ #update POOL_str
                my $npc =  $cache[$curPos]->[0]->[12] + $cache[$curPos]->[0]->[26] ;
                $cache[$curPos]->[0]->[$reposi] =~ s/$rmv\[\d+?\]/$rmv\[$npc\]/;
                $cache[$curPos+$ici]->[0]->[$reposi]  = "NA";
                next;
             }
             $cache[$curPos]->[0]->[$reposi] += $cache[$curPos+$ici]->[0]->[$reposi];
             $cache[$curPos+$ici]->[0]->[$reposi] = "NA";
           }
	}
	#last case is future position, in which current is left alone.
  }

  }
			
  #add high quality data filter.
  $cache[$curPos] = hqDataFlag($cache[$curPos]); #add hq data flags

  #add indel marker to line of current position.
  $l = join(",",@{$cache[$curPos]->[0]}).",".$cache[$curPos]->[1].",".$cache[$curPos]->[2].",".$cache[$curPos]->[3];
  
  if($fitted==0 && ($dhqCnt <= $trainNum || ( $rhqCnt <= $trainNum && $dnaOnly!=1) ) ){

    #print out initial data
    printOut("#INIT:".$l."\n")

  }elsif($fitted==0){
    #fit distribution to initial data
    fitDist();
  
    $fitted=1;    
    return(2);
  }else{
 
    $cache[$curPos] = evalDat($cache[$curPos]); #per row.    
    $l .= ",".$cache[$curPos]->[4].",".$cache[$curPos]->[5].",".$cache[$curPos]->[6]; #add p-values.

    if( 
	($cache[$curPos]->[4] ne "NA" &&  $cache[$curPos]->[4]  < $pvCut) ||
	($cache[$curPos]->[5] ne "NA" &&  $cache[$curPos]->[5]  < $pvCut) ||
	$verboseOut ){ #only print sites with above significance unless verbose out.
	    printOut($l."\n");
    }

  }
  $pln++;

  #remove first element, keeps at maximum $indelShadow size.
  if(scalar(@cache)>$indelShadow){
    shift @cache;
  }
}


sub printAllCache{
  for(my $curPos=0;$curPos<=$#cache;$curPos++){
#    my $l = $cache[$curPos]->[0].",".$cache[$curPos]->[1]."\n";
    my $l = join(",",@{$cache[$curPos]->[0]}).",".$cache[$curPos]->[1]."\n";
    printOut($l);
  }
}

### @@@ end check
sub lineComp {
  my $a = shift;
  my $b = shift;
  if (($a->[0] == $b->[0] && $a->[1]< $b->[1]) || $a->[0] < $b->[0]){
    return -1; #a behind b
  }
  if (($b->[0] == $a->[0] && $b->[1]< $a->[1]) || $b->[0] < $a->[0]){
    return 1; #b behind a
  }
  if ($a->[0] == $b->[0] && $a->[1] == $b->[1] ){
    return 2;
  }

}


sub openSNP{
	if($snpFile ne "-1" && $snpFile ne "blank"){
		open(SNP,$snpFile) || die "cannot open SNP $snpFile.";
		my $tmph = <SNP>;
		@snpf = split /\t/, <SNP>;

	}
}



####################
# end functions
####################


open(OF,">${resultsDir}/unceqr_proc.all.csv") || die "cannot access ${resultsDir}/unceqr_proc.all file";



###
# print all options to LOG.
###
	printOut($header,1);
	printOut("# \n# \n",1);
	printOut(("#"x20)."  \n",1);
	printOut("# \n# \n",1);
	printOut("#startDate: ".`date`,1);

	printOut("#tumorRNA $rnaf\n",1);
	printOut("#tumorDNA $tdnaf\n",1);
	printOut("#normalDNA $ndnaf\n",1);
	printOut("#normalDNAchr $ndnafChr\n",1);
	printOut("#tumorRNAchr $rnafChr\n",1);
	printOut("#tumorDNAchr $tdnafChr\n",1);
	printOut("#fastaWchr $fastaWchr\n",1);
	printOut("#fastaWoChr $fastaWoChr\n",1);
	printOut("#resultsDir $resultsDir\n",1);
	printOut("#regionsToQuery $posf\n",1);
	printOut("#dense $dense\n",1);
	printOut("#mainProc $mainProc\n",1);
	printOut("#snpFile $snpFile\n",1);
	printOut("#st $st\n",1);
	printOut("#trimEnd $trimEnd\n",1);
	printOut("#minTMapQ $minTMapQ\n",1);
	printOut("#minTBaseQ $minTBaseQ\n",1);
	printOut("#flagFilter $flagFilter\n",1);
	printOut("#maxIH $maxIH\n",1);
	printOut("#maxNM $maxNM\n",1);
	printOut("#maxDepth $maxDepth\n",1);
	printOut("#maxHomopolymer $maxHomopolymer\n",1);	
	printOut("#normalLoose $normalLoose\n",1);
	printOut("#preMinNormCov $preMinNormCov\n",1);
	printOut("#preTumorCov $preMinTumorCov\n",1);
	printOut("#minInNormFrac $minInNormFrac\n",1);
	printOut("#regionRestrict $regionRestrict\n",1);
	printOut("#minNormCnt $minNormCnt\n",1);
	printOut("#maxNormPlural $maxNormPlural\n",1);
	printOut("#minDnaTumorCnt $minDnaTumorCnt\n",1);
	printOut("#maxDnaBias $maxDnaBias\n",1);
	printOut("#maxDnaStruckProp $maxDnaStruckProp\n",1);
	printOut("#minRnaTumorCnt $minRnaTumorCnt\n",1);
	printOut("#maxRnaBias: $maxRnaBias\n",1);
	printOut("#maxRnaStruckProp: $maxRnaStruckProp\n",1);
	printOut("#maxDnaTumorPluralProp: $maxDnaTumorPluralProp\n",1);
	printOut("#maxRnaTumorPluralProp: $maxRnaTumorPluralProp\n",1);
	printOut("#trainNum: $trainNum\n",1);
	printOut("#dnaOnly: $dnaOnly\n",1);
	printOut("#rnaOnly: $rnaOnly\n",1);
	printOut("#priorMutFile: $maf\n",1);
        printOut("#endMutFile: $mafEnd\n",1);
	printOut("#indelShadow: $cacheLen\n",1);
	printOut("#indel merge distance: $posDist\n",1);
	printOut("#normIndelShadowFrac: $normIndelShadowFrac\n",1);
	printOut("#medStart: $medStart\n",1);
	printOut("#sampleId: $sampId\n",1);
	printOut("#verboseOut: $verboseOut\n",1);
	printOut("#pvCut: $pvCut\n",1);
	printOut("#unceqrver: $unceqrver\n",1);
	printOut("#mode: $mode\n",1);
	printOut("#\n#\n",1);



my $samtoolsGenome = 0;



my @cn = ("chr","pos","universal","normCnt","normPlural","normChars","normAllChars","norm_str",
"DNA_refCnt", "DNA_refCntNoIndel",  "DNA_nonRefCnt", "DNA_nonRefCntNoIndel","DNA_maxVal", "DNA_majNonRef", "DNA_nonRefPlural","DNA_str", "DNA_refPlus", "DNA_refNeg",  "DNA_varPlus", "DNA_varNeg", "DNA_struckProp","DNA_medDist",
"RNA_refCnt", "RNA_refCntNoIndel",  "RNA_nonRefCnt", "RNA_nonRefCntNoIndel","RNA_maxVal","RNA_majNonRef", "RNA_nonRefPlural","RNA_str",  "RNA_refPlus", "RNA_refNeg",  "RNA_varPlus", "RNA_varNeg","RNA_struckProp","RNA_medDist",
"POOL_refCnt", "POOL_refCntNoIndel",  "POOL_nonRefCnt", "POOL_nonRefCntNoIndel","POOL_maxVal", "POOL_majNonRef", "POOL_nonRefPlural","POOL_str",  "POOL_refPlus", "POOL_refNeg",  "POOL_varPlus", "POOL_varNeg", "POOL_struckProp","POOL_medDist",
"isSNP","featureId","indelShadow","DNA_hq","RNA_hq","p.value.dna","p.value.rna","p.value.meta");


for(my $i=0;$i<scalar(@cn);$i++){
  $revCol{$cn[$i]} = $i;
}



printOut( join(",",@cn)."\n" );

my $regArg;

my $goptRna = ($rnafChr eq "Y") ? " -f $fastaWchr " : " -f $fastaWoChr ";
my $goptTDna= ($tdnafChr eq "Y") ? " -f $fastaWchr " : " -f $fastaWoChr ";
my $goptNDna= ($ndnafChr eq "Y") ? " -f $fastaWchr " : " -f $fastaWoChr ";

my $filterOpt= " -B -A -q $minTMapQ -Q $minTBaseQ -T $trimEnd -H $maxIH -N $maxNM -d $maxDepth -K $flagFilter -Z $maxHomopolymer ";

#note: B disable BAQ computaton, A print anomalous pairs.


my $nfilterOpt = ($normalLoose eq "Y") ? " -B -A -q 0 -Q 0 -T 0 -H 10000 -N 1000 -d $maxDepth -K 0 -Z 100 "   : $filterOpt;

openSNP();

###
#filtered regions; regions to be filtered after the fact; must sorted in chr, pos order.
##

my (%chrMin,%chrMax);



	my @echrs;
	my $curre=0;
	if($posf ne ""){
	  @echrs=();
	  open(PF,$posf);
	 # my $tmp = <PF>; #header line
	  while(<PF>){
	    chomp;
	    my @f = split /\t/;
	    $f[0] =~ s/chr//;
	    $f[1]++; #because BED is zero based and SAM is one based.
	    $f[2];

	    if( !exists($chrMin{$f[0]}) ){
	    $chrMin{$f[0]} = $f[1] ;
	    $chrMax{$f[0]} = $f[2]  ;
	    }else{
	    $chrMin{$f[0]} = ($chrMin{$f[0]} > $f[1]) ? $f[1]:$chrMin{$f[0]}  ;
	    $chrMax{$f[0]} = ($chrMax{$f[0]} < $f[2]) ? $f[2]:$chrMax{$f[0]}  ;
	    }
	    @echrs[scalar(@echrs)] = \@f;
	  }
	  close PF;
	  printOut "#positions read: ".scalar(@echrs)."\n";
	}




#loop
my $regCmd;
if ($dense eq "Y"){
#  $regArg = " -l $posf ";

  #fitting:
  fl: foreach my $chrl ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y") { #
   if(!exists($chrMin{$chrl}) || !exists($chrMax{$chrl})){
     next;
   }
    printOut "#chrf $chrl\n";
    printOut "#".`date`."";
    $regCmd = " grep -P '^$chrl\\t' ";
    $regArg = " -r $chrl:$chrMin{$chrl}-$chrMax{$chrl} "; #this sets both region and BED; region takes precedence; reason for this is to terminate the region search after chromo
    my $ret = queryBam("f");
    if ($ret == 2){
      last fl;
    }
  }

  $curre=0; #reset
  close SNP;
  openSNP();
  foreach my $chrl ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y") { #
   if(!exists($chrMin{$chrl}) || !exists($chrMax{$chrl})){
     next;
   }
    printOut "#chrp $chrl\n";
    $regCmd = " grep -P '^$chrl\\t' ";
    $regArg = " -r $chrl:$chrMin{$chrl}-$chrMax{$chrl} "; #this sets both region and BED; region takes precedence; reason for this is to terminate the region search after chromo
    my $ret = queryBam("p");
  }




}elsif($maf ne "" && $maf ne "-1"){
  open(PF,$maf ) || die "cannot open $maf ";
  my $pcnt;
  PFL: while(<PF>){
    chomp;
    my @f = split /\t/;
    if ($f[15] !~ /$sampId/ ){next;} #lines with sample idenfitier anywhere in line are taken.
    $f[5] = $f[5]-10; #one based extra padding
    $f[6] = $f[6]+10; #one based extra padding
    $f[5] = ($f[5] < 1 ) ? 1 : $f[5]; 
    $regArg = " -r $f[4]:$f[5]-$f[6] ";
    $regCmd = " cat ";
    my $ret = queryBam("all");
    if ($ret == 2){
      last PFL;
    }
  }
  close PF;
  $curre=0; #reset
  close SNP;
  openSNP();
  if($fitted==0){
    fitDist();
    $fitted=1;
  }

  open(PF,$maf ) || die "cannot open $maf ";
  my $pcnt;
  PFL: while(<PF>){
    chomp;
    my @f = split /\t/;
    if ($f[15] !~ /$sampId/ ){next;} #lines with sample idenfifier anywhere in line are taken.
    $f[5] = $f[5]-10; #one based extra padding
    $f[6] = $f[6]+10; #one based extra padding
    $f[5] = ($f[5] < 1 ) ? 1 : $f[5]; 
    $regArg = " -r $f[4]:$f[5]-$f[6] ";
    $regCmd = " cat ";
    my $ret = queryBam("all");
  }
  close PF;


}else{ #non-dense positions
  open(PF,$posf) || die "cannot open $posf";
  my $pcnt=0;
  DFL: while(<PF>){
    chomp;
    my @f = split /\t/;
    $f[1]++; #bed is zero based;
    $f[1] -= $regPadd;
    $f[2] += $regPadd;
    $regArg = " -r $f[0]:$f[1]-$f[2] ";
    $regCmd = " cat ";
    my $ret = queryBam("all");
    if($ret == 2){
	last DFL;
    }
  }
  close PF;

  open(PF,$posf) || die "cannot open $posf";
  my $pcnt=0;
  DFL: while(<PF>){
    chomp;
    my @f = split /\t/;
    $f[1]++; #bed is zero based;
    $f[1] -= $regPadd;
    $f[2] += $regPadd;
    $regArg = " -r $f[0]:$f[1]-$f[2] ";
    $regCmd = " cat ";
    my $ret = queryBam("all");
  }
  close PF;
}


sub queryBam{
	my ($sampleb) = @_; #sample base position;

	my ($tcmd,$ncmd,$rcmd);
	#regArg comes with chr by default

	my $sedc;
	my $tregArg = $regArg;
	my $oregArg = $regArg;


	#prepare region control for mpileup: region and list; region is bound if both.
	if($rnafChr eq "Y"){
		$sedc = " sed 's/^/chr/' "; #prepend chr
		$regArg =~ s/\-r /\-r chr/;
	}else{
		$sedc = " cat ";
	}
        $sedc .= ($dense eq "Y") ? ' | perl -e "while(<>){ my (\\$chr,\\$p1,\\$p2) = \\$_ =~ /(\\w+)\\t(\\d+)\\t+(\\d+)/; \\$regPadd = '.$regPadd.';  \\$p1 -= \\$regPadd; \\$p2 += \\$regPadd; print \\"\\$chr\\t\\$p1\\t\\$p2\\n\\";   }"' : "";

	$rcmd = ($dense eq "Y") ? "$regCmd $posf | $sedc | $st mpileup $filterOpt $goptRna $regArg -l - $rnaf": "$st mpileup $filterOpt $regArg $goptRna $rnaf ";
	$rcmd  = ($rnaf eq "blank") ? "cat /dev/null" : $rcmd;
	printOut("#RNAcmd $rcmd\n");

	if($rnaf ne "blank"){
		open (RF, "-|", "$rcmd") || die "cannot open rna";
	}


	if($tdnafChr eq "Y"){
		$sedc = " sed 's/^/chr/' "; #prepend chr
		$tregArg =~ s/\-r /\-r chr/;
	}else{
		$sedc = " cat ";
	}
        $sedc .= ($dense eq "Y") ? ' | perl -e "while(<>){ my (\\$chr,\\$p1,\\$p2) = \\$_ =~ /(\\w+)\\t(\\d+)\\t+(\\d+)/; \\$regPadd = '.$regPadd.';  \\$p1 -= \\$regPadd; \\$p2 += \\$regPadd; print \\"\\$chr\\t\\$p1\\t\\$p2\\n\\";   }"' : "";

        $tcmd = ($dense eq "Y") ? "$regCmd $posf | $sedc | $st mpileup $filterOpt $goptTDna $tregArg -l - $tdnaf": "$st mpileup $filterOpt $tregArg $goptTDna $tdnaf ";
	$tcmd  = ($tdnaf eq "blank") ? "cat /dev/null" : $tcmd;
	printOut("#TCMD: $tcmd\n");
	if($tdnaf ne "blank" ){
	  open (TF, "-|", "$tcmd") || die "cannot open TF";
	}

	#normal and tumor come from same for now.
	$tregArg=$oregArg;
	if($ndnafChr eq "Y"){
		$sedc = " sed 's/^/chr/' "; #prepend chr
		$tregArg =~ s/\-r /\-r chr/;
	}else{
		$sedc = " cat ";
	}
        $sedc .= ($dense eq "Y") ? ' | perl -e "while(<>){ my (\\$chr,\\$p1,\\$p2) = \\$_ =~ /(\\w+)\\t(\\d+)\\t+(\\d+)/; \\$regPadd = '.$regPadd.';  \\$p1 -= \\$regPadd; \\$p2 += \\$regPadd; print \\"\\$chr\\t\\$p1\\t\\$p2\\n\\";   }"' : "";


	$ncmd = ($dense eq "Y") ? "$regCmd $posf | $sedc | $st mpileup $nfilterOpt $goptNDna $tregArg -l - $ndnaf": "$st mpileup $nfilterOpt $tregArg $goptNDna $ndnaf ";
	$ncmd  = ($ndnaf eq "blank") ? "cat /dev/null" : $ncmd;
	printOut("#normalCMD $ncmd\n");
	#die;
	if($ndnaf ne "blank"){
	  open (NF, "-|", $ncmd) || die "cannot open normal";
	}

	###
	#PROCESS LINES
	###
	my ($rline,$nline,$tline)=(0,0,0);
	my ($rinc,$ninc,$tinc)=(1,1,1);
	#my ($reof,$neof,$teof)=0;
	my @nRef;

	my $cond=1;

	my @rf = (10**10**10,10**10**10);
	my @nf = (10**10**10,10**10**10);
	my @tf = (10**10**10,10**10**10);

	POSW: while(1){
		
	if($rinc && $rnaf != -1 && !eof(RF)  ){
	  $rinc=0;
	  $rline=<RF>;
	  if(!eof(RF)){
  	    chomp $rline;
	    @rf = split /\t/, $rline;
	    $ry=1;
	  }else{
  	    $rf[0] = $rf[1] = 10**10**10;

	    
  	  }
	}
	
	if($tinc && $tdnaf != -1 && !eof(TF) ){
	  $tinc=0;
	  $tline=<TF>;
 	  if(!eof(TF)){
  	    chomp $tline;
  	    @tf = split /\t/, $tline;
	    $dy=1;
  	  }else{
  	    $tf[0] = $tf[1] = 10**10**10;

	    
  	  }
	}	
	
	if($ninc && $ndnaf != -1 && !eof(NF)  ){
	  $ninc=0;
	  $nline=<NF>;
  	  if(!eof(NF)){
  	    chomp $nline;
  	    @nf = split /\t/, $nline;
	    $chri = $nf[0];
	    $pos = $nf[1];
            @nRef = getRef(@nf);
 	    $ny=1;
  	  }else{
  	    $nf[0] = $nf[1] = 10**10**10;
	    @nRef=();

	    

  	  }
	}	
	

	#useful test for line sync, printOut "#lines $rf[0],$rf[1],$rinc ".eof(RF)."\t$tf[0],$tf[1],$tinc ".eof(TF)."\t$nf[0],$nf[1],$ninc ".eof(NF)."\n";

	if(eof(RF) && eof(NF) && eof(TF)){
	  last;
	}	
	
	###
	#compare fields based on common lines, increment lagging lines forward
	###

	my @tDnaCnt;
	my @tRnaCnt;
	my @poolCnt;

	my $ur;

	my $isSNP=0;
	
	$rf[0] =~ s/chr//;
	if($regChrLess == "N"){
		$tf[0] =~ s/chr//;
		$nf[0] =~ s/chr//;
	}
	

	if ( lineComp(\@rf,\@tf)==2 && lineComp(\@rf,\@nf)==2 ){    
            @nRef = getRef(@nf);
	    @tDnaCnt = countRef(\@tf, \@nRef);
	    @tRnaCnt = countRef(\@rf, \@nRef);
	    @poolCnt = countRef(mergeArr(\@tf,\@rf),\@nRef);
	    $pos=$rf[1];
	    $chri=$rf[0];
	    $tinc=$rinc=$ninc=1;
	    #$ur = getRefFASTA($chri,$pos);
	    $ur=$rf[2];
	    $ur =~ tr/[atcg]/[ATCG]/;
 	  
	#SINGLE CASES BELOW
	}elsif( lineComp(\@nf,\@rf)==-1 && lineComp(\@nf,\@tf)==-1){ #just normal
	  $ninc=1;	
  	  $pos=$nf[1];
	  $chri=$nf[0];
	  $ur = $nf[2];
	  $ur =~ tr/[atcg]/[ATCG]/; 	  

	}elsif( lineComp(\@rf,\@nf)==-1 && lineComp(\@rf,\@tf)==-1  ){ #just RNA
	  $pos=$rf[1];
	  $chri=$rf[0];
	  #$ur = getRefFASTA($chri,$pos);
	  $ur=$rf[2];
	  $ur =~ tr/[atcg]/[ATCG]/;
	  my %tmp; $tmp{$ur}=-1;
	  addSNP(\%tmp);
	  @nRef = ('',\%tmp); # universal ref
	  @poolCnt = @tRnaCnt = countRef(\@rf, \@nRef);
	  $rinc=1;	

	}elsif( lineComp(\@tf,\@rf)==-1 && lineComp(\@tf,\@nf)==-1 ){ #just tumor
	  $pos=$tf[1];
	  $chri=$tf[0];
	  $ur = $tf[2];
	  $ur =~ tr/[atcg]/[ATCG]/;
	  my %tmp; $tmp{$ur}=-1;
	  addSNP(\%tmp);
	  @nRef = ('',\%tmp); # universal ref
          @poolCnt = @tDnaCnt = countRef(\@tf, \@nRef);		  
	  $tinc=1;

	#PAIR CASES BELOW
	}elsif(lineComp(\@nf,\@rf)==-1 && lineComp(\@nf,\@tf)==2 ){ #tumor and normal
          @nRef = getRef(@nf);
          @poolCnt = @tDnaCnt = countRef(\@tf, \@nRef);	
	  $ninc=1;
	  $tinc=1;
    	  $pos=$nf[1];
	  $chri=$nf[0];
	  $ur = $nf[2];
          $ur =~ tr/[atcg]/[ATCG]/;  

	}elsif( lineComp(\@nf,\@rf)==2 && lineComp(\@nf,\@tf)==-1 ){ #normal and RNA
          @nRef = getRef(@nf);
	  @poolCnt = @tRnaCnt = countRef(\@rf, \@nRef);
	  $ninc=1;
	  $rinc=1;
          $pos=$nf[1]; 
	  $chri=$nf[0];
	  $ur = $nf[2];
    	  $ur =~ tr/[atcg]/[ATCG]/;

	}elsif(lineComp(\@rf,\@nf)==-1 && lineComp(\@rf,\@tf)==2){ #tumor and RNA, no normal.
 	  $pos=$rf[1];
	  $chri=$rf[0];
	  $ur = $rf[2];
	  $ur =~ tr/[atcg]/[ATCG]/;
	  my %tmp; $tmp{$ur}++;
	  addSNP(\%tmp);
	  @nRef = ('',\%tmp); # universal ref
          @tDnaCnt = countRef(\@tf, \@nRef);	
          @tRnaCnt = countRef(\@rf, \@nRef);
	  @poolCnt = countRef(mergeArr(\@tf,\@rf),\@nRef);
	  $tinc=1;
	  $rinc=1;

	}else{
	    die "trapped. case not specified. @nRef";
	}
	


	if ($sampleb eq "f"){
		#base sample
		my $r = rand();
		if ($r > 0.25  ){ #(1.414257e-08 * 10) ){
		  next POSW;
		}
	}	




	###
	# Region restriction.
	###

	#this is the region list for whole chromosome search
	if($regionRestrict eq "Y"){
		my $chriNoChr = $chri;
		$chri =~ s/chr//;
		#fast forward regions to CHR or to LP; assumes chr are in order
		while(
				($echrs[$curre]->[0] ne $chriNoChr) ||
				($echrs[$curre]->[0] eq $chriNoChr && $echrs[$curre]->[2] + 1 + $regPadd < $pos ) #pass region if right is less than position, not in region yet.
				&& $curre<scalar(@echrs)-1 ){  

			#test printOut "#advance reg: $chriNoChr $echrs[$curre]->[0] $pos $echrs[$curre]->[1] $echrs[$curre]->[2]\n";
			$curre++;	
		}

		#if regions are out, have advanced beyond the size of the array.
		if($curre == scalar(@echrs) ){
			last POSW;
			
		}

		 #fast forward ALN if less than left position of region.
		if($echrs[$curre]->[0] eq $chriNoChr && $pos < $echrs[$curre]->[1] - 1 - $regPadd ){  #should be lt not lte
			#test  printOut "#advance aln: $chriNoChr $echrs[$curre]->[0] $pos $echrs[$curre]->[1] \n";
		          next POSW;
		}
	}	

	    my $normChari = join("-",keys(%{$nRef[1]}));
	if( ($tDnaCnt[0] + $tDnaCnt[2]  >=$preMinTumorCov || $tRnaCnt[0] + $tRnaCnt[2] >= $preMinTumorCov || $poolCnt[0] + $poolCnt[1] >= $preMinTumorCov || $normChari =~ /ins|del/)
		&& $nRef[0] >= $preMinNormCov){ #only print if one of data types has coverage.
	    my $ot;
	    $ot .= "$chri,$pos,$ur,";
	    #my $normAchar = join("-",keys(%{$nRef[3]}));

	    $ot .= join(",",($nRef[0],$nRef[2],$normChari,$nRef[3],$nRef[4])).",";
	    $ot .= join(",",(@tDnaCnt[0..12],$tDnaCnt[15]) ).",";
            $ot .= join(",",(@tRnaCnt[0..12],$tRnaCnt[15]) ).",";
	    $ot .= join(",",(@poolCnt[0..12],$poolCnt[15]) ).",";

	    #SNP pos
	    $ot .= "$isSNP";

	    #region
	    $ot .=  ",$curre\n"; #featureId
	    my $rt;
	    if($cacheLen>0){
		    $rt = printOutCache($ot,$nRef[7]); #this includes normIndelFrac
	    }else{
		    $ot =~ s/\n//;
		    $rt = printOutCache($ot,"0"); #zero for indelShadow		
	    }
       		 #test for line sync, printOut "#synch: $rf[0],$rf[1],$rinc ".eof(RF)."\t$tf[0],$tf[1],$tinc ".eof(TF)."\t$nf[0],$nf[1],$ninc ".eof(NF)."\n";
        
	   if($rt ==2){
		return(2);
	   }
	 }
	}
}
if($cacheLen>0){
printOut("empty Cache\n");
printAllCache();
}

printOut("#dtn: $dtn");
printOut("#rtn: $rtn");
printOut("#mtn: $mtn");


printOut("#unceqr_proc fileMD5SUM: .".`md5sum $0`);
printOut("#unceqr_proc pwd: .".`pwd`);
printOut("#unceqr_proc uname: .".`uname -a`);
printOut("#unceqr_proc finish date".`date`);

close SNP;
close OF;

1;

