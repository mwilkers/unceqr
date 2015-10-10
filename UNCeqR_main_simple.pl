use Getopt::Long qw(GetOptionsFromArray);
Getopt::Long::Configure("pass_through");

use strict;
use vars qw($rnaf $tdnaf $ndnaf $preMinTumorCov $preMinNormCov $posf $startChr $endChr $maf $sampId $tumorLoose $normalLoose $filterTumor $testReg $regionRestrict $minTMapQ $minTBaseQ $maxNM $maxIH $trimEnd $flagFilter $st $minInNormFrac $filterScript $maxDepth $header $mainProc $resultsDir $snpFile $fastaWchr $fastaWoChr $tdnafChr $ndnafChr $rnafChr $minNormCnt $maxNormPlural $minDnaTumorCnt $maxDnaBias $maxDnaStruckProp $minRnaTumorCnt $maxRnaBias $maxRnaStruckProp $trainNum $maxDnaTumorPluralProp $maxRnaTumorPluralProp $dnaOnly $rnaOnly $localPath $laneFile $maxHomopolymer $mafEnd $annFlag $verboseOut $R $Rscript $R_LIBS_USER $perlLib $annovarPath $analysisPath $medStart $indelShadow $sampleFile $curRecord $pvCut);

require 'UNCeqR_conf.pl';



my $rootDir = $sampId;
$resultsDir .= "/$rootDir";

`mkdir $resultsDir`;

open(STDERR,">${resultsDir}/unceqr_stderror") || die "cannot access ${resultsDir}/unceqr_stderror  log";
print STDERR "#starting UNCeqR_proc\n";

  require "UNCeqR_proc.pl";


#convert test.csv to VCF file
my $pv = 0.001;
my $rcmd = "perl -I$localPath $localPath/UNCeqR_2vcf.pl $resultsDir/unceqr_proc.all.csv $pv >& $resultsDir/log.vcf";
print `$rcmd`;


