#DESCRIPTION: transforms unceqr_test.csv output to VCF format.

print STDERR "#NOTE: UNCeqR_2vcf.pl transforms each line of unceqr_test.csv to vcf rows.  There is the rare possibility, after filtering for significance, that two indels may overlap the same region and could affect compatible representation in the VCF file.  This is not supported and users interested in this should edit/curate their VCF as desired.\n";


my $fn = shift;
my $cut = shift;

open(OF,$fn) || die "cannot open $fn\n";

my $head;
W: while(<OF>){
 $head = $_;
 if( $_ !~ /^#/){
   last W;
 }
}

chomp $head;
my @cn = split /\,/, $head;
my %revCol;
for(my $i=0;$i<scalar(@cn);$i++){
  $revCol{$cn[$i]} = $i;
}



my $l=0;

open(VCFOUT,">$fn.vcf") || die "cannot open $fn.vcf\n";

print VCFOUT join("\t",("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"))."\n";
OFL: while(<OF>){
  if(/^#/){
    next;
  }
  chomp;
  my @f = split /\,/;

  #limit to significant calls from any model.
  if(!( 

	($f[$revCol{"p.value.dna"}] ne "NA" && $f[$revCol{"p.value.dna"}]<$cut) || 
	($f[$revCol{"p.value.rna"}] ne "NA" && $f[$revCol{"p.value.rna"}]<$cut) || 
	($f[$revCol{"p.value.meta"}] ne "NA" && $f[$revCol{"p.value.meta"}]<$cut) )){
    next;
  }

  my $refAllele = $f[$revCol{"universal"}]; #use reference allele.  not germline.
#  my $varAllele = $f[41]; #POOL majNonRef.

  my $varAllele = ($f[$revCol{"p.value.dna"}] ne "NA" && ($f[$revCol{"p.value.rna"}] eq "NA" || $f[$revCol{"p.value.dna"}] <= $f[$revCol{"p.value.rna"}]) ) ? $f[$revCol{"DNA_majNonRef"}] : $f[$revCol{"RNA_majNonRef"}]; #take alt allele based on significance.


  if($varAllele =~ /ins/){
    $varAllele =~ s/ins//;
    $varAllele = "$refAllele$varAllele";
  }

  if($varAllele =~ /del/){
    my ($tmp) = $varAllele =~ /del(.+)/;
    $varAllele=$refAllele;
    $refAllele .= $tmp;
  }

  $refAllele =~ tr/a-z/A-Z/;
  $varAllele =~ tr/a-z/A-Z/;

  print VCFOUT join("\t",($f[0],$f[1],".",$refAllele,$varAllele,".",".","."))."\n";
  $l++;

}

close VCFOUT;
close OF;


exit();






