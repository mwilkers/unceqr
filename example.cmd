perl UNCeqR_main_simple.pl -tumorDNA samtools-0.1.18_unceqr/examples/ex1.bam -tumorRNA samtools-0.1.18_unceqr/examples/ex1.bam -sampleId test -resultsDir ./ -normalDNAchr N -tumorDNAchr N -tumorRNAchr N -fastaWoChr samtools-0.1.18_unceqr/examples/ex1.fa -regionsToQuery example_reg.txt -verboseOut 0 -dense N -minNormCnt 0 -preMinNormCov 0 -minRnaTumorCnt 1 -preMinTumorCov 1 -minDnaTumorCnt 1 -trainNum 500 -mode denovo -indelShadow 0