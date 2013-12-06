#!/usr/bin/perl -w

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use Math::Random;
use Math::Round;
use Storable;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use List::Util qw (sum shuffle min max);
use threads;
use threads::shared;
use Statistics::Descriptive;
use URI::Escape;
use Cwd 'abs_path';
#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
#	Description
#		This is a Perl script finds the transcription start sites by comparing the read5End pileups of two types of libraries: TEX treated (TEX) and Standard (STD). 
#
#	Input
#
#		--TEXRead5EndPileupIndxPathAndSize=		path and integer[compulsory]; path of the TEX treated libraries index.hsh.pls and its lib size, in format of path,size; could be multiple;
#		--STDRead5EndPileupIndxPathAndSize=		path and integer [compulsory]; path of the standard libraries index.hsh.pls and its lib size, in format of path,size; could be multiple;
#		--TEXBaseCompositionPileupIndxPath=		path [compulsory]; path of the based composition pileup index.hsh.pls of the TEX treated libraries for checking start sites;
#		--fullReadDefineNATPileupIndxPath=		path [compulsory]; path of the full read pileup index.hsh.pls which will be used to define genes with or without NAT ;
#		--trnsfrgInfoHshStorablePath=			path [compulsory]; path of the trnsfrgInfoHsh generated in cisNATFinder;
#		--ratioCutoffPct=						percentage [95]; percentage cutoff for “ratio” used to distinguish exon and tss distributions
#		--countCutoffPct=						percentage [50]; percentage cutoff for “count” used to distinguish exon and tss distributions
#		--fastaPath=							file path [compulsory]; the path fasta file contains the genome sequence, for generating blank perl storables;
#		--gffPath=								path[compulsory]; path of the reference GFF for gene annotation;
#		--outDir=								output directory
#
#	v0.2
#		[Tue 17 Sep 2013 15:11:30 CEST] trnsfrgInfoHshStorablePath compatible now with transfragDiscoverer
#		[Tue 17 Sep 2013 16:52:47 CEST] removed TSSi analyses
#
#	v0.3
#		[Wed 18 Sep 2013 20:23:32 CEST] won’t print individual histogram, 
#
#	v0.4
#		[Sat 21 Sep 2013 14:27:33 CEST] will report a final noise index based on the ratio of finalized site in exon and putative TSS;
#		[Sat 21 Sep 2013 19:41:34 CEST] ratioCutoffPct and countCutoffPct move to input options
#
#	#-----finalPooleSTDTrack
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSFinder/v0.4/TSSFinder_v0.4.pl
#	--gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff
#	--TEXRead5EndPileupIndxPath=/Volumes/C_Analysis/NGS/results/E011_TEX5Map_basic_min35nt/E011_merged_TEX5Map_A_B_test_extTreat_with3_NM/mergedBam/BAMToReadEndPerlStorable/countMode.5.offset.0.baseComp.no/cntgCovPls/index.hsh.pls,1
#	--STDRead5EndPileupIndxPath=/Volumes/C_Analysis/NGS/results/E011_TEX5Map_basic_min35nt/E011_optTEX5Map_TEXminusTAPminus/mergedBam/BAMToReadEndPerlStorable/countMode.full.offset.0.baseComp.no/cntgCovPls/index.hsh.pls,1
#	--TEXBaseCompositionPileupIndxPath=/Volumes/C_Analysis/NGS/results/E011_TEX5Map_basic_min35nt/E011_merged_TEX5Map_A_B_test_extTreat_with3_NM/mergedBam/BAMToReadEndPerlStorable/countMode.5.offset.0.baseComp.yes/cntgCovPls/index.hsh.pls
#	--fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/fasta/genome.sorted.fa
#	--trnsfrgInfoHshStorablePath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/transfragDiscoverer/v0.6/transfragDiscoverer/storable/trnsfrgInfoHsh.pls
#	--fullReadDefineNATPileupIndxPath=/Volumes/C_Analysis/NGS/results/E011_TEX5Map_basic_min35nt/E011_optTEX5Map_TEXminusTAPminus/mergedBam/BAMToReadEndPerlStorable/countMode.full.offset.0.baseComp.no/cntgCovPls/index.hsh.pls
#	--ratioCutoffPct=95
#	--countCutoffPct=50
#	--outDir=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSFinder/v0.4/finalPooleSTDTrack/
#
#
#	#-----polIIPausing
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSFinder/v0.4/TSSFinder_v0.4.pl
#	--gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff
#	--TEXRead5EndPileupIndxPath=/Volumes/C_Analysis/NGS/results/E023_TSS35To50nt_basic_small/E023_TSS35To50nt_HM1_heatShockT0H/filterBamWith32NtLonger/BAMToReadEndPerlStorable/countMode.5.offset.0.baseComp.no/cntgCovPls/index.hsh.pls,1
#	--STDRead5EndPileupIndxPath=/Volumes/C_Analysis/NGS/results/E011_TEX5Map_basic_min35nt/E011_optTEX5Map_TEXminusTAPminus/mergedBam/BAMToReadEndPerlStorable/countMode.full.offset.0.baseComp.no/cntgCovPls/index.hsh.pls,1
#	--TEXBaseCompositionPileupIndxPath=/Volumes/C_Analysis/NGS/results/E023_TSS35To50nt_basic_small/E023_TSS35To50nt_HM1_heatShockT0H/filterBamWith32NtLonger/BAMToReadEndPerlStorable/countMode.5.offset.0.baseComp.yes/cntgCovPls/index.hsh.pls
#	--fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/fasta/genome.sorted.fa
#	--trnsfrgInfoHshStorablePath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/transfragDiscoverer/v0.6/transfragDiscoverer/storable/trnsfrgInfoHsh.pls
#	--fullReadDefineNATPileupIndxPath=/Volumes/C_Analysis/NGS/results/E011_TEX5Map_basic_min35nt/E011_optTEX5Map_TEXminusTAPminus/mergedBam/BAMToReadEndPerlStorable/countMode.full.offset.0.baseComp.no/cntgCovPls/index.hsh.pls
#	--outDir=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSFinder/v0.3/polIIPausing/
#
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	[2013-10-20 17:06]	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSFinder/v0.4/TSSFinder_v0.4.pl --gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff --TEXRead5EndPileupIndxPathAndSize=/Volumes/C_Analysis/NGS/results/E011_TEX5Map_basic_min35nt/E011_merged_TEX5Map_A_B_test_extTreat_with3_NM/mergedBam/BAMToReadEndPerlStorable/countMode.5.offset.0.baseComp.no/cntgCovPls/index.hsh.pls,32800863 --STDRead5EndPileupIndxPathAndSize=/Volumes/C_Analysis/NGS/results/E011_TEX5Map_basic_min35nt/E011_optTEX5Map_TEXminusTAPminus/mergedBam/BAMToReadEndPerlStorable/countMode.full.offset.0.baseComp.no/cntgCovPls/index.hsh.pls,32800863 --TEXBaseCompositionPileupIndxPath=/Volumes/C_Analysis/NGS/results/E011_TEX5Map_basic_min35nt/E011_merged_TEX5Map_A_B_test_extTreat_with3_NM/mergedBam/BAMToReadEndPerlStorable/countMode.5.offset.0.baseComp.yes/cntgCovPls/index.hsh.pls --fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/fasta/genome.sorted.fa --trnsfrgInfoHshStorablePath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/transfragDiscoverer/v0.6/transfragDiscoverer/storable/trnsfrgInfoHsh.pls --fullReadDefineNATPileupIndxPath=/Volumes/C_Analysis/NGS/results/E011_TEX5Map_basic_min35nt/E011_optTEX5Map_TEXminusTAPminus/mergedBam/BAMToReadEndPerlStorable/countMode.full.offset.0.baseComp.no/cntgCovPls/index.hsh.pls --ratioCutoffPct=99.9 --countCutoffPct=90 --outDir=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSFinder/v0.4/TEXTest/
#
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSFinder/v0.4/TSSFinder_v0.4.pl
#	--gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff
#	--TEXRead5EndPileupIndxPathAndSize=/Volumes/C_Analysis/NGS/results/E011_TEX5Map_basic_min35nt/E011_merged_TEX5Map_A_B_test_extTreat_with3_NM/mergedBam/BAMToReadEndPerlStorable/countMode.5.offset.0.baseComp.no/cntgCovPls/index.hsh.pls,32800863
#	--STDRead5EndPileupIndxPathAndSize=/Volumes/C_Analysis/NGS/results/E011_TEX5Map_basic_min35nt/E011_optTEX5Map_TEXminusTAPminus/mergedBam/BAMToReadEndPerlStorable/countMode.full.offset.0.baseComp.no/cntgCovPls/index.hsh.pls,32800863
#	--TEXBaseCompositionPileupIndxPath=/Volumes/C_Analysis/NGS/results/E011_TEX5Map_basic_min35nt/E011_merged_TEX5Map_A_B_test_extTreat_with3_NM/mergedBam/BAMToReadEndPerlStorable/countMode.5.offset.0.baseComp.yes/cntgCovPls/index.hsh.pls
#	--fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/fasta/genome.sorted.fa
#	--trnsfrgInfoHshStorablePath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/transfragDiscoverer/v0.6/transfragDiscoverer/storable/trnsfrgInfoHsh.pls
#	--fullReadDefineNATPileupIndxPath=/Volumes/C_Analysis/NGS/results/E011_TEX5Map_basic_min35nt/E011_optTEX5Map_TEXminusTAPminus/mergedBam/BAMToReadEndPerlStorable/countMode.full.offset.0.baseComp.no/cntgCovPls/index.hsh.pls
#	--ratioCutoffPct=99.9
#	--countCutoffPct=90
#	--outDir=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSFinder/v0.4/TEXTest/
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
my $scriptDirPath = dirname(rel2abs($0));
my $globalTmpLogPath = "$scriptDirPath/tmp.log.txt";
open TMPLOG, ">", $globalTmpLogPath;
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|3694, readParameters|4232
#	secondaryDependOnSub: currentTime|1071
#
#<section ID="startingTasks" num="0">
#-----print start log
&printCMDLogOrFinishMessage("CMDLog");#->3694

#-----read parameters
my ($TEXRead5EndPileupIndxPathAndSizeAry_ref, $STDRead5EndPileupIndxPathAndSizeAry_ref, $TEXBaseCompositionPileupIndxPath, $fullReadDefineNATPileupIndxPath, $trnsfrgInfoHshStorablePath, $fastaPath, $countCutoffPct, $ratioCutoffPct, $gffPath, $outDir) = &readParameters();#->4232
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineHardCodedParam
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedParam" num="1">
my $forcePoolLibRd5End = 'no';
my $forceRunGoldenmRNASetTSSi = 'no';
my $forceRunMemeAndDreme = 'no';
my $forceDefineEnd3NAT = 'yes';

my $copyIGVTrack = 'no'; #---copy all IGV tracks into one single folder, for easier migration of all tracks to another compute for local displays;
my $fontSize = 10; #---fontSize for IGV
my $maxPctDistBtwATGTrnsfrgEnd = 75; #---maximum percentile of mRNA sorted by distance between ATG and transfrag will be used as he golden set for bona fide coding mRNA TSS
my $exonEndTrim = 50; #----nt to be trim from exon fro defining a negetive set in TEX vs STD ratio training, to avoid possible TSS in exon;
my $minOriginalCount = 3; #----minium count in the original TEX ans STD to be considered as non-zero
my $goldenmRNASetSizeLimit = 10000; #-----limit the golden geneset size

my $TSSSearchRngHsh_ref = {};
$TSSSearchRngHsh_ref->{'ATG'}{'up'} = 150;
$TSSSearchRngHsh_ref->{'ATG'}{'dn'} = 150;
$TSSSearchRngHsh_ref->{'TAA'}{'up'} = 150;
$TSSSearchRngHsh_ref->{'TAA'}{'dn'} = 150;

my $motifEValueLimitHsh_ref = {};
$motifEValueLimitHsh_ref->{'meme'} = 1e-50;
$motifEValueLimitHsh_ref->{'dreme'} = 1e-100;
my $maxHitPVal = 0.01; #----maximum p-value as va valid MAST hit

#my $memeDremeToRunHsh_ref = {'meme'=>1};
#my $memeDremeToRunHsh_ref = {'dreme'=>1};
#my $memeDremeToRunHsh_ref = {'meme'=>1, 'dreme'=>1};
my $memeDremeToRunHsh_ref = {};

my $trimBaseComPlotRngHsh_ref = {};
$trimBaseComPlotRngHsh_ref->{'start'} = 0; #-----starting positon of the full seq in subtr in &plotBaseCompositionAroundmRNAReferencePoint
$trimBaseComPlotRngHsh_ref->{'length'} = 200; #-----length of full seq in subtr in &plotBaseCompositionAroundmRNAReferencePoint

#----fix the PutativeTSSReg in noise score calculation, for better comparison across different parameters. use “auto” to noy fix it
my $fixPutativeTSSRegHsh_ref = {
	'lowerLimit' => -20,
	'upperLimit' => 0,
};

#---range of the gene end to be counted measured by the percentile，appear in TSSATGSenseHist.pdf and TSSTAAAntisenseHist.pdf
my $validTSSLimitHsh_ref = {};
$validTSSLimitHsh_ref->{'ATG'}{'sense'}{'lowerPctLimit'} = 5;
$validTSSLimitHsh_ref->{'ATG'}{'sense'}{'upperPctLimit'} = 95;
$validTSSLimitHsh_ref->{'TAA'}{'antisense'}{'lowerPctLimit'} = 5;
$validTSSLimitHsh_ref->{'TAA'}{'antisense'}{'upperPctLimit'} = 95;

my $maxThread = 10; #---number of threads to be spawned

my $resultDirTag = "MC$minOriginalCount.PD$maxPctDistBtwATGTrnsfrgEnd.CR$ratioCutoffPct.CC$countCutoffPct";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_defineHardCodedPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedPath" num="2">
my $IGVGenomeID = "EHI_v3.0";
my $IGVTrackToCopyAry_ref = ();
my $hardCodedIGVPathHsh_ref = {};
$hardCodedIGVPathHsh_ref->{'GCPath'} = "/Volumes/A_MPro2TB/NGS/IGVInfo/infoTracks/EHI/GC.Repetitiveness/win.100_step.10.GC.tdf";
$hardCodedIGVPathHsh_ref->{'reptvPath'} = "/Volumes/A_MPro2TB/NGS/IGVInfo/infoTracks/EHI/GC.Repetitiveness/win.100_step.10.Reptv.tdf";
$hardCodedIGVPathHsh_ref->{'plusCovTDFPath'} = "/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/artefactCorrectedPileup/EHI_Standard/modelFit/W200.S100.I1.R0.10000000/modelFitData.CL0.99/correctedWig.upr/corrected.plus.wig.tdf";
$hardCodedIGVPathHsh_ref->{'minusCovTDFPath'} = "/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/artefactCorrectedPileup/EHI_Standard/modelFit/W200.S100.I1.R0.10000000/modelFitData.CL0.99/correctedWig.upr/corrected.minus.wig.tdf";
$hardCodedIGVPathHsh_ref->{'polyABamPath'} = "/Volumes/A_MPro2TB/NGS/analyses/EHI_polyA/polyATailDiscoverer/N.yes.D10.8.U10.8.T6.6.M2.B5.Y80.P70.L20.A5/BR5BA20BU60BD60R2P3A10G30C0B200S200K6.6P50.50.0.50/filterSam/all.filter.bam";
$hardCodedIGVPathHsh_ref->{'nonAnnotatedPolymerFreqPath'} = "/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/GFFHandle/GFFSeqExtractor/v0.8/EHI_v2_allFearues.forPileupCounter/polymerFreq/nonAnnotated.polymerFreq.txt";
$hardCodedIGVPathHsh_ref->{'gffPath'} = $gffPath;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_defineOutDirPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutDirPath" num="3">
my @mkDirAry;
my $libStorableDir = "$outDir/libData/storable/"; push @mkDirAry, $libStorableDir; 
my $libWigDir = "$outDir/libData/wig/"; push @mkDirAry, $libWigDir;
my $originalRd5EndStorableDir = "$libStorableDir/originalRd5EndPileup/"; push @mkDirAry, $originalRd5EndStorableDir;
my $scaledRd5EndStorableDir = "$libStorableDir/scaledRd5EndPileup/"; push @mkDirAry, $scaledRd5EndStorableDir;
my $resultStorableDir = "$outDir/result/$resultDirTag/storable/"; push @mkDirAry, $resultStorableDir;
my $resultLogDir = "$outDir/result/$resultDirTag/log/"; push @mkDirAry, $resultLogDir;
my $resultGFFDir = "$outDir/result/$resultDirTag/GFF/"; push @mkDirAry, $resultGFFDir;
my $resultWigDir = "$outDir/result/$resultDirTag/wig/"; push @mkDirAry, $resultWigDir;
my $resultXMLDir = "$outDir/result/$resultDirTag/XML/"; push @mkDirAry, $resultXMLDir;
my $resultDremeDir = "$outDir/result/$resultDirTag/dreme/"; push @mkDirAry, $resultDremeDir;
my $resultMemeDir = "$outDir/result/$resultDirTag/meme/"; push @mkDirAry, $resultMemeDir;
my $resultMastDir = "$outDir/result/$resultDirTag/mast/"; push @mkDirAry, $resultMastDir;
my $resultFastaDir = "$outDir/result/$resultDirTag/fasta/"; push @mkDirAry, $resultFastaDir;
my $resultIGVTrackCopyDir = "$outDir/result/$resultDirTag/IGVTrackCopy/"; push @mkDirAry, $resultIGVTrackCopyDir;

my $ggplotDirHsh_ref = {};
my @ggplotFileTypeAry = qw /dat pdf R log/;
foreach my $fileType (@ggplotFileTypeAry) {$ggplotDirHsh_ref->{$fileType} = "$outDir/result/$resultDirTag/ggplot/$fileType"; push @mkDirAry, $ggplotDirHsh_ref->{$fileType};}

my $weblogoDirHsh_ref = {};
my @weblogFileTypeAry = qw /pdf fasta/;
foreach my $fileType (@weblogFileTypeAry) {$weblogoDirHsh_ref->{$fileType} = "$outDir/result/$resultDirTag/weblogo/$fileType"; push @mkDirAry, $weblogoDirHsh_ref->{$fileType};}

my $filterTEXCovStorableDir = "$resultStorableDir/filterTEXCov/"; push @mkDirAry, $filterTEXCovStorableDir;

system ("mkdir -pm 777 $_") foreach @mkDirAry;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_defineOutFilePath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutFilePath" num="4">
my $XMLPath = "$resultXMLDir/TSSFinder.igv.xml";
my $trnsfrgGFFPathHsh_ref = {};
foreach my $plusOrMinus ('plus', 'minus') {
	$trnsfrgGFFPathHsh_ref->{$plusOrMinus} = "$resultGFFDir/trnsfrg.$plusOrMinus.gff";
}

my $wigglePathHsh_ref = {};
my @TSSItemAry = qw /TSS_goldenmRNA_TEX_original/;
foreach my $item (@TSSItemAry) {
	foreach my $plusOrMinus ('plus', 'minus') {
		$wigglePathHsh_ref->{$item}{$plusOrMinus} = "$resultWigDir/$item.$plusOrMinus.wig";
	}
}

foreach my $plusOrMinus ('plus', 'minus') {
	$wigglePathHsh_ref->{'filterTEX'}{$plusOrMinus} = "$resultWigDir/filterTEX.$plusOrMinus.wig";
}
my $TSSvsExonCutoffHshPath = "$resultStorableDir/TSSvsExonCutoffHsh.pls";
my $geneBasedTSSLogPath = "$resultLogDir/geneBasedTSSLog.xls";
my $mRNANATInfoHshPath = "$resultStorableDir/mRNANATInfoHsh.pls";
my $mRNANATInfoLogPath = "$resultLogDir/mRNANATInfoHsh.xls";
my $strndEndValidTSSInfoHshPath = "$resultStorableDir/strndEndValidTSSInfoHsh.pls";
my $geneBasedTSSInfoHshPath = "$resultStorableDir/geneBasedTSSInfoHsh.pls";
my $mRNARefPtHshPath = "$resultStorableDir/mRNARefPtHsh.pls";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 5_processInputData
#	primaryDependOnSub: checkGeneInfo|718, defineLibTypeInfo|1166, getAndPrintTrnsfrg|1672, getCtgryGeneInfo|1997, getIndivCntgCovPlsPath|2146, poolLibRd5EndCov|3584, readGFF_oneRNAPerGene|4075, readMultiFasta|4174, zipUnzipCntgCovInPlsPathHsh|4603
#	secondaryDependOnSub: createEmptyGenomeCovPerlStorable|1010, currentTime|1071, getIndivCntgCovPlsPath|2146, printGFF_oneRNAPerGene_chooseStrnd_filterAry|3727, storeAndScaleRd5EndFromBothLibraries|4452, warningToProceed|4583
#
#<section ID="processInputData" num="5">
#----------Define library information
my ($libTypeInfoHsh_ref) = &defineLibTypeInfo($STDRead5EndPileupIndxPathAndSizeAry_ref, $TEXRead5EndPileupIndxPathAndSizeAry_ref, $libWigDir);#->1166

#----------Read fasta
my ($fastaHsh_ref) = &readMultiFasta($fastaPath);#->4174

#----------Read Gff
my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);#->4075
&checkGeneInfo($geneInfoHsh_ref);#->718

#----------Get bfmRNA and mRNA ranges
my @mRNAAry = qw/bfmRNA/;
my ($mRNAInfoHsh_ref, $mRNAByCntgHsh_ref)= &getCtgryGeneInfo($geneInfoHsh_ref, \@mRNAAry);#->1997

#----------Retrieve the transfrags
my ($trnsfrgInfoHsh_ref) = &getAndPrintTrnsfrg($trnsfrgInfoHshStorablePath, $trnsfrgGFFPathHsh_ref);#->1672

#----------Pool 5End Read Cov from all libraries into one set of storable
my ($rd5EndPlsInfoHsh_ref) = &poolLibRd5EndCov($fastaHsh_ref, $originalRd5EndStorableDir, $scaledRd5EndStorableDir, $libTypeInfoHsh_ref, $forcePoolLibRd5End, $maxThread);#->3584

#----------Get the paths of the TEXBaseCompositionPileup
my ($TEXBaseCompositionPileupPathHsh_ref) = &getIndivCntgCovPlsPath($TEXBaseCompositionPileupIndxPath);#->2146

#----------Get the paths of the fullReadDefineNATPileup
my ($fullReadDefineNATPileupPathHsh_ref) = &getIndivCntgCovPlsPath($fullReadDefineNATPileupIndxPath);#->2146
&zipUnzipCntgCovInPlsPathHsh('unzip', $fullReadDefineNATPileupPathHsh_ref);#->4603
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 6_trainningSTDTEXRatio
#	primaryDependOnSub: trainingUsingExonVsTssData|4548
#	secondaryDependOnSub: analyzeSTDTEXRd5EndRatio|489, getRd5EndCountForExonAndTSSOfGoldenmRNASet|2451, predictBonaFideTSSinGoldenmRNASet|3628, reportStatus|4274
#
#<section ID="trainningSTDTEXRatio" num="6">
#---check overlapping between mRNA and trnsfrg, and define a set of gold standard genes, then find TSS
my ($TSSvsExonCutoffHsh_ref) = &trainingUsingExonVsTssData($forceRunGoldenmRNASetTSSi, $trnsfrgInfoHsh_ref, $mRNAInfoHsh_ref, $ggplotDirHsh_ref, $maxPctDistBtwATGTrnsfrgEnd, $libTypeInfoHsh_ref, $mRNAByCntgHsh_ref, $wigglePathHsh_ref, $goldenmRNASetSizeLimit, $exonEndTrim, $rd5EndPlsInfoHsh_ref, $minOriginalCount, $TSSvsExonCutoffHshPath, $ratioCutoffPct, $countCutoffPct, $maxThread);#->4548
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 7_filterTEXTrackAfterTrainning
#	primaryDependOnSub: filterTEXSiteBasedOnTrainingCutoff|1275, printWiggleSingleTrackFromCntgCovPlsPathHsh|3997
#	secondaryDependOnSub: checkRunningThreadAndWaitToJoin|981, createEmptyGenomeCovPerlStorable|1010, generateThreadHshWithRandomItem|1645, getIndivCntgCovPlsPath|2146, reportStatus|4274
#
#<section ID="filterTEXTrackAfterTrainning" num="7">
my ($filterTEXCovPlsPathHsh_ref)= &filterTEXSiteBasedOnTrainingCutoff($TSSvsExonCutoffHsh_ref, $rd5EndPlsInfoHsh_ref, $libTypeInfoHsh_ref, $fastaHsh_ref, $filterTEXCovStorableDir, $minOriginalCount, $maxThread);#->1275
&printWiggleSingleTrackFromCntgCovPlsPathHsh($filterTEXCovPlsPathHsh_ref, 0, $wigglePathHsh_ref->{'filterTEX'}{'plus'}, 'no');#->3997
&printWiggleSingleTrackFromCntgCovPlsPathHsh($filterTEXCovPlsPathHsh_ref, 1, $wigglePathHsh_ref->{'filterTEX'}{'minus'}, 'no');#->3997
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 8_investigateGeneTSSRelativePos
#	primaryDependOnSub: calculateNoiseIndex|625, getGeneWithValidTSSAtGeneEnd|2096, getPutativeTSSRegionBasedOnPeakValidSite|2364, investigateTSSRelativeToATGAndTAA|3132, plotTSSRelativeToATGAndTAAHistogram|3494, printGeneBasedTSSLog|3798
#	secondaryDependOnSub: getCoverageOfItemRngType_multiStrand|1897, getMaxCovAndPosAroundGeneEdge|2296, ggplotHistogram|2780, printBothFHAndStdout|3671, printMaxCovAndPosAroundGeneEdgeInfo|3955, reportStatus|4274
#
#<section ID="investigateGeneTSSRelativePos" num="8">
#----Calculate the noiseIndex 
my ($putativeTSSRegHsh_ref, $validTSSDstrbtnSenseHsh_ref) = &getPutativeTSSRegionBasedOnPeakValidSite($mRNAInfoHsh_ref, $mRNAByCntgHsh_ref, $filterTEXCovPlsPathHsh_ref, $TSSSearchRngHsh_ref, $ggplotDirHsh_ref, $resultLogDir, $maxThread, $fixPutativeTSSRegHsh_ref);#->2364
my ($noiseIndexHsh_ref) = &calculateNoiseIndex($putativeTSSRegHsh_ref, $validTSSDstrbtnSenseHsh_ref, $TSSSearchRngHsh_ref, $resultStorableDir, $resultLogDir);#->625

#----plot the TSS distribution around ATG and TAA
my ($geneBasedTSSInfoHsh_ref, $TSSmRNAEndFreqHsh_ref) = &investigateTSSRelativeToATGAndTAA($filterTEXCovPlsPathHsh_ref, $mRNAInfoHsh_ref, $TSSSearchRngHsh_ref, $geneBasedTSSInfoHshPath);#->3132
my ($strndEndValidTSSInfoHsh_ref) = &getGeneWithValidTSSAtGeneEnd($geneBasedTSSInfoHsh_ref, $TSSmRNAEndFreqHsh_ref, $strndEndValidTSSInfoHshPath, $validTSSLimitHsh_ref);#->2096
&plotTSSRelativeToATGAndTAAHistogram($TSSmRNAEndFreqHsh_ref, $ggplotDirHsh_ref, $strndEndValidTSSInfoHsh_ref, $validTSSLimitHsh_ref);#->3494
&printGeneBasedTSSLog($geneBasedTSSLogPath, $geneBasedTSSInfoHsh_ref, $strndEndValidTSSInfoHsh_ref);#->3798
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 9_analyzeNucleotideComposition
#	primaryDependOnSub: defineDegenerateNucleotideHsh|1089, getBaseAtTSSAndExon|1699, getSequenceAroundmRNAReferencePoint|2642, getSiteBaseMismatchAndProportion|2702, plotBaseCompositionAroundmRNAReferencePoint|3302, plotWeblogoAroundTSS|3546, printmRNARefPtLog|4035
#	secondaryDependOnSub: calculateBaseCompositionInAlignments|581, createWeblogo|1045, currentTime|1071, ggplotXYLinesMultipleSamples|2961, reportStatus|4274
#
#<section ID="analyzeNucleotideComposition" num="9">
my ($degenNtHsh_ref) = &defineDegenerateNucleotideHsh();#->1089
my ($siteBaseStringWithReadHsh_ref, $siteBaseGenomeProportionHsh_ref, $mRNARefPtHsh_ref) = &getBaseAtTSSAndExon($strndEndValidTSSInfoHsh_ref, $mRNAInfoHsh_ref, $TEXBaseCompositionPileupPathHsh_ref, $rd5EndPlsInfoHsh_ref, $libTypeInfoHsh_ref, $fastaHsh_ref, $exonEndTrim, $validTSSLimitHsh_ref, $mRNARefPtHshPath);#->1699
&printmRNARefPtLog($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $resultLogDir);#->4035
my ($siteBaseWithReadProportionHsh_ref, $siteBaseWithReadMismatchHsh_ref) = &getSiteBaseMismatchAndProportion($siteBaseStringWithReadHsh_ref, $siteBaseGenomeProportionHsh_ref);#->2702
my ($seqAroundSiteHsh_ref) = &getSequenceAroundmRNAReferencePoint($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $fastaHsh_ref, $resultFastaDir);#->2642
&plotWeblogoAroundTSS($seqAroundSiteHsh_ref, $weblogoDirHsh_ref); #->3546
&plotBaseCompositionAroundmRNAReferencePoint($seqAroundSiteHsh_ref, $ggplotDirHsh_ref, $trimBaseComPlotRngHsh_ref);#->3302
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 10_analyzeTSSMotif
#	primaryDependOnSub: plotMotifOccurenceWithMASTMatch|3361, searchMotifAroundmRNAReferencePoint|4318
#	secondaryDependOnSub: checkRunningThreadAndWaitToJoin|981, currentTime|1071, generateDREMECmd|1527, generateMEMECmd|1581, getDREMEMotif|2034, getMASTLogPostionalData|2179, getMEMEMotif|2229, ggplotXYLinesMultipleSamples|2961, reportStatus|4274, runMAST|4297
#
#<section ID="analyzeTSSMotif" num="10">
my ($dremeXMLPathHsh_ref, $memeXMLPathHsh_ref) = &searchMotifAroundmRNAReferencePoint($seqAroundSiteHsh_ref, $resultDremeDir, $resultMemeDir, $forceRunMemeAndDreme, $hardCodedIGVPathHsh_ref, $memeDremeToRunHsh_ref);#->4318
&plotMotifOccurenceWithMASTMatch($dremeXMLPathHsh_ref, $memeXMLPathHsh_ref, $ggplotDirHsh_ref, $seqAroundSiteHsh_ref, $motifEValueLimitHsh_ref, $resultMastDir, $degenNtHsh_ref, $hardCodedIGVPathHsh_ref, $maxHitPVal);#->3361
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 11_outputWiggleXML
#	primaryDependOnSub: gatherIGVXMLTrack|1460, outputOriginalAndScaledRd5EndWiggle|3234, printIGVXML|3848
#	secondaryDependOnSub: checkRunningThreadAndWaitToJoin|981, currentTime|1071, printWiggleSingleTrackFromCntgCovPlsPathHsh|3997, reportStatus|4274
#
#<section ID="outputWiggleXML" num="11">
#-----print the original and called wiggles
&outputOriginalAndScaledRd5EndWiggle($libTypeInfoHsh_ref, $rd5EndPlsInfoHsh_ref, $forcePoolLibRd5End);#->3234

#-----output the XML
my ($IGVTrackInfoHsh_ref) = &gatherIGVXMLTrack($libTypeInfoHsh_ref, $trnsfrgGFFPathHsh_ref, $wigglePathHsh_ref);#->1460
&printIGVXML($fontSize, $XMLPath, $hardCodedIGVPathHsh_ref, $IGVGenomeID, $IGVTrackInfoHsh_ref, $copyIGVTrack, $resultIGVTrackCopyDir);#->3848
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 12_finishingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|3694
#	secondaryDependOnSub: currentTime|1071
#
#<section ID="finishingTasks" num="12">
#-----open the XML in finder
#system ("open $outDir/result/$resultDirTag/ggplot/pdf/");

#-----print the log
&printCMDLogOrFinishMessage("finishMessage");#->3694
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

#====================================================================================================================================================#
#List of subroutines by category
#
#	XML [n=1]:
#		printIGVXML
#
#	alignment [n=1]:
#		calculateBaseCompositionInAlignments
#
#	baseComposition [n=2]:
#		calculateBaseCompositionInAlignments, defineDegenerateNucleotideHsh
#
#	coverage [n=5]:
#		getCoverageOfItemRngType_multiStrand, getMaxCovAndPosAroundGeneEdge, getPutativeTSSRegionBasedOnPeakValidSite
#		getRd5EndSurroundingmRNAATG, printMaxCovAndPosAroundGeneEdgeInfo
#
#	fasta [n=1]:
#		readMultiFasta
#
#	general [n=8]:
#		checkGeneInfo, currentTime, getCtgryGeneInfo
#		printBothFHAndStdout, printCMDLogOrFinishMessage, readMultiFasta
#		readParameters, reportStatus
#
#	gff [n=6]:
#		checkGeneInfo, findDistanceBetweenATGAndTrnsfrgEnd, getPutativeTSSRegionBasedOnPeakValidSite
#		getRd5EndSurroundingmRNAATG, printGFF_oneRNAPerGene_chooseStrnd_filterAry, readGFF_oneRNAPerGene
#
#	ggplot [n=5]:
#		ggplotHistogram, ggplotTwoSampleHistogram, ggplotXYErrorBarPlot
#		ggplotXYLinesMultipleSamples, ggplotXYScatterPlot
#
#	math [n=1]:
#		calculateXYPairRatio
#
#	multithread [n=3]:
#		checkRunningThreadAndWaitToJoin, generateThreadHshWithRandomCntg, generateThreadHshWithRandomItem
#
#	plotInR [n=5]:
#		ggplotHistogram, ggplotTwoSampleHistogram, ggplotXYErrorBarPlot
#		ggplotXYLinesMultipleSamples, ggplotXYScatterPlot
#
#	range [n=1]:
#		checkOverlapAndProximity_withMargin
#
#	reporting [n=3]:
#		currentTime, printBothFHAndStdout, printMaxCovAndPosAroundGeneEdgeInfo
#
#	specific [n=13]:
#		analyzeSTDTEXRd5EndRatio, calculateNoiseIndex, defineGoldenmRNASet
#		defineLibTypeInfo, downSampleGoldenATGRd5EndDataHshForTesting, filterTEXSiteBasedOnTrainingCutoff
#		findDistanceBetweenATGAndTrnsfrgEnd, gatherIGVXMLTrack, getPutativeTSSRegionBasedOnPeakValidSite
#		getRd5EndSurroundingmRNAATG, plotWeblogoAroundTSS, printIGVXML
#		printmRNARefPtLog
#
#	storable [n=3]:
#		createEmptyGenomeCovPerlStorable, getIndivCntgCovPlsPath, zipUnzipCntgCovInPlsPathHsh
#
#	thirdPartyApp [n=1]:
#		plotWeblogoAroundTSS
#
#	thridPartyApp [n=6]:
#		createWeblogo, generateDREMECmd, generateMEMECmd
#		getMASTLogPostionalData, getMEMEMotif, runMAST
#
#	unassigned [n=22]:
#		getAndPrintTrnsfrg, getBaseAtTSSAndExon, getDREMEMotif
#		getGeneWithValidTSSAtGeneEnd, getRd5EndCountForExonAndTSSOfGoldenmRNASet, getSequenceAroundmRNAReferencePoint
#		getSiteBaseMismatchAndProportion, identifiyPeakSiteWithinRegion, intervalizeXYPairsForErrorBarPlot
#		investigateTSSRelativeToATGAndTAA, outputOriginalAndScaledRd5EndWiggle, outputTSSResultWiggle
#		plotBaseCompositionAroundmRNAReferencePoint, plotMotifOccurenceWithMASTMatch, plotTSSRelativeToATGAndTAAHistogram
#		poolLibRd5EndCov, predictBonaFideTSSinGoldenmRNASet, printGeneBasedTSSLog
#		searchMotifAroundmRNAReferencePoint, storeAndScaleRd5EndFromBothLibraries, trainingUsingExonVsTssData
#		warningToProceed
#
#	wiggle [n=1]:
#		printWiggleSingleTrackFromCntgCovPlsPathHsh
#
#====================================================================================================================================================#

sub analyzeSTDTEXRd5EndRatio {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: calculateXYPairRatio|685, ggplotTwoSampleHistogram|2848, reportStatus|4274
#	appearInSub: trainingUsingExonVsTssData|4548
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_trainningSTDTEXRatio|308
#	input: $countCutoffPct, $ggplotDirHsh_ref, $minOriginalCount, $ratioCutoffPct, $rd5EndTSSExonCountHsh_ref
#	output: $TSSvsExonCutoffHsh_ref
#	toCall: my ($TSSvsExonCutoffHsh_ref) = &analyzeSTDTEXRd5EndRatio($rd5EndTSSExonCountHsh_ref, $ggplotDirHsh_ref, $minOriginalCount, $ratioCutoffPct, $countCutoffPct);
#	calledInLine: 4574
#....................................................................................................................................................#

	my ($rd5EndTSSExonCountHsh_ref, $ggplotDirHsh_ref, $minOriginalCount, $ratioCutoffPct, $countCutoffPct) = @_;
	
	my $countBinWidth = 0.2;
	my $ratioBinWidth = 0.2;
	my $dataPtMax = 2000;
	my $ratioExtraArg = ' + xlim(c(-20, 20))';
	my $countExtraArg = ' + xlim(c(0, 20))';

	my $TSSvsExonCutoffHsh_ref = {};

	my $tmpPlotItemHsh_ref = {};
	$tmpPlotItemHsh_ref->{'tss'}{'ratioHist'} = 'StdTexRatioHist_tss';
	$tmpPlotItemHsh_ref->{'exon'}{'ratioHist'} = 'StdTexRatioHist_exon';
	$tmpPlotItemHsh_ref->{'tss'}{'TEXCountZeroSTDHist'} = 'TexCountHist_tss';
	$tmpPlotItemHsh_ref->{'exon'}{'TEXCountZeroSTDHist'} = 'TexCountHist_exon';

	my $twoSamplePlotAryHsh_ref = {};

	my $sub_ggplotDir = "STD.vs.TEX.ratio";
	system "mkdir -pm 777 $ggplotDirHsh_ref->{$_}/$sub_ggplotDir/" foreach (keys %{$ggplotDirHsh_ref});
	
	foreach my $tssOrExon (keys %{$rd5EndTSSExonCountHsh_ref}) {

		&reportStatus("Analyzing STD TEX Rd5End Ratio for $tssOrExon in histogram", 10, "\n");#->4274

		#-----choose the input value in linear or log2 scale here
		my (undef, $YToXRatioAry_ref, undef, undef) = &calculateXYPairRatio(\@{$rd5EndTSSExonCountHsh_ref->{$tssOrExon}{'original'}{'linear'}}, $minOriginalCount);#->685
		$twoSamplePlotAryHsh_ref->{'ratio'}{$tssOrExon}{'plotAry_ref'} = $YToXRatioAry_ref;

		#-----choose the input value in linear or log2 scale here
		my (undef, undef, undef, $YValXZeroAry_ref) = &calculateXYPairRatio(\@{$rd5EndTSSExonCountHsh_ref->{$tssOrExon}{'original'}{'linear'}}, $minOriginalCount);#->685
		$twoSamplePlotAryHsh_ref->{'TEXCountZeroSTD'}{$tssOrExon}{'plotAry_ref'} = $YValXZeroAry_ref;
		
		my $tmpTwoSamplePlotInfoHsh_ref = {};
		$tmpTwoSamplePlotInfoHsh_ref->{'ratio'}{'plotItem'} = 'StdTexRatioHist_both';
		$tmpTwoSamplePlotInfoHsh_ref->{'ratio'}{'binWidth'} = $ratioBinWidth;
		$tmpTwoSamplePlotInfoHsh_ref->{'ratio'}{'extraArg'} = $ratioExtraArg;
		$tmpTwoSamplePlotInfoHsh_ref->{'ratio'}{'cutoffPct'} = $ratioCutoffPct;

		$tmpTwoSamplePlotInfoHsh_ref->{'TEXCountZeroSTD'}{'plotItem'} = 'TexCountHist_both';
		$tmpTwoSamplePlotInfoHsh_ref->{'TEXCountZeroSTD'}{'binWidth'} = $countBinWidth;
		$tmpTwoSamplePlotInfoHsh_ref->{'TEXCountZeroSTD'}{'extraArg'} = $countExtraArg;
		$tmpTwoSamplePlotInfoHsh_ref->{'TEXCountZeroSTD'}{'cutoffPct'} = $countCutoffPct;

		foreach my $ratioOrTEXCountZeroSTD (keys %{$twoSamplePlotAryHsh_ref}) {
			my $plotAryHsh_ref = {};
			my $cutoffPct = $tmpTwoSamplePlotInfoHsh_ref->{$ratioOrTEXCountZeroSTD}{'cutoffPct'};
			foreach my $tssOrExon (%{$twoSamplePlotAryHsh_ref->{$ratioOrTEXCountZeroSTD}}) {
				my $plotAry_ref = $twoSamplePlotAryHsh_ref->{$ratioOrTEXCountZeroSTD}{$tssOrExon}{'plotAry_ref'};
				$plotAryHsh_ref->{$tssOrExon} = $plotAry_ref;
			}

			#----calculate the value limit of cutoffPct for adding a vertical line in the histogram
			my $valueStatObj = Statistics::Descriptive::Full->new();
			$valueStatObj->add_data(@{$plotAryHsh_ref->{'exon'}});
			my $tmpCutoffValueHsh_ref = {};
			$tmpCutoffValueHsh_ref->{'linear'} = sprintf "%.10f", $valueStatObj->percentile($cutoffPct);
			$tmpCutoffValueHsh_ref->{'log2'} = -20;
			$tmpCutoffValueHsh_ref->{'log2'} = sprintf "%.10f", log($tmpCutoffValueHsh_ref->{'linear'})/log(2) if $tmpCutoffValueHsh_ref->{'linear'} > 0;
			$TSSvsExonCutoffHsh_ref->{$ratioOrTEXCountZeroSTD} = $tmpCutoffValueHsh_ref->{'linear'};
			
			my $item = $tmpTwoSamplePlotInfoHsh_ref->{$ratioOrTEXCountZeroSTD}{'plotItem'};
			my $dataPath = "$ggplotDirHsh_ref->{'dat'}/$sub_ggplotDir/$item.dat";
			my $pdfPath = "$ggplotDirHsh_ref->{'pdf'}/$sub_ggplotDir/$item.pdf";
			my $RScriptPath = "$ggplotDirHsh_ref->{'R'}/$sub_ggplotDir/$item.R";
			my $logPath = "$ggplotDirHsh_ref->{'log'}/$sub_ggplotDir/$item.log";
			my $leftXAxisPercentileLimit = 'min';
			my $rightXAxisPercentileLimit = 'max';
			my $xAxis = "$ratioOrTEXCountZeroSTD";
			my $binWidth = $tmpTwoSamplePlotInfoHsh_ref->{$ratioOrTEXCountZeroSTD}{'binWidth'};
			my $log2OrLinear = 'log2';
			my $extraArg = $tmpTwoSamplePlotInfoHsh_ref->{$ratioOrTEXCountZeroSTD}{'extraArg'};
			$extraArg .= " + geom_vline(xintercept=c($tmpCutoffValueHsh_ref->{$log2OrLinear}), linetype=\"dotted\") + annotate(\"text\", x=$tmpCutoffValueHsh_ref->{$log2OrLinear}, y=0, label=\"$cutoffPct\%\=$tmpCutoffValueHsh_ref->{$log2OrLinear}\", vjust=-0.2, hjust=-0.1, angle=90)";
			&ggplotTwoSampleHistogram($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftXAxisPercentileLimit, $rightXAxisPercentileLimit, $xAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear);#->2848
		}
	}
	
	return ($TSSvsExonCutoffHsh_ref);
}
sub calculateBaseCompositionInAlignments {
#....................................................................................................................................................#
#	subroutineCategory: alignment, baseComposition
#	dependOnSub: reportStatus|4274
#	appearInSub: plotBaseCompositionAroundmRNAReferencePoint|3302
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 9_analyzeNucleotideComposition|349
#	input: $seqAlignHsh_ref
#	output: $baseCompByBaseHsh_ref
#	toCall: my ($baseCompByBaseHsh_ref) = &calculateBaseCompositionInAlignments($seqAlignHsh_ref);
#	calledInLine: 3342
#....................................................................................................................................................#

	my ($seqAlignHsh_ref) = @_;
	
	my $baseCountByPosHsh_ref = {};
	my $baseCompByBaseHsh_ref = {};
	my $validSeqNum = 0;
	my $tmpLengthHsh_ref = {};
	
	foreach my $seqName (keys %{$seqAlignHsh_ref}) {
		next if $seqAlignHsh_ref->{$seqName} =~ m/[^ATGCatgc]/;
		$validSeqNum++;
		my @seqAry = split //, $seqAlignHsh_ref->{$seqName};
		$tmpLengthHsh_ref->{@seqAry}++;
		for my $pos (0..$#seqAry) {
			my $base = $seqAry[$pos];
			$base =~ tr/atgc/ATGC/;
			$baseCountByPosHsh_ref->{$pos}{$base}++;
		}
	}
	
	my $lengthNum = keys %{$tmpLengthHsh_ref};
	&reportStatus("WARNING: Length of the sequences in the alignment is not uniform", 10, "\n") if $lengthNum > 1;#->4274
	foreach my $base (qw/A T G C/) {
		foreach my $pos (sort {$a <=> $b} keys %{$baseCountByPosHsh_ref}) {
			$baseCountByPosHsh_ref->{$pos}{$base} = 0 if not $baseCountByPosHsh_ref->{$pos}{$base};
			$baseCompByBaseHsh_ref->{$base}{$pos} = $baseCountByPosHsh_ref->{$pos}{$base}/$validSeqNum;
		}
	}
	
	return $baseCompByBaseHsh_ref;
	
}
sub calculateNoiseIndex {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: printBothFHAndStdout|3671, reportStatus|4274
#	appearInSub: >none
#	primaryAppearInSection: 8_investigateGeneTSSRelativePos|331
#	secondaryAppearInSection: >none
#	input: $TSSSearchRngHsh_ref, $putativeTSSRegHsh_ref, $resultLogDir, $resultStorableDir, $validTSSDstrbtnSenseHsh_ref
#	output: $noiseIndexHsh_ref
#	toCall: my ($noiseIndexHsh_ref) = &calculateNoiseIndex($putativeTSSRegHsh_ref, $validTSSDstrbtnSenseHsh_ref, $TSSSearchRngHsh_ref, $resultStorableDir, $resultLogDir);
#	calledInLine: 338
#....................................................................................................................................................#
	my ($putativeTSSRegHsh_ref, $validTSSDstrbtnSenseHsh_ref, $TSSSearchRngHsh_ref, $resultStorableDir, $resultLogDir) = @_;
	
	my $exonEndTrim = 10;
	my $margin = $TSSSearchRngHsh_ref->{'ATG'}{'up'};
	my %tmpRngBoundHsh = ();
	$tmpRngBoundHsh{'exon'}{'start'} = $margin + $exonEndTrim;
	$tmpRngBoundHsh{'TSS'}{'start'} = $margin + $putativeTSSRegHsh_ref->{'lowerLimit'}; #--e.g. 150 + (-15)
	$tmpRngBoundHsh{'TSS'}{'end'} = $margin + $putativeTSSRegHsh_ref->{'upperLimit'};;
	
	my %siteCountHsh = ();
	
	&reportStatus("Calculating noise index", 10, "\n");#->4274
	
	foreach my $mRNAID (keys %{$validTSSDstrbtnSenseHsh_ref}) {

		my $covAry_ref = $validTSSDstrbtnSenseHsh_ref->{$mRNAID}{'s'};
		$tmpRngBoundHsh{'exon'}{'end'} = $#{$covAry_ref} - $margin - $exonEndTrim;
		
		foreach my $exonOrTSS (keys %tmpRngBoundHsh) {
			my $regAry_ref = [@$covAry_ref[$tmpRngBoundHsh{$exonOrTSS}{'start'}..$tmpRngBoundHsh{$exonOrTSS}{'end'}]];
			foreach (@{$regAry_ref}) {
				$siteCountHsh{$exonOrTSS}{'valid'}++ if $_ > 0;
				$siteCountHsh{$exonOrTSS}{'total'}++;
			}
		}
	}
	
	my $exon_ValidSitePerNt = sprintf "%.10f", $siteCountHsh{'exon'}{'valid'}/$siteCountHsh{'exon'}{'total'};
	my $TSS_ValidSitePerNt = sprintf "%.10f", $siteCountHsh{'TSS'}{'valid'}/$siteCountHsh{'TSS'}{'total'};
	my $noiseIndex = sprintf "%.10f", $exon_ValidSitePerNt/$TSS_ValidSitePerNt;

	my $noiseIndexFH;
	open $noiseIndexFH, ">", "$resultLogDir/noise.index.txt";
	&printBothFHAndStdout("exon_ValidSitePerNt = $exon_ValidSitePerNt", 10, $noiseIndexFH);#->3671
	&printBothFHAndStdout("TSS_ValidSitePerNt = $TSS_ValidSitePerNt", 10, $noiseIndexFH);#->3671
	&printBothFHAndStdout("noiseIndex = $noiseIndex", 10, $noiseIndexFH);#->3671
	close $noiseIndexFH;
	
	my $noiseIndexHsh_ref = {
		"noiseIndex" => $noiseIndex,
		"exon_ValidSitePerNt" => $exon_ValidSitePerNt,
		"TSS_ValidSitePerNt" => $TSS_ValidSitePerNt,
	};
	
	store($noiseIndexHsh_ref, "$resultStorableDir/noiseIndexHsh.pls");

	return ($noiseIndexHsh_ref);
}
sub calculateXYPairRatio {
#....................................................................................................................................................#
#	subroutineCategory: math
#	dependOnSub: >none
#	appearInSub: analyzeSTDTEXRd5EndRatio|489
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $XYPairAry_ref, $minValAsNonZero
#	output: $XToYRatioAry_ref, $XValYZeroAry_ref, $YToXRatioAry_ref, $YValXZeroAry_ref
#	toCall: my ($XToYRatioAry_ref, $YToXRatioAry_ref, $XValYZeroAry_ref, $YValXZeroAry_ref) = &calculateXYPairRatio($XYPairAry_ref, $minValAsNonZero);
#	calledInLine: 527, 531
#....................................................................................................................................................#
	
	my ($XYPairAry_ref, $minValAsNonZero) = @_;
	
	my $XToYRatioAry_ref = ();
	my $YToXRatioAry_ref = ();
	my $YValXZeroAry_ref = ();
	my $XValYZeroAry_ref = ();
	
	foreach my $XYPair (@{$XYPairAry_ref}) {
		my ($XVal, $YVal) = split /,/, $XYPair;
		if ($YVal >= $minValAsNonZero and $XVal >= $minValAsNonZero) {
			push @{$XValYZeroAry_ref}, $XVal/$YVal;
			push @{$YToXRatioAry_ref}, $YVal/$XVal;
		} elsif ($YVal < $minValAsNonZero and $XVal >= $minValAsNonZero) {
			push @{$XValYZeroAry_ref}, $XVal;
		} elsif ($YVal >= $minValAsNonZero and $XVal < $minValAsNonZero) {
			push @{$YValXZeroAry_ref}, $YVal;
		}
	}
	return ($XToYRatioAry_ref, $YToXRatioAry_ref, $XValYZeroAry_ref, $YValXZeroAry_ref);
}
sub checkGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: currentTime|1071
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|273
#	secondaryAppearInSection: >none
#	input: $geneInfoHsh_ref
#	output: none
#	toCall: &checkGeneInfo($geneInfoHsh_ref);
#	calledInLine: 286
#....................................................................................................................................................#
	
	my ($geneInfoHsh_ref) = @_;
	
	print "[".&currentTime()."] Checking gene categories.\n";#->1071
	my $ctrgyCountHsh_ref = {};
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		$ctrgyCountHsh_ref->{$ctgry}++;
	}
	
	foreach my $ctgry (sort keys %{$ctrgyCountHsh_ref}) {
		print "[".&currentTime()."] Item in $ctgry = $ctrgyCountHsh_ref->{$ctgry}\n";#->1071
	}
}
sub checkOverlapAndProximity_withMargin {
#....................................................................................................................................................#
#	subroutineCategory: range
#	dependOnSub: currentTime|1071, generateThreadHshWithRandomCntg|1618, reportStatus|4274
#	appearInSub: findDistanceBetweenATGAndTrnsfrgEnd|1386
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $checkPrxmty, $maxThread, $qryInfoHsh_ref, $qryMargin, $qryRngType, $refInfoHsh_ref, $refMargin, $refRngType, $reportExactMatch
#	output: $hitAndPrxmtyByQryHsh_ref, $hitAndPrxmtyByRefHsh_ref
#	toCall: my ($hitAndPrxmtyByRefHsh_ref, $hitAndPrxmtyByQryHsh_ref) = &checkOverlapAndProximity_withMargin($refInfoHsh_ref, $qryInfoHsh_ref, $checkPrxmty, $reportExactMatch, $maxThread, $refRngType, $qryRngType, $refMargin, $qryMargin);
#	calledInLine: 1407
#....................................................................................................................................................#
	#---incoming variables
	my ($refInfoHsh_ref, $qryInfoHsh_ref, $checkPrxmty, $reportExactMatch, $maxThread, $refRngType, $qryRngType, $refMargin, $qryMargin) = @_;

	#---outgoing variables
	my $refGeneNumTotal = 0;
	
	#---make a tmpHsh to contain all cntgs the have either ref and qry
	my $tmpCntgHsh_ref = {};
	my $refCntgHsh_ref = {};
	my $qryCntgHsh_ref = {};
	
	foreach my $refGeneID (keys %{$refInfoHsh_ref}) {
		if ($refInfoHsh_ref->{$refGeneID}{$refRngType}) {
			@{$refInfoHsh_ref->{$refGeneID}{$refRngType}} = sort {$a <=> $b} @{$refInfoHsh_ref->{$refGeneID}{$refRngType}};
			$refGeneNumTotal++;
			$tmpCntgHsh_ref->{$refInfoHsh_ref->{$refGeneID}{'cntg'}}++;
			$refCntgHsh_ref->{$refInfoHsh_ref->{$refGeneID}{'cntg'}}{$refGeneID}++;
		}
	}
	
	foreach my $qryGeneID (keys %{$qryInfoHsh_ref}) {
		if ($qryInfoHsh_ref->{$qryGeneID}{$qryRngType}) {
			@{$qryInfoHsh_ref->{$qryGeneID}{$qryRngType}} = sort {$a <=> $b} @{$qryInfoHsh_ref->{$qryGeneID}{$qryRngType}};
			$tmpCntgHsh_ref->{$qryInfoHsh_ref->{$qryGeneID}{'cntg'}}++;
			$qryCntgHsh_ref->{$qryInfoHsh_ref->{$qryGeneID}{'cntg'}}{$qryGeneID}++;
		}
	}
	
	my $totalCntgNum = keys %{$tmpCntgHsh_ref};
	my @cntgAry = keys %{$tmpCntgHsh_ref};

	my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomCntg($maxThread, \@cntgAry);#->1618
	my $refGeneNumProc :shared = 0;
	my %threadHsh = ();

	foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
		my $cntgAry_ref = \@{$randCntgInThreadHsh_ref->{$threadNum}};
		my $cntgNum = @{$randCntgInThreadHsh_ref->{$threadNum}};
		&reportStatus("$cntgNum cntgs spawned to thread $threadNum", 20, "\n");#->4274

		#---spawn a new thread
		($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781
		
			sub {
				my ($cntgAry_ref) = @_;
	
				my $hitAndPrxmtyByRefHsh_InThr_ref = {};
				my $hitAndPrxmtyByQryHsh_InThr_ref = {};

				foreach my $cntg (@{$cntgAry_ref}) {

					#---update the on-screen progress
					#&reportStatus("Finding overlaping on $cntg", 20, "\r");#->4274
		
					my $tmpPrxmtyByRefHsh_ref = {};
					my $tmpPrxmtyByQryHsh_ref = {};

					if ((exists $qryCntgHsh_ref->{$cntg}) and (exists $refCntgHsh_ref->{$cntg})) {#---if there are both ref and qry can both be found on cntg
						foreach my $refGeneID (keys %{$refCntgHsh_ref->{$cntg}}) {#--- all ftur on the $strnd of $cntg of refGff
							my ($refStart, $refEnd) = ($refInfoHsh_ref->{$refGeneID}{$refRngType}->[0]-$refMargin, $refInfoHsh_ref->{$refGeneID}{$refRngType}->[-1]+$refMargin);
				
							$refGeneNumProc++;
				
							&reportStatus("$refGeneNumProc of $refGeneNumTotal reference genes checked", 50, "\r");#->4274

							foreach my $qryGeneID (keys %{$qryCntgHsh_ref->{$cntg}}) {#--- all ftur on the $strnd of $cntg of QryGtf
	
								my $samestrnd = "no";
								$samestrnd = "yes" if ($refInfoHsh_ref->{$refGeneID}{'strnd'} eq $qryInfoHsh_ref->{$qryGeneID}{'strnd'});
								my ($qryStart, $qryEnd) = ($qryInfoHsh_ref->{$qryGeneID}{$qryRngType}->[0]-$qryMargin, $qryInfoHsh_ref->{$qryGeneID}{$qryRngType}->[-1]+$qryMargin);

								my $scene;
								my $ovrlpSize;

								if (($refStart == $qryStart) && ($refEnd == $qryEnd)) {#---scene 0
									$scene = 'exactMatch';
									$ovrlpSize = $qryEnd - $qryStart;
						
								} elsif (($refStart<=$qryStart)&&($refEnd>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 1
									$scene = 'overlapTail';
									$ovrlpSize = $refEnd - $qryStart;

								} elsif (($refStart>=$qryStart)&&($refStart<=$qryEnd)&&($refEnd>=$qryEnd)) {#---scene 2
									$scene = 'overlapHead';
									$ovrlpSize = $qryEnd - $refStart;

								} elsif (($refStart<=$qryStart)&&($refEnd>=$qryEnd)) {#---scene 3
									$scene = 'cover';
									$ovrlpSize = $qryEnd - $qryStart;

								} elsif (($refStart>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 4
									$scene = 'within';
									$ovrlpSize = $refEnd - $refStart;

								#------Proximity with ref's tail proximal to qry's head
								} elsif (($refEnd<=$qryStart)&&($refEnd<$qryEnd)) {#---scene 5 ---> ref Tail, qry Head

									$scene = 'prxmtyTail';

									if ($checkPrxmty eq "yes") {
										my $tmpPrmxty = $qryStart - $refEnd;
										$tmpPrxmtyByRefHsh_ref->{'XS'}{$refGeneID}{"T"}{$qryGeneID} = $tmpPrmxty;
										$tmpPrxmtyByQryHsh_ref->{'XS'}{$qryGeneID}{"H"}{$refGeneID} = $tmpPrmxty;

										if ($samestrnd eq "yes") {
											$tmpPrxmtyByRefHsh_ref->{'SS'}{$refGeneID}{"T"}{$qryGeneID} = $tmpPrmxty;
											$tmpPrxmtyByQryHsh_ref->{'SS'}{$qryGeneID}{"H"}{$refGeneID} = $tmpPrmxty;
										}
									}

								#------Proximity with ref's head proximal to qry's tail
								} elsif (($refStart>=$qryEnd)&&($refStart>$qryStart)) {#---scene 6 ---> ref Head, qry Tail

									$scene = 'prxmtyHead';

									if ($checkPrxmty eq "yes") {
										my $tmpPrmxty = $refStart - $qryEnd;
										$tmpPrxmtyByRefHsh_ref->{'XS'}{$refGeneID}{"H"}{$qryGeneID} = $tmpPrmxty;
										$tmpPrxmtyByQryHsh_ref->{'XS'}{$qryGeneID}{"T"}{$refGeneID} = $tmpPrmxty;

										if ($samestrnd eq "yes") {
											$tmpPrxmtyByRefHsh_ref->{'SS'}{$refGeneID}{"H"}{$qryGeneID} = $tmpPrmxty;
											$tmpPrxmtyByQryHsh_ref->{'SS'}{$qryGeneID}{"T"}{$refGeneID} = $tmpPrmxty;
										}
									}

								} else {#---BUG! possibly other scene?
									#print "[".&currentTime()."] refStart=$refStart; refEnd=$refEnd; qryStart=$qryStart; qryEnd=$qryEnd\n";#->1071
									die "Unexpected overlapping scene between $refGeneID and $qryGeneID. It's a Bug. Program qutting.\n";
								}
					
								if ($scene ne 'prxmtyTail' and $scene ne 'prxmtyHead' and not ($reportExactMatch eq 'no' and $scene eq 'exactMatch')) {

									@{$hitAndPrxmtyByRefHsh_InThr_ref->{'XS'}{'hit'}{$refGeneID}{$qryGeneID}} = ($scene, $ovrlpSize);
									@{$hitAndPrxmtyByQryHsh_InThr_ref->{'XS'}{'hit'}{$qryGeneID}{$refGeneID}} = ($scene, $ovrlpSize);

									if ($samestrnd eq "yes") {
										@{$hitAndPrxmtyByRefHsh_InThr_ref->{'SS'}{'hit'}{$refGeneID}{$qryGeneID}} = ($scene, $ovrlpSize);
										@{$hitAndPrxmtyByQryHsh_InThr_ref->{'SS'}{'hit'}{$qryGeneID}{$refGeneID}} = ($scene, $ovrlpSize);
									}
								}
							}
						}
					}

					#---find the closest proximity for all refs
					if ($checkPrxmty eq "yes") {
						my $refQryRefHsh_ref = {};

						$refQryRefHsh_ref->{'ref'}{'tmpPrxmtyHsh_ref'} = $tmpPrxmtyByRefHsh_ref;
						$refQryRefHsh_ref->{'ref'}{'cntgHsh_ref'} = $refCntgHsh_ref;
						$refQryRefHsh_ref->{'ref'}{'hitAndPrxmtyHsh_ref'} = $hitAndPrxmtyByRefHsh_InThr_ref;

						$refQryRefHsh_ref->{'qry'}{'tmpPrxmtyHsh_ref'} = $tmpPrxmtyByQryHsh_ref;
						$refQryRefHsh_ref->{'qry'}{'cntgHsh_ref'} = $qryCntgHsh_ref;
						$refQryRefHsh_ref->{'qry'}{'hitAndPrxmtyHsh_ref'} = $hitAndPrxmtyByQryHsh_InThr_ref;
			
						foreach my $refOrQry ('ref', 'qry') {

							my $cntgHsh_ref = $refQryRefHsh_ref->{$refOrQry}{'cntgHsh_ref'};
							my $tmpPrxmtyHsh_ref = $refQryRefHsh_ref->{$refOrQry}{'tmpPrxmtyHsh_ref'};
							my $hitAndPrxmtyHsh_ref = $refQryRefHsh_ref->{$refOrQry}{'hitAndPrxmtyHsh_ref'};

							foreach my $ftur (keys %{$cntgHsh_ref->{$cntg}}) {
								foreach my $XSOrSS ('XS', 'SS') {
									foreach my $HOrT ('H', 'T') {
										$tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{"edge"} = -999 if (not exists $tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT});
										foreach my $otherFtur (sort {$tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{$a} <=> $tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{$b}} keys %{$tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}}) {
											@{$hitAndPrxmtyHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT}} = ($tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{$otherFtur}, $otherFtur);
											last; #---sample the smallest only
										}
									}
								}
							}
						}
					}
				}
				
				return ($hitAndPrxmtyByRefHsh_InThr_ref, $hitAndPrxmtyByQryHsh_InThr_ref);
			}
			,($cntgAry_ref)
		);
	}
	
	
	my %tmpTransferThrDataHsh = ();
	$tmpTransferThrDataHsh{'ref'}{'all'} = {};
	$tmpTransferThrDataHsh{'qry'}{'all'} = {};
	
	while (keys %threadHsh) {
		foreach my $threadNum (keys %threadHsh) {
			my $thr = $threadHsh{$threadNum};
			if (not $thr->is_running()) {
				($tmpTransferThrDataHsh{'ref'}{'thr'}, $tmpTransferThrDataHsh{'qry'}{'thr'}) = $thr->join;
				foreach my $refOrQry (keys %tmpTransferThrDataHsh) {
					my ($allHsh_ref, $thrHsh_ref) = ($tmpTransferThrDataHsh{$refOrQry}{'all'}, $tmpTransferThrDataHsh{$refOrQry}{'thr'});
					foreach my $XSOrSS ('XS', 'SS') {
						foreach my $ftur (keys %{$thrHsh_ref->{$XSOrSS}{'hit'}}) {
							foreach my $hitftur (keys %{$thrHsh_ref->{$XSOrSS}{'hit'}{$ftur}}) {
								@{$allHsh_ref->{$XSOrSS}{'hit'}{$ftur}{$hitftur}} = @{$thrHsh_ref->{$XSOrSS}{'hit'}{$ftur}{$hitftur}};
							}
						}
						
						if ($thrHsh_ref->{$XSOrSS}{'prxmty'}) {
							foreach my $ftur (keys %{$thrHsh_ref->{$XSOrSS}{'prxmty'}}) {
								foreach my $HOrT (keys %{$thrHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}}) {
									@{$allHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT}} = @{$thrHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT}} if $thrHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT};
								}
							}
						}
					}
				}
				delete $threadHsh{$threadNum};
			}
		}
		sleep 1;
	}

	print "\n";

	my $hitAndPrxmtyByRefHsh_ref = $tmpTransferThrDataHsh{'ref'}{'all'};
	my $hitAndPrxmtyByQryHsh_ref = $tmpTransferThrDataHsh{'qry'}{'all'};

	return ($hitAndPrxmtyByRefHsh_ref, $hitAndPrxmtyByQryHsh_ref);
}
sub checkRunningThreadAndWaitToJoin {
#....................................................................................................................................................#
#	subroutineCategory: multithread
#	dependOnSub: reportStatus|4274
#	appearInSub: filterTEXSiteBasedOnTrainingCutoff|1275, outputOriginalAndScaledRd5EndWiggle|3234, searchMotifAroundmRNAReferencePoint|4318, storeAndScaleRd5EndFromBothLibraries|4452
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_analyzeTSSMotif|365, 11_outputWiggleXML|376, 7_filterTEXTrackAfterTrainning|319
#	input: $sleepTime, $verbose
#	output: none
#	toCall: &checkRunningThreadAndWaitToJoin($verbose, $sleepTime);
#	calledInLine: 1374, 3264, 4447, 4544
#....................................................................................................................................................#
	
	my ($verbose, $sleepTime) = @_;
	
	my @runningThrAry = threads->list(threads::running);
	my @joinableThrAry = threads->list(threads::joinable);
	while (@runningThrAry or @joinableThrAry) {
		@runningThrAry = threads->list(threads::running);
		@joinableThrAry = threads->list(threads::joinable);
		foreach my $joinableThr (@joinableThrAry) {
			$joinableThr->detach() if not $joinableThr->is_running();
		}
		my $numThreadRunning = scalar @runningThrAry;
		&reportStatus("The last $numThreadRunning threads are still running", 20, "\r") if $verbose eq 'yes';#->4274

		sleep $sleepTime;
	}
}
sub createEmptyGenomeCovPerlStorable {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: currentTime|1071
#	appearInSub: filterTEXSiteBasedOnTrainingCutoff|1275, poolLibRd5EndCov|3584
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_processInputData|273, 7_filterTEXTrackAfterTrainning|319
#	input: $cntgCovPlsDir, $fastaHsh_ref
#	output: $cntgCovIdxHshPath
#	toCall: my ($cntgCovIdxHshPath) = &createEmptyGenomeCovPerlStorable($cntgCovPlsDir, $fastaHsh_ref);
#	calledInLine: 1300, 3618
#....................................................................................................................................................#

	my ($cntgCovPlsDir, $fastaHsh_ref) = @_;
	
	my $cntgCovPlsIdxHsh_ref = {};
	foreach my $cntg (keys %{$fastaHsh_ref}) {
		print "[".&currentTime()."] Creating empty storabe for $cntg         \r";#->1071
		my $cntgLen = length($fastaHsh_ref->{$cntg});
		my $cntgCovAry_ref = ();
		foreach (1..$cntgLen) {
			push @{$cntgCovAry_ref}, undef;
		}
		my $cntgCovPlsName = "$cntg.ary.pls";
		my $cntgCovPlsPath = "$cntgCovPlsDir/$cntgCovPlsName";
		$cntgCovPlsIdxHsh_ref->{$cntg} = $cntgCovPlsName;
		store($cntgCovAry_ref, "$cntgCovPlsPath");
	}

	my $cntgCovIdxHshPath = "$cntgCovPlsDir/index.hsh.pls";
	store($cntgCovPlsIdxHsh_ref, "$cntgCovIdxHshPath");
	print "\n";
	
	return $cntgCovIdxHshPath;
}
sub createWeblogo {
#....................................................................................................................................................#
#	subroutineCategory: thridPartyApp
#	dependOnSub: >none
#	appearInSub: plotWeblogoAroundTSS|3546
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 9_analyzeNucleotideComposition|349
#	input: $fastaPath, $pdfPath, $seqAlignHsh_ref, $seqType, $title
#	output: none
#	toCall: &createWeblogo($seqAlignHsh_ref, $pdfPath, $fastaPath, $seqType, $title);
#	calledInLine: 3578
#....................................................................................................................................................#

	my ($seqAlignHsh_ref, $pdfPath, $fastaPath, $seqType, $title) = @_;
	
	#---print the fasta
	open (FASTA, ">", $fastaPath);
	foreach my $seqName (keys %{$seqAlignHsh_ref}) {
		print FASTA ">$seqName\n";
		print FASTA "$seqAlignHsh_ref->{$seqName}\n";
	}
	close FASTA;
	
	system ("weblogo --fin $fastaPath --datatype fasta --format pdf --fout $pdfPath --sequence-type $seqType --title \"$title\" --color-scheme classic");
	
}
sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general, reporting
#	dependOnSub: >none
#	appearInSub: checkGeneInfo|718, checkOverlapAndProximity_withMargin|744, createEmptyGenomeCovPerlStorable|1010, defineGoldenmRNASet|1123, defineLibTypeInfo|1166, gatherIGVXMLTrack|1460, getAndPrintTrnsfrg|1672, getBaseAtTSSAndExon|1699, getCtgryGeneInfo|1997, getIndivCntgCovPlsPath|2146, getSequenceAroundmRNAReferencePoint|2642, getSiteBaseMismatchAndProportion|2702, plotMotifOccurenceWithMASTMatch|3361, printBothFHAndStdout|3671, printCMDLogOrFinishMessage|3694, readGFF_oneRNAPerGene|4075, readMultiFasta|4174, reportStatus|4274, warningToProceed|4583, zipUnzipCntgCovInPlsPathHsh|4603
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|114, 10_analyzeTSSMotif|365, 11_outputWiggleXML|376, 12_finishingTasks|391, 5_processInputData|273, 9_analyzeNucleotideComposition|349
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 732, 740, 882, 1026, 1139, 1161, 1201, 1220, 1479, 1686, 1837, 1891, 2016, 2029, 2174, 2660, 2724, 3402, 3688, 3689, 3714, 3717, 3722, 4095, 4192, 4291, 4596, 4598, 4618
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub defineDegenerateNucleotideHsh {
#....................................................................................................................................................#
#	subroutineCategory: baseComposition
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 9_analyzeNucleotideComposition|349
#	secondaryAppearInSection: >none
#	input: none
#	output: $degenNtHsh_ref
#	toCall: my ($degenNtHsh_ref) = &defineDegenerateNucleotideHsh();
#	calledInLine: 354
#....................................................................................................................................................#

	my $degenNtHsh_ref = {};
	
	$degenNtHsh_ref->{'A'} = 'A';
	$degenNtHsh_ref->{'C'} = 'C';
	$degenNtHsh_ref->{'G'} = 'G';
	$degenNtHsh_ref->{'T'} = 'T';
	$degenNtHsh_ref->{'U'} = 'U';
	$degenNtHsh_ref->{'R'} = 'AG';
	$degenNtHsh_ref->{'Y'} = 'CT';
	$degenNtHsh_ref->{'S'} = 'CG';
	$degenNtHsh_ref->{'W'} = 'AT';
	$degenNtHsh_ref->{'K'} = 'GT';
	$degenNtHsh_ref->{'M'} = 'AC';
	$degenNtHsh_ref->{'B'} = 'CGT';
	$degenNtHsh_ref->{'D'} = 'AGT';
	$degenNtHsh_ref->{'H'} = 'ACT';
	$degenNtHsh_ref->{'V'} = 'ACG';
	$degenNtHsh_ref->{'N'} = 'ACGT';
	
	return $degenNtHsh_ref;
}
sub defineGoldenmRNASet {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: currentTime|1071
#	appearInSub: predictBonaFideTSSinGoldenmRNASet|3628
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $distBtwATGTrnsfrgEndHsh_ref, $goldenmRNASetSizeLimit, $maxPctDistBtwATGTrnsfrgEnd
#	output: $distBtwATGTrnsfrgEndCutoff, $goldenmRNAListHsh_ref
#	toCall: my ($goldenmRNAListHsh_ref, $distBtwATGTrnsfrgEndCutoff) = &defineGoldenmRNASet($distBtwATGTrnsfrgEndHsh_ref, $maxPctDistBtwATGTrnsfrgEnd, $goldenmRNASetSizeLimit);
#	calledInLine: 3648
#....................................................................................................................................................#

	my ($distBtwATGTrnsfrgEndHsh_ref, $maxPctDistBtwATGTrnsfrgEnd, $goldenmRNASetSizeLimit) = @_;
	
	my $goldenmRNAListHsh_ref = {};
	
	print "[".&currentTime()."] Defining the golden mRNA set                         \n";#->1071
	my $distBtwATGTrnsfrgEndStatObj = Statistics::Descriptive::Full->new();
	foreach my $mRNAID (keys %{$distBtwATGTrnsfrgEndHsh_ref}) {
		$distBtwATGTrnsfrgEndStatObj->add_data($distBtwATGTrnsfrgEndHsh_ref->{$mRNAID});
	}

	my $distBtwATGTrnsfrgEndCutoff = $distBtwATGTrnsfrgEndStatObj->percentile($maxPctDistBtwATGTrnsfrgEnd);
	
	my @tmpListAry;
	
	foreach my $mRNAID (keys %{$distBtwATGTrnsfrgEndHsh_ref}) {
		push @tmpListAry, $mRNAID if $distBtwATGTrnsfrgEndHsh_ref->{$mRNAID} < $distBtwATGTrnsfrgEndCutoff; 
	}
	
	@tmpListAry = shuffle(@tmpListAry);
	my $goldenSetmRNASize = $goldenmRNASetSizeLimit;
	$goldenSetmRNASize = $#tmpListAry if $goldenSetmRNASize >= @tmpListAry;
	
	for my $i (0..$goldenSetmRNASize) {
		$goldenmRNAListHsh_ref->{$tmpListAry[$i]}++;
	}
	
	print "[".&currentTime()."] At percentile $maxPctDistBtwATGTrnsfrgEnd,  distBtwATGTrnsfrgEndCutoff=$distBtwATGTrnsfrgEndCutoff, with $goldenSetmRNASize gene selected                        \n";#->1071

	return ($goldenmRNAListHsh_ref, $distBtwATGTrnsfrgEndCutoff);
}
sub defineLibTypeInfo {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: currentTime|1071
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|273
#	secondaryAppearInSection: >none
#	input: $STDRead5EndPileupIndxPathAndSizeAry_ref, $TEXRead5EndPileupIndxPathAndSizeAry_ref, $resultWigDir
#	output: $libTypeInfoHsh_ref
#	toCall: my ($libTypeInfoHsh_ref) = &defineLibTypeInfo($STDRead5EndPileupIndxPathAndSizeAry_ref, $TEXRead5EndPileupIndxPathAndSizeAry_ref, $resultWigDir);
#	calledInLine: 279
#....................................................................................................................................................#

	my ($STDRead5EndPileupIndxPathAndSizeAry_ref, $TEXRead5EndPileupIndxPathAndSizeAry_ref, $resultWigDir) = @_;

	my $libTypeInfoHsh_ref = {};
	
	$libTypeInfoHsh_ref->{'STD'}{'read5EndPileupIndxPathAndSizeAry_ref'} = $STDRead5EndPileupIndxPathAndSizeAry_ref;
	$libTypeInfoHsh_ref->{'STD'}{'plus'}{'aryIndex'} = 0; #---the index on the pooled cntgCovAry
	$libTypeInfoHsh_ref->{'STD'}{'minus'}{'aryIndex'} = 1; #---the index on the pooled cntgCovAry

	$libTypeInfoHsh_ref->{'TEX'}{'read5EndPileupIndxPathAndSizeAry_ref'} = $TEXRead5EndPileupIndxPathAndSizeAry_ref;
	$libTypeInfoHsh_ref->{'TEX'}{'plus'}{'aryIndex'} = 2; #---the index on the pooled cntgCovAry
	$libTypeInfoHsh_ref->{'TEX'}{'minus'}{'aryIndex'} = 3; #---the index on the pooled cntgCovAry

	
	foreach my $TEXOrSTD ('TEX', 'STD') {

		#----get the size and path
		my $read5EndPileupIndxPathAry_ref = (); 
		my $sizeAry_ref = (); 
		
		foreach my $pathAndSize (@{$libTypeInfoHsh_ref->{$TEXOrSTD}{'read5EndPileupIndxPathAndSizeAry_ref'}}) {
			my ($path, $size) = split /,/, $pathAndSize;
			push @{$sizeAry_ref}, $size;
			push @{$read5EndPileupIndxPathAry_ref}, $path;
			print "[".&currentTime()."] $TEXOrSTD library size = $size         \n";#->1071
		}

		$libTypeInfoHsh_ref->{$TEXOrSTD}{'read5EndPileupIndxPathAry_ref'} = $read5EndPileupIndxPathAry_ref;
		$libTypeInfoHsh_ref->{$TEXOrSTD}{'sizeAry_ref'} = $sizeAry_ref;
	
		#----create the paths for wiggle
		foreach my $originalOrScaled ('original', 'scaled') {
			foreach my $plusOrMinus ('plus', 'minus') {
				$libTypeInfoHsh_ref->{$TEXOrSTD}{'wig'}{$originalOrScaled}{$plusOrMinus} = "$resultWigDir/$TEXOrSTD.$originalOrScaled.$plusOrMinus.wig";
			}
		}
	}
	
	#---calculate the scaling factor
	my $meanSize = (sum(@{$libTypeInfoHsh_ref->{'TEX'}{'sizeAry_ref'}}) + sum(@{$libTypeInfoHsh_ref->{'STD'}{'sizeAry_ref'}}))/2;
	
	foreach my $TEXOrSTD ('TEX', 'STD') {
		$libTypeInfoHsh_ref->{$TEXOrSTD}{'scalingFactor'} = $meanSize/sum(@{$libTypeInfoHsh_ref->{$TEXOrSTD}{'sizeAry_ref'}});
		print "[".&currentTime()."] $TEXOrSTD library scaling factor = $libTypeInfoHsh_ref->{$TEXOrSTD}{'scalingFactor'}         \n";#->1071
	}
	
	return ($libTypeInfoHsh_ref);

}
sub downSampleGoldenATGRd5EndDataHshForTesting {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: >none
#	appearInSub: predictBonaFideTSSinGoldenmRNASet|3628
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $goldenATGRd5EndDataHsh_ref
#	output: $filteredGoldenATGRd5EndDataHsh_ref
#	toCall: my ($filteredGoldenATGRd5EndDataHsh_ref) = &downSampleGoldenATGRd5EndDataHshForTesting($goldenATGRd5EndDataHsh_ref);
#	calledInLine: 3658
#....................................................................................................................................................#

	my ($goldenATGRd5EndDataHsh_ref) = @_;
	
	my $filterHsh_ref = {};
	
	#$filterHsh_ref->{'rngMin,rngMax'} = number;

	$filterHsh_ref->{'5-100'}{'num'} = 30;
	$filterHsh_ref->{'101-1000'}{'num'} = 30;
	$filterHsh_ref->{'1001-2000'}{'num'} = 30;
	$filterHsh_ref->{'2001-9999999'}{'num'} = 30;
	
	#-----push the mRNAID into different ranges of peak site cov
	foreach my $mRNAID (keys %{$goldenATGRd5EndDataHsh_ref}) {
		foreach my $pos (sort {$goldenATGRd5EndDataHsh_ref->{$mRNAID}{$b} <=> $goldenATGRd5EndDataHsh_ref->{$mRNAID}{$a}} keys %{$goldenATGRd5EndDataHsh_ref->{$mRNAID}}) {
			my $peakCov = $goldenATGRd5EndDataHsh_ref->{$mRNAID}{$pos};
			foreach my $rng (keys %{$filterHsh_ref}) {
				my ($rngMin, $rngMax) = split /-/, $rng;
				push @{$filterHsh_ref->{$rng}{'mRNAIDAry'}}, $mRNAID if $peakCov >= $rngMin and $peakCov <= $rngMax;
			}
			last;
		}
	}
	
	my $filteredGoldenATGRd5EndDataHsh_ref = {};
	foreach my $rng (keys %{$filterHsh_ref}) {
		my $storedNum = 0;
		foreach my $mRNAID (@{$filterHsh_ref->{$rng}{'mRNAIDAry'}}) {
			last if ($storedNum >= $filterHsh_ref->{$rng}{'num'});
			%{$filteredGoldenATGRd5EndDataHsh_ref->{$mRNAID}} = %{$goldenATGRd5EndDataHsh_ref->{$mRNAID}};
			$storedNum++;
		}
	}
	
	return $filteredGoldenATGRd5EndDataHsh_ref;
}
sub filterTEXSiteBasedOnTrainingCutoff {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: checkRunningThreadAndWaitToJoin|981, createEmptyGenomeCovPerlStorable|1010, generateThreadHshWithRandomItem|1645, getIndivCntgCovPlsPath|2146, reportStatus|4274
#	appearInSub: >none
#	primaryAppearInSection: 7_filterTEXTrackAfterTrainning|319
#	secondaryAppearInSection: >none
#	input: $TSSvsExonCutoffHsh_ref, $fastaHsh_ref, $filterTEXCovStorableDir, $libTypeInfoHsh_ref, $maxThread, $minOriginalCount, $rd5EndPlsInfoHsh_ref
#	output: $filterTEXCovPlsPathHsh_ref
#	toCall: my ($filterTEXCovPlsPathHsh_ref) = &filterTEXSiteBasedOnTrainingCutoff($TSSvsExonCutoffHsh_ref, $rd5EndPlsInfoHsh_ref, $libTypeInfoHsh_ref, $fastaHsh_ref, $filterTEXCovStorableDir, $minOriginalCount, $maxThread);
#	calledInLine: 324
#....................................................................................................................................................#

	my ($TSSvsExonCutoffHsh_ref, $rd5EndPlsInfoHsh_ref, $libTypeInfoHsh_ref, $fastaHsh_ref, $filterTEXCovStorableDir, $minOriginalCount, $maxThread) = @_;
	
	my $filterTEXCntgCovIdxHshPath = "$filterTEXCovStorableDir/index.hsh.pls";
	my $filterTEXCovPlsPathHsh_ref;
	
	if (-s $filterTEXCntgCovIdxHshPath) {
		
		&reportStatus("filterTEXCntgCovIdxHshPath found. skip filtering", 10, "\n");#->4274
		$filterTEXCovPlsPathHsh_ref = &getIndivCntgCovPlsPath($filterTEXCntgCovIdxHshPath);#->2146
		
	} else {

		#---generate empty storable for filterTEXC
		$filterTEXCntgCovIdxHshPath = &createEmptyGenomeCovPerlStorable($filterTEXCovStorableDir, $fastaHsh_ref);#->1010
		$filterTEXCovPlsPathHsh_ref = &getIndivCntgCovPlsPath($filterTEXCntgCovIdxHshPath);#->2146
	
		my $ratioCutoff = $TSSvsExonCutoffHsh_ref->{'ratio'};
		my $countCutoff = $TSSvsExonCutoffHsh_ref->{'TEXCountZeroSTD'};

		my $itemAry_ref = [(keys %{$rd5EndPlsInfoHsh_ref->{'original'}{'covPlsPathHsh_ref'}})];
		my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomItem($maxThread, $itemAry_ref);#->1645
		my $cntgProc :shared = 0;
		my $filterSite :shared = 0;
		my %threadHsh = ();
		my $totalRatioSite :shared = 0;
		my $totalTEXCountZeroSTDSite :shared = 0;
		my $validRatioSite :shared = 0;
		my $validTEXCountZeroSTDSite :shared = 0;
		
		foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
			my $cntgAry_ref = $randCntgInThreadHsh_ref->{$threadNum};
			($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781
				sub {

					my ($cntgAry_ref) = @_;

					foreach my $cntg (@{$cntgAry_ref}) {
						
						$cntgProc++;
						
						&reportStatus("Filtering TEX Cov with ratio > $ratioCutoff or count > $countCutoff in $cntgProc th cntg. $filterSite sites passed", 10, "\r");#->4274
			
						my $originalCntgCovAry_ref = retrieve($rd5EndPlsInfoHsh_ref->{'original'}{'covPlsPathHsh_ref'}->{$cntg});
						my $filterTEXCntgCovAry_ref = retrieve($filterTEXCovPlsPathHsh_ref->{$cntg});
		
						foreach my $i (0..$#{$originalCntgCovAry_ref}) {
							if ($originalCntgCovAry_ref->[$i]) {
								my $tmpFilterTEXCovHsh_ref = {};
								my @tmpOriginalCntgCovAry = split /,/, $originalCntgCovAry_ref->[$i];
								foreach my $strndInWord ('plus', 'minus') {
									my $TEXCov = $tmpOriginalCntgCovAry[$libTypeInfoHsh_ref->{'TEX'}{$strndInWord}{'aryIndex'}];
									my $STDCov = $tmpOriginalCntgCovAry[$libTypeInfoHsh_ref->{'STD'}{$strndInWord}{'aryIndex'}];
									$tmpFilterTEXCovHsh_ref->{$strndInWord} = 0;
									if ($TEXCov >= $minOriginalCount) {
						
										if ($STDCov >= $minOriginalCount) {#----case 1: with both TEX and STD
											$totalRatioSite++;
											if ($TEXCov >= $TSSvsExonCutoffHsh_ref->{'ratio'}*$STDCov) {
												$tmpFilterTEXCovHsh_ref->{$strndInWord} = $TEXCov;
												$filterSite++;
												$validRatioSite++;
											}
										} else {#----case 2: with only TEX
											$totalTEXCountZeroSTDSite++;
											if ($TEXCov >= $TSSvsExonCutoffHsh_ref->{'TEXCountZeroSTD'}) {
												$tmpFilterTEXCovHsh_ref->{$strndInWord} = $TEXCov;
												$filterSite++;
												$validTEXCountZeroSTDSite++;
											}
										}
						
									}
								}
								if ($tmpFilterTEXCovHsh_ref->{'plus'} > 0 or $tmpFilterTEXCovHsh_ref->{'minus'} > 0) {
									$filterTEXCntgCovAry_ref->[$i] = join ",", ($tmpFilterTEXCovHsh_ref->{'plus'}, $tmpFilterTEXCovHsh_ref->{'minus'});
								}
							}
						}
						store($filterTEXCntgCovAry_ref, $filterTEXCovPlsPathHsh_ref->{$cntg});
					}
					
					return ();
				}#---end of sub
				,($cntgAry_ref)
			); #---end of threads->new
		}#---end of foreach my $threadNum

		&checkRunningThreadAndWaitToJoin('no', 1);#->981
	
		print "\n";
		&reportStatus("ratioSite = $validRatioSite / $totalRatioSite", 10, "\n");#->4274
		&reportStatus("TEXCountZeroSTDSite = $validTEXCountZeroSTDSite / $totalTEXCountZeroSTDSite", 10, "\n");#->4274
		
	}
	

	return ($filterTEXCovPlsPathHsh_ref);
}
sub findDistanceBetweenATGAndTrnsfrgEnd {
#....................................................................................................................................................#
#	subroutineCategory: specific, gff
#	dependOnSub: checkOverlapAndProximity_withMargin|744, ggplotHistogram|2780, reportStatus|4274
#	appearInSub: predictBonaFideTSSinGoldenmRNASet|3628
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $ggplotDirHsh_ref, $mRNAInfoHsh_ref, $maxThread, $trnsfrgInfoHsh_ref
#	output: $distBtwATGTrnsfrgEndHsh_ref
#	toCall: my ($distBtwATGTrnsfrgEndHsh_ref) = &findDistanceBetweenATGAndTrnsfrgEnd($trnsfrgInfoHsh_ref, $mRNAInfoHsh_ref, $ggplotDirHsh_ref, $maxThread);
#	calledInLine: 3645
#....................................................................................................................................................#

	my ($trnsfrgInfoHsh_ref, $mRNAInfoHsh_ref, $ggplotDirHsh_ref, $maxThread) = @_;

	#----check the overlapping distance
	my $checkPrxmty = 'no';
	my $reportExactMatch = 'yes';
	my $refRngType = 'geneRng';
	my $qryRngType = 'geneRng';
	my $refMargin = 0;
	my $qryMargin = 0;
	my ($hitAndPrxmtyToTrnsfrgBymRNAHsh_ref, $hitAndPrxmtyTomRNAByTrnsfrgHsh_ref) = &checkOverlapAndProximity_withMargin($mRNAInfoHsh_ref, $trnsfrgInfoHsh_ref, $checkPrxmty, $reportExactMatch, 2, $refRngType, $qryRngType, $refMargin, $qryMargin);#->744
	my $distBtwATGTrnsfrgEndHsh_ref = {};
	my $maxATGUpStrmDist = 100;
	
	#----check the overlapping distance
	&reportStatus("Calculate distance between ATG and trnsfrg", 10, "\n");#->4274
	foreach my $mRNAID (sort keys %{$hitAndPrxmtyToTrnsfrgBymRNAHsh_ref->{'SS'}{'hit'}}) {
		my $strnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
		my ($RNAStart, $RNAEnd) = @{$mRNAInfoHsh_ref->{$mRNAID}{'RNARng'}};
		my @tmpTrnsfrgBoundAry = ();
		foreach my $trnsfrgID (sort keys %{$hitAndPrxmtyToTrnsfrgBymRNAHsh_ref->{'SS'}{'hit'}{$mRNAID}}) {
			my ($trnsfrgStart, $trnsfrgEnd) = @{$trnsfrgInfoHsh_ref->{$trnsfrgID}{'RNARng'}};
			push @tmpTrnsfrgBoundAry, ($trnsfrgStart, $trnsfrgEnd);
		}
		
		@tmpTrnsfrgBoundAry = sort {$a <=> $b} @tmpTrnsfrgBoundAry;
		my $distBtwATGTrnsfrgEnd;
		if ($strnd eq '+') {
			$distBtwATGTrnsfrgEnd = $RNAStart - $tmpTrnsfrgBoundAry[0];
		} elsif ($strnd eq '-') {
			$distBtwATGTrnsfrgEnd = $tmpTrnsfrgBoundAry[-1] - $RNAEnd;
		} else {
			die "GFF strand error.\n";
		}
		$distBtwATGTrnsfrgEndHsh_ref->{$mRNAID} = $distBtwATGTrnsfrgEnd if ($distBtwATGTrnsfrgEnd > -1*$maxATGUpStrmDist);
	}
	
	#----plot the distribution
	my @plotAry = map {$distBtwATGTrnsfrgEndHsh_ref->{$_}} keys %{$distBtwATGTrnsfrgEndHsh_ref};

	my $sub_ggplotDir = "distBtwATGTrnsfrgEndHist";
	system "mkdir -pm 777 $ggplotDirHsh_ref->{$_}/$sub_ggplotDir/" foreach (keys %{$ggplotDirHsh_ref});
	
	my $item = 'distBtwATGTrnsfrgEndHist';
	my $dataPath = "$ggplotDirHsh_ref->{'dat'}/$sub_ggplotDir/$item.dat";
	my $pdfPath = "$ggplotDirHsh_ref->{'pdf'}/$sub_ggplotDir/$item.pdf";
	my $RScriptPath = "$ggplotDirHsh_ref->{'R'}/$sub_ggplotDir/$item.R";
	my $logPath = "$ggplotDirHsh_ref->{'log'}/$sub_ggplotDir/$item.log";
	my $leftXAxisPercentileLimit = 'min';
	my $rightXAxisPercentileLimit = 95;
	my $binWidth = 1;
	my $xAxis = 'distBtwATGTrnsfrgEnd';
	my $dataPtMax = 99999999;
	my $extraArg = '';
	my $log2OrLinear = 'linear';
	my $height = 8;
	my $width = 16;
	
	&ggplotHistogram(\@plotAry, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftXAxisPercentileLimit, $rightXAxisPercentileLimit, $xAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear, $height, $width);#->2780
	
	return ($distBtwATGTrnsfrgEndHsh_ref);
}
sub gatherIGVXMLTrack {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: currentTime|1071
#	appearInSub: >none
#	primaryAppearInSection: 11_outputWiggleXML|376
#	secondaryAppearInSection: >none
#	input: $libTypeInfoHsh_ref, $trnsfrgGFFPathHsh_ref, $wigglePathHsh_ref
#	output: $IGVTrackInfoHsh_ref
#	toCall: my ($IGVTrackInfoHsh_ref) = &gatherIGVXMLTrack($libTypeInfoHsh_ref, $trnsfrgGFFPathHsh_ref, $wigglePathHsh_ref);
#	calledInLine: 385
#....................................................................................................................................................#
	
	my ($libTypeInfoHsh_ref, $trnsfrgGFFPathHsh_ref, $wigglePathHsh_ref) = @_;
	
	my $IGVTrackInfoHsh_ref = {};
	my $colorHsh_ref = {};
	$colorHsh_ref->{'plus'} = '255,153,153';
	$colorHsh_ref->{'minus'} = '153,153,255';
	
	print "[".&currentTime()."] Gathering track for IGV.                               \n";#->1071
	
	foreach my $TEXOrSTD (sort keys %{$libTypeInfoHsh_ref}) {
		#foreach my $originalOrScaled (sort keys %{$libTypeInfoHsh_ref->{$TEXOrSTD}{'wig'}}) {
		foreach my $originalOrScaled ('original') {
			foreach my $plusOrMinus (sort keys %{$libTypeInfoHsh_ref->{$TEXOrSTD}{'wig'}{$originalOrScaled}}) {
				my $path = $libTypeInfoHsh_ref->{$TEXOrSTD}{'wig'}{$originalOrScaled}{$plusOrMinus};
				$IGVTrackInfoHsh_ref->{$path}{'name'} = "$TEXOrSTD.$originalOrScaled.$plusOrMinus";
				$IGVTrackInfoHsh_ref->{$path}{'displayMode'} = "SQUISHED";
				$IGVTrackInfoHsh_ref->{$path}{'DataRange_type'} = "LOG";
				$IGVTrackInfoHsh_ref->{$path}{'autoScale'} = "true";
				$IGVTrackInfoHsh_ref->{$path}{'color'} = $colorHsh_ref->{$plusOrMinus};
				$IGVTrackInfoHsh_ref->{$path}{'windowFunction'} = 'none';
				$IGVTrackInfoHsh_ref->{$path}{'renderer'} = "BAR_CHART";
				$IGVTrackInfoHsh_ref->{$path}{'height'} = 20;
			}
		}
	}
	
	foreach my $plusOrMinus (keys %{$trnsfrgGFFPathHsh_ref}) {
		my $path = $trnsfrgGFFPathHsh_ref->{$plusOrMinus};
		$IGVTrackInfoHsh_ref->{$path}{'name'} = "trnsfrg.$plusOrMinus";
		$IGVTrackInfoHsh_ref->{$path}{'displayMode'} = "SQUISHED";
		$IGVTrackInfoHsh_ref->{$path}{'DataRange_type'} = "LINEAR";
		$IGVTrackInfoHsh_ref->{$path}{'autoScale'} = "true";
		$IGVTrackInfoHsh_ref->{$path}{'color'} = $colorHsh_ref->{$plusOrMinus};
		$IGVTrackInfoHsh_ref->{$path}{'windowFunction'} = 'none';
		$IGVTrackInfoHsh_ref->{$path}{'renderer'} = "GENE_TRACK";
		$IGVTrackInfoHsh_ref->{$path}{'height'} = 20;
	}
	
	foreach my $item (keys %{$wigglePathHsh_ref}) {
		foreach my $plusOrMinus (keys %{$wigglePathHsh_ref->{$item}}) {
			my $path = $wigglePathHsh_ref->{$item}{$plusOrMinus};
			$IGVTrackInfoHsh_ref->{$path}{'name'} = "$item.$plusOrMinus";
			$IGVTrackInfoHsh_ref->{$path}{'displayMode'} = "SQUISHED";
			$IGVTrackInfoHsh_ref->{$path}{'DataRange_type'} = "LOG";
			$IGVTrackInfoHsh_ref->{$path}{'autoScale'} = "true";
			$IGVTrackInfoHsh_ref->{$path}{'color'} = $colorHsh_ref->{$plusOrMinus};
			$IGVTrackInfoHsh_ref->{$path}{'windowFunction'} = 'none';
			$IGVTrackInfoHsh_ref->{$path}{'renderer'} = "BAR_CHART";
			$IGVTrackInfoHsh_ref->{$path}{'height'} = 20;
		}
	}
	
	return $IGVTrackInfoHsh_ref;
}
sub generateDREMECmd {
#....................................................................................................................................................#
#	subroutineCategory: thridPartyApp
#	dependOnSub: >none
#	appearInSub: searchMotifAroundmRNAReferencePoint|4318
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_analyzeTSSMotif|365
#	input: $dremeMode, $dremeOutDir, $maxk, $minE, $mink, $negLogPath, $negSeqPath, $negativeAlignHsh_ref, $posSeqPath, $positiveAlignHsh_ref, $shuffleLogPath
#	output: $dremeOutputHsh_ref
#	toCall: my ($dremeOutputHsh_ref) = &generateDREMECmd($positiveAlignHsh_ref, $negativeAlignHsh_ref, $posSeqPath, $negSeqPath, $dremeOutDir, $shuffleLogPath, $negLogPath, $maxk, $mink, $minE, $dremeMode);
#	calledInLine: 4401
#....................................................................................................................................................#

	my ($positiveAlignHsh_ref, $negativeAlignHsh_ref, $posSeqPath, $negSeqPath, $dremeOutDir, $shuffleLogPath, $negLogPath, $maxk, $mink, $minE, $dremeMode) = @_;
	
	my $tmpDataInfoHsh_ref = {};
	$tmpDataInfoHsh_ref->{'pos'}{'hsh_ref'} = $positiveAlignHsh_ref;
	$tmpDataInfoHsh_ref->{'pos'}{'path'} = $posSeqPath;
	$tmpDataInfoHsh_ref->{'neg'}{'hsh_ref'} = $negativeAlignHsh_ref;
	$tmpDataInfoHsh_ref->{'neg'}{'path'} = $negSeqPath;
	
	foreach my $posOrNeg (keys %{$tmpDataInfoHsh_ref}) {
		my $seqPath = $tmpDataInfoHsh_ref->{$posOrNeg}{'path'};
		my $seqAlignHsh_ref = $tmpDataInfoHsh_ref->{$posOrNeg}{'hsh_ref'};
		open (SEQ, ">", $seqPath);
		foreach my $seqName (keys %{$seqAlignHsh_ref}) {
			next if $seqAlignHsh_ref->{$seqName} =~ m/[^ATGCatgc]/;
			print SEQ ">".$seqName."\n";
			print SEQ $seqAlignHsh_ref->{$seqName}."\n";
		}
		close SEQ;
	}
	
	my $dremeOutputHsh_ref = {};
	my @shuffleOrNegativeAry = ();
	push @shuffleOrNegativeAry, 'shuffle' if $dremeMode eq 'shuffle' or $dremeMode eq 'both';
	push @shuffleOrNegativeAry, 'negative' if $dremeMode eq 'negative' or $dremeMode eq 'both';
	
	foreach my $shuffleOrNegative (@shuffleOrNegativeAry) {
		my $dremeOutSubDir = "$dremeOutDir/$shuffleOrNegative/"; system ("mkdir -pm 777 $dremeOutSubDir");
		$dremeOutputHsh_ref->{$shuffleOrNegative}{'dir'} = $dremeOutSubDir;
		my $cmd = "dreme -norc -oc $dremeOutSubDir -p $tmpDataInfoHsh_ref->{'pos'}{'path'} -maxk $maxk -mink $mink -e $minE";
		my $log = $shuffleLogPath;
		if ($shuffleOrNegative eq 'negative') {
			$log = $negLogPath;
			$cmd .= " -n $tmpDataInfoHsh_ref->{'neg'}{'path'} "
		}
		$cmd .= " &>$log";
		$dremeOutputHsh_ref->{$shuffleOrNegative}{'cmd'} = $cmd;
		$dremeOutputHsh_ref->{$shuffleOrNegative}{'xml'} = "$dremeOutSubDir/dreme.xml";
	}
	
	return ($dremeOutputHsh_ref);
}
sub generateMEMECmd {
#....................................................................................................................................................#
#	subroutineCategory: thridPartyApp
#	dependOnSub: >none
#	appearInSub: searchMotifAroundmRNAReferencePoint|4318
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_analyzeTSSMotif|365
#	input: $logPath, $maxk, $memeOutDir, $memeSeqLimit, $mink, $nonAnnotatedPolymerFreqPath, $seqAlignHsh_ref, $seqPath
#	output: $XMLPath, $cmd
#	toCall: my ($cmd, $XMLPath) = &generateMEMECmd($seqAlignHsh_ref, $seqPath, $memeOutDir, $memeSeqLimit, $mink, $maxk, $logPath, $nonAnnotatedPolymerFreqPath);
#	calledInLine: 4431
#....................................................................................................................................................#

	my ($seqAlignHsh_ref, $seqPath, $memeOutDir, $memeSeqLimit, $mink, $maxk, $logPath, $nonAnnotatedPolymerFreqPath) = @_;
	
	#---generate a random name list for sampling limited seq
	my @randomSeqNameAry = ();
	foreach my $seqName (keys %{$seqAlignHsh_ref}) {
		push @randomSeqNameAry, $seqName;
	}
	@randomSeqNameAry = shuffle(@randomSeqNameAry);
	$memeSeqLimit = $#randomSeqNameAry if $memeSeqLimit > $#randomSeqNameAry;
	
	open (SEQ, ">", $seqPath);
	foreach my $i (0..$memeSeqLimit) {
		my $seqName = $randomSeqNameAry[$i];
		print SEQ ">".$seqName."\n";
		print SEQ $seqAlignHsh_ref->{$seqName}."\n";
	}
	close SEQ;
	
	my $cmd = "meme $seqPath -p 10 -oc $memeOutDir -bfile $nonAnnotatedPolymerFreqPath -dna -minw $mink -maxw $maxk -nmotifs 3 &>$logPath";
	my $XMLPath = "$memeOutDir/meme.xml";
	
	return ($cmd, $XMLPath);
	
}
sub generateThreadHshWithRandomCntg {
#....................................................................................................................................................#
#	subroutineCategory: multithread
#	dependOnSub: >none
#	appearInSub: checkOverlapAndProximity_withMargin|744
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $cntgAry_ref, $threadToSpawn
#	output: $randCntgInThreadHsh_ref
#	toCall: my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomCntg($threadToSpawn, $cntgAry_ref);
#	calledInLine: 786
#....................................................................................................................................................#

	my ($threadToSpawn, $cntgAry_ref) = @_;

	my @shuffleCntgAry = shuffle(@{$cntgAry_ref});
	my $threadNum = 1;
	my $randCntgInThreadHsh_ref = {};
	foreach my $cntg (@{$cntgAry_ref}) {
		$threadNum = 1 if $threadNum > $threadToSpawn;
		push @{$randCntgInThreadHsh_ref->{$threadNum}}, $cntg;
		$threadNum++;
	}
	
	return $randCntgInThreadHsh_ref;

}
sub generateThreadHshWithRandomItem {
#....................................................................................................................................................#
#	subroutineCategory: multithread
#	dependOnSub: >none
#	appearInSub: filterTEXSiteBasedOnTrainingCutoff|1275, getCoverageOfItemRngType_multiStrand|1897, getRd5EndCountForExonAndTSSOfGoldenmRNASet|2451, storeAndScaleRd5EndFromBothLibraries|4452
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 7_filterTEXTrackAfterTrainning|319
#	input: $itemAry_ref, $maxThread
#	output: $randItemForThrHsh_ref
#	toCall: my ($randItemForThrHsh_ref) = &generateThreadHshWithRandomItem($maxThread, $itemAry_ref);
#	calledInLine: 1307, 1911, 2475, 4483
#....................................................................................................................................................#

	my ($maxThread, $itemAry_ref) = @_;

	my @shuffleItemAry = shuffle(@{$itemAry_ref});
	my $threadNum = 1;
	my $randItemForThrHsh_ref = {};
	foreach my $item (@shuffleItemAry) {
		$threadNum = 1 if $threadNum > $maxThread;
		push @{$randItemForThrHsh_ref->{$threadNum}}, $item;
		$threadNum++;
	}
	
	return ($randItemForThrHsh_ref);

}
sub getAndPrintTrnsfrg {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: currentTime|1071, printGFF_oneRNAPerGene_chooseStrnd_filterAry|3727
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|273
#	secondaryAppearInSection: >none
#	input: $trnsfrgGFFPathHsh_ref, $trnsfrgInfoHshStorablePath
#	output: $trnsfrgInfoHsh_ref
#	toCall: my ($trnsfrgInfoHsh_ref) = &getAndPrintTrnsfrg($trnsfrgInfoHshStorablePath, $trnsfrgGFFPathHsh_ref);
#	calledInLine: 293
#....................................................................................................................................................#

	my ($trnsfrgInfoHshStorablePath, $trnsfrgGFFPathHsh_ref) = @_;

	print "[".&currentTime()."] Retrieving trnsfrgInfoHshStorable                        \n";#->1071
	my $trnsfrgInfoHsh_ref = retrieve($trnsfrgInfoHshStorablePath);
	
	my $gffGeneLineOnly = 'yes';
	foreach my $strndInWord (keys %{$trnsfrgGFFPathHsh_ref}) {
		my $outGFFPath = $trnsfrgGFFPathHsh_ref->{$strndInWord};
		&printGFF_oneRNAPerGene_chooseStrnd_filterAry($trnsfrgInfoHsh_ref, $outGFFPath, $gffGeneLineOnly, $strndInWord, []);#->3727
	}
	
	return $trnsfrgInfoHsh_ref;

}
sub getBaseAtTSSAndExon {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: currentTime|1071, reportStatus|4274
#	appearInSub: >none
#	primaryAppearInSection: 9_analyzeNucleotideComposition|349
#	secondaryAppearInSection: >none
#	input: $TEXBaseCompositionPileupPathHsh_ref, $exonEndTrim, $fastaHsh_ref, $libTypeInfoHsh_ref, $mRNAInfoHsh_ref, $mRNARefPtHshPath, $rd5EndPlsInfoHsh_ref, $strndEndValidTSSInfoHsh_ref, $validTSSLimitHsh_ref
#	output: $mRNARefPtHsh_ref, $siteBaseGenomeProportionHsh_ref, $siteBaseStringWithReadHsh_ref
#	toCall: my ($siteBaseStringWithReadHsh_ref, $siteBaseGenomeProportionHsh_ref, $mRNARefPtHsh_ref) = &getBaseAtTSSAndExon($strndEndValidTSSInfoHsh_ref, $mRNAInfoHsh_ref, $TEXBaseCompositionPileupPathHsh_ref, $rd5EndPlsInfoHsh_ref, $libTypeInfoHsh_ref, $fastaHsh_ref, $exonEndTrim, $validTSSLimitHsh_ref, $mRNARefPtHshPath);
#	calledInLine: 355
#....................................................................................................................................................#

	my ($strndEndValidTSSInfoHsh_ref, $mRNAInfoHsh_ref, $TEXBaseCompositionPileupPathHsh_ref, $rd5EndPlsInfoHsh_ref, $libTypeInfoHsh_ref, $fastaHsh_ref, $exonEndTrim, $validTSSLimitHsh_ref, $mRNARefPtHshPath) = @_;
	
	my $randomExonSitePerGene = 50; #---number of random exon per gene to be sampled
	my $numberOfTSSPerGene = 1; #---number of TSS per gene to be sampled sorted by cov
	my $minCountToSample = 5; #----minium read count of the a site to be sampled 
	my $NAT5EndRegReltvPosLowerLimit = $validTSSLimitHsh_ref->{'TAA'}{'antisense'}{'lowerValLimit'};
	my $NAT5EndRegReltvPosUpperLimit = $validTSSLimitHsh_ref->{'TAA'}{'antisense'}{'upperValLimit'};
	my $mRNAUTR5ReltvPosLowerLimit = $validTSSLimitHsh_ref->{'ATG'}{'sense'}{'lowerValLimit'};
	my $mRNARefPtHsh_ref = {}; #----collect the mRNA_TSS, NAT_TSS, ATG and TAA as reference point, only the genes with NAT_TSS with collect TAA and only gene with mRNA_TSS with collect ATG
	
	my $siteTypeHsh_ref = {};
	$siteTypeHsh_ref->{'ATG'}{'sense'} = "mRNA_TSS";
	$siteTypeHsh_ref->{'TAA'}{'antisense'} = "NAT_TSS";
	
	&reportStatus("Picking valid TSS postion for base scanning", 40, "\n");#->4274
	
	#---get mRNA with TSS at ATG or TAA and the site with highest rd5Ends
	my $siteByCntgBymRNAHsh_ref = {};
	my $allTSSSiteHash_ref = {};
	foreach my $ATGOrTAA (keys %{$strndEndValidTSSInfoHsh_ref}) {
		foreach my $senseOrAntisense (keys %{$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}}) {
			my $siteType = $siteTypeHsh_ref->{$ATGOrTAA}{$senseOrAntisense};
			foreach my $mRNAID (keys %{$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}}) {
				my $cntg = $mRNAInfoHsh_ref->{$mRNAID}{'cntg'};

				#----sort the $rltvPos by cov
				my $storedTSSNum = 0;
				foreach my $rltvPos (sort {$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{$mRNAID}{$b}{'cov'} <=> $strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{$mRNAID}{$a}{'cov'}} keys %{$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{$mRNAID}}) {
					my $absPos = $strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{$mRNAID}{$rltvPos}{'absPos'};
					
					#---get only the TSS peak as reference point
					$mRNARefPtHsh_ref->{$siteType}{$mRNAID} = $absPos if not defined $mRNARefPtHsh_ref->{$siteType}{$mRNAID};
					if ($storedTSSNum < $numberOfTSSPerGene) {#--take the highest cov site
						push @{$siteByCntgBymRNAHsh_ref->{$cntg}{$mRNAID}{$siteType}}, $absPos;
						
						$storedTSSNum++;
					}
					$allTSSSiteHash_ref->{$cntg}{$mRNAID}{$absPos} = $siteType;
				}
				
				if ($siteType eq 'NAT_TSS') {
					if ($mRNAInfoHsh_ref->{$mRNAID}{'strnd'} eq '-') {
						$mRNARefPtHsh_ref->{'NAT_TSS'}{$mRNAID}--;
					} else {
						$mRNARefPtHsh_ref->{'NAT_TSS'}{$mRNAID}++;
					}
				}
			}
		}
	}

	&reportStatus("Picking random exon postion for base scanning", 40, "\n");#->4274

	#---get the exon sites of the genes with TSS
	foreach my $cntg (keys %{$siteByCntgBymRNAHsh_ref}) {
		foreach my $mRNAID (keys %{$siteByCntgBymRNAHsh_ref->{$cntg}}) {

			my $strnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
			my ($RNAStart, $RNAEnd) = @{$mRNAInfoHsh_ref->{$mRNAID}{'RNARng'}};
			my ($ATGPos, $TAAPos) = ($RNAStart, $RNAEnd);
			if ($strnd eq '-') {
				($ATGPos, $TAAPos) = ($TAAPos, $ATGPos);
				$TAAPos--; #---get the outside position as the reference
			} else {
				$TAAPos++; #---get the outside position as the reference
			}

			$mRNARefPtHsh_ref->{'mRNA_ATG'}{$mRNAID} = $ATGPos if $mRNARefPtHsh_ref->{'mRNA_TSS'}{$mRNAID};
			$mRNARefPtHsh_ref->{'mRNA_TAA'}{$mRNAID} = $TAAPos if $mRNARefPtHsh_ref->{'NAT_TSS'}{$mRNAID};

			#---get the exon position
			my @exonPosAry = ();

			@{$mRNAInfoHsh_ref->{$mRNAID}{'exonRng'}} = sort {$a <=> $b} @{$mRNAInfoHsh_ref->{$mRNAID}{'exonRng'}};
			for (my $i=0; $i < $#{$mRNAInfoHsh_ref->{$mRNAID}{'exonRng'}}; $i += 2) {
				my ($exonStart, $exonEnd) = ($mRNAInfoHsh_ref->{$mRNAID}{'exonRng'}->[$i], $mRNAInfoHsh_ref->{$mRNAID}{'exonRng'}->[$i+1]);
				push @exonPosAry, ($exonStart+$exonEndTrim..$exonEnd-$exonEndTrim);
			}

			my @shuffleIndexAry = shuffle(0..$#exonPosAry);
			foreach my $shuffledIndex (@shuffleIndexAry) {
				push @{$siteByCntgBymRNAHsh_ref->{$cntg}{$mRNAID}{'mRNA_exon'}}, $exonPosAry[$shuffledIndex];
				last if @{$siteByCntgBymRNAHsh_ref->{$cntg}{$mRNAID}{'mRNA_exon'}} >= $randomExonSitePerGene;
			}
			
			#----get the mRNA UTR5 and NAT end 5 reg sites
			my $tmpRngAryHsh_ref = {};
			if ($strnd eq '+') {
				@{$tmpRngAryHsh_ref->{'mRNA_UTR5'}} = ($RNAStart+$mRNAUTR5ReltvPosLowerLimit..$RNAStart);#----$mRNAUTR5ReltvPosLowerLimit is negetive
				@{$tmpRngAryHsh_ref->{'NAT_end5Reg'}} = ($RNAEnd+$NAT5EndRegReltvPosLowerLimit..$RNAEnd+$NAT5EndRegReltvPosUpperLimit);
			} else {
				@{$tmpRngAryHsh_ref->{'mRNA_UTR5'}} = ($RNAStart..$RNAStart-$mRNAUTR5ReltvPosLowerLimit);#----$mRNAUTR5ReltvPosLowerLimit is negetive
				@{$tmpRngAryHsh_ref->{'NAT_end5Reg'}} = ($RNAEnd-$NAT5EndRegReltvPosUpperLimit..$RNAEnd-$NAT5EndRegReltvPosLowerLimit);#----$NAT5EndRegReltvPosLowerLimit is negetive
			}
			
			#----get the site only if the it is not a TSS
			foreach my $mRNA_UTR5OrNAT_end5Reg (keys %{$tmpRngAryHsh_ref}) {
				foreach my $pos (@{$tmpRngAryHsh_ref->{$mRNA_UTR5OrNAT_end5Reg}}) {
					if (not $allTSSSiteHash_ref->{$cntg}{$mRNAID}{$pos}) {
						push @{$siteByCntgBymRNAHsh_ref->{$cntg}{$mRNAID}{$mRNA_UTR5OrNAT_end5Reg}}, $pos;
					}
				}
			}
		}
	}
	
	#---set up a hash to translate the strnd of the gene to the strnd of the site
	my $tmpStrndHsh_ref = {};
	$tmpStrndHsh_ref->{'+'}{'mRNA_UTR5'} = 'plus';
	$tmpStrndHsh_ref->{'+'}{'mRNA_TSS'} = 'plus';
	$tmpStrndHsh_ref->{'+'}{'mRNA_exon'} = 'plus';
	$tmpStrndHsh_ref->{'+'}{'NAT_TSS'} = 'minus';
	$tmpStrndHsh_ref->{'+'}{'NAT_end5Reg'} = 'minus';

	$tmpStrndHsh_ref->{'-'}{'mRNA_UTR5'} = 'minus';
	$tmpStrndHsh_ref->{'-'}{'mRNA_TSS'} = 'minus';
	$tmpStrndHsh_ref->{'-'}{'mRNA_exon'} = 'minus';
	$tmpStrndHsh_ref->{'-'}{'NAT_TSS'} = 'plus';
	$tmpStrndHsh_ref->{'-'}{'NAT_end5Reg'} = 'plus';
	
	my $siteBaseStringWithReadHsh_ref = {};
	my $siteBaseGenomeProportionHsh_ref = {};
	my $siteTypeCountHsh_ref = {};
	
	foreach my $cntg (keys %{$siteByCntgBymRNAHsh_ref}) {

		print "[".&currentTime()."] Scanning base composition of TSS/Exon on $cntg                          \r";#->1071

		my $originalCntgCovAry_ref = retrieve($rd5EndPlsInfoHsh_ref->{'original'}{'covPlsPathHsh_ref'}->{$cntg});
		my $baseComAry_ref = retrieve($TEXBaseCompositionPileupPathHsh_ref->{$cntg});
		my $cntgSeq = $fastaHsh_ref->{$cntg};
		my $cntgLen = length $cntgSeq;
		foreach my $mRNAID (keys %{$siteByCntgBymRNAHsh_ref->{$cntg}}) {
			my $geneStrnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
			foreach my $siteType (keys %{$siteByCntgBymRNAHsh_ref->{$cntg}{$mRNAID}}) {
				my $siteStrnd = $tmpStrndHsh_ref->{$geneStrnd}{$siteType};
				my $covAryIndex = $libTypeInfoHsh_ref->{'TEX'}{$siteStrnd}{'aryIndex'};
				foreach my $pos (@{$siteByCntgBymRNAHsh_ref->{$cntg}{$mRNAID}{$siteType}}) {
					#---set the index to pos - 1
					my $i = $pos - 1;
					next if $i > $cntgLen;
					my $genomeBase = substr $cntgSeq, $i, 1;
					$genomeBase =~ tr/ACGTacgt/TGCAtgca/ if $siteStrnd eq 'minus';

					#---collect the genome base with or without read
					$siteBaseGenomeProportionHsh_ref->{$siteType}{'count'}{$genomeBase}++;

					#----check whether there's read5End
					if ($originalCntgCovAry_ref->[$i]) {

						#---get the cov of the site in the TEX lib
						my @tmpCovAry = split /,/, $originalCntgCovAry_ref->[$i];
						my $cov = $tmpCovAry[$covAryIndex];
		
						if ($cov > $minCountToSample) {
							#----get the base combination in TEX library
							my $tmpBaseComHsh_ref = {};
							($tmpBaseComHsh_ref->{'plus'}, $tmpBaseComHsh_ref->{'minus'}) = split /,/, $baseComAry_ref->[$i];
							
							#----get the base in the genome sequence
							if ($genomeBase ne 'N') {
								$siteTypeCountHsh_ref->{$siteType}++;
								$siteBaseStringWithReadHsh_ref->{$siteType}{$cntg}{$pos}{'genomeBase'} = $genomeBase;
								$siteBaseStringWithReadHsh_ref->{$siteType}{$cntg}{$pos}{'readBase'} = $tmpBaseComHsh_ref->{$siteStrnd};
								#print TMPLOG join "\t", ($siteType, $cntg, $pos, $genomeBase, $tmpBaseComHsh_ref->{$siteStrnd}."\n"); #---debug
							}
						}
					}
				}
			}
		}
	}

	print "\n";
	
	store ($mRNARefPtHsh_ref, $mRNARefPtHshPath);
	
	#---report the count
	foreach my $siteType (keys %{$siteTypeCountHsh_ref}) {
		my $count = $siteTypeCountHsh_ref->{$siteType};
		print "[".&currentTime()."] $siteType collected = $count                          \n";#->1071
	}
	
	return ($siteBaseStringWithReadHsh_ref, $siteBaseGenomeProportionHsh_ref, $mRNARefPtHsh_ref);
}
sub getCoverageOfItemRngType_multiStrand {
#....................................................................................................................................................#
#	subroutineCategory: coverage
#	dependOnSub: generateThreadHshWithRandomItem|1645, reportStatus|4274
#	appearInSub: getPutativeTSSRegionBasedOnPeakValidSite|2364
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 8_investigateGeneTSSRelativePos|331
#	input: $dirtnAry_ref, $itemByCntgHsh_ref, $itemInfoHsh_ref, $margin, $maxThread, $pileupStorablePathHsh_ref, $rngType
#	output: $covHsh_ref
#	toCall: my ($covHsh_ref) = &getCoverageOfItemRngType_multiStrand($pileupStorablePathHsh_ref, $itemInfoHsh_ref, $itemByCntgHsh_ref, $maxThread, $rngType, $margin, $dirtnAry_ref);
#	calledInLine: 2384
#....................................................................................................................................................#
	my ($pileupStorablePathHsh_ref, $itemInfoHsh_ref, $itemByCntgHsh_ref, $maxThread, $rngType, $margin, $dirtnAry_ref) = @_;
	
	my $cntgAry_ref = [(keys %{$pileupStorablePathHsh_ref})];
	my $randCntgInThreadHsh_ref = &generateThreadHshWithRandomItem($maxThread, $cntgAry_ref);#->1645
	my $cntgProc :shared = 0;
	my %threadHsh = ();

	my %strndCovToGetHsh = ();

	$strndCovToGetHsh{'+'}{'s'} = '+';
	$strndCovToGetHsh{'+'}{'a'} = '-';
	$strndCovToGetHsh{'-'}{'s'} = '-';
	$strndCovToGetHsh{'-'}{'a'} = '+';

	foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
		my $cntgAry_ref = \@{$randCntgInThreadHsh_ref->{$threadNum}};
		my $cntgNum = @{$randCntgInThreadHsh_ref->{$threadNum}};
		&reportStatus("$cntgNum cntgs spawned to thread $threadNum", 20, "\r");#->4274

		#---spawn a new thread
		($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781
		
			sub {
				my ($cntgAry_ref) = @_;
				my $covInThrHsh_ref = {};
				foreach my $cntg (@{$cntgAry_ref}) {
					$cntgProc++;
					my $cntgCovAry_ref = retrieve($pileupStorablePathHsh_ref->{$cntg});
					next if not exists $itemByCntgHsh_ref->{$cntg};
					foreach my $itemID (keys %{$itemByCntgHsh_ref->{$cntg}}) {
						next if not $itemInfoHsh_ref->{$itemID}{$rngType};

						my @rng = sort {$a <=> $b} @{$itemInfoHsh_ref->{$itemID}{$rngType}};
						unshift @rng, ($rng[0]-1-$margin, $rng[0]-1);
						push @rng, ($rng[-1]+1, $rng[-1]+1+$margin);
	
						my $posCovAry_ref = ();
						for (my $i=0; $i < $#rng; $i += 2) {
							foreach my $j ($rng[$i]-1..$rng[$i+1]-1) {
								my %tmpCovHsh = ('+'=>0, '-'=>0);
								($tmpCovHsh{'+'}, $tmpCovHsh{'-'}) = split /,/, $cntgCovAry_ref->[$j] if ($cntgCovAry_ref->[$j]);
								foreach my $dirtn (@{$dirtnAry_ref}) {
									push @{$posCovAry_ref->{$dirtn}}, $tmpCovHsh{$strndCovToGetHsh{$itemInfoHsh_ref->{$itemID}{'strnd'}}{$dirtn}};
								}
							}
						}
	
						foreach my $dirtn (@{$dirtnAry_ref}) {
							if ($itemInfoHsh_ref->{$itemID}{'strnd'} eq '-') {
								@{$covInThrHsh_ref->{$itemID}{$dirtn}} = reverse @{$posCovAry_ref->{$dirtn}};
							} else {
								@{$covInThrHsh_ref->{$itemID}{$dirtn}} = @{$posCovAry_ref->{$dirtn}};
							}
						}
					}
					
					&reportStatus("Finished counting items on $cntgProc cntg", 20, "\r");#->4274
				}
				return ($covInThrHsh_ref);
			}
			,($cntgAry_ref)
		);
	}
	
	my $covHsh_ref = {};
	
	while (keys %threadHsh) {
		foreach my $threadNum (keys %threadHsh) {
			my $thr = $threadHsh{$threadNum};
			if (not $thr->is_running()) {
				my ($covInThrHsh_ref) = $thr->join;
				foreach my $itemID (keys %{$covInThrHsh_ref}) {
					foreach my $dirtn (keys %{$covInThrHsh_ref->{$itemID}}) {
						$covHsh_ref->{$itemID}{$dirtn} = $covInThrHsh_ref->{$itemID}{$dirtn};
					}
				}
				delete $threadHsh{$threadNum};
			}
		}
		sleep 1;
	}

	my $numItem = scalar(keys %{$covHsh_ref});
	my $dirtnStr = join " & ", @{$dirtnAry_ref};
	&reportStatus("$rngType coverage of $numItem items in $dirtnStr direction were stored", 10, "\n");#->4274
	
	return ($covHsh_ref);
}
sub getCtgryGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|1071
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|273
#	secondaryAppearInSection: >none
#	input: $ctgryAry_ref, $geneInfoHsh_ref
#	output: $geneCtgryByCntgHsh_ref, $geneCtgryInfoHsh_ref
#	toCall: my ($geneCtgryInfoHsh_ref, $geneCtgryByCntgHsh_ref) = &getCtgryGeneInfo($geneInfoHsh_ref, $ctgryAry_ref);
#	calledInLine: 290
#....................................................................................................................................................#

	my ($geneInfoHsh_ref, $ctgryAry_ref) = @_;
	
	my $geneCtgryInfoHsh_ref = {};
	my $geneCtgryByCntgHsh_ref = {};
	
	my $ctgryStr = join ",", @{$ctgryAry_ref};

	print "[".&currentTime()."] Filtering GFF on cgtry $ctgryStr.\n";#->1071
	
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		my $cntg = $geneInfoHsh_ref->{$geneID}{'cntg'};
		if (grep /^$ctgry$/, @{$ctgryAry_ref}) {
			%{$geneCtgryInfoHsh_ref->{$geneID}} = %{$geneInfoHsh_ref->{$geneID}};
			$geneCtgryByCntgHsh_ref->{$cntg}{$geneID}++;
		}
	}
	
	my $numGene = keys %{$geneCtgryInfoHsh_ref};
	
	print "[".&currentTime()."] $numGene gene filtered on cgtry $ctgryStr.\n";#->1071
	
	return $geneCtgryInfoHsh_ref, $geneCtgryByCntgHsh_ref;
}
sub getDREMEMotif {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: plotMotifOccurenceWithMASTMatch|3361
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_analyzeTSSMotif|365
#	input: $degenNtHsh_ref, $dremeXMLToRead
#	output: $motifInfoHsh_ref
#	toCall: my ($motifInfoHsh_ref) = &getDREMEMotif($dremeXMLToRead, $degenNtHsh_ref);
#	calledInLine: 3449
#....................................................................................................................................................#

	my ($dremeXMLToRead, $degenNtHsh_ref) = @_;
	
	my $motifInfoHsh_ref;
	my $motifNum = 0;
	open (DREMEXML, "$dremeXMLToRead");
	while (my $theLine = <DREMEXML>) {
		chomp $theLine;
		if ($theLine =~ m/<motif id=/) {
			$motifNum++;
			my @theLineSplt = split / |\>/, $theLine;
			my ($motifSeq, $evalue);
			foreach my $arg (@theLineSplt) {
				$arg =~ s/\"//g;
				if ($arg =~ m/^seq=/) {$motifSeq = substr ($arg, index ($arg, "=")+1);}
				elsif ($arg =~ m/^evalue=/) {$evalue = substr ($arg, index ($arg, "=")+1);}
			}
			
			#---generate the regular expression of the degenerate bases
			my @motifSeqAry = split //, $motifSeq;
			my $regexString = '';
			foreach my $degenNt (@motifSeqAry) {
				my @baseComAry = split //, $degenNtHsh_ref->{$degenNt};
				my $baseComWithComma = join ",", @baseComAry;
				$regexString .= "{$baseComWithComma}";
			}

			$motifInfoHsh_ref->{$motifSeq}{'regexString'} = $regexString;
			$motifInfoHsh_ref->{$motifSeq}{"evalue"} = $evalue;
			$motifInfoHsh_ref->{$motifSeq}{"motifNum"} = $motifNum;
			@{$motifInfoHsh_ref->{$motifSeq}{"matchSeq"}} = glob $regexString;

		}
=pod
		#---obsolete, will use regexStr and glob to generate matchSeq
		if ($theLine =~ m/<match seq=/) {
			my @theLineSplt = split / |\>/, $theLine;
			foreach my $arg (@theLineSplt) {
				if ($arg =~ m/^seq=/) {
					$arg =~ s/\"//g;
					my $matchSeq = substr ($arg, index ($arg, "=")+1); 
					push @{$motifInfoHsh_ref->{$motifSeq}{"matchSeq"}}, $matchSeq;
				}
			}
		}
=cut
	}
	close (DREMEXML);
	return $motifInfoHsh_ref;
}
sub getGeneWithValidTSSAtGeneEnd {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 8_investigateGeneTSSRelativePos|331
#	secondaryAppearInSection: >none
#	input: $TSSmRNAEndFreqHsh_ref, $geneBasedTSSInfoHsh_ref, $strndEndValidTSSInfoHshPath, $validTSSLimitHsh_ref
#	output: $strndEndValidTSSInfoHsh_ref
#	toCall: my ($strndEndValidTSSInfoHsh_ref) = &getGeneWithValidTSSAtGeneEnd($geneBasedTSSInfoHsh_ref, $TSSmRNAEndFreqHsh_ref, $strndEndValidTSSInfoHshPath, $validTSSLimitHsh_ref);
#	calledInLine: 342
#....................................................................................................................................................#

	my ($geneBasedTSSInfoHsh_ref, $TSSmRNAEndFreqHsh_ref, $strndEndValidTSSInfoHshPath, $validTSSLimitHsh_ref) = @_;
	#----set the Pct limit and count the number of genes with valid TSS
	
	my $strndEndValidTSSInfoHsh_ref = {};
	foreach my $ATGOrTAA (keys %{$validTSSLimitHsh_ref}) {
		foreach my $senseOrAntisense (keys %{$validTSSLimitHsh_ref->{$ATGOrTAA}}) {

			#----add the pct line to define TSS to be stored
			my $valueStatObj = Statistics::Descriptive::Full->new();
			$valueStatObj->add_data(@{$TSSmRNAEndFreqHsh_ref->{$ATGOrTAA}{$senseOrAntisense}});
	
			my $lowerPctLimit = $validTSSLimitHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{'lowerPctLimit'};
			my $upperPctLimit = $validTSSLimitHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{'upperPctLimit'};
	
			my $lowerValLimit = sprintf "%.3f", $valueStatObj->percentile($lowerPctLimit);
			my $upperValLimit = sprintf "%.3f", $valueStatObj->percentile($upperPctLimit);
			
			$validTSSLimitHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{'lowerValLimit'} = $lowerValLimit;
			$validTSSLimitHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{'upperValLimit'} = $upperValLimit;
			
			foreach my $mRNAID (keys %{$geneBasedTSSInfoHsh_ref}) {
				foreach my $i (0..$#{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{'rltvPos'}}) {
					my $rltvPos = ${$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{'rltvPos'}}[$i];
					if ($rltvPos <= $upperValLimit and $rltvPos >= $lowerValLimit) {
						$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{$mRNAID}{$rltvPos}{'cov'} = ${$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{'cov'}}[$i];
						$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{$mRNAID}{$rltvPos}{'absPos'} = ${$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{'absPos'}}[$i];
					}
				}
			}
		}
	}
	
	store ($strndEndValidTSSInfoHsh_ref, $strndEndValidTSSInfoHshPath);
	
	return ($strndEndValidTSSInfoHsh_ref);
	
}
sub getIndivCntgCovPlsPath {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: currentTime|1071
#	appearInSub: filterTEXSiteBasedOnTrainingCutoff|1275, poolLibRd5EndCov|3584, storeAndScaleRd5EndFromBothLibraries|4452
#	primaryAppearInSection: 5_processInputData|273
#	secondaryAppearInSection: 5_processInputData|273, 7_filterTEXTrackAfterTrainning|319
#	input: $cntgCovPlsIndexPath
#	output: $cntgCovInPlsPathHsh_ref
#	toCall: my ($cntgCovInPlsPathHsh_ref) = &getIndivCntgCovPlsPath($cntgCovPlsIndexPath);
#	calledInLine: 299, 302, 1295, 1301, 3611, 3619, 4470
#....................................................................................................................................................#
	
	my ($cntgCovPlsIndexPath) = @_;

	$cntgCovPlsIndexPath =~ s/\.gz$//;

	my $cntgCovInPlsPathHsh_ref = {};
	
	system ("gzip -fd $cntgCovPlsIndexPath.gz") if -s "$cntgCovPlsIndexPath.gz";
	my %plsIndexHsh = %{retrieve($cntgCovPlsIndexPath)};
	
	my (undef, $cntgCovStroableDir, undef) = fileparse($cntgCovPlsIndexPath, qr/\.[^.]*/);
	foreach my $cntg (keys %plsIndexHsh) {
		my $cntgCovPlsPath = "$cntgCovStroableDir/$plsIndexHsh{$cntg}";
		die "cntgCovPlsPath $cntgCovPlsPath is invalid\n" if ((not -s $cntgCovPlsPath) and (not -s $cntgCovPlsPath.".gz"));
		$cntgCovInPlsPathHsh_ref->{$cntg} = $cntgCovPlsPath;
	}
	my $numCntg = keys %{$cntgCovInPlsPathHsh_ref};
	print "[".&currentTime()."] pls path of $numCntg contig stored.\n";#->1071
	
	return $cntgCovInPlsPathHsh_ref;
}
sub getMASTLogPostionalData {
#....................................................................................................................................................#
#	subroutineCategory: thridPartyApp
#	dependOnSub: >none
#	appearInSub: plotMotifOccurenceWithMASTMatch|3361
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_analyzeTSSMotif|365
#	input: $mastHitLog, $maxHitPVal, $maxMotifEVal, $maxPos, $motifInfoHsh_ref, $totalSeqNum
#	output: $motifPctHsh_ref
#	toCall: my ($motifPctHsh_ref) = &getMASTLogPostionalData($mastHitLog, $motifInfoHsh_ref, $maxHitPVal, $maxMotifEVal, $maxPos, $totalSeqNum);
#	calledInLine: 3421, 3465
#....................................................................................................................................................#
	my ($mastHitLog, $motifInfoHsh_ref, $maxHitPVal, $maxMotifEVal, $maxPos, $totalSeqNum) = @_;
	
	my $motifPctHsh_ref = {};
	my $motifInSeqCountHsh_ref = {};

	my %motifCountHsh = ();
	
	#---get the motif seq and reference to their number and in MAST the motif is represented as +\d
	my %motifNumToMotifSeqCoversionHsh = ();
	foreach my $motifSeq (keys %{$motifInfoHsh_ref}) {
		my $motifNum = $motifInfoHsh_ref->{$motifSeq}{"motifNum"};
		my $evalue = $motifInfoHsh_ref->{$motifSeq}{"evalue"};
		if ($evalue <= $maxMotifEVal) {
			$motifNumToMotifSeqCoversionHsh{'+'.$motifNum} = $motifSeq;
		}
	}
	
	open MASTLOG, "<", $mastHitLog;
	while (<MASTLOG>) {
		next if $_ =~ m/^#/;
		chomp;
		my ($sequence_name, $motif, $hit_start, $hit_end, $score, $hit_pValue) = split / +/;
		if ($hit_pValue <= $maxHitPVal and exists $motifNumToMotifSeqCoversionHsh{$motif}) {
			$motifCountHsh{$motifNumToMotifSeqCoversionHsh{$motif}}{$hit_start}++;
			push @{$motifInSeqCountHsh_ref->{$motifNumToMotifSeqCoversionHsh{$motif}}{$sequence_name}}, $hit_start;
		}
	}
	close MASTLOG;
	
	foreach my $motifSeq (keys %motifCountHsh) {
		foreach my $pos (1..$maxPos) {
			$motifPctHsh_ref->{$motifSeq}{$pos} = 0;
			$motifPctHsh_ref->{$motifSeq}{$pos} = 100*$motifCountHsh{$motifSeq}{$pos}/$totalSeqNum if $motifCountHsh{$motifSeq}{$pos};
		}
	}

	return ($motifPctHsh_ref);
}
sub getMEMEMotif {
#....................................................................................................................................................#
#	subroutineCategory: thridPartyApp
#	dependOnSub: >none
#	appearInSub: plotMotifOccurenceWithMASTMatch|3361
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_analyzeTSSMotif|365
#	input: $MEMEXMLToRead, $degenNtHsh_ref
#	output: $motifInfoHsh_ref
#	toCall: my ($motifInfoHsh_ref) = &getMEMEMotif($MEMEXMLToRead, $degenNtHsh_ref);
#	calledInLine: 3407
#....................................................................................................................................................#

	my ($MEMEXMLToRead, $degenNtHsh_ref) = @_;
	
	my $motifInfoHsh_ref = {};
	my $regexEvalueHsh_ref = {};
	my %revDegenNtHsh = reverse %{$degenNtHsh_ref}; #----reverse the key value, now baseCombination is the key and the degenerate letter is the value

	my $motifID;

	open (MEMEXML, "$MEMEXMLToRead");
	while (my $theLine = <MEMEXML>) {
		chomp $theLine;
		if ($theLine =~ m/<motif id="([^\"]+)".+name="([^\"]+)".+e_value="([^\"]+)"/) {
			$motifID = $1;
			$regexEvalueHsh_ref->{$motifID}{"motifNum"} = $2;
			$regexEvalueHsh_ref->{$motifID}{"evalue"} = $3;
		}
		
		$regexEvalueHsh_ref->{$motifID}{"regex"} = $1 if ($theLine =~ m/^([\[\]ATGC]+)$/);
	}
	
	foreach my $motifID (keys %{$regexEvalueHsh_ref}) {
		my $matchedSeqHsh_ref = {};
		my $regex = $regexEvalueHsh_ref->{$motifID}{"regex"};
		my @regexAry = split //, $regex;
		my $position = 0;
		my $inBracket = 'no';
		foreach my $i (0..$#regexAry) {
			$inBracket = 'yes' if ($regexAry[$i] eq '[');
			$inBracket = 'no' if ($regexAry[$i] eq ']');
			if ($regexAry[$i] =~ m/A|T|G|C/) {
				push @{$matchedSeqHsh_ref->{$position}}, $regexAry[$i];
			}
			$position++ if $inBracket eq 'no';
		}

		my $regexString = '';
		my $motifSeq = '';
		foreach my $position (sort {$a <=> $b} keys %{$matchedSeqHsh_ref}) {
			my $baseComWithComma = join ",", sort {$a cmp $b} @{$matchedSeqHsh_ref->{$position}};
			my $baseComNoComma = join '', sort {$a cmp $b} @{$matchedSeqHsh_ref->{$position}};
			$regexString .= "{$baseComWithComma}";
			$motifSeq .= $revDegenNtHsh{$baseComNoComma};#---convert to degenerate base
		}
	
		$motifInfoHsh_ref->{$motifSeq}{'regexString'} = $regexString;
		$motifInfoHsh_ref->{$motifSeq}{'evalue'} = $regexEvalueHsh_ref->{$motifID}{"evalue"};
		$motifInfoHsh_ref->{$motifSeq}{'motifNum'} = $regexEvalueHsh_ref->{$motifID}{"motifNum"};
		@{$motifInfoHsh_ref->{$motifSeq}{"matchSeq"}} = glob $regexString;
		
	}
	close (MEMEXML);

	return $motifInfoHsh_ref;
}
sub getMaxCovAndPosAroundGeneEdge {
#....................................................................................................................................................#
#	subroutineCategory: coverage
#	dependOnSub: reportStatus|4274
#	appearInSub: getPutativeTSSRegionBasedOnPeakValidSite|2364
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 8_investigateGeneTSSRelativePos|331
#	input: $covHsh_ref, $dirtnAry_ref, $headOrTailAry_ref, $insideSrchDist, $refPtFromAryEndDist
#	output: $maxPosCovByItemHsh_ref, $maxPosDstrbtnHsh_ref, $validItemCountHsh_ref
#	toCall: my ($maxPosCovByItemHsh_ref, $maxPosDstrbtnHsh_ref, $validItemCountHsh_ref) = &getMaxCovAndPosAroundGeneEdge($covHsh_ref, $headOrTailAry_ref, $dirtnAry_ref, $refPtFromAryEndDist, $insideSrchDist);
#	calledInLine: 2391
#....................................................................................................................................................#
	my ($covHsh_ref, $headOrTailAry_ref, $dirtnAry_ref, $refPtFromAryEndDist, $insideSrchDist) = @_;
	
	#---$refPtFromAryEndDist refers to the distance between the reference point, e.g. start and end of CDSRng, to the end of the covAry, it usually equals to the "margin" used in &getCoverageOfItemRngType_multiStrand
	
	my $maxPosCovByItemHsh_ref = {};
	my $maxPosDstrbtnHsh_ref = {};
	my $validItemCountHsh_ref = {};

	foreach my $itemID (keys %{$covHsh_ref}) {
		foreach my $dirtn (@{$dirtnAry_ref}) {
			if (not $covHsh_ref->{$itemID}{$dirtn}) {
				&reportStatus("warning: $itemID does not have coverage in $dirtn direction", 10, "\n");#->4274
				next;
			}
			my $covAry_ref = $covHsh_ref->{$itemID}{$dirtn};
			my $srchAry_ref = [];
		
			foreach my $headOrTail (@{$headOrTailAry_ref}) {

				if ($headOrTail eq 'head') {
					$srchAry_ref = [@$covAry_ref[0..$refPtFromAryEndDist+$insideSrchDist-1]];
				} elsif ($headOrTail eq 'tail') {
					$srchAry_ref = [@$covAry_ref[$#{$covAry_ref}-$refPtFromAryEndDist-$insideSrchDist..$#{$covAry_ref}]];
					@{$srchAry_ref} = reverse @{$srchAry_ref};
				} else {
					die "DEBUG: headOrTail = $headOrTail. Invalid\n";
				}
		
				my $maxCov = max(@{$srchAry_ref}); #---http://stackoverflow.com/questions/2700302/how-do-i-get-a-slice-from-an-array-reference
				my $maxPos = -9999;

				if ($maxCov > 0) {
					my $maxIdx = 0;
					$srchAry_ref->[$maxIdx] > $srchAry_ref->[$_] or $maxIdx = $_ for 1 .. $#{$srchAry_ref};

					if ($headOrTail eq 'head') {
						$maxPos = $maxIdx - $refPtFromAryEndDist;
					} elsif ($headOrTail eq 'tail') {
						$maxPos = $refPtFromAryEndDist - $maxIdx;
					} else {
						die "DEBUG: headOrTail = $headOrTail. Invalid\n";
					}
					push @{$maxPosDstrbtnHsh_ref->{$dirtn}{$headOrTail}}, $maxPos;
					$validItemCountHsh_ref->{$dirtn}{$headOrTail}{'valid'}++;
				}

				$validItemCountHsh_ref->{$dirtn}{$headOrTail}{'total'}++;
				$maxPosCovByItemHsh_ref->{$itemID}{$dirtn}{$headOrTail}{'maxCov'} = $maxCov;
				$maxPosCovByItemHsh_ref->{$itemID}{$dirtn}{$headOrTail}{'maxPos'} = $maxPos;
				
			}
		}
	}

	return ($maxPosCovByItemHsh_ref, $maxPosDstrbtnHsh_ref, $validItemCountHsh_ref);
}
sub getPutativeTSSRegionBasedOnPeakValidSite {
#....................................................................................................................................................#
#	subroutineCategory: specific, coverage, gff
#	dependOnSub: getCoverageOfItemRngType_multiStrand|1897, getMaxCovAndPosAroundGeneEdge|2296, ggplotHistogram|2780, printMaxCovAndPosAroundGeneEdgeInfo|3955
#	appearInSub: >none
#	primaryAppearInSection: 8_investigateGeneTSSRelativePos|331
#	secondaryAppearInSection: >none
#	input: $TSSSearchRngHsh_ref, $filterTEXCovPlsPathHsh_ref, $fixPutativeTSSRegHsh_ref, $ggplotDirHsh_ref, $mRNAByCntgHsh_ref, $mRNAInfoHsh_ref, $maxThread, $resultLogDir
#	output: $putativeTSSRegHsh_ref, $validTSSDstrbtnSenseHsh_ref
#	toCall: my ($putativeTSSRegHsh_ref, $validTSSDstrbtnSenseHsh_ref) = &getPutativeTSSRegionBasedOnPeakValidSite($mRNAInfoHsh_ref, $mRNAByCntgHsh_ref, $filterTEXCovPlsPathHsh_ref, $TSSSearchRngHsh_ref, $ggplotDirHsh_ref, $resultLogDir, $maxThread, $fixPutativeTSSRegHsh_ref);
#	calledInLine: 337
#....................................................................................................................................................#
	my ($mRNAInfoHsh_ref, $mRNAByCntgHsh_ref, $filterTEXCovPlsPathHsh_ref, $TSSSearchRngHsh_ref, $ggplotDirHsh_ref, $resultLogDir, $maxThread, $fixPutativeTSSRegHsh_ref) = @_;
	
	#---get the coverage
	my $pileupStorablePathHsh_ref = $filterTEXCovPlsPathHsh_ref;
	my $itemInfoHsh_ref = $mRNAInfoHsh_ref;
	my $itemByCntgHsh_ref = $mRNAByCntgHsh_ref;
	my $rngType = 'CDSRng';
	my $margin = $TSSSearchRngHsh_ref->{'ATG'}{'up'};
	my $dirtnAry_ref = [qw/s/];
	my ($covHsh_ref) = &getCoverageOfItemRngType_multiStrand($pileupStorablePathHsh_ref, $itemInfoHsh_ref, $itemByCntgHsh_ref, $maxThread, $rngType, $margin, $dirtnAry_ref);#->1897
	
	
	#---get the peak position
	my $headOrTailAry_ref = [qw/head/];
	my $refPtFromAryEndDist = $margin;
	my $insideSrchDist = 0;
	my ($maxPosCovByItemHsh_ref, $maxPosDstrbtnHsh_ref, $validItemCountHsh_ref) = &getMaxCovAndPosAroundGeneEdge($covHsh_ref, $headOrTailAry_ref, $dirtnAry_ref, $refPtFromAryEndDist, $insideSrchDist);#->2296
	
	#---print the log for double check
	&printMaxCovAndPosAroundGeneEdgeInfo($maxPosCovByItemHsh_ref, $resultLogDir, "peak.valid.TSS.UTR5");#->3955
	
	#---plot the distribution
	my $sub_ggplotDir = "validTSSDstrbtn";
	system "mkdir -pm 777 $ggplotDirHsh_ref->{$_}/$sub_ggplotDir/" foreach (keys %{$ggplotDirHsh_ref});

	my $plotAry_ref = $maxPosDstrbtnHsh_ref->{'s'}{'head'};

	my $valueStatObj = Statistics::Descriptive::Full->new();
	$valueStatObj->add_data(@{$plotAry_ref});

	my $item = "peakSiteOnlyInUTR5";
	my $dataPath = "$ggplotDirHsh_ref->{'dat'}/$sub_ggplotDir/$item.dat";
	my $pdfPath = "$ggplotDirHsh_ref->{'pdf'}/$sub_ggplotDir/$item.pdf";
	my $RScriptPath = "$ggplotDirHsh_ref->{'R'}/$sub_ggplotDir/$item.R";
	my $logPath = "$ggplotDirHsh_ref->{'log'}/$sub_ggplotDir/$item.log";
	my $leftXAxisPercentileLimit = 1;
	my $rightXAxisPercentileLimit = 'max';
	my $binWidth = 1;
	my $xAxis = "peak.valid.TSS.within.$margin.upstream.ATG";
	my $dataPtMax = 99999999;
	my $log2OrLinear = 'linear';

	my $lowerPctLimit = 5;
	my $upperPctLimit = 99;
	my $lowerValLimit = $valueStatObj->percentile($lowerPctLimit);
	my $upperValLimit = $valueStatObj->percentile($upperPctLimit);

	if ($fixPutativeTSSRegHsh_ref->{'lowerLimit'} ne 'auto') {
		$lowerValLimit = $fixPutativeTSSRegHsh_ref->{'lowerLimit'};
		$lowerPctLimit = 'fixed';
	}

	if ($fixPutativeTSSRegHsh_ref->{'upperLimit'} ne 'auto') {
		$upperValLimit = $fixPutativeTSSRegHsh_ref->{'upperLimit'};
		$upperPctLimit = 'fixed';
	}
	
	my $extraArg = '';
	$extraArg .= " + geom_vline(xintercept=c($lowerValLimit), linetype=\"dotted\") + annotate(\"text\", x=$lowerValLimit, y=0, label=\"$lowerPctLimit\%\=$lowerValLimit\", vjust=-0.2, hjust=-0.1, angle=90)";
	$extraArg .= " + geom_vline(xintercept=c($upperValLimit), linetype=\"dotted\") + annotate(\"text\", x=$upperValLimit, y=0, label=\"$upperPctLimit\%\=$upperValLimit\", vjust=-0.2, hjust=-0.1, angle=90)";
	my $height = 8;
	my $width = 16;

	my ($plotValueAry_ref) = &ggplotHistogram($plotAry_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftXAxisPercentileLimit, $rightXAxisPercentileLimit, $xAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear, $height, $width);#->2780
	
	my $putativeTSSRegHsh_ref = {
		'lowerLimit' => $lowerValLimit,
		'upperLimit' => $upperValLimit,
	};
	
	#---save to latter use
	my $validTSSDstrbtnSenseHsh_ref = $covHsh_ref;
	
	return ($putativeTSSRegHsh_ref, $validTSSDstrbtnSenseHsh_ref);
}
sub getRd5EndCountForExonAndTSSOfGoldenmRNASet {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: generateThreadHshWithRandomItem|1645, reportStatus|4274
#	appearInSub: trainingUsingExonVsTssData|4548
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_trainningSTDTEXRatio|308
#	input: $exonEndTrim, $goldenmRNATSSResultHsh_ref, $libTypeInfoHsh_ref, $mRNAInfoHsh_ref, $maxThread, $minOriginalCount, $rd5EndPlsInfoHsh_ref
#	output: $rd5EndTSSExonCountHsh_ref
#	toCall: my ($rd5EndTSSExonCountHsh_ref) = &getRd5EndCountForExonAndTSSOfGoldenmRNASet($mRNAInfoHsh_ref, $goldenmRNATSSResultHsh_ref, $exonEndTrim, $libTypeInfoHsh_ref, $rd5EndPlsInfoHsh_ref, $minOriginalCount, $maxThread);
#	calledInLine: 4573
#....................................................................................................................................................#

	my ($mRNAInfoHsh_ref, $goldenmRNATSSResultHsh_ref, $exonEndTrim, $libTypeInfoHsh_ref, $rd5EndPlsInfoHsh_ref, $minOriginalCount, $maxThread) = @_;
	
	#---generate mRNAID by cntg hash
	my $mRNAIDByCntgHsh_ref = {};
	foreach my $mRNAID (keys %{$goldenmRNATSSResultHsh_ref}) {
		my $cntg = $mRNAInfoHsh_ref->{$mRNAID}{'cntg'};
		push @{$mRNAIDByCntgHsh_ref->{$cntg}}, $mRNAID;
	}

	my $exonPosHsh_ref = {};#---to contain the exon positions and 

	my $itemAry_ref = [(keys %{$mRNAIDByCntgHsh_ref})];
	my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomItem($maxThread, $itemAry_ref);#->1645
	my $cntgProc :shared = 0;
	my %threadHsh = ();
	
	foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
		my $cntgAry_ref = $randCntgInThreadHsh_ref->{$threadNum};
		($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781
			sub {

				my ($cntgAry_ref) = @_;
	
				my $rd5EndTSSExonCountHsh_inThr_ref = {};
				
				#---get the exon positions of the goldenmRNA set
				foreach my $cntg (@{$cntgAry_ref}) {
					
					$cntgProc++;
					
					&reportStatus("Getting TEX vs STD Read5End training data on $cntgProc th cntg", 10, "\r");#->4274

					my $cntgCovAryHsh_ref = {};
					foreach my $originalOrScaled (keys %{$rd5EndPlsInfoHsh_ref}) {
						$cntgCovAryHsh_ref->{$originalOrScaled} = retrieve($rd5EndPlsInfoHsh_ref->{$originalOrScaled}{'covPlsPathHsh_ref'}->{$cntg});
					}

					foreach my $mRNAID (@{$mRNAIDByCntgHsh_ref->{$cntg}}) {

						my $strnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
						my $strndInWord = 'plus'; $strndInWord = 'minus' if $strnd eq '-';
						my $posAryHsh_ref = {};
			
						#---get the tss pos
						foreach my $pos (keys %{$goldenmRNATSSResultHsh_ref->{$mRNAID}{'tss'}}) {
							push @{$posAryHsh_ref->{'tss'}}, $pos;
						}
		
						#---get the exon pos
						@{$mRNAInfoHsh_ref->{$mRNAID}{'exonRng'}} = sort {$a <=> $b} @{$mRNAInfoHsh_ref->{$mRNAID}{'exonRng'}};
						for (my $i=0; $i < $#{$mRNAInfoHsh_ref->{$mRNAID}{'exonRng'}}; $i += 2) {
							my ($exonStart, $exonEnd) = ($mRNAInfoHsh_ref->{$mRNAID}{'exonRng'}->[$i], $mRNAInfoHsh_ref->{$mRNAID}{'exonRng'}->[$i+1]);
							foreach my $pos ($exonStart+$exonEndTrim..$exonEnd-$exonEndTrim) {
								push @{$posAryHsh_ref->{'exon'}}, $pos;
							}
						}
			
						foreach my $tssOrExon(keys %{$posAryHsh_ref}) {
							foreach my $pos (@{$posAryHsh_ref->{$tssOrExon}}) {
								if ($cntgCovAryHsh_ref->{'original'}->[$pos-1]) {
									my $tmpCountHsh_ref = {};
									my $tmpRd5EndAryHsh_ref = {};
						
									#----get original, scaled of TEX ans STD in linear and log2
									foreach my $originalOrScaled (keys %{$cntgCovAryHsh_ref}) {
										@{$tmpRd5EndAryHsh_ref->{$originalOrScaled}} = split /,/, $cntgCovAryHsh_ref->{$originalOrScaled}->[$pos-1];
										foreach my $TEXOrSTD (keys %{$libTypeInfoHsh_ref}) {
											$tmpCountHsh_ref->{$originalOrScaled}{'linear'}{$TEXOrSTD} = @{$tmpRd5EndAryHsh_ref->{$originalOrScaled}}[$libTypeInfoHsh_ref->{$TEXOrSTD}{$strndInWord}{'aryIndex'}];
											if ($tmpCountHsh_ref->{$originalOrScaled}{'linear'}{$TEXOrSTD} > 0) {
												$tmpCountHsh_ref->{$originalOrScaled}{'log2'}{$TEXOrSTD} = log($tmpCountHsh_ref->{$originalOrScaled}{'linear'}{$TEXOrSTD})/log(2);
											} else {
												$tmpCountHsh_ref->{$originalOrScaled}{'log2'}{$TEXOrSTD} = 0;
											}
										}
									}

									if ($tmpCountHsh_ref->{'original'}{'linear'}{'TEX'} > $minOriginalCount or $tmpCountHsh_ref->{'original'}{'linear'}{'STD'} > $minOriginalCount) {
										foreach my $originalOrScaled (keys %{$tmpCountHsh_ref}) {
											foreach my $linearOrLog2 (keys %{$tmpCountHsh_ref->{$originalOrScaled}}) {
												my $covStr = $tmpCountHsh_ref->{$originalOrScaled}{$linearOrLog2}{'STD'}.",".$tmpCountHsh_ref->{$originalOrScaled}{$linearOrLog2}{'TEX'};
												push @{$rd5EndTSSExonCountHsh_inThr_ref->{$tssOrExon}{$originalOrScaled}{$linearOrLog2}}, $covStr;
												#print TMPLOG "$tssOrExon\t$strnd\t$pos\t$covStr\t$originalOrScaled\t$linearOrLog2\n";
											}
										}
									}
								}
							}
						}
					} #---end of foreach my $mRNAID (@{$mRNAIDByCntgHsh_ref->{$cntg}}) {
				}#---end of foreach $cntg (keys %{$mRNAIDByCntgHsh_ref}) {
				return ($rd5EndTSSExonCountHsh_inThr_ref);
			}#---end of sub
			,($cntgAry_ref)
		); #---end of threads->new
	}#---end of foreach my $threadNum
	
	my $rd5EndTSSExonCountHsh_ref = {};

	while (keys %threadHsh) {
		foreach my $threadNum (keys %threadHsh) {
			my $thr = $threadHsh{$threadNum};
			if (not $thr->is_running()) {
				my ($rd5EndTSSExonCountHsh_inThr_ref) = $thr->join;
				foreach my $tssOrExon (keys %{$rd5EndTSSExonCountHsh_inThr_ref}) {
					foreach my $originalOrScaled (keys %{$rd5EndTSSExonCountHsh_inThr_ref->{$tssOrExon}}) {
						foreach my $linearOrLog2 (keys %{$rd5EndTSSExonCountHsh_inThr_ref->{$tssOrExon}{$originalOrScaled}}) {
							foreach my $covStr (@{$rd5EndTSSExonCountHsh_inThr_ref->{$tssOrExon}{$originalOrScaled}{$linearOrLog2}}) {
								push @{$rd5EndTSSExonCountHsh_ref->{$tssOrExon}{$originalOrScaled}{$linearOrLog2}}, $covStr;
							}
						}
					}
				}
				delete $threadHsh{$threadNum};
			}
		}
		sleep 1;
	}

	return ($rd5EndTSSExonCountHsh_ref);
}
sub getRd5EndSurroundingmRNAATG {
#....................................................................................................................................................#
#	subroutineCategory: specific, gff, coverage
#	dependOnSub: reportStatus|4274
#	appearInSub: predictBonaFideTSSinGoldenmRNASet|3628
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $ATGDnStrmSrchRegSize, $ATGUpStrmSrchRegSize, $TEXOrSTD, $libTypeInfoHsh_ref, $mRNAByCntgHsh_ref, $mRNAFilterListHsh_ref, $mRNAInfoHsh_ref, $originalOrScaled, $rd5EndPlsInfoHsh_ref
#	output: $ATGRd5EndDataHsh_ref
#	toCall: my ($ATGRd5EndDataHsh_ref) = &getRd5EndSurroundingmRNAATG($mRNAInfoHsh_ref, $mRNAByCntgHsh_ref, $rd5EndPlsInfoHsh_ref, $ATGUpStrmSrchRegSize, $ATGDnStrmSrchRegSize, $libTypeInfoHsh_ref, $mRNAFilterListHsh_ref, $TEXOrSTD, $originalOrScaled);
#	calledInLine: 3655
#....................................................................................................................................................#
	
	my ($mRNAInfoHsh_ref, $mRNAByCntgHsh_ref, $rd5EndPlsInfoHsh_ref, $ATGUpStrmSrchRegSize, $ATGDnStrmSrchRegSize, $libTypeInfoHsh_ref, $mRNAFilterListHsh_ref, $TEXOrSTD, $originalOrScaled) = @_;

	my $ATGRd5EndDataHsh_ref = {};
	
	foreach my $cntg (sort keys %{$mRNAByCntgHsh_ref}) {
		my $tmpCntgAry_refHsh_ref = {};
		foreach my $originalOrScaled (keys %{$rd5EndPlsInfoHsh_ref}) {
			$tmpCntgAry_refHsh_ref->{$originalOrScaled} = retrieve($rd5EndPlsInfoHsh_ref->{$originalOrScaled}{'covPlsPathHsh_ref'}->{$cntg});
		}
		
		foreach my $mRNAID (sort keys %{$mRNAByCntgHsh_ref->{$cntg}}) {
			
			#---skip the mRNA if it is not in the filter
			next if not $mRNAFilterListHsh_ref->{$mRNAID};
			
			&reportStatus("Getting read 5End around ATG of $mRNAID", 10, "\r");#->4274

			my $tmpATGRd5EndDataHsh_ref = {};
			
			my $strnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
			
			my @tmpBoundAry = sort {$a <=> $b} @{$mRNAInfoHsh_ref->{$mRNAID}{'CDSRng'}};
			my @srchRngAry = ();
			my $wordStrnd;
			if ($strnd eq '+') {
				@srchRngAry = ($tmpBoundAry[0]-$ATGUpStrmSrchRegSize..$tmpBoundAry[0]+$ATGDnStrmSrchRegSize);
				$wordStrnd = 'plus';
			} else {
				@srchRngAry = ($tmpBoundAry[-1]-$ATGDnStrmSrchRegSize..$tmpBoundAry[-1]+$ATGUpStrmSrchRegSize);
				$wordStrnd = 'minus';
			}
			
			foreach my $pos (@srchRngAry) {
				my $i = $pos-1;
				if (defined $tmpCntgAry_refHsh_ref->{$originalOrScaled}->[$i]) {
					my @tmpRd5EndAry = split /,/, $tmpCntgAry_refHsh_ref->{$originalOrScaled}->[$i];
					my $cov = $tmpRd5EndAry[$libTypeInfoHsh_ref->{$TEXOrSTD}{$wordStrnd}{'aryIndex'}];
					$ATGRd5EndDataHsh_ref->{$mRNAID}{$pos} = $cov if $cov > 0;
				}
			}
		}
	}
	
	return ($ATGRd5EndDataHsh_ref);
}
sub getSequenceAroundmRNAReferencePoint {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: currentTime|1071
#	appearInSub: >none
#	primaryAppearInSection: 9_analyzeNucleotideComposition|349
#	secondaryAppearInSection: >none
#	input: $fastaHsh_ref, $mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $resultFastaDir
#	output: $seqAroundSiteHsh_ref
#	toCall: my ($seqAroundSiteHsh_ref) = &getSequenceAroundmRNAReferencePoint($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $fastaHsh_ref, $resultFastaDir);
#	calledInLine: 358
#....................................................................................................................................................#

	my ($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $fastaHsh_ref, $resultFastaDir) = @_;
	
	my $siteSeqRng = 100;
	my $seqAroundSiteHsh_ref = {};
	
	foreach my $siteType (keys %{$mRNARefPtHsh_ref}) {
		print "[".&currentTime()."] Getting sequences around $siteType                          \n";#->1071
		my $fastaPath = "$resultFastaDir/$siteType.surround.$siteSeqRng.nt.fasta";
		open (FASTA, ">", $fastaPath);
		foreach my $mRNAID (keys %{$mRNARefPtHsh_ref->{$siteType}}) {
			if ($mRNAInfoHsh_ref->{$mRNAID}{'strnd'} eq '+') {

				$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'} = substr $fastaHsh_ref->{$mRNAInfoHsh_ref->{$mRNAID}{'cntg'}}, $mRNARefPtHsh_ref->{$siteType}{$mRNAID}-$siteSeqRng-1, $siteSeqRng;
				$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'} = substr $fastaHsh_ref->{$mRNAInfoHsh_ref->{$mRNAID}{'cntg'}}, $mRNARefPtHsh_ref->{$siteType}{$mRNAID}-1, $siteSeqRng;

			} else {
				
				$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'} = substr $fastaHsh_ref->{$mRNAInfoHsh_ref->{$mRNAID}{'cntg'}}, $mRNARefPtHsh_ref->{$siteType}{$mRNAID}-1+1, $siteSeqRng;
				$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'} = substr $fastaHsh_ref->{$mRNAInfoHsh_ref->{$mRNAID}{'cntg'}}, $mRNARefPtHsh_ref->{$siteType}{$mRNAID}-$siteSeqRng-1+1, $siteSeqRng;
				
				foreach my $upStrmOrDnStrm (keys %{$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}}) {
					$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{$upStrmOrDnStrm} = reverse $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{$upStrmOrDnStrm};
					$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{$upStrmOrDnStrm} =~ tr/ACGTacgt/TGCAtgca/;
				}
			}
			
			#---reverse complement the seq if the purpose is to investigate NAT
			if ($siteType eq 'mRNA_TAA' or $siteType eq 'NAT_TSS') {
				$seqAroundSiteHsh_ref->{$siteType."_originalDirtn"}{$mRNAID}{'upStrm'} = $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'};
				$seqAroundSiteHsh_ref->{$siteType."_originalDirtn"}{$mRNAID}{'dnStrm'} = $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'};
				($seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'}, $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'}) = ($seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'}, $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'});
				foreach my $upStrmOrDnStrm (keys %{$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}}) {
					$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{$upStrmOrDnStrm} = reverse $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{$upStrmOrDnStrm};
					$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{$upStrmOrDnStrm} =~ tr/ACGTacgt/TGCAtgca/;
				}
			}
			
			print FASTA ">".$mRNAID."\n";
			print FASTA $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'}.$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'}."\n";
			
			#print TMPLOG join "\t", ($siteType, $mRNAInfoHsh_ref->{$mRNAID}{'strnd'}, $mRNARefPtHsh_ref->{$siteType}{$mRNAID}, $mRNAID, $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'}, $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'}."\n");
		}
		close FASTA;
	}
	
	return $seqAroundSiteHsh_ref;
}
sub getSiteBaseMismatchAndProportion {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: currentTime|1071
#	appearInSub: >none
#	primaryAppearInSection: 9_analyzeNucleotideComposition|349
#	secondaryAppearInSection: >none
#	input: $siteBaseGenomeProportionHsh_ref, $siteBaseStringWithReadHsh_ref
#	output: $siteBaseWithReadMismatchHsh_ref, $siteBaseWithReadProportionHsh_ref
#	toCall: my ($siteBaseWithReadProportionHsh_ref, $siteBaseWithReadMismatchHsh_ref) = &getSiteBaseMismatchAndProportion($siteBaseStringWithReadHsh_ref, $siteBaseGenomeProportionHsh_ref);
#	calledInLine: 357
#....................................................................................................................................................#

	my ($siteBaseStringWithReadHsh_ref, $siteBaseGenomeProportionHsh_ref) = @_;
	
	my $minMismatchProportion = 0.6;

	my $siteBaseWithReadProportionHsh_ref = {};
	my $siteBaseWithReadMismatchHsh_ref = {};
	
	#---get the compositon
	foreach my $siteType (sort keys %{$siteBaseStringWithReadHsh_ref}) {

		print "[".&currentTime()."] Analyzing $siteType sites                          \n";#->1071

		#---get the count
		foreach my $cntg (sort keys %{$siteBaseStringWithReadHsh_ref->{$siteType}}) {
			foreach my $pos (sort keys %{$siteBaseStringWithReadHsh_ref->{$siteType}{$cntg}}) {
				my $genomeBase = $siteBaseStringWithReadHsh_ref->{$siteType}{$cntg}{$pos}{'genomeBase'};
				my $readBaseHsh_ref = {};
				($readBaseHsh_ref->{'A'}, $readBaseHsh_ref->{'T'}, $readBaseHsh_ref->{'G'}, $readBaseHsh_ref->{'C'}) = split /:/, $siteBaseStringWithReadHsh_ref->{$siteType}{$cntg}{$pos}{'readBase'};

				#---get the total count for calculating mismatch proportion
				my $totalReadCount = 0;
				foreach my $readBase (qw /A T G C/) {
					$totalReadCount += $readBaseHsh_ref->{$readBase};
				}
				$siteBaseWithReadProportionHsh_ref->{$siteType}{'count'}{'genomeBase'}{$genomeBase}++;

				#----with only get the base with highest count
				foreach my $readBase (sort {$readBaseHsh_ref->{$b} <=> $readBaseHsh_ref->{$a}} keys %{$readBaseHsh_ref}) {
					$siteBaseWithReadProportionHsh_ref->{$siteType}{'count'}{'readBase'}{$readBase}++;

					#---record the mismatch, if it is a mismatch (of course), and proportion > minMismatchProportion
					if ($readBase ne $genomeBase and $readBaseHsh_ref->{$readBase}/$totalReadCount > $minMismatchProportion) {
						$siteBaseWithReadMismatchHsh_ref->{$siteType}{'mismatch'}{$genomeBase}{$readBase}++;
					}
					last;
				}
			}
		}
		
		#---get the compositon with read
		foreach my $genomeBaseOrReadBase (sort keys %{$siteBaseWithReadProportionHsh_ref->{$siteType}{'count'}}) {
			my $totalCount = 0;
			foreach my $base (sort keys %{$siteBaseWithReadProportionHsh_ref->{$siteType}{'count'}{$genomeBaseOrReadBase}}) {
				$totalCount += $siteBaseWithReadProportionHsh_ref->{$siteType}{'count'}{$genomeBaseOrReadBase}{$base};
			}
			foreach my $base (sort keys %{$siteBaseWithReadProportionHsh_ref->{$siteType}{'count'}{$genomeBaseOrReadBase}}) {
				my $proportion = sprintf "%.5f", $siteBaseWithReadProportionHsh_ref->{$siteType}{'count'}{$genomeBaseOrReadBase}{$base}/$totalCount;
				$siteBaseWithReadProportionHsh_ref->{$siteType}{'proportion'}{$genomeBaseOrReadBase}{$base} = $proportion;
				print TMPLOG join "\t", ($siteType, $genomeBaseOrReadBase, $base, $proportion."\n");#---debug
			}
		}
		
		my $totalGenomeCount = 0;
		foreach my $base (sort keys %{$siteBaseGenomeProportionHsh_ref->{$siteType}{'count'}}) {
			$totalGenomeCount += $siteBaseGenomeProportionHsh_ref->{$siteType}{'count'}{$base};
		}
		foreach my $base (sort keys %{$siteBaseGenomeProportionHsh_ref->{$siteType}{'count'}}) {
			my $proportion = sprintf "%.5f", $siteBaseGenomeProportionHsh_ref->{$siteType}{'count'}{$base}/$totalGenomeCount;
			$siteBaseGenomeProportionHsh_ref->{$siteType}{'proportion'}{$base} = $proportion;
			print TMPLOG join "\t", ($siteType, 'withOrWithRead', $base, $proportion."\n");#---debug
		}
	}
	
	return $siteBaseWithReadProportionHsh_ref, $siteBaseWithReadMismatchHsh_ref;
}
sub ggplotHistogram {
#....................................................................................................................................................#
#	subroutineCategory: ggplot, plotInR
#	dependOnSub: >none
#	appearInSub: findDistanceBetweenATGAndTrnsfrgEnd|1386, getPutativeTSSRegionBasedOnPeakValidSite|2364, plotTSSRelativeToATGAndTAAHistogram|3494
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 8_investigateGeneTSSRelativePos|331
#	input: $RScriptPath, $binWidth, $dataPath, $dataPtMax, $extraArg, $height, $leftXAxisPercentileLimit, $log2OrLinear, $logPath, $pdfPath, $plotAry_ref, $rightXAxisPercentileLimit, $width, $xAxis
#	output: $plotValueAry_ref
#	toCall: my ($plotValueAry_ref) = &ggplotHistogram($plotAry_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftXAxisPercentileLimit, $rightXAxisPercentileLimit, $xAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear, $height, $width);
#	calledInLine: 1455, 2438, 3541
#....................................................................................................................................................#
	
	my ($plotAry_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftXAxisPercentileLimit, $rightXAxisPercentileLimit, $xAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear, $height, $width) = @_;
	
	my $valueStatObj = Statistics::Descriptive::Full->new();
	$valueStatObj->add_data(@{$plotAry_ref});

	my $leftXAxisLimitValue;
	if ($leftXAxisPercentileLimit eq 'min') {
		$leftXAxisLimitValue = $valueStatObj->min();
	} else {
		$leftXAxisLimitValue = $valueStatObj->percentile($leftXAxisPercentileLimit);
	}

	my $rightXAxisLimitValue;
	if ($rightXAxisPercentileLimit eq 'max') {
		$rightXAxisLimitValue = $valueStatObj->max();
	} else {
		$rightXAxisLimitValue = $valueStatObj->percentile($rightXAxisPercentileLimit);
	}
	
	#---trim the end values
	my @trimmedAry = ();
	foreach my $value (@{$plotAry_ref}) {
		my $transformedValue = $value;
		$transformedValue = log($value)/log(2) if $log2OrLinear eq 'log2';
		push @trimmedAry, $transformedValue if $value <= $rightXAxisLimitValue and $value >= $leftXAxisLimitValue;
	}
	
	$dataPtMax = @trimmedAry if $dataPtMax > @trimmedAry;

	#---down sample the data point number
	my @shuffleIndexAry = shuffle(0..$#trimmedAry);
	my $plotValueAry_ref = ();
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA $xAxis."\n";
	foreach my $i (0..$dataPtMax-1) {
		my $shuffleValue = $trimmedAry[$shuffleIndexAry[$i]];
		push @{$plotValueAry_ref}, $shuffleValue;
		print PLOTDATA $shuffleValue."\n";
		
	}
	close PLOTDATA;
	
	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "ggplot(dataFrame, aes(x=$xAxis)) + ggtitle(\"Distribution of $xAxis $log2OrLinear scale [n=$dataPtMax]\") + geom_histogram(binwidth=$binWidth, aes(y = ..density.., fill = ..count..)) + geom_density() $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\", height=$height, width=$width)\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath ");

	return $plotValueAry_ref;
	
}
sub ggplotTwoSampleHistogram {
#....................................................................................................................................................#
#	subroutineCategory: ggplot, plotInR
#	dependOnSub: >none
#	appearInSub: analyzeSTDTEXRd5EndRatio|489
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $RScriptPath, $binWidth, $dataPath, $dataPtMax, $extraArg, $leftXAxisPercentileLimit, $log2OrLinear, $logPath, $pdfPath, $plotAryHsh_ref, $rightXAxisPercentileLimit, $xAxis
#	output: none
#	toCall: &ggplotTwoSampleHistogram($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftXAxisPercentileLimit, $rightXAxisPercentileLimit, $xAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear);
#	calledInLine: 574
#....................................................................................................................................................#

	my ($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftXAxisPercentileLimit, $rightXAxisPercentileLimit, $xAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA "sample\t$xAxis\n";

	foreach my $sample (keys %{$plotAryHsh_ref}) {
		
		my $plotAry_ref = $plotAryHsh_ref->{$sample};
	
		my $valueStatObj = Statistics::Descriptive::Full->new();
		$valueStatObj->add_data(@{$plotAry_ref});

		my $leftXAxisLimitValue;
		if ($leftXAxisPercentileLimit eq 'min') {
			$leftXAxisLimitValue = $valueStatObj->min();
		} else {
			$leftXAxisLimitValue = $valueStatObj->percentile($leftXAxisPercentileLimit);
		}

		my $rightXAxisLimitValue;
		if ($rightXAxisPercentileLimit eq 'max') {
			$rightXAxisLimitValue = $valueStatObj->max();
		} else {
			$rightXAxisLimitValue = $valueStatObj->percentile($rightXAxisPercentileLimit);
		}
	
		#---trim the end values
		my @trimmedAry = ();
		foreach my $value (@{$plotAry_ref}) {
			my $transformedValue = $value;
			$transformedValue = log($value)/log(2) if $log2OrLinear eq 'log2';
			push @trimmedAry, $transformedValue if $value <= $rightXAxisLimitValue and $value >= $leftXAxisLimitValue;
		}
		
		my $indivDataPtMax = $dataPtMax;
		$indivDataPtMax = @trimmedAry if $indivDataPtMax > @trimmedAry;

		#---down sample the data point number
		my @shuffleIndexAry = shuffle(0..$#trimmedAry);
		foreach my $i (0..$indivDataPtMax-1) {
			my $shuffleValue = $trimmedAry[$shuffleIndexAry[$i]];
			print PLOTDATA "$sample\t$shuffleValue\n";
		}
	
	}
	close PLOTDATA;

	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "ggplot(dataFrame, aes($xAxis, fill = sample)) + ggtitle(\"Distribution of $xAxis $log2OrLinear scale\") + geom_density(alpha = 0.2) $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\")\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath ");

}
sub ggplotXYErrorBarPlot {
#....................................................................................................................................................#
#	subroutineCategory: ggplot, plotInR
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $RScriptPath, $XAxis, $XYErrorBarPlotHsh_ref, $YAxis, $dataPath, $logPath, $pdfPath, $plotDataPoint
#	output: none
#	toCall: &ggplotXYErrorBarPlot($XAxis, $YAxis, $XYErrorBarPlotHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $plotDataPoint);
#	calledInLine: none
#....................................................................................................................................................#

	my ($XAxis, $YAxis, $XYErrorBarPlotHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $plotDataPoint) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA join "", ((join "\t", ($XAxis, $YAxis)), "\n");
	foreach my $XVal (sort {$a <=> $b} keys %{$XYErrorBarPlotHsh_ref}) {
		my $YValAry_ref = $XYErrorBarPlotHsh_ref->{$XVal};
		foreach my $YVal (sort {$a <=> $b} @{$YValAry_ref}) {
			print PLOTDATA join "", ((join "\t", ($XVal, $YVal)), "\n");
		}
	}
	close PLOTDATA;
	
	my $geom_point = ''; 
	$geom_point = '+ geom_point()' if $plotDataPoint eq 'yes';
	
	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "stderr <- function(x){sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))}"."\n";
	print R "lowStderr <- function(x){return(mean(x)-stderr(x))}"."\n";
	print R "highStderr <- function(x){return(mean(x)+stderr(x))}"."\n";
	print R "ggplot(dataFrame, aes(x=$XAxis, y=$YAxis)) + stat_smooth(method = \"lm\", fill = \"grey\", formula = y ~ x, level = 0.999, fullrange=TRUE) + stat_smooth(method = \"lm\", fill = \"darkgrey\", formula = y ~ x, level = 0.99, fullrange=TRUE) + stat_smooth(method = \"lm\", fill = \"blue\", formula = y ~ x, level = 0.95, fullrange=TRUE) + stat_smooth(method = \"lm\", fill = \"purple\", formula = y ~ x, level = 0.90, fullrange=TRUE) + stat_smooth(method = \"lm\", fill = \"red\", formula = y ~ x, level = 0.80, fullrange=TRUE) + stat_smooth(method = \"lm\", fill = \"yellow\", formula = y ~ x, level = 0.60, fullrange=TRUE) $geom_point + stat_summary(fun.y=mean, geom=\"point\", size=2) + stat_summary(fun.ymin=lowStderr, fun.ymax=highStderr, geom=\"errorbar\")"."\n";
	print R "ggsave(file=\"$pdfPath\")\n";
	print R "fitModel <- lm($YAxis ~ $XAxis, data = dataFrame)"."\n";
	print R "summary(fitModel)"."\n";
	print R "confint(fitModel)"."\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath ");
}
sub ggplotXYLinesMultipleSamples {
#....................................................................................................................................................#
#	subroutineCategory: ggplot, plotInR
#	dependOnSub: >none
#	appearInSub: plotBaseCompositionAroundmRNAReferencePoint|3302, plotMotifOccurenceWithMASTMatch|3361
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_analyzeTSSMotif|365, 9_analyzeNucleotideComposition|349
#	input: $RScriptPath, $XAXis, $YAxis, $YVariable, $dataPath, $extraArg, $height, $logPath, $pdfPath, $plotDataHsh_ref, $width
#	output: none
#	toCall: &ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $YVariable, $extraArg, $height, $width);
#	calledInLine: 3356, 3439, 3483
#....................................................................................................................................................#

	my ($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $YVariable, $extraArg, $height, $width) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA join "", (join "\t", ($YVariable, $YAxis, $XAXis)), "\n";
	foreach my $YCategory (sort keys %{$plotDataHsh_ref}) {
		foreach my $XVal (sort {$a <=> $b} keys %{$plotDataHsh_ref->{$YCategory}}) {
			my $YVal = $plotDataHsh_ref->{$YCategory}{$XVal};
			print PLOTDATA join "", (join "\t", ($YCategory, $YVal, $XVal)), "\n";
		}
	}
	close PLOTDATA;

	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "ggplot(dataFrame, aes(x=$XAXis, y=$YAxis, colour=$YVariable)) + geom_line() $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\", height=$height, width=$width)\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath");

}
sub ggplotXYScatterPlot {
#....................................................................................................................................................#
#	subroutineCategory: ggplot, plotInR
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $RScriptPath, $XAxis, $XYPairScatterPlotAry_ref, $YAxis, $dataPath, $logPath, $pdfPath
#	output: none
#	toCall: &ggplotXYScatterPlot($XAxis, $YAxis, $XYPairScatterPlotAry_ref, $dataPath, $pdfPath, $RScriptPath, $logPath);
#	calledInLine: none
#....................................................................................................................................................#

	my ($XAxis, $YAxis, $XYPairScatterPlotAry_ref, $dataPath, $pdfPath, $RScriptPath, $logPath) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA join "", ((join "\t", ($XAxis, $YAxis)), "\n");
	foreach my $YXPair (@{$XYPairScatterPlotAry_ref}) {
		my ($XVal, $YVal) = split /,/, $YXPair;
		print PLOTDATA join "", ((join "\t", ($XVal, $YVal)), "\n");
	}
	close PLOTDATA;
	
	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "ggplot(dataFrame, aes(x=$XAxis, y=$YAxis)) + stat_smooth(method = \"lm\", fill = \"grey\", formula = y ~ x, level = 0.999, fullrange=TRUE) + stat_smooth(method = \"lm\", fill = \"darkgrey\", formula = y ~ x, level = 0.99, fullrange=TRUE) + stat_smooth(method = \"lm\", fill = \"blue\", formula = y ~ x, level = 0.95, fullrange=TRUE) + stat_smooth(method = \"lm\", fill = \"purple\", formula = y ~ x, level = 0.90, fullrange=TRUE) + stat_smooth(method = \"lm\", fill = \"red\", formula = y ~ x, level = 0.80, fullrange=TRUE) + stat_smooth(method = \"lm\", fill = \"yellow\", formula = y ~ x, level = 0.60, fullrange=TRUE) + geom_point()"."\n";
	print R "ggsave(file=\"$pdfPath\")\n";
	print R "fitModel <- lm($YAxis ~ $XAxis, data = dataFrame)"."\n";
	print R "summary(fitModel)"."\n";
	print R "confint(fitModel)"."\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath ");
}
sub identifiyPeakSiteWithinRegion {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|4274
#	appearInSub: predictBonaFideTSSinGoldenmRNASet|3628
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $goldenATGRd5EndDataHsh_ref, $mRNAInfoHsh_ref
#	output: $goldenmRNATSSResultHsh_ref
#	toCall: my ($goldenmRNATSSResultHsh_ref) = &identifiyPeakSiteWithinRegion($goldenATGRd5EndDataHsh_ref, $mRNAInfoHsh_ref);
#	calledInLine: 3660
#....................................................................................................................................................#

	my ($goldenATGRd5EndDataHsh_ref, $mRNAInfoHsh_ref) = @_;
	
	&reportStatus("Identify peak sites as TSS", 10, "\n");#->4274
	
	my $goldenmRNATSSResultHsh_ref;

	foreach my $mRNAID (keys %{$goldenATGRd5EndDataHsh_ref}) {
		$goldenmRNATSSResultHsh_ref->{$mRNAID}{'cntg'} = $mRNAInfoHsh_ref->{$mRNAID}{'cntg'};
		$goldenmRNATSSResultHsh_ref->{$mRNAID}{'strnd'} = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
		foreach my $pos (sort {$goldenATGRd5EndDataHsh_ref->{$mRNAID}{$b} <=> $goldenATGRd5EndDataHsh_ref->{$mRNAID}{$a}} keys %{$goldenATGRd5EndDataHsh_ref->{$mRNAID}}) {
			$goldenmRNATSSResultHsh_ref->{$mRNAID}{'tss'}{$pos} = $goldenATGRd5EndDataHsh_ref->{$mRNAID}{$pos};
			last;
		}
	}

	return ($goldenmRNATSSResultHsh_ref);
}
sub intervalizeXYPairsForErrorBarPlot {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $XValueInterval, $XYPairAry_ref, $maxIntervalSample, $minIntervalSample
#	output: $XYErrorBarPlotHsh_ref
#	toCall: my ($XYErrorBarPlotHsh_ref) = &intervalizeXYPairsForErrorBarPlot($XYPairAry_ref, $XValueInterval, $maxIntervalSample, $minIntervalSample);
#	calledInLine: none
#....................................................................................................................................................#

	my ($XYPairAry_ref, $XValueInterval, $maxIntervalSample, $minIntervalSample) = @_;
	
	#---construct a hash for sorting by XVal
	my $pairID = 0;
	my $tmpValuePairByIDHsh_ref = {};
	my $XValMin = 999999999999;
	foreach my $valuePair (@{$XYPairAry_ref}) {
		my ($XVal, $YVal) = split /,/, $valuePair;
		$tmpValuePairByIDHsh_ref->{$pairID}{'X'} = $XVal;
		$tmpValuePairByIDHsh_ref->{$pairID}{'Y'} = $YVal;
		$XValMin = $XVal if ($XVal < $XValMin);
		$pairID++;
	}
	
	#---divide the pairs into intervals
	my $intervalLimit = $XValMin + $XValueInterval;
	my $intervalNumber = 0;
	my $intervalValuesHsh_ref = {};
	foreach my $pairID (sort {$tmpValuePairByIDHsh_ref->{$a}{'X'} <=> $tmpValuePairByIDHsh_ref->{$b}{'X'}} keys %{$tmpValuePairByIDHsh_ref}) {
		my $XVal = $tmpValuePairByIDHsh_ref->{$pairID}{'X'};
		my $YVal = $tmpValuePairByIDHsh_ref->{$pairID}{'Y'};
		if ($XVal <= $intervalLimit) {#----small than interval limit
			push @{$intervalValuesHsh_ref->{$intervalNumber}{'X'}}, $XVal;
			push @{$intervalValuesHsh_ref->{$intervalNumber}{'Y'}}, $YVal;
		} else {#----greater the interval limit, increse the limit
			$intervalNumber++;
			$intervalLimit += $XValueInterval;
		}
	}
	
	#---down sample the pairs and mean the X
	my $XYErrorBarPlotHsh_ref;
	
	foreach my $intervalNumber (keys %{$intervalValuesHsh_ref}) {
		my $pairNum = $#{$intervalValuesHsh_ref->{$intervalNumber}{'X'}};
		next if $pairNum < $minIntervalSample;
		my $pairToSample = $maxIntervalSample;
		$pairToSample = $pairNum if $maxIntervalSample > $pairNum;
		
		my @shuffleIndexAry = shuffle(0..$pairNum);
		my $XValAry_ref = ();
		my $YValAry_ref = ();
		foreach my $i (1..$pairToSample) {
			my $shuffledIndex = $shuffleIndexAry[$i];
			push @{$XValAry_ref}, ${$intervalValuesHsh_ref->{$intervalNumber}{'X'}}[$shuffledIndex];
			push @{$YValAry_ref}, ${$intervalValuesHsh_ref->{$intervalNumber}{'Y'}}[$shuffledIndex];
		}
		
		my $XValMean = sum(@{$XValAry_ref})/@{$XValAry_ref};
		
		$XYErrorBarPlotHsh_ref->{$XValMean} = $YValAry_ref;
	}
	
	$intervalValuesHsh_ref = {};
	$tmpValuePairByIDHsh_ref = {};

	return $XYErrorBarPlotHsh_ref;
}
sub investigateTSSRelativeToATGAndTAA {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|4274
#	appearInSub: >none
#	primaryAppearInSection: 8_investigateGeneTSSRelativePos|331
#	secondaryAppearInSection: >none
#	input: $TSSSearchRngHsh_ref, $filterTEXCovPlsPathHsh_ref, $geneBasedTSSInfoHshPath, $mRNAInfoHsh_ref
#	output: $TSSmRNAEndFreqHsh_ref, $geneBasedTSSInfoHsh_ref
#	toCall: my ($geneBasedTSSInfoHsh_ref, $TSSmRNAEndFreqHsh_ref) = &investigateTSSRelativeToATGAndTAA($filterTEXCovPlsPathHsh_ref, $mRNAInfoHsh_ref, $TSSSearchRngHsh_ref, $geneBasedTSSInfoHshPath);
#	calledInLine: 341
#....................................................................................................................................................#

	my ($filterTEXCovPlsPathHsh_ref, $mRNAInfoHsh_ref, $TSSSearchRngHsh_ref, $geneBasedTSSInfoHshPath) = @_;
	
	#---generate mRNAID by cntg hash
	my $mRNAIDByCntgHsh_ref = {};
	foreach my $mRNAID (keys %{$mRNAInfoHsh_ref}) {
		my $cntg = $mRNAInfoHsh_ref->{$mRNAID}{'cntg'};
		push @{$mRNAIDByCntgHsh_ref->{$cntg}}, $mRNAID;
	}

	my $TSSmRNAEndFreqHsh_ref = {};
	my $geneBasedTSSInfoHsh_ref = {};

	#---get the exon positions of the goldenmRNA set
	foreach my $cntg (keys %{$mRNAIDByCntgHsh_ref}) {
	
		&reportStatus("Getting TSS frequency around ATG and TAA", 40, "\r");#->4274

		my $filterTEXCntgCovAry_ref = retrieve($filterTEXCovPlsPathHsh_ref->{$cntg});

		foreach my $mRNAID (@{$mRNAIDByCntgHsh_ref->{$cntg}}) {
			my $strnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
			
			#---assume $RNAStart = ATG and RNAEnd = TAA
			my ($RNAStart, $RNAEnd) = @{$mRNAInfoHsh_ref->{$mRNAID}{'RNARng'}};
			my $tmpPosSrchAryHsh_ref = {};
			
			foreach my $rltvPos (-1*$TSSSearchRngHsh_ref->{'ATG'}{'up'}..$TSSSearchRngHsh_ref->{'ATG'}{'dn'}) {
				if ($strnd eq '+') {
					$tmpPosSrchAryHsh_ref->{'ATG'}{$rltvPos} = $RNAStart+$rltvPos;
				} else {
					$tmpPosSrchAryHsh_ref->{'ATG'}{$rltvPos} = $RNAEnd-$rltvPos;
				}
			}

			foreach my $rltvPos (-1*$TSSSearchRngHsh_ref->{'TAA'}{'up'}..$TSSSearchRngHsh_ref->{'TAA'}{'dn'}) {
				if ($strnd eq '+') {
					$tmpPosSrchAryHsh_ref->{'TAA'}{$rltvPos} = $RNAEnd+$rltvPos;
				} else {
					$tmpPosSrchAryHsh_ref->{'TAA'}{$rltvPos} = $RNAStart-$rltvPos;
				}
			}
			
			foreach my $ATGOrTAA (keys %{$tmpPosSrchAryHsh_ref}) {
				
				foreach my $senseOrAntisense ('sense', 'antisense') {
					foreach my $rltvPosOrCov ('rltvPos', 'cov') {
						@{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{$rltvPosOrCov}} = ();
					}
				}
	
				foreach my $rltvPos (keys %{$tmpPosSrchAryHsh_ref->{$ATGOrTAA}}) {
					my $absPos = $tmpPosSrchAryHsh_ref->{$ATGOrTAA}{$rltvPos};

					#---pos > 0 to ensure the search is not out of bound
					if ($absPos > 0 and $filterTEXCntgCovAry_ref->[$absPos-1]) {
						my $i = $absPos-1;
						my ($senseCov, $antisenseCov) = split /,/, $filterTEXCntgCovAry_ref->[$i];
						($antisenseCov, $senseCov) = ($senseCov, $antisenseCov) if $strnd eq '-';

						if ($senseCov > 0) {
							push @{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{'sense'}{'rltvPos'}}, $rltvPos;
							push @{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{'sense'}{'cov'}}, $senseCov;
							push @{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{'sense'}{'absPos'}}, $absPos;
							push @{$TSSmRNAEndFreqHsh_ref->{$ATGOrTAA}{'sense'}}, $rltvPos;
						}

						if ($antisenseCov > 0) {
							push @{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{'antisense'}{'rltvPos'}}, $rltvPos;
							push @{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{'antisense'}{'cov'}}, $antisenseCov;
							push @{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{'antisense'}{'absPos'}}, $absPos;
							push @{$TSSmRNAEndFreqHsh_ref->{$ATGOrTAA}{'antisense'}}, $rltvPos;
						}
					}
				}
				
				foreach my $senseOrAntisense ('sense', 'antisense') {
					if (@{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{'rltvPos'}} == 0) {
						@{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{'rltvPos'}} = (-9999);
						@{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{'cov'}} = (-9999);
					}
				}
			}
		}
	}
	
	store ($geneBasedTSSInfoHsh_ref, $geneBasedTSSInfoHshPath);
	
	return $geneBasedTSSInfoHsh_ref, $TSSmRNAEndFreqHsh_ref;
}
sub outputOriginalAndScaledRd5EndWiggle {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: checkRunningThreadAndWaitToJoin|981, printWiggleSingleTrackFromCntgCovPlsPathHsh|3997, reportStatus|4274
#	appearInSub: >none
#	primaryAppearInSection: 11_outputWiggleXML|376
#	secondaryAppearInSection: >none
#	input: $forcePoolLibRd5End, $libTypeInfoHsh_ref, $rd5EndPlsInfoHsh_ref
#	output: none
#	toCall: &outputOriginalAndScaledRd5EndWiggle($libTypeInfoHsh_ref, $rd5EndPlsInfoHsh_ref, $forcePoolLibRd5End);
#	calledInLine: 382
#....................................................................................................................................................#

	my ($libTypeInfoHsh_ref, $rd5EndPlsInfoHsh_ref, $forcePoolLibRd5End) = @_;
	
	foreach my $TEXOrSTD (sort keys %{$libTypeInfoHsh_ref}) {
		foreach my $originalOrScaled (sort keys %{$libTypeInfoHsh_ref->{$TEXOrSTD}{'wig'}}) {
			my $covPlsPathHsh_ref = $rd5EndPlsInfoHsh_ref->{'original'}{'covPlsPathHsh_ref'};
			foreach my $plusOrMinus (sort keys %{$libTypeInfoHsh_ref->{$TEXOrSTD}{'wig'}{$originalOrScaled}}) {
				my $wigPath = $libTypeInfoHsh_ref->{$TEXOrSTD}{'wig'}{$originalOrScaled}{$plusOrMinus};
				if ($forcePoolLibRd5End eq 'yes' or not -s $wigPath) {
					&reportStatus("Output wiggle to $TEXOrSTD $originalOrScaled $plusOrMinus in multiple threads", 10, "\n");#->4274
					my $aryIndex = $libTypeInfoHsh_ref->{$TEXOrSTD}{$plusOrMinus}{'aryIndex'};
					threads->create(\&printWiggleSingleTrackFromCntgCovPlsPathHsh, ($covPlsPathHsh_ref, $aryIndex, $wigPath, 'no'));#->3997
				} else {
					&reportStatus("Wiggle to $TEXOrSTD $originalOrScaled $plusOrMinus found. skip outputting", 10, "\n");#->4274
				}
			}
		}
	}

	&checkRunningThreadAndWaitToJoin('no', 1);#->981
}
sub outputTSSResultWiggle {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: predictBonaFideTSSinGoldenmRNASet|3628
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $goldenTSSResultHsh_ref, $goldenTSSWigPathHsh_ref
#	output: none
#	toCall: &outputTSSResultWiggle($goldenTSSWigPathHsh_ref, $goldenTSSResultHsh_ref);
#	calledInLine: 3666
#....................................................................................................................................................#

	my ($goldenTSSWigPathHsh_ref, $goldenTSSResultHsh_ref) = @_;
	
	my $posByCntgHsh_ref = {};
	foreach my $regionID (keys %{$goldenTSSResultHsh_ref}) {
		my $cntg = $goldenTSSResultHsh_ref->{$regionID}{'cntg'};
		my $strnd = $goldenTSSResultHsh_ref->{$regionID}{'strnd'};
		foreach my $pos (keys %{$goldenTSSResultHsh_ref->{$regionID}{'tss'}}) {
			$posByCntgHsh_ref->{$strnd}{$cntg}{$pos} = $goldenTSSResultHsh_ref->{$regionID}{'tss'}{$pos};
		}
	}

	foreach my $strnd ('+','-') {
		open WIGGLE, ">", $goldenTSSWigPathHsh_ref->{$strnd};
		foreach my $cntg (sort keys %{$posByCntgHsh_ref->{$strnd}}) {
			print WIGGLE "variableStep chrom=$cntg span=1\n";
			foreach my $pos (sort keys %{$posByCntgHsh_ref->{$strnd}{$cntg}}) {
				print WIGGLE join '', ((join "\t", ($pos, $posByCntgHsh_ref->{$strnd}{$cntg}{$pos})), "\n");
			}
		}
		close WIGGLE;
	}
}
sub plotBaseCompositionAroundmRNAReferencePoint {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: calculateBaseCompositionInAlignments|581, ggplotXYLinesMultipleSamples|2961, reportStatus|4274
#	appearInSub: >none
#	primaryAppearInSection: 9_analyzeNucleotideComposition|349
#	secondaryAppearInSection: >none
#	input: $ggplotDirHsh_ref, $seqAroundSiteHsh_ref, $trimBaseComPlotRngHsh_ref
#	output: none
#	toCall: &plotBaseCompositionAroundmRNAReferencePoint($seqAroundSiteHsh_ref, $ggplotDirHsh_ref, $trimBaseComPlotRngHsh_ref);
#	calledInLine: 360
#....................................................................................................................................................#

	my ($seqAroundSiteHsh_ref, $ggplotDirHsh_ref, $trimBaseComPlotRngHsh_ref) = @_;
	
	my $tmpItemHsh_ref = {};
	
	$tmpItemHsh_ref->{'mRNA_TSS'} = 'baseComp_mRNA_TSS';
	$tmpItemHsh_ref->{'NAT_TSS'} = 'baseComp_NAT_TSS';
	$tmpItemHsh_ref->{'mRNA_TAA'} = 'baseComp_mRNA_TAA';
	$tmpItemHsh_ref->{'mRNA_ATG'} = 'baseComp_mRNA_ATG';
	$tmpItemHsh_ref->{'NAT_TSS_originalDirtn'} = 'baseComp_NAT_TSS_originalDirtn';
	$tmpItemHsh_ref->{'mRNA_TAA_originalDirtn'} = 'baseComp_mRNA_TAA_originalDirtn';

	my $sub_ggplotDir = "baseCompstn";
	system "mkdir -pm 777 $ggplotDirHsh_ref->{$_}/$sub_ggplotDir/" foreach (keys %{$ggplotDirHsh_ref});
	
	foreach my $siteType (sort keys %{$seqAroundSiteHsh_ref}) {
		
		&reportStatus("Plotting sequence compostion around $siteType", 10, "\n");#->4274

		my $item = $tmpItemHsh_ref->{$siteType};
		
		my $seqAlignHsh_ref = {};
		foreach my $mRNAID (sort keys %{$seqAroundSiteHsh_ref->{$siteType}}) {
			my $fullSeq = $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'}.$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'};
			$seqAlignHsh_ref->{$mRNAID} = substr $fullSeq, $trimBaseComPlotRngHsh_ref->{'start'}, $trimBaseComPlotRngHsh_ref->{'length'};
		}


		my $seqLength = length ($seqAlignHsh_ref->{(keys %{$seqAlignHsh_ref})[rand keys %{$seqAlignHsh_ref}]}); #---from my $random_value = $hash{(keys %hash)[rand keys %hash]} in http://stackoverflow.com/questions/8547642/select-a-random-hash-key
		my $baseCompByBaseHsh_ref = &calculateBaseCompositionInAlignments($seqAlignHsh_ref);#->581
		
		{	
			my $plotDataHsh_ref = $baseCompByBaseHsh_ref;
			my $dataPath = "$ggplotDirHsh_ref->{'dat'}/$sub_ggplotDir/$item.dat";
			my $pdfPath = "$ggplotDirHsh_ref->{'pdf'}/$sub_ggplotDir/$item.pdf";
			my $RScriptPath = "$ggplotDirHsh_ref->{'R'}/$sub_ggplotDir/$item.R";
			my $logPath = "$ggplotDirHsh_ref->{'log'}/$sub_ggplotDir/$item.log";
			my $XAXis = 'relativePositon';
			my $YAxis = 'proportion';
			my $YVariable = 'base';
			my $extraArg = "+ ylim(0,1) + scale_x_continuous(breaks=seq(0, $seqLength, by=5))";
			my $height = 6;
			my $width = 14;
			&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $YVariable, $extraArg, $height, $width);#->2961
		}
	}
}
sub plotMotifOccurenceWithMASTMatch {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: currentTime|1071, getDREMEMotif|2034, getMASTLogPostionalData|2179, getMEMEMotif|2229, ggplotXYLinesMultipleSamples|2961, runMAST|4297
#	appearInSub: >none
#	primaryAppearInSection: 10_analyzeTSSMotif|365
#	secondaryAppearInSection: >none
#	input: $degenNtHsh_ref, $dremeXMLPathHsh_ref, $ggplotDirHsh_ref, $hardCodedIGVPathHsh_ref, $maxHitPVal, $memeXMLPathHsh_ref, $motifEValueLimitHsh_ref, $resultMastDir, $seqAroundSiteHsh_ref
#	output: 
#	toCall: &plotMotifOccurenceWithMASTMatch($dremeXMLPathHsh_ref, $memeXMLPathHsh_ref, $ggplotDirHsh_ref, $seqAroundSiteHsh_ref, $motifEValueLimitHsh_ref, $resultMastDir, $degenNtHsh_ref, $hardCodedIGVPathHsh_ref, $maxHitPVal);
#	calledInLine: 371
#....................................................................................................................................................#
	
	my ($dremeXMLPathHsh_ref, $memeXMLPathHsh_ref, $ggplotDirHsh_ref, $seqAroundSiteHsh_ref, $motifEValueLimitHsh_ref, $resultMastDir, $degenNtHsh_ref, $hardCodedIGVPathHsh_ref, $maxHitPVal) = @_;
	
	my $nonAnnotatedPolymerFreqPath = $hardCodedIGVPathHsh_ref->{'nonAnnotatedPolymerFreqPath'};
	my $fullSeqPathHsh_ref = {};
	my $mastFastaDir = "$resultMastDir/inputFasta/";
	system ("mkdir -pm 777 $mastFastaDir");
	
	my %totalSeqNumHsh;
	my %seqLengthHsh;
	
	foreach my $seqSiteType (sort keys %{$seqAroundSiteHsh_ref}) {
		$fullSeqPathHsh_ref->{$seqSiteType} = "$mastFastaDir/$seqSiteType.fasta";
		open FASTA, ">", $fullSeqPathHsh_ref->{$seqSiteType};
		$totalSeqNumHsh{$seqSiteType} = keys %{$seqAroundSiteHsh_ref->{$seqSiteType}};
		foreach my $mRNAID (sort keys %{$seqAroundSiteHsh_ref->{$seqSiteType}}) {
			print FASTA ">$mRNAID\n";
			print FASTA $seqAroundSiteHsh_ref->{$seqSiteType}{$mRNAID}{'upStrm'}.$seqAroundSiteHsh_ref->{$seqSiteType}{$mRNAID}{'dnStrm'}."\n";
			$seqLengthHsh{$seqSiteType} = length ($seqAroundSiteHsh_ref->{$seqSiteType}{$mRNAID}{'upStrm'}.$seqAroundSiteHsh_ref->{$seqSiteType}{$mRNAID}{'dnStrm'});
		}
		close FASTA;
	}

	foreach my $sub_ggplotDir (qw/meme dreme/) {
		system "mkdir -pm 777 $ggplotDirHsh_ref->{$_}/$sub_ggplotDir/" foreach (keys %{$ggplotDirHsh_ref});
	}

	#-----plot the Dreme and meme motif
	foreach my $siteType (sort keys %{$dremeXMLPathHsh_ref}) {

		print "[".&currentTime()."] Running mast for dreme motif occurrence for $siteType                          \n";#->1071
		
		foreach my $region (sort keys %{$dremeXMLPathHsh_ref->{$siteType}}) {
			{#-----plot the meme motif distribution
				my $memeXMLPath = $memeXMLPathHsh_ref->{$siteType}{$region};
				my ($motifInfoHsh_ref) = &getMEMEMotif($memeXMLPath, $degenNtHsh_ref);#->2229
				
				foreach my $seqSiteType (sort keys %{$fullSeqPathHsh_ref}) {
					my @mkDirAry;
					my $mastParentDir = "$resultMastDir/meme/$siteType/$region/$seqSiteType/"; push @mkDirAry, $mastParentDir;
					my $mastOutDir = "$mastParentDir/out"; push @mkDirAry, $mastOutDir;
					foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}
					my $XMLPath = $memeXMLPath;
					my $fastaPath = $fullSeqPathHsh_ref->{$seqSiteType};
					my ($mastHitLog) = &runMAST($XMLPath, $mastOutDir, $fastaPath, $nonAnnotatedPolymerFreqPath);#->4297

					my $maxMotifEVal = $motifEValueLimitHsh_ref->{'meme'};
					my $maxPos = $seqLengthHsh{$seqSiteType};
					my $totalSeqNum = $totalSeqNumHsh{$seqSiteType};
					my ($motifPctHsh_ref) = &getMASTLogPostionalData($mastHitLog, $motifInfoHsh_ref, $maxHitPVal, $maxMotifEVal, $maxPos, $totalSeqNum);#->2179
					
					{
						my $nameTag = "$siteType.$region.$seqSiteType.meme.mast.p$maxHitPVal";
						my $plotDataHsh_ref = $motifPctHsh_ref;
						my $item = $nameTag;
						my $sub_ggplotDir = 'meme';
						my $dataPath = "$ggplotDirHsh_ref->{'dat'}/$sub_ggplotDir/$item.dat";
						my $pdfPath = "$ggplotDirHsh_ref->{'pdf'}/$sub_ggplotDir/$item.pdf";
						my $RScriptPath = "$ggplotDirHsh_ref->{'R'}/$sub_ggplotDir/$item.R";
						my $logPath = "$ggplotDirHsh_ref->{'log'}/$sub_ggplotDir/$item.log";
						my $XAXis = 'relativePositon';
						my $YAxis = 'proportion';
						my $YVariable = 'motif';
						my $seqLength = $seqLengthHsh{$seqSiteType};
						my $extraArg = " + scale_x_continuous(breaks=seq(0, $seqLength, by=10))";
						my $height = 6;
						my $width = 14;
						&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $YVariable, $extraArg, $height, $width);#->2961
					}
				}
			}

			{
				#-----plot the dreme motif distribution
				foreach my $shuffleOrNegative (sort keys %{$dremeXMLPathHsh_ref->{$siteType}{$region}}) {

					my $dremeXMLPath = $dremeXMLPathHsh_ref->{$siteType}{$region}{$shuffleOrNegative};
					my $motifInfoHsh_ref = &getDREMEMotif($dremeXMLPath, $degenNtHsh_ref);#->2034

					foreach my $seqSiteType (sort keys %{$fullSeqPathHsh_ref}) {
						my @mkDirAry;
						my $mastParentDir = "$resultMastDir/dreme/$siteType/$shuffleOrNegative/$region/$seqSiteType"; push @mkDirAry, $mastParentDir;
						my $mastOutDir = "$mastParentDir/out"; push @mkDirAry, $mastOutDir;
						foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}
					
						my $mastTag = "meme.$siteType.$region.$seqSiteType";
						my $XMLPath = $dremeXMLPath;
						my $fastaPath = $fullSeqPathHsh_ref->{$seqSiteType};
						my ($mastHitLog) = &runMAST($XMLPath, $mastOutDir, $fastaPath, $nonAnnotatedPolymerFreqPath);#->4297

						my $maxMotifEVal = $motifEValueLimitHsh_ref->{'dreme'};
						my $maxPos = $seqLengthHsh{$seqSiteType};
						my $totalSeqNum = $totalSeqNumHsh{$seqSiteType};
						my ($motifPctHsh_ref) = &getMASTLogPostionalData($mastHitLog, $motifInfoHsh_ref, $maxHitPVal, $maxMotifEVal, $maxPos, $totalSeqNum);#->2179

						{
							my $nameTag = "$siteType.$region.$seqSiteType.dreme.$shuffleOrNegative.mast.p$maxHitPVal";
							my $plotDataHsh_ref = $motifPctHsh_ref;
							my $item = $nameTag;
							my $sub_ggplotDir = 'dreme';
							my $dataPath = "$ggplotDirHsh_ref->{'dat'}/$sub_ggplotDir/$item.dat";
							my $pdfPath = "$ggplotDirHsh_ref->{'pdf'}/$sub_ggplotDir/$item.pdf";
							my $RScriptPath = "$ggplotDirHsh_ref->{'R'}/$sub_ggplotDir/$item.R";
							my $logPath = "$ggplotDirHsh_ref->{'log'}/$sub_ggplotDir/$item.log";
							my $XAXis = 'relativePositon';
							my $YAxis = 'proportion';
							my $YVariable = 'motif';
							my $seqLength = $seqLengthHsh{$seqSiteType};
							my $extraArg = " + scale_x_continuous(breaks=seq(0, $seqLength, by=10))";
							my $height = 6;
							my $width = 14;
							&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $YVariable, $extraArg, $height, $width);#->2961
						}
					}
				}
			}
		}
	}

	return ();
}
sub plotTSSRelativeToATGAndTAAHistogram {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: ggplotHistogram|2780
#	appearInSub: >none
#	primaryAppearInSection: 8_investigateGeneTSSRelativePos|331
#	secondaryAppearInSection: >none
#	input: $TSSmRNAEndFreqHsh_ref, $ggplotDirHsh_ref, $strndEndValidTSSInfoHsh_ref, $validTSSLimitHsh_ref
#	output: none
#	toCall: &plotTSSRelativeToATGAndTAAHistogram($TSSmRNAEndFreqHsh_ref, $ggplotDirHsh_ref, $strndEndValidTSSInfoHsh_ref, $validTSSLimitHsh_ref);
#	calledInLine: 343
#....................................................................................................................................................#

	my ($TSSmRNAEndFreqHsh_ref, $ggplotDirHsh_ref, $strndEndValidTSSInfoHsh_ref, $validTSSLimitHsh_ref) = @_;

	my $tmpPlotItemHsh_ref = {};
	$tmpPlotItemHsh_ref->{'ATG'}{'sense'} = 'TSSATGSenseHist';
	$tmpPlotItemHsh_ref->{'TAA'}{'antisense'} = 'TSSTAAAntisenseHist';

	my $sub_ggplotDir = "validTSSDstrbtn";
	system "mkdir -pm 777 $ggplotDirHsh_ref->{$_}/$sub_ggplotDir/" foreach (keys %{$ggplotDirHsh_ref});

	foreach my $ATGOrTAA (keys %{$strndEndValidTSSInfoHsh_ref}) {
		foreach my $senseOrAntisense (keys %{$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}}) {
			
			my $geneWithTSSWithinRng = keys %{$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}};

			my $item = $tmpPlotItemHsh_ref->{$ATGOrTAA}{$senseOrAntisense};
			my $plotAry_ref = $TSSmRNAEndFreqHsh_ref->{$ATGOrTAA}{$senseOrAntisense};
			my $dataPath = "$ggplotDirHsh_ref->{'dat'}/$sub_ggplotDir/$item.dat";
			my $pdfPath = "$ggplotDirHsh_ref->{'pdf'}/$sub_ggplotDir/$item.pdf";
			my $RScriptPath = "$ggplotDirHsh_ref->{'R'}/$sub_ggplotDir/$item.R";
			my $logPath = "$ggplotDirHsh_ref->{'log'}/$sub_ggplotDir/$item.log";
			my $leftXAxisPercentileLimit = 'min';
			my $rightXAxisPercentileLimit = 'max';
			my $binWidth = 1;
			my $xAxis = "TSS.around.$ATGOrTAA.$senseOrAntisense.$geneWithTSSWithinRng.genes.within.rng";
			my $dataPtMax = 99999999;
			my $log2OrLinear = 'linear';
			my $extraArg = '';
			my $lowerValLimit = $validTSSLimitHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{'lowerValLimit'};
			my $upperValLimit = $validTSSLimitHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{'upperValLimit'};
			my $lowerPctLimit = $validTSSLimitHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{'lowerPctLimit'};
			my $upperPctLimit = $validTSSLimitHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{'upperPctLimit'};
			$extraArg .= " + geom_vline(xintercept=c($lowerValLimit), linetype=\"dotted\") + annotate(\"text\", x=$lowerValLimit, y=0, label=\"$lowerPctLimit\%\=$lowerValLimit\", vjust=-0.2, hjust=-0.1, angle=90)";
			$extraArg .= " + geom_vline(xintercept=c($upperValLimit), linetype=\"dotted\") + annotate(\"text\", x=$upperValLimit, y=0, label=\"$upperPctLimit\%\=$upperValLimit\", vjust=-0.2, hjust=-0.1, angle=90)";
			my $height = 8;
			my $width = 16;
			&ggplotHistogram($plotAry_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftXAxisPercentileLimit, $rightXAxisPercentileLimit, $xAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear, $height, $width);#->2780
		}
	}
}
sub plotWeblogoAroundTSS {
#....................................................................................................................................................#
#	subroutineCategory: specific, thirdPartyApp
#	dependOnSub: createWeblogo|1045
#	appearInSub: >none
#	primaryAppearInSection: 9_analyzeNucleotideComposition|349
#	secondaryAppearInSection: >none
#	input: $seqAroundSiteHsh_ref, $weblogoDirHsh_ref
#	output: none
#	toCall: &plotWeblogoAroundTSS($seqAroundSiteHsh_ref, $weblogoDirHsh_ref);
#	calledInLine: 359
#....................................................................................................................................................#
	
	my ($seqAroundSiteHsh_ref, $weblogoDirHsh_ref) = @_;
	
	my $upStrmRng = 5;
	my $dnStrmRng = 5;
	
	foreach my $siteType (keys %{$seqAroundSiteHsh_ref}) {
		my $seqAlignHsh_ref = {};
		foreach my $mRNAID (keys %{$seqAroundSiteHsh_ref->{$siteType}}) {
			my $upStrmSeq = substr $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'}, -1*$upStrmRng, $upStrmRng;
			my $dnStrmRng = substr $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'}, 0, $dnStrmRng;
			$seqAlignHsh_ref->{$mRNAID} = $upStrmSeq.$dnStrmRng;
		}
		
		{
			my $numSeq = keys %{$seqAroundSiteHsh_ref->{$siteType}};
			my $nameTag = "around.$siteType";
			my $pdfPath = "$weblogoDirHsh_ref->{pdf}/$nameTag.pdf";
			my $fastaPath = "$weblogoDirHsh_ref->{fasta}/$nameTag.fasta";
			my $title = "$siteType\_-$upStrmRng..+$dnStrmRng\[n=$numSeq\]";
			my $seqType = 'dna';
			&createWeblogo($seqAlignHsh_ref, $pdfPath, $fastaPath, $seqType, $title);#->1045
		}
	}
	
}
sub poolLibRd5EndCov {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: createEmptyGenomeCovPerlStorable|1010, getIndivCntgCovPlsPath|2146, storeAndScaleRd5EndFromBothLibraries|4452, warningToProceed|4583
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|273
#	secondaryAppearInSection: >none
#	input: $fastaHsh_ref, $forcePoolLibRd5End, $libTypeInfoHsh_ref, $maxThread, $originalRd5EndStorableDir, $scaledRd5EndStorableDir
#	output: $rd5EndPlsInfoHsh_ref
#	toCall: my ($rd5EndPlsInfoHsh_ref) = &poolLibRd5EndCov($fastaHsh_ref, $originalRd5EndStorableDir, $scaledRd5EndStorableDir, $libTypeInfoHsh_ref, $forcePoolLibRd5End, $maxThread);
#	calledInLine: 296
#....................................................................................................................................................#

	my ($fastaHsh_ref, $originalRd5EndStorableDir, $scaledRd5EndStorableDir, $libTypeInfoHsh_ref, $forcePoolLibRd5End, $maxThread) = @_;
	
	my $rd5EndPlsInfoHsh_ref = {};
	
	$rd5EndPlsInfoHsh_ref->{'original'}{'rd5EndIdxHshPath'} = "$originalRd5EndStorableDir/index.hsh.pls";
	$rd5EndPlsInfoHsh_ref->{'scaled'}{'rd5EndIdxHshPath'} = "$scaledRd5EndStorableDir/index.hsh.pls";
	$rd5EndPlsInfoHsh_ref->{'original'}{'storableDir'} = "$originalRd5EndStorableDir";
	$rd5EndPlsInfoHsh_ref->{'scaled'}{'storableDir'} = "$scaledRd5EndStorableDir";

	if ((-s $rd5EndPlsInfoHsh_ref->{'original'}{'rd5EndIdxHshPath'} or -s $rd5EndPlsInfoHsh_ref->{'original'}{'rd5EndIdxHshPath'}.".gz") 
		and (-s $rd5EndPlsInfoHsh_ref->{'scaled'}{'rd5EndIdxHshPath'} or -s $rd5EndPlsInfoHsh_ref->{'scaled'}{'rd5EndIdxHshPath'}.".gz") 
		and $forcePoolLibRd5End eq 'no') {
		my $message = "originalRd5EndIdxHshPath and scaledRd5EndIdxHshPath found, will be retrieved in a few seconds......\r";
		&warningToProceed($message);#->4583
		foreach my $originalOrScaled (keys %{$rd5EndPlsInfoHsh_ref}) {
			$rd5EndPlsInfoHsh_ref->{$originalOrScaled}{'covPlsPathHsh_ref'} = &getIndivCntgCovPlsPath($rd5EndPlsInfoHsh_ref->{$originalOrScaled}{'rd5EndIdxHshPath'});#->2146
		}
		
	} else {

		foreach my $originalOrScaled (keys %{$rd5EndPlsInfoHsh_ref}) {
			#----------Create empty stroable
			&createEmptyGenomeCovPerlStorable($rd5EndPlsInfoHsh_ref->{$originalOrScaled}{'storableDir'}, $fastaHsh_ref);#->1010
			$rd5EndPlsInfoHsh_ref->{$originalOrScaled}{'covPlsPathHsh_ref'} = &getIndivCntgCovPlsPath($rd5EndPlsInfoHsh_ref->{$originalOrScaled}{'rd5EndIdxHshPath'});#->2146
		}

		&storeAndScaleRd5EndFromBothLibraries($rd5EndPlsInfoHsh_ref, $libTypeInfoHsh_ref, $fastaHsh_ref, $maxThread);#->4452
	}
	
	return ($rd5EndPlsInfoHsh_ref);
}
sub predictBonaFideTSSinGoldenmRNASet {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: defineGoldenmRNASet|1123, downSampleGoldenATGRd5EndDataHshForTesting|1227, findDistanceBetweenATGAndTrnsfrgEnd|1386, getRd5EndSurroundingmRNAATG|2584, identifiyPeakSiteWithinRegion|3031, outputTSSResultWiggle|3267
#	appearInSub: trainingUsingExonVsTssData|4548
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_trainningSTDTEXRatio|308
#	input: $forceRunGoldenmRNASetTSSi, $ggplotDirHsh_ref, $goldenmRNASetSizeLimit, $libTypeInfoHsh_ref, $mRNAByCntgHsh_ref, $mRNAInfoHsh_ref, $maxPctDistBtwATGTrnsfrgEnd, $maxThread, $rd5EndPlsInfoHsh_ref, $trnsfrgInfoHsh_ref, $wigglePathHsh_ref
#	output: $goldenmRNATSSResultHsh_ref
#	toCall: my ($goldenmRNATSSResultHsh_ref) = &predictBonaFideTSSinGoldenmRNASet($forceRunGoldenmRNASetTSSi, $trnsfrgInfoHsh_ref, $mRNAInfoHsh_ref, $ggplotDirHsh_ref, $maxPctDistBtwATGTrnsfrgEnd, $libTypeInfoHsh_ref, $mRNAByCntgHsh_ref, $rd5EndPlsInfoHsh_ref, $wigglePathHsh_ref, $goldenmRNASetSizeLimit, $maxThread);
#	calledInLine: 4572
#....................................................................................................................................................#
	
	my ($forceRunGoldenmRNASetTSSi, $trnsfrgInfoHsh_ref, $mRNAInfoHsh_ref, $ggplotDirHsh_ref, $maxPctDistBtwATGTrnsfrgEnd, $libTypeInfoHsh_ref, $mRNAByCntgHsh_ref, $rd5EndPlsInfoHsh_ref, $wigglePathHsh_ref, $goldenmRNASetSizeLimit, $maxThread) = @_;

	my $upStrmBufferSize = 20; #----extra size added to the distBtwATGTrnsfrgEndCutoff
	
	#-----find the distance between 
	my ($distBtwATGTrnsfrgEndHsh_ref) = &findDistanceBetweenATGAndTrnsfrgEnd($trnsfrgInfoHsh_ref, $mRNAInfoHsh_ref, $ggplotDirHsh_ref, $maxThread);#->1386
	
	#-----define a golden set of 
	my ($goldenmRNAListHsh_ref, $distBtwATGTrnsfrgEndCutoff) = &defineGoldenmRNASet($distBtwATGTrnsfrgEndHsh_ref, $maxPctDistBtwATGTrnsfrgEnd, $goldenmRNASetSizeLimit);#->1123

	#-----get the rd5EndCount around the mRNA ATG
	my $TEXOrSTD = 'TEX';
	my $originalOrScaled = 'original';
	my $ATGUpStrmSrchRegSize = $distBtwATGTrnsfrgEndCutoff+$upStrmBufferSize;
	my $ATGDnStrmSrchRegSize = 0;
	my ($goldenATGRd5EndDataHsh_ref) = &getRd5EndSurroundingmRNAATG($mRNAInfoHsh_ref, $mRNAByCntgHsh_ref, $rd5EndPlsInfoHsh_ref, $ATGUpStrmSrchRegSize, $ATGDnStrmSrchRegSize, $libTypeInfoHsh_ref, $goldenmRNAListHsh_ref, $TEXOrSTD, $originalOrScaled);#->2584

	#----to down sample
	#$goldenATGRd5EndDataHsh_ref = &downSampleGoldenATGRd5EndDataHshForTesting($goldenATGRd5EndDataHsh_ref);#->1227

	my ($goldenmRNATSSResultHsh_ref) = &identifiyPeakSiteWithinRegion($goldenATGRd5EndDataHsh_ref, $mRNAInfoHsh_ref);#->3031
	
	my $TSSItem = 'TSS_goldenmRNA_TEX_original';
	my $goldenTSSWigPathHsh_ref = {};
	$goldenTSSWigPathHsh_ref->{'+'} = $wigglePathHsh_ref->{$TSSItem}{'plus'};
	$goldenTSSWigPathHsh_ref->{'-'} = $wigglePathHsh_ref->{$TSSItem}{'minus'};
	&outputTSSResultWiggle($goldenTSSWigPathHsh_ref, $goldenmRNATSSResultHsh_ref);#->3267
	
	return ($goldenmRNATSSResultHsh_ref);
}
sub printBothFHAndStdout {
#....................................................................................................................................................#
#	subroutineCategory: general, reporting
#	dependOnSub: currentTime|1071
#	appearInSub: calculateNoiseIndex|625
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 8_investigateGeneTSSRelativePos|331
#	input: $FH, $message, $numTrailingSpace
#	output: 
#	toCall: &printBothFHAndStdout($message, $numTrailingSpace, $FH);
#	calledInLine: 669, 670, 671
#....................................................................................................................................................#
	
	my ($message, $numTrailingSpace, $FH) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);

	print "[".&currentTime()."] ".$message.$trailingSpaces."\n";#->1071
	print {$FH} "[".&currentTime()."] ".$message."\n";#->1071

	return ();
}
sub printCMDLogOrFinishMessage {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|1071
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|114, 12_finishingTasks|391
#	secondaryAppearInSection: >none
#	input: $CMDLogOrFinishMessage
#	output: none
#	toCall: &printCMDLogOrFinishMessage($CMDLogOrFinishMessage);
#	calledInLine: 120, 400
#....................................................................................................................................................#

	my ($CMDLogOrFinishMessage) = @_;
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		print CMDLOG "[".&currentTime()."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";#->1071
		close CMDLOG;
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->1071
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->1071
		print "=========================================================================\n\n";
	}
}
sub printGFF_oneRNAPerGene_chooseStrnd_filterAry {
#....................................................................................................................................................#
#	subroutineCategory: gff
#	dependOnSub: >none
#	appearInSub: getAndPrintTrnsfrg|1672
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_processInputData|273
#	input: $filterAry_ref, $geneInfoHsh_ref, $gffGeneLineOnly, $outGFFPath, $strndInWord
#	output: none
#	toCall: &printGFF_oneRNAPerGene_chooseStrnd_filterAry($geneInfoHsh_ref, $outGFFPath, $gffGeneLineOnly, $strndInWord, $filterAry_ref);
#	calledInLine: 1692
#....................................................................................................................................................#

	my ($geneInfoHsh_ref, $outGFFPath, $gffGeneLineOnly, $strndInWord, $filterAry_ref) = @_;

	#------gffGeneLineOnly: print only the gene line, which is perfect for trnsfrgs
	#------strndInWord: plus, minus or both
		
	$outGFFPath =~ s/\/+/\//g;
	
	my $strndOut = '+';
	$strndOut = '-' if $strndInWord eq 'minus';
	
	open (GFFOUT, ">", "$outGFFPath");
	print GFFOUT "##gff-version\t3\n";

	foreach my $geneID (sort {$a cmp $b} keys %{$geneInfoHsh_ref}) {
		
		my $toPrint = 'yes';
		foreach my $filter (@{$filterAry_ref}) {
			$toPrint = 'no' if ($geneInfoHsh_ref->{$geneID}{$filter}) eq 'no';
		}
		
		if ($toPrint eq 'yes') {
		
			my $strnd = $geneInfoHsh_ref->{$geneID}{"strnd"};

			next if $strndOut ne $strnd and $strndInWord ne 'both';
		
			my $cntg = $geneInfoHsh_ref->{$geneID}{"cntg"};
			my $ctgry = $geneInfoHsh_ref->{$geneID}{"ctgry"};
			my $description = $geneInfoHsh_ref->{$geneID}{"description"};
		
			my ($geneStart, $geneEnd) = sort {$a <=> $b} @{$geneInfoHsh_ref->{$geneID}{'geneRng'}};
			print GFFOUT join "", (join "\t", ($cntg, 'BCP', 'gene', $geneStart, $geneEnd, ".", $strnd, ".", "ID=$geneID;Name=$geneID;description=$description")), "\n";

			if ($gffGeneLineOnly eq 'no') {#-----will print the RNA and exon lines also, aim to avoid display annoyance on IGV
				my $RNAID = $geneInfoHsh_ref->{$geneID}{"RNAID"};
				my ($RNAStart, $RNAEnd) = sort {$a <=> $b} @{$geneInfoHsh_ref->{$geneID}{'RNARng'}};
				print GFFOUT join "", (join "\t", ($cntg, 'BCP', $ctgry, $RNAStart, $RNAEnd, ".", $strnd, ".", "ID=$RNAID;Name=$geneID;Parent=$geneID;description=$description;")), "\n";
			
				foreach my $rngType ('exonRng', 'CDSRng') {
					my $gffRng = 'exon';
					$gffRng = 'CDS' if $rngType eq 'CDSRng';
				
					if ($geneInfoHsh_ref->{$geneID}{$rngType}) {
						my @rngAry = sort {$a <=> $b} @{$geneInfoHsh_ref->{$geneID}{$rngType}};
						my $num = 1;
						for (my $i=0; $i < $#rngAry; $i += 2) {
							$num++;
							my ($start, $end) = ($rngAry[$i], $rngAry[$i+1]);
							my $ID = "$rngType\_$num\_$RNAID";
							print GFFOUT join "", (join "\t", ($cntg, "BCP", $gffRng, $start, $end, ".", $strnd, ".", "ID=$ID;Name=$ID;Parent=$RNAID;description=.;")), "\n";
						}
					}
				} #---end of foreach my $rngType ('exonRng', 'CDSRng') 
			} #---end of if ($gffGeneLineOnly eq 'no') 
		} #---end of if ($toPrint eq 'yes')
	}
	close GFFOUT;
}
sub printGeneBasedTSSLog {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 8_investigateGeneTSSRelativePos|331
#	secondaryAppearInSection: >none
#	input: $geneBasedTSSInfoHsh_ref, $geneBasedTSSLogPath, $strndEndValidTSSInfoHsh_ref
#	output: none
#	toCall: &printGeneBasedTSSLog($geneBasedTSSLogPath, $geneBasedTSSInfoHsh_ref, $strndEndValidTSSInfoHsh_ref);
#	calledInLine: 344
#....................................................................................................................................................#

	my ($geneBasedTSSLogPath, $geneBasedTSSInfoHsh_ref, $strndEndValidTSSInfoHsh_ref) = @_;
	
	open (GENETSSLOG, ">", $geneBasedTSSLogPath);
	foreach my $mRNAID (sort keys %{$geneBasedTSSInfoHsh_ref}) {
		my @outputAry;
		push @outputAry, 'mRNAID';
		foreach my $ATGOrTAA (sort keys %{$strndEndValidTSSInfoHsh_ref}) {
			foreach my $senseOrAntisense (sort keys %{$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}}) {
				foreach my $rltvPosOrCov (sort keys %{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}}) {
					push @outputAry, $ATGOrTAA."_".$senseOrAntisense."_".$rltvPosOrCov;
				}
				push @outputAry, $ATGOrTAA."_".$senseOrAntisense."_withinRng";
			}
		}
		print GENETSSLOG join "", ((join "\t", @outputAry), "\n");
		last;
	}

	foreach my $mRNAID (sort keys %{$geneBasedTSSInfoHsh_ref}) {
		my @outputAry;
		push @outputAry, $mRNAID;
		foreach my $ATGOrTAA (sort keys %{$strndEndValidTSSInfoHsh_ref}) {
			foreach my $senseOrAntisense (sort keys %{$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}}) {
				foreach my $rltvPosOrCov (sort keys %{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}}) {
					my $str = join ",", @{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{$rltvPosOrCov}};
					push @outputAry, $str;
				}
				my $withinRng = 'no';
				$withinRng = 'yes' if defined $strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{$mRNAID};
				push @outputAry, $withinRng;
			}
		}
		print GENETSSLOG join "", ((join "\t", @outputAry), "\n");
	}
	close GENETSSLOG;

}
sub printIGVXML {
#....................................................................................................................................................#
#	subroutineCategory: specific, XML
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 11_outputWiggleXML|376
#	secondaryAppearInSection: >none
#	input: $IGVGenomeID, $IGVTrackInfoHsh_ref, $XMLPath, $copyIGVTrack, $fontSize, $hardCodedIGVPathHsh_ref, $resultIGVTrackCopyDir
#	output: none
#	toCall: &printIGVXML($fontSize, $XMLPath, $hardCodedIGVPathHsh_ref, $IGVGenomeID, $IGVTrackInfoHsh_ref, $copyIGVTrack, $resultIGVTrackCopyDir);
#	calledInLine: 386
#....................................................................................................................................................#

	my ($fontSize, $XMLPath, $hardCodedIGVPathHsh_ref, $IGVGenomeID, $IGVTrackInfoHsh_ref, $copyIGVTrack, $resultIGVTrackCopyDir) = @_;
	
	open (XML, ">", $XMLPath);
	printf XML "%0s", "<\?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"\?>\n";
	printf XML "%0s", "<Session genome=\"$IGVGenomeID\" version=\"5\">\n";
	printf XML "%4s", "<Resources>\n";
	
	if ($copyIGVTrack eq 'yes') {
		foreach my $track (keys %{$hardCodedIGVPathHsh_ref}) {
			my $originalPath = $hardCodedIGVPathHsh_ref->{$track};
			die "$originalPath for IGV not found\n" if not -s $originalPath;
			my ($name, $dir, $extension) = fileparse($originalPath, qr/\.[^.]*/);
			my $newPath = "$resultIGVTrackCopyDir/$name$extension";
			system ("cp -f $originalPath $resultIGVTrackCopyDir");
			system ("cp -f $originalPath.bai $resultIGVTrackCopyDir") if $extension eq '.bam';
			$hardCodedIGVPathHsh_ref->{$track} = $newPath;
		}
	
		foreach my $originalPath (keys %{$IGVTrackInfoHsh_ref}) {
			die "$originalPath for IGV not found\n" if not -s $originalPath;
			my ($name, $dir, $extension) = fileparse($originalPath, qr/\.[^.]*/);
			system ("cp -f $originalPath $resultIGVTrackCopyDir");
			my $newPath = "$resultIGVTrackCopyDir/$name$extension";
			%{$IGVTrackInfoHsh_ref->{$newPath}} = %{$IGVTrackInfoHsh_ref->{$originalPath}};
			delete $IGVTrackInfoHsh_ref->{$originalPath};
		}
	}

	foreach my $track (keys %{$hardCodedIGVPathHsh_ref}) {
		my $path = $hardCodedIGVPathHsh_ref->{$track};
		die "$path for IGV not found\n" if not -s $path;
		printf XML "%8s", "<Resource path=\"$path\"\/>\n";
	}

	foreach my $path (keys %{$IGVTrackInfoHsh_ref}) {
		die "$path for IGV not found\n" if not -s $path;
		printf XML "%8s", "<Resource path=\"$path\"\/>\n";
	}

	printf XML "%4s", "</Resources>\n";
	printf XML "%4s", "<Panel height=\"900\" name=\"DataPanel\" width=\"1200\">\n";

	printf XML "%8s", "<Track colorScale=\"ContinuousColorScale;20.0;10.0;25.0;35.0;153,153,255;255,255,255;255,153,153\" fontSize=\"$fontSize\" id=\"$hardCodedIGVPathHsh_ref->{'GCPath'}\"  name=\"GC\" normalize=\"false\" renderer=\"HEATMAP\" height=\"20\" sortable=\"true\" visible=\"true\" windowFunction=\"none\">\n";
	printf XML "%12s", "<DataRange type=\"LINEAR\"/>\n";
	printf XML "%8s", "</Track>\n";

	printf XML "%8s", "<Track colorScale=\"ContinuousColorScale;2.0;1.0;3.0;15.0;255,255,255;153,153,255;255,153,153\" fontSize=\"$fontSize\" id=\"$hardCodedIGVPathHsh_ref->{'reptvPath'}\"  name=\"Repetitiveness\" normalize=\"false\" renderer=\"HEATMAP\" height=\"20\" sortable=\"true\" visible=\"true\" windowFunction=\"none\">\n";
	printf XML "%12s", "<DataRange type=\"LINEAR\"/>\n";
	printf XML "%8s", "</Track>\n";

	printf XML "%8s", "<Track autoScale=\"false\" clazz=\"org.broad.igv.track.FeatureTrack\" colorScale=\"ContinuousColorScale;0.0;17.0;255,255,255;0,0,178\" displayMode=\"SQUISHED\"  fontSize=\"$fontSize\" height=\"20\" id=\"EHI_v13_genes\" name=\"bonaFide mRNA\" renderer=\"BASIC_FEATURE\" sortable=\"false\" visible=\"true\" windowFunction=\"count\">\n";
	printf XML "%12s", "<DataRange type=\"LINEAR\"/>\n";
	printf XML "%8s", "</Track>\n";

	printf XML "%8s", "<Track autoScale=\"false\" clazz=\"org.broad.igv.track.FeatureTrack\" colorScale=\"ContinuousColorScale;0.0;40.0;255,255,255;0,0,178\" displayMode=\"COLLAPSED\"  fontSize=\"$fontSize\" height=\"20\" id=\"$hardCodedIGVPathHsh_ref->{'gffPath'}\"  name=\"All Genomic Features\" renderer=\"GENE_TRACK\" sortable=\"false\" visible=\"true\" windowFunction=\"count\">\n";
	printf XML "%12s", "<DataRange type=\"LINEAR\"/>\n";
	printf XML "%8s", "</Track>\n";

	printf XML "%8s", "<Track autoScale=\"true\" color=\"255,153,153\" displayMode=\"COLLAPSED\" fontSize=\"$fontSize\" height=\"100\" id=\"$hardCodedIGVPathHsh_ref->{'plusCovTDFPath'}\"  name=\"Corrected Plus Cov\" normalize=\"false\" renderer=\"BAR_CHART\" sortable=\"true\" visible=\"true\" windowFunction=\"mean\">\n";
	printf XML "%12s", "<DataRange type=\"LOG\"/>\n";
	printf XML "%8s", "</Track>\n";

	printf XML "%8s", "<Track autoScale=\"true\" color=\"153,153,255\" displayMode=\"COLLAPSED\" fontSize=\"$fontSize\" height=\"100\" id=\"$hardCodedIGVPathHsh_ref->{'minusCovTDFPath'}\"  name=\"Corrected Minus Cov\" normalize=\"false\" renderer=\"BAR_CHART\" sortable=\"true\" visible=\"true\" windowFunction=\"mean\">\n";
	printf XML "%12s", "<DataRange type=\"LOG\"/>\n";
	printf XML "%8s", "</Track>\n";

	foreach my $path (sort keys %{$IGVTrackInfoHsh_ref}) {
		my $name = $IGVTrackInfoHsh_ref->{$path}{'name'};
		my $displayMode = $IGVTrackInfoHsh_ref->{$path}{'displayMode'};
		my $DataRange_type = $IGVTrackInfoHsh_ref->{$path}{'DataRange_type'};
		my $autoScale = $IGVTrackInfoHsh_ref->{$path}{'autoScale'};
		my $color = $IGVTrackInfoHsh_ref->{$path}{'color'};
		my $windowFunction = $IGVTrackInfoHsh_ref->{$path}{'windowFunction'};
		my $renderer = $IGVTrackInfoHsh_ref->{$path}{'renderer'};
		my $height = $IGVTrackInfoHsh_ref->{$path}{'height'};

		printf XML "%8s", "<Track autoScale=\"$autoScale\" clazz=\"org.broad.igv.track.FeatureTrack\" displayMode=\"$displayMode\" color=\"$color\" fontSize=\"$fontSize\" height=\"$height\" id=\"$path\" name=\"$name\" renderer=\"$renderer\" sortable=\"false\" visible=\"true\" windowFunction=\"$windowFunction\">\n";
		printf XML "%12s", "<DataRange type=\"$DataRange_type\"/>\n";
		printf XML "%8s", "</Track>\n";
	}
	
	printf XML "%4s", "</Panel>\n";
	printf XML "%4s", "<Panel height=\"300\" width=\"1200\">\n";
	printf XML "%8s", "<Track altColor=\"0,0,178\" autoScale=\"false\" color=\"0,0,178\" displayMode=\"EXPANDED\" fontSize=\"$fontSize\" id=\"$hardCodedIGVPathHsh_ref->{'polyABamPath'}\" name=\"polyAReads\" showSpliceJunctions=\"false\" sortable=\"true\" visible=\"true\">\n";
	printf XML "%12s", "<RenderOptions colorByTag=\"\" colorOption=\"READ_STRAND\" flagUnmappedPairs=\"false\" groupByTag=\"\" maxInsertSize=\"1000\" minInsertSize=\"50\" shadeBasesOption=\"QUALITY\" shadeCenters=\"true\" showAllBases=\"false\" sortByTag=\"\"/>\n";
	printf XML "%8s", "</Track>\n";
	printf XML "%4s", "</Panel>\n";
	printf XML "%4s", "<PanelLayout dividerFractions=\"0.6637168141592921,0.9494310998735778\"/>\n";

	printf XML "%0s", "</Session>\n";

	close XML;

}
sub printMaxCovAndPosAroundGeneEdgeInfo {
#....................................................................................................................................................#
#	subroutineCategory: coverage, reporting
#	dependOnSub: >none
#	appearInSub: getPutativeTSSRegionBasedOnPeakValidSite|2364
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 8_investigateGeneTSSRelativePos|331
#	input: $maxPosCovByItemHsh_ref, $nameTag, $resultLogDir
#	output: 
#	toCall: &printMaxCovAndPosAroundGeneEdgeInfo($maxPosCovByItemHsh_ref, $resultLogDir, $nameTag);
#	calledInLine: 2394
#....................................................................................................................................................#
	my ($maxPosCovByItemHsh_ref, $resultLogDir, $nameTag) = @_;
	
	open MAXCOVLOG, ">", "$resultLogDir/$nameTag.maxPosCov.xls";

	foreach my $itemID (sort keys %{$maxPosCovByItemHsh_ref}) {
		my @outputAry = ("itemID");
		foreach my $dirtn (sort keys %{$maxPosCovByItemHsh_ref->{$itemID}}) {
			foreach my $headOrTail (sort keys %{$maxPosCovByItemHsh_ref->{$itemID}{$dirtn}}) {
				push @outputAry, ($dirtn.".".$headOrTail.".maxCov", $dirtn.".".$headOrTail.".maxPos");
			}
		}
		print MAXCOVLOG join "", (join "\t", (@outputAry)), "\n";
		last;
	}

	foreach my $itemID (sort keys %{$maxPosCovByItemHsh_ref}) {
		my @outputAry = ($itemID);
		foreach my $dirtn (sort keys %{$maxPosCovByItemHsh_ref->{$itemID}}) {
			foreach my $headOrTail (sort keys %{$maxPosCovByItemHsh_ref->{$itemID}{$dirtn}}) {
				my $maxCov = $maxPosCovByItemHsh_ref->{$itemID}{$dirtn}{$headOrTail}{'maxCov'};
				my $maxPos = $maxPosCovByItemHsh_ref->{$itemID}{$dirtn}{$headOrTail}{'maxPos'};
				push @outputAry, ($maxCov, $maxPos);
			}
		}
		print MAXCOVLOG join "", (join "\t", (@outputAry)), "\n";
	}
	close MAXCOVLOG;
	
	return ();
}
sub printWiggleSingleTrackFromCntgCovPlsPathHsh {
#....................................................................................................................................................#
#	subroutineCategory: wiggle
#	dependOnSub: >none
#	appearInSub: outputOriginalAndScaledRd5EndWiggle|3234
#	primaryAppearInSection: 7_filterTEXTrackAfterTrainning|319
#	secondaryAppearInSection: 11_outputWiggleXML|376
#	input: $aryIndex, $cntgCovPlsPathHsh_ref, $gzip, $wigPath
#	output: none
#	toCall: &printWiggleSingleTrackFromCntgCovPlsPathHsh($cntgCovPlsPathHsh_ref, $aryIndex, $wigPath, $gzip);
#	calledInLine: 325, 326, 3256
#....................................................................................................................................................#
	
	my ($cntgCovPlsPathHsh_ref, $aryIndex, $wigPath, $gzip) = @_;
	
	open (WIGGLE, ">", $wigPath);
	
	foreach my $cntg (sort keys %{$cntgCovPlsPathHsh_ref}) {

		print WIGGLE "variableStep chrom=$cntg span=1\n";

		my $cntgCovPlsPath = "$cntgCovPlsPathHsh_ref->{$cntg}";
 		system ("gzip -df $cntgCovPlsPath.gz") if (-s "$cntgCovPlsPath.gz" and $gzip eq 'yes');
		my $cntgCovAry_ref = retrieve($cntgCovPlsPath);
		system ("gzip -f $cntgCovPlsPath") if (-s $cntgCovPlsPath and $gzip eq 'yes');
		for my $i (0..$#{$cntgCovAry_ref}) {
			if ($cntgCovAry_ref->[$i]) {
				my @tmpCovAry = split /,/, $cntgCovAry_ref->[$i];
				my $cov = $tmpCovAry[$aryIndex];
				if ($cov > 0) {
					my $pos = $i + 1;
					print WIGGLE join '', ((join "\t", ($pos, $cov)), "\n");
				}
			}
		}
	}
	close WIGGLE;
}
sub printmRNARefPtLog {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 9_analyzeNucleotideComposition|349
#	secondaryAppearInSection: >none
#	input: $mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $resultLogDir
#	output: 
#	toCall: &printmRNARefPtLog($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $resultLogDir);
#	calledInLine: 356
#....................................................................................................................................................#
	my ($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $resultLogDir) = @_;
	
	open GENEREFPTLOG, ">". "$resultLogDir/mRNA.ref.point.log.txt";
	my @headerAry = qw/mRNAID description locationTag strnd CDSLength/;
	foreach my $siteType (keys %{$mRNARefPtHsh_ref}) {
		push @headerAry, $siteType;
	}
	print GENEREFPTLOG join "", (join "\t", (@headerAry)), "\n";

	foreach my $mRNAID (sort keys %{$mRNAInfoHsh_ref}) {
		my $description = $mRNAInfoHsh_ref->{$mRNAID}{'description'};
		my $cntg = $mRNAInfoHsh_ref->{$mRNAID}{'cntg'};
		my $strnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
		my @CDSRngAry = sort {$a <=> $b} @{$mRNAInfoHsh_ref->{$mRNAID}{'CDSRng'}};
		my $locationTag = $cntg.":".$CDSRngAry[-1]."-".$CDSRngAry[0];
		my $CDSLength = $CDSRngAry[-1] - $CDSRngAry[0] + 1;
		my @outputAry = ($mRNAID, $description, $locationTag, $strnd, $CDSLength);
		foreach my $siteType (keys %{$mRNARefPtHsh_ref}) {
			my $pos = 'none';
			$pos = $mRNARefPtHsh_ref->{$siteType}{$mRNAID} if $mRNARefPtHsh_ref->{$siteType}{$mRNAID};
			push @outputAry, $pos;
		}
		print GENEREFPTLOG join "", (join "\t", (@outputAry)), "\n";
	}
	close GENEREFPTLOG;

	return ();
}
sub readGFF_oneRNAPerGene {
#....................................................................................................................................................#
#	subroutineCategory: gff
#	dependOnSub: currentTime|1071
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|273
#	secondaryAppearInSection: >none
#	input: $gffPath
#	output: $geneInfoHsh_ref
#	toCall: my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);
#	calledInLine: 285
#....................................................................................................................................................#

	my ($gffPath) = @_;

	my $geneInfoHsh_ref = {};
	
	#---read the gff
	my $geneByRNAHsh_ref = {};

	open (GFF, $gffPath);
	print "[".&currentTime()."] Reading: $gffPath\n";#->1071
	while (my $theLine = <GFF>) {

		chomp $theLine;
		
		last if $theLine =~ m/^##FASTA/;
		
		if ($theLine !~ m/^\#|^\@/ and $theLine !~ m/\tsupercontig\t/) {

			my ($seq, undef, $geneCategory, $featureStart, $featureEnd, undef, $geneStrd, undef, $dscrptns) = split (/\t/, $theLine);
			
			#----assigne all non -/+ will be treated as plus
			$geneStrd = "+" if (($geneStrd ne "-") and ($geneStrd ne "+"));
			
			my @dscrptnsSplt = split /;/, $dscrptns;
			my ($unqID, $parent);
			my $geneName = "unknown";
			foreach my $theDscptn (@dscrptnsSplt) {
				if ($theDscptn =~ m/^ID=/) {$unqID = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^Parent=/) {$parent = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^description=/) {$geneName = substr ($theDscptn, index ($theDscptn, "=")+1);}
			}

			if ($geneCategory eq "gene") {#---gene
				
				my $geneID = $unqID;
				
				$geneInfoHsh_ref->{$geneID}{'strnd'} = $geneStrd;
				$geneInfoHsh_ref->{$geneID}{'cntg'} = $seq;
				$geneInfoHsh_ref->{$geneID}{'description'} = uri_unescape($geneName);
				$geneInfoHsh_ref->{$geneID}{'description'} =~ s/\+/ /g;
				@{$geneInfoHsh_ref->{$geneID}{'geneRng'}} = ($featureStart, $featureEnd);

			} elsif ($geneCategory eq "CDS") {#---Only for coding genes
				
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'CDSRng'}}, ($featureStart, $featureEnd);
				
			} elsif ($geneCategory eq "exon") {#---exon, may be exons of alternative transcripts, wiull sort out later
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'exonRng'}}, ($featureStart, $featureEnd);
				
			} else {#---can be tRNA, rRNA, mRNA, repRNA, ncRNA
				my $RNAID = $unqID;
				my $geneID = $parent;
				$geneByRNAHsh_ref->{$RNAID} = $geneID;
				$geneInfoHsh_ref->{$geneID}{'ctgry'} = $geneCategory;
				@{$geneInfoHsh_ref->{$geneID}{'RNARng'}} = ($featureStart, $featureEnd);
				$geneInfoHsh_ref->{$geneID}{'RNAID'} = $RNAID;
			}
		}#---end of if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
	}#---end of while (my $theLine = <INFILE>)
	close GFF;
	
	#---get the UTR if any
	my $minUTRLength = 10;
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {

		if (exists $geneInfoHsh_ref->{$geneID}{'CDSRng'}) {
			my $exonMin = min(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $exonMax = max(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $CDSMin = min(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});
			my $CDSMax = max(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});

			if ($geneInfoHsh_ref->{$geneID}{'strnd'} eq '+') {
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			} else {
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			}
		}
	}
	
	return ($geneInfoHsh_ref);
}
sub readMultiFasta {
#....................................................................................................................................................#
#	subroutineCategory: fasta, general
#	dependOnSub: currentTime|1071
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|273
#	secondaryAppearInSection: >none
#	input: $fastaPath
#	output: $fastaHsh_ref
#	toCall: my ($fastaHsh_ref) = &readMultiFasta($fastaPath);
#	calledInLine: 282
#....................................................................................................................................................#

	my ($fastaPath) = @_;

	my ($seq, $seqName);
	my $fastaHsh_ref = {};
	my $i = 0;

	print "[".&currentTime()."] Reading: $fastaPath.\n";#->1071
	
	open (INFILE, $fastaPath);
	chomp (my $curntLine = <INFILE>); #get the first line
	while (my $nextLine = <INFILE>) {
		chomp $nextLine;
		
		#---Only two types of line in current line, the header or seq
		if ($curntLine =~ m/^>/) {#-- header line
			my @theLineSplt = split (/\|/, $curntLine);
			$seqName = $theLineSplt[0]; #---get the first tag
			$seqName =~ s/ //g; #---remove space
			$seqName =~ s/>//g; #---remove space
		} else {#--seq line
			$seq = $seq.$curntLine;
		}
		
		#---check if next line has a > or that's the end of file
		if ($nextLine =~ m/^>/) {
			$seq =~ tr/a-z/A-Z/;
			$fastaHsh_ref->{$seqName} = $seq;
			$seq = "";

			#---ad hoc limit
			#my $cntgNum = keys %{$fastaHsh_ref}; last if $cntgNum > 100;

		} elsif (eof(INFILE)) {#---this is the last line
			$seq =~ tr/a-z/A-Z/;
			$seq = $seq.$nextLine;
			$fastaHsh_ref->{$seqName} = $seq;
		}
		
		#---next line becomes current line
		$curntLine = $nextLine;
	}

	close INFILE;
	return ($fastaHsh_ref);
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|114
#	secondaryAppearInSection: >none
#	input: none
#	output: $STDRead5EndPileupIndxPathAndSizeAry_ref, $TEXBaseCompositionPileupIndxPath, $TEXRead5EndPileupIndxPathAndSizeAry_ref, $countCutoffPct, $fastaPath, $fullReadDefineNATPileupIndxPath, $gffPath, $outDir, $ratioCutoffPct, $trnsfrgInfoHshStorablePath
#	toCall: my ($TEXRead5EndPileupIndxPathAndSizeAry_ref, $STDRead5EndPileupIndxPathAndSizeAry_ref, $TEXBaseCompositionPileupIndxPath, $fullReadDefineNATPileupIndxPath, $trnsfrgInfoHshStorablePath, $fastaPath, $countCutoffPct, $ratioCutoffPct, $gffPath, $outDir) = &readParameters();
#	calledInLine: 123
#....................................................................................................................................................#
	
	my ($TEXRead5EndPileupIndxPathAndSizeAry_ref, $STDRead5EndPileupIndxPathAndSizeAry_ref, $TEXBaseCompositionPileupIndxPath, $fullReadDefineNATPileupIndxPath, $trnsfrgInfoHshStorablePath, $fastaPath, $gffPath, $countCutoffPct, $ratioCutoffPct, $outDir);

	my $dirPath = dirname(rel2abs($0));
	$outDir = "$dirPath/TSSFinder/";
	$countCutoffPct = 95;
	$ratioCutoffPct = 50;
	
	GetOptions 	("TEXRead5EndPileupIndxPathAndSize=s" => \@{$TEXRead5EndPileupIndxPathAndSizeAry_ref},
				 "STDRead5EndPileupIndxPathAndSize=s" => \@{$STDRead5EndPileupIndxPathAndSizeAry_ref},
				 "TEXBaseCompositionPileupIndxPath=s" => \$TEXBaseCompositionPileupIndxPath,
				 "fullReadDefineNATPileupIndxPath=s" => \$fullReadDefineNATPileupIndxPath,
				 "trnsfrgInfoHshStorablePath=s"  => \$trnsfrgInfoHshStorablePath,
				 "countCutoffPct:f"  => \$countCutoffPct,
				 "ratioCutoffPct:f"  => \$ratioCutoffPct,
				 "fastaPath=s"  => \$fastaPath,
				 "gffPath=s"  => \$gffPath,
				 "outDir:s"  => \$outDir)

	or die	("Error in command line arguments\n");
	
	#---check file
	foreach my $fileToCheck ($gffPath, $fastaPath, $trnsfrgInfoHshStorablePath, $TEXBaseCompositionPileupIndxPath, $fullReadDefineNATPileupIndxPath) {
		die "Can't read $fileToCheck" if (not -s $fileToCheck and not -s "$fileToCheck.gz");
	}

	system "mkdir -p -m 777 $outDir/";

	return($TEXRead5EndPileupIndxPathAndSizeAry_ref, $STDRead5EndPileupIndxPathAndSizeAry_ref, $TEXBaseCompositionPileupIndxPath, $fullReadDefineNATPileupIndxPath, $trnsfrgInfoHshStorablePath, $fastaPath, $countCutoffPct, $ratioCutoffPct, $gffPath, $outDir);
}
sub reportStatus {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|1071
#	appearInSub: analyzeSTDTEXRd5EndRatio|489, calculateBaseCompositionInAlignments|581, calculateNoiseIndex|625, checkOverlapAndProximity_withMargin|744, checkRunningThreadAndWaitToJoin|981, filterTEXSiteBasedOnTrainingCutoff|1275, findDistanceBetweenATGAndTrnsfrgEnd|1386, getBaseAtTSSAndExon|1699, getCoverageOfItemRngType_multiStrand|1897, getMaxCovAndPosAroundGeneEdge|2296, getRd5EndCountForExonAndTSSOfGoldenmRNASet|2451, getRd5EndSurroundingmRNAATG|2584, identifiyPeakSiteWithinRegion|3031, investigateTSSRelativeToATGAndTAA|3132, outputOriginalAndScaledRd5EndWiggle|3234, plotBaseCompositionAroundmRNAReferencePoint|3302, searchMotifAroundmRNAReferencePoint|4318, storeAndScaleRd5EndFromBothLibraries|4452, trainingUsingExonVsTssData|4548
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_analyzeTSSMotif|365, 11_outputWiggleXML|376, 6_trainningSTDTEXRatio|308, 7_filterTEXTrackAfterTrainning|319, 8_investigateGeneTSSRelativePos|331, 9_analyzeNucleotideComposition|349
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 524, 613, 647, 793, 807, 818, 1004, 1294, 1327, 1377, 1378, 1412, 1725, 1762, 1925, 1964, 1992, 2318, 2493, 2611, 3045, 3159, 3254, 3258, 3330, 4409, 4412, 4416, 4436, 4439, 4442, 4475, 4498, 4567
#....................................................................................................................................................#
	
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->1071

	return ();
	
}
sub runMAST {
#....................................................................................................................................................#
#	subroutineCategory: thridPartyApp
#	dependOnSub: >none
#	appearInSub: plotMotifOccurenceWithMASTMatch|3361
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_analyzeTSSMotif|365
#	input: $XMLPath, $fastaPath, $mastOutDir, $nonAnnotatedPolymerFreqPath
#	output: $mastHitLog
#	toCall: my ($mastHitLog) = &runMAST($XMLPath, $mastOutDir, $fastaPath, $nonAnnotatedPolymerFreqPath);
#	calledInLine: 3416, 3460
#....................................................................................................................................................#
	my ($XMLPath, $mastOutDir, $fastaPath, $nonAnnotatedPolymerFreqPath) = @_;
	
	my $mastHitLog = "$mastOutDir/mast.hit.log";
	my $mastErrorLog = "$mastOutDir/error.log";
	my $mastCmd = "mast $XMLPath $fastaPath -oc $mastOutDir -norc -ev 1e+10 -mt 10 -hit_list >$mastHitLog 2>$mastErrorLog";
	system ("$mastCmd");

	return ($mastHitLog);
}
sub searchMotifAroundmRNAReferencePoint {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: checkRunningThreadAndWaitToJoin|981, generateDREMECmd|1527, generateMEMECmd|1581, reportStatus|4274
#	appearInSub: >none
#	primaryAppearInSection: 10_analyzeTSSMotif|365
#	secondaryAppearInSection: >none
#	input: $forceRunMemeAndDreme, $hardCodedIGVPathHsh_ref, $memeDremeToRunHsh_ref, $resultDremeDir, $resultMemeDir, $seqAroundSiteHsh_ref
#	output: $dremeXMLPathHsh_ref, $memeXMLPathHsh_ref
#	toCall: my ($dremeXMLPathHsh_ref, $memeXMLPathHsh_ref) = &searchMotifAroundmRNAReferencePoint($seqAroundSiteHsh_ref, $resultDremeDir, $resultMemeDir, $forceRunMemeAndDreme, $hardCodedIGVPathHsh_ref, $memeDremeToRunHsh_ref);
#	calledInLine: 370
#....................................................................................................................................................#

	my ($seqAroundSiteHsh_ref, $resultDremeDir, $resultMemeDir, $forceRunMemeAndDreme, $hardCodedIGVPathHsh_ref, $memeDremeToRunHsh_ref) = @_;
	
	my $dremeXMLPathHsh_ref = {};
	my $memeXMLPathHsh_ref = {};

	my $minE = 0.001;
	my $memeSeqLimit = 1000;

	#---search sites only in mRNA_TSS
	#foreach my $siteType (qw /mRNA_TSS NAT_TSS/) {
	foreach my $siteType (qw /mRNA_TSS/) {

		my $positiveAlignHsh_ref = {};
		my $negativeAlignHsh_ref = {};
		my $tmpPosNegRngHsh_ref = {};
		my $tmpMaxMinKHsh_ref = {};
		
		#---manually define the range to be searched, 65 to 100 is the best after all

		$tmpPosNegRngHsh_ref->{'dnCore'}{'positive'}{'startPos'} = 120;
		$tmpPosNegRngHsh_ref->{'dnCore'}{'positive'}{'length'} = 20;
		$tmpPosNegRngHsh_ref->{'dnCore'}{'negative'}{'startPos'} = 0;
		$tmpPosNegRngHsh_ref->{'dnCore'}{'negative'}{'length'} = 20;
		$tmpMaxMinKHsh_ref->{'dnCore'}{'mink'} = 5;
		$tmpMaxMinKHsh_ref->{'dnCore'}{'maxk'} = 5;

		$tmpPosNegRngHsh_ref->{'TATA'}{'positive'}{'startPos'} = 40;
		$tmpPosNegRngHsh_ref->{'TATA'}{'positive'}{'length'} = 30;
		$tmpPosNegRngHsh_ref->{'TATA'}{'negative'}{'startPos'} = 170;
		$tmpPosNegRngHsh_ref->{'TATA'}{'negative'}{'length'} = 30;
		$tmpMaxMinKHsh_ref->{'TATA'}{'mink'} = 8;
		$tmpMaxMinKHsh_ref->{'TATA'}{'maxk'} = 8;

		$tmpPosNegRngHsh_ref->{'upCore'}{'positive'}{'startPos'} = 80;
		$tmpPosNegRngHsh_ref->{'upCore'}{'positive'}{'length'} = 15;
		$tmpPosNegRngHsh_ref->{'upCore'}{'negative'}{'startPos'} = 180;
		$tmpPosNegRngHsh_ref->{'upCore'}{'negative'}{'length'} = 15;
		$tmpMaxMinKHsh_ref->{'upCore'}{'mink'} = 6;
		$tmpMaxMinKHsh_ref->{'upCore'}{'maxk'} = 6;

		#$tmpPosNegRngHsh_ref->{'initiator'}{'positive'}{'startPos'} = 96;
		#$tmpPosNegRngHsh_ref->{'initiator'}{'positive'}{'length'} = 8;
		#$tmpPosNegRngHsh_ref->{'initiator'}{'negative'}{'startPos'} = 190;
		#$tmpPosNegRngHsh_ref->{'initiator'}{'negative'}{'length'} = 8;
		#$tmpMaxMinKHsh_ref->{'initiator'}{'mink'} = 6;
		#$tmpMaxMinKHsh_ref->{'initiator'}{'maxk'} = 6;

		foreach my $region (keys %{$tmpPosNegRngHsh_ref}) {
			my $mink = $tmpMaxMinKHsh_ref->{$region}{'mink'};
			my $maxk = $tmpMaxMinKHsh_ref->{$region}{'maxk'};
			
			foreach my $mRNAID (sort keys %{$seqAroundSiteHsh_ref->{$siteType}}) {
				my $fullSeq = $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'}.$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'};
				$positiveAlignHsh_ref->{$mRNAID} = substr $fullSeq, $tmpPosNegRngHsh_ref->{$region}{'positive'}{'startPos'}, $tmpPosNegRngHsh_ref->{$region}{'positive'}{'length'};
				$negativeAlignHsh_ref->{$mRNAID} = substr $fullSeq, $tmpPosNegRngHsh_ref->{$region}{'negative'}{'startPos'}, $tmpPosNegRngHsh_ref->{$region}{'negative'}{'length'};
			}

			if ($memeDremeToRunHsh_ref->{'dreme'}) {
				my @mkDirAry;
				my $dremeParentDir = "$resultDremeDir/$siteType/$region/"; push @mkDirAry, $dremeParentDir;
				my $dremeFastaDir = "$dremeParentDir/fasta"; push @mkDirAry, $dremeFastaDir;
				my $dremeLogDir = "$dremeParentDir/log"; push @mkDirAry, $dremeLogDir;
				my $dremeOutDir = "$dremeParentDir/out"; push @mkDirAry, $dremeLogDir;
				my $dremeCmdDir = "$dremeParentDir/cmd"; push @mkDirAry, $dremeCmdDir;
				foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}

				my $posSeqPath = "$dremeFastaDir/positive.fasta";
				my $negSeqPath = "$dremeFastaDir/negative.fasta";
				my $shuffleLogPath = "$dremeLogDir/shuffle.log.txt";
				my $negLogPath = "$dremeLogDir/negative.log.txt";
				my $dremeMode = 'negative'; #---shuffle or negative or both
				my $dremeOutputHsh_ref = &generateDREMECmd($positiveAlignHsh_ref, $negativeAlignHsh_ref, $posSeqPath, $negSeqPath, $dremeOutDir, $shuffleLogPath, $negLogPath, $maxk, $mink, $minE, $dremeMode);#->1527
				
				foreach my $shuffleOrNegative (keys %{$dremeOutputHsh_ref}) {

					my $dremeXMLPath = $dremeOutputHsh_ref->{$shuffleOrNegative}{'xml'};
					$dremeXMLPathHsh_ref->{$siteType}{$region}{$shuffleOrNegative} = $dremeXMLPath;
				
					if (not (-s $dremeXMLPath) or $forceRunMemeAndDreme eq 'yes') {
						&reportStatus("Issuing a Dreme thread to search for motifs around $region of $siteType", 10, "\n");#->4274
						threads->create(sub{system "echo \"$dremeOutputHsh_ref->{$shuffleOrNegative}{'cmd'}\" >$dremeCmdDir/dreme.$shuffleOrNegative.cmd; $dremeOutputHsh_ref->{$shuffleOrNegative}{'cmd'} ";});
					} else{
						&reportStatus("skipping dreme. XML found for $region $shuffleOrNegative $siteType", 10, "\n");#->4274
					}
				}
			} else {
				&reportStatus("skipping dreme for $siteType", 10, "\n");#->4274
			}
			
			if ($memeDremeToRunHsh_ref->{'meme'}) {
				my $seqAlignHsh_ref = $positiveAlignHsh_ref;
				my @mkDirAry;
				my $memeParentDir = "$resultMemeDir/$siteType/$region/"; push @mkDirAry, $memeParentDir;
				my $memeFastaDir = "$memeParentDir/fasta"; push @mkDirAry, $memeFastaDir;
				my $memeLogDir = "$memeParentDir/log"; push @mkDirAry, $memeLogDir;
				my $memeOutDir = "$memeParentDir/out"; push @mkDirAry, $memeOutDir;
				my $memeCmdDir = "$memeParentDir/out"; push @mkDirAry, $memeCmdDir;
				foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}
				my $seqPath = "$memeFastaDir/positive.fasta";
				my $logPath = "$memeLogDir/meme.log.txt";
				my $nonAnnotatedPolymerFreqPath = $hardCodedIGVPathHsh_ref->{'nonAnnotatedPolymerFreqPath'};
				my ($cmd, $XMLPath) = &generateMEMECmd($seqAlignHsh_ref, $seqPath, $memeOutDir, $memeSeqLimit, $mink, $maxk, $logPath, $nonAnnotatedPolymerFreqPath);#->1581
				
				$memeXMLPathHsh_ref->{$siteType}{$region} = $XMLPath;
				
				if (not (-s $XMLPath) or $forceRunMemeAndDreme eq 'yes' or $memeDremeToRunHsh_ref->{'meme'}) {
					&reportStatus("Issuing a meme thread to search for motifs around $region of $siteType", 10, "\n");#->4274
					threads->create(sub{system "echo \"$cmd\" >$memeCmdDir/meme.cmd; $cmd";});
				} else{
					&reportStatus("skipping meme. XML found for $region $siteType", 10, "\n");#->4274
				}
			} else {
				&reportStatus("skipping meme for $siteType", 10, "\n");#->4274
			}
		}
	}
	
	&checkRunningThreadAndWaitToJoin('yes', 1);#->981
	
	return ($dremeXMLPathHsh_ref, $memeXMLPathHsh_ref);
}
sub storeAndScaleRd5EndFromBothLibraries {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: checkRunningThreadAndWaitToJoin|981, generateThreadHshWithRandomItem|1645, getIndivCntgCovPlsPath|2146, reportStatus|4274, zipUnzipCntgCovInPlsPathHsh|4603
#	appearInSub: poolLibRd5EndCov|3584
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_processInputData|273
#	input: $fastaHsh_ref, $libTypeInfoHsh_ref, $maxThread, $rd5EndPlsInfoHsh_ref
#	output: none
#	toCall: &storeAndScaleRd5EndFromBothLibraries($rd5EndPlsInfoHsh_ref, $libTypeInfoHsh_ref, $fastaHsh_ref, $maxThread);
#	calledInLine: 3622
#....................................................................................................................................................#
	
	my ($rd5EndPlsInfoHsh_ref, $libTypeInfoHsh_ref, $fastaHsh_ref, $maxThread) = @_;
	
	my $tmpCntgCovInPlsPathHsh_refAry_ref = ();
	
	foreach my $TEXOrSTD (keys %{$libTypeInfoHsh_ref}) {
		foreach my $read5EndPileupIndxPath (@{$libTypeInfoHsh_ref->{$TEXOrSTD}{'read5EndPileupIndxPathAry_ref'}}) {
			my $cntgCovInPlsPathHsh_ref = &getIndivCntgCovPlsPath($read5EndPileupIndxPath);#->2146
			push @{$tmpCntgCovInPlsPathHsh_refAry_ref}, $cntgCovInPlsPathHsh_ref;
			push @{$libTypeInfoHsh_ref->{$TEXOrSTD}{'cntgCovInPlsPathHsh_ref'}}, $cntgCovInPlsPathHsh_ref;
		}
		my $numLib = @{$libTypeInfoHsh_ref->{$TEXOrSTD}{'cntgCovInPlsPathHsh_ref'}};
		&reportStatus("$numLib libraries store for $TEXOrSTD", 10, "\n");#->4274
	}

	#-----unzip the original lib Ary
	foreach my $cntgCovInPlsPathHsh_ref (@{$tmpCntgCovInPlsPathHsh_refAry_ref}) {&zipUnzipCntgCovInPlsPathHsh('unzip', $cntgCovInPlsPathHsh_ref);}#->4603

	#-----multithread pooling 
	my $itemAry_ref = [(keys %{$rd5EndPlsInfoHsh_ref->{'original'}{'covPlsPathHsh_ref'}})];
	my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomItem($maxThread, $itemAry_ref);#->1645
	my $cntgProc :shared = 0;
	my %threadHsh = ();

	foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
		my $cntgAry_ref = $randCntgInThreadHsh_ref->{$threadNum};
		($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781
			sub {

				my ($cntgAry_ref) = @_;

				foreach my $cntg (@{$cntgAry_ref}) {

					$cntgProc++;
					
					&reportStatus("Rd5End of $cntgProc cntgs were pooled and scaled", 10, "\r");#->4274

					my $originalRd5EndCntgCovAry_ref = retrieve($rd5EndPlsInfoHsh_ref->{'original'}{'covPlsPathHsh_ref'}->{$cntg});
					my $scaledRd5EndCntgCovAry_ref = retrieve($rd5EndPlsInfoHsh_ref->{'scaled'}{'covPlsPathHsh_ref'}->{$cntg});

					foreach my $TEXOrSTD (keys %{$libTypeInfoHsh_ref}) {

						my $scalingFactor = $libTypeInfoHsh_ref->{$TEXOrSTD}{'scalingFactor'};
						my $plusAryIndex = $libTypeInfoHsh_ref->{$TEXOrSTD}{'plus'}{'aryIndex'};
						my $minusAryIndex = $libTypeInfoHsh_ref->{$TEXOrSTD}{'minus'}{'aryIndex'};
			
						#---retieve the cntg cov ary for each of the TEXOrSTD libs
						foreach my $cntgCovInPlsPathHsh_ref (@{$libTypeInfoHsh_ref->{$TEXOrSTD}{'cntgCovInPlsPathHsh_ref'}}) {
							my $cntgCovPlsPath = $cntgCovInPlsPathHsh_ref->{$cntg};
							my $cntgCovAry_ref = retrieve($cntgCovPlsPath);
							for my $i (0..$#{$cntgCovAry_ref}) {
								if ($cntgCovAry_ref->[$i]) {
									my ($plusOriginalCov, $minusOriginalCov) = split /,/, $cntgCovAry_ref->[$i];

									my @tmpOriginalCovAry = (0,0,0,0);
									@tmpOriginalCovAry = split /,/, $originalRd5EndCntgCovAry_ref->[$i] if ($originalRd5EndCntgCovAry_ref->[$i]);
									$tmpOriginalCovAry[$plusAryIndex] += $plusOriginalCov;
									$tmpOriginalCovAry[$minusAryIndex] += $minusOriginalCov;
									$originalRd5EndCntgCovAry_ref->[$i] = join ",", @tmpOriginalCovAry;
						
									my $plusScaledCov = sprintf "%.2f", $plusOriginalCov*$scalingFactor;
									my $minusScaledCov = sprintf "%.2f", $minusOriginalCov*$scalingFactor;
									my @tmpScaledCovAry = (0,0,0,0);
									@tmpScaledCovAry = split /,/, $scaledRd5EndCntgCovAry_ref->[$i] if ($scaledRd5EndCntgCovAry_ref->[$i]);
									$tmpScaledCovAry[$plusAryIndex] += $plusScaledCov;
									$tmpScaledCovAry[$minusAryIndex] += $minusScaledCov;
									$scaledRd5EndCntgCovAry_ref->[$i] = join ",", @tmpScaledCovAry;
								}
							}
						}
					}
		
					store($originalRd5EndCntgCovAry_ref, $rd5EndPlsInfoHsh_ref->{'original'}{'covPlsPathHsh_ref'}->{$cntg});
					store($scaledRd5EndCntgCovAry_ref, $rd5EndPlsInfoHsh_ref->{'scaled'}{'covPlsPathHsh_ref'}->{$cntg});
				}
				
			}
			,($cntgAry_ref)
		); #---end of threads->new
	}#---end of foreach my $threadNum

	&checkRunningThreadAndWaitToJoin('no', 1);#->981

}
sub trainingUsingExonVsTssData {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: analyzeSTDTEXRd5EndRatio|489, getRd5EndCountForExonAndTSSOfGoldenmRNASet|2451, predictBonaFideTSSinGoldenmRNASet|3628, reportStatus|4274
#	appearInSub: >none
#	primaryAppearInSection: 6_trainningSTDTEXRatio|308
#	secondaryAppearInSection: >none
#	input: $TSSvsExonCutoffHshPath, $countCutoffPct, $exonEndTrim, $forceRunGoldenmRNASetTSSi, $ggplotDirHsh_ref, $goldenmRNASetSizeLimit, $libTypeInfoHsh_ref, $mRNAByCntgHsh_ref, $mRNAInfoHsh_ref, $maxPctDistBtwATGTrnsfrgEnd, $maxThread, $minOriginalCount, $ratioCutoffPct, $rd5EndPlsInfoHsh_ref, $trnsfrgInfoHsh_ref, $wigglePathHsh_ref
#	output: $TSSvsExonCutoffHsh_ref
#	toCall: my ($TSSvsExonCutoffHsh_ref) = &trainingUsingExonVsTssData($forceRunGoldenmRNASetTSSi, $trnsfrgInfoHsh_ref, $mRNAInfoHsh_ref, $ggplotDirHsh_ref, $maxPctDistBtwATGTrnsfrgEnd, $libTypeInfoHsh_ref, $mRNAByCntgHsh_ref, $wigglePathHsh_ref, $goldenmRNASetSizeLimit, $exonEndTrim, $rd5EndPlsInfoHsh_ref, $minOriginalCount, $TSSvsExonCutoffHshPath, $ratioCutoffPct, $countCutoffPct, $maxThread);
#	calledInLine: 314
#....................................................................................................................................................#

	my ($forceRunGoldenmRNASetTSSi, $trnsfrgInfoHsh_ref, $mRNAInfoHsh_ref, $ggplotDirHsh_ref, $maxPctDistBtwATGTrnsfrgEnd, $libTypeInfoHsh_ref, $mRNAByCntgHsh_ref, $wigglePathHsh_ref, $goldenmRNASetSizeLimit, $exonEndTrim, $rd5EndPlsInfoHsh_ref, $minOriginalCount, $TSSvsExonCutoffHshPath, $ratioCutoffPct, $countCutoffPct, $maxThread) = @_;
	
	my $TSSvsExonCutoffHsh_ref = {};
	
	#---find the cutoff hash or do the training
	if (-s $TSSvsExonCutoffHshPath) {

		&reportStatus("TSSvsExonCutoffHsh found and will be retrieved. Skipping training", 30, "\n");#->4274
		($TSSvsExonCutoffHsh_ref) = retrieve($TSSvsExonCutoffHshPath);

	} else {

		my ($goldenmRNATSSResultHsh_ref) = &predictBonaFideTSSinGoldenmRNASet($forceRunGoldenmRNASetTSSi, $trnsfrgInfoHsh_ref, $mRNAInfoHsh_ref, $ggplotDirHsh_ref, $maxPctDistBtwATGTrnsfrgEnd, $libTypeInfoHsh_ref, $mRNAByCntgHsh_ref, $rd5EndPlsInfoHsh_ref, $wigglePathHsh_ref, $goldenmRNASetSizeLimit, $maxThread);#->3628
		my $rd5EndTSSExonCountHsh_ref = &getRd5EndCountForExonAndTSSOfGoldenmRNASet($mRNAInfoHsh_ref, $goldenmRNATSSResultHsh_ref, $exonEndTrim, $libTypeInfoHsh_ref, $rd5EndPlsInfoHsh_ref, $minOriginalCount, $maxThread);#->2451
		($TSSvsExonCutoffHsh_ref) = &analyzeSTDTEXRd5EndRatio($rd5EndTSSExonCountHsh_ref, $ggplotDirHsh_ref, $minOriginalCount, $ratioCutoffPct, $countCutoffPct);#->489
		store($TSSvsExonCutoffHsh_ref, $TSSvsExonCutoffHshPath);

	}

	return $TSSvsExonCutoffHsh_ref;
	
}
sub warningToProceed {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: currentTime|1071
#	appearInSub: poolLibRd5EndCov|3584
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_processInputData|273
#	input: $message
#	output: none
#	toCall: &warningToProceed($message);
#	calledInLine: 3609
#....................................................................................................................................................#

	my ($message) = @_;
	print "[".&currentTime()."] $message\n";#->1071
	sleep 1;
	#print "[".&currentTime()."] press \"n\" to exit or press anything else to proceed:";#->1071
	#chomp (my $stdIn = <STDIN>);
	#die "exiting....." if $stdIn eq 'n';
}
sub zipUnzipCntgCovInPlsPathHsh {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: currentTime|1071
#	appearInSub: storeAndScaleRd5EndFromBothLibraries|4452
#	primaryAppearInSection: 5_processInputData|273
#	secondaryAppearInSection: >none
#	input: $cntgCovInPlsPathHsh_ref, $zipUnzip
#	output: none
#	toCall: &zipUnzipCntgCovInPlsPathHsh($zipUnzip, $cntgCovInPlsPathHsh_ref);
#	calledInLine: 303, 4479
#....................................................................................................................................................#

	my ($zipUnzip, $cntgCovInPlsPathHsh_ref) = @_;
	
	foreach my $cntg (sort keys %{$cntgCovInPlsPathHsh_ref}) {
		print "[".&currentTime()."] $zipUnzip cntg ary.                \r";#->1071
		my $cntgCovPlsPath = "$cntgCovInPlsPathHsh_ref->{$cntg}";
		if ($zipUnzip eq 'unzip') {
			system ("gzip -df $cntgCovPlsPath.gz") if (-s "$cntgCovPlsPath.gz");
		} else {
			system ("gzip -f $cntgCovPlsPath") if (-s "$cntgCovPlsPath");
		}
	}
	print "\n";
}

exit;
