#!/usr/bin/perl -w
#######################################################
# Author :  Jainy Thomas
# date   :  Jan 2017
# email  :  jainythomas1@gmail.com
# Pupose :  to genotype the locus of a novel TE insertion
#           
#####################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Cwd;
use Bio::SearchIO; 
use Bio::SeqIO;
use Bio::DB::Fasta;
use MIME::Lite;
use Data::Dumper;
use File::Copy;
use List::MoreUtils qw(uniq);
use List::MoreUtils qw(firstidx);


my $version = "14.0";
my $scriptname = "processbam_extract_GM_scores.pl";
my $changelog = "
#   - v1.0 = 16 Jan 2017 
#	- v2.0 = 24 Jan 2017
#				modified the script so that so that path of the output can be given at the commandline
#				using the picardtools version 2.8.1 and java 1.8
#	- v2.1 = 1  Feb 2017
#				included concatenate function for contigs and singlets
#	- v2.5 = 1  Feb 2017
#				the output will be directed to a folder with region within each ltr
#				included the verbose option
#   - v2.6 = 13 Feb 2017
#				can be selected in the script only IGV run or extracted reads needs to be done or both	
#   - v4.0 = 22 March 2017
#               extract genomic sequences and do blast with assembled contigs, rewrite the contigs in the orientation and align with genomic sequence   
#				To do: select blast can be performed with genomic sequence or consensus
#   - v5.0 = 19 April 2017
#               genomic sequence will be extracted for a shorter flanking from the breakpoint 
#				reversing the query and database options in blast
#	- v6.0 = 3 	May 2017
#				used subroutine renaming the fasta headers with filename. It helps to create one blast output file. will be easier to recognise the query sequences
#	- v6.2 = 8 	May 2017
#				trying to extract the discordant reads, but the current methods extract other unmapped reads as well that are in proper pair
#				going to try a different methods that extract only discordant reads (added notes on different stuff)
#   - v6.3 = 10 May 2017
#				Modify the script to use a single methods out of the above 
#				To do : add the blast option for TEs , accept the TE sequences (fasta format) and add a column to the file (TE name)
#	- v7.0 = 11 May 2017
#   - v8.0 = 20 June 2017
#     			changed name from Run_samtools_picard_cap3_selection_blast_genotyping.v7.0.pl  to find_genotypes.v8.0.pl
#				Rewrite the script the files are read by picard tools only once
#				
#   - 8.2 = 6 July 2017
# 				retried to see if there is any difference in acqusition of discordant reads No difference identified
#				fixed bug. script failed if no discordant reads were identified for an individual
#   - 8.3 = 7 July 2017
#				removed the different methods for obtaining the discordant reads
#	- 8.5 =18 July 2017
#				identify and soft clipped reads remove it for genome blasting
#   - 8.7 = 21 July 2017
#				introduced UrQt to trim reads, but works only locally the script terminates abruptly
#   - 8.8 = 26	July 2017
#               copying the fastq files locally for the UrQt to process
# 				modifying the script so that script can be used only for IGV purposes
#	- 8.9 = 28 July 2017
#				removed UrQt as the script get killed without any notice due to microcuts in the program
#	-8.9_1 =23 August 2017
#			    introduced reconstructTEfolder for concatenating allreads from all individual for TE reconstruction	
#   -8.9_2 = 26 September 2017
#				deprecated contigTE blast and disread TEblast as clement is doing with the the concatenated read TE
#				modified the path with tetype to genomeloc 
#	-11.0 = 2 Novemeber 2017
#				deprecate extracting only soft clipped reads
#	-12.0 = 25 January 2018
#				introducing splitting mapped reads to R1, R2 and unpaired reads for downstream pipelines (orientTE and extract TE)
#				also need to split discordant reads to corresponding R1, R2 and Up files
#	-13.0 = 31 March 2018
#				find_genotypes.v13.0_1KGP.pl current name
#				introduced mappability scores to the script that will help us to identify regions with low mapability and hence low confident insertions 		
#	-14.0 = 04-05-18
#				changed name to processbam_extract_GM_scores.pl	
#	-15.0 = 05 April 18
#				changing the software locations as commandline arguments for ease of integrating into the pipeline
#				all the optional arguments are changing to mandatory options except for the mappability option		
#	-v15.5 = 7 September 2018
#				added discoassembly folder and its contents, 
#	-v16.0 = 18 Septemeber 2018
#				tried different methods to extract soft clipped reads
#	-v17.0 = 18 September 2018
#				finalised one and removed others
#	-v18.0 = 18 January 2019
#				changed how reads are identified from array for one subroutie splitdismatepairs
#				introduced quality cut off while acquiring reads
#				To do: picard tools some times generated duplicate entry for reads with truncated reads that will cause read pairs to get generated in non-equal way, so need to fix that issue
#	-v19.0 = 26 February 2019
#				it is not the issue with the picard tools,but a bug in the some bam files in the script, fixed that			
#		
\n";

my $usage = "\nUsage [$version]: 
    perl $scriptname -t <table> -f <split files> -g <pathtogenomefile> -bl <bam location> -pt <path to picard tools> -sq <path of seqtk> -bu <path of bamutils>[-p <path of the outputdirectory processbamout>] -bt <path to bedtools bin directory> [-v] [-c] [-h] [-s] 
	
	
	
    MANDATORY ARGUMENT:	
    -t,--table 			(STRING) file contain accession information first column needs to be the IDs, second column BAMIDs
    -f,--file  			(STRING) file containing accession information output from splitfile/or make list
    -g,--genome       	(STRING) path to the genome file
    -bl,--bam location 	(STRING) location of bam files
    
    -pt,--picardtools   (STRING) path to picard tools
    OPTIONAL ARGUMENTS:
    -p,--path         	(STRING) output directory name (path)
                            	 Default = <current working directory>
    
    -c,--chlog  		(BOOL)   Print updates
    -v,--v      		(BOOL)   Print version if only option
    -s,--verbose		(BOOL)   The script will talk to you
    -h,--help>  		(BOOL)   Print this usage\n\n";
   


#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($file,$path,$table,$GENOME,$bamlocation,$seqtkpro,$bamUtilpro,$picardtools,$verbose,$help,$v,$chlog,$bedtoolsdir);
GetOptions ('f=s' => \$file,
            'p=s' => \$path,
            'g=s' => \$GENOME,
            't=s' => \$table,
            'bl=s'=> \$bamlocation,
            'pt=s'=> \$picardtools,
            'sq=s'=> \$seqtkpro,
            'bu=s'=> \$bamUtilpro, 
            'bt=s'=> \$bedtoolsdir,	
            'c'   => \$chlog, 
            'h'   => \$help,
            's'   => \$verbose, 
            'v'   => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script $scriptname version $version\n\n" if ((! $table) && (! $file) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $table) || (! $file) ||  ($help));
my $cwd = getcwd();
#die $usage if (($mappability eq "yes") && ((! $mysqldb) || (! $user) ||  ($password))) ;
$path = $cwd if (!$path) ;
#my $mapscores = "$path/file.mappabilityscores.txt";

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
#neccessary tools #added to path
#my $seqtkpro = "/home/jainy/software/seqtk";#vader server
#my $picardpro = "/home/jainy/software/picard-2.9.2";#vader server
#my $bamlocation = "/vbod2/cgoubert/Correct_Genotypes/1KGP_bams";#vader server
#my $bamUtilpro = "/home/jainy/software/bamUtil";#Yodaserver installed github version
#my $bedtoolspro = "/home/jainy/software/bedtools2/bin";#vaderserver,Yodaserver
#my $mysqltable = "hg19wgEncodeCrgMapabilityAlign100mer_index";



my %bamfile = ();
my $bamid ;
my @dbIDs;
my $db;
my $extractgenomicseqout;
my $uniqueid;
my $dispath;
my $individual;
my $genomeloc;
my @alldisreads =();
my %hashindividual =();
my %discordantmatelist =();
my $discordmatefile;
my $indiallreads;
my $dismateseqfilepath;
my $dismateIDlistpath;
my %cordinatescore;
my %allindimappingscores;
my $dbh;

#loading bam file ids for the individuals
%bamfile = load_file ($table);

open (my $fh, "<", $file) or confess "\n ERROR (main): could not open to read $file $!\n";
	while(<$fh>) {
		chomp (my $line = $_);	#NA06985    14_85161024 14 85161024 AluYc1 + null
		my @col = split(/\s+/,$line);
		$individual = $col[0];
		my $chr = $col[2] ;
		my $start = $col[3] - 250;
		my $end = $col[3] + 250;
		my $genomeloc = $chr.":".$start."-".$end;
		#my $m5start = $col[3] - 500;
		#my $m5end = $col[3] + 500;
		#$genomeloc =  $chr.":".$m5start."-".$m5end;
		#$cordinates500{$genomeloc} = 1;
		#$cordinates250{$genomeloc2} = 1;
		
		#my $tetype = $col[3];
		&find_bamid();
		$uniqueid = $individual.".".$genomeloc;
		
		#for viewing in IGV- step1
		make_path  ("$path/IGV/$genomeloc");
		unless (-e "$path/IGV/$genomeloc/$uniqueid.bam") {
			system("samtools view -b -o $path/IGV/$genomeloc/$uniqueid.bam $bamlocation/$bamid $genomeloc") == 0 or die ("unable to run command on $uniqueid \n");
		}
		print STDERR " 	Extracting 	$uniqueid using samtools done\n" if ($verbose);
		system ("samtools index -b $path/IGV/$genomeloc/$uniqueid.bam ") == 0 or die ("unable to create index file of $uniqueid.bam \n");
		print STDERR " indexing the bamfile done \n" if ($verbose) ;
		
		##Extracting discordant reads
		make_path  ("$path/Discordantreads/$individual");
		system ("samtools view -q 20 -b -F 3854 $bamlocation/$bamid $genomeloc > $path/Discordantreads/$individual/$uniqueid.discordantF3854.outfile.bam") == 0 or die ("unable to extract readsby F 3854 flag $!");
		system ("samtools bam2fq $path/Discordantreads/$individual/$uniqueid.discordantF3854.outfile.bam | $seqtkpro/seqtk seq -A -q20 > $path/Discordantreads/$individual/$uniqueid.dismapped.reads.fasta") == 0 or die ("unable to convert discordant bam file  to fasta $uniqueid \n");
		#copy to make discassembly
		make_path ("$path/Discoassembly/$genomeloc");
		#copy("$path/Discordantreads/$individual/$uniqueid.dismapped.reads.fasta", "$path/Discoassembly/$genomeloc/$uniqueid.dismapped.reads.fasta") or die "Copy failed:$!";
		###########################################
		#Trying to see which method extract 
		
		make_path ("$path/splitreads/$genomeloc");
		#extracting split reads(without -M during bwa mem)
		
		system ("samtools view -q 20 -b -f 2048 $bamlocation/$bamid $genomeloc > $path/splitreads/$genomeloc/$uniqueid.splitF2048.outfile.bam") == 0 or die ("unable to extract readsby F 3854 flag $!");
		system ("samtools bam2fq $path/splitreads/$genomeloc/$uniqueid.splitF2048.outfile.bam | $seqtkpro/seqtk seq -A -q20 > $path/splitreads/$genomeloc/$uniqueid.splitF2048.reads.fasta") == 0 or die ("unable to convert discordant bam file  to fasta $uniqueid \n");
		
		###########################################
		#identifying mates of discordant reads
		$dispath = "$path/Discordantreads/$individual";
		my $mappedisreads = "$uniqueid.dismapped.reads.fasta";
		%discordantmatelist = &load_readIDs ($dispath,$mappedisreads);
		#print discordant read mates to a file
		make_path  ("$path/Discordantreads/$file/dismateIDLists");	
		$discordmatefile = "$path/Discordantreads/$file/dismateIDLists/$uniqueid.discordantmatesreadIDlist.txt";
		&print_hash(%discordantmatelist);
		#loading the all the discordant reads that needs to be extracted to a hash
		&collectreads_indi();
	}
#print Dumper %discordantmatelist, "\n";	
#print Dumper %hashindividual, "\n";			
close $fh;
#print Dumper %allindimappingscores,"\n";



make_path  ("$path/Allreadsforbam/$file");
$indiallreads = "$path/Allreadsforbam/$file"; 
&printreads_indi();
&extractreads_bam();
make_path ("$path/Discordantreads/$file/discordantmatesonly");
$dismateseqfilepath = "$path/Discordantreads/$file/discordantmatesonly"; 
$dismateIDlistpath = "$path/Discordantreads/$file/dismateIDLists";
&extractdiscordmates();


# index the genome and connect to the fasta file
my $reindex;
my $indexfile = "$GENOME.index";
	if (-e $indexfile) {
		$reindex = 0;
		print "\t Genome previously indexed - Skipping indexing...\n\t (if you want to reindex the genome, delete $indexfile)\n";
	} else {
		$reindex = 1;
		print "\t  Genome not indexed - Indexing $GENOME...\n";
	}
$db = Bio::DB::Fasta->new( $GENOME, -reindex=>$reindex) or die "\t    ERROR - Failed to create Bio::DB::Fasta object from $GENOME $!\n";
#create list of the ID of the genome file  
@dbIDs = $db->get_all_ids();

#my $tetype;
#my $gextract;
my $allreadrenamedfile;
my $disrenamedfile;
#my $softreads;
#my %softlist;
my $allreadnosoftrenamedfile;

open ($fh, "<", $file) or confess "\n ERROR (main): could not open to read $file $!\n";
	while(<$fh>) {
		chomp (my $line = $_);
		my @col = split(/\s+/,$line);
		my $numcol = @col;
		print "the number of colums is $numcol\n";
		$individual = $col[0];
		my $chr = $col[2] ;
		#my $gstart = $col[2] - 20;
		#my $gend = $col[2] + 20;
		#$gextract = $chr."_".$gstart."-".$gend;
		my $start = $col[3] - 250;
		my $end = $col[3] + 250;
		#my $start2 = $col[2] - 500;
		#my $end2 = $col[2] + 500;
		#$genomeloc =  $chr.":".$start2."-".$end2;
		$genomeloc = $chr.":".$start."-".$end;
		#$tetype = $col[3];
		&find_bamid();
		$uniqueid = $individual.".".$genomeloc;
		my $epath = "$path/ExtractedReads/$genomeloc";
		my $dmatefile;
		make_path  ("$path/Discosplitassembly/$genomeloc") ;		
		#Extracting reads for the assembling the reads
	
		make_path  ("$path/ExtractedReads/$genomeloc") ;		
		#extracting only the mapped reads
		unless (-e "$path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.Up.fasta") {
			system("samtools view -q 20 -b -F 4 $bamlocation/$bamid $genomeloc > $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.bam ") == 0 or die ("unable to extract bam at $uniqueid \n");
			#to extract softclipped reads
			#my $cmd = q(awk '$6 ~ /S/ {print">"$1; print $10}');#to capture all reads that are softclipped
			# to capture soft clipped reads and add the read pair info to each read#optical duplicates,supplementaryalignments,not primaryalignment and qcfailed reads are not taken into account here see table /Users/jainy/Desktop/Projects/scripts/myscripts/samtoolsflag.txt
			my $cmd = q(awk '$6 ~ /S/ && $2 <= 127 && $2 >= 64 {print">"$1"/1\n"$10}; $6 ~ /S/ && $2 <= 191 && $2 >= 128 {print">"$1"/2\n"$10}');
			system ("samtools view $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.bam | $cmd > $path/splitreads/$genomeloc/$uniqueid.onlymappedreadIDs.bam.softclippedreads.fasta") == 0 or die ("unable to extract softclipped read names at $uniqueid \n");
			
			# to extract soft clipped reads in fasta example
			#samtools view B3.2:62211518-62212018.onlymappedreadIDs.bam | awk '$6 ~ /S/ {print">"$1; print $10}' > B3.2:62211518-62212018.softclippedreads.fasta
			system ("samtools bam2fq $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.bam | $seqtkpro/seqtk seq -A -q20  > $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.fasta") == 0 or die ("unable to convert to fasta $uniqueid \n");
			#sort the bam file on name
			system ("samtools sort -n -o $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.bam $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.bam") == 0 or die ("unable to name sort bam file $uniqueid \n");
			#Extract mapped reads as R1, R2 and Unparied from the bam
			system ("$bamUtilpro/bam bam2FastQ --in $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.bam --firstOut $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDsBU.namesorted.R1.fastq --secondOut $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDsBU.namesorted.R2.fastq --unpairedOut $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.Up.fastq --noReverseComp") == 0 or die ("unable to extract reads using bamUtil $uniqueid \n");
			unlink ("$path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDsBU.namesorted.R1.fastq","$path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDsBU.namesorted.R2.fastq") or warn "Could not unlink : $!";
			#extracting the R1 and R2 from bedtools as there is a bug with bamutil for R1 and R2
			system ("$bedtoolsdir/bedtools bamtofastq -i $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.bam -fq $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.R1.fastq -fq2 $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.R2.fastq") == 0 or die ("unable to extract reads using bamUtil $uniqueid \n");
			#changing fastq to fasta
			system ("$seqtkpro/seqtk seq -A -q20 $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.R1.fastq > $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.R1.fasta") == 0 or die ("unable to run seqtk on R1 $uniqueid \n");
			system ("$seqtkpro/seqtk seq -A -q20 $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.R2.fastq > $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.R2.fasta") == 0 or die ("unable to run seqtk on R2 $uniqueid \n");
			system ("$seqtkpro/seqtk seq -A -q20 $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.Up.fastq > $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.Up.fasta") == 0 or die ("unable to run seqtk on Up $uniqueid \n");
		}
		
		if (-e "$path/Discordantreads/$file/discordantmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta") {
			unless (-z "$path/Discordantreads/$file/discordantmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta") {
				
				#copy the corresponding the discordant mate reads
				copy("$path/Discordantreads/$file/discordantmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta", "$path/ExtractedReads/$genomeloc/$uniqueid.discordantmatesreadIDlist.txt.fasta") or die "Copy failed:$!";

				$dmatefile = "$uniqueid.discordantmatesreadIDlist.txt.fasta";
				&splitdismatetopairs($epath,$dmatefile);
				#concatenate mapped and discordant reads that is split based on read pair info 
				system ("cat $epath/$dmatefile.r1 $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.R1.fasta > $path/ExtractedReads/$genomeloc/$uniqueid.allreads.namesorted.D1.R1.fasta") == 0 or die ("unable to concatenate discordantsplitreads $uniqueid \n");
				system ("cat $epath/$dmatefile.r2 $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.R2.fasta > $path/ExtractedReads/$genomeloc/$uniqueid.allreads.namesorted.D2.R2.fasta") == 0 or die ("unable to concatenate discordantsplitreads $uniqueid \n");

				#copy("$path/Discordantreads/$file/discordantmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta", "$path/Discoassembly/$genomeloc/$uniqueid.discordantmatesreadIDlist.txt.fasta") or die "Copy failed:$!"; 
				
				#print "concatenate discordant reads and its mates\n";
				system ("cat $path/Discordantreads/$individual/$uniqueid.dismapped.reads.fasta $path/Discordantreads/$file/discordantmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta > $path/Discoassembly/$genomeloc/$uniqueid.disreadsmates.fasta") == 0 or die ("unable to concatenate allreads for assembly5 $uniqueid $! \n");
				
				
				#concatenate mapped and discordant reads
				
				#system ("cat $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.fasta $path/ExtractedReads/$genomeloc/$uniqueid.discordantmatesreadIDlist.txt.fasta > $path/ExtractedReads/$genomeloc/$uniqueid.allpairedendread.fasta") == 0 or die ("unable to concatenate allreads for assembly $uniqueid $! \n");
			}
		} 
# 		else {
# 			#copy("$path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.fasta", "$path/ExtractedReads/$genomeloc/$uniqueid.allpairedendread.fasta") or die "Copy failed:$!"; 
# 		}
		#add the option if the prediction is present (will help to get the cases where TEs present and will increase the chance of having a good TE assembly)
		#copy the all reads file to a folder where assembly for the TEsequence is made when the prediction is positive
		#concatenate discordants reads mates and split reads
		my $discmatereads = "$path/Discoassembly/$genomeloc/$uniqueid.disreadsmates.fasta";
		my $softclipedrds = "$path/splitreads/$genomeloc/$uniqueid.onlymappedreadIDs.bam.softclippedreads.fasta";
		my $supplereads = "$path/splitreads/$genomeloc/$uniqueid.splitF2048.reads.fasta";
		my $discmatesize = -s $discmatereads if (-e $discmatereads);
		my $softsize = -s $softclipedrds if (-e $softclipedrds );
		my $supplesize = -s $supplereads if (-e $supplereads);
		$supplesize = 0 if (! -e $supplereads);
		$softsize = 0 if (! -e $softclipedrds);
		$discmatesize = 0 if (! -e $discmatereads);
		
		if (($softsize > 0) && ($supplesize > 0) && ($discmatesize > 0)) {
			system qq(cat $softclipedrds $supplereads $discmatereads >> $path/Discosplitassembly/$genomeloc/$uniqueid.supsoft.disreadsmates.fasta) == 0 or die ("unable to concatenate allreads for assembly1 $uniqueid $! \n");
		
		} elsif (($softsize > 0) && ($supplesize == 0) && ($discmatesize > 0)) {
		
			system ("cat $softclipedrds $path/Discoassembly/$genomeloc/$uniqueid.disreadsmates.fasta > $path/Discosplitassembly/$genomeloc/$uniqueid.supsoft.disreadsmates.fasta ") == 0 or die ("unable to concatenate softreads $uniqueid $! \n");
		} elsif  (($softsize == 0) && ($supplesize > 0) && ($discmatesize > 0)) {
			system ("cat $supplereads $path/Discoassembly/$genomeloc/$uniqueid.disreadsmates.fasta > $path/Discosplitassembly/$genomeloc/$uniqueid.supsoft.disreadsmates.fasta ") == 0 or die ("unable to concatenate allreads for assembly2 $uniqueid $! \n");
		} elsif (($softsize == 0) && ($supplesize == 0) && ($discmatesize > 0)) {
			copy("$path/Discoassembly/$genomeloc/$uniqueid.disreadsmates.fasta", "$path/Discosplitassembly/$genomeloc/$uniqueid.disreadsmates.fasta") or die "Copy failed:$!"; 
		
		} elsif (($softsize > 0) && ($supplesize == 0) && ($discmatesize == 0)){
			copy("$softclipedrds", "$path/Discosplitassembly/$genomeloc/$uniqueid.softonly.disreadsmates.fasta") or die ("unable to copy $uniqueid $! \n");
		
		} elsif (($softsize > 0) && ($supplesize > 0) && ($discmatesize == 0)) {
		
			system ("cat $softclipedrds $supplereads > $path/Discosplitassembly/$genomeloc/$uniqueid.supsoftonly.disreadsmates.fasta ") == 0 or die ("unable to concatenate allreads for assembly4 $uniqueid $! \n");
		}
				 
		
		make_path  ("$path/orientTE/$genomeloc");
		if (-e "$path/ExtractedReads/$genomeloc/$uniqueid.allreads.namesorted.D1.R1.fasta") {
			copy("$path/ExtractedReads/$genomeloc/$uniqueid.allreads.namesorted.D1.R1.fasta", "$path/orientTE/$genomeloc/$uniqueid.allreads.namesorted.D1.R1.fasta") or die "Copy failed D1. R1:$!";
		} else {
			copy("$path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.R1.fasta", "$path/orientTE/$genomeloc/$uniqueid.allreads.namesorted.R1.fasta") or die "Copy failed R1:$!";
		}
		if (-e "$path/ExtractedReads/$genomeloc/$uniqueid.allreads.namesorted.D2.R2.fasta") {
			copy("$path/ExtractedReads/$genomeloc/$uniqueid.allreads.namesorted.D2.R2.fasta", "$path/orientTE/$genomeloc/$uniqueid.allreads.namesorted.D2.R2.fasta") or die "Copy failed D1. R1:$!";
		} else {
			copy("$path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.R2.fasta", "$path/orientTE/$genomeloc/$uniqueid.allreads.namesorted.R2.fasta") or die "Copy failed R1:$!";
		}
		if (-e "$path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.D.Up.fasta") {
			copy("$path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.D.Up.fasta", "$path/orientTE/$genomeloc/$uniqueid.allreads.namesorted.D.Up.fasta") or die "Copy failed D. Up:$!";
		} elsif (-e "$path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.Up.fasta") {
			copy("$path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.Up.fasta", "$path/orientTE/$genomeloc/$uniqueid.allreads.namesorted.Up.fasta") or die "Copy failed Up:$!";
		}
		#extracting the genomic sequence for the orient TE script

		make_path  ("$path/ExtractGenomicsequences");
		$extractgenomicseqout = "$path/ExtractGenomicsequences/$genomeloc.extract.seq.fa";
		unless (-e $extractgenomicseqout) {
			my $out = Bio::SeqIO->newFh(-format => 'Fasta',
										-file   => ">$extractgenomicseqout") 
										or die "\t    ERROR - Failed to create SeqIO FH object from $extractgenomicseqout $!\n";
			my $log = "$extractgenomicseqout.log";
			open (LOG, ">>$log") or die "\t    ERROR - can not create log file $log $!\n";
			my $newId = join('_',$chr,$start,$end);	
				if ("@dbIDs" =~ m/(\S*)($chr)(\S*)/) {
					my $subSeq = $db->seq($chr,$start,$end);	# extract target sequence
					my $seqobj = Bio::Seq->new( -display_id => $newId,
													   -seq => $subSeq); # create object with target
					print $out $seqobj;	# print it out (in fasta format)
				} else {
					print LOG "$chr was not found in $GENOME\n";
				}
	}
	close LOG;
	}
close $fh;
&email;
exit;
#-----------------------------------------------------------------------------
#----------------------------------- SUB -------------------------------------
#-----------------------------------------------------------------------------
sub load_file {
	my ($file1) = @_;
	my %sbamfile;
	open (my $th, "<", $file1) or confess "\n ERROR (main): could not open to read $file1 $!\n";
		while (my $data = <$th>) { #reading the table
			chomp $data; #storing the table values to the $data and removing the extraline
			my @namebam = split(/\s+/,$data); # splitting the data based on tab and storing into the arrray
			my $name = $namebam[0];
			$sbamfile{$name} = $namebam[1]; #loading to the hash
		}
	return (%sbamfile);
	close $th;	
}
sub find_bamid {
	if (exists ($bamfile { $individual } ))  {
		$bamid = $bamfile{$individual};
		print STDERR " the file now analysing is $bamid \n" if ($verbose);
	}
	else {
		die ("bamid cannot be successfully captured! Please check the input files $! \n");
	}
}
sub load_readIDs {
	my ($dpath,$mappedid) = @_;
	open (my $idh,"<","$dpath/$mappedid") || die ("unable to open $mappedid $! \n ");
	my %idlist;
	while (<$idh>) {
		 my $head = $_;
		 chomp $head;
		if ($head =~ /^\>(.*)\/(\d)/) {
			$head = "$1\/2" if ($2 == 1);
			$head = "$1\/1" if ($2 == 2);
			$idlist{$head} =1;
		} elsif ($head =~ /^\>(.*)/)  {
		 	$head = "$1\/0";
		 	$idlist{$head} =1;
		} else { 
			next;
		}
	}
	return (%idlist);
	close $idh;
}
sub print_hash {
	my %hashtoprint = @_;
	open (my $hp,">","$discordmatefile") || die ("failed to open file to write discordant mate list $!\n");
	foreach my $mate (sort keys %hashtoprint) {
		print $hp "$mate\n";
	}
	close $hp;
}
sub collectreads_indi {#collect all reads for an individual
	@alldisreads =();
	foreach my $read (sort keys %discordantmatelist) {
		$read = substr $read, 0,-2;#remove /1or2 from the file
		push (@{$hashindividual{$individual}},$read);
		(@{$hashindividual{$individual}}) = uniq (@{$hashindividual{$individual}});#for making the array unique
	}
	return (%hashindividual);
}
sub printreads_indi {
	foreach my $indi (sort keys %hashindividual) {
		open (my $ih, ">","$indiallreads/$indi.allreadIDs.txt" ) or die ("cannot write file $indi.allreadIDs.txt $!\n");
		foreach my $reads (@{$hashindividual{$indi}}) {
			print $ih "$reads\n";
		}
		close $ih;
	}
}
sub extractreads_bam {
	my @indifiles = `ls $path/Allreadsforbam/$file`;
	my $bam_id;
	foreach my $indivi (@indifiles) {
		my @allreadsfile_name = split (/\./,$indivi);
		my $allreadindi = $allreadsfile_name[0];
		if (exists ($bamfile{$allreadindi})) {
			$bam_id = $bamfile{$allreadindi};
		} else {
			print STDERR "bam_id cannot be identified for $allreadindi\n";
		}
		#extract reads using picard tools
		unless (-e "$path/Allreadsforbam/$file/$allreadindi.allreadIDs.bam") {
			system ("java -jar $picardtools/picard.jar FilterSamReads INPUT=$bamlocation/$bam_id VALIDATION_STRINGENCY=LENIENT FILTER=includeReadList READ_LIST_FILE=$path/Allreadsforbam/$file/$allreadindi.allreadIDs.txt WRITE_READS_FILES=false OUTPUT=$path/Allreadsforbam/$file/$allreadindi.allreadIDs.bam") == 0 or die ("unable to run picard tools in $file on $allreadindi \n");
		}
			system ("samtools bam2fq $path/Allreadsforbam/$file/$allreadindi.allreadIDs.bam | $seqtkpro/seqtk seq -A -q20 > $path/Allreadsforbam/$file/$allreadindi.allreadIDs.fasta") == 0 or die ("unable to convert to fasta $allreadindi \n");
			#system ("bedtools bamtofastq -i $path/Allreadsforbam/$file/$allreadindi.allreadIDs.bam -fq $path/Allreadsforbam/$file/$allreadindi.allreadIDs.fq" ) == 0 or die ("unable to convert to fastq $allreadindi \n");#didnt add mate information, so did not work so switched to previous version
			#system ("$seqtkpro/seqtk seq -A -q20 $path/Allreadsforbam/$file/$allreadindi.allreadIDs.fq > $path/Allreadsforbam/$file/$allreadindi.allreadIDs.fasta") == 0 or die ("unable to convert to fasta $allreadindi \n");

		
	}
}
sub extractdiscordmates {
	my @dismates = `ls $dismateIDlistpath`;
	foreach my $matefile (@dismates) {
		chomp $matefile;
		my @matefilename = split (/\./,$matefile);
		my $indivifilename = $matefilename[0];
		my %mates = ();
		open (my $mh, "<", "$dismateIDlistpath/$matefile") or confess "\n ERROR (main): could not open to read $matefile $!\n";
		while (my $dataline = <$mh>) { 
			chomp $dataline; 
			$mates{$dataline} = 1; 
		}
		#print Dumper %mates, "\n";
		
		if (-e "$path/Allreadsforbam/$file/$indivifilename.allreadIDs.fasta") {
			my $readio_obj = Bio::SeqIO->new(-file 	 => "$path/Allreadsforbam/$file/$indivifilename.allreadIDs.fasta", 
											 -format => 'fasta') 
										 or die "\t    ERROR - Failed to create SeqIO FH object from $indivifilename.allreadIDs.fasta $!\n";  
			my $outreadio_obj = Bio::SeqIO->new(-file   => ">$dismateseqfilepath/$matefile.fasta",
												-format => 'fasta') 
											 or die "\t    ERROR - Failed to create SeqIO FH object from $dismateseqfilepath/$matefile.fasta $!\n";  
			while (my $seq = $readio_obj->next_seq() ){
				my $header = $seq->display_id;
				if (exists $mates{$header}) {
					$outreadio_obj->write_seq($seq);
				} else {
					next;
				}
			}
		}
	}
}
sub renameseq_filename {
	my ($contigfile,$fpath) = @_;
	make_path  ("$fpath/Renamedcontigs");
	open (my $bhout, ">","$fpath/Renamedcontigs/$uniqueid.rename.fasta") or die "\n ERROR (main): could not open to read  $!\n";
	open (my $bh, "<", "$fpath/$contigfile") or confess "\n ERROR (main): could not open to read $contigfile $!\n";
		while(my $dataline = <$bh>) {
			chomp($dataline);
			#print STDERR "$line\n";
				if ($dataline =~ m/^\>\w+\d+/) {
				$dataline =~ s/^\>(\w+\d+)/\>$uniqueid\.$1/;
				print $bhout "$dataline\n";
				}
				else {
				print $bhout "$dataline\n";
				}
		}
	close $bh;
	close $bhout;
}
sub splitdismatetopairs {
	my ($dmpath,$dmfile) = @_;
	my %dmlistR1 = ();
	my %dmlistR2 = ();
	my %upaddlist = ();
	my @seqArrayR1 = ();
	my @seqArrayR2 = ();
	my @seqArrayUp = ();
	my %listR1 = ();
	my %listR2 = ();
	my %dup =();
	
	my $rdioobj1 = Bio::SeqIO->new(-file 	 => "$path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.R1.fasta", 
									 -format => 'fasta') 
								 or die "\t    ERROR - Failed to create SeqIO FH object from $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.R1.fasta $!\n";  
	while (my $seq = $rdioobj1->next_seq() ){#loading the reads in the discordant reads into a hash and also into separate arrays
		my $header = $seq->display_id;
		if ($header =~ /^(.*)\/(\d)/) {
			if ($2 == 1) {
				$listR1{$1}=1;
				#push (@seqArrayR1,$seq);
			} elsif ($2 == 2) {
				$listR2{$1}=1;
				#push (@seqArrayR2,$seq);
			}
		} 
	}	
	
	my $rdioobj2 = Bio::SeqIO->new(-file 	 => "$path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.R2.fasta", 
									 -format => 'fasta') 
								 or die "\t    ERROR - Failed to create SeqIO FH object from $path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.R2.fasta $!\n";  
	while (my $seq = $rdioobj2->next_seq() ){#loading the reads in the discordant reads into a hash and also into separate arrays
		my $header = $seq->display_id;
		if ($header =~ /^(.*)\/(\d)/) {
			if ($2 == 1) {
				$listR1{$1}=1;
				#push (@seqArrayR1,$seq);
			} elsif ($2 == 2) {
				$listR2{$1}=1;
				#push (@seqArrayR2,$seq);
			}
		} 
	}									 
								 
	my $readio_obj = Bio::SeqIO->new(-file 	 => "$dmpath/$dmfile", 
									 -format => 'fasta') 
								 or die "\t    ERROR - Failed to create SeqIO FH object from $dmfile $!\n";  
	while (my $seq = $readio_obj->next_seq() ){#loading the reads in the discordant reads into a hash and also into separate arrays
		my $header = $seq->display_id;
		if ($header =~ /^(.*)\/(\d)/) {
			if ($2 == 1) {
				if (exists ($listR1{$1})) {
					next;
				} elsif (exists ($dmlistR1{$1})) {
					delete ($dmlistR1{$1});
					$dup{$1} =1;
					next;
					
				} else {
					$dmlistR1{$1}=1;
					push (@seqArrayR1,$seq);
				}
				
			} elsif ($2 == 2) {
				if (exists ($listR2{$1})) {
					next;
				} elsif (exists ($dmlistR2{$1})) {
					delete ($dmlistR2{$1});
					$dup{$1} =1;
					next;
					
				} else {
					$dmlistR2{$1}=1;
					push (@seqArrayR2,$seq);
				}
				
			}
		} 
	}
	
	
	my $DUpfile = "$uniqueid.onlymappedreadIDs.namesorted.D.Up.fasta";
	my $Upreadio_obj = Bio::SeqIO->new(-file 	 => "$path/ExtractedReads/$genomeloc/$uniqueid.onlymappedreadIDs.namesorted.Up.fasta", 
									 -format => 'fasta') 
								 or die "\t    ERROR - Failed to create SeqIO FH object from $dmfile $!\n";  
	while (my $seq = $Upreadio_obj->next_seq() ){
		my $header = $seq->display_id;
		if (exists $dup{$header}) {
			push (@seqArrayUp, $seq);
			next;
		}
		if (exists ($dmlistR1{$header})) {#check if the read in unpaired file is present in the discordant read1 list has remove it discordant list
			delete ($dmlistR1{$header});
			$seq->display_id("$header/2");#modify display_id add /1 or /2
			push (@seqArrayR2,$seq);
		} elsif (exists ($dmlistR2{$header})) {
			delete ($dmlistR2{$header});
			$seq->display_id("$header/1");
			push (@seqArrayR1,$seq);
		} else {
			push (@seqArrayUp, $seq);
			next;
		}
	}
	
	 
	if (%dmlistR1) {
		for (my $i = 0; $i <= ($#seqArrayR1); $i++) {
			#print "last index is $#seqArrayR1\n";
			my $readname = $seqArrayR1[$i]->display_id;
			#print "$readname is the index is $i\n";
			$readname =~ /^(.*)\/(\d)/;
			if (exists $dmlistR1{$1}){
				#print $seqArrayR1[$i]->display_id,"is the read name with index $i\n";
				my $upseq = splice (@seqArrayR1,$i,1);
				#print "the index is $i\n";
				#$read = substr $read, 0,-2;
				$upseq->display_id("$1");
				#print $upseq->display_id,"\n";
				push (@seqArrayUp, $upseq);
				#$i= -1;
				$i = $i-1;
			} else {
				next;
			}
		}	
	}
	
	if (%dmlistR2) {
		for (my $i = 0; $i <= $#seqArrayR2; $i++) {
			#print "last index is $#seqArrayR2\n";
			my $readname = $seqArrayR2[$i]->display_id;
			#print "$readname is the index is $i\n";
			$readname =~ /^(.*)\/(\d)/;
			if (exists $dmlistR2{$1}){
				#print $seqArrayR2[$i]->display_id,"\n";
				#print "the index is $i\n";
				$readname = substr $readname, 0,-2;
				my $upseq = splice (@seqArrayR2,$i,1);
				$upseq->display_id("$readname");
				#print $upseq->display_id,"\n";
				push (@seqArrayUp, $upseq);
				#$i=-1;
				$i = $i-1;
			} else {
				next;
			}
			#print $seqArrayR1[$i]->seq;
		}	
	}
	
	if (%dup){
		for (my $i = 0; $i <= $#seqArrayR1; $i++) {
			#print "last index is $#seqArrayR2\n";
			my $readname = $seqArrayR1[$i]->display_id;
			#print "$readname is the index is $i\n";
			$readname =~ /^(.*)\/(\d)/;
			if (exists $dup{$1}){
				#print $seqArrayR2[$i]->display_id,"\n";
				#print "the index is $i\n";
				$readname = substr $readname, 0,-2;
				my $upseq = splice (@seqArrayR1,$i,1);
				$upseq->display_id("$readname");
				#print $upseq->display_id,"\n";
				push (@seqArrayUp, $upseq);
				#$i=-1;
				$i = $i-1;
			} else {
				next;
			}
			#print $seqArrayR1[$i]->seq;
		}	
	}
	if (%dup) {
		for (my $i = 0; $i <= $#seqArrayR2; $i++) {
			#print "last index is $#seqArrayR2\n";
			my $readname = $seqArrayR2[$i]->display_id;
			#print "$readname is the index is $i\n";
			$readname =~ /^(.*)\/(\d)/;
			if (exists $dup{$1}){
				#print $seqArrayR2[$i]->display_id,"\n";
				#print "the index is $i\n";
				$readname = substr $readname, 0,-2;
				my $upseq = splice (@seqArrayR2,$i,1);
				$upseq->display_id("$readname");
				#print $upseq->display_id,"\n";
				push (@seqArrayUp, $upseq);
				#$i=-1;
				$i = $i-1;
			} else {
				next;
			}
			#print $seqArrayR1[$i]->seq;
		}	
	}
	
	#instead of writing a new file open the R1 file and check if its already present or not and then write a new file all the R1 reads, R2 reads, then remove concatenation step.
	
	
	&sortfasta_header($dmpath,"$dmfile.r1",@seqArrayR1);
	&sortfasta_header($dmpath,"$dmfile.r2",@seqArrayR2);
	&sortfasta_header($dmpath,"$DUpfile",@seqArrayUp);
	#print STDERR "$DUpfile printed in $dmpath \n";
}

sub sortfasta_header {
	my ($spath,$filetosort,@seqarray) = @_;
		
	my $sortedReadio_obj = Bio::SeqIO->new(-file   => ">$spath/$filetosort",
										      -format => 'fasta') 
									         or die "\t    ERROR - Failed to create SeqIO FH object from $spath/$filetosort $!\n";  
	@seqarray = sort { ($a->display_id cmp $b->display_id) } @seqarray;
	foreach my $seque (@seqarray) {
		$sortedReadio_obj->write_seq($seque);	
	}
}
sub email {
	my $to = 'jainythomas1@gmail.com';
	my $cc = 'jainyt@genetics.utah.edu';
	my $from = 'jainy@vader.genetics.utah.edu';
	my $subject = 'genotyping';
	my $message = "genotyping $file DONE..";

	my $msg = MIME::Lite->new(
							From     => $from,
							To       => $to,
							Cc       => $cc,
							Subject  => $subject,
							Data     => $message
							);	 
	$msg->send;
	print "Email Sent Successfully\n";
}