#!/usr/bin/perl -w
#######################################################
# Author :  Jainy Thomas
# date   :  November 2017
# email  :  jainythomas1@gmail.com
# Pupose :  concatenate all reads from a locus from all individuals,assemble and blast against TE and genomics seq, find orientation of the TE (plus or minus)
#			and extract TE sequences with TSD and without TSDs
#           
#####################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Cwd;
use Data::Dumper;
use Bio::SeqIO;
use Bio::DB::Fasta;
use File::Copy;
use List::MoreUtils qw(uniq);
use feature 'fc';
use MIME::Lite;


my $version = "12.6";
my $scriptname = "orientTE_extractTE.pl";
my $changelog = "
#   - v1.0 = 3 November 2017 
#	- v2.5 = 8 November 2017 
#			   	extract full length TE from the assembled sequences and prints to a file (didnt finish it)
#   - v3.0 = 25 January 2018
#				changed the assembly method from cap3 to other next gene seq assemblers 
#	- v4.0 = 31 January 2018
#				extraction of TEsequences also made possible
#	- v4.1 = 31 January 2018
#				reverse compliment the flank for if the TE is in the minus orientation when blasted
#	- v4.2 = 2 February 2018
#				changed how the TEsequences read by the script, from the table and from folder with fasta sequences	
#	- v4.5 = 5 February 2018
#				generates a output table with prediction and if a full lengthTE was able to extract (yes or no), output the fasta 
#				sequence with both TSDs in one file and the others in another file (header locus,queryTE,orientation, TE)
#	- v5.0 = 4 April 2018
#				modifies the output table with TSD and name of the sequence
#	- v6.0 = 2 August 2018
# 				modified the input table format for TE and position column 1 chr_brkp col 2 is TE
#	- v6.1 = 6 September 2018
#				fixed a bug when spade the output is  contigs.fasta and added the option of passing threads to spade , 
#	- v7.0 = 8 September 2018
#				making a separate assembly with discordant reads and its mate
#				To do. Add split reads to it
#   - v8.0 = 2 October 2018
#				Added the option of possibility of doing split reads
#				introduced a method to identify split reads and add the split reads
#				fixed bug when tried to split split reads as well (added teh subroutinsplitdismatesoftreads)
#				spade failed for some loci for all reads ( dont know) and produced contigs.fasta instead of scaffold.fasta
#				added the option to run minia if spade produced contigs.fasta 
#	- v9.0 = 3 October 2018
#				comparing three outputs from three different methods and generating a final output
#				To do: add length of the sequence as a criteria to pick the partial sequence
#
#	- v10.0 = 4 October 2018
#				picking the longest one from the partial and full lengths
#	- v11.0 = 13 November 2018
#				introducing check for strand printing only if all predictions point to one strand
#				add length parameter to extract TEseq
#	- v11.1 = 5 December 2018
#				changed from 55 to 45 kmer length of minia
#				modified for DC, if the folder exists in the Assembled_DC it wont rerun the assembly
#				for now assembly features are cancelled for SPDC
#				introduce sort reverse for all read assembly
#	- v11.3	= 22 December 2018
#				assembly features are turned on but only if the directory is absent
#				minimum TSD is changed to 7
#				and sequence extracted within TE is limited to 4
#	- v12.0 = 11 February 2019
#				redo the assembly if the query cov is less than particular amount, to make sure assembly is not hindering established for discosplitassmebly
#				also introduced the features where a list of regions can be provided to do only particular set of loci
#	- v12.1 = 14 February 2019
#				introduced similar type of need based assembly to disc assembly and orientTE
#				also planning to introduce only the second and third set of assembly only for loci that did not work
#	- v12.3 = 21 February 2019
#				going to remove dipspades as it may improve time. Since I dont see contig made from dipSPA in test I did, In addition, some times there is no output created and stalls to go others			
#
#	- v12.4 = 22 February 2019
#				changing the how the assembly is getting performed
#	- v12.5 = 22 February 2019
#				opening the output from split assembly and performs the second assembly only on the failed ones
#	- v12.6 = 7 March 2019
#				added the option to compare the output from files instead of the need to restart the whole assembly again, also modified, to minimum requirement to flanking 25 bps for detection of orientation 
#				
#				
\n";

my $usage = "\nUsage [$version]: 
    perl $scriptname  -d <directory orientTE> -g <directory> -t <TEdirectory>  -l <list of TEinsertion table> -dc <Discassembly>-ds <Discosplitassembly> -a yes / no [-fl dir list to examine or -dr] [-p <path of the outputdirectory>][-o <output file>] [-v] [-c] [-h] 
	
    MANDATORY ARGUMENTS:	
    -d,--dirallreads    (STRING) input directory with all reads for all loci <path of the directory needs to be given>
	-g,--genomedir  (STRING) input directory with the corresponding reference genome sequence   
    -t,--TEdir	  	(STRING) input directory with TE sequences (separate fasta files with location name)
    -l,--list		(STRING) input file containing the information on TE insertion
    -dc,--disdir	(STRING) input directory (discoassembly)
    -bp,--blastn   (STRING) location of blastn
    -cp,--CAP3     (STRING) location of cap3 assembler
    -mn,--minia    (STRING) location of minia
    -sp,--spade    (STRING) location of spade
    -ds,--displit   (STRING) input directory discsplitdir
    -fl,--listfile  (STRING)	input file containing directories
    -dr --directores (BOOL) if directories are used instead of list of loci
    -a,--auto (STRING) automatically goes to next step and does next assembly only for those that does not have a complete assembly yes or no
    -oc      (BOOL) if you just wanna compare the three output and generate a new summary file
    OPTIONAL ARGUMENTS:
    -p,--path   (STRING) output directory name (path)
                         Default = <current working directory>
    -o,--output (STRING) output file with strand prediction
    -te,--TEout (STRING) output file with Assembled TE sequences only if nearly full length
    -c,--chlog  (BOOL)   Print updates
    -v,--v      (BOOL)   Print version if only option
    -s,--verbose(BOOL)   The script will talk to you
    -h,--help>  (BOOL)   Print this usage\n\n";

#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($rdir,$gdir,$TEdir,$onlycompare,$listte,$auto,$listfile,$folder,$teout,$discdir,$discsplidir,$miniadir,$CAP3dir,$BLASTdir,$spadedir,$cpus,$path,$fineout,$verbose,$help,$v,$chlog);
GetOptions ('d=s' => \$rdir,
			'g=s' => \$gdir,
			't=s' => \$TEdir,
			'l=s' => \$listte,
            'p=s' => \$path,
            'o=s' => \$fineout,
           'te=s' => \$teout,
           'dc=s' => \$discdir,
           'ds=s' => \$discsplidir,
           'bp=s' => \$BLASTdir,
           'cp=s' => \$CAP3dir,
           'mn=s' => \$miniadir,
           'sp=s' => \$spadedir,
           'cu=s' => \$cpus,
           'dr'   => \$folder,
           'a=s'  => \$auto,
           'fl=s' => \$listfile,
           'oc'	  => \$onlycompare,		
            'c'   => \$chlog, 
            'h'   => \$help,
            's'   => \$verbose, 
            'v'   => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script $scriptname version $version\n\n" if ((! $help)  && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $TEdir)  || ($help));
my $cwd = getcwd();
$path = $cwd if (!$path) ;
$cpus = 3 if (! $cpus);
#make_path("$path");
$fineout = "genomeloc.strand.prediction.$version.txt" if (! $fineout);
$teout = "Assembled_TEsequences.$version.txt" if (! $teout);
#die "\n -d $rdir does not exist?\n\n"  if (! -e $rdir);
#die "\n -g $gdir does not exist?\n\n"  if (! -e $gdir);
#die "\n -t $TEdir does not exist?\n\n"  if (! -e $TEdir);
#$rdir = $1 if ($rdir =~ /^(.*)\/$/);
#$gdir = $1 if ($gdir =~ /^(.*)\/$/);
#$TEdir = $1 if ($TEdir =~ /^(.*)\/$/);

#my $logfile = "$path/$genomloc.log";
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
#change to commandline
#my $CAP3pro = "/home/jainy/software/CAP3";#vader server
#my $BLASTpro = "/home/jainy/software/ncbi-blast-2.6.0+/bin";
#my $miniapro = "/home/jainy/software/minia-v2.0.7-Source/build/bin";
#my $SPAdepro = "/home/jainy/software/SPAdes-3.11.1-Linux/bin";#vaderserver,Yodaserver

#-d "$path/Assembled_TEreads"?die "$path/Assembled_TEreads already exist? please delete/rename the existing folder\n":make_path ("$path/Assembled_TEreads");
my $genomloc;
my $directory;
my @zerofiles;
my $R1out;
my %TEinfo;
my $contigmaxname; 
my $contigreadnum;
my $totalreadnum; 
my $outpath;
my %strandprediction;
my %extractcontiginfo;
my %alldirectory;
my %discpredictionout;
my %allreadpredictionout;
my %splitdiscpredictionout;
my @allreadfullTE;
my @allreadpartTE;
my @discsplitfullTE;
my @discsplitpartTE;
my @discfullTE;
my @discpartTE;
my %fullTEallread;
my %partTEallread;
my %fullTEdiscsplit;
my %partTEdiscsplit;
my %fullTEdisc;
my %partTEdisc;
my @fullTE_allmethods;
my @partTE_allmethods;



#load TE info
%TEinfo = &load_file($listte) if ($listte);
#print Dumper %TEinfo;
#concatenate all the files in the directory 
my $left;
my $right;
my @matchingtsdlist;
my @full_lengTEs ;
my @partial_teseq;
my $dspath;
my $dreadmatefile;
my $DUpfile;
my $DmatefileR1;
my $DmatefileR2;
my $slgedcdir;
my @unassdispl;
my @unassdisc;




if ($discsplidir) {
	print "#####################starting Discordantreadsplitassmbly#####################\n";
	my @discorspli;
	@discorspli = &readlocifromfile($listfile) if ($listfile);
	@discorspli = `ls $discsplidir` if ($folder);
	
	foreach $slgedcdir (@discorspli) {
		#my $present = 0;
		my $success = 0;
		chomp ($slgedcdir);
		$alldirectory{$slgedcdir} =1;
		
		unless (-e "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.fasta") {
			make_path ("$path/Assembled_SPDCreads/$slgedcdir");
			my @disremafiles = glob ("'$discsplidir/$slgedcdir/*disreadsmates.fasta'");
			#unlink glob "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreads.fasta_k45_ma3.*" if (-e "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreads.fasta_k45_ma3.contigs.fa");
			#unlink glob "$path/Assembled_SPDCreads/$slgedcdir/Renamed_Assembledseq/*.fasta.*" if (-e "$path/Assembled_SPDCreads/$slgedcdir/Renamed_Assembledseq/$slgedcdir.rename.fasta");
			#unlink("$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.extract.seq.fa.tabular.blast.out") if (-e "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.extract.seq.fa.tabular.blast.out");

			my $discout = "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.fasta";
			&concatenatefiles($discout,@disremafiles);
		}
		$dspath = "$path/Assembled_SPDCreads/$slgedcdir";
		$dreadmatefile = "$slgedcdir.concate.disreadsmates.fasta";
		
		$DUpfile = "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.Up.fasta";
		$DmatefileR1 = "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.R1.fasta";
		$DmatefileR2 = "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.R2.fasta"; 
		my $readcount = &splitsoftdismatetopairs($dspath,$dreadmatefile);
		#print " the total readcount is $readcount\n";
	
		#check rom here
		my $queryte = &find_TEseq($slgedcdir);
		print STDERR "the query file is $queryte.fasta for the locus $slgedcdir\n";
	
		if (($readcount > 1) && ($readcount < 100)) {
			#unless ($present ==1) {
			system ("$CAP3dir/cap3 $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.fasta > $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.asmbl.fasta") == 0 or warn ("unable to assemble fasta $directory \n");
			#}
			if (-e "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.fasta.cap.contigs") {
				my $assembledfile = "$slgedcdir.concate.disreadsmates.fasta.cap.contigs";
				&renameseq_filename($assembledfile,$dspath,$slgedcdir,"cap");
			}
			
			&blast_file ($TEdir,"$queryte.fasta",$slgedcdir,$dspath,"cap");	#blast with the query TEsequence
			my $blastfile = "$dspath/$slgedcdir.$queryte.fasta.cap.tabular.blast.out";
			my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","cap");
			$qst = 0 if (not defined $qst);
			$qen = 0 if (not defined $qen);
			&blast_file ($gdir,"$slgedcdir.extract.seq.fa",$slgedcdir,$dspath,"cap");
			my $blastgenomefile = "$dspath/$slgedcdir.extract.seq.fa.cap.tabular.blast.out";
			my $genomeblast = &parse_blast($blastgenomefile,"95","cap");
			&compare_twoHOH($slgedcdir,$TEblast,$genomeblast);
		
		} elsif ($readcount >= 100) {
# 			unless ($present ==1) {
# 				system ("$spadedir/spades.py -1 $DmatefileR1 -2 $DmatefileR2 -s $DUpfile --careful --only-assembler -t $cpus -o $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreadsSPAdeout") == 0 or 
# 				#system ("$spadedir/dipspades.py -1 $DmatefileR1 -2 $DmatefileR2 -s $DUpfile -t $cpus --only-assembler --expect-rearrangements -o $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreadsdipSPAdeout") == 0 or 
# 				system ("$spadedir/spades.py --s1 $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.fasta --careful --only-assembler -t $cpus -o $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreadssingleSPAdeout") == 0 or 
# 				system ("$miniadir/minia -in $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreads.fasta_k45_ma3") == 0 or warn ("unable to assemble disfasta $slgedcdir \n");
# 			}
			unless (-e "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreadsSPAdeout/scaffolds.fasta" ) {
				system ("$spadedir/spades.py -1 $DmatefileR1 -2 $DmatefileR2 -s $DUpfile --careful --only-assembler -t $cpus -o $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreadsSPAdeout")== 0 or warn ("unable to assemble disfasta $slgedcdir using spades \n");
			}	
			if (-e "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreadsSPAdeout/scaffolds.fasta") {
				my %longestassem;
				#Rename scaffolds.fasta 
				#print "checking of spadesscaffpld";
				copy("$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreadsSPAdeout/scaffolds.fasta", "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreads.scaffolds.fasta") or die "Copy failed scaffolds.fasta $slgedcdir:$!";
				#rename the fasta sequence with its filename
				my $assembledfile = "$slgedcdir.disreads.scaffolds.fasta";
				&renameseq_filename($assembledfile,$dspath,$slgedcdir,"SPAscaf");
				&blast_file ($TEdir,"$queryte.fasta",$slgedcdir,$dspath,"SPAscaf");	#blast with the query TEsequence
				my $blastfile = "$dspath/$slgedcdir.$queryte.fasta.SPAscaf.tabular.blast.out";#$qloc.$type.tabular.blast.out
				my $blastfilesize = -s $blastfile;
				my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","SPAscaf") if ($blastfilesize > 0);
				$qst = 0 if (not defined $qst);
				$qen = 0 if (not defined $qen);
				print "$qst,$qen for $slgedcdir SPAscaf\n";
				my $TEblastSPAscaf = $TEblast;
				print Dumper $TEblastSPAscaf,"\n";
				$longestassem{'SPAscaf'} = ($qen - $qst) ;
				print "I was here station1\n";
				print "status of success is $success\n";
				if (($blastfilesize == 0) ||((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280))|| (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))   ) {
					system ("$spadedir/spades.py --s1 $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.fasta --careful --only-assembler -t $cpus -o $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreadssingleSPAdeout") == 0 or warn ("unable to assemble disfasta $slgedcdir \n");
					#system ("$miniadir/minia -in $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreads.fasta_k45_ma3") == 0 or warn ("unable to assemble disfasta $slgedcdir \n");
					if (-e "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreadssingleSPAdeout/scaffolds.fasta") {
						copy("$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreadssingleSPAdeout/scaffolds.fasta",	"$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreads.single_scaffolds.fasta") or die "Copy failed single_scaffolds.fasta $slgedcdir $!";
						my $assembledfile = "$slgedcdir.disreads.single_scaffolds.fasta";
						&renameseq_filename($assembledfile,$dspath,$slgedcdir,"singSPA");	
						&blast_file ($TEdir,"$queryte.fasta",$slgedcdir,$dspath,"singSPA");
						my $blastfile = "$dspath/$slgedcdir.$queryte.fasta.singSPA.tabular.blast.out";#$qloc.$type.tabular.blast.out
						my $blastfilesize = -s $blastfile;
						my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","singSPA") if ($blastfilesize > 0);
						$qst = 0 if (not defined $qst);
						$qen = 0 if (not defined $qen);
						print "$qst,$qen for $slgedcdir singSPA\n";
						my $TEblastsingSPA = $TEblast;
						print Dumper $TEblastsingSPA,"\n";
						$longestassem{'singSPA'} = ($qen - $qst) ;
						print "I was here station2\n";
						print "status of success is $success\n";
						if (($blastfilesize == 0) || ((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280))|| (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))   ) {
							system ("$miniadir/minia -in $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreads.fasta_k45_ma3") == 0 or warn ("unable to assemble disfasta $slgedcdir \n");	
							if (-e "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreads.fasta_k45_ma3.contigs.fa") { 
								my $assembledfile = "$slgedcdir.disreads.fasta_k45_ma3.contigs.fa";
								&renameseq_filename($assembledfile,$dspath,$slgedcdir,"mini");
								&blast_file ($TEdir,"$queryte.fasta",$slgedcdir,$dspath,"mini");
								my $blastfile = "$dspath/$slgedcdir.$queryte.fasta.mini.tabular.blast.out";#$qloc.$type.tabular.blast.out
								my $blastfilesize = -s $blastfile;
								my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","mini");
								$qst = 0 if (not defined $qst);
								$qen = 0 if (not defined $qen);
								print "$qst,$qen for $slgedcdir mini\n";
								my $TEblastmini = $TEblast;
								print Dumper $TEblastmini,"\n";
								$longestassem{'mini'} = ($qen - $qst) ;
								print "I was here station3\n";
								print "status of success is $success\n";
								if (($blastfilesize == 0) || ((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280))|| (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))   ) {
									my $longest = &findmaxvaluehash(%longestassem);
									&blast_file ($gdir,"$slgedcdir.extract.seq.fa",$slgedcdir,$dspath,$longest);
									my $blastgenomefile = "$dspath/$slgedcdir.extract.seq.fa.$longest.tabular.blast.out";
									my $genomeblast = &parse_blast($blastgenomefile,"95",$longest);
									&compare_twoHOH($slgedcdir,$TEblastmini,$genomeblast) if ($longest eq "mini");
									&compare_twoHOH($slgedcdir,$TEblastsingSPA,$genomeblast) if ($longest eq "singSPA");
									#&compare_twoHOH($slgedcdir,$TEblastdipSPA,$genomeblast) if ($longest eq "dipSPA");
									&compare_twoHOH($slgedcdir,$TEblastSPAscaf,$genomeblast) if ($longest eq "SPAscaf");
									print Dumper %$genomeblast, "is the genome blast\n";
									$success = 1;
									print "I was here station4\n";
									print "status of success is $success\n";
								} else {
									&blast_file ($gdir,"$slgedcdir.extract.seq.fa",$slgedcdir,$dspath,"mini");
									my $blastgenomefile = "$dspath/$slgedcdir.extract.seq.fa.mini.tabular.blast.out";
									my $genomeblast = &parse_blast($blastgenomefile,"95","mini");
									print Dumper %$genomeblast, "is the genome blast\n";
									&compare_twoHOH($slgedcdir,$TEblastmini,$genomeblast);
									$success = 1;
									print "I was here station5\n";
									print "status of success is $success\n";
								}
							}
						} else {
							&blast_file ($gdir,"$slgedcdir.extract.seq.fa",$slgedcdir,$dspath,"singSPA");
							my $blastgenomefile = "$dspath/$slgedcdir.extract.seq.fa.singSPA.tabular.blast.out";
							my $genomeblast = &parse_blast($blastgenomefile,"95","singSPA");
							print Dumper %$genomeblast, "is the genome blast\n";
							&compare_twoHOH($slgedcdir,$TEblastsingSPA,$genomeblast);
							$success = 1;
							print "I was here station6\n";
							print "status of success is $success\n";#was here
						}
					}
				} else {
					#parse the blastoutput for max score
					&blast_file ($gdir,"$slgedcdir.extract.seq.fa",$slgedcdir,$dspath,"SPAscaf");
					my $blastgenomefile = "$dspath/$slgedcdir.extract.seq.fa.SPAscaf.tabular.blast.out";
					my $genomeblast = &parse_blast($blastgenomefile,"95","SPAscaf");
					print Dumper %$genomeblast, "is the genome blast\n";
					&compare_twoHOH($slgedcdir,$TEblastSPAscaf,$genomeblast);
					$success = 1;
					print "I was here station7\n";
					print "status of success is $success\n";
				}  
				
			} 
			
			
			if ((! -e "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreadssingleSPAdeout/scaffolds.fasta" ) && ($success == 0)) {
				my %longestassem;
				system ("$spadedir/spades.py --s1 $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.fasta --careful --only-assembler -t $cpus -o $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreadssingleSPAdeout") == 0 or warn ("unable to assemble disfasta $slgedcdir using singspades \n");
				if (-e "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreadssingleSPAdeout/scaffolds.fasta") {
					copy("$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreadssingleSPAdeout/scaffolds.fasta",	"$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreads.single_scaffolds.fasta") or die "Copy failed single_scaffolds.fasta $slgedcdir $!";
					my $assembledfile = "$slgedcdir.disreads.single_scaffolds.fasta";
					&renameseq_filename($assembledfile,$dspath,$slgedcdir,"singSPA");	
					&blast_file ($TEdir,"$queryte.fasta",$slgedcdir,$dspath,"singSPA");
					my $blastfile = "$dspath/$slgedcdir.$queryte.fasta.singSPA.tabular.blast.out";#$qloc.$type.tabular.blast.out
					my $blastfilesize = -s $blastfile;
					my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","singSPA") if ($blastfilesize > 0);
					$qst = 0 if (not defined $qst);
					$qen = 0 if (not defined $qen);
					my $TEblastsingSPA = $TEblast;
					$longestassem{'singSPA'} = ($qen - $qst) ;
					print "I was here station8\n";
					print "status of success is $success\n";
				
					if (($blastfilesize == 0) || ((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280)) || (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))  ) {
						system ("$miniadir/minia -in $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreads.fasta_k45_ma3") == 0 or warn ("unable to assemble disfasta $slgedcdir \n");	
						if (-e "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreads.fasta_k45_ma3.contigs.fa") { 
							my $assembledfile = "$slgedcdir.disreads.fasta_k45_ma3.contigs.fa";
							&renameseq_filename($assembledfile,$dspath,$slgedcdir,"mini");
							&blast_file ($TEdir,"$queryte.fasta",$slgedcdir,$dspath,"mini");
							my $blastfile = "$dspath/$slgedcdir.$queryte.fasta.mini.tabular.blast.out";#$qloc.$type.tabular.blast.out
							my $blastfilesize = -s $blastfile;
							my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","mini");
							$qst = 0 if (not defined $qst);
							$qen = 0 if (not defined $qen);
							my $TEblastmini = $TEblast;
							$longestassem{'mini'} = ($qen - $qst) ;
							print "I was here station9\n";
							print "status of success is $success\n";
							if (($blastfilesize == 0) || ((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280)) || (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))  ) {
								my $longest = &findmaxvaluehash(%longestassem);
								&blast_file ($gdir,"$slgedcdir.extract.seq.fa",$slgedcdir,$dspath,$longest);
								my $blastgenomefile = "$dspath/$slgedcdir.extract.seq.fa.$longest.tabular.blast.out";
								my $genomeblast = &parse_blast($blastgenomefile,"95",$longest);
								#print "discsplit genomeblast\n";
								print Dumper %$genomeblast, "is the genome blast\n";
								#predict the orientation of TE insertion
								&compare_twoHOH($slgedcdir,$TEblastmini,$genomeblast) if ($longest eq "mini");
								&compare_twoHOH($slgedcdir,$TEblastsingSPA,$genomeblast) if ($longest eq "singSPA");
								$success = 1;
								print "I was here station10\n";
								print "status of success is $success\n";
							} else {
								&blast_file ($gdir,"$slgedcdir.extract.seq.fa",$slgedcdir,$dspath,"mini");
								my $blastgenomefile = "$dspath/$slgedcdir.extract.seq.fa.mini.tabular.blast.out";
								my $genomeblast = &parse_blast($blastgenomefile,"95","mini");
								#print "discsplit genomeblast\n";
								print Dumper %$genomeblast, "is the genome blast\n";
								#predict the orientation of TE insertion
								&compare_twoHOH($slgedcdir,$TEblastmini,$genomeblast);
								$success = 1;
								print "I was here station11\n";
								print "status of success is $success\n";
							}
						}
					} else {
						&blast_file ($gdir,"$slgedcdir.extract.seq.fa",$slgedcdir,$dspath,"singSPA");
						my $blastgenomefile = "$dspath/$slgedcdir.extract.seq.fa.singSPA.tabular.blast.out";
						my $genomeblast = &parse_blast($blastgenomefile,"95","singSPA");
						#print "discsplit genomeblast\n";
						print Dumper %$genomeblast, "is the genome blast\n";
						#predict the orientation of TE insertion
						&compare_twoHOH($slgedcdir,$TEblastsingSPA,$genomeblast);
						$success = 1;
						print "I was here station12\n";
						print "status of success is $success\n";
					}
				} 
			}
			
			if (($success == 0) && (! -e "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreads.fasta_k45_ma3.contigs.fa")) {
				my %longestassem;
				system ("$miniadir/minia -in $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.concate.disreadsmates.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreads.fasta_k45_ma3") == 0 or warn ("unable to assemble disfasta $slgedcdir \n");
				if (-e "$path/Assembled_SPDCreads/$slgedcdir/$slgedcdir.disreads.fasta_k45_ma3.contigs.fa") { 
					my $assembledfile = "$slgedcdir.disreads.fasta_k45_ma3.contigs.fa";
					&renameseq_filename($assembledfile,$dspath,$slgedcdir,"mini");
					&blast_file ($TEdir,"$queryte.fasta",$slgedcdir,$dspath,"mini");
					my $blastfile = "$dspath/$slgedcdir.$queryte.fasta.mini.tabular.blast.out";#$qloc.$type.tabular.blast.out
					my $blastfilesize = -s $blastfile;
					my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","mini") ;
					$qst = 0 if (not defined $qst);
					$qen = 0 if (not defined $qen);
					my $TEblastmini = $TEblast;
					$longestassem{'mini'} = ($qen - $qst) ;
					print "I was here station13\n";
					print "status of success is $success\n";
				
					if (($blastfilesize == 0) || ((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280)) || (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))  ) {
						my $longest = &findmaxvaluehash(%longestassem);
						&blast_file ($gdir,"$slgedcdir.extract.seq.fa",$slgedcdir,$dspath,$longest);
						my $blastgenomefile = "$dspath/$slgedcdir.extract.seq.fa.$longest.tabular.blast.out";
						my $genomeblast = &parse_blast($blastgenomefile,"95",$longest);
						#print "discsplit genomeblast\n";
						print Dumper %$genomeblast, "is the genome blast\n";
						#predict the orientation of TE insertion
						&compare_twoHOH($slgedcdir,$TEblastmini,$genomeblast) if ($longest eq "mini");
						print "I was here station14\n";
						print "status of success is $success\n";
					} else {
						&blast_file ($gdir,"$slgedcdir.extract.seq.fa",$slgedcdir,$dspath,"mini");
						my $blastgenomefile = "$dspath/$slgedcdir.extract.seq.fa.mini.tabular.blast.out";
						my $genomeblast = &parse_blast($blastgenomefile,"95","mini");
						#print "discsplit genomeblast\n";
						print Dumper %$genomeblast, "is the genome blast\n";
						#predict the orientation of TE insertion
						&compare_twoHOH($slgedcdir,$TEblastmini,$genomeblast);
						print "I was here station15\n";
						print "status of success is $success\n";
					
					}
				}
			}
			
			 
			
			if ($readcount==0) {
				push (@zerofiles,"$slgedcdir");
			}
		}
	
	}

	print Dumper %extractcontiginfo;
	#print Dumper %alldirectory;

	print  "DONE with prediction of orientation.... \n";

	print  "\n\nstarting to extract TEsequences.... \n";

	$left = 5;
	$right = 3;
	@matchingtsdlist=();
	@full_lengTEs=();
	@partial_teseq=();
	#my $teseq;
	&loadseq_tohash("Assembled_SPDCreads");
	print Dumper %extractcontiginfo;
	#find the contigs with TEsequence and flanking sequence and extract the flanking sequence
	foreach my $locus (sort keys %extractcontiginfo) {
		print "locus is $locus\n";#added they type of assembly so should be processed similarly
		my ($typeasm,$genomeloc,$rest) = split (/\./,$locus,3);
		my ($leftkmer);
		my ($rightkmer);
		my $leftflankseq;
		my $rightflankseq;
		my $start;
		my $end;
		my $strand;
		my $seq;
		my $qTE;
		my $newlocus;
		my $teseq;
		my $testart;
		my $teend;
		my $tsdseq;
		my $tsdsucc;
		my $contlen;
		foreach my $info (keys %{$extractcontiginfo{$locus}}) {
			#print "locus is $info\n";
			$start = $extractcontiginfo{$locus}{$info} if ($info eq 'sstart');
			$end = $extractcontiginfo{$locus}{$info} if ($info eq 'send');
			$strand = $extractcontiginfo{$locus}{$info} if ($info eq 'strand');
			#print STDERR " strand is $strand\n"; 
			$seq = $extractcontiginfo{$locus}{$info} if ($info eq 'seq');
			$qTE = $extractcontiginfo{$locus}{$info} if ($info eq 'qid');
			$contlen = $extractcontiginfo{$locus}{$info} if ($info eq 'slen');
			
		}
		$newlocus = $locus.".".$start.".".$end.".".$strand if ($strand eq 'plus');
		$newlocus = $locus.".".$end.".".$start.".".$strand if ($strand eq 'minus');
		my $seq_obj = Bio::Seq->new( -display_id => $newlocus, -seq => $seq); # create object with seq
		my $lengthseq = $seq_obj->length;#length of the contig
		#print "the length of the $lengthseq\n";
		if ($strand eq "plus") {# here I am extracting flanking sequence of the TE
			my $leftstart = $start - 40;
			my $leftend = $start + 4;
			my $rightend = $end + 40;
			my $rightstart = $end - 4;
			#print "for a plus $leftstart\t$leftend\t$rightstart\t$rightend\n";
			$leftstart = 1 if ($leftstart <= 0);
			$rightend = $contlen if ($rightend > $contlen);
			$leftflankseq = $seq_obj->subseq($leftstart,$leftend) ;
			$rightflankseq = $seq_obj->subseq($rightstart,$rightend);
		} elsif ($strand eq "minus") {
			my $leftstart = $end - 40;
			my $leftend = $end + 4;
			my $rightend = $start + 40;
			my $rightstart = $start - 4;
			$rightend = $contlen if ($rightend >= $contlen);
			#print "for a minus $leftstart\t$leftend\t$rightstart\t$rightend\n";
			#print STDERR"contig length is $contlen \n";
			$leftstart = 1 if ($leftstart <= 0);
			my $leflankseq = $seq_obj->subseq($leftstart,$leftend);
			my $riflankseq = $seq_obj->subseq($rightstart,$rightend);
			$rightflankseq = &revcom_seq($leflankseq);
			$leftflankseq = &revcom_seq($riflankseq);
		}
		#print "$locus left seq is $leftflankseq\n";
		#print "$locus right seq is $rightflankseq\n";
		#identifies TSD sequence
		$leftkmer = &generate_kmer_pos($leftflankseq,$left,$newlocus) if defined ($leftflankseq);
		$rightkmer = &generate_kmer_pos($rightflankseq,$right,$newlocus) if defined ($rightflankseq);
		#print Dumper %$leftkmer;
		#print Dumper %$rightkmer;
		&compare_kmers($leftkmer,$rightkmer);
		my $tsdinfo = $matchingtsdlist[0];
		@matchingtsdlist = (); 
		#print "tsd for $newlocus is $tsdinfo\n";
		if (defined ($tsdinfo)) {
			($testart,$teend,$tsdseq,$tsdsucc) = &findcord_bothtsds($tsdinfo);
			#print "$locus output from findcordbothtsd is $testart, $teend\n";
			$teend = $lengthseq if  ($teend > $lengthseq);
			$teseq = $seq_obj->subseq($testart,$teend);
			$teseq = &revcom_seq($teseq) if ($strand eq 'minus');
			#my $lenseq = ($teend-$testart) +1;
			my $lenseq = length($teseq);
			if ($tsdsucc eq "yes") {
			
				my $newID = join (".", $locus,$tsdseq,$qTE,$lenseq);
				&checkprediction($genomeloc,$tsdsucc,$newID);
				my $seqout_obj = Bio::Seq->new( -display_id => $newID, -seq => $teseq);
				push (@full_lengTEs,$seqout_obj);
			
			} elsif ($tsdsucc eq "no") {
			
				my $newid =  join (".", $locus,$qTE,"noTSD",$lenseq);
				&checkprediction($genomeloc,$tsdsucc,$newid);
				my $seqout_obj = Bio::Seq->new( -display_id => $newid, -seq => $teseq);
				push (@partial_teseq,$seqout_obj);
			}
		} else {
			$teseq = $seq_obj->subseq($start,$end) if ($strand eq 'plus');
			$teseq = $seq_obj->subseq($end,$start) if ($strand eq 'minus');
			$teseq = &revcom_seq($teseq) if ($strand eq 'minus');
			my $lenseq = length($teseq);
			#my $lenseq = ($start - $end) +1 if ($strand eq 'minus') ;
			#$lenseq = ($end - $start) +1 if ($strand eq 'plus');
			my $new_id =  join (".", $locus,$qTE,"noTSD",$lenseq);
			my $seqout_obj = Bio::Seq->new( -display_id => $new_id, -seq => $teseq);
			push (@partial_teseq,$seqout_obj);
			&checkprediction($genomeloc,"no",$new_id);
		}
	}


	#compare two hashes of alldirectory strand prediction
	#print the strandprediction hash
	@discsplitfullTE = @full_lengTEs;
	@discsplitpartTE = @partial_teseq;
	&printarray_fasta ("$path/$teout.full_len.fasta","displ",@full_lengTEs);
	&printarray_fasta ("$path/$teout.partial.fasta","displ",@partial_teseq);
	&find_no_te_cand();
	#print Dumper %strandprediction;
	%splitdiscpredictionout = %strandprediction;
	&print_keyvalue_HOH("$path/$fineout","displ",\%strandprediction);
	
	print "------------DONE with discordant split only assembly---------------------\n";
	
	if (($auto eq "yes") && ($listfile) && ($discdir)) {
		print "------------extracting candidates for with discordant read only assembly---------------------\n";
		
		open (my $gh, "<", "$path/$fineout.displ") or confess "\n ERROR (main): could not open to read $path/$fineout.displ $!\n";
			while (my $data = <$gh>) { #reading the table
				chomp $data; #storing the table values to the $data and removing the extraline
				my @name = split(/\s+/,$data); # splitting the data based on tab and storing into the arrray
				my $dirc = $name[0];#first column is genome location
				my $suc = $name[1];
				push @unassdispl,$dirc if ($suc eq "no");
				
			}
		close ($gh);
		print Dumper @unassdispl,"unassembled ones from Discosplitassembly\n";
	}

}


if ($discdir) {
	print "################################Discordantreadsmatesassembly###############################\n";
	my $gedcdir;
	my $dpath;
	my @discor;
	@discor = `ls $discdir` if ($folder);
	@discor = &readlocifromfile($listfile) if (($listfile) && ($auto eq "no"));
	@discor = @unassdispl if (($listfile) && ($auto eq "yes"));
	%strandprediction=();
	%extractcontiginfo=();
	%alldirectory=();

	foreach $gedcdir (@discor) {
		#my $present = 0;
		my $success = 0;
		chomp ($gedcdir);
		$alldirectory{$gedcdir} =1;
		
		unless (-e "$path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.fasta") {
			make_path ("$path/Assembled_DCreads/$gedcdir");
			my @disremafiles = glob ("'$discdir/$gedcdir/*disreadsmates.fasta'");
			my $discout = "$path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.fasta";
			#unlink glob "$path/Assembled_DCreads/$gedcdir/$gedcdir.disreads.fasta_k45_ma3.*" if (-e "$path/Assembled_DCreads/$gedcdir/$gedcdir.disreads.fasta_k45_ma3.contigs.fa");
			#unlink glob "$path/Assembled_DCreads/$gedcdir/Renamed_Assembledseq/*.fasta.*" if (-e "$path/Assembled_DCreads/$gedcdir/Renamed_Assembledseq/$gedcdir.rename.fasta");
			#unlink("$path/Assembled_DCreads/$gedcdir/$gedcdir.extract.seq.fa.tabular.blast.out") if (-e "$path/Assembled_DCreads/$gedcdir/$gedcdir.extract.seq.fa.tabular.blast.out");
			&concatenatefiles($discout,@disremafiles);
		}	
		$dpath = "$path/Assembled_DCreads/$gedcdir";
		$dreadmatefile = "$gedcdir.concate.disreadsmates.fasta";
		$DUpfile = "$path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.Up.fasta";
		$DmatefileR1 = "$path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.R1.fasta";
		$DmatefileR2 = "$path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.R2.fasta"; 
		#my $readcount = &splitdismatetopairs($dpath,$dreadmatefile);
		my $readcount = &splitsoftdismatetopairs($dpath,$dreadmatefile);
		my $queryte = &find_TEseq($gedcdir);
		print STDERR "the query file is $queryte.fasta for the locus $gedcdir\n";
		if (($readcount > 1) && ($readcount < 100)) {
			#unless ($present ==1 ){
			system ("$CAP3dir/cap3 $path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.fasta > $path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.asmbl.fasta") == 0 or warn ("unable to assemble fasta $directory \n");
			#}	
			if (-e "$path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.fasta.cap.contigs") {
				my $assembledfile = "$gedcdir.concate.disreadsmates.fasta.cap.contigs";
				
				&renameseq_filename($assembledfile,$dpath,$gedcdir,"cap");
			} else {
				print "No assembly detected for $gedcdir\n"
			}
			
			
			&blast_file ($TEdir,"$queryte.fasta",$gedcdir,$dpath,"cap");	#blast with the query TEsequence
			#print STDERR "the query is $TEdir/$queryte.fasta\n";
			#print STDERR "the disdirectory is $gedcdir\n";
			my $blastfile = "$dpath/$gedcdir.$queryte.fasta.cap.tabular.blast.out";
			my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","cap");
			$qst = 0 if (not defined $qst);
			$qen = 0 if (not defined $qen);
			#print Dumper %$TEblast;
			#parse the blastoutput for max score
			&blast_file ($gdir,"$gedcdir.extract.seq.fa",$gedcdir,$dpath,"cap");
			my $blastgenomefile = "$dpath/$gedcdir.extract.seq.fa.cap.tabular.blast.out";
			my $genomeblast = &parse_blast($blastgenomefile,"95","cap");
			#print "the genome blast results for disc\n";
			#print Dumper %$genomeblast;
			#predict the orientation of TE insertion
			&compare_twoHOH($gedcdir,$TEblast,$genomeblast);
			
	
		} elsif ($readcount >= 100) {
			#unless ($present ==1 ){
				#system ("$spadedir/spades.py -1 $DmatefileR1 -2 $DmatefileR2 -s $DUpfile --careful --only-assembler -t $cpus -o $path/Assembled_DCreads/$gedcdir/$gedcdir.disreadsSPAdeout") == 0 or 
				#system ("$spadedir/dipspades.py -1 $DmatefileR1 -2 $DmatefileR2 -s $DUpfile -t $cpus --only-assembler --expect-rearrangements -o $path/Assembled_DCreads/$gedcdir/$gedcdir.disreadsdipSPAdeout") == 0 or 
				#system ("$spadedir/spades.py --s1 $path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.fasta --careful --only-assembler -t $cpus -o $path/Assembled_DCreads/$gedcdir/$gedcdir.disreadssingleSPAdeout") == 0 or 
				#system ("$miniadir/minia -in $path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_DCreads/$gedcdir/$gedcdir.disreads.fasta_k45_ma3") == 0 or warn ("unable to assemble disfasta $gedcdir \n");
			#}
			unless (-e "$path/Assembled_DCreads/$gedcdir/$gedcdir.disreadsSPAdeout/scaffolds.fasta"){
				system ("$spadedir/spades.py -1 $DmatefileR1 -2 $DmatefileR2 -s $DUpfile --careful --only-assembler -t $cpus -o $path/Assembled_DCreads/$gedcdir/$gedcdir.disreadsSPAdeout") == 0 or warn ("unable to assemble disfasta $gedcdir using spades \n");
			}
			if (-e "$path/Assembled_DCreads/$gedcdir/$gedcdir.disreadsSPAdeout/scaffolds.fasta") {
				my %longestassem;
				#Rename scaffolds.fasta 
				#print "checking of spadesscaffpld";
				copy("$path/Assembled_DCreads/$gedcdir/$gedcdir.disreadsSPAdeout/scaffolds.fasta", "$path/Assembled_DCreads/$gedcdir/$gedcdir.disreads.scaffolds.fasta") or die "Copy failed scaffolds.fasta $gedcdir:$!";
				#rename the fasta sequence with its filename
				my $assembledfile = "$gedcdir.disreads.scaffolds.fasta";
				&renameseq_filename($assembledfile,$dpath,$gedcdir,"SPAscaf");
				&blast_file ($TEdir,"$queryte.fasta",$gedcdir,$dpath,"SPAscaf");
				my $blastfile = "$dpath/$gedcdir.$queryte.fasta.SPAscaf.tabular.blast.out";
				my $blastfilesize = -s $blastfile;
				my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","SPAscaf");	#blast with the query TEsequence
				$qst = 0 if (not defined $qst);
				$qen = 0 if (not defined $qen);
				print "$qst,$qen for $gedcdir SPAscaf\n";
				my $TEblastSPAscaf = $TEblast;
				print Dumper $TEblastSPAscaf,"\n";
				$longestassem{'SPAscaf'} = ($qen - $qst) ;
				if (($blastfilesize == 0) ||((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280))|| (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))   ) {
					system ("$spadedir/spades.py --s1 $path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.fasta --careful --only-assembler -t $cpus -o $path/Assembled_DCreads/$gedcdir/$gedcdir.disreadssingleSPAdeout") == 0 or warn ("unable to assemble disfasta $gedcdir \n");
					#system ("$miniadir/minia -in $path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_DCreads/$gedcdir/$gedcdir.disreads.fasta_k45_ma3") == 0 or 
					if (-e "$path/Assembled_DCreads/$gedcdir/$gedcdir.disreadssingleSPAdeout/scaffolds.fasta") {
						copy("$path/Assembled_DCreads/$gedcdir/$gedcdir.disreadssingleSPAdeout/scaffolds.fasta",	"$path/Assembled_DCreads/$gedcdir/$gedcdir.disreads.single_scaffolds.fasta") or die "Copy failed single_scaffolds.fasta $gedcdir $!";
						my $assembledfile = "$gedcdir.disreads.single_scaffolds.fasta";
						&renameseq_filename($assembledfile,$dpath,$gedcdir,"singSPA");	
						&blast_file ($TEdir,"$queryte.fasta",$gedcdir,$dpath,"singSPA");
						my $blastfile = "$dpath/$gedcdir.$queryte.fasta.singSPA.tabular.blast.out";#$qloc.$type.tabular.blast.out
						$qst = 0 if (not defined $qst);
						$qen = 0 if (not defined $qen);
						my $blastfilesize = -s $blastfile;
						my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","singSPA") if ($blastfilesize > 0);
						print "$qst,$qen for $gedcdir singSPA\n";
						my $TEblastsingSPA = $TEblast;
						print Dumper $TEblastsingSPA,"\n";
						$longestassem{'singSPA'} = ($qen - $qst) ;
						if (($blastfilesize == 0) || ((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280))  ) {
							system ("$miniadir/minia -in $path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_DCreads/$gedcdir/$gedcdir.disreads.fasta_k45_ma3") == 0 or warn ("unable to assemble disfasta $gedcdir \n");
							if (-e "$path/Assembled_DCreads/$gedcdir/$gedcdir.disreads.fasta_k45_ma3.contigs.fa") { 
								my $assembledfile = "$gedcdir.disreads.fasta_k45_ma3.contigs.fa";
								&renameseq_filename($assembledfile,$dpath,$gedcdir,"mini");
								&blast_file ($TEdir,"$queryte.fasta",$gedcdir,$dpath,"mini");
								my $blastfile = "$dpath/$gedcdir.$queryte.fasta.mini.tabular.blast.out";#$qloc.$type.tabular.blast.out
								my $blastfilesize = -s $blastfile;
								my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","mini");
								$qst = 0 if (not defined $qst);
								$qen = 0 if (not defined $qen);
								print "$qst,$qen for $gedcdir mini\n";
								my $TEblastmini = $TEblast;
								print Dumper $TEblastmini,"\n";
								$longestassem{'mini'} = ($qen - $qst) ;
								if (($blastfilesize == 0) || ((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280))|| (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))   ) {
									my $longest = &findmaxvaluehash(%longestassem);
									&blast_file ($gdir,"$gedcdir.extract.seq.fa",$gedcdir,$dpath,$longest);
									my $blastgenomefile = "$dpath/$gedcdir.extract.seq.fa.$longest.tabular.blast.out";
									my $genomeblast = &parse_blast($blastgenomefile,"95",$longest);
									#print "discsplit genomeblast\n";
									#print Dumper %$genomeblast;
									#predict the orientation of TE insertion
									&compare_twoHOH($gedcdir,$TEblastmini,$genomeblast) if ($longest eq "mini");
									&compare_twoHOH($gedcdir,$TEblastsingSPA,$genomeblast) if ($longest eq "singSPA");
									#&compare_twoHOH($gedcdir,$TEblastdipSPA,$genomeblast) if ($longest eq "dipSPA");
									&compare_twoHOH($gedcdir,$TEblastSPAscaf,$genomeblast) if ($longest eq "SPAscaf");
									$success = 1;
								} else {
									&blast_file ($gdir,"$gedcdir.extract.seq.fa",$gedcdir,$dpath,"mini");
									my $blastgenomefile = "$dpath/$gedcdir.extract.seq.fa.mini.tabular.blast.out";
									my $genomeblast = &parse_blast($blastgenomefile,"95","mini");
									#print "discsplit genomeblast\n";
									#print Dumper %$genomeblast;
									#predict the orientation of TE insertion
									&compare_twoHOH($gedcdir,$TEblastmini,$genomeblast);
									$success = 1;
								}
							}
						}else {
							&blast_file ($gdir,"$gedcdir.extract.seq.fa",$gedcdir,$dpath,"singSPA");
							my $blastgenomefile = "$dpath/$gedcdir.extract.seq.fa.singSPA.tabular.blast.out";
							my $genomeblast = &parse_blast($blastgenomefile,"95","singSPA");
							#print "discsplit genomeblast\n";
							#print Dumper %$genomeblast;
							#predict the orientation of TE insertion
							&compare_twoHOH($gedcdir,$TEblastsingSPA,$genomeblast);
							$success = 1;
						}
						
					}
				} else {
					#parse the blastoutput for max score
					&blast_file ($gdir,"$gedcdir.extract.seq.fa",$gedcdir,$dpath,"SPAscaf");
					my $blastgenomefile = "$dpath/$gedcdir.extract.seq.fa.SPAscaf.tabular.blast.out";
					my $genomeblast = &parse_blast($blastgenomefile,"95","SPAscaf");
					#print "discsplit genomeblast\n";
					#print Dumper %$genomeblast;
					#predict the orientation of TE insertion
					&compare_twoHOH($gedcdir,$TEblastSPAscaf,$genomeblast);
					$success = 1;
				}
				
			} 
			if ((! -e "$path/Assembled_DCreads/$gedcdir/$gedcdir.disreadssingleSPAdeout/scaffolds.fasta") &&($success == 0 )){
				my %longestassem;
				system ("$spadedir/spades.py --s1 $path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.fasta --careful --only-assembler -t $cpus -o $path/Assembled_DCreads/$gedcdir/$gedcdir.disreadssingleSPAdeout") == 0 or warn ("unable to assemble singSPA $gedcdir \n");	
				if (-e "$path/Assembled_DCreads/$gedcdir/$gedcdir.disreadssingleSPAdeout/scaffolds.fasta") {
					copy("$path/Assembled_DCreads/$gedcdir/$gedcdir.disreadssingleSPAdeout/scaffolds.fasta",	"$path/Assembled_DCreads/$gedcdir/$gedcdir.disreads.single_scaffolds.fasta") or die "Copy failed single_scaffolds.fasta $gedcdir $!";
					my $assembledfile = "$gedcdir.disreads.single_scaffolds.fasta";
					&renameseq_filename($assembledfile,$dpath,$gedcdir,"singSPA");	
					&blast_file ($TEdir,"$queryte.fasta",$gedcdir,$dpath,"singSPA");
					my $blastfile = "$dpath/$gedcdir.$queryte.fasta.singSPA.tabular.blast.out";#$qloc.$type.tabular.blast.out
					my $blastfilesize = -s $blastfile;
					my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","singSPA") if ($blastfilesize > 0);
					$qst = 0 if (not defined $qst);
					$qen = 0 if (not defined $qen);
					my $TEblastsingSPA = $TEblast;
					$longestassem{'singSPA'} = ($qen - $qst) ;
				
					if (($blastfilesize == 0) || ((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280)) || (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))  ) {
						system ("$miniadir/minia -in $path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_DCreads/$gedcdir/$gedcdir.disreads.fasta_k45_ma3") == 0 or warn ("unable to assemble disfasta $gedcdir \n");	
						if (-e "$path/Assembled_DCreads/$gedcdir/$gedcdir.disreads.fasta_k45_ma3.contigs.fa") { 
							my $assembledfile = "$gedcdir.disreads.fasta_k45_ma3.contigs.fa";
							&renameseq_filename($assembledfile,$dpath,$gedcdir,"mini");
							&blast_file ($TEdir,"$queryte.fasta",$gedcdir,$dpath,"mini");
							my $blastfile = "$dpath/$gedcdir.$queryte.fasta.mini.tabular.blast.out";#$qloc.$type.tabular.blast.out
							my $blastfilesize = -s $blastfile;
							my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","mini");
							$qst = 0 if (not defined $qst);
							$qen = 0 if (not defined $qen);
							my $TEblastmini = $TEblast;
							$longestassem{'mini'} = ($qen - $qst) ;
						
							if (($blastfilesize == 0) || ((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280)) || (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))  ) {
								my $longest = &findmaxvaluehash(%longestassem);
								&blast_file ($gdir,"$gedcdir.extract.seq.fa",$gedcdir,$dpath,$longest);
								my $blastgenomefile = "$dpath/$gedcdir.extract.seq.fa.$longest.tabular.blast.out";
								my $genomeblast = &parse_blast($blastgenomefile,"95",$longest);
								#print "discsplit genomeblast\n";
								#print Dumper %$genomeblast;
								#predict the orientation of TE insertion
								&compare_twoHOH($gedcdir,$TEblastmini,$genomeblast) if ($longest eq "mini");
								&compare_twoHOH($gedcdir,$TEblastsingSPA,$genomeblast) if ($longest eq "singSPA");
								$success = 1;
							} else {
								&blast_file ($gdir,"$gedcdir.extract.seq.fa",$gedcdir,$dpath,"mini");
								my $blastgenomefile = "$dpath/$gedcdir.extract.seq.fa.mini.tabular.blast.out";
								my $genomeblast = &parse_blast($blastgenomefile,"95","mini");
								#print "discsplit genomeblast\n";
								#print Dumper %$genomeblast;
								#predict the orientation of TE insertion
								&compare_twoHOH($gedcdir,$TEblastmini,$genomeblast);
								$success = 1;
							}
						}
					} else {
						&blast_file ($gdir,"$gedcdir.extract.seq.fa",$gedcdir,$dpath,"singSPA");
						my $blastgenomefile = "$dpath/$gedcdir.extract.seq.fa.singSPA.tabular.blast.out";
						my $genomeblast = &parse_blast($blastgenomefile,"95","singSPA");
						#print "discsplit genomeblast\n";
						#print Dumper %$genomeblast;
						#predict the orientation of TE insertion
						&compare_twoHOH($gedcdir,$TEblastsingSPA,$genomeblast);
						$success = 1;
					}
				}
			}
			if ( (! -e "$path/Assembled_DCreads/$gedcdir/$gedcdir.disreads.fasta_k45_ma3.contigs.fa")  && ($success == 0)) {
				system ("$miniadir/minia -in $path/Assembled_DCreads/$gedcdir/$gedcdir.concate.disreadsmates.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_DCreads/$gedcdir/$gedcdir.disreads.fasta_k45_ma3") == 0 or warn ("unable to assemble disfasta $gedcdir \n");	
				my %longestassem;
				if (-e "$path/Assembled_DCreads/$gedcdir/$gedcdir.disreads.fasta_k45_ma3.contigs.fa") { 
					my $assembledfile = "$gedcdir.disreads.fasta_k45_ma3.contigs.fa";
					&renameseq_filename($assembledfile,$dpath,$gedcdir,"mini");
					&blast_file ($TEdir,"$queryte.fasta",$gedcdir,$dpath,"mini");
					my $blastfile = "$dpath/$gedcdir.$queryte.fasta.mini.tabular.blast.out";#$qloc.$type.tabular.blast.out
					my $blastfilesize = -s $blastfile;
					my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","mini") ;
					my $TEblastmini = $TEblast;
					$longestassem{'mini'} = ($qen - $qst) ;
				
					if (($blastfilesize == 0) || ((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280)) || (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))  ) {
						my $longest = &findmaxvaluehash(%longestassem);
						&blast_file ($gdir,"$gedcdir.extract.seq.fa",$gedcdir,$dpath,$longest);
						my $blastgenomefile = "$dpath/$gedcdir.extract.seq.fa.$longest.tabular.blast.out";
						my $genomeblast = &parse_blast($blastgenomefile,"95",$longest);
						#print "discsplit genomeblast\n";
						#print Dumper %$genomeblast;
						#predict the orientation of TE insertion
						&compare_twoHOH($gedcdir,$TEblastmini,$genomeblast) if ($longest eq "mini");
					} else {
						&blast_file ($gdir,"$gedcdir.extract.seq.fa",$gedcdir,$dpath,"mini");
						my $blastgenomefile = "$dpath/$gedcdir.extract.seq.fa.mini.tabular.blast.out";
						my $genomeblast = &parse_blast($blastgenomefile,"95","mini");
						#print "discsplit genomeblast\n";
						#print Dumper %$genomeblast;
						#predict the orientation of TE insertion
						&compare_twoHOH($gedcdir,$TEblastmini,$genomeblast);
					}
				}
			}
			

		} elsif ($readcount==0) {
			push (@zerofiles,"$gedcdir");
		}
	}

	#print Dumper %extractcontiginfo;
	#print Dumper %alldirectory;

	print  "DONE with prediction of orientation.... \n";

	print  "\n\nstarting to extract TEsequences.... \n";

	$left = 5;
	$right = 3;
	@matchingtsdlist =();
	@full_lengTEs=();
	@partial_teseq=();
	#my $teseq;
	&loadseq_tohash("Assembled_DCreads");
	#print Dumper %extractcontiginfo;
	#find the contigs with TEsequence and flanking sequence and extract the flanking sequence
	foreach my $locus (sort keys %extractcontiginfo) {
		#print "locus is $locus\n";
		my ($typeasm,$genomeloc,$rest) = split (/\./,$locus,3);
		my ($leftkmer);
		my ($rightkmer);
		my $leftflankseq;
		my $rightflankseq;
		my $start;
		my $end;
		my $strand;
		my $seq;
		my $qTE;
		my $newlocus;
		my $teseq;
		my $testart;
		my $teend;
		my $tsdseq;
		my $tsdsucc;
		my $contlen;
		foreach my $info (keys %{$extractcontiginfo{$locus}}) {
			#print "locus is $info\n";
			$start = $extractcontiginfo{$locus}{$info} if ($info eq 'sstart');
			$end = $extractcontiginfo{$locus}{$info} if ($info eq 'send');
			$strand = $extractcontiginfo{$locus}{$info} if ($info eq 'strand');
			#print STDERR " strand is $strand\n"; 
			$seq = $extractcontiginfo{$locus}{$info} if ($info eq 'seq');
			$qTE = $extractcontiginfo{$locus}{$info} if ($info eq 'qid');
			$contlen = $extractcontiginfo{$locus}{$info} if ($info eq 'slen');
		
		}
		$newlocus = $locus.".".$start.".".$end.".".$strand if ($strand eq 'plus');
		$newlocus = $locus.".".$end.".".$start.".".$strand if ($strand eq 'minus');
		my $seq_obj = Bio::Seq->new( -display_id => $newlocus, -seq => $seq); # create object with seq
		my $lengthseq = $seq_obj->length;#length of the contig
		if ($strand eq "plus") {# here I am extracting flanking sequence of the TE
			my $leftstart = $start - 40;
			my $leftend = $start + 4;
			my $rightend = $end + 40;
			my $rightstart = $end - 4;
			#print "for a plus $leftstart\t$leftend\t$rightstart\t$rightend\n";
			$leftstart = 1 if ($leftstart <= 0);
			$rightend = $contlen if ($rightend > $contlen);
			$leftflankseq = $seq_obj->subseq($leftstart,$leftend) ;
			$rightflankseq = $seq_obj->subseq($rightstart,$rightend);
		} elsif ($strand eq "minus") {
			my $leftstart = $end - 40;
			my $leftend = $end + 4;
			my $rightend = $start + 40;
			my $rightstart = $start - 4;
		
			#print "for a minus $leftstart\t$leftend\t$rightstart\t$rightend\n";
			$leftstart = 1 if ($leftstart <= 0);
			$rightend = $contlen if ($rightend > $contlen);
			my $leflankseq = $seq_obj->subseq($leftstart,$leftend);
			my $riflankseq = $seq_obj->subseq($rightstart,$rightend);
			$rightflankseq = &revcom_seq($leflankseq);
			$leftflankseq = &revcom_seq($riflankseq);
		}
		#print "$locus left seq is $leftflankseq\n";
		#print "$locus right seq is $rightflankseq\n";
		#identifies TSD sequence
		$leftkmer = &generate_kmer_pos($leftflankseq,$left,$newlocus) if defined ($leftflankseq);
		$rightkmer = &generate_kmer_pos($rightflankseq,$right,$newlocus) if defined ($rightflankseq);
		#print Dumper %$leftkmer;
		#print Dumper %$rightkmer;
		&compare_kmers($leftkmer,$rightkmer);
		my $tsdinfo = $matchingtsdlist[0];
		@matchingtsdlist = (); 
		#print "tsd for $newlocus is $tsdinfo\n";
		if (defined ($tsdinfo)) {
			($testart,$teend,$tsdseq,$tsdsucc) = &findcord_bothtsds($tsdinfo);
			#print "$locus output from findcordbothtsd is $testart, $teend\n";
			$teend = $lengthseq if ($teend > $lengthseq);
			$teseq = $seq_obj->subseq($testart,$teend);
			$teseq = &revcom_seq($teseq) if ($strand eq 'minus');
			my $lenseq = length($teseq);
			if ($tsdsucc eq "yes") {
			
				my $newID = join (".", $locus,$tsdseq,$qTE,$lenseq);
				&checkprediction($genomeloc,$tsdsucc,$newID);
				my $seqout_obj = Bio::Seq->new( -display_id => $newID, -seq => $teseq);
				push (@full_lengTEs,$seqout_obj);
			
			} elsif ($tsdsucc eq "no") {
			
				my $newid =  join (".", $locus,$qTE,"noTSD",$lenseq);
				&checkprediction($genomeloc,$tsdsucc,$newid);
				my $seqout_obj = Bio::Seq->new( -display_id => $newid, -seq => $teseq);
				push (@partial_teseq,$seqout_obj);
			}
		} else {
			$teseq = $seq_obj->subseq($start,$end) if ($strand eq 'plus');
			$teseq = $seq_obj->subseq($end,$start) if ($strand eq 'minus');
			$teseq = &revcom_seq($teseq) if ($strand eq 'minus');
			my $lenseq = length($teseq);
			my $new_id =  join (".", $locus,$qTE,"noTSD",$lenseq);
			my $seqout_obj = Bio::Seq->new( -display_id => $new_id, -seq => $teseq);
			push (@partial_teseq,$seqout_obj);
			&checkprediction($genomeloc,"no",$new_id);
		}
	}


	#compare two hashes of alldirectory strand prediction
	#print the strandprediction hash
	
	@discfullTE = @full_lengTEs;
	@discpartTE = @partial_teseq;
	&printarray_fasta ("$path/$teout.full_len.fasta","disc",@full_lengTEs);
	&printarray_fasta ("$path/$teout.partial.fasta","disc",@partial_teseq);
	&find_no_te_cand();
	#print Dumper %strandprediction;
	%discpredictionout = %strandprediction;
	&print_keyvalue_HOH("$path/$fineout","disc",\%strandprediction);
	
	
	print "------------DONE with discordant only assembly---------------------\n";
	if (($auto eq "yes") && ($listfile) && ($rdir)) {
		print "------------extracting candidates for with discordant read only assembly---------------------\n";
		
		open (my $dh, "<", "$path/$fineout.disc") or confess "\n ERROR (main): could not open to read $path/$fineout.disc $!\n";
			while (my $data = <$dh>) { #reading the table
				chomp $data; #storing the table values to the $data and removing the extraline
				my @name = split(/\s+/,$data); # splitting the data based on tab and storing into the arrray
				my $dirc = $name[0];#first column is genome location
				my $suc = $name[1];
				push @unassdisc,$dirc if ($suc eq "no");
			}
		close ($dh);
		print Dumper @unassdisc,"unassembled ones from Discassembly\n";
	}
	
}



if ($rdir) {
	print "###############################allreadsassembly###############################\n";
	%strandprediction=();
	%extractcontiginfo=();
	%alldirectory=();
	my @dir;
	@dir =`ls $rdir`;
	@dir = &readlocifromfile($listfile) if ($listfile);
	@dir = @unassdisc if (($listfile) && ($auto));
	
	foreach $directory (sort {$a cmp $b} @dir) {
		#my $present = 0;
		my $success = 0;
		chomp ($directory);
		$alldirectory{$directory} =1;
		#my %longestassem;
		print STDERR "the directory now analysing is $directory\n";
		$outpath = "$path/Assembled_TEreads/$directory";
		my $R1out = "$path/Assembled_TEreads/$directory/$directory.concatenated.R1.fasta";
		my $R2out = "$path/Assembled_TEreads/$directory/$directory.concatenated.R2.fasta";
		my $Upout = "$path/Assembled_TEreads/$directory/$directory.concatenated.Up.fasta";
		unless (-e "$path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta") {
			make_path ("$path/Assembled_TEreads/$directory") ;
			my @R1files = glob ("'$rdir/$directory/*R1.fasta'");
			my @R2files = glob ("'$rdir/$directory/*R2.fasta'");
			my @UPfiles = glob ("'$rdir/$directory/*Up.fasta'");
			&concatenatefiles($R1out,@R1files);
			&concatenatefiles($R2out,@R2files);
			&concatenatefiles($Upout,@UPfiles);
			system ("cat $path/Assembled_TEreads/$directory/$directory.concatenated.R1.fasta $path/Assembled_TEreads/$directory/$directory.concatenated.R2.fasta $path/Assembled_TEreads/$directory/$directory.concatenated.Up.fasta > $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta") == 0 or die ("unable to concatenate discordantsplitreads $directory \n");
		}
		
		my $R1filesize = (-s "$path/Assembled_TEreads/$directory/$directory.concatenated.R1.fasta");
		my $R2filesize = (-s "$path/Assembled_TEreads/$directory/$directory.concatenated.R2.fasta");
		
		my $queryte = &find_TEseq($directory);
		print STDERR "the query file is $queryte.fasta for the locus $directory\n";
		#assemble the concatenated sequence
	
		# print STDERR "assemble the concatenated sequence\n";
		if (($R1filesize > 0) || ($R1filesize > 0)) {
			#unless ($present == 1) {
				#Assemble only the mapped reads		   
				#system ("$spadedir/spades.py -1 $R1out -2 $R2out -s $Upout --careful --only-assembler -t $cpus -o $path/Assembled_TEreads/$directory/$directory.allreadsSPAdeout") == 0 or 
				#system ("$spadedir/dipspades.py -1 $R1out -2 $R2out -s $Upout -t $cpus --only-assembler --expect-rearrangements -o $path/Assembled_TEreads/$directory/$directory.allreadsdipSPAdeout") == 0 or 
				#system ("$spadedir/spades.py --s1 $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta -t $cpus --only-assembler --careful $path/Assembled_TEreads/$directory/$directory.allreadsingleSPAdeout") == 0 or 
				#system ("$miniadir/minia -in $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta_k45_ma3") == 0 or warn ("unable to assemble fasta $directory \n");
				#system ("$CAP3dir/cap3 $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta > $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.asmbl.fasta") == 0 or warn ("unable to assemble fasta $directory \n");
			#}
			unless (-e "$path/Assembled_TEreads/$directory/$directory.allreadsSPAdeout/scaffolds.fasta") {
				system ("$spadedir/spades.py -1 $R1out -2 $R2out -s $Upout --careful --only-assembler -t $cpus -o $path/Assembled_TEreads/$directory/$directory.allreadsSPAdeout") == 0 or warn ("unable to assemble fasta with spades $directory \n");
			}
			if (-e "$path/Assembled_TEreads/$directory/$directory.allreadsSPAdeout/scaffolds.fasta") {
				my %longestassem;
				#Rename scaffolds.fasta 
				copy("$path/Assembled_TEreads/$directory/$directory.allreadsSPAdeout/scaffolds.fasta", "$path/Assembled_TEreads/$directory/$directory.allreads.scaffolds.fasta") or die "Copy failed scaffolds.fasta $directory:$!";
				#rename the fasta sequence with its filename
				my $assembledfile = "$directory.allreads.scaffolds.fasta";
				&renameseq_filename($assembledfile,$outpath,$directory,"SPAscaf");
				&blast_file ($TEdir,"$queryte.fasta",$directory,$outpath,"SPAscaf");	#blast with the query TEsequence
				my $blastfile = "$outpath/$directory.$queryte.fasta.SPAscaf.tabular.blast.out";#$qloc.$type.tabular.blast.out
				my $blastfilesize = -s $blastfile;
				my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","SPAscaf") if ($blastfilesize > 0);
				print "$qst,$qen for $directory SPAscaf\n";
				my $TEblastSPAscaf = $TEblast;
				print Dumper $TEblastSPAscaf,"\n";
				$longestassem{'SPAscaf'} = ($qen - $qst) if ($blastfilesize > 0);
				if (($blastfilesize == 0) ||((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280))|| (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))   ) {
					system ("$spadedir/spades.py --s1 $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta --careful --only-assembler -t $cpus -o $path/Assembled_TEreads/$directory/$directory.allreadsingleSPAdeout") == 0 or warn ("unable to assemble disfasta $directory \n");
					#system ("$miniadir/minia -in $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta_k45_ma3") == 0 or 
					if (-e "$path/Assembled_TEreads/$directory/$directory.allreadsingleSPAdeout/scaffolds.fasta") {
						copy("$path/Assembled_TEreads/$directory/$directory.allreadsingleSPAdeout/scaffolds.fasta",	"$path/Assembled_TEreads/$directory/$directory.allreads.single_scaffolds.fasta") or die "Copy failed single_scaffolds.fasta $directory $!";
						my $assembledfile = "$directory.allreads.single_scaffolds.fasta";
						&renameseq_filename($assembledfile,$outpath,$directory,"singSPA");	
						&blast_file ($TEdir,"$queryte.fasta",$directory,$outpath,"singSPA");
						my $blastfile = "$outpath/$directory.$queryte.fasta.singSPA.tabular.blast.out";#$qloc.$type.tabular.blast.out
						my $blastfilesize = -s $blastfile;
						my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","singSPA") if ($blastfilesize > 0);
						$qst = 0 if (not defined $qst);
						$qen = 0 if (not defined $qen);
						print "$qst,$qen for $directory singSPA\n";
						my $TEblastsingSPA = $TEblast;
						print Dumper $TEblastsingSPA,"\n";
						$longestassem{'singSPA'} = ($qen - $qst) ;
						#$longestassem{$directory}{'singSPA'}{'qend'} = $qen if ($blastfilesize > 0);
						if (($blastfilesize == 0) || ((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280)) || (($blastfilesize > 0) && (not defined $qst) && (not defined $qen)) ) {
							system ("$miniadir/minia -in $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta_k45_ma3") == 0 or warn ("unable to assemble disfasta $directory \n");	
							if (-e "$path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta_k45_ma3.contigs.fa") { 
								my $assembledfile = "$directory.concatenated.allreads.fasta_k45_ma3.contigs.fa";
								&renameseq_filename($assembledfile,$outpath,$directory,"mini");
								&blast_file ($TEdir,"$queryte.fasta",$directory,$outpath,"mini");
								my $blastfile = "$outpath/$directory.$queryte.fasta.mini.tabular.blast.out";#$qloc.$type.tabular.blast.out
								my $blastfilesize = -s $blastfile;
								my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","mini");
								$qst = 0 if (not defined $qst);
								$qen = 0 if (not defined $qen);
								print "$qst,$qen for $directory mini\n";
								my $TEblastmini = $TEblast;
								print Dumper $TEblastmini,"\n";
								$longestassem{'mini'} = ($qen - $qst) ;
								
								if (($blastfilesize == 0) || ((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280))|| (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))   ) {
									my $longest = &findmaxvaluehash(%longestassem);
									&blast_file ($gdir,"$directory.extract.seq.fa",$directory,$outpath,$longest);
									my $blastgenomefile = "$outpath/$directory.extract.seq.fa.$longest.tabular.blast.out";
									my $genomeblast = &parse_blast($blastgenomefile,"95",$longest);
									#print "discsplit genomeblast\n";
									#print Dumper %$genomeblast;
									#predict the orientation of TE insertion
									&compare_twoHOH($directory,$TEblastmini,$genomeblast) if ($longest eq "mini");
									&compare_twoHOH($directory,$TEblastsingSPA,$genomeblast) if ($longest eq "singSPA");
									#&compare_twoHOH($directory,$TEblastdipSPA,$genomeblast) if ($longest eq "dipSPA");
									&compare_twoHOH($directory,$TEblastSPAscaf,$genomeblast) if ($longest eq "SPAscaf");
									$success = 1;
								} else {
									&blast_file ($gdir,"$directory.extract.seq.fa",$directory,$outpath,"mini");
									my $blastgenomefile = "$outpath/$directory.extract.seq.fa.mini.tabular.blast.out";
									my $genomeblast = &parse_blast($blastgenomefile,"95","mini");
									#print "discsplit genomeblast\n";
									#print Dumper %$genomeblast;
									#predict the orientation of TE insertion
									&compare_twoHOH($directory,$TEblastmini,$genomeblast);
									$success = 1;
								}
							}
						} else {
							&blast_file ($gdir,"$directory.extract.seq.fa",$directory,$outpath,"singSPA");
							my $blastgenomefile = "$outpath/$directory.extract.seq.fa.singSPA.tabular.blast.out";
							my $genomeblast = &parse_blast($blastgenomefile,"95","singSPA");
							#print "discsplit genomeblast\n";
							#print Dumper %$genomeblast;
							#predict the orientation of TE insertion
							&compare_twoHOH($directory,$TEblastsingSPA,$genomeblast);
							$success = 1;
						}
					}	
				} else {
					#parse the blastoutput for max score
					&blast_file ($gdir,"$directory.extract.seq.fa",$directory,$outpath,"SPAscaf");
					my $blastgenomefile = "$outpath/$directory.extract.seq.fa.SPAscaf.tabular.blast.out";
					my $genomeblast = &parse_blast($blastgenomefile,"95","SPAscaf");
					#print "discsplit genomeblast\n";
					#print Dumper %$genomeblast;
					#predict the orientation of TE insertion
					&compare_twoHOH($directory,$TEblastSPAscaf,$genomeblast);
					$success = 1;
				}  
			}  
 			if ((! -e "$path/Assembled_TEreads/$directory/$directory.allreadsingleSPAdeout/scaffolds.fasta") && ($success == 0 )) {
				my %longestassem;
				system ("$spadedir/spades.py --s1 $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta --careful --only-assembler -t $cpus -o $path/Assembled_TEreads/$directory/$directory.allreadsingleSPAdeout") == 0 or warn ("unable to assemble singSPA $directory \n");	
 				if (-e "$path/Assembled_TEreads/$directory/$directory.allreadsingleSPAdeout/scaffolds.fasta") {
					copy("$path/Assembled_TEreads/$directory/$directory.allreadsingleSPAdeout/scaffolds.fasta",	"$path/Assembled_TEreads/$directory/$directory.allreads.single_scaffolds.fasta") or die "Copy failed single_scaffolds.fasta $directory $!";
					my $assembledfile = "$directory.allreads.single_scaffolds.fasta";
					&renameseq_filename($assembledfile,$outpath,$directory,"singSPA");	
					&blast_file ($TEdir,"$queryte.fasta",$directory,$outpath,"singSPA");
					my $blastfile = "$outpath/$directory.$queryte.fasta.singSPA.tabular.blast.out";#$qloc.$type.tabular.blast.out
					my $blastfilesize = -s $blastfile;
					my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","singSPA") if ($blastfilesize > 0);
					$qst = 0 if (not defined $qst);
					$qen = 0 if (not defined $qen);
					my $TEblastsingSPA = $TEblast;
					$longestassem{'singSPA'} = ($qen - $qst) ;
				
					if (($blastfilesize == 0) || ((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280)) || (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))  ) {
						system ("$miniadir/minia -in $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta_k45_ma3") == 0 or warn ("unable to assemble disfasta $directory \n");	
						if (-e "$path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta_k45_ma3.contigs.fa") { 
							my $assembledfile = "$directory.concatenated.allreads.fasta_k45_ma3.contigs.fa";
							&renameseq_filename($assembledfile,$outpath,$directory,"mini");
							&blast_file ($TEdir,"$queryte.fasta",$directory,$outpath,"mini");
							my $blastfile = "$outpath/$directory.$queryte.fasta.mini.tabular.blast.out";#$qloc.$type.tabular.blast.out
							my $blastfilesize = -s $blastfile;
							my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","mini");
							$qst = 0 if (not defined $qst);
							$qen = 0 if (not defined $qen);
							my $TEblastmini = $TEblast;
							$longestassem{'mini'} = ($qen - $qst) ;
						
							if (($blastfilesize == 0) || ((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280)) || (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))  ) {
								my $longest = &findmaxvaluehash(%longestassem);
								&blast_file ($gdir,"$directory.extract.seq.fa",$directory,$outpath,$longest);
								my $blastgenomefile = "$outpath/$directory.extract.seq.fa.$longest.tabular.blast.out";
								my $genomeblast = &parse_blast($blastgenomefile,"95",$longest);
								#print "discsplit genomeblast\n";
								#print Dumper %$genomeblast;
								#predict the orientation of TE insertion
								&compare_twoHOH($directory,$TEblastmini,$genomeblast) if ($longest eq "mini");
								&compare_twoHOH($directory,$TEblastsingSPA,$genomeblast) if ($longest eq "singSPA");
								$success=1;
							} else {
								&blast_file ($gdir,"$directory.extract.seq.fa",$directory,$outpath,"mini");
								my $blastgenomefile = "$outpath/$directory.extract.seq.fa.mini.tabular.blast.out";
								my $genomeblast = &parse_blast($blastgenomefile,"95","mini");
								#print "discsplit genomeblast\n";
								#print Dumper %$genomeblast;
								#predict the orientation of TE insertion
								&compare_twoHOH($directory,$TEblastmini,$genomeblast);
								$success=1;
							}
						}
					} else {
						&blast_file ($gdir,"$directory.extract.seq.fa",$directory,$outpath,"singSPA");
						my $blastgenomefile = "$outpath/$directory.extract.seq.fa.singSPA.tabular.blast.out";
						my $genomeblast = &parse_blast($blastgenomefile,"95","singSPA");
						#print "discsplit genomeblast\n";
						#print Dumper %$genomeblast;
						#predict the orientation of TE insertion
						&compare_twoHOH($directory,$TEblastsingSPA,$genomeblast);
						$success=1;
					}
				} 
 			}
 			if ((! -e "$path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta_k45_ma3.contigs.fa") && ($success == 0)){
 				my %longestassem;
 				system ("$miniadir/minia -in $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta -kmer-size 45 -abundance-min 3 -out $path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta_k45_ma3") == 0 or warn ("unable to assemble disfasta $directory \n");	
				if (-e "$path/Assembled_TEreads/$directory/$directory.concatenated.allreads.fasta_k45_ma3.contigs.fa") { 
					my $assembledfile = "$directory.concatenated.allreads.fasta_k45_ma3.contigs.fa";
					&renameseq_filename($assembledfile,$outpath,$directory,"mini");
					&blast_file ($TEdir,"$queryte.fasta",$directory,$outpath,"mini");
					my $blastfile = "$outpath/$directory.$queryte.fasta.mini.tabular.blast.out";#$qloc.$type.tabular.blast.out
					my $blastfilesize = -s $blastfile;
					my ($qst,$qen,$TEblast) = &parse_blast($blastfile,"90","mini") ;
					$qst = 0 if (not defined $qst);
					$qen = 0 if (not defined $qen);
					my $TEblastmini = $TEblast;
					$longestassem{'mini'} = ($qen - $qst) ;
				
					if (($blastfilesize == 0) || ((defined $qst) && ($qst > 4)) || ((defined $qen) && ($qen < 280)) || (($blastfilesize > 0) && (not defined $qst) && (not defined $qen))  ) {
						my $longest = &findmaxvaluehash(%longestassem);
						&blast_file ($gdir,"$directory.extract.seq.fa",$directory,$outpath,$longest);
						my $blastgenomefile = "$outpath/$directory.extract.seq.fa.$longest.tabular.blast.out";
						my $genomeblast = &parse_blast($blastgenomefile,"95",$longest);
						#print "discsplit genomeblast\n";
						#print Dumper %$genomeblast;
						#predict the orientation of TE insertion
						&compare_twoHOH($directory,$TEblastmini,$genomeblast) if ($longest eq "mini");
						$success=1;
					} else {
						&blast_file ($gdir,"$directory.extract.seq.fa",$directory,$outpath,"mini");
						my $blastgenomefile = "$outpath/$directory.extract.seq.fa.mini.tabular.blast.out";
						my $genomeblast = &parse_blast($blastgenomefile,"95","mini");
						#print "discsplit genomeblast\n";
						#print Dumper %$genomeblast;
						#predict the orientation of TE insertion
						&compare_twoHOH($directory,$TEblastmini,$genomeblast);
						$success=1;
					}
				} 
 			}
			
		} else {
			push (@zerofiles,$R1out);
		} 
		
	}

	#print Dumper %extractcontiginfo;
	#print Dumper %alldirectory;

	print  "DONE with prediction of orientation.... \n";

	print  "\n\nstarting to extract TEsequences.... \n";

	$left = 5;
	$right = 3;
	@matchingtsdlist=();
	@full_lengTEs=();
	@partial_teseq=();
	#my $teseq;
	print "loading sequences\n";
	&loadseq_tohash("Assembled_TEreads");
	#print Dumper %extractcontiginfo;
	#find the contigs with TEsequence and flanking sequence and extract the flanking sequence
	foreach my $locus (sort keys %extractcontiginfo) {
		print "locus is $locus\n";
		my ($typeasm,$genomeloc,$rest) = split (/\./,$locus,3);
		print "genomelocus is $genomeloc\n";
		my ($leftkmer);
		my ($rightkmer);
		my $leftflankseq;
		my $rightflankseq;
		my $start;
		my $end;
		my $strand;
		my $seq;
		my $qTE;
		my $newlocus;
		my $teseq;
		my $testart;
		my $teend;
		my $tsdseq;
		my $tsdsucc;
		my $contlen;
		foreach my $info (keys %{$extractcontiginfo{$locus}}) {
			#print "locus is $info\n";
			$start = $extractcontiginfo{$locus}{$info} if ($info eq 'sstart');
			$end = $extractcontiginfo{$locus}{$info} if ($info eq 'send');
			$strand = $extractcontiginfo{$locus}{$info} if ($info eq 'strand');
			#print STDERR " strand is $strand\n"; 
			$seq = $extractcontiginfo{$locus}{$info} if ($info eq 'seq');
			$qTE = $extractcontiginfo{$locus}{$info} if ($info eq 'qid');
			$contlen = $extractcontiginfo{$locus}{$info} if ($info eq 'slen');
		}
		$newlocus = $locus.".".$start.".".$end.".".$strand if ($strand eq 'plus');
		$newlocus = $locus.".".$end.".".$start.".".$strand if ($strand eq 'minus');
		my $seq_obj = Bio::Seq->new( -display_id => $newlocus, -seq => $seq); # create object with seq
		my $lengthseq = $seq_obj->length;#length of the contig
		if ($strand eq "plus") {# here I am extracting flanking sequence of the TE
			my $leftstart = $start - 40;
			my $leftend = $start + 4;
			my $rightend = $end + 40;
			my $rightstart = $end - 4;
			#print "for a plus $leftstart\t$leftend\t$rightstart\t$rightend\n";
			$leftstart = 1 if ($leftstart <= 0);
			$rightend = $contlen if ($rightend > $contlen);
			$leftflankseq = $seq_obj->subseq($leftstart,$leftend) ;
			$rightflankseq = $seq_obj->subseq($rightstart,$rightend);
		} elsif ($strand eq "minus") {
			my $leftstart = $end - 40;
			my $leftend = $end + 4;
			my $rightend = $start + 40;
			my $rightstart = $start - 4;
		
			#print "for a minus $leftstart\t$leftend\t$rightstart\t$rightend\n";
			$leftstart = 1 if ($leftstart <= 0);
			$rightend = $contlen if ($rightend > $contlen);
			my $leflankseq = $seq_obj->subseq($leftstart,$leftend);
			my $riflankseq = $seq_obj->subseq($rightstart,$rightend);
			$rightflankseq = &revcom_seq($leflankseq);
			$leftflankseq = &revcom_seq($riflankseq);
		}
		#print "$locus left seq is $leftflankseq\n";
		#print "$locus right seq is $rightflankseq\n";
		#identifies TSD sequence
		$leftkmer = &generate_kmer_pos($leftflankseq,$left,$newlocus) if defined ($leftflankseq);
		$rightkmer = &generate_kmer_pos($rightflankseq,$right,$newlocus) if defined ($rightflankseq);
		#print Dumper %$leftkmer;
		#print Dumper %$rightkmer;
		&compare_kmers($leftkmer,$rightkmer);
		my $tsdinfo = $matchingtsdlist[0];
		@matchingtsdlist = (); 
		#print "tsd for $newlocus is $tsdinfo\n";
		if (defined ($tsdinfo)) {
			($testart,$teend,$tsdseq,$tsdsucc) = &findcord_bothtsds($tsdinfo);
			#print "$locus output from findcordbothtsd is $testart, $teend\n";
			$teend = $lengthseq if ($teend > $lengthseq);
			$teseq = $seq_obj->subseq($testart,$teend);
			$teseq = &revcom_seq($teseq) if ($strand eq 'minus');
			my $lenseq = length($teseq);
			if ($tsdsucc eq "yes") {
				my $newID = join (".", $locus,$tsdseq,$qTE,$lenseq);
				&checkprediction($genomeloc,$tsdsucc,$newID);
				my $seqout_obj = Bio::Seq->new( -display_id => $newID, -seq => $teseq);
				push (@full_lengTEs,$seqout_obj);
			} elsif ($tsdsucc eq "no") {
				my $newid =  join (".", $locus,$qTE,"noTSD",$lenseq);
				&checkprediction($genomeloc,$tsdsucc,$newid);
				my $seqout_obj = Bio::Seq->new( -display_id => $newid, -seq => $teseq);
				push (@partial_teseq,$seqout_obj);
			}
		} else {
			$teseq = $seq_obj->subseq($start,$end) if ($strand eq 'plus');
			$teseq = $seq_obj->subseq($end,$start) if ($strand eq 'minus');
			$teseq = &revcom_seq($teseq) if ($strand eq 'minus');
			my $lenseq = length($teseq);
			my $new_id =  join (".", $locus,$qTE,"noTSD",$lenseq);
			my $seqout_obj = Bio::Seq->new( -display_id => $new_id, -seq => $teseq);
			push (@partial_teseq,$seqout_obj);
			&checkprediction($genomeloc,"no",$new_id);
		}
	}


	#compare two hashes of alldirectory strand prediction
	#print the strandprediction hash
	@allreadfullTE = @full_lengTEs;#do that opening the file if need to do 
	@allreadpartTE = @partial_teseq;
	&printarray_fasta ("$path/$teout.full_len.fasta","all",@full_lengTEs);
	&printarray_fasta ("$path/$teout.partial.fasta","all",@partial_teseq);
	&find_no_te_cand();
	#print Dumper %strandprediction;
	%allreadpredictionout= %strandprediction;
	&print_keyvalue_HOH("$path/$fineout","all",\%strandprediction);
	
	
	print "------------DONE all read assembly---------------------\n";
}

if ($onlycompare) {
	print "------------starting to compare the output files generated ---------------------\n";
	%splitdiscpredictionout = &load_file_hash("$path/$fineout","displ");
	#print "$path/$fineout.displ\n";
	
	%discpredictionout = &load_file_hash("$path/$fineout","disc");
	
	%allreadpredictionout = &load_file_hash("$path/$fineout","all");
	
	@discsplitfullTE = &load_seq_array("$path/$teout.full_len.fasta","displ");
	@discfullTE = &load_seq_array("$path/$teout.full_len.fasta","disc");
	@allreadfullTE = &load_seq_array("$path/$teout.full_len.fasta","all");
	@discsplitpartTE = &load_seq_array("$path/$teout.partial.fasta","displ");
	@discpartTE = &load_seq_array("$path/$teout.partial.fasta","disc");
	@allreadpartTE = &load_seq_array("$path/$teout.partial.fasta","all");
}

# print Dumper %allreadpredictionout,"all predictionout\n";
# print Dumper %discpredictionout,"disc predictionout\n";
# print Dumper %splitdiscpredictionout,"splitpredictionout\n";

print "checking before passing to subroutine\n";

print "comparing all three outputs for the final output\n";

my ($compout,$fullout,$partout) = &comparethreehashes(\%splitdiscpredictionout,\%discpredictionout,\%allreadpredictionout);

# print Dumper %$compout,"compareout\n";
# 
# print Dumper %$fullout,"fullout\n";
# 
# print Dumper %$partout,"partial\n";
# 
# 
# print Dumper %fullTEallread,"fullTEallread\n" ;
# print Dumper %partTEallread,"partTEallread\n";
# print Dumper %fullTEdiscsplit,"fullTEdiscsplit\n";
# print Dumper %partTEdiscsplit,"partTEdiscsplit\n";
# print Dumper %fullTEdisc,"fullTEdisc\n";
# print Dumper %partTEdisc,"partTEdisc\n";

&extractfromarrary("full",\@allreadfullTE,\%fullTEallread) if (%fullTEallread);
&extractfromarrary("part",\@allreadpartTE,\%partTEallread) if (%partTEallread);
&extractfromarrary("full",\@discsplitfullTE,\%fullTEdiscsplit) if (%fullTEdiscsplit);
&extractfromarrary("part",\@discsplitpartTE,\%partTEdiscsplit) if (%partTEdiscsplit);
&extractfromarrary("full",\@discfullTE,\%fullTEdisc) if (%fullTEdisc);
&extractfromarrary("part",\@discpartTE,\%partTEdisc) if (%partTEdisc);
&printarray_fasta ("$path/$teout.full_len.fasta","summary",@fullTE_allmethods);
&printarray_fasta ("$path/$teout.partial.fasta","summary",@partTE_allmethods);
&printkeyvalue("$path/$fineout","summary",$compout);

print "------------DONE with the orientTEextractTE script---------------------\n";
#&email;
exit;
#-----------------------------------------------------------------------------
#----------------------------------- SUB -------------------------------------
#-----------------------------------------------------------------------------
sub load_file_hash {
	my $fipath = shift;
	my $typef = shift;
	my %hs=();
	my $filetopen = $fipath.".".$typef;
	print "$filetopen\n";
	open (my $dh, "<", "$filetopen") or confess "\n ERROR (sub): could not open to read $filetopen $!\n";
		while (my $data = <$dh>) {
			chomp $data;
			#print "$data\n";
			my @cols = split(/\s+/,$data);
			#print "@cols\n";
			my $geloc= $cols[0];
			#print "$geloc\n";
			$hs{$geloc}{'fulllengTE'} = $cols[1];
			$hs{$geloc}{'seq_id'} = $cols[3];
			$hs{$geloc}{'strand'} = $cols[2];
		
		}
		print Dumper%hs,"within sub\n"; 
	
	close ($dh);
	return (%hs);
}
sub load_seq_array{
	my $fipath = shift;
	my $typef = shift;
	my $filetopen = $fipath.".".$typef;
	my @seqArray = ();
	my $readio_obj = Bio::SeqIO->new(-file 	=> "$filetopen", 
									-format     => "fasta" 
									) or die "\t    ERROR - Failed to create SeqIO FH object from $filetopen $!\n";  
		while (my $seque = $readio_obj->next_seq() ){
			push @seqArray,$seque;
		}
		
	return (@seqArray);
}
sub checkprediction {
	my ($gloca,$fTEstatus,$ID) = @_;
	if (exists ($strandprediction{$gloca}))  {
		$strandprediction{$gloca}{'fulllengTE'} = $fTEstatus;
		$strandprediction{$gloca}{'seq_id'} = $ID;
	} else {
		$strandprediction{$gloca}{'fulllengTE'} = $fTEstatus;
		$strandprediction{$gloca}{'strand'} = '?';
		$strandprediction{$gloca}{'seq_id'} = $ID;
	}
}
sub load_file {
	my ($file1) = @_;
	my %stefile;
	open (my $th, "<", $file1) or confess "\n ERROR (sub): could not open to read $file1 $!\n";
		while (my $data = <$th>) { #reading the table
			chomp $data; #storing the table values to the $data and removing the extraline
			my @name = split(/\s+/,$data); # splitting the data based on tab and storing into the arrray
			my $dirc = $name[0];#first column is genome location
			$stefile{$dirc} = $name[1]; #name of the TE inserted is loaded 
		}
	return (%stefile);
	close $th;	
}
sub blast_file {
	my ($qpath,$qloc,$directory,$dbpath,$type) = @_;
	#unless (-e "$dbpath/Renamed_Assembledseq/$directory.$type.rename.fasta.nhr") {
	system ("$BLASTdir/makeblastdb -in $dbpath/Renamed_Assembledseq/$directory.$type.rename.fasta -input_type fasta -dbtype nucl") == 0 or die ("unable to makeblastdb on $dbpath/Renamed_Assembledseq/$directory.$type.rename.fasta \n");
	#}
	if ($qloc eq "$directory.extract.seq.fa") {
		system ("$BLASTdir/blastn -db $dbpath/Renamed_Assembledseq/$directory.$type.rename.fasta -query $qpath/$qloc -evalue 0.001 -outfmt \"6 qseqid sseqid pident qlen length slen sstrand qstart qend sstart send qcovs mismatch gapopen evalue bitscore \" -out $dbpath/$qloc.$type.tabular.blast.out");
	} else {
		system ("$BLASTdir/blastn -db $dbpath/Renamed_Assembledseq/$directory.$type.rename.fasta -query $qpath/$qloc -evalue 0.001 -outfmt \"6 qseqid sseqid pident qlen length slen sstrand qstart qend sstart send qcovs mismatch gapopen evalue bitscore \" -out $dbpath/$directory.$qloc.$type.tabular.blast.out");
	}
}
sub findmaxvaluehash {
	my (%hashdata) = @_;
	my $keymax;
	while (my ($key, $value) = each %hashdata) {
		if (not defined $keymax) {
			$keymax = $key;
			next;
		}
		if ($hashdata{$keymax} < $value) {
			$keymax = $key;
		}
	}
	return $keymax;
}	
sub concatenatefiles {
	my ($outputfile, @files) = @_;
	foreach $genomloc (@files){
		chomp ($genomloc);
		#print STDERR "the file now analysing is $genomloc\n";
		#system ("cat $rdir/$directory/$genomloc >> $out") == 0 or die ("unable to concatenate $genomloc \n");
		system ("cat $genomloc >> $outputfile") == 0 or die ("unable to concatenate $genomloc \n");#when using glob path is included so no need to give path
	}
}
sub find_TEseq {
	my ($gloc) = @_;
	my ($c,$starts,$ends) = split (/[:-]/,$gloc);
	my $brkp = $starts +250;
	my $gloca = $c."_".$brkp;
	#print STDERR "the sequence TE is $gloc";
	my $TEseqfile;
	if (exists ($TEinfo{$gloca}))  {
		$TEseqfile = $TEinfo{$gloca};
	}
	else {
		die  "Not able to get the TEseq file \n";
	}
	return ($TEseqfile);
}
sub renameseq_filename {
	my ($contigfile,$fpath,$directory,$type) = @_;
	#print STDERR "$fpath and $contigfile\n";
	make_path  ("$fpath/Renamed_Assembledseq");
	open (my $bhout, ">","$fpath/Renamed_Assembledseq/$directory.$type.rename.fasta") or die "\n ERROR (main): could not open to write $fpath/Renamed_Assembledseq/$directory.$type.rename.fasta $!\n";
	open (my $bh, "<", "$fpath/$contigfile") or confess "\n ERROR (main): could not open to read $contigfile $!\n";
		while(my $dataline = <$bh>) {
			chomp($dataline);
			#print STDERR "$line\n";
				if ($dataline =~ m/^\>\w+\d+/) {
					my $header = $dataline;
					$header =~s/\>//g;
					$header =~ s/\./\_/g;
					$dataline = ">".$directory.".".$header;#dipspades contigs had a different format for fasta header
					#$dataline =~ s/^\>(\w+\d+)/\>$directory\.$1/;
					print $bhout "$dataline\n";
				}
				else {
					print $bhout "$dataline\n";
				}
		}
	close $bh;
	close $bhout;
}
sub parse_blast {
	#parse blast tabular output unique hit with highest blast score
	my $blast = shift;
	my $idlimit = shift;
	my $assembtype = shift;
	my %d = ();
	my $questar;
	my $quend;
	open (my $bl, "<", $blast) or confess "\n ERROR (sub parse_blast): could not open to read $blast $!\n";
	LINE: while (my $line = <$bl>) {
		chomp $line;
		next LINE if ($line !~ /\w/);#  0      1        2    3     4         5     6      7       8      9   10   11    12       13      14     15
 		my @col = split(/\t/,$line); # qseqid sseqid pident qlen aln_length slen sstrand qstart qend sstart send qcovs mismatch gapopen evalue bitscore
		my $contig = $assembtype.".".$col[1];
		my $score = $col[15];
		my $strand = $col[6];
		my $pident = $col[2];
		my $querycov = $col[11];
		my $qstart = $col[7];
		my $qend = $col[8];
		my $qlen = $col[3];
		my $percequerycov = ($querycov/$qlen)*100 ;
		my $qcutoff;
		if ($qlen >= 500) {
			$qcutoff =4.8;#minimum 25bps homology 
		} else {
			$qcutoff =20;
		}
		
		#if 	(($pident > 90) && ($querycov > 40))  {	#need to modify if the stringency has to be increased
		if 	(($pident > $idlimit) && ($percequerycov > $qcutoff))   {	
			
			if (exists ($d{$contig})) {
				my $currentscore = $d{$contig}{'score'};
				#print "current score is $currentscore\n";
				#print "score is $score\n";
				if ($currentscore > $score){
					next;
				} else {
					$d{$contig}{'score'} = $score;
					$d{$contig}{'strand'} = $strand;
					$d{$contig}{'sstart'}  = $col[9];
					$d{$contig}{'send'}  = $col[10];
					$d{$contig}{'qid'} = $col[0];
					$d{$contig}{'slen'} = $col[5];
					$d{$contig}{'qstart'} = $qstart;
					$d{$contig}{'qend'} = $qend; 
				}
			} else {
				$d{$contig}{'score'}  = $score;
				$d{$contig}{'strand'}  = $strand;
				$d{$contig}{'sstart'}  = $col[9];
				$d{$contig}{'send'}  = $col[10];
				$d{$contig}{'qid'} = $col[0];
				$d{$contig}{'slen'} = $col[5];
				$d{$contig}{'qstart'} = $qstart;
				$d{$contig}{'qend'} = $qend; 
			
			}
		}
	}
	#my $noofdisread = keys (%d);
	foreach my $contigkey (sort {$d{$b}->{'score'} <=> $d{$a}->{'score'} } keys (%d)) {
		#print "$contigkey\n";
		foreach my $features (keys %{$d{$contigkey}}) {#derefence second level
		#print "$features\n";
			$questar = $d{$contigkey}->{$features} if ($features eq 'qstart');
			$quend = $d{$contigkey}->{$features} if ($features eq 'qend');
				#print STDERR "TEstrand is $TEstrand and gstrand is $genomestrand \n";
		}
		last;
	}
	return ($questar,$quend,\%d);
}	
sub compare_twoHOH {
	my ($genomloc,$TEhash,$ghash) = @_;
	my %TEhash = %$TEhash;
	my %ghash =%$ghash;
	my $match;
	print Dumper %ghash,"\n";
	
	foreach my $contigkey (sort {$TEhash{$b}->{'score'} <=> $TEhash{$a}->{'score'} } keys (%TEhash)) {#sort on values of HoH
		print STDERR "key in TEH is $contigkey\n";
		my $TEstrand;
		my $genomestrand;
		if (exists ($ghash{$contigkey})) {
			print STDERR " $contigkey exists in ghash\n";
			#%extractcontiginfo = %{ $TEhash{$contigkey} };# copying only next level of keys and values
			$extractcontiginfo{$contigkey} = $TEhash{$contigkey};
			$match = 1;
			foreach my $features (keys %{$TEhash{$contigkey}}) {#derefence second level
				$TEstrand = $TEhash{$contigkey}->{$features} if ($features eq 'strand');
				$genomestrand = $ghash{$contigkey}->{$features} if ($features eq 'strand');
				print STDERR "TEstrand is $TEstrand and gstrand is $genomestrand \n";
			}
			if ((($TEstrand eq 'plus') && ($genomestrand eq 'minus')) || (($TEstrand eq 'minus') && ($genomestrand eq 'plus'))) {
				$strandprediction{$genomloc}{'strand'} = '-';	
			} elsif ((($TEstrand eq 'plus') && ($genomestrand eq 'plus')) || (($TEstrand eq 'minus') && ($genomestrand eq 'minus'))) {	
				$strandprediction{$genomloc}{'strand'} = '+';	
			}
			last;		
		} else {
			next;
		}
	}
	unless ($match) {
		foreach my $contigkey (sort {$TEhash{$b}->{'score'} <=> $TEhash{$a}->{'score'} } keys (%TEhash)) {
			$extractcontiginfo{$contigkey} = $TEhash{$contigkey};
			last;#you will get the one largest score
		}
	}
	return (\%strandprediction,\%extractcontiginfo);
}
sub find_no_te_cand {
	foreach my $locus (sort keys %alldirectory) {
		if (exists ($strandprediction{$locus})) {
			next;
		} else {
			$strandprediction{$locus}{'fulllengTE'} = 'no';
			$strandprediction{$locus}{'strand'} = '?';
			$strandprediction{$locus}{'seq_id'} = 'noTE';
		}	
	}
	return (\%strandprediction);
}
sub print_keyvalue_HOH {
	my ($fl,$type,$hashtoprint) = @_;
	my %hashtoprint =%$hashtoprint;
	my $filepath ="$fl.$type"; 
	open (my $kp,">","$filepath") || die ("failed to open file to write $filepath $!\n");
	foreach my $mate (sort keys %hashtoprint) {
		#splitting to get TSD separate
		if ($hashtoprint{$mate}{'fulllengTE'} eq "yes") {
			my $info = $hashtoprint{$mate}{'seq_id'};
			my @features = split (/\./, $info);
			my $tsdfeat = $features[1];
			print $kp "$mate\t$hashtoprint{$mate}{'fulllengTE'}\t$hashtoprint{$mate}{'strand'}\t$hashtoprint{$mate}{'seq_id'}\t$tsdfeat\n";
		} else {
			print $kp "$mate\t$hashtoprint{$mate}{'fulllengTE'}\t$hashtoprint{$mate}{'strand'}\t$hashtoprint{$mate}{'seq_id'}\n";
		}	
	}
	close $kp;
}
sub printkeyvalue {
	my ($fl,$type,$hashtoprint) = @_;
	my %hashtoprint =%$hashtoprint;
	my $filepaths ="$fl.$type"; 
	open (my $jp,">","$filepaths") || die ("failed to open file to write $filepaths $!\n");
	foreach my $loci (sort keys %hashtoprint) {
		#splitting to get TSD separate
		my $locidetails = $hashtoprint{$loci};
		my @locdet = split (/\t/, $locidetails);
		if ($locdet[0] =~ /yes/) {
			print $jp "$loci\t$locidetails";
		} else {
			print $jp "$loci\t$locdet[0]\t$locdet[1]\t$locdet[2]";
		}
	}
	close $jp;
}
sub comparethreehashes {
	#my ($hash1,$hash2,$hash3) = @_;
	my ($hash1) = shift;#allread changed to #discsplit
	my ($hash2) = shift;#discread
	my ($hash3) = shift;#discsplit changed to #allread
	my %fulllengthIDs;
	my %partialIDs;

	my %comparedoutput;
	my %hash1 = %$hash1;
	my %hash2 = %$hash2;
	my %hash3 = %$hash3;

	# print Dumper %hash1,"\n";
	# print Dumper %hash2,"\n";
	# print Dumper %hash3,"\n";


		foreach my $mate (sort keys %hash1) {
			print "the locus analysing is $mate\n";
			my $info;
			my @features;
			my $tsdfeat;
			my $hash1strand;
			my $hash1fte;
			my $hash1seqid;
			my $hash1length;
			my $hash2strand;
			my $hash2length;
			my $hash2fte;
			my $hash2seqid;
			my $hash3strand;
			my $hash3fte;
			my $hash3seqid;
			my $hash3length;
			my $hashseqid;
			my $hashstrand;
			my $details;
			my $detailspart1;
			my $fleng;
			my $partle;
			my $telength;
			my $identifier;
			my $telength1;
			my $telength2;
			my $telength3;
			my $hashlength;
			my $hash_seqid;
			my $strand1;
			my $strand2;
			my $strand3;
			#checking hash1 alread
			if (($hash1{$mate}{'fulllengTE'} eq "yes") && (($hash1{$mate}{'strand'} eq "+") ||($hash1{$mate}{'strand'} eq "-")) ) {
				#print "checking in hash1\n";
				$info = $hash1{$mate}{'seq_id'};
				
				@features = split (/\./, $info);
				$tsdfeat = $features[3];
				$telength1 = $features[-1];
				$telength = $telength1;
				$strand1 = $hash1{$mate}{'strand'};
				$details = "$hash1{$mate}{'fulllengTE'}\t$hash1{$mate}{'strand'}\t$hash1{$mate}{'seq_id'}\t$tsdfeat\n";
				$identifier = 1;
				$hash_seqid = $hash1{$mate}{'seq_id'};
			} else {
				$hash1fte = $hash1{$mate}{'fulllengTE'};
				$hash1strand = $hash1{$mate}{'strand'};
				$hash1seqid = $hash1{$mate}{'seq_id'};
				my @nameinfo = split (/\./, $hash1seqid);
				$hash1length =  $nameinfo[-1];
				$strand1 = $hash1strand if ($hash1strand ne '?');
				#$identifier = 1;
				print "$hash1seqid is the info in comparethreehashes hash1\n";
				#print "going to check in hash2 and hash3\n";
				#$hash1tsd = $tsdfeat if ($hash1fte eq "yes");
				#if (exists $hash2{$mate}) {
			}
			#checking hash2 discread
			if (($hash2{$mate}{'fulllengTE'} eq "yes") && (($hash2{$mate}{'strand'} eq "+") ||($hash2{$mate}{'strand'} eq "-"))) {
				#print "checking in hash2\n";
				$info = $hash2{$mate}{'seq_id'};
				@features = split (/\./, $info);
				$tsdfeat = $features[3];
				$telength2 = $features[-1];
				$strand2 = $hash2{$mate}{'strand'};
				$identifier = 2;
				#$details = "$hash2{$mate}{'fulllengTE'}\t$hash2{$mate}{'strand'}\t$hash2{$mate}{'seq_id'}\t$tsdfeat\n";
				if (($telength) && ($telength2 > $telength))  {
				
					$details = "$hash2{$mate}{'fulllengTE'}\t$hash2{$mate}{'strand'}\t$hash2{$mate}{'seq_id'}\t$tsdfeat\n";
					$identifier = 2 ;
					$telength = $telength2;
					$hash_seqid = $hash2{$mate}{'seq_id'};
					#print "checking in hash2 gt hash1\n";
				#} elsif (($telength) && ($telength2 <= $telength)) {
					#next;
					#print "checking in hash2 ltoreq hash1\n";
					#$details = "$hash1{$mate}{'fulllengTE'}\t$hash1{$mate}{'strand'}\t$hash1{$mate}{'seq_id'}\t$tsdfeat\n";
					#$identifier = 1 ;
				} elsif (! defined $telength) {
					#print "checking in hash1 notdefined hash1\n";
					$details = "$hash2{$mate}{'fulllengTE'}\t$hash2{$mate}{'strand'}\t$hash2{$mate}{'seq_id'}\t$tsdfeat\n";
					$identifier =2;
					$telength = $telength2;
					$hash_seqid = $hash2{$mate}{'seq_id'};
				}
				
			} else {
				$hash2fte = $hash2{$mate}{'fulllengTE'};
				$hash2strand = $hash2{$mate}{'strand'};
				$hash2seqid = $hash2{$mate}{'seq_id'} ;
				$strand2 = $hash2strand if ($hash2strand ne '?');
				my @nameinfo = split (/\./, $hash2seqid);
				$hash2length =  $nameinfo[-1];
				#$identifier =2;
				#print " notdefined in hash2\n";
				#$hash2tsd = $tsdfeat if ($hash2fte eq "yes");
			}
			#print "going to check in  hash3\n";
			if (($hash3{$mate}{'fulllengTE'} eq "yes") && (($hash3{$mate}{'strand'} eq "+") ||($hash3{$mate}{'strand'} eq "-"))) {
				#print "checking in hash3\n";
				$info = $hash3{$mate}{'seq_id'};
				@features = split (/\./, $info);
				$tsdfeat = $features[3];
				$telength3 = $features[-1];
				$strand3 = $hash3{$mate}{'strand'};
				$identifier =3;
				#$details = "$hash3{$mate}{'fulllengTE'}\t$hash3{$mate}{'strand'}\t$hash3{$mate}{'seq_id'}\t$tsdfeat\n";
				if (($telength) && ($telength3 > $telength))  {
					$details = "$hash3{$mate}{'fulllengTE'}\t$hash3{$mate}{'strand'}\t$hash3{$mate}{'seq_id'}\t$tsdfeat\n";
					$identifier =3;
					$telength = $telength3;
					$hash_seqid = $hash3{$mate}{'seq_id'};
					#print "checking in hash3 gt hash2\n";
				#} elsif (($telength) && ($telength3 <= $telength)) {
					#$details = $details;
					#next;
					#$details = "$hash2{$mate}{'fulllengTE'}\t$hash2{$mate}{'strand'}\t$hash2{$mate}{'seq_id'}\t$tsdfeat\n";
					#$identifier =2;
				} elsif (!defined $telength) {
					$details = "$hash3{$mate}{'fulllengTE'}\t$hash3{$mate}{'strand'}\t$hash3{$mate}{'seq_id'}\t$tsdfeat\n";
					$identifier =3;
					$telength = $telength3;
					$hash_seqid = $hash3{$mate}{'seq_id'};
				}
			
			} else {
				$hash3fte = $hash3{$mate}{'fulllengTE'};
				$hash3strand = $hash3{$mate}{'strand'};
				$hash3seqid = $hash3{$mate}{'seq_id'} ;
				my @nameinfo = split (/\./, $hash3seqid);
				$hash3length =  $nameinfo[-1];
				$strand3 = $hash3strand if ($hash3strand ne '?');
				#$identifier =3;
				#$hash3tsd = $tsdfeat if ($hash3fte eq "yes");
			}
			
			if ($details) {
				my $finalstrand;
				print "identifier is in station 1 is for $mate $identifier\n ";
				$fullTEdiscsplit{$hash_seqid} = 1 if ($identifier ==1);#here is the problem not the correct id is populated
				$fullTEdisc{$hash_seqid} = 1 if ($identifier ==2);
				$fullTEallread{$hash_seqid} = 1 if ($identifier ==3);
				if ((($strand1) && ($strand2)) && ($strand1 ne $strand2)) {
					$finalstrand = "?";
					
				}
				if ((($strand2) && ($strand3)) && ($strand2 ne $strand3)) {
					$finalstrand = "?";
				}
				
				if ((($strand1) && ($strand3)) && ($strand1 ne $strand3)) {
					$finalstrand = "?";
				}
				if ($finalstrand eq '?') {
					my @moddetails = split(/\t/, $details);
					$moddetails[1] = "?";
					$details = "$moddetails[0]\t$moddetails[1]\t$moddetails[2]\t$moddetails[3]\n";
				}
				
			}	
		
			#print "$details after first round of check\n";
			unless ($details) {
			print "identifier is in station 2 is for $mate $identifier\n ";
				if ((($hash1fte eq "yes") || ($hash2fte eq "yes") || ($hash3fte eq "yes")) && (($hash1strand ne "?") || ($hash2strand ne "?" ) || ($hash3strand ne "?"))) {
					if ($hash1fte eq "yes") {
						$hashseqid = $hash1seqid ;
						$hashlength = $hash1length ;
						$identifier =1;
					} 
					if ($hash2fte eq "yes") {
						if (($hashlength) && ($hash2length > $hashlength)) {
							$hashseqid = $hash2seqid ;
							$hashlength = $hash2length;
							$identifier =2;
						#} elsif (($hashlength) && ($hash2length <= $hashlength)) {
							#next;
							#print "do nothing\n";
						} elsif (!defined $hashlength) {
							$hashseqid = $hash2seqid ;
							$hashlength = $hash2length;
							$identifier = 2;
						}
					 }
					 if ($hash3fte eq "yes") {
						if (($hashlength) && ($hash3length > $hashlength)) {
							$hashseqid = $hash3seqid ;
							$hashlength = $hash3length;
							$identifier = 3;
						#} elsif (($hashlength) && ($hash3length <= $hashlength)) {
							#next;
							#print "do nothing\n";
						} elsif (!defined $hashlength) {
							$hashseqid = $hash3seqid ;
							$hashlength = $hash3length;
							$identifier = 3;
						}
					 }
					my @feat = split (/\./, $hashseqid);
					my $tsdf = $feat[3];
					$hashstrand = $hash1strand if ($hash1strand ne "?");
					$hashstrand = $hash2strand if ($hash2strand ne "?");
					$hashstrand = $hash3strand if ($hash3strand ne "?");
					
					if ((($hash1strand ne "?") && ($hash2strand ne "?"))  && ($hash1strand eq $hash2strand)) {
						$hashstrand = $hash1strand;
					}
					
					if ((($hash2strand ne "?") && ($hash3strand ne "?"))  && ($hash2strand eq $hash3strand)) {
						$hashstrand = $hash2strand;
					}
					
					if ((($hash1strand ne "?") && ($hash3strand ne "?"))  && ($hash1strand eq $hash3strand)) {
						$hashstrand = $hash3strand;
					}
					$hashstrand ="?" if (not defined $hashstrand);
					$detailspart1 = "yes\t$hashstrand\t$hashseqid\t$tsdf\n";
					$fleng = "yes\t$hashstrand\t$hashseqid\t$tsdf\n";
					#print "$detailspart1 not found first round\n";
					if ($hashseqid) {
						$fullTEdiscsplit{$hashseqid} = 1 if ($identifier ==1);
						$fullTEdisc{$hashseqid} = 1 if ($identifier ==2);
						$fullTEallread{$hashseqid} = 1 if ($identifier ==3);
					 }
				} elsif ((($hash1fte eq "yes") || ($hash2fte eq "yes") || ($hash3fte eq "yes")) && (($hash1strand eq "?") && ($hash2strand eq "?" ) && ($hash3strand eq "?"))) {
				
					if ($hash1fte eq "yes") {
						$hashseqid = $hash1seqid;
						$hashlength = $hash1length;
						$identifier = 1;
					} 
				
					if ($hash2fte eq "yes") {
						if (($hashlength) && ($hash2length > $hashlength)) {
							$hashseqid = $hash2seqid;
							$hashlength = $hash2length;
							$identifier = 2;
						#} elsif (($hashlength) && ($hash2length <= $hashlength)) {
							#next;
							#print "do nothing\n";
						} elsif (!defined$hashlength) {
							$hashseqid = $hash2seqid;
							$hashlength = $hash2length;
							$identifier = 2;
						}
					 }
				 
					 if ($hash3fte eq "yes") {
						if (($hashlength) && ($hash3length > $hashlength)) {
							$hashseqid = $hash3seqid ;
							$hashlength = $hash3length;
							$identifier = 3;
						#} elsif (($hashlength) && ($hash3length <= $hashlength)) {
							#next;
							#print "do nothing\n";
						} elsif (!defined$hashlength) {
							$hashseqid = $hash3seqid;
							$hashlength = $hash3length;
							$identifier = 3;
						}
					 }
		
					my @feat = split (/\./, $hashseqid);
					my $tsdf = $feat[3];
					$detailspart1 = "yes\t?\t$hashseqid\t$tsdf\n";
				
					#print "$detailspart1 detailspart1\n";
					$fleng = "yes\t$hashstrand\t$hashseqid\t$tsdf\n";
					if ($hashseqid) {
						$fullTEdiscsplit{$hashseqid} = 1 if ($identifier ==1);
						$fullTEdisc{$hashseqid} = 1 if ($identifier ==2);
						$fullTEallread{$hashseqid} = 1 if ($identifier ==3);
					}
				 } elsif ((($hash1fte eq "no") && ($hash2fte eq "no") && ($hash3fte eq "no")) && (($hash1strand ne "?") || ($hash2strand ne "?" ) || ($hash3strand ne "?"))) {
				
					if ($hash1seqid ne "noTE") {
						$hashseqid = $hash1seqid;
						$hashlength = $hash1length;
						$identifier = 1;
					}
					if ($hash2seqid ne "noTE") {
						if (($hashlength) && ($hash2length > $hashlength)) {
							$hashseqid = $hash2seqid;
							$hashlength = $hash2length;
							$identifier = 2;
						#} elsif (($hashlength) && ($hash2length <= $hashlength)) {
							#next;
							#print "do nothing\n";
				
						} elsif (!defined$hashlength) {
							$hashseqid = $hash2seqid;
							$hashlength = $hash2length;
							$identifier = 2;
						}
					}
					if ($hash3seqid ne "noTE") {
						if (($hashlength) && ($hash3length > $hashlength)) {
							$hashseqid = $hash3seqid ;
							$hashlength = $hash3length;
							$identifier = 3;
						#} elsif (($hashlength) && ($hash3length <= $hashlength)) {
							#next;
							#print "do nothing\n";
						} elsif (!defined $hashlength) {
							$hashseqid = $hash3seqid;
							$hashlength = $hash3length;
							$identifier = 3;
						}
					} 
					$hashstrand = $hash1strand if ($hash1strand ne "?");
					$hashstrand = $hash2strand if ($hash2strand ne "?");
					$hashstrand = $hash3strand if ($hash3strand ne "?");
					#if you want majority then change ne to eq and hashstrand eq one of it.
					if ((($hash1strand ne "?") && ($hash2strand ne "?"))  && ($hash1strand ne $hash2strand)) {
						$hashstrand = "?";
					}
					
					if ((($hash2strand ne "?") && ($hash3strand ne "?"))  && ($hash2strand ne $hash3strand)) {
						$hashstrand = "?";
					}
					
					if ((($hash1strand ne "?") && ($hash3strand ne "?"))  && ($hash1strand ne $hash3strand)) {
						$hashstrand = "?";
					}
					
					# if ((($hash1strand ne "?") && ($hash2strand ne "?") && ($hash3strand ne "?")) && ($hash1strand ne $hash2strand))  {
# 						$hashstrand = "?";
# 					}
# 					if ((($hash1strand ne "?") && ($hash2strand ne "?") && ($hash3strand ne "?")) && ($hash2strand ne $hash3strand))  {
# 						$hashstrand = "?";
# 					}
# 					
# 					if ((($hash1strand ne "?") && ($hash2strand ne "?") && ($hash3strand ne "?")) && ($hash1strand ne $hash3strand))  {
# 						$hashstrand = "?";
# 					}
					
					$detailspart1 = "no\t$hashstrand\t$hashseqid\n";
				

					$partle = "no\t$hashstrand\t$hashseqid\n";
					if ($hashseqid) {
						$partTEdiscsplit{$hashseqid} = 1 if ($identifier ==1);
						$partTEdisc{$hashseqid} = 1 if ($identifier ==2);
						$partTEallread{$hashseqid} = 1 if ($identifier ==3);
					}
				
				
				 } elsif ((($hash1fte eq "no") && ($hash2fte eq "no") && ($hash3fte eq "no")) && (($hash1strand eq "?") && ($hash2strand eq "?" ) && ($hash3strand eq "?"))) {
				
							
					if ($hash1seqid ne "noTE") {
						$hashseqid = $hash1seqid;
						$hashlength = $hash1length;
						$identifier = 1;
					}
					if ($hash2seqid ne "noTE") {
						if (($hashlength) && ($hash2length > $hashlength)) {
							$hashseqid = $hash2seqid;
							$hashlength = $hash2length;
							$identifier = 2;
						} elsif (!defined $hashlength) {
							$hashseqid = $hash2seqid;
							$hashlength = $hash2length;
							$identifier = 2;
						}
					}
					if ($hash3seqid ne "noTE") {
						if (($hashlength) && ($hash3length > $hashlength)) {
							$hashseqid = $hash3seqid ;
							$hashlength = $hash3length;
							$identifier = 3;
						} elsif (!defined $hashlength) {
							$hashseqid = $hash3seqid;
							$hashlength = $hash3length;
							$identifier = 3;
						}
					}
					$detailspart1 = "no\t?\t$hashseqid\n";
					#print "$detailspart1 detailspart1 allq\n";
					$partle = "no\t?\t$hashseqid\n";
					if ($hashseqid) {
						$partTEdiscsplit{$hashseqid} = 1 if ($identifier == 1);
						$partTEdisc{$hashseqid} = 1 if ($identifier == 2);
						$partTEallread{$hashseqid} = 1 if ($identifier == 3);
					}
				 }
			}
			$comparedoutput{$mate} = $details if ($details);
			$comparedoutput{$mate} = $detailspart1 if ($detailspart1);
			$fulllengthIDs{$mate} = $details if ($details);
			$fulllengthIDs{$mate} = $details if ($fleng);
			$partialIDs{$mate} = $partle if ($partle);
		}
		return (\%comparedoutput,\%fulllengthIDs,\%partialIDs);
}
sub extractfromarrary {
	my $typeoutput = shift;
	my $seqarray = shift;
	my $seqlist= shift;
	my @seqarray = @$seqarray;
	my %seqlist = %$seqlist;
		for (my $l = 0; $l <= ($#seqarray); $l++) {
			my $TEname = $seqarray[$l]->display_id;
			if (exists $seqlist{$TEname}){
				#print $seqArrayR1[$i]->display_id,"is the read name with index $i\n";
				my $TEseq = splice (@seqarray,$l,1);
				push (@fullTE_allmethods, $TEseq) if ($typeoutput eq "full");
				push (@partTE_allmethods, $TEseq) if ($typeoutput eq "part");
				$l = $l-1;
			} else {
				next;
			}
		}	
}
sub printarray_fasta {
	my ($f,$type,@seqarray) = @_;
	
	my $filetoprint	= "$f.$type";	
	my $seqarray_obj = Bio::SeqIO->new(-file   => ">$filetoprint",
										      -format => 'fasta') 
									         or die "\t    ERROR - Failed to create SeqIO FH object from $filetoprint  $!\n";  
	foreach my $seque (@seqarray) {
		$seqarray_obj->write_seq($seque);	
	}
}
sub loadseq_tohash {
	my ($assemtype)=@_;
	
	foreach  my $locus (keys %extractcontiginfo) {
	
		my ($typeasb,$dirinfo,$contigname) = split(/\./,$locus,3); 
		#print "type ofasse is $typeasb, locus is $locus, dirname is $dirinfo and contigname is $contigname\n";
		my $readio_obj = Bio::SeqIO->new(-file 	=> "$path/$assemtype/$dirinfo/Renamed_Assembledseq/$dirinfo.$typeasb.rename.fasta", 
									-format     => "fasta" 
									) or die "\t    ERROR - Failed to create SeqIO FH object from $path/$assemtype/$dirinfo/Renamed_Assembledseq/$dirinfo.$typeasb.rename.fasta $!\n";  
		while (my $seque = $readio_obj->next_seq() ){
			my $header = $seque->display_id;
			my $sequence = $seque->seq;
			my $originallocus = $dirinfo.".".$contigname;
			if ($header =~ /$originallocus$/) {
				$extractcontiginfo{$locus}{'seq'} = $sequence;
			} else {
				next;
			}
		}
	}
}
sub generate_kmer_pos {
	my $kmersize = 9;#changed to 9 to make it more stric
	my ($seq,$direction,$genomeloc) = @_;
	my $lenseq = length $seq;	
	#print "the length of the sequence $lenseq\n";
	#print "the length of the left right $lenleft\n";    
    my %kmers;    
    for (my $j = $kmersize; $j <= 26; $j++) {#index starts from 0
    	my $klen = $j;
		if ($direction == 3) {
			for (my $i = 0; ($i + $j) <= $lenseq; $i++) {
				$kmers{$genomeloc}{$klen}{$i} = substr($seq, $i, $j);
			}
		}
		if ($direction == 5) {# index starts from the kmerlength number
			for (my $i = $j; $i <= $lenseq; $i++) {#to extract kmers from the end to the beginning#check the lenght 
				#my $s = ($i+$j); 
				#print "i is $i\n";
				#print "j is $j\n";
				#print "s is $s\n";
				$kmers{$genomeloc}{$klen}{$i} = substr($seq, -$i, $j ); #This is the logic used #-4,4 #-5,4 #-6,4
			}
		}
	}
    return (\%kmers);
}
sub compare_kmers {
	my ($l, $r) = @_;
	my %lefthash = %$l;
	my %righthash = %$r;
	my %tsdlist;
	my @potentialtsdlist;
	for my $loc (keys %lefthash) {
		for my $klen ( sort { $a <=> $b } keys %{$lefthash{$loc}} ) {#sort keys numerically
			#print "$klen is the first level key\n"; 
			if (exists ($righthash{$loc}{$klen})) {
				for my $lindex (sort keys $lefthash{$loc}{$klen}) {
					#print "$lindex is the second level key for secondhash\n"; 
					for my $rindex (sort keys $righthash{$loc}{$klen}) {
					     my $lefttsd = $lefthash{$loc}{$klen}{$lindex};
					     my $righttsd = $righthash{$loc}{$klen}{$rindex};
						 #if ( $lefthash{$loc}{$klen}{$lindex} eq $righthash{$loc}{$klen}{$rindex} ) {
						 if ( fc($lefttsd) eq fc($righttsd) ) {#to check equality irrespective of the case
							 #print "for $klen is the matching TSD  $lefthash{$loc}{$klen}{$lindex} is at $lindex, $righthash{$loc}{$klen}{$rindex} is at $rindex\n";
							 my $tsdinfo = "$loc.$klen.$lefthash{$loc}{$klen}{$lindex}.$lindex.$rindex";
							 push (@potentialtsdlist,$tsdinfo);
							 #print "$tsdinfo\n";
						 } #elsif (amatch ($lefthash{$g},["i 9%"],$righthash{$m}) ) {
							#print "the approximate matching TSD $lefthash{$g} is at $g, $righthash{$m} at  $m\n";
						 #}
					}			
				}
				#print "@potentialtsdlist\n";
				my $count = @potentialtsdlist;
				if ($count == 1) {
					push (@matchingtsdlist,$potentialtsdlist[0]);
					@matchingtsdlist = uniq @matchingtsdlist ;
					#print STDERR "@matchingtsdlist\n";
				} else {
					@potentialtsdlist = ();
				}		
			} else {
				#print STDERR "not identified in righthash\n"
			}
		}
	}
	return (@matchingtsdlist);
}
sub findcord_bothtsds {
	my ($element) = @_;
	my $genomloc;
	my $tsd1start;
	my $tsd1end;
	my $tsd2start;
	my $tsd2end;
	my $tsdsuccess;
	#addedtype of assembly	
	          	#		0					1					   2   				3  4   5    6		     7  	8 9	
	print "$element\n";#mini.8:128356417-128356917.NODE_2_length_316_cov_59_472803.96.239.minus.15.AGCCAGGCAATTTTT.39.6#renamed to remove . within contig name
	my @details = split (/\./,$element);
	
	$genomloc = $details[1];
	# my $seqinfo = $details[1];
# 	if ($seqinfo =~ /^NODE/) {
# 		my @seqdetails = split (/\_/,$seqinfo);
# 		my $seqlen = $seqdetails[3];
# 		
# 	} elsif ($seqinfo =~ /^\d+/) {
# 	
# 	}
	
	my $te_start = $details[3];
	my $te_end = $details[4];
	my $tstrand = $details[5];
	my $tkmerlen = $details[6];
	my $tsd = $details[7];
	my $lefindex = $details[8];
	my $rigindex = $details[9];
	#my $gID = $genomloc.".".$details[1].".".$details[2];
	my $diff = ($lefindex - $tkmerlen);
	if ($diff < 10) {
		if ($tstrand eq 'plus') {
			#my $gchr = $tchr;
			$tsd1start = (($te_start+4) - $lefindex)+1;
			$tsd1end = ($tsd1start+$tkmerlen);
			$tsd2start = ($te_end-4) + $rigindex;
			$tsd2end = ($tsd2start + $tkmerlen)-1;
			$tsdsuccess = "yes";
			#$gloc = $gchr.":".$tsd1start."-".$tsd1end."\t".$gchr.":".$tsd2start."-".$tsd2end;
			#$TEcor_bothtsd{$gloc} = $tsd;
			#$exactmatch{$gID} = 'exactmatch';
		} elsif ($tstrand eq 'minus') {#leftflank and rightflank is reverse complemented to identify the cordinates need to be modified
			#my $gchr = $tchr;
			
			$tsd1end = (($te_start+4) - $rigindex)+1;#  the sequence is reverse complemented so it is polyTend
			$tsd1start = ($tsd1end - $tkmerlen);
			$tsd2end = (($te_end-4) + $lefindex);
			$tsd2start = ($tsd2end - $tkmerlen)-1 ;#+1 and -1 is to get the TE cordinates after and before TSD
			$tsdsuccess = "yes";
			#$gloc = $gchr.":".$tsd1start."-".$tsd1end."\t".$gchr.":".$tsd2start."-".$tsd2end;
			#my $revtsd = &revcom_seq($tsd);      
			#$TEcor_bothtsd{$gloc} = $revtsd;
			#$exactmatch{$gID} = 'exactmatch';
		}
	} else {
		#print STDERR "TSD with exact match cannot be identified for $element\n";
		$tsdsuccess = "no";
		$tsd1start = $te_start;
		$tsd2end  = $te_end;
	}
	return ($tsd1start,$tsd2end,$tsd,$tsdsuccess);
}
sub revcom_seq {
	my ($inputseq) = @_; 
	my $seqobj = Bio::Seq->new(-seq => "$inputseq", #to obtain reverse compliment of tsd when the strand is minus as I am extracting the reverse compliment of the sequence 
                         	   -alphabet => 'dna' );
	$seqobj = $seqobj->revcom; #reverse compliment the sequence
	my $revtsd = $seqobj->seq;
	return ($revtsd);
}
sub splitsoftdismatetopairs {
	my ($dmpath,$dmfile) = @_;
	my %dmlistR1 = ();
	my %dmlistR2 = ();
	my $totalreadcount;
	my @seqArrayR1 = ();
	my @seqArrayR2 = ();
	my @seqArrayUp = ();
	my $duplicateR1 = 0;
	my $duplicateR2 = 0;
	my $readio_obj = Bio::SeqIO->new(-file 	 => "$dmpath/$dmfile", 
									 -format => 'fasta') 
								 or die "\t    ERROR - Failed to create SeqIO FH object from $dmfile $!\n";  
	while (my $seq = $readio_obj->next_seq() ){#loading the reads in the discordant reads into a hash and also into separate arrays
		#chomp $seq;
		my $header = $seq->display_id;
		if (defined ($header =~ /^(.*)\/(\d)/)) {
			if ($2 == 1) {
				if (exists $dmlistR1{$1}) {#trying to remove duplicate reads because merging with soft clipped reads can create duplicate reads
					$duplicateR1 +=1; 
					next;
				} else {
					$dmlistR1{$1}=1;
					push (@seqArrayR1,$seq) if (defined $seq);
				}
			} elsif ($2 == 2) {
				if (exists $dmlistR2{$1}) {
					$duplicateR2 +=1; 
					next;
				} else {
					$dmlistR2{$1}=1;
					push (@seqArrayR2,$seq) if (defined $seq);
				}
			}
		} 
	}
	#print "Duplicates are $duplicateR1 and $duplicateR1\n";
	#print "$nume\n";
	$totalreadcount = (keys %dmlistR1) + (keys %dmlistR2);
	
	foreach my $record (keys %dmlistR1) {
		if (exists ($dmlistR2{$record})) {#check if the read in unpaired file is present in the discordant read1 list has remove it discordant list
			delete ($dmlistR1{$record});
			delete ($dmlistR2{$record});
		}  else {
			next;
		}
	}
	foreach my $record (keys %dmlistR2) {
		if (exists ($dmlistR1{$record})) {#check if the read in unpaired file is present in the discordant read1 list has remove it discordant list
			delete ($dmlistR1{$record});
			delete ($dmlistR2{$record});
		}  else {
			next;
		}
	}
	
	#print Dumper %dmlistR1,"\n";
	#print Dumper %dmlistR2,"\n";
	#print Dumper @seqArrayR1,"\n";
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
	#print Dumper @seqArrayUp,"\n";
	
	&sortfasta_header($DmatefileR1,@seqArrayR1);
	&sortfasta_header($DmatefileR2,@seqArrayR2);
	&sortfasta_header($DUpfile,@seqArrayUp);
	#print STDERR "$DUpfile printed in $dmpath \n";
	return $totalreadcount;
}
sub splitdismatetopairs {
	my ($dmpath,$dmfile) = @_;
	my %dmlistR1 = ();
	my %dmlistR2 = ();
	my $totalreadcount;
	my @seqArrayR1 = ();
	my @seqArrayR2 = ();
	my @seqArrayUp = ();
	my $readio_obj = Bio::SeqIO->new(-file 	 => "$dmpath/$dmfile", 
									 -format => 'fasta') 
								 or die "\t    ERROR - Failed to create SeqIO FH object from $dmfile $!\n";  
	while (my $seq = $readio_obj->next_seq() ){#loading the reads in the discordant reads into a hash and also into separate arrays
		my $header = $seq->display_id;
		if ($header =~ /^(.*)\/(\d)/) {
			if ($2 == 1) {
				$dmlistR1{$1}=1;
				push (@seqArrayR1,$seq);
			} elsif ($2 == 2) {
				$dmlistR2{$1}=1;
				push (@seqArrayR2,$seq);
			}
		} 
	}
	
	$totalreadcount = (keys %dmlistR1) + (keys %dmlistR2);
	
	foreach my $record (keys %dmlistR1) {
		if (exists ($dmlistR2{$record})) {#check if the read in unpaired file is present in the discordant read1 list has remove it discordant list
			delete ($dmlistR1{$record});
			delete ($dmlistR2{$record});
		}  else {
			next;
		}
	}
	foreach my $record (keys %dmlistR2) {
		if (exists ($dmlistR1{$record})) {#check if the read in unpaired file is present in the discordant read1 list has remove it discordant list
			delete ($dmlistR1{$record});
			delete ($dmlistR2{$record});
		}  else {
			next;
		}
	}
	
	
	 
	if (%dmlistR1) {# puttting the rest of the rest of the sequences in the discordant mate file to unpaired sequences
		foreach my $read (keys %dmlistR1){
			my $ind = firstidx { $_->display_id eq "$read/1" } @seqArrayR1;
			
			my $upseq = splice (@seqArrayR1,$ind,1);
			
			$read = substr $read, 0,-2;
			$upseq->display_id("$read");
			push (@seqArrayUp, $upseq);
		}
	} 
	
	if (%dmlistR2) {
		foreach my $read (keys %dmlistR2){
			my $ind = firstidx { $_->display_id eq "$read/2" } @seqArrayR2;
			
			my $upseq = splice (@seqArrayR2,$ind,1);
			
			$read = substr $read, 0,-2;
			$upseq->display_id("$read");
			
			push (@seqArrayUp, $upseq);
		}
	
	}
	
	
	&sortfasta_header($DmatefileR1,@seqArrayR1);
	&sortfasta_header($DmatefileR2,@seqArrayR2);
	&sortfasta_header($DUpfile,@seqArrayUp);
	#print STDERR "$DUpfile printed in $dmpath \n";
	return $totalreadcount;
}
sub sortfasta_header {
	my ($filetosort,@seqarray) = @_;
		
	my $sortedReadio_obj = Bio::SeqIO->new(-file   => ">$filetosort",
										      -format => 'fasta') 
									         or die "\t    ERROR - Failed to create SeqIO FH object from $filetosort $!\n";  
	@seqarray = sort { ($a->display_id cmp $b->display_id) } @seqarray;
	foreach my $seque (@seqarray) {
		$sortedReadio_obj->write_seq($seque);	
	}
}
sub readlocifromfile {
	my ($locilist) = @_;
	my @locinames;
	open (my $lh, "<", $locilist) or confess "\n ERROR (main): could not open to read $locilist $!\n";
		while (my $data = <$lh>) { #reading the table
			chomp $data; #storing the table values to the $data and removing the extraline
			push @locinames,$data if ($data !~ /^$/); #  storing into the arrray
		}
	return (@locinames);
	close $lh;	
}
sub email {
	my $to = 'jainythomas1@gmail.com';
	my $cc = 'jainyt@genetics.utah.edu';
	my $from = 'jainy@freya.genetics.utah.edu';
	my $subject = 'orientTE_extractTE';
	my $message = "see output in $path..";

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