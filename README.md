# Non-reference-Alu_Assembly
 Assemble non-reference Alu insertions when breakpoints are provided


Usage:
Step 1:
Extract reads from a 250bp window of the breakpoint including mates of the discordant reads

	perl 01_processbam_forreadextract_v19.0.pl -t /path/to/BAMID.txt -f inputfile.txt -p /path/to/outputfolder -bl /path/to/BAMfiles -g /path/to/referencegenome.fa -pt /path/to/picardtools -sq /path/to/seqtk -bu /path/to/bamUtil -bt /path/to/bedtools2/bin


Sample BAMID file

	LP6005441-DNA_A01	LP6005441-DNA_A01.38.sorted.bam
	LP6005441-DNA_A03	LP6005441-DNA_A03.38.sorted.bam
	LP6005441-DNA_A04	LP6005441-DNA_A04.38.sorted.bam
	LP6005441-DNA_A05	LP6005441-DNA_A05.38.sorted.bam	
	
Sample input file (without header):
	BAMID				uniqueid		chr		breakpoint TE	strand(optional) TSD(optional)
	
	LP6005441-DNA_A01	chr1_1564206	chr1	1564206	AluY	+	AAGAATACTC
	LP6005441-DNA_A03	chr1_1564206	chr1	1564206	AluY	+	AAGAATACTC
	LP6005441-DNA_A04	chr1_1564206	chr1	1564206	AluY	+	AAGAATACTC
	LP6005441-DNA_A05	chr1_1564206	chr1	1564206	AluY	+	AAGAATACTC
	LP6005441-DNA_A01	chr1_2564285	chr1	2564285	AluY	+	AGAAGCGGTG
	LP6005441-DNA_A03	chr1_2564285	chr1	2564285	AluY	+	AGAAGCGGTG
	LP6005441-DNA_A04	chr1_2564285	chr1	2564285	AluY	+	AGAAGCGGTG
	LP6005441-DNA_A05	chr1_2564285	chr1	2564285	AluY	+	AGAAGCGGTG

"Parallel"(https://www.gnu.org/software/parallel/) can be used to total decrease time taken for step1. The input file can be split into multiple files and move to a directory (Lets call it "inputdir"). While splitting the input file, have all the candidates from all individual in a single file.  Make a list of input files (call it "list_of_infputfiles"). Now cd inputdir and launch the following command. -j can be changed as per the number of the input files and number of CPUs available.

	cat ../list_of_inputfiles.txt | nohup /usr/bin/parallel -j 4 --results /path/to/stderr 'perl 01_processbam_forreadextract_v18.0.pl -t /path/to/BAMID.txt -f {} -p /path/to/outputfolder -bl /path/to/BAMfiles -g /path/to/referencegenome.fa -pt /path/to/picardtools -sq /path/to/seqtk -bu /path/to/bamUtil -bt /path/to/bedtools2/bin' &

Step 2:
denovo assembly of reads, identify the Alu sequence and extract the assembled sequence using the output from step 1
Create a folder named e.g TEsequences and copy consensus AluY sequences (see sample_files). The consensus AluY sequences is used to identify homologous AluY sequences from the assembled contigs


	perl 04_orientTE_extractTE_v12.6.pl -g /path/to/outputfolder/ExtractGenomicsequences -t /path/to/TEsequences -ds /path/to/outputfolder/Discosplitassembly -dc /path/to/outputfolder/Discoassembly -d /path/to/outputfolder/orientTE -fl /path/to/list_locations.txt -l /path/to/position_TE.txt -bp /path/to/localBlast -cp /path/to/CAP3 -sp /path/to/SPAdes-3.11.1-Linux/bin -mn /path/to/minia-v2.0.7-Source/build/bin -cu numberofcpus -p /path/to/output 

sample list_locations.txt (used only if the locations the in the list are to be analysed)

	chr10:116369944-116370444
	chr10:13533871-13534371
	chr10:23579915-23580415
	chr10:26701998-26702498

Sample position_TE.txt

	chr1_1564206    AluY
	chr1_2564285    AluY
	chr1_2622360    AluY
	chr1_3046212    AluY
	chr1_3095323    AluY
	

