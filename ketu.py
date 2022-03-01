#!/usr/bin/python

"""
All encompassing NGS pipeline in Python go mokas!~
"""

import re
import os
#import compileall

BASEDIR = os.path.dirname(os.path.realpath(__file__));
os.chdir(BASEDIR);
os.environ['BASEDIR'] = os.path.dirname(os.path.realpath(__file__));
#compileall.compile_dir(BASEDIR)

intro="""\
Welcome to Ketu Computational Genomics Toolkit v0.6.29
#                                                 
#      _/    _/  _/_/_/_/  _/_/_/_/_/  _/    _/   
#     _/  _/    _/            _/      _/    _/    
#    _/_/      _/_/_/        _/      _/    _/     
#   _/  _/    _/            _/      _/    _/      
#  _/    _/  _/_/_/_/      _/        _/_/         
#                                                 
#                                              
You will interact in a Q/A format via the terminal, just like classic role-playing games:
"""
print intro;

#read-processor
def proc():
	print "\nKetu Can process assembeled reads to yeild more accurate analsis v0.6.29";

	message = """\
Ketu can perform the following Processing steps
Please select one to continue:

[1] Remove Duplicate Reads
[2] Local Realignment
[3] Base Quality Recalibration
""";
	print message;

	choice = raw_input("\n:");
	if choice == "1":
		print "You have decided to Remove Duplicate Reads\n";

	else:
		if choice == "2":
			print "You have decided to perform Local Realignment\n";

		else:
			if choice == "3":
				print "You have decided to perform Base Quality Recalibration\n";

	input_variable = raw_input("\nTo begin Processing Please drag & drop a BAM file: \n")
	os.environ['input_variable'] = input_variable;
	print "Let's take a quick peek inside our BAM file.."
	os.system("samtools view $input_variable | head -n 5 > del.tmp");
	os.system("more del.tmp");

	print "\n The file looks okay to me.\n"
	os.system("rm del.tmp");

	input_reference = raw_input("Please drag & drop a Reference Genome File:\n")
	os.environ['input_reference'] = input_reference;
	print "Loading Reference Genome from:\n"; print input_reference; print "\nLoading...\n";

	if choice == "1":
		print "Removing reads with identical external coordinates, keeping read with highest mapping quality:\n";
		os.system("java -Xmx4g -jar $BASEDIR/MarkDuplicates.jar INPUT=$input_variable OUTPUT=$input_variable.rmdup.bam ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true");
		os.system("echo Your new BAM file with Duplicate Reads Removed can be found at $input_variable.rmdup.bam");

	else:
		if choice == "2":
			print "Realigning Reads around local regions...";
			os.system("java -Xmx4g -jar $BASEDIR/GenomeAnalysisTK.jar -I $input_variable -R $input_reference -T RealignerTargetCreator -o $BASEDIR/IndelRealigner.intervals");

			os.system("java -Xmx4g -jar $BASEDIR/GenomeAnalysisTK.jar -I $input_variable -R $input_reference -T IndelRealigner -targetIntervals $BASEDIR/IndelRealigner.intervals -o $input_variable.realn.bam");
			os.system("echo Your new BAM file, Realigned around Local Regions, can be found at $input_variable.realn.bam");

			#os.system("rm $BASEDIR/IndelRealigner.intervals");

		else:
			if choice =="3":
				print "Recalibrating Base Qualities based on known human variants";
				input_db = raw_input("Please drag & drop a known database of variants, e.g. dbSNP, Mills gold standard:\n")
				os.environ['input_db'] = input_db;

				os.system("java -Xmx4g -jar $BASEDIR/GenomeAnalysisTK.jar -T BaseRecalibrator -R $input_reference -I $input_variable -knownSites $input_db -o $BASEDIR/recal_data.grp");
				os.system("java -Xmx4g -jar $BASEDIR/GenomeAnalysisTK.jar -T PrintReads -R $input_reference -BQSR $BASEDIR/recal_data.grp -I $input_variable -o $input_variable.recal.bam");

				os.system("java -Xmx4g -jar $BASEDIR/GenomeAnalysisTK.jar -T BaseRecalibrator -R $input_reference -I $input_variable.recal.bam -knownSites $input_db --plot_pdf_file $input_variable.pdf -BQSR $BASEDIR/recal_data.grp -o $BASEDIR/recal_data.after.grp");

				os.system("echo Your new BAM file with Base Quality Recalibration can be found at $input_variable.recal.bam\nA comparison plot of uncalibrated/recalibrated files can be found at $input_variable.pdf");


#variant discover
def disc():
	print "\nKetu can Detect Variants via UnifiedGenotyper & mPileup v0.6.28";

	choice = raw_input("\nPlease select an algorithm to continue:\n[1] UnifiedGenotyper\n[2] mPileup\n")

	if choice == "1":
		print "\nYou have picked UnifiedGenotyper to detect variants.\n";
		detect_algo = "unifiedgeno";

	else:
		if choice == "2":
			print "You have picked mPileup to detect variants.\n";
			detect_algo = "mpileup";

	#file input to run detection algorithms on
	input_variable = raw_input("\nTo begin Discovering variants Please drag & drop a BAM file:\n")
	os.environ['input_variable'] = input_variable;

	#print "\nLet's take a peak inside:\n";
	#os.system("$BASEDIR/samtools view $input_variable | head -n 5 > tmp1.tmp");
	#os.system("more tmp1.tmp");
	#os.system("rm tmp1.tmp");

	message = """\
Checking BAM file...

...

Looks okay to me.

Initiating Discovery modules...

...\n"""

	print message;

	input_reference = raw_input("Please drag & drop a Reference Genome File To Continue Variant Detection:\n")
	os.environ['input_reference'] = input_reference;

	message = """\
Loading...

All parameters are GO\n"""

	print message;

	if detect_algo == "unifiedgeno":
		os.system("java -Xmx4g -jar $BASEDIR/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $input_reference -I $input_variable -o $input_variable.raw.vcf");

	else:
		if detect_algo == "mpileup":
			os.system("$BASEDIR/samtools mpileup -uf $input_reference $input_variable | $BASEDIR/bcftools/bcftools view -bvcg - > $input_variable.raw.bcf");
			os.system("$BASEDIR/bcftools/bcftools view $input_variable.raw.view > $input_variable.raw.vcf");

	os.system("echo Your new VCF file can be found at $input_variable.raw.vcf");


#tail method
def tail( f, window=20 ):
    BUFSIZ = 1024
    f.seek(0, 2)
    bytes = f.tell()
    size = window
    block = -1
    data = []
    while size > 0 and bytes > 0:
        if (bytes - BUFSIZ > 0):
            # Seek back one whole BUFSIZ
            f.seek(block*BUFSIZ, 2)
            # read BUFFER
            data.append(f.read(BUFSIZ))
        else:
            # file too small, start from begining
            f.seek(0,0)
            # only read what was not read
            data.append(f.read(bytes))
        linesFound = data[-1].count('\n')
        size -= linesFound
        bytes -= BUFSIZ
        block -= 1
    return '\n'.join(''.join(data).splitlines()[-window:])


#assembly-controller
def asmb():
	print "\nKetu Can utilize two different Assembly Algorithms v0.6.28";
	choice = raw_input("\nPlease Select an Option to Continue: \n [1] BWA\n [2] Bowtie\n");

	if choice == "1":
		print "\nYou have picked Burrows Wheeler Aligner (BWA) to continue.\n";
		alignment_algo = "bwa";
	else:
		if choice == "2":
			print "You have picked the Bowtie Aligner to continue.\n";
			alignment_algo = "bowtie";

	#build reference index
	print "To begin Genome Assembly you will need an Reference Genome Index for the specified Alignment Algorithm";

	build_index = raw_input("\nWould you like to build an Index ?\n[1] Yes\n[2] No\n");

	if build_index == "1":
		if alignment_algo == "bwa":
			print "Preparing to Build Alignment Index of Reference Genome for BWA\n";

		else:
			if alignment_algo == "bowtie":
				print "Preparing to Build Alignment Index of Reference Genome for Bowtie\n";
	else:
		if build_index == "2":
			print "You have decided not to build an Alignment Index.\n";
			print "Please make sure you have an Alignment Index as it is necessary for the upcoming processes.\n";

	if build_index == "1":
		reference = raw_input("\nPlease drag & drop Reference Genome File (recommend using Broad Inst. UCSC hg19):\n");
		os.environ['reference'] = reference;

		if alignment_algo == "bwa":

			print "Warning this process will take upto several hours...warm the kettle...\n";
			os.system("$BASEDIR/bwa index -a bwtsw $reference");

		else:
			if alignment_algo == "bowtie":

				print "Warning this process will take upto several hours...warm the kettle...\n";
				os.system("$BASEDIR/bowtie-build $reference hg19");
	# pick and view reference
	reference_dir = raw_input("\nPlease drag & drop your Alinger Indexed Reference File, making sure the Index files are in the same Directory as the Reference File: \n");
	os.environ['reference_dir'] = reference_dir;

	input_file = open(reference_dir, 'r')
	print "\nLet's take a peak inside:\n";

	print tail(input_file, window=10);print "\nLooks okay to me.\n...\n";

	#begin assembly
	if alignment_algo == "bwa":
		print "To begin pair-end read alignment with BWA please drag & drop two FASTQ files\n";

		read1 = raw_input("\nFile 1:\n");
		os.environ['read1'] = read1;

		read2 = raw_input("\nFile 2:\n");
		os.environ['read2'] = read2;

		os.system("$BASEDIR/bwa aln $reference_dir $read1 > $read1.sai");
		os.system("$BASEDIR/bwa aln $reference_dir $read2 > $read2.sai");

		os.system("$BASEDIR/bwa sampe $reference_dir $read1.sai $read2.sai $read1.fastq $read2.fastq > $read1.sam");
		os.system("$BASEDIR/samtools view -bS -o $read1.bam $read1.sam");

		os.system("echo Your final aligned file for this pair of raw reads can be found under $read1.bam");

	else:
		if alignment_algo == "bowtie":
			read = raw_input("To begin single-end read alignment with Bowtie please drag & drop a FASTQ file\n");
			os.environ['read'] = read;

			os.system("$BASEDIR/bowtie -S $reference_dir $read $read.sam");
			os.system("$BASEDIR/samtools view -bS -o $read.bam $read.sam");
			os.system("$BASEDIR/samtools sort $read.bam $read.sorted");

			os.system("echo Your final aligned file for this raw read can be found under $read.sorted.bam");

	print "\nGenome Assembly complete!\n";


#annotator
def anno():
	print "\nKetu Can add Variant Annotations & Functional Annotations"

	choice = raw_input("\nPlease Select an Option to Continue: \n [1] Variant Annotation\n [2] Functional Annotation\n")

	if choice == "1":
		print "\nYou have picked Variant Annotation, which will assign rsIDs to variants.\n";
		anno_algo = "variantanno";

	else:
		if choice == "2":
			print "You have picked Functional Annotation, which will predict functional effects of variants.\n";
			anno_algo = "functionalanno";

	#file input to annotate
	input_variable = raw_input("\nTo begin any Annotation Please drag & drop a raw VCF file:\n")
	os.environ['input_variable'] = input_variable;

	input_file = open(input_variable, 'r')
	print "\nLet's take a peak inside:\n";
	print tail(input_file, window=3);print "\nLooks okay to me.\n...\n";

	#determine need for reference files
	if anno_algo == "variantanno":
		input_reference = raw_input("Please drag & drop a Reference Genome File:\n")
		os.environ['input_reference'] = input_reference;
		print "Loading Reference Genome from:\n"; print input_reference; print "\nLoading...\n";

		input_vcf = raw_input("Please drag & drop SNP Database VCF File:\n")
		os.environ['input_vcf'] = input_vcf;
		print "Loading Reference Genome from:\n"; print input_vcf; print "\nLoading...\n";

	else:
		if anno_algo == "functionalanno":
			print "Functional Annotation utilizes the SnpEff package\n";

			snpeffchoice = raw_input("Please select a SnpEff Reference library to use: \n [1] HG19\n [2] GRCh37\n")
			if snpeffchoice == "1":
				hg_dl = raw_input("Do you already have HG19 downloaded for SnpEff? \n [y] Yes\n [n] No\n")

				if hg_dl == "n":
					print os.path.dirname(os.path.realpath(__file__));
					os.system("java -Xmx4g -jar $BASEDIR/snpEff.jar download hg19 -c $BASEDIR/snpEff.config");
				else:
					if hg_dl == "y":
						print "\nLoading HG19 Database...\n";
			else:
				if snpeffchoice == "2":
					grch_dl = raw_input("Do you already have GRCh37 downaded for SnpEff? \n [y] Yes\n [n] No\n")

					if grch_dl == "n":
						print os.path.dirname(os.path.realpath(__file__));
						os.system("java -Xmx4g -jar $BASEDIR/snpEff.jar download GRCh37.64 -c $BASEDIR/snpEff.config");
					else:
						if grch_dl == "y":
							print "\nLoading GRCh37 Database...\n";

	print "\nInitiating java modules...\nAll parameters are GO\nRunning Annotator...\n";

	#executing annotation
	if anno_algo == "variantanno":
		os.system("java -Xmx4g -jar $BASEDIR/GenomeAnalysisTK.jar -T VariantAnnotator -R $input_reference --variant $input_variable --dbsnp $input_vcf -L $input_variable --alwaysAppendDbsnpId -o $input_variable.anno.vcf");

	else:
		if anno_algo == "functionalanno":
			if snpeffchoice == "1":
				os.system("java -Xmx4g -jar $BASEDIR/snpEff.jar eff -v -onlyCoding true -i vcf -o vcf hg19 $input_variable -c $BASEDIR/snpEff.config > $input_variable.anno.vcf");
			else:
				if snpeffchoice == "2":
					os.system("java -Xmx4g -jar $BASEDIR/snpEff.jar eff -v -onlyCoding true -i vcf -o vcf GRCh37.64 $input_variable -c $BASEDIR/snpEff.config > $input_variable.anno.vcf");

	print "Annotation complete...\n";
	os.system("echo Your new Annotated file can be found at $input_variable.anno.vcf");

#filtration
def filt():
	print "\nKetu Can help create custom filters for Genomic Variants v0.6.28"

	input_variable = raw_input("\nTo begin Filtering variants Please drag & drop a VCF file:\n")
	os.environ['input_variable'] = input_variable;

	input_file = open(input_variable, 'r')
	print "\nLet's take a peak inside:\n";

	print tail(input_file, window=3);print "\nLooks okay to me.\n...\n";

	print "\nInitiating Filtering modules..\n...\n";

	input_reference = raw_input("Please drag & drop a Reference Genome File To Continue Variant Detection:\n")
	os.environ['input_reference'] = input_reference;

	os.system("echo Loading Reference Genome from: $input_reference");

	message = """\
Loading...

Now we're ready to begin Filtering, as a starting point, filtering by Depth of Read is quite helpful.
\n"""

	print message;

	low_vcf = raw_input("What is the Lowest Depth of Read you'd like to be visible?\n");
	os.environ['low_vcf'] = low_vcf;
	os.system("echo Tagging all variants with Depth of Read below $low_vcf as 'Low'");

	message = """\
\nKetu allows for 4 levels of stringencies when it comes to filtering,
please enter the Deapth of Read for the following filter stringencies-
\n"""

	print message;

	mid_vcf = raw_input("Mid: \n");
	os.environ['mid_vcf'] = mid_vcf;
	high_vcf = raw_input("High:\n");
	os.environ['high_vcf'] = high_vcf;
	final_vcf = raw_input("Final Filter:\n");
	os.environ['final_vcf'] = final_vcf;

	print "\nThe Final filter allows you to assign it a custom name\n";
	final_name = raw_input("What would you like to call it? (i.e. 'Xtreme')\n");
	os.environ['final_name'] = final_name;

	os.system("""
        echo Review Depth of Read Filter Assignments:
        echo Low below $low_vcf
        echo Mid up to $mid_vcf
        echo High up to $high_vcf
        echo $final_name above $final_vcf
        """);
	print "\nLoading Filter Stringencies...\nApplying Filters...\n";

	os.system("""java -Xmx4g -jar $BASEDIR/GenomeAnalysisTK.jar -T VariantFiltration -R $input_reference --variant $input_variable -o $input_variable.filt.vcf\
        --filterExpression " DP < $low_vcf " --filterName "Low" \
        --filterExpression " DP > $low_vcf && DP < $mid_vcf " --filterName "Mid" \
        --filterExpression " DP > $mid_vcf && DP < $high_vcf " --filterName "High" \
        --filterExpression " DP > $final_vcf " --filterName "$final_name"\
        """
		  );

	print "\nFiltering complete...\n";
	os.system("echo Your new Filtered file is located at $input_variable.filt.vcf");

	print "\nWould you like to search through your Variants according to the new Filters?";

	message = """\
\nWhich Filter Would you like to be visibile ?
        [1] Low
        [2] Mid
        [3] High
        [4] """;os.system("echo $final_name");

	choice = raw_input("\n [1] Yes\n [2] No\n");

	if choice == "1":
		print message;os.system("echo $final_name");
		visible = raw_input();
		if visible == "1":
			grepme = "'Low'";
		else:
			if visible == "2":
				grepme = "'Mid'";
			else:
				if visible == "3":
					grepme = "'High'";
				else:
					if visible == "4":
						grepme = final_name;

	os.environ['grepme'] = grepme;

	print "\nIt helps to look at a single Chromosome worth of variants at a time.\n";
	chrom = raw_input("\nWhich Chromosome would you like to view?\n");
	os.environ['chrom'] = chrom;

	os.system("echo Searching Chromosome $chrom and filtering by $grepme...");

	os.system("grep %s $input_variable.filt.vcf > del.tmp" % (grepme));
	os.system("grep chr%s del.tmp > del2.tmp" % (chrom));
	os.system("less del2.tmp");
	os.system("rm del2.tmp");
	os.system("rm del.tmp");



#master-controller
print "Please select an option (i.e. '3' for Variant Discovery):\n"
choice_msg="""\

[1] Genome Alignment - From FASTQ short-reads to BAM files
  - Using either BWA or Bowtie algortihms

[2] Process/Recalibrate - Adjust assembeled BAM files according to best practices
  - Header adjustments, mate-pair fixes, etc

[3] Variant Discvoery, Annotation & Filtration

"""
choice = raw_input(choice_msg);

if choice == "1":
        print "\nYou have picked [1] Genome Alignment\n";
        input_option = "assembly";

else:
        if choice == "2":
                print "\nYou have picked [2] Process/Recalibraten\n";
                input_option = "process";
        else:
            if choice =="3":
                print "\nYou have picked [3] Variant Discvoery, Annotation & Filtration\n";
                input_option = "variant";

if input_option == "assembly":
	asmb();
else:
    if input_option == "process":
        proc();
    else:
        if input_option == "variant":
            variant_msg="""\
Initiating Variant-related Moldules..
        ...

        Please select an option:
        [1] Discvoery
        [2] Annotation
        [3] Filtration
        [4] Disease Liklihood
        """
            variant_menu = raw_input(variant_msg);
            if variant_menu == "1":
                print "\nYou have decided to Discover Variants\n";
                disc();
            else:
                if variant_menu == "2":
                    print "\nYou have decided to Annotate Variants\n";
                    anno();
                else:
                    if variant_menu == "3":
                        print "\nYou have decided to Filter Variants\n";
                        filt();
		    else:
			if variant_menu == "4":
			    print "\nYou have decided to calculate Disease Likelihood\n";
			    os.system("python DeadNode_ketu.py");

print "\nThank you for using Ketu Computational Genomics Toolkit, goodbye :)\n";
