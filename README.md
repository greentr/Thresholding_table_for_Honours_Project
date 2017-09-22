# Filtering-Script
README:
Variant filtering and metric producing script for UNDR ROVER VCF files by Thomas Green for 2017 Honours project.
Made to work with multiple VCFs produced from UNDR ROVER after a sequencing run using a cell line DNA

Required to use:
	Files from github:
		folder named: TP-FP-outputs
		add_sample_count.py made by Khalid Mahmood
		filt.sh
		Sens_PPV.py
		Blank text file named: Sens_PPV-statistics.txt
	Python 2.7
	Python tools:
		cyvcf2
		pybedtools
		pysam
	Command line tools:
		vt
	VCF file containing truth set named TRUTHSET.vcf
	Folder containing VCFs

To run:
	Use command in terminal:
		bash filt.sh <name of folder with VCFs> <min read pairs depth> <min % of read pairs variant appears in>

Produces: 
	Folder within VCF folder that contains VCFs filtered for read pairs and percentage
	VCF files in folder 'TP-FP-outputs' that list true positives and false positives that appear with a sample count of how many samples they appear in
	Information in the text file 'Sens_PPV-statistics.txt' detailing the sensitivity and positive predictive value with the given conditions

# Filtering-Script
