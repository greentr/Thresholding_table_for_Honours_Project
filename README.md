Sens_PPV_table README

This is a python script for creating a tsv table listing the sensitivity and positive predictive values resulting from different levels of filtering for depth and proportion thresholds using sample DNA and vcfs created from UNDR ROVER.


Usage:

python Sens_PPV_table.py -c <folder containing sample vcfs> -t <link to vcf containing known true variants>


Output:

Folders within condition folder for each filtering condition

Filtered vcfs within above folders

tsv file containing a tabulated list of sensitivities and positive predictive values for all filtering conditions


Currently 242 filtering conditions are used.

Depth thresholds of 0 to 50 with intervals of 5, 50 to 100 with intervals of 10 and 100 to 400 with intervals of 50.

Proportion thresholds of 0 to 50 with intervals of 5

This can be altered by changing the numpy aranges in the python script

Looks like:

pct_options = numpy.arange(0, 51, 5)

np1 = numpy.arange(0, 50, 5)

np2 = numpy.arange(50, 100, 10)

np3 = numpy.arange(100, 401, 50)

np_options = numpy.concatenate((np1, np2, np3))


pct is proportion threshold

np is depth threshold


Python 2.7 modules used:

cyvcf2 0.7.2

numpy 1.13.1

pybedtools 0.7.10

pysam 0.9.0
