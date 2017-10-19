import sys
import os
import argparse
import datetime
from os import listdir
from os.path import isfile, join
import numpy
import pybedtools
from cyvcf2 import VCF, Writer

current_dir = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--condition", help="link to the folder containing the sample vcfs for a given condition")
parser.add_argument("-t", "--truthset", help="define the vcf containing a list of known true variants")
args = parser.parse_args()

if args.condition:
    condition_full = args.condition
    condition = os.path.basename(os.path.normpath(condition_full))
else:
    print "A condition must be provided using '-c' or '--condition' followed by a link to the folder containing the sample vcfs for a given condition"

# check condition exists
if os.path.exists(condition_full):
    print "Condition is %s" % condition
    pass
else:
    print "Error: Condition does not exist"
    quit()

# define truth set
if args.truthset:
    truth_set = args.truthset
else:
    print "A truth set must be provided using '-t' or '--truthset' followed by a link to the vcf contain a list of true variants"
    quit()

truthsetcyvcf2 = VCF(truth_set)
truthsetnumber = 0
for truthsetvariant in truthsetcyvcf2:
    truthsetnumber += 1
print ('Truth set has %s variants') % truthsetnumber

a = pybedtools.BedTool(truth_set)

# define samples
og_vcfs = [files for files in listdir(condition_full) if (isfile(join(condition_full, files)) and (files.endswith('.vcf')))]
sample_number = len(og_vcfs)
print ("There are %s samples") % sample_number

completed_number = 0

# make tsv file for table
if not os.path.exists(current_dir + '/' + 'Sens_PPV_table_outputs'):
            os.makedirs(current_dir + '/' + 'Sens_PPV_table_outputs')
now = datetime.datetime.now()
f = open("Sens_PPV_table_outputs/Output_%s_%s.tsv" % (condition, now.strftime("%Y-%m-%d_%H:%M")), "w")
f.write("Proportion \t")        
f.write("Depth \t")
f.write("Sensitivity \t")
f.write("PPV \n")
f.close()

# define filtering conditions
pct_options = numpy.arange(0, 51, 5)
np1 = numpy.arange(0, 50, 5)
np2 = numpy.arange(50, 100, 10)
np3 = numpy.arange(100, 401, 50)
np_options = numpy.concatenate((np1, np2, np3))
total_filter_number = (len(pct_options) * len(np_options))
print "Total number of filtering conditions is %s" % total_filter_number

# make progess meter
progress = round(((float(completed_number)/float(total_filter_number))*100),2)
print "\r", "%s percent complete |" % progress,

# begin loop for each condition
for pct in pct_options:
    for np in np_options:
        os.chdir(current_dir)
        
        # make folder for filtered vcfs
        if not os.path.exists(current_dir + '/' + condition + '/' + 'filtered_NP' + str(np) + '_PCT' + str(pct)):
            os.makedirs(current_dir + '/' + condition + '/' + 'filtered_NP' + str(np) + '_PCT' + str(pct))
        mypath = current_dir + '/' + condition + '/' + 'filtered_NP' + str(np) + '_PCT' + str(pct)
        
        # filter for depth and proportion thresholds as well as removing variants next to homopolymeric runs
        for vcf in og_vcfs:
            os.chdir(condition_full)
            cyvcf_vcf = VCF(vcf)
            new_filtered_vcf = mypath + '/' + "filt_%s" % vcf
            w = Writer(new_filtered_vcf, cyvcf_vcf)
            for v in cyvcf_vcf:
                if v.INFO["NP"] >= np and v.INFO["PCT"] >= pct:
                    try:
                        if v.INFO["HRUN"] is 'A' or 'C' or 'G' or 'T':
                            pass
                        else:
                            pass
                    except KeyError:
                        w.write_record(v)
            w.close()
        
        print "\r", "%s percent complete /" % progress, # update progress meter
        
        # define filtered vcfs
        variant_file_list = [files for files in listdir(mypath) if (isfile(join(mypath, files)) and (files.endswith('.vcf')))]  
        os.chdir(current_dir)
        
        # make metrics dictionary
        metrics = {}
        for vcf in variant_file_list:
            key = os.path.basename(vcf)
            metrics[key] = {}
        
        # intersect between truth set and filtered vcfs and add count for true positives
        for vcf in variant_file_list:
            print "\r", "%s percent complete -" % progress, #update progress meter
            os.chdir(mypath)
            b = pybedtools.BedTool(vcf)
            a_and_b = b.intersect(a, c=True, header=True, stream=True)
            true_positives = 0
            false_positives = 0
            
            # count the numbers to determine true positive and false positive count
            for line in a_and_b:
                if line[8] == '0':
                    false_positives += 1
                if line[8] == '1':
                    true_positives += 1
            print "\r", "%s percent complete \\" % progress, # update progress meter
            
            # account for bedtools' errors in intersecting
            if true_positives > truthsetnumber:
                diff = true_positives - truthsetnumber
                false_positives += diff
                true_positives -= diff
                
            #calculate other statistical values
            false_negatives = truthsetnumber - true_positives
            total_positives = true_positives + false_negatives
            TP = float(true_positives)
            P = float(total_positives)
            FP = float(false_positives)
            sensitivity = TP / P
            try:
                PPV = TP / (TP + FP)
            except ZeroDivisionError: # fails if no positive calls are made on a single sample
                pass
            
            # add relevant metrics to dictionary
            access_key = os.path.basename(vcf)
            metrics[access_key]['sensitivity'] = sensitivity
            metrics[access_key]['Positive_Predictive_Value'] = PPV

        # prepare metrics for output
        x = sum(metrics[sens]['sensitivity'] for sens in metrics) / len(metrics)
        y = sum(metrics[PosPred]['Positive_Predictive_Value'] for PosPred in metrics) / len(metrics)

        # add new line to table (tsv) containing data from this condition
        os.chdir(current_dir)
        f = open("Sens_PPV_table_outputs/Output_%s_%s.tsv" % (condition, now.strftime("%Y-%m-%d_%H:%M")), "a")
        f.write("%s \t" % pct)        
        f.write("%s \t" % np)
        f.write("%s \t" % x)
        f.write("%s \n" % y)
        f.close()
        
        # increase count for progress meter
        completed_number += 1
        progress = round(((float(completed_number)/float(total_filter_number))*100),2)
        print "\r", "%s percent complete |" % progress, # update progress meter

# announce completion
print "\n"
print "done"