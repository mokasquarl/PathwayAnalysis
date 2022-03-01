#!/usr/bin/python

"""
Definitive Desired Diagnosis (d^3)
from KEGG mapped NGS variants ~go mokas!
"""

import re
import os
import sys
import subprocess
import math
#import compileall

#setting current working directory
BASEDIR = os.path.dirname(os.path.realpath(__file__));
os.chdir(BASEDIR);
os.environ['BASEDIR'] = os.path.dirname(os.path.realpath(__file__));

intro="""\
#
#           _/  _/_/_/
#      _/_/_/        _/
#   _/    _/    _/_/
#  _/    _/        _/
#   _/_/_/  _/_/_/
#
#

Definitive Desired Diagnosis (d^3)
"""
print intro;
#all pathways present in a pedigree multi-vcf
pedi_hid=""
pedi_hid=open("hid.ded","r");
n=0

input_vcf = raw_input("\nPlease input VCF file for trio:")
os.environ['input_vcf'] = input_vcf;

pre_file=""
pre_file=open(input_vcf,"r");

sub_vcf=""
sub_vcf=open("sub_vcf.del","w");

for line in pre_file:
    if not line.startswith("#"):
        elem1=line.split("\t");
        mom= elem1[9];
        dad=elem1[10];
        child=elem1[11];

        elemPLm=elem1[9].split(":");
        elemPLd=elem1[10].split(":");
        elemPLc=elem1[11].split(":");
        momPL=elemPLm[0];
        dadPL=elemPLd[0];
        childPL=elemPLc[0];
        childGT=elemPLc[2].split(",");

        if (childPL == '1/1') and (childPL != momPL) and (childPL != dadPL):
            sub_vcf.write(line);

os.environ['vcf_to_use'] = "sub_vcf.del";

result_file=""
result_file=open("lastDeadNode.results", "w");

#begin scraping KEGG db
for line in bfile:
    elems=line.split(" ");
    print elems[0];
    os.environ['hsaid'] = elems[0];
    os.system("curl http://www.genome.jp/dbget-bin/get_linkdb?-t+disease+path:$hsaid | grep 'ds:H[0-9]*' | cut -d '>' -f 2,3 > 2.del");
   #get Disease ID (HID) from hsa_id list in "hid.ded"

    cfile=open("2.del","r");
    for line in cfile:
        elems2=line.split("</a>               ");
        print elems2[0],elems2[1];
        os.environ['hid'] = elems2[0];
        os.environ['hid1'] = elems2[1];
        os.system("curl http://www.genome.jp/dbget-bin/get_linkdb?-t+genes+disease:$hid | grep 'hsa:' | cut -d '>' -f 3 | cut -d ',' -f 1 | grep -o '\w[A-Z]\w*' | grep -v 'EC' > 3.del");

        os.system("more 3.del | tr '\n' '\|' | sed '$s/.$//' > 4.del");
       # find all genes linked to a given disease "HID"

        dfile=open("4.del","r");
        for line in dfile:
            elems3=line.split();
           # os.system("grep -E %s $input_vcf | grep -v 'PASS' | grep -i -v 'intron' | grep -v 'NONE' | grep -v -c 'SNPEFF_EFFECT=SYNONYMOUS_CODING'" % (elems3[0]));
            cmd="grep -v 'SNPEFF_EFFECT=SYNONYMOUS_CODING' $vcf_to_use | grep -v 'PASS' | grep -i -v 'intron' | grep -v 'NONE' | grep -v 'random' | grep -v 'Un_*' | grep -E '" + elems3[0] + "'";
            #cmd="tee >(grep -v 'SNPEFF_EFFECT=SYNONYMOUS_CODING' $vcf_to_use | grep -v 'PASS' | grep -i -v 'intron' | grep -v 'NONE' | grep 'FQ=-' | grep -E '" + elems3[0] + "'";

            print elems[0],os.system(str(cmd));
            print "--:--:-- --:--:-- --:--:----:--:-- --:--:-- --:--:----:--:-- --:--:-- --:--:----:--:-- --:--:-- --:--:--";
           # clean up

os.system("rm 1.del");
os.system("rm 2.del");
os.system("rm 3.del");
os.system("rm 4.del");

#begin second part of DEADNODE

input_results=""
input_results=open("lastDeadNode.results", "r");
os.environ['input_results'] = input_results;

#taking results from 1st 1/2, extracting and sorting disease IDS
os.system("grep -o 'H[0-9][0-9][0-9][0-9][0-9][0-9]*' $input_results | sort | uniq  > 1.del");

result_count=""
result_count=open("1.del","r");

start_marker = ""
end_marker = '--:--:-- --:--:-- --:--:----:--:-- --:--:-- --:--:----:--:-- --:--:-- --:--:----:--:-- --:--:-- --:--:--'

for line in result_count:
    elem1=line.split();
    start_marker = elem1[0];
    #set HID from 1.del as start marker

    var_count = 0;

    #count the number of variants/disease
    with open(input_results) as inf:
        showline = False
        for line in inf:
            if start_marker in line:
                showline = True
            if showline:
                var_count = var_count+1;
            if end_marker in line:
                showline = False
    if var_count > 5:
        var_count = var_count - 5;


#looking through disease ranked variants determine Functional stats
        disease_id = elem1[0];
        variant_store=""
        variant_store=open("2.del","w");
        start_marker = disease_id
    	mis_count = 0;
        non_count = 0;
        gene_print=""

        with open(input_results) as inf:
            showline = False
            for line in inf:
                if start_marker in line:
                    showline = True
                if end_marker in line:
                    showline = False
                if showline and not line.startswith("H") and not line.startswith("h") and not line.startswith("0") and not line.startswith("\n") and not line.startswith("#") and not line.startswith("256"):
                    elem2=line.split("\t");
                    elem3=elem2[7].split(";");
                    delem_count=""
                    delem_count = 0;
		    print line;
                    #adding up mutation types and printing gene names

                    while delem_count < 15:
                        delem_count = delem_count + 1;

                        if elem3[delem_count] is not None:
                            gene_name = re.compile("SNPEFF_GENE_NAME=[^\s]*");
                            mis_sense = re.compile("SNPEFF_FUNCTIONAL_CLASS=MISSENSE");
                            non_sense = re.compile("SNPEFF_FUNCTIONAL_CLASS=NONSENSE");

                            if gene_name.match(elem3[delem_count]):
                                elem4=elem3[delem_count].split("=");
                                gene_print=elem4[1];
                            if mis_sense.match(elem3[delem_count]):
                                mis_count = mis_count + 1;
                            if non_sense.match(elem3[delem_count]):
                                non_count = non_count + 1;

                            if non_count > 0 and non_sense.match(elem3[delem_count]):
                                print gene_print;

        #scoring diseases
        if mis_count !=0 or non_count !=0:
	    cancer_retin=""
	    cancer_retin = 1;

            cache_file = "kegg.cache";
            start_cache = elem1[0];
            end_cache = "#";
	    score_count=3.14159;
            with open(cache_file) as inf:
                showline = False
                for line in inf:
                    if start_cache in line:
                        showline = True
                    if showline:
                        print line;
			if line.startswith("#"):
			    elem5=line.split("#");
			    matchObj = re.search( r'cancer', line, re.M|re.I)
			    if matchObj:
			        cancer_retin = 0.1;
			    matchObj2 = re.search( r'retin[^\s]*', line, re.M|re.I)
			    if matchObj2:
				cancer_retin = 0.1;
                    if end_cache in line:
                        showline = False
	    print "Nonsense Mutations: ", non_count;
            print "Missense Mutations: ", mis_count;
	    score_count = (non_count * cancer_retin) + (mis_count*0.1*cancer_retin);
	    print "Disease Rank Score: ", score_count;
	    print end_marker;

os.system("rm 1.del");
