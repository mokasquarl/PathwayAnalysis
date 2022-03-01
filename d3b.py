#!/usr/bin/python

"""
Definitive Desired Diagnosis (d^3)
no REST calls from KEGG mapped NGS variants ~go mokas!
"""

import re
import os
import sys
import subprocess
import math
#import compileall
import pdb
import fileinput

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
n=0

input_vcf = "/Users/mokas/Desktop/codebase/max/working/KEGG/Good_Trio_preKG.vcf";
os.environ['input_vcf'] = input_vcf;

pre_file=""
pre_file=open(input_vcf,"r");

f=open('VB_polyphen_query.txt','w+');

with open('/Users/mokas/Desktop/codebase/max/working/Eff/M_XY_hg19_Eff_anno.vcf', mode='r') as polyphen:
	for line in polyphen:
		if not line.startswith("#"):
			element=line.split("\t");
			chrosom=str(element[0]).rstrip();
			positon=element[1];
			referen=str(element[3]).rstrip();
			alterne=str(element[4]).rstrip();
			f.write(chrosom+':'+positon+' '+referen+'/'+alterne+'\n');
	f.close();
pdb.set_trace();

#sub_vcf=""
#sub_vcf=open("sub_vcf.del","r");

pos_to_ccds=""
pos_to_ccds=open("ccdsKgMap.txt","r");

ccds_to_hsa=""
ccds_to_hsa=open("keggPathway.txt","r");

uc_id=""
os.environ['uc_id'] = uc_id;

hsa_name=""
os.environ['hsa_name'] = hsa_name;
vcf_fmt=""
def writeAgain():
    f = open('tmpFmtln.del', 'r');
    lines = f.readlines();
    f.close();
    f = open('tmpFmtln.del', 'w');
    for line in lines:
        if line != vcf_fmt:
            f.write(line);
    f.close();

fmt_file=open('tmpFmtln.del', 'w+');
f=open('tmpCP.del', 'w+');

for line in ccds_to_hsa:
    element=line.split("\t");
    uc_id = element[0];
    os.environ['uc_id'] = uc_id;
    duh_result=os.system("grep $uc_id ./ccdsKgMap.txt > tmpCCDS.del");
    doh_result=os.system("cat tmpCCDS.del | cut -f3,4,5 > tmpChrPos.del");
    doh_chr=os.popen('cat tmpChrPos.del | cut -f1').read();
    print ('this CCDS' + doh_chr);
    doh_poA=os.popen('cat tmpChrPos.del | cut -f2').read();
    doh_poB=os.popen('cat tmpChrPos.del | cut -f3').read();
    chr_doh=str(doh_chr).rstrip();
    try:
        poA_doh=int(doh_poA);
        poB_doh=int(doh_poB);
    except ValueError:
        pass

    hsa_name=str(element[2]).rstrip();
    os.environ['hsa_name'] = hsa_name;
    hsa_grep=os.system("grep $hsa_name ./keggMapDesc.txt > tmpDesc.del");
    hsa_desc=os.popen('cat tmpDesc.del | cut -f2').read();
    desc_hsa=str(hsa_desc).replace(" ","").rstrip();
    os.environ['desc_hsa'] = desc_hsa;
    print ('with '+chr_doh);

    with open('/Users/mokas/Desktop/codebase/max/working/KEGG/Good_Trio_preKG.vcf', mode='r') as inf:
        for line in inf:
            if not line.startswith("#"):
                elemo=line.split("\t");
                chro=str(elemo[0]).rstrip();
                poso=int(elemo[1]);
                if chr_doh == chro and poso >= poA_doh and poso <= poB_doh:
                    vcf_line=line;
                    os.environ['vcf_line'] = vcf_line;
                    if fmt_file.tell() == 0:
                        fmt_file.write(str(vcf_line).replace(";VDB", ";KP="+hsa_name+":"+desc_hsa+";VDB"));
                    elif fmt_file.tell() != 0:
                        prevFile=str(vcf_line).replace(";VDB", ";KP="+hsa_name+":"+desc_hsa+";VDB");
                        subElem=prevFile.split("\t");
                        print subElem[0];
                        subChr=str(subElem[0]).rstrip();
                        try:
                            subPos=int(subElem[1]);
                        except ValueError:
                            pass
                        os.environ['subChr'] = subChr;
                        os.environ['subPos'] = str(subPos);
                        isPresent=os.popen('grep $subChr ./tmpFmtln.del | grep -c $subPos').read();
                        try:
                            presentNum=int(isPresent);
                        except ValueError:
                            pass
                        print hsa_name;
                        print presentNum;
                        if presentNum == 0:
                            fmt_file.write(str(vcf_line).replace(";VDB", ";KP="+hsa_name+":"+desc_hsa+";VDB"));
                        elif presentNum != 0:
                            lineMod = os.popen('grep $subChr ./tmpFmtln.del | grep $subPos').read();
                            os.environ['lineMod'] = str(lineMod);
                            os.popen("grep -v '$lineMod' ./tmpFmtln.del>./tmpCP.del");
                            f.write(lineMod);
                            os.popen("cp ./tmpCP.del ./tmpFmtln.del");
#                    for line in fileinput.input("tmpFmtln.del", inplace = 1):
#                        sub_elem=line.split("\t");
#                        sub_chr=str(sub_elem[0]).rstrip();
#                        sub_pos=int(sub_elem[1]);
#                        if chr_doh == sub_chr and poso == sub_pos:
#                            print line.replace(";VDB", "|"+hsa_name+":"+desc_hsa+";VDB"),
#                        if line.find(vcf_line) == -1 and chr_doh != sub_chr and poso != sub_pos:
#                            fmt_file.write(str(vcf_line).replace(";VDB", ";KP="+hsa_name+":"+desc_hsa+";VDB"));

                    if chr_doh == 'chr18':
                        pdb.set_trace();
#                    os.system("echo $fmt_line >> tmpFmtln.del");

for line in pre_file:
    if not line.startswith("#"):
        elem1=line.split("\t");
        chrom= elem1[0];
        pos=elem1[1];
        print chrom;
        print pos;

os.environ['vcf_to_use'] = "sub_vcf.del";

result_file=""
result_file=open("lastDeadNode.results", "w");

bfile=open("hid.ded","r");
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
