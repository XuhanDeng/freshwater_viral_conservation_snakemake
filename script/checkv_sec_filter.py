import os
import pandas
import argparse


import pandas as pd

parser = argparse.ArgumentParser(description="virsoter2,max_score>=0.8;dvf:score>=0.8,p<0.05;genomad != provirus",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i1","--checkv",help="quality_summary.tsv")
parser.add_argument("-i2","--vs2",help="virsorter_result;final-viral-score.tsv")
parser.add_argument("-i3","--genomad",help="virsorter_result;_virus_summary.tsv")
parser.add_argument("-o","--output",help="the out put dir you want store file")
args = parser.parse_args()
checkv_in = pd.read_csv(args.checkv,sep="\t")

vs2_in = pd.read_csv(args.vs2,sep="\t")

genomad_in = pd.read_csv(args.genomad,sep="\t")


#vs2
vs2_in["name"] = vs2_in["seqname"].str.split("|").str[0]
vs2_in["full/partial"] = vs2_in["seqname"].str.split("|").str[2]

vs2_filter = vs2_in[(vs2_in["hallmark"]>0)& (vs2_in["full/partial"]=="full")]

vs2_list = vs2_filter.loc[:,["name"]]
#genomad
genomad_in["name"] = genomad_in["seq_name"]
genomad_filter = genomad_in[(genomad_in["topology"]!="Provirus") & (genomad_in["n_hallmarks"]>0)]
genomad_list = genomad_filter.loc[:,["name"]]


#checkv

checkv_in["name"]=checkv_in["contig_id"]
checkv_filter = checkv_in[checkv_in.loc[:,"viral_genes"]>0]
checkv_list = checkv_filter.loc[:,["name"]]

##merge_index_table
list_merge = pd.DataFrame()

list_merge["name"] = pd.concat([vs2_list.iloc[:,0],genomad_list.iloc[:,0],checkv_list.iloc[:,0]],ignore_index=True)

list_merge = list_merge.drop_duplicates("name")

base_name = os.path.basename(args.genomad).split("_scaffolds_")[0]

list_merge.to_csv(os.path.join(args.output,base_name+"_checkv_extract.txt"),index=False,header=None)