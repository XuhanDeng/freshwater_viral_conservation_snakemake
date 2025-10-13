import os
import pandas
import argparse


import pandas as pd

parser = argparse.ArgumentParser(description="virsoter2,max_score>=0.8;dvf:score>=0.8,p<0.05;genomad != provirus",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i1","--dvf",help="deepvirfinder_outputtable;dvfpred.txt")
parser.add_argument("-i2","--vs2",help="virsorter_result;final-viral-score.tsv")
parser.add_argument("-i3","--genomad",help="virsorter_result;_virus_summary.tsv")
parser.add_argument("-o","--output",help="the out put dir you want store file")
args = parser.parse_args()
dvf_in = pd.read_csv(args.dvf,sep="\t")

vs2_in = pd.read_csv(args.vs2,sep="\t")

genomad_in = pd.read_csv(args.genomad,sep="\t")


#dvf_filter
dvf_in["tag1"] = "dvf"
dvf_filter = dvf_in[(dvf_in["score"]>=0.8) & (dvf_in["pvalue"]<0.05)]
dvf_list = dvf_filter.loc[:,["name","tag1"]]

#vs2
vs2_in["tag2"] = "vs2"
vs2_in["name"] = vs2_in["seqname"].str.split("|").str[0]
vs2_in["full/partial"] = vs2_in["seqname"].str.split("|").str[2]

vs2_filter = vs2_in[(vs2_in["max_score"]>=0.8) & (vs2_in["full/partial"]=="full")]

vs2_list = vs2_filter.loc[:,["name","tag2"]]
#genomad
genomad_in["tag3"] = "genomad"
genomad_in["name"] = genomad_in["seq_name"]
genomad_filter = genomad_in[genomad_in["topology"]!="Provirus"]
genomad_list=genomad_filter.loc[:,["name","tag3"]]


##merge_index_table
merge_index = pd.DataFrame()

merge_index["name"]=pd.concat([dvf_list.loc[:,"name"],genomad_list.loc[:,"name"],vs2_list.loc[:,"name"]],ignore_index=True)

merge_unique = merge_index.drop_duplicates("name")
# 第一次合并 merge_unique 和 dvf_list，基于 'name' 列，使用左连接
merge_temp = pd.merge(merge_unique, dvf_list, how="left", on="name")

# 第二次合并 merge_temp 和 genomad_list，基于 'name' 列，使用左连接
merge_temp = pd.merge(merge_temp, genomad_list, how="left", on="name")

# 第三次合并 merge_temp 和 vs2_list，基于 'name' 列，使用左连接
merge_add_label = pd.merge(merge_temp, vs2_list, how="left", on="name")

extract_list = merge_add_label.loc[:,"name"]

base_name = os.path.basename(args.dvf).split("_scaffolds_")[0]

# Create subdirectories
csv_dir = os.path.join(args.output, "csv")
list_dir = os.path.join(args.output, "list")
os.makedirs(csv_dir, exist_ok=True)
os.makedirs(list_dir, exist_ok=True)

# Save outputs
merge_add_label.to_csv(os.path.join(csv_dir, base_name+"_merged_results.csv"), index=False)
extract_list.to_csv(os.path.join(list_dir, base_name+"_merge3_list.txt"), index=False, header=None)