#!/usr/bin/env python
# coding: utf-8

#initialize module
import sys, getopt
import math
import numpy as np
import gzip

#initialize variables
#matrix of cells x regions
cell_ids = {}
#dictionary of cell barcodes with total counts per chr
cell_counts = {}
chrom_step = 1000000
chrom_max = 0
file_init = ''
line1 = ''
min_filter = 10000
output_file = ''
#chrom sizes file
chrom_sizefile = ''
#initialization
try:
   opts, args = getopt.getopt(sys.argv[1:],"hi:o:b:f:g:",["help","ifile=","ofile=","binsize=","frags=","genome="])
   #opts, args = getopt.getopt(sys.argv[1:],"h",["help"])
   print("process_fragment_file",len(opts))
except getopt.GetoptError as err:
   print('process_fragment_file.py -i <inputfile.tsv.gz> -o <outputfile.tsv> -b <binsize> -f <minfrags> -g <genomefile.chrom.sizes>')
   print(str(err),"process_fragment_file")
   sys.exit(2)
for opt, arg in opts:
   print("process_fragment_file",opt)
   if opt in ("-h","help"):
    print('process_fragment_file.py  -i <inputfile.tsv.gz> -o <outputfile.tsv> -b <binsize> -f <minfrags> -g <genomefile>')
    sys.exit()
   elif opt in ("-i", "--ifile"):
     file_init = arg
   elif opt in ("-o", "--ofile"):
     output_file = arg
   elif opt in ("-b", "--binsize"):
     chrom_step = int(arg)
   elif opt in ("-f", "--frags"):
     min_filter = int(arg)
   elif opt in ("-g", "--genome"):
     chrom_sizefile = arg

if len(opts)<5:
   print('process_fragment_file.py  -i <inputfile.tsv.gz> -o <outputfile.tsv> -b <binsize> -f <minfrags> -g <genomefile>')
   sys.exit(2)

#load sizes
chrom_sizes = {}
def init_chrom_sizes():
    chrom_size = {}
    with open(chrom_sizefile) as reader:
        line1 = reader.readline().rstrip()
        while line1 != '':
          line_split = line1.split('\t')
          chrom_size[line_split[0]]=line_split[1]
          line1 = reader.readline().rstrip()
    return chrom_size

chrom_sizes = init_chrom_sizes()
print(chrom_sizes)
print("Initializing chromosomes")
#initialize chromosome bins
def init_chrom_list_all():
    chrom_min = 0
    chrom_total_list = {}
    for chrom in chrom_sizes:
        #chrom_list_temp = {}
        chrom_max1 = int(chrom_sizes[chrom])
        chrom_list_temp = np.zeros(int(chrom_max1 / chrom_step)+1)
#            chrom_list_temp[i*chrom_step] = 0
        chrom_total_list[chrom] = chrom_list_temp
    return chrom_total_list

blank_list_all = init_chrom_list_all()


print("Reading file")
#read file
lines_read = 0
with gzip.open(file_init,'rt') as reader:
    for line1 in reader:
      #print(line1 + " " + str(len(cell_ids)))
      #remove header from new version of cellranger-atac
      if (line1.find('#')!=-1):
        print("commented header line detected...skipping")
        continue 
      line_split = line1.split('\t')
      readLength=int(line_split[2])-int(line_split[1])-1
      #test if chromosome is in the list
      if (line_split[0] in chrom_sizes.keys()):
        #print(line_split[0] + "\n")
        if (line_split[3] in cell_ids):
          setIndex = int(math.floor(int(line_split[1]) / chrom_step))
          cell_ids[line_split[3]][line_split[0]][setIndex] += readLength
          cell_counts[line_split[3]] += 1
        else:
          cell_ids[line_split[3]]=init_chrom_list_all()
          setIndex = int(math.floor(int(line_split[1]) / chrom_step))
          cell_ids[line_split[3]][line_split[0]][setIndex] += readLength
          cell_counts[line_split[3]]=1
      else:
        #abnormal chromosome - will not include in final results
        print("warning: chromosome " + line_split[0] + " found in fragment file but not in chromosome file - skipping fragment")
      lines_read+=1
      if lines_read%100000==0:
        print("Reading line " + str(lines_read) + "\n")
print("File loaded into memory")


#how many pass filter    
filter_passing = 0

#minimum number of fragments per barcode - quality filter

print("Writing output")
with open(output_file,'w') as writer:
    #initialize header
    writer.write("Cell_id")
    for chrom in chrom_sizes:
         for i in range(0,len(blank_list_all[chrom])):
              writer.write("\t" + chrom + "_" + str(i*chrom_step))
    writer.write("\n")
    for cell in cell_ids:
        if (cell_counts[cell] > min_filter):
           chr_list = cell_ids[cell]
           filter_passing += 1
           writer.write(cell)
           for chrom in chrom_sizes:
               writer.write('\t')
               writer.write('\t'.join(map(str, chr_list[chrom])))
                 
           writer.write("\n")
    
print("Done writing output")


