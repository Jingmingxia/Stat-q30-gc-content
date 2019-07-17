#!/usr/bin/env python
#encoding: utf-8
#by jmxia
import sys,os

################################################
# The script is to stat the q20, q30, gc_content and output fasta 
# Usage: python $0 input_fq output_fa
###############################################

def stat_qual(qstr):
    q20 = 0
    q30 = 0
    for q in qstr:
        qual = ord(q) - 33   # ord function return the value of the ASCII
        if qual >= 30:
            q30 += 1
            q20 += 1
        elif qual >= 20:
            q20 += 1
    return q20, q30

def read_fq(files):
  dict_seq = {}
  dict_qual = {}
  n = 0
  for line in open(files,'r'):
    line = line.strip()
    if not line:
      continue
    n += 1
    if line.startswith("@") and n != 4:
      lines = line.strip()
      scaf = lines
      seq = ""
      qual = ""
      n = 1
    elif n == 2:
      seq = line
      dict_seq[scaf] = seq
    elif n == 3:
      plus = line
    elif n == 4:
      qual = line
      dict_qual[scaf] = qual
  return dict_seq,dict_qual

fq_file = sys.argv[1]
out_fasta = sys.argv[2]
dict_seq,dict_qual = read_fq(fq_file)

a,t,c,g,n,dna = 0,0,0,0,0,0
with open(out_fasta,'w') as fh:
  for key in dict_seq:
    k = key.replace("@",">")
    fh.write("{0}\n{1}\n".format(k,dict_seq[key]))
    a += dict_seq[key].count("A")
    t += dict_seq[key].count("T")
    c += dict_seq[key].count("C")
    g += dict_seq[key].count("G")
    n += dict_seq[key].count("N")
    dna += len(dict_seq[key])

q20_count = 0
q30_count = 0
reads = 0
for key in dict_qual:
  reads += 1
  q20,q30 = stat_qual(dict_qual[key])
  q20_count += q20
  q30_count += q30

ap = round(float(a)*100/float(dna),2)
tp = round(float(t)*100/float(dna),2)
cp = round(float(c)*100/float(dna),2)
gp = round(float(g)*100/float(dna),2)
np = round(float(n)*100/float(dna),2)
gc = g + c
gcp = round(float(gc)*100/float(dna),2)
q30 = round(float(q30_count)*100/float(dna),2)
q20 = round(float(q20_count)*100/float(dna),2)
print("""Base\tNumber\tPercentage
A:\t{0}\t{1}
T:\t{2}\t{3}
C\t{4}\t{5}
G\t{6}\t{7}
N\t{8}\t{9}
GC\t{10}\t{11}
Total\t{12}\t{13}
""".format(a,ap,t,tp,c,cp,g,gp,n,np,gc,gcp,dna,100.00))
print("""Total_reads:{0}\n
Quality_type\tBases\tPercentage
q20:\t{1}\t{2}
q30:\t{3}\t{4}
""".format(reads,q20_count,q20,q30_count,q30))

