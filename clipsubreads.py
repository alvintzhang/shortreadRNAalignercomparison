#Author: Alvin Zhang
#Last changes: 6/13/24
#Parsing CIGAR string + creating query/reference position table (two arrays) that also keep track of splices

import pysam
import sys


bamfile = pysam.AlignmentFile("/Users/AlvinZhang2026/bio_data/several_reads2.bam", "rb")
seq = ""
SoH = 0
cigar = ""

#just to get a basic read

for read in bamfile.fetch():
    print("read name: ",read.query_name)
    print("reference name: ",read.reference_name)
    print("first position on ref: ",read.reference_start)
    print("first position on read: ",read.query_alignment_start)
    print("read sequence: ",read.query_sequence)

    seq = read.query_sequence
    cigar = read.cigarstring

    print("cigar: ",read.cigarstring)

    for (op, le) in read.cigartuples:
        if op == 0 or op == 7 or op == 8:
            print("match: ",le)
        elif op == 1:
            print("insertion: ",le)
        elif op == 2:
            print("deletion: ",le)
        elif op == 3:
            print("intron: ",le)
        elif op == 4:
            print("soft clip: ",le)
            SoH = 0
        elif op == 5:
            SoH = 1
            print("hard clip: ",le)

#print(seq)
#print(SoH)

query_positions = []
reference_positions = []

for read in bamfile.fetch():
    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
            query_positions.append(query_pos)
            reference_positions.append(ref_pos)

print(query_positions) #correct but no splice
print(reference_positions) #correct but no splice

#splicing and deletion are the same ~ you can't tell

#print(SoH)









query_positions2 = []
reference_positions2 = []

splices = [] #holds last match and first match (beginning/end of splice)
#splices holds positions of the query column/row of the table

queryindex = 0
referenceindex = 0

#for base in read.query_sequence:
    #query_positions2.append(read.query_alignment_start+queryindex)
    #reference_positions2.append(read.reference_start+referenceindex)

for read in bamfile.fetch():

    if(SoH==0):
        queryindex = read.query_alignment_start
    else:
        queryindex = 0
    referenceindex = read.reference_start

    for op, le in read.cigartuples:
        if(op == 0):
            for i in range(le):
                query_positions2.append(queryindex)
                reference_positions2.append(referenceindex)
                queryindex+=1
                referenceindex+=1
        elif(op == 1):
            queryindex+=le
            query_positions2.append(queryindex)
        elif(op == 2):
            referenceindex+=le
            reference_positions2.append(referenceindex)
        elif(op == 3):
            splices.append([queryindex,queryindex+1])
            for x in range(2):
                query_positions2.append(queryindex)
                queryindex+=1
            reference_positions2.append(referenceindex)
            referenceindex+=le
            reference_positions2.append(referenceindex)











print(query_positions2)
print(reference_positions2)
print(splices)

bamfile.close()



