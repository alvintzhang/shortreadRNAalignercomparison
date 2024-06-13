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
            if(SoH == 0):
                SoH = 1
        elif op == 5:
            if(SoH == 0):
                SoH = 2
            print("hard clip: ",le)



query_positions = []
reference_positions = []

for read in bamfile.fetch():
    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
            query_positions.append(query_pos)
            reference_positions.append(ref_pos)

print(query_positions) #correct but no splice
print(reference_positions) #correct but no splice




query_positions2 = []
reference_positions2 = []

splices = [] #holds last match and first match (beginning/end of splice)
#splices holds positions of the query column/row of the table

queryindex = 0
referenceindex = 0

for read in bamfile.fetch():

    if(SoH==1):
        queryindex = read.query_alignment_start
    elif(SoH==2):
        queryindex = 0
    else:
        print("No clipping Error")
        sys.exit(1)
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
        elif(op == 2):
            referenceindex+=le
        elif(op == 3):
            splices.append([queryindex-1,queryindex])
            referenceindex+=le



dict = {}
print(query_positions2)
print(reference_positions2)
print(splices)
assert len(query_positions2)==len(reference_positions2)
for i in range(len(query_positions2)):
    dict[query_positions2[i]] = reference_positions2[i]





for x in splices:
    print("{}:{},{}:{}".format(x[0],dict[x[0]],x[1], dict[x[1]])) #printing out the splices and ref positions
bamfile.close()



