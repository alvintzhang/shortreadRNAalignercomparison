#Author: Alvin Zhang
#Last changes: 6/13/24
#Parsing CIGAR string + creating query/reference position table (two arrays) that also keep track of splices

import pysam
import sys

bamfile = pysam.AlignmentFile("/Users/AlvinZhang2026/bio_data/several_reads2.bam", "rb")


#parseCigar method which reads the given bam file
#returns the sequence, cigar string, and whether the left clip is soft/hard
def parseCigar(givenBamfile):
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

    return [seq, cigar, SoH]


# Creates two arrays that hold the query position and reference positions of the read
# Cannot differentiate between deletions and splices
def QRtable1(givenBamfile):
    query_positions = []
    reference_positions = []

    for read in givenBamfile.fetch():
        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
            query_positions.append(query_pos)
            reference_positions.append(ref_pos)

    print(query_positions) #correct but no splice
    print(reference_positions) #correct but no splice

QRtable1(bamfile)

# Just like QRtable1 but also returns an array with the splices
# Therefore able to differentiate between all operations returned by cigartuples
def QRtable2(givenBamfile, SoH):
    query_positions2 = []
    reference_positions2 = []

    splices = [] # holds last match and first match (beginning/end of splice)
    # splices holds positions of the query column/row of the table

    queryindex = 0
    referenceindex = 0

    for read in givenBamfile.fetch():

        if(SoH==1):
            queryindex = read.query_alignment_start
        elif(SoH==2):
            queryindex = 0
        else:
            print("No clipping Error")
            sys.exit(1)
        referenceindex = read.reference_start

        for op, le in read.cigartuples:
            if(op == 0): #match or mismatch
                for i in range(le):
                    query_positions2.append(queryindex)
                    reference_positions2.append(referenceindex)
                    queryindex+=1
                    referenceindex+=1
            elif(op == 1): #insertion
                queryindex+=le
            elif(op == 2): #deletion
                referenceindex+=le
            elif(op == 3): #splice
                splices.append([queryindex-1,queryindex])
                referenceindex+=le

    array = [query_positions2, reference_positions2, splices]
    return array

def findRefSplicePos(qp, rp, s):
    dict = {}
    assert len(query_positions2)==len(reference_positions2)
    for i in range(len(qp)):
        dict[qp[i]] = rp[i]



    for x in splices:
        print("{}:{},{}:{}".format(x[0],dict[x[0]],x[1], dict[x[1]])) #printing out the splices and ref positions

    return dict

parseCigarArray = parseCigar(bamfile) #calling parseCigar
SoftorHard = parseCigarArray[2]
QRtable2Array = QRtable2(bamfile, SoftorHard)

query_positions2 = QRtable2Array[0]
reference_positions2 = QRtable2Array[1]
splices = QRtable2Array[2]

print(query_positions2)
print(reference_positions2)
print(splices)

findRefSplicePos(query_positions2, reference_positions2, splices)



bamfile.close()