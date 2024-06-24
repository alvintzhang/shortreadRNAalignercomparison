#Author: Alvin Zhang
#Last changes: 6/14/24
#Processes a BAM file to extract and parse CIGAR strings, create position tables, and identify splices.

import pysam
import sys

bamfile = pysam.AlignmentFile("/Users/AlvinZhang2026/bio_data/several_reads2.bam", "rb")



#parseCigar method which reads the given BAM file
#returns the sequences, CIGAR strings, and whether the left clip is soft or hard (SoH variable)
def parseCigar(givenBamfile):
    seq = [] #all the reads in the BAM file
    SoH = [] #soft or hard of left clip for each read
    cigar = [] #array of strings with CIGAR strings of all the potential reads

    for read in bamfile.fetch(): #if there is more than one read in the file
        if(read.is_secondary):
            continue
        print("read name: ",read.query_name)
        print("reference name: ",read.reference_name)
        print("first position on ref: ",read.reference_start)
        print("first position on read: ",read.query_alignment_start)
        print("read sequence: ",read.query_sequence)
        seq.append(read.query_sequence)
        cigar.append(read.cigarstring)

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
                if not SoH:
                    SoH.append(1)
            elif op == 5:
                if not SoH:
                    SoH.append(2)
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

    array = [query_positions, reference_positions]
    return array

# Just like QRtable1 but also returns an array with the splices
# Therefore able to differentiate between all operations returned by cigartuples
def QRtable2(givenBamfile, SoH):
    query_positions2 = []
    reference_positions2 = []

    splices = [] # holds last match and first match (beginning/end of splice)
    # splices holds positions of the query column/row of the table
    SoHindex = 0
    queryindex = 0
    referenceindex = 0

    for read in givenBamfile.fetch():

        if SoH==1 :
            queryindex = read.query_alignment_start
        elif SoH==2 :
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
    assert len(qp)==len(rp)
    for i in range(len(qp)):
        dict[qp[i]] = rp[i]

    qrDASplicePos = []
    for x in s:
        qrDASplicePos.append("{}:{},{}:{}".format(x[0],dict[x[0]],x[1], dict[x[1]])) #printing out the splices and ref positions

    return [dict, qrDASplicePos]

def extract_subreads(query_seq, query_start):
    # Extract the subread from the query sequence
    subread = query_seq[query_start-1:query_start-1+150]  # Adjusting for 0-based indexing in Python
    return subread



#def positiontoBase(parseCigarArray, query_positions2):
 #   qbase = []
  #  seq = parseCigarArray[0]
   # for x in query_positions2:
    #    qbase.append(seq[x])

    #return qbase

def findRandomSubread(table):



#TESTS
def testFunctions(givenBamfile):
    parseCigarArray = parseCigar(givenBamfile) #calling parseCigar
    SoftorHardArray = (parseCigarArray[2])
    SoftorHard = SoftorHardArray[0]
    QRtable2Array = QRtable2(givenBamfile, SoftorHard)

    correctData = QRtable1(givenBamfile)
    for x in correctData:
        print(x)

    query_positions2 = QRtable2Array[0]
    reference_positions2 = QRtable2Array[1]
    splices = QRtable2Array[2]

    print(query_positions2)
    print(reference_positions2)
    print(splices)

    findRefSplicePos(query_positions2, reference_positions2, splices)
    arr = findRefSplicePos(query_positions2, reference_positions2, splices)
    dict = arr[0]
    a = []

    for x in splices:
        b = [x[0], dict[x[0]]]
        c = [x[1], dict[x[1]]]

        a.append([b,c])

    for x in a:
        for y in x:
            print(y)
        print("")


   # print(positiontoBase(parseCigarArray, query_positions2))


testFunctions(bamfile)


bamfile.close()