import pysam
import sys

bamfile = pysam.AlignmentFile("/Users/AlvinZhang2026/bio_data/several_reads2.bam", "rb")


# parseCigar method which reads the given BAM file
# returns the sequences, CIGAR strings, and whether the left clip is soft or hard (SoH variable)
def parseCigar(givenBamfile):
    seq = []  # all the reads in the BAM file
    SoH = []  # soft or hard of left clip for each read
    cigar = []  # array of strings with CIGAR strings of all the potential reads

    for read in givenBamfile.fetch():  # if there is more than one read in the file
        if read.is_secondary:
            continue
        print("read name: ", read.query_name)
        print("reference name: ", read.reference_name)
        print("first position on ref: ", read.reference_start)
        print("first position on read: ", read.query_alignment_start)
        print("read sequence: ", read.query_sequence)
        seq.append(read.query_sequence)
        cigar.append(read.cigarstring)

        print("cigar: ", read.cigarstring)

        for (op, le) in read.cigartuples:
            if op == 0 or op == 7 or op == 8:
                print("match: ", le)
            elif op == 1:
                print("insertion: ", le)
            elif op == 2:
                print("deletion: ", le)
            elif op == 3:
                print("intron: ", le)
            elif op == 4:
                print("soft clip: ", le)
                if not SoH:
                    SoH.append(1)
            elif op == 5:
                if not SoH:
                    SoH.append(2)
                print("hard clip: ", le)

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

    splices = []  # holds last match and first match (beginning/end of splice)
    # splices holds positions of the query column/row of the table
    SoHindex = 0
    queryindex = 0
    referenceindex = 0

    for read in givenBamfile.fetch():

        if SoH == 1:
            queryindex = read.query_alignment_start
        elif SoH == 2:
            queryindex = 0
        else:
            print("No clipping Error")
            sys.exit(1)
        referenceindex = read.reference_start

        for op, le in read.cigartuples:
            if op == 0:  # match or mismatch
                for i in range(le):
                    query_positions2.append(queryindex)
                    reference_positions2.append(referenceindex)
                    queryindex += 1
                    referenceindex += 1
            elif op == 1:  # insertion
                queryindex += le
            elif op == 2:  # deletion
                referenceindex += le
            elif op == 3:  # splice
                splices.append([queryindex - 1, queryindex])
                referenceindex += le

    array = [query_positions2, reference_positions2, splices]
    return array


def findRefSplicePos(qp, rp, s):
    dict = {}
    assert len(qp) == len(rp)
    for i in range(len(qp)):
        dict[qp[i]] = rp[i]

    qrDASplicePos = []
    for x in s:
        qrDASplicePos.append(
            "{}:{},{}:{}".format(x[0], dict[x[0]], x[1], dict[x[1]]))  # printing out the splices and ref positions

    return [dict, qrDASplicePos]


def extract_subreads(query_seq, start_positions, length):
    # Extract the subreads from the query sequences
    subreads = []
    for start_pos in start_positions:
        subread = query_seq[start_pos:start_pos + length]  # Extract subread of given length
        subreads.append(subread)
    return subreads


def compare_subreads_with_rna(subreads, short_read_rna):
    # Compare extracted subreads with short-read RNA sequences
    matches = []
    for subread in subreads:
        if subread in short_read_rna:
            matches.append(subread)
    return matches


def testFunctions(givenBamfile):
    parseCigarArray = parseCigar(givenBamfile)  # calling parseCigar
    seqs = parseCigarArray[0]
    SoftorHardArray = parseCigarArray[2]
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

        a.append([b, c])

    for x in a:
        for y in x:
            print(y)
        print("")

    # Example short-read RNA sequence for comparison (replace with actual data)
    short_read_rna = seqs[0]

    # Extract subreads (example: first 10 positions)
    length_of_subread = 150
    subreads = extract_subreads(seqs[0], query_positions2[:10], length_of_subread)

    # Compare subreads with short-read RNA sequence
    matches = compare_subreads_with_rna(subreads, short_read_rna)
    print("Matches with short-read RNA:", matches)


testFunctions(bamfile)

bamfile.close()
