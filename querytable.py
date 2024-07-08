import pysam
import sys

bamfile = pysam.AlignmentFile("/Users/AlvinZhang2026/bio_data/subread3.bam", "rb")


# This code processes and analyzes sequencing data in a BAM file to understand sequencing alignment, interpreting how
# reads map to the reference genome through long read rna sequencing (Minimap2). It then identifies splice sites and
# creates a table with the query positions and reference positions to extract 150 base subreads from. Using short read rna
# sequencing on the subreads, the code then compares the short read data with the original long read data to tests for
# accuracy and precision

# Enhancing accuracy and precision, this code allows for more thorough studying of gene expression and regulation.

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


# Because our splices array contains the index of where the splices are based on the query index
# We want to convert from the query index to the reference positions in order to find the actual bases where the splices are
# This function also prints out the results
def findRefSplicePos(qp, rp, s):
    dict = {}
    assert len(qp) == len(rp)
    for i in range(len(qp)):
        dict[qp[i]] = rp[i]

    # printing results
    qrDASplicePos = []
    for x in s:
        qrDASplicePos.append(
            "{}:{},{}:{}".format(x[0], dict[x[0]], x[1], dict[x[1]]))  # printing out the splices and ref positions
    # returns a dict with the both the query positions and the reference positions
    return [dict, qrDASplicePos]


# Extracts all subreads of the given length from the query sequences
# Returns a list with the subreads
def extract_subreads(query_seq, start_positions, length):
    # Extract the subreads from the query sequences
    subreads = []
    for start_pos in start_positions:
        if (start_pos + length <= len(query_seq)):
            subread = query_seq[start_pos:start_pos + length]  # Extract subread of given length
            subreads.append(subread)

    # returns the list of extracted subreads
    return subreads


# comparison method
def compare_subreads(subreads, short_read_rna):
    # Compare extracted subreads with short-read RNA sequences
    matches = []
    for subread in subreads:
        if subread in short_read_rna:
            matches.append(subread)
    return matches


# converts subreads to fasta file to run short read (STAR/HISAT2/Minimap2) RNA seq on
def subreads_to_fasta(subreads, output_file):
    with open(output_file, 'w') as fasta_file:
        for i, subread in enumerate(subreads):
            # Write the header line with a unique identifier
            fasta_file.write(f">subread_{i + 1}\n")
            # Write the sequence data
            fasta_file.write(f"{subread}\n")

# This function tests all the functions written above and prints out their results
def testFunctions(givenBamfile):
    parseCigarArray = parseCigar(givenBamfile)  # calling parseCigar
    seqs = parseCigarArray[0]
    SoftorHardArray = parseCigarArray[2]
    SoftorHard = SoftorHardArray[0]
    QRtable2Array = QRtable2(givenBamfile, SoftorHard)

    # first dataset
    correctData = QRtable1(givenBamfile)
    for x in correctData:
        print(x)

    query_positions2 = QRtable2Array[0]
    reference_positions2 = QRtable2Array[1]
    splices = QRtable2Array[2]

    # second dataset that also contains the positions of splices within the sequences
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

    # Minimap2 long read example sequence to run
    read_rna = seqs[0]

    # Extracts our appropriate subreads
    length_of_subread = 150  # 150 base pairs
    subreads = extract_subreads(seqs[0], query_positions2, length_of_subread)

    # Make sure that the subreads match back up with the original long read sequence results
    matches = compare_subreads(subreads, read_rna)
    print("Matches with short-read RNA:", matches)

    output_file1 = 'new_subreads.fasta'

    # fasta output converting from subreads in BAM format
    output_file2 = "/Users/AlvinZhang2026/bio_data/new_subreads.fasta"
    subreads_to_fasta(matches, output_file2)


testFunctions(bamfile)

bamfile.close()

import pysam
import sys

bamfile = pysam.AlignmentFile("/Users/AlvinZhang2026/bio_data/subread3.bam", "rb")


# This code processes and analyzes sequencing data in a BAM file to understand sequencing alignment, interpreting how
# reads map to the reference genome through long read rna sequencing (Minimap2). It then identifies splice sites and
# creates a table with the query positions and reference positions to extract 150 base subreads from. Using short read rna
# sequencing on the subreads, the code then compares the short read data with the original long read data to tests for
# accuracy and precision

# Enhancing accuracy and precision, this code allows for more thorough studying of gene expression and regulation.

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


# Because our splices array contains the index of where the splices are based on the query index
# We want to convert from the query index to the reference positions in order to find the actual bases where the splices are
# This function also prints out the results
def findRefSplicePos(qp, rp, s):
    dict = {}
    assert len(qp) == len(rp)
    for i in range(len(qp)):
        dict[qp[i]] = rp[i]

    # printing results
    qrDASplicePos = []
    for x in s:
        qrDASplicePos.append(
            "{}:{},{}:{}".format(x[0], dict[x[0]], x[1], dict[x[1]]))  # printing out the splices and ref positions
    # returns a dict with the both the query positions and the reference positions
    return [dict, qrDASplicePos]


# Extracts all subreads of the given length from the query sequences
# Returns a list with the subreads
def extract_subreads(query_seq, start_positions, length):
    # Extract the subreads from the query sequences
    subreads = []
    for start_pos in start_positions:
        if (start_pos + length <= len(query_seq)):
            subread = query_seq[start_pos:start_pos + length]  # Extract subread of given length
            subreads.append(subread)

    # returns the list of extracted subreads
    return subreads


# comparison method
def compare_subreads(subreads, longread):
    # Compare extracted subreads with the long-read sequence
    matches = []
    for subread in subreads:
        if subread in longread:
            matches.append(subread)
    return matches


# converts subreads to fasta file to run short read (STAR/HISAT2/Minimap2) RNA seq on
def subreads_to_fasta(subreads, output_file):
    with open(output_file, 'w') as fasta_file:
        for i, subread in enumerate(subreads):
            # Write the header line with a unique identifier
            fasta_file.write(f">subread_{i + 1}\n")
            # Write the sequence data
            fasta_file.write(f"{subread}\n")

# This function tests all the functions written above and prints out their results
def testFunctions(givenBamfile):
    parseCigarArray = parseCigar(givenBamfile)  # calling parseCigar
    seqs = parseCigarArray[0]
    SoftorHardArray = parseCigarArray[2]
    SoftorHard = SoftorHardArray[0]
    QRtable2Array = QRtable2(givenBamfile, SoftorHard)

    # first dataset
    correctData = QRtable1(givenBamfile)
    for x in correctData:
        print(x)

    query_positions2 = QRtable2Array[0]
    reference_positions2 = QRtable2Array[1]
    splices = QRtable2Array[2]

    # second dataset that also contains the positions of splices within the sequences
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

    # Minimap2 long read example sequence to run
    read_rna = seqs[0]

    # Extracts our appropriate subreads
    length_of_subread = 150  # 150 base pairs
    subreads = extract_subreads(seqs[0], query_positions2, length_of_subread)

    # Make sure that the subreads match back up with the original long read sequence results
    matches = compare_subreads(subreads, read_rna)
    print("Matches with short-read RNA:", matches)

    output_file1 = 'new_subreads.fasta'

    # fasta output converting from subreads in BAM format
    output_file2 = "/Users/AlvinZhang2026/bio_data/new_subreads.fasta"
    subreads_to_fasta(matches, output_file2)

    arr1 = []
    arr2 = []

    for x in range(10):
        arr1.append(query_positions2[x])
    for x in range(10):
        arr2.append(reference_positions2[x])

    print(arr1)
    print(arr2)

testFunctions(bamfile)

bamfile.close()

