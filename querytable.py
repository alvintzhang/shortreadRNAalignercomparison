import pysam
import sys
import random

bamfile = pysam.AlignmentFile("/Users/AlvinZhang2026/bio_data/several_reads2.bam", "rb")

# Enhancing accuracy and precision, this code allows for more thorough studying of gene expression and regulation.

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

def QRtable1(givenBamfile):
    query_positions = []
    reference_positions = []

    for read in givenBamfile.fetch():
        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
            query_positions.append(query_pos)
            reference_positions.append(ref_pos)

    array = [query_positions, reference_positions]
    return array

def QRtable2(givenBamfile, SoH):
    query_positions2 = []
    reference_positions2 = []

    splices = []  # holds last match and first match (beginning/end of splice)
    splices2 = []
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
                splices2.append([[queryindex - 1, referenceindex],[queryindex, referenceindex+le]])
                query_positions2.append(queryindex)
                reference_positions2.append(referenceindex)
                referenceindex += le

    array = [query_positions2, reference_positions2, splices, splices2]
    return array

def findRefSplicePos(qp, rp, s):
    dict = {}
    assert len(qp) == len(rp)
    for i in range(len(qp)):
        dict[qp[i]] = rp[i]

    qrDASplicePos = []
    for x in s:
        qrDASplicePos.append(
            "{}:{},{}:{}".format(x[0], dict[x[0]], x[1], dict[x[1]]))
    return [dict, qrDASplicePos]

def extract_subreads(query_seq, start_positions, length):
    subreads = []
    used_positions = set()

    for start_pos in start_positions:
        if start_pos not in used_positions and (start_pos + length <= len(query_seq)):
            subread = query_seq[start_pos:start_pos + length]
            subreads.append(subread)
            used_positions.update(range(start_pos, start_pos + length))

    return subreads

def extract_random_subreads(query_seq, query_positions, ref_positions, length):
    subreads = []
    num_elements = len(query_seq) // length

    for _ in range(num_elements):
        start_pos = random.randint(0, len(query_seq) - length)
        subread = query_seq[start_pos:start_pos + length]
        end_pos = start_pos + length - 1
        ref_start = ref_positions[start_pos]
        ref_end = ref_positions[end_pos]
        subreads.append((subread, start_pos, end_pos, ref_start, ref_end))

    return subreads

def compare_subreads(subreads, short_read_rna):
    matches = []
    for subread in subreads:
        if subread in short_read_rna:
            matches.append(subread)
    return matches

def subreads_to_fasta(subreads, output_file):
    with open(output_file, 'w') as fasta_file:
        for i, subread in enumerate(subreads):
            fasta_file.write(f">subread_{i + 1}\n")
            fasta_file.write(f"{subread}\n")

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

    print(QRtable2Array[3])

    print(query_positions2)
    print(reference_positions2)

    read_rna = seqs[0]

    length_of_subread = 150  # 150 base pairs
    subreads = extract_subreads(seqs[0], query_positions2, length_of_subread)

    matches = compare_subreads(subreads, read_rna)
    print("Subreads:", matches)

    output_file1 = 'new_subreads.fasta'

    output_file2 = "/Users/AlvinZhang2026/bio_data/new_subreads.fasta"
    subreads_to_fasta(matches, output_file2)

    subreads2 = extract_random_subreads(seqs[0], query_positions2, reference_positions2, length_of_subread)
    print("Random Subreads:")
    for i, (subread, start_pos, end_pos, ref_start, ref_end) in enumerate(subreads2):
        print(f"Subread {i + 1}: {subread} (Query Start: {start_pos}, Query End: {end_pos}, Ref Start: {ref_start}, Ref End: {ref_end})")

    output_file3 = "/Users/AlvinZhang2026/bio_data/random_subreads.fasta"

    subreadseqs = [subread[0] for subread in subreads2]
    subreads_to_fasta(subreadseqs, output_file3)

    matches2 = compare_subreads(subreadseqs, read_rna)

testFunctions(bamfile)

bamfile.close()
