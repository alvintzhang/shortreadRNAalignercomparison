import pysam
import os
import random

# Open the BAM file
bamfile = pysam.AlignmentFile("/Users/AlvinZhang2026/hg002_revio_grch38_minimap2_juncbed.chr20.part.bam", "rb")

# Parses the CIGAR strings from the BAM file, returning sequences, CIGAR strings, and clip types (soft or hard).
def parseCigar(givenBamfile):
    seq = []  # List to store sequences and their query names
    cigar = []  # List to store CIGAR strings
    SoH = []  # Soft or hard clip indicator

    for read in givenBamfile.fetch():
        if read.is_secondary:
            continue
        seq.append((read.query_name, read.query_sequence))
        cigar.append(read.cigarstring)

        if read.cigartuples:
            for (op, le) in read.cigartuples:
                if op == 4 and not SoH:
                    SoH.append(1)  # Soft clip
                elif op == 5 and not SoH:
                    SoH.append(2)  # Hard clip

    return [seq, cigar, SoH]

# Creates arrays of query and reference positions, unable to differentiate between deletions and splices.
def QRtable1(givenBamfile):
    query_positions = []
    reference_positions = []

    for read in givenBamfile.fetch():
        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
            query_positions.append(query_pos)
            reference_positions.append(ref_pos)

    return [query_positions, reference_positions]

# Creates arrays of query and reference positions, including splice information, allowing differentiation between all operations.
def QRtable2(givenBamfile, SoH):
    query_positions2 = []
    reference_positions2 = []
    splices = []  # Positions of splices (query index)
    splices2 = []  # Detailed splice information (query and reference positions)

    for read in givenBamfile.fetch():
        queryindex = read.query_alignment_start if SoH == 1 else 0
        referenceindex = read.reference_start

        for op, le in read.cigartuples:
            if op == 0:  # Match
                for i in range(le):
                    query_positions2.append(queryindex)
                    reference_positions2.append(referenceindex)
                    queryindex += 1
                    referenceindex += 1
            elif op == 1:  # Insertion
                queryindex += le
            elif op == 2:  # Deletion
                referenceindex += le
            elif op == 3:  # Skip (splice)
                splices.append([queryindex - 1, queryindex])
                splices2.append([[queryindex - 1, referenceindex], [queryindex, referenceindex + le]])
                query_positions2.append(queryindex)
                reference_positions2.append(referenceindex)
                referenceindex += le

    return [query_positions2, reference_positions2, splices, splices2]

# Converts query positions of splices to reference positions, returning a dictionary and formatted results.
def findRefSplicePos(qp, rp, s):
    ref_dict = {}
    assert len(qp) == len(rp)
    for i in range(len(qp)):
        ref_dict[qp[i]] = rp[i]

    qrDASplicePos = []
    for x in s:
        qrDASplicePos.append(f"{x[0]}:{ref_dict[x[0]]},{x[1]}:{ref_dict[x[1]]}")
    return [ref_dict, qrDASplicePos]

# Extracts subreads from the query sequences, returning a list with detailed subread information.
def extract_subreads(query_seq, start_positions, length, long_read_id, reference_positions):
    subreads = []
    used_positions = set()

    for i, start_pos in enumerate(start_positions):
        if start_pos is not None and start_pos not in used_positions and (start_pos + length <= len(query_seq)):
            subread = query_seq[start_pos:start_pos + length]
            ref_start = reference_positions[i] if i < len(reference_positions) else None
            ref_end = reference_positions[i + length - 1] if i + length - 1 < len(reference_positions) else None
            subread_id = f"{long_read_id}_subread_{i + 1}_start_{start_pos}_end_{start_pos + length - 1}"
            subreads.append((subread, start_pos, start_pos + length - 1, ref_start, ref_end, long_read_id, subread_id))
            used_positions.update(range(start_pos, start_pos + length))

    return subreads

# Compares extracted subreads with short-read RNA sequences, returning matching subreads.
def compare_subreads(subreads, short_read_rna):
    matches = []
    for subread in subreads:
        if subread[0] in short_read_rna:
            matches.append(subread)
    return matches

# Converts subreads to a FASTA file format for further RNA sequencing analysis.
def subreads_to_fasta(subreads, output_file):
    with open(output_file, 'w') as fasta_file:
        for subread, start, end, ref_start, ref_end, long_read_id, subread_id in subreads:
            fasta_file.write(f">{subread_id}\n")
            fasta_file.write(f"{subread}\n")

# Converts subreads to a BAM file format for further analysis.
def subreads_to_bam(subreads, output_file):
    temp_bam_file = output_file + ".temp.bam"
    with pysam.AlignmentFile(temp_bam_file, "wb", header=bamfile.header) as bam_out:
        for i, (subread, start, end, ref_start, ref_end, long_read_id, subread_id) in enumerate(subreads):
            read = pysam.AlignedSegment(bamfile.header)
            read.query_name = subread_id
            read.query_sequence = subread
            read.query_qualities = None
            read.reference_start = ref_start if ref_start is not None else start
            read.flag = 0
            read.cigar = [(0, len(subread))]
            bam_out.write(read)

    os.rename(temp_bam_file, output_file)

# Generates random subreads for a given query sequence
def generate_random_subreads(query_seq, subread_length):
    num_subreads = len(query_seq) // subread_length  # Number of subreads based on the sequence length
    random_subreads = []
    for _ in range(num_subreads):
        if len(query_seq) <= subread_length:
            break
        start = random.randint(0, len(query_seq) - subread_length)
        subread = query_seq[start:start + subread_length]
        random_subreads.append((subread, start, start + subread_length - 1))
    return random_subreads

# Tests all the functions, printing their results.
def testFunctions(givenBamfile):
    parseCigarArray = parseCigar(givenBamfile)
    seqs = parseCigarArray[0]
    cigars = parseCigarArray[1]
    SoHArray = parseCigarArray[2]
    SoH = SoHArray[0] if SoHArray else 0
    QRtable2Array = QRtable2(givenBamfile, SoH)

    # Collect all subreads and random subreads for all long reads
    all_subreads = []
    random_subreads_list = []
    long_read_ids = []  # List to store long read sequence IDs
    long_read_cigars = []  # List to store long read CIGAR strings

    for (read_id, read_seq), cigar in zip(seqs, cigars):
        long_read_ids.append(read_id)  # Add the long read sequence ID to the list
        long_read_cigars.append(cigar)  # Add the long read CIGAR string to the list
        print(f"Processing long read ID: {read_id}")
        print(f"Long read CIGAR string: {cigar}")  # Print the CIGAR string

        length_of_subread = 150

        # Generate random subreads for the current long read
        random_subreads = generate_random_subreads(read_seq, length_of_subread)

        for i, (subread, start, end) in enumerate(random_subreads):
            subread_id = f"{read_id}_random_subread_{i + 1}_start_{start}_end_{end}"
            random_subreads_list.append((subread, start, end, subread_id))
            print(f"Random Subread ID: {subread_id}, Subread: {subread}, Start: {start}, End: {end}")

        all_subreads.extend(random_subreads)  # Collect random subreads from all reads

    # Print all long read sequence IDs
    print("Long Read Sequence IDs:")
    for long_read_id in long_read_ids:
        print(long_read_id)

    # Output random subreads to FASTA
    output_file3 = "/Users/AlvinZhang2026/bio_data/extractedrandomsubreads1.fasta"
    subreads_to_fasta([(subread, start, end, None, None, read_id, subread_id) for subread, start, end, subread_id in random_subreads_list], output_file3)

# Run the test function
testFunctions(bamfile)
