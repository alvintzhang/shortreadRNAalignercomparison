import pysam
import os
import sys

# Open the BAM file
bamfile = pysam.AlignmentFile("/Users/AlvinZhang2026/bio_data/several_reads.bam", "rb")

def parseCigar(givenBamfile):
    seq = []
    cigar = []
    SoH = []

    for read in givenBamfile.fetch():
        if read.is_secondary:
            continue
        seq.append((read.query_name, read.query_sequence))
        cigar.append(read.cigarstring)

        for (op, le) in read.cigartuples:
            if op == 4 and not SoH:
                SoH.append(1)
            elif op == 5 and not SoH:
                SoH.append(2)

    return [seq, cigar, SoH]

def QRtable1(givenBamfile):
    query_positions = []
    reference_positions = []

    for read in givenBamfile.fetch():
        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
            query_positions.append(query_pos)
            reference_positions.append(ref_pos)

    return [query_positions, reference_positions]

def QRtable2(givenBamfile, SoH):
    query_positions2 = []
    reference_positions2 = []
    splices = []
    splices2 = []

    for read in givenBamfile.fetch():
        if SoH == 1:
            queryindex = read.query_alignment_start
        elif SoH == 2:
            queryindex = 0
        else:
            sys.exit(1)
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
            elif op == 3:  # Skip
                splices.append([queryindex - 1, queryindex])
                splices2.append([[queryindex - 1, referenceindex], [queryindex, referenceindex + le]])
                query_positions2.append(queryindex)
                reference_positions2.append(referenceindex)
                referenceindex += le

    return [query_positions2, reference_positions2, splices, splices2]

def findRefSplicePos(qp, rp, s):
    ref_dict = {}
    assert len(qp) == len(rp)
    for i in range(len(qp)):
        ref_dict[qp[i]] = rp[i]

    qrDASplicePos = []
    for x in s:
        qrDASplicePos.append(
            "{}:{},{}:{}".format(x[0], ref_dict[x[0]], x[1], ref_dict[x[1]]))
    return [ref_dict, qrDASplicePos]

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

def compare_subreads(subreads, short_read_rna):
    matches = []
    for subread in subreads:
        if subread[0] in short_read_rna:
            matches.append(subread)
    return matches

def subreads_to_fasta(subreads, output_file):
    with open(output_file, 'w') as fasta_file:
        for subread, start, end, ref_start, ref_end, long_read_id, subread_id in subreads:
            fasta_file.write(f">{subread_id}\n")
            fasta_file.write(f"{subread}\n")

def subreads_to_bam(subreads, output_file):
    temp_bam_file = output_file + ".temp.bam"
    with pysam.AlignmentFile(temp_bam_file, "wb", header=bamfile.header) as bam_out:
        for i, (subread, start, end, ref_start, ref_end, long_read_id, subread_id) in enumerate(subreads):
            read = pysam.AlignedSegment(bamfile.header)  # Correct class name
            read.query_name = subread_id
            read.query_sequence = subread
            read.query_qualities = None
            read.reference_start = ref_start if ref_start is not None else start  # Using ref_start if available
            read.flag = 0
            read.cigar = [(0, len(subread))]
            bam_out.write(read)

    os.rename(temp_bam_file, output_file)

def testFunctions(givenBamfile):
    parseCigarArray = parseCigar(givenBamfile)
    seqs = parseCigarArray[0]
    SoHArray = parseCigarArray[2]
    SoH = SoHArray[0]
    QRtable2Array = QRtable2(givenBamfile, SoH)

    correctData = QRtable1(givenBamfile)
    for x in correctData:
        print(x)

    query_positions2 = QRtable2Array[0]
    reference_positions2 = QRtable2Array[1]
    splices = QRtable2Array[2]

    print(QRtable2Array[3])

    all_subreads = []  # Collecting subreads for all long reads
    for read_id, read_seq in seqs:
        print(f"Processing long read ID: {read_id}")
        length_of_subread = 150
        subreads = extract_subreads(read_seq, query_positions2, length_of_subread, read_id, reference_positions2)

        print(f"Subreads from {read_id}:")
        for subread, start, end, ref_start, ref_end, long_read_id, subread_id in subreads:
            print(f"Subread ID: {subread_id}, Subread: {subread}, Start: {start}, End: {end}")

        matches = compare_subreads(subreads, read_seq)
        print("Matched Subreads:", matches)

        all_subreads.extend(subreads)  # Collect subreads from all reads

    output_file2 = "/Users/AlvinZhang2026/bio_data/new_subreads.fasta"
    subreads_to_fasta(all_subreads, output_file2)

    output_file_bam = "/Users/AlvinZhang2026/bio_data/new_subreads.bam"
    subreads_to_bam(all_subreads, output_file_bam)

    # Generating random subreads for demonstration purposes
    print("Random Subreads:")
    for i, (subread, start_pos, end_pos, ref_start, ref_end, long_read_id, subread_id) in enumerate(all_subreads):
        subread_id = f"{long_read_id}_random_subread_{i + 1}_start_{start_pos}_end_{end_pos}"
        print(f"Subread ID: {subread_id}, Subread {i + 1}: {subread} (Query Start: {start_pos}, Query End: {end_pos}, Ref Start: {ref_start}, Ref End: {ref_end})")

    output_file3 = "/Users/AlvinZhang2026/bio_data/random_subreads.fasta"
    subreads_to_fasta([(subread, start, end, ref_start, ref_end, long_read_id, f"{long_read_id}_random_subread_{i + 1}_start_{start}_end_{end}") for i, (subread, start, end, ref_start, ref_end, long_read_id, subread_id) in enumerate(all_subreads)], output_file3)

# Run the test function
testFunctions(bamfile)

# Close the BAM file
bamfile.close()
