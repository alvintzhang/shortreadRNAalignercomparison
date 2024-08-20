import pysam

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
        seq.append((read.query_name, read.query_sequence))  # Store query name and sequence
        cigar.append(read.cigarstring)  # Store CIGAR string
        if read.cigartuples:
            for (op, le) in read.cigartuples:
                if op == 4 and not SoH:  # Check for soft clip
                    SoH.append(1)  # Soft clip indicator
                elif op == 5 and not SoH:  # Check for hard clip
                    SoH.append(2)  # Hard clip indicator

    return [seq, cigar, SoH]  # Return sequences, CIGAR strings, and clip information

# Creates arrays of query and reference positions, including splice information, allowing differentiation between all operations.
def QRtable2(givenBamfile, SoH):
    query_positions2 = []  # List to store query positions
    reference_positions2 = []  # List to store reference positions
    splices = []  # Positions of splices (query index)
    splices2 = []  # Detailed splice information (query and reference positions)
    cigars = []  # List to store CIGAR strings

    for read in givenBamfile.fetch():
        queryindex = read.query_alignment_start if SoH == 1 else 0  # Adjust query index based on clip type
        referenceindex = read.reference_start  # Initialize reference index
        cigar_string = []  # Temporary storage for CIGAR operations

        for op, le in read.cigartuples:
            if op == 0:  # Match
                cigar_string.append(f"{le}M")
                for i in range(le):
                    query_positions2.append(queryindex)  # Store query position
                    reference_positions2.append(referenceindex)  # Store reference position
                    queryindex += 1  # Increment query index
                    referenceindex += 1  # Increment reference index
            elif op == 1:  # Insertion
                cigar_string.append(f"{le}I")
                queryindex += le  # Skip inserted bases
            elif op == 2:  # Deletion
                cigar_string.append(f"{le}D")
                referenceindex += le  # Skip deleted bases
            elif op == 3:  # Skip (splice)
                cigar_string.append(f"{le}N")
                splices.append([queryindex - 1, queryindex])  # Store splice positions
                splices2.append([[queryindex - 1, referenceindex], [queryindex, referenceindex + le]])
                query_positions2.append(queryindex)
                reference_positions2.append(referenceindex)
                referenceindex += le  # Increment reference index by splice length

        cigars.append(''.join(cigar_string))  # Store the CIGAR string for the read

    return [query_positions2, reference_positions2, splices, splices2, cigars]  # Return positions, splice info, and CIGAR strings

# Generates the CIGAR string for a subread based on its query and reference positions
def generate_subread_cigar(query_positions, reference_positions, subread_start, subread_end):
    subread_cigar = []
    last_op = None
    count = 0

    for i in range(subread_start, subread_end + 1):
        if i >= len(query_positions):
            break
        if last_op is None:
            last_op = "M"  # Assume the first operation is a match
            count = 1
        elif query_positions[i] == query_positions[i - 1] + 1 and reference_positions[i] == reference_positions[i - 1] + 1:
            if last_op == "M":
                count += 1
            else:
                subread_cigar.append(f"{count}{last_op}")
                last_op = "M"
                count = 1
        elif query_positions[i] == query_positions[i - 1] + 1:
            if last_op == "I":
                count += 1
            else:
                subread_cigar.append(f"{count}{last_op}")
                last_op = "I"
                count = 1
        elif reference_positions[i] == reference_positions[i - 1] + 1:
            if last_op == "D":
                count += 1
            else:
                subread_cigar.append(f"{count}{last_op}")
                last_op = "D"
                count = 1
        else:
            if last_op is not None:
                subread_cigar.append(f"{count}{last_op}")
            count = 1
            last_op = "N"  # Assume it's a splice

    if count > 0 and last_op is not None:
        subread_cigar.append(f"{count}{last_op}")

    return ''.join(subread_cigar)

# Extracts subreads from the query sequences, returning a list with detailed subread information.
def extract_subreads(query_seq, start_positions, length, long_read_id, query_positions, reference_positions):
    subreads = []  # List to store subreads
    used_positions = set()  # Track positions that have been used

    for i, start_pos in enumerate(start_positions):
        if start_pos is not None and start_pos not in used_positions and (start_pos + length <= len(query_seq)):
            subread = query_seq[start_pos:start_pos + length]  # Extract subread
            query_start = query_positions[start_pos] if start_pos < len(query_positions) else None  # Query start position
            query_end = query_positions[start_pos + length - 1] if start_pos + length - 1 < len(query_positions) else None  # Query end position
            subread_cigar = generate_subread_cigar(query_positions, reference_positions, start_pos, start_pos + length - 1)
            subread_id = f"{long_read_id}_random_subread_{i + 1}_start_{query_start}_end_{query_end}_cigar_{subread_cigar}"  # Generate subread ID with CIGAR
            subreads.append((subread, query_start, query_end, long_read_id, subread_id))  # Store subread information
            used_positions.update(range(start_pos, start_pos + length))  # Mark positions as used

    return subreads  # Return list of subreads

# Converts subreads to a FASTA file format for further RNA sequencing analysis.
def subreads_to_fasta(subreads, output_file):
    with open(output_file, 'w') as fasta_file:
        for subread, query_start, query_end, long_read_id, subread_id in subreads:
            fasta_file.write(f">{subread_id}\n")  # Write header with subread ID
            fasta_file.write(f"{subread}\n")  # Write subread sequence

# Tests all the functions, printing their results.
def testFunctions(givenBamfile):
    parseCigarArray = parseCigar(givenBamfile)  # Parse CIGAR strings
    seqs = parseCigarArray[0]  # Extract sequences
    SoHArray = parseCigarArray[2]  # Extract soft/hard clip information
    SoH = SoHArray[0] if SoHArray else 0  # Determine soft/hard clip type
    QRtable2Array = QRtable2(givenBamfile, SoH)  # Generate query and reference position tables
    query_positions = QRtable2Array[0]  # Extract query positions
    reference_positions = QRtable2Array[1]  # Extract reference positions

    # Collect all subreads for all long reads
    all_subreads = []  # List to store all subreads
    long_read_ids = []  # List to store long read sequence IDs

    for read_id, read_seq in seqs:
        long_read_ids.append(read_id)  # Add the long read sequence ID to the list
        print(f"Processing long read ID: {read_id}")
        length_of_subread = 150  # Define length of subreads

        # Extract subreads for the current long read
        subreads = extract_subreads(read_seq, query_positions, length_of_subread, read_id, query_positions, reference_positions)

        for subread, query_start, query_end, long_read_id, subread_id in subreads:
            print(f"Subread ID: {subread_id}, Sequence: {subread}")
            all_subreads.append((subread, query_start, query_end, long_read_id, subread_id))  # Add subreads to the list

    # Write subreads to a FASTA file
    subreads_to_fasta(all_subreads, "/Users/AlvinZhang2026/bio_data/subreads_output.fasta")

# Run the test function
testFunctions(bamfile)
