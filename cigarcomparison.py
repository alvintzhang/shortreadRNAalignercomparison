import pysam
import re
from tabulate import tabulate

'''This code processes short-read RNA-seq alignment data and compares the CIGAR string of the aligned subreads against the "original" CIGAR strings included in the sequence IDs. It then produces a short summary that gets printed out while also recording the raw data into a TSV file for viewing on excel or google sheets'''

def read_and_process_reads(long_read_bam, subread_bam):
    long_read_cigar_dict = {}  # Dictionary to store long read CIGAR strings and sequences
    subread_cigar_list = []  # List to store subread CIGAR strings and relevant info

    print("Opening BAM files...")
    long_read_file = pysam.AlignmentFile(long_read_bam, "rb")  # Open long read BAM file
    subread_file = pysam.AlignmentFile(subread_bam, "rb")  # Open subread BAM file

    print("Processing long reads...")
    # Process long reads and extract CIGAR strings and sequences
    for read in long_read_file.fetch():
        if read.cigarstring:
            base_seq_id = read.query_name.replace('/', '=')
            long_read_cigar_dict[base_seq_id] = (read.cigarstring, read.query_sequence)
        else:
            print(f"Warning: Read {read.query_name} has no CIGAR string")

    print(f"Long reads dictionary: {len(long_read_cigar_dict)} entries")

    print("Processing subreads...")
    # Process subreads and extract CIGAR strings, start and stop positions
    for read in subread_file.fetch():
        full_seq_id = read.query_name
        base_seq_id = full_seq_id.split('random_subread')[0].rstrip('_')
        base_seq_id = base_seq_id.replace('=', '/')
        cigar = read.cigarstring

        # Extract start and stop positions from the sequence ID
        start_stop_match = re.search(r'_start_(\d+)_end_(\d+)$', full_seq_id)
        if start_stop_match:
            start, stop = map(int, start_stop_match.groups())
        else:
            start, stop = 0, 0

        if cigar:
            subread_cigar_list.append((base_seq_id, cigar, start, stop, full_seq_id))
        else:
            print(f"Warning: Subread {full_seq_id} has no CIGAR string")

    print(f"Subread list: {len(subread_cigar_list)} entries")

    return long_read_cigar_dict, subread_cigar_list

def parse_cigar(cigar_string):
    """Parse CIGAR string into a list of (length, operation) tuples"""
    return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)]

def get_cigar_name(op):
    """Map CIGAR operation to descriptive name"""
    return {
        'M': 'match',
        'I': 'insertion',
        'D': 'deletion',
        'N': 'splice',
        'S': 'soft_clip',
        'H': 'hard_clip',
        'P': 'padding',
        '=': 'match',
        'X': 'mismatch'
    }.get(op, 'unknown')

def generate_expected_cigar_string(long_cigar, start, stop):
    """Generate expected CIGAR string for the given start and stop positions"""
    expected_cigar = []
    current_position = 0
    total_matches = 0

    for length, op in long_cigar:
        if current_position >= stop:
            break
        if current_position + length <= start:
            current_position += length
            continue

        overlap_start = max(start, current_position)
        overlap_end = min(stop, current_position + length)
        overlap_length = overlap_end - overlap_start

        if overlap_length > 0:
            expected_cigar.append((overlap_length, op))
            if op == 'M':
                total_matches += overlap_length

        current_position += length

    # Ensure the CIGAR string is exactly 150 base pairs long
    if total_matches < 150:
        remaining_matches = 150 - total_matches
        for i, (length, op) in enumerate(expected_cigar):
            if op == 'M':
                if length >= remaining_matches:
                    expected_cigar[i] = (length + remaining_matches, 'M')
                    remaining_matches = 0
                    break
                else:
                    remaining_matches -= length

        if remaining_matches > 0:
            expected_cigar.append((remaining_matches, 'M'))

    expected_cigar_string = ''.join([f"{length}{op}" for length, op in expected_cigar])
    return expected_cigar_string

def calculate_accuracy_precision(summary, shared_summary):
    """Calculate accuracy and precision based on summary and shared counts"""
    accurate_bases = sum(shared_summary[op] for op in ['match', 'insertion', 'deletion', 'splice', 'mismatch'])
    total_detected_bases = sum(summary[op] for op in ['match', 'insertion', 'deletion', 'splice', 'mismatch'])

    accuracy = accurate_bases / total_detected_bases if total_detected_bases > 0 else 0
    precision = accurate_bases / (accurate_bases + summary['unmatched_sub']) if (accurate_bases + summary['unmatched_sub']) > 0 else 0

    return accuracy, precision

def calculate_splice_junction_metrics(long_cigar, sub_cigar):
    """Calculate splice junction accuracy"""
    long_splices = [(length, op) for length, op in long_cigar if op == 'N']
    sub_splices = [(length, op) for length, op in sub_cigar if op == 'N']

    # Calculate intersection length (L1) and unique lengths (L2, L3)
    intersection_length = sum(min(l1, l2) for (l1, _), (l2, _) in zip(long_splices, sub_splices))
    long_only_length = sum(l1 for l1, op in long_splices if op == 'N' and (l1, op) not in sub_splices)
    sub_only_length = sum(l2 for l2, op in sub_splices if op == 'N' and (l2, op) not in long_splices)

    # L2 + L3 + L1 = union
    splice_accuracy = intersection_length / (intersection_length + long_only_length + sub_only_length) if (intersection_length + long_only_length + sub_only_length) > 0 else 0

    return splice_accuracy

def calculate_approximate_coordinate_accuracy(long_cigar, sub_cigar):
    """Calculate approximate coordinate accuracy based on CIGAR strings"""
    long_coords = []
    sub_coords = []
    current_position = 0

    # Track the position for long read CIGAR operations
    for length, op in long_cigar:
        if op in ['M', 'D', 'N']:
            long_coords.append((current_position, current_position + length))
        current_position += length

    current_position = 0
    # Track the position for subread CIGAR operations
    for length, op in sub_cigar:
        if op in ['M', 'D', 'N']:
            sub_coords.append((current_position, current_position + length))
        current_position += length

    intersection_length = 0

    # Calculate the overlap between long read and subread coordinates
    for long_start, long_end in long_coords:
        for sub_start, sub_end in sub_coords:
            overlap_start = max(long_start, sub_start)
            overlap_end = min(long_end, sub_end)
            if overlap_start < overlap_end:
                intersection_length += overlap_end - overlap_start

    # Calculate the union length
    union_length = (sum(end - start for start, end in long_coords) +
                    sum(end - start for start, end in sub_coords) - intersection_length)

    # Calculate the approximate coordinate accuracy
    coord_accuracy = intersection_length / union_length if union_length > 0 else 0

    return coord_accuracy

def align_subread_to_longread(long_read_cigar, long_read_sequence, subread_cigar, start, stop):
    """Align subread to long read and calculate accuracy, precision, and other metrics"""
    long_cigar = parse_cigar(long_read_cigar)
    sub_cigar = parse_cigar(subread_cigar)

    expected_cigar_string = generate_expected_cigar_string(long_cigar, start, stop)

    summary = {
        'match': 0,
        'insertion': 0,
        'deletion': 0,
        'splice': 0,
        'mismatch': 0,
        'soft_clip': 0,
        'unmatched_sub': 0
    }

    shared_summary = {
        'match': 0,
        'insertion': 0,
        'deletion': 0,
        'splice': 0,
        'mismatch': 0,
        'soft_clip': 0
    }

    expected_cigar_parsed = parse_cigar(expected_cigar_string)
    sub_cigar_parsed = parse_cigar(subread_cigar)

    expected_index = 0
    sub_index = 0

    while expected_index < len(expected_cigar_parsed) and sub_index < len(sub_cigar_parsed):
        expected_len, expected_op = expected_cigar_parsed[expected_index]
        sub_len, sub_op = sub_cigar_parsed[sub_index]

        if expected_op == sub_op:
            common_len = min(expected_len, sub_len)
            if get_cigar_name(expected_op) in shared_summary:
                shared_summary[get_cigar_name(expected_op)] += common_len

            expected_len -= common_len
            sub_len -= common_len

            if expected_len == 0:
                expected_index += 1
            else:
                expected_cigar_parsed[expected_index] = (expected_len, expected_op)

            if sub_len == 0:
                sub_index += 1
            else:
                sub_cigar_parsed[sub_index] = (sub_len, sub_op)
        else:
            expected_index += 1
            sub_index += 1

    # Calculate accuracy, precision, and other metrics
    accuracy, precision = calculate_accuracy_precision(summary, shared_summary)
    splice_accuracy = calculate_splice_junction_metrics(long_cigar, sub_cigar_parsed)
    coord_accuracy = calculate_approximate_coordinate_accuracy(long_cigar, sub_cigar_parsed)

    return accuracy, precision, splice_accuracy, coord_accuracy, expected_cigar_string

def main():
    long_read_bam = '/Users/AlvinZhang2026/hg002_revio_grch38_minimap2_juncbed.chr20.part.bam'
    subread_bam = '/Users/AlvinZhang2026/STARsortedoutput16.bam'

    long_read_cigar_dict, subread_cigar_list = read_and_process_reads(long_read_bam, subread_bam)

    results = []
    for subread_info in subread_cigar_list:
        base_seq_id, sub_cigar, start, stop, full_seq_id = subread_info
        long_read_cigar = long_read_cigar_dict.get(base_seq_id, (None, None))

        if long_read_cigar[0]:
            accuracy, precision, splice_accuracy, coord_accuracy, expected_cigar_string = align_subread_to_longread(
                long_read_cigar[0], long_read_cigar[1], sub_cigar, start, stop)

            results.append([
                full_seq_id,
                accuracy,
                precision,
                splice_accuracy,
                coord_accuracy,
                expected_cigar_string,
                sub_cigar
            ])

    # Print results in a table format
    print(tabulate(results, headers=['Subread ID', 'Accuracy', 'Precision', 'Splice Accuracy', 'Coordinate Accuracy', 'Expected CIGAR', 'Subread CIGAR'], tablefmt='grid'))

if __name__ == "__main__":
    main()
