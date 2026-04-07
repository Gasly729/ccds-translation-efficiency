#!/usr/bin/env python3
"""
Generate transcript-relative regions.bed for C. elegans from GFF3 annotation.

Reads a GFF3, builds exon structures for each mRNA transcript, and maps
CDS/UTR features from genomic to transcript-relative coordinates.

Output: BED6 format (transcript_id, start, end, region_type, 0, +)
where coordinates are 0-based transcript-relative.
"""

import sys
from collections import defaultdict

def parse_parent(attrs):
    """Extract Parent ID from GFF3 attributes."""
    for field in attrs.split(';'):
        if field.startswith('Parent='):
            return field[7:]  # Remove 'Parent='
    return None

def parse_id(attrs):
    """Extract ID from GFF3 attributes."""
    for field in attrs.split(';'):
        if field.startswith('ID='):
            return field[3:]
    return None

def genomic_to_transcript(genomic_start, genomic_end, exons, strand):
    """
    Map genomic interval [genomic_start, genomic_end) to transcript coordinates.
    Exons should be sorted by genomic position.
    Returns (tx_start, tx_end) in 0-based transcript coordinates.
    """
    tx_offset = 0
    tx_start = None
    tx_end = None
    
    for ex_start, ex_end in exons:
        ex_len = ex_end - ex_start
        # Overlap between feature and exon
        overlap_start = max(genomic_start, ex_start)
        overlap_end = min(genomic_end, ex_end)
        
        if overlap_start < overlap_end:
            if tx_start is None:
                tx_start = tx_offset + (overlap_start - ex_start)
            tx_end = tx_offset + (overlap_end - ex_start)
        
        tx_offset += ex_len
    
    if strand == '-':
        total_len = sum(e[1] - e[0] for e in exons)
        if tx_start is not None and tx_end is not None:
            new_start = total_len - tx_end
            new_end = total_len - tx_start
            tx_start, tx_end = new_start, new_end
    
    return tx_start, tx_end

def main(gff3_file, output_bed, lengths_file):
    # Step 1: Parse GFF3
    transcripts = {}  # transcript_id -> {chrom, strand, exons: [], features: []}
    
    print("Parsing GFF3...", file=sys.stderr)
    with open(gff3_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            chrom, source, ftype, start, end, score, strand, phase, attrs = parts
            start = int(start) - 1  # Convert to 0-based
            end = int(end)
            
            if ftype == 'mRNA':
                tid = parse_id(attrs)
                if tid:
                    # Remove prefix like 'transcript:'
                    transcripts[tid] = {
                        'chrom': chrom, 'strand': strand,
                        'exons': [], 'features': []
                    }
            
            elif ftype == 'exon':
                parent = parse_parent(attrs)
                if parent and parent in transcripts:
                    transcripts[parent]['exons'].append((start, end))
            
            elif ftype in ('CDS', 'five_prime_UTR', 'three_prime_UTR'):
                parent = parse_parent(attrs)
                if parent and parent in transcripts:
                    region_name = ftype
                    if ftype == 'five_prime_UTR':
                        region_name = 'UTR5'
                    elif ftype == 'three_prime_UTR':
                        region_name = 'UTR3'
                    transcripts[parent]['features'].append(
                        (start, end, region_name)
                    )
    
    print(f"Parsed {len(transcripts)} transcripts", file=sys.stderr)
    
    # Step 2: Read lengths file to get valid transcript IDs
    valid_ids = set()
    with open(lengths_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts:
                valid_ids.add(parts[0])
    
    print(f"Valid transcript IDs from lengths: {len(valid_ids)}", file=sys.stderr)
    
    # Step 3: Generate regions
    written = 0
    skipped_no_cds = 0
    with open(output_bed, 'w') as out:
        for tid, info in transcripts.items():
            # Check if this transcript ID (with prefix) is in our lengths file
            full_id = tid  # e.g., "transcript:Y74C9A.3.1"
            if full_id not in valid_ids:
                continue
            
            if not info['exons']:
                continue
            
            # Sort exons by position
            info['exons'].sort()
            strand = info['strand']
            
            # Check for CDS
            has_cds = any(f[2] == 'CDS' for f in info['features'])
            if not has_cds:
                skipped_no_cds += 1
                continue
            
            # Map features to transcript coordinates
            regions = []
            for gstart, gend, rtype in info['features']:
                tx_start, tx_end = genomic_to_transcript(
                    gstart, gend, info['exons'], strand
                )
                if tx_start is not None and tx_end is not None and tx_start < tx_end:
                    regions.append((tx_start, tx_end, rtype))
            
            # Sort by start position
            regions.sort()
            
            for tx_start, tx_end, rtype in regions:
                out.write(f"{full_id}\t{tx_start}\t{tx_end}\t{rtype}\t0\t+\n")
                written += 1
    
    print(f"Written {written} regions, skipped {skipped_no_cds} transcripts without CDS", 
          file=sys.stderr)

if __name__ == '__main__':
    gff3 = sys.argv[1]
    bed_out = sys.argv[2]
    lengths = sys.argv[3]
    main(gff3, bed_out, lengths)
