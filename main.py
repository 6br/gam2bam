#/usr/bin/env python3
#pip install pystream-protobuf
import stream 
import vg_pb2
import pysam
import sys
import gzip
import logging
#from numba import jit

gamfile = sys.argv[1]
tmpfilename = sys.argv[2] 
headerfile = sys.argv[3]
mate_pair = sys.argv[4] == "True"
samheader = pysam.AlignmentFile(headerfile, "rb")

logging.info(samheader.references)
logging.info(samheader.lengths)
ref_len = { k: v for (k,v) in zip(samheader.references, samheader.lengths) }

#@jit
def is_seq_valid(seq):
    """Returns True if sequence is DNA otherwise False"""
    valid_bases = ['A', 'T', 'G', 'C']
    for base in seq:
        if base not in valid_bases:
            return False
    return True

alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):    
    return "".join(complement.get(base, base) for base in reversed(list(seq)))
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq) 
    logger.debug(bases)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    logger.debug(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

#@jit
def consume_cigar(cigars):
    return sum([i[1] for i in cigars if i[0] != 2])

#@jit
def consume_ref(cigars):
    return sum([i[1] for i in cigars if i[0] != 1])

#@jit
def to_cigar(edit, extended):
   if edit.from_length == edit.to_length:
       if extended and not edit.sequence:
           return (7, edit.from_length)
       elif extended:
           return (8, edit.from_length)
       else:
           return (0, edit.from_length)
   elif edit.from_length:
       return (2, edit.from_length) # DEL
   elif edit.to_length: 
       return (1, edit.to_length) #INS

with stream.open(gamfile, 'rb') as istream:
    with pysam.AlignmentFile(tmpfilename, "wb", template=samheader, threads=4) as outf:
        for data in istream:
            m = vg_pb2.Alignment()
            m.ParseFromString(data)
            read_pos = 0
            logger.debug(m)
            for s in m.path.mapping:
                cigar = tuple([ to_cigar(i, False) for i in s.edit])
                read_end = read_pos + consume_cigar(cigar)
                a = pysam.AlignedSegment(header=samheader.header)
                a.query_sequence = m.sequence[read_pos:read_end] # e.g. "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
                if s.position.is_reverse:
                    cigar = tuple(reversed(cigar))
                    a.query_sequence = reverse_complement(m.sequence[read_pos:read_end])
                a.flag = 16 if s.position.is_reverse else 0
                #if s.position.is_reverse:
                    #cigar = tuple(sorted(cigar))
                if mate_pair and m.fragment_next.name:
                    a.query_name = m.name.rstrip() + "_1"
                elif mate_pair and m.fragment_prev.name:
                    a.query_name = m.name.rstrip() + "_2"
                else:
                    a.query_name = m.name #e.g. "read_28833_29006_6945"
                a.reference_name = str(s.position.node_id)

                a.reference_start = s.position.offset if not s.position.is_reverse else ref_len[str(s.position.node_id)] - consume_ref(cigar)  - s.position.offset
                a.mapping_quality = m.mapping_quality #e.g. 20
                a.cigar = cigar #e.g. ((0,10), (2,1), (0,25))
                a.tags = () #e.g. (("NM", 1), ("RG", "L1"))
                outf.write(a)
                read_pos = read_end
