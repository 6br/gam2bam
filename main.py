import stream #pip install pystream-protobuf
import vg_pb2
import pysam
import sys
import gzip
from numba import jit

header = { 'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1575, 'SN': 'chr1'},
                   {'LN': 1584, 'SN': 'chr2'}] }

gamfile = sys.argv[1]
tmpfilename = sys.argv[2] 
headerfile = sys.argv[3]

samheader = pysam.AlignmentFile(headerfile, "rb")
print(samheader.references)
print(samheader.lengths)
ref_len = { k: v for (k,v) in zip(samheader.references, samheader.lengths) }
print(ref_len)

@jit
def is_seq_valid(seq):
    """Returns True if sequence is DNA otherwise False"""
    valid_bases = ['A', 'T', 'G', 'C']
    for base in seq:
        if base not in valid_bases:
            return False
    return True

alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

@jit
def reverse_complement(seq):    
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq) 
    #print(bases)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    #print(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

@jit
def consume_cigar(cigars):
    return sum([i[1] for i in cigars if i[0] != 2])

@jit
def consume_ref(cigars):
    return sum([i[1] for i in cigars if i[0] != 1])

@jit
def to_cigar(edit):
   if edit.from_length == edit.to_length and not edit.sequence:
       return (0, edit.from_length) # (7, edit.from_length)
   elif edit.from_length == edit.to_length:
       return (0, edit.from_length) # (8. edit.to_length)
   elif edit.from_length:
       return (2, edit.from_length) # DEL
   elif edit.to_length: 
       return (1, edit.to_length) #INS

#f = gzip.open(gamfile, 'rb')
with stream.open(gamfile, 'rb') as istream:
    with pysam.AlignmentFile(tmpfilename, "wb", template=samheader) as outf:
        for data in istream:
            m = vg_pb2.Alignment()
            m.ParseFromString(data)
        # work with message

#for m in stream.parse(gamfile, vg_pb2.Alignment):
            read_pos = 0
            print(m)
            for s in m.path.mapping:
                cigar = tuple([ to_cigar(i) for i in s.edit])
                read_end = read_pos + consume_cigar(cigar)
                a = pysam.AlignedSegment(header=samheader.header)
                a.query_name = m.name #"read_28833_29006_6945"
                a.query_sequence = m.sequence[read_pos:read_end] #"AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
                if s.position.is_reverse:
                    pass
                    #cigar = tuple(sorted(cigar))
                    a.query_sequence = reverse_complement(m.sequence[read_pos:read_end])
                a.flag = 16 if s.position.is_reverse else 0
                #if s.position.is_reverse:
                    #cigar = tuple(sorted(cigar))
                a.reference_name = str(s.position.node_id)
                a.reference_start = s.position.offset if not s.position.is_reverse else ref_len[s.position.node_id] - consume_ref(cigar)  - s.position.offset
                a.mapping_quality = m.mapping_quality #20
                a.cigar = cigar #((0,10), (2,1), (0,25))
        #a.next_reference_id = 0
        #a.next_reference_start=199
        #a.template_length=167
        #a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
                a.tags = () #tuple([ to_cigar(i) for i in s.edit]) #(("NM", 1),
                  #("RG", "L1"))
                outf.write(a)
            read_end = read_pos

