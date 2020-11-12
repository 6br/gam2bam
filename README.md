# gam2bam

Convert GAM format to BAM format (regarding each node in GAM format as each reference sequence in BAM format)

## usage

Input: input.vg and input.gam

```bash
pip install pystream-protobuf
vg view input.vg > input.gfa
awk '/^S/{print ">"$2"\n"$3}' input.gfa | fold > input.fa
samtools faidx input.fa
touch void.sam
samtools view -bt input.fa.fai void.sam > header.bam
python3 main.py input.gam output.bam header.bam   
```

