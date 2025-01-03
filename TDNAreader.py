#!/home/hongb/anaconda3/bin/python3
import subprocess as sp
import sys
import os
import argparse
import time
import multiprocessing as mp
from functools import partial
import pandas as pd
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description = 
    '''
    T-DNAreader is a bioinformatics tool designed to identify T-DNA insertion sites (TISs) in plant genomes using sequencing data.
    Contact address: miso5103@snu.ac.kr
    ''', formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument('-b', dest = 'bam', required=True, help = 
    '''Input BAM or SAM file(s).
Multiple files should be separated by ',' (e.g., rep1.bam,rep2.bam,rep3.bam).
For paired-end reads, pairs of BAM files should be separated by ':' (e.g., rep1_1.bam:rep1_2.bam,rep2_1.bam:rep2_2.bam).
Strongly recommend using all biological replicates.
    ''')
    parser.add_argument('--paired', dest = 'paired', action='store_true', help = 
    '''The input BAM files are from paired-end sequencing.
**Note: it is required to align paired-end reads separately to the plant genome (treating them as single-end reads).
Paired BAM files should be delimited by ':' in the input.
    ''')
    parser.add_argument('-tdna', dest = 'tdna', help = 
    '''Path to the T-DNA fasta file containing sequences from LB to RB.
A Bowtie2 index for this fasta file is required and must be created prior to running T-DNAreader.
Bowtie2 index can be generated using 'bowtie2-build tdna.fa'
Sequences outside the LB-RB regions (e.g., vector backbone sequences) should be excluded.
    ''')
    parser.add_argument('-o', dest = 'outdir', default = './', help = 
    '''Path to the output directory where all output files will be written.
Default = './'
    ''')
    parser.add_argument('-p', dest = 'prefix', required=True, help = 
    '''Prefix for naming output files.
    ''')
    parser.add_argument('-bL', dest = 'bowtieL', default = 10, type=int, help = 
    '''Length of seed substring (in bp) used in Bowtie2 alignment.
Must be >3 and <=l1. 
Default = 10
    ''')
    parser.add_argument('-d', dest = 'd', default=200, type=int, help = 
    '''Maximum distance from T-DNA borders to distinguish between border regions and intergnal regions.
Reads mapped within this distance from T-DNA borders used the threshold specified by -l1.
Reads mapped beyond this distance use the threshold specified by l2.
Default = 200 (100<=d<=300 is recommended).
    ''')
    parser.add_argument('-l1', dest = 'l1', default = 18, type=int, help = 
    '''Minimum length of sequences matched to T-DNA junctions (reads mapped close to T-DNA borders).
Higher values increase stringecy, potentially reducing false positives but may miss some true TISs.
Default = 18 (15<=l1<=300 is recommended). 
    ''')
    parser.add_argument('-l2', dest = 'l2', default = 30, type=int, help = 
    '''Minimum length of sequences matched with T-DNA internal regions (reads mapped far from T-DNA borders).
Higher values increase stringecy, potentially reducing false positives but may miss some true TISs.
Default = 30 (30<=l2<=50 is recommended)
    ''')
    parser.add_argument('-m', dest = 'm', default=1, type=int, help = 
    '''Maximum number of mismatches allowed in T-DNA alignment.
Default = 1 (bp). (0<=m<=3 is recommended)
    ''')
    parser.add_argument('-t', dest = 'thread', default = 1, type = int, help = 
    '''Number of threads to use. 
Default = 1. (t = {number of samples} is recommended more most efficient analysis)
    ''')
    parser.add_argument('-bed', dest='bed', default = False, help = 
    '''(Optional) A bed file containing genomic information. 
BED6 format is required (6 columns).
Gene names or identifiers should be in column 4.
If specified, overlapping genes with TISs will be identified.
    ''')
    parser.add_argument('-bl', dest = 'blacklist', default = False, help = 
    '''(Optional) A bed file specifying regions to exclude from analysis.
If the inserted T-DNA of transgenic plants contain endogenous genomic sequences, 
it is recommended to specify the positions of those sequences to avoid misinterpretation.
    ''')
    parser.add_argument('--tmp', dest = 'keeptmp', default = False, action = 'store_true', help = 
    '''(Optional) Keep intermediate files generated by T-DNAreader''')
    return parser.parse_args()

def MaketmpDir(outdir, prefix):
    if os.path.isdir(f'{outdir}') == False:
        os.mkdir(f'{outdir}')
    if os.path.isdir(f'{outdir}/tmp') == False:
        os.mkdir(f'{outdir}/tmp')

def getFirstGap(cigar):
    return int(cigar.split('S')[0])
def getLastGap(cigar):
    return int(cigar.split('M')[-1].split('S')[0])

def getSoftClippedSeqFas(outdir, l1, inprefix, outprefix):
    infile = f'{inprefix}.SoftClipped.rmdups.sam'
    outfile = f'{outdir}/tmp/{outprefix}.SoftClipped.fastq'
    with open(infile) as i1, open(outfile,'w') as o1:
        n=1
        for line in i1:
            if line[0] == "@":
                continue
            splits = line.rstrip().split('\t')
            readid, cigar, seq, seqQ = splits[0], splits[5], splits[9], splits[10]
            seq_first, seq_last, seqQ_first, seqQ_last = False, False, False, False

            if cigar[-1] == 'S':
                s2 = getLastGap(cigar)
                if s2 >= l1:
                    seq_last = seq[-s2:]
                    seqQ_last = seqQ[-s2:]
                
                if 'S' in cigar[:-1]:
                    s1 = getFirstGap(cigar)
                    if s1 >= l1:
                        seq_first = seq[:s1]
                        seqQ_first = seqQ[:s1]      
            else:
                s1 = getFirstGap(cigar)
                if s1 >= l1:
                    seq_first = seq[:s1]
                    seqQ_first = seqQ[:s1]

            if seq_first:
                o1.write(f'@{readid}_L {n} length={s1}' + '\n' + seq_first + '\n' + f'+{readid}_L {n} length={s1}' + '\n' + seqQ_first + '\n')
                n += 1
            if seq_last:
                o1.write(f'@{readid}_R {n} length={s2}' + '\n' + seq_last + '\n' + f'+{readid}_R {n} length={s2}' + '\n' + seqQ_last + '\n')
                n += 1

def getSoftClippedSeqFasParallel(reps, prefixlist, outdir, l1, t):

    pool = mp.Pool(processes=t)
    func_parallel = partial(getSoftClippedSeqFas, outdir, l1)
    pool.starmap(func_parallel, zip(reps, prefixlist))
    pool.close()
    pool.join()

## Split local and global 
def filterSAM(outdir, l1:int, d:int, l2:int, m:int, prefix):
    infile = f'{outdir}/tmp/{prefix}.SoftClipped.tdna.sam'
    outfile = f'{outdir}/tmp/{prefix}.SoftClipped.tdna.filtered.txt'

    with open(infile) as i1, open(outfile, 'w') as o1:
        tdnaSizeDic = dict()
        for line in i1:
            if line[:3] == "@SQ":
                tdna = line.rstrip().split('\t')[1].split(':')[1]
                size = int(line.rstrip().split(':')[-1])
                tdnaSizeDic[tdna] = size
            elif line[0] == "@":
                pass
            else:
                splits = line.rstrip().split('\t')
                readid, flag, tdna, pos, cigar, seq, md = splits[0], splits[1], splits[2], splits[3], splits[5], splits[9], splits[-2]
                border = readid.split('_')[1]
                readid = readid.split('_')[0]

                size = tdnaSizeDic[tdna]
                mismatch = 0
                ## Clip cigar ##
                if cigar[-1] == "S":
                    cigar = 'M'.join(cigar.split('M')[:-1]) + 'M'
                if 'S' in cigar:
                    cigar = cigar.split('S')[1]
                ################

                ## Calculate mismatches (including in/del) ##
                if 'I' in cigar:
                    N_in = int(cigar.split('I')[0].split('M')[1])
                    mismatch = mismatch + N_in
                mismatch = mismatch + sum([md.split(':')[-1].count(i) for i in ['A', 'T', 'G', 'C']])
                if mismatch > m:
                    continue
                #############################################

                ## Position cutoff (near T-DNA border) + matched basepair cutoff ##
                if int(pos) < d or int(pos) > (size - d):
                    matchedList = md.split(':')[-1].replace('A','C').replace('T','C').replace('G','C').replace('^','').split('C')
                    matchedLength = sum([int(i) for i in matchedList if i != ''])
                    if matchedLength >= l1:
                        if int(pos) < d:
                            o1.write('\t'.join([readid, border, flag, tdna, pos, cigar, seq, str(matchedLength)]) + '\n')
                        else:
                            o1.write('\t'.join([readid, border, flag, tdna, f'{pos}({str((int(pos)-size))})', cigar, seq, str(matchedLength)]) + '\n')
                ###################################################################
                ## Only matched basepair cutoff ##
                else:
                    matchedList = md.split(':')[-1].replace('A','C').replace('T','C').replace('G','C').replace('^','').split('C')
                    matchedLength = sum([int(i) for i in matchedList if i != ''])
                    if matchedLength >= l2:
                        o1.write('\t'.join([readid, border, flag, tdna, pos, cigar, seq, str(matchedLength)]) + '\n')
                ###################################

def filterSAMParallel(prefixlist, outdir, l1, d, l2, m, t):
    pool = mp.Pool(processes=t)
    func_parallel = partial(filterSAM, outdir, l1, d, l2, m)
    pool.map(func_parallel, prefixlist)
    pool.close()
    pool.join()

def MapEndPositionParallel(prefixlist, outdir, t):
    pool = mp.Pool(processes=t)
    func_parallel = partial(MapEndPosition, outdir)
    pool.map(func_parallel, prefixlist)
    pool.close()
    pool.join()

def MapEndPosition(outdir, prefix):
    infile = f'{outdir}/tmp/{prefix}.tdnaSites.txt'
    outfile = f'{outdir}/tmp/{prefix}.tdnaSites.addEnd.txt'
    with open(infile) as i1, open(outfile,'w') as o1:
        startCount = dict()
        endCount = dict()
        idlist = []
        for line in i1:
            splits = line.rstrip().split('\t')
            chrom, start, readid, cigar = splits[0], splits[1], splits[3], splits[5]
            end = str(getMappedPosition(int(start), cigar))
            dup_id = ':'.join([chrom, start, end, cigar])
            ## Remove remaining duplicate reads
            if dup_id not in idlist:
                idlist.append(dup_id)
                splits[2] = end
                o1.write('\t'.join(splits) + '\n')

def getMappedPosition(start, cigar):
    if cigar[-1] == "S":
        cigar = 'M'.join(cigar.split('M')[:-1]) + 'M'
    if 'S' in cigar:
        cigar = cigar.split('S')[1]
    cigar = cigar.replace('N','M').replace('D','M')
    Isplits = cigar.split('I')
    Msplits = []
    for i in Isplits:
        while 1:
            if i[-1] != 'M':
                i = i[:-1]
            else:
                Msplits.append(i)
                break    
    extension = sum([int(i) for i in ''.join(Msplits)[:-1].split('M')])
    return (start + extension - 1)

def reverseComplement(seq):
    bpDic = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N' : 'N'}
    r_seq = ''
    for bp in seq[::-1]:
        r_seq = r_seq+bpDic[bp]
    return r_seq

def FilterNoise(outdir, prefixlist, layout, outfile):
    dflist = []
    for prefix in prefixlist:
        genomic = pd.read_csv(f'{outdir}/tmp/{prefix}.tdnaSites.addEnd.txt', sep='\t', header=None, names=['chrom', 'start', 'end', 'readid', 'flag', 'cigar', 'seq'])
        tdna = pd.read_csv(f'{outdir}/tmp/{prefix}.SoftClipped.tdna.filtered.txt', sep='\t', header=None, names=['readid', 'border', 'tdna_flag', 'tdna', 'tdna_pos', 'tdna_cigar', 'tdna_seq', 'matchedL'])
        merged = genomic.merge(tdna, on='readid')[['chrom', 'start', 'end', 'readid', 'flag', 'cigar', 'seq', 'tdna', 'border', 'tdna_flag', 'tdna_pos', 'tdna_cigar', 'matchedL', 'tdna_seq']]
        merged['prefix'] = prefix
        merged['FL'] = merged.cigar.str.split('S').str[0]
        merged['RL'] = merged.cigar.str.split('M').str[-1].str.split('S').str[0]
        merged['FL'] = np.where(merged.FL.str.contains('M'), 0, merged['FL'])
        merged['FL'] = merged['FL'].astype(int)
        merged['RL'] = merged['RL'].replace(to_replace='', value=0)
        merged['RL'] = merged['RL'].astype(int)        
        merged = merged[merged[['FL', 'RL']].max(axis=1) >= merged.tdna_seq.str.len()]
        dflist.append(merged.iloc[:,:-2])
    mergeddf = pd.concat(dflist, ignore_index=True).sort_values(by=['chrom', 'start', 'end', 'cigar', 'readid'])

    if layout == True:
        mergeddf['sample'] = mergeddf['prefix'].str.split('-').str[:-1].apply(lambda x : "-".join(x))
        mergeddf = mergeddf.drop_duplicates(subset=['chrom', 'start', 'cigar', 'sample'], keep='first')
        mergeddf = mergeddf.drop_duplicates(subset=['chrom', 'end', 'cigar', 'sample'], keep='first')
        mergeddf = mergeddf.drop_duplicates(subset=['chrom', 'start', 'seq', 'sample'], keep='first')
        mergeddf = mergeddf.drop_duplicates(subset=['chrom', 'end', 'seq', 'sample'], keep='first')
        mergeddf = mergeddf.iloc[:,:-1]
    else:
        mergeddf = mergeddf.drop_duplicates(subset=['chrom', 'start', 'cigar', 'prefix'], keep='first')
        mergeddf = mergeddf.drop_duplicates(subset=['chrom', 'end', 'cigar', 'prefix'], keep='first')
        mergeddf = mergeddf.drop_duplicates(subset=['chrom', 'start', 'seq', 'prefix'], keep='first')
        mergeddf = mergeddf.drop_duplicates(subset=['chrom', 'end', 'seq', 'prefix'], keep='first')
    startPos = dict()
    endPos = dict()
    idlist = list()
    if len(mergeddf) > 0:
        for idx, row in mergeddf.iterrows():
            chrom, start, end, readid, lr = row[0], row[1], row[2], row[3], row[8]
            if readid not in idlist:
                if lr == 'L':
                    if f'{chrom}:{start}' in startPos:
                        startPos[f'{chrom}:{start}'] += 1
                    else:
                        startPos[f'{chrom}:{start}'] = 1
                elif lr == 'R':
                    if f'{chrom}:{end}' in endPos:
                        endPos[f'{chrom}:{end}'] += 1
                    else:
                        endPos[f'{chrom}:{end}'] = 1
                idlist.append(readid)
    else:
        print("## No T-DNA insertion identified. ")
        return True
    filtereddf = []
    LBs = [keys for keys in startPos if startPos[keys] > 1]
    RBs = [keys for keys in endPos if endPos[keys] > 1]                       

    for regions in LBs:
        chrom, start = regions.split(':')
        start = int(start)
        LB_df = mergeddf.query('(chrom == @chrom) and (start == @start)').reset_index(drop=True)
        LB_df['end'] = start
        filtereddf.append(LB_df)
    for regions in RBs:
        chrom, end = regions.split(':')
        end = int(end)
        RB_df = mergeddf.query('(chrom == @chrom) and (end == @end)').reset_index(drop=True)
        RB_df['start'] = end
        filtereddf.append(RB_df)
    
    try:
        outdf = pd.concat(filtereddf, ignore_index=True)
    except ValueError:
        return 0
    
    seqlist = []
    for idx, row in outdf.iterrows():
        seq, tseq, t_direction, lr = row[6], row[13], row[9], row[8]
        if t_direction == 0:
            expr = '>>>'
        elif t_direction == 16:
            tseq = reverseComplement(tseq)
            expr = '<<<'

        pos = seq.find(tseq)
        seq2 = seq[:pos] + expr + seq[pos:(pos+len(tseq))] + expr + seq[(pos+len(tseq)):]
        seqlist.append(seq2)
    outdf.loc[:,'seq'] = seqlist
    outdf = outdf.sort_values(by = ['chrom', 'start', 'end'])
    outdf.to_csv(f'{outdir}/tmp/{outfile}.TDNA.tmp', sep='\t', index=False)

def Cleartmp(outdir, reps, prefix, prefixlist):
    sp.run(f'rm -r {outdir}/tmp', shell=True)
    for rep in reps:
        sp.run(f'rm {rep}.SoftClipped.sam', shell=True)
        sp.run(f'rm {rep}.SoftClipped.rmdups.sam', shell=True)
    for p in prefixlist:
        sp.run(f'rm {p}.SoftClipped.rmdups.sam', shell=True)
if __name__ == '__main__':
    args = parse_args()
    MaketmpDir(args.outdir, args.prefix)
    if not os.path.isfile(args.tdna):
        print(f'### Check T-DNA.fa : {args.tdna}')

    ### Parse bamfiles ###
    reps = args.bam.split(',')
    n_sample = len(reps)
    if args.paired == True:
        allreps = []
        for rep in reps:
            for fr in rep.split(':'):
                allreps.append(fr)
        reps = allreps
    reps_join = ' '.join(reps)
    
    prefixlist = []
    for i in range(1, n_sample+1):
        if args.paired == True:
            prefixlist.append(f'{args.prefix}_{i}-F')
            prefixlist.append(f'{args.prefix}_{i}-R')
        else:
            prefixlist.append(f'{args.prefix}_{i}')
    prefix_join = ' '.join(prefixlist)
    
    print('### Sample name ###')
    for s, p in zip(reps, prefixlist):
        print(f'{s} : {p}')
    
    ### Extracting hybrid reads (T-DNA : Genome) ###
    
    # Filtering only soft clipped reads
    sp.run(f"parallel -j {args.thread} \"samtools head {{1}} > {{1}}.SoftClipped.sam ; samtools view {{1}} | awk '\$6 ~/S/' >> {{1}}.SoftClipped.sam ; \
           samtools markdup -r {{1}}.SoftClipped.sam {{1}}.SoftClipped.rmdups.sam\" ::: $(echo {reps_join})", shell=True)

    # Clipping unmatched sequences
    getSoftClippedSeqFasParallel(reps, prefixlist, args.outdir, args.l1, args.thread)
    sp.run(f"parallel -j {args.thread} \"trim_galore -q 10 --phred33 -o {args.outdir}/tmp/ --basename {{1}}.SoftClipped -j 1 --stringency 3 \
            --length {args.l1} --illumina {args.outdir}/tmp/{{1}}.SoftClipped.fastq --suppress_warn\" ::: $(echo {prefix_join})", shell=True)

    # Mapping local -> global with no overlaps
    sp.run(f"parallel -j {args.thread} \"bowtie2 -x {args.tdna} -U {args.outdir}/tmp/{{1}}.SoftClipped_trimmed.fq -L {args.bowtieL} --local |\
            awk '\$3 != \\\"*\\\"' | tee {args.outdir}/tmp/{{1}}.SoftClipped.tdna.sam | cut -f 1 | sort -k1,1 | uniq > \
            {args.outdir}/tmp/{{1}}.SoftClipped.tdna.readid \" ::: $(echo {prefix_join})", shell=True)
    sp.run(f"parallel -j {args.thread} \"bowtie2 -x {args.tdna} -U {args.outdir}/tmp/{{1}}.SoftClipped_trimmed.fq -L {args.bowtieL} | awk '\$3 != \\\"*\\\"' \
            | grep -wvf {args.outdir}/tmp/{{1}}.SoftClipped.tdna.readid >> {args.outdir}/tmp/{{1}}.SoftClipped.tdna.sam >> \
            {args.outdir}/tmp/{{1}}.SoftClipped.tdna.sam ; samtools sort {args.outdir}/tmp/{{1}}.SoftClipped.tdna.sam -o {args.outdir}/tmp/{{1}}.SoftClipped.tdna.sam \" \
            ::: $(echo {prefix_join})", shell=True)

    # Flitering T-DNA mapped reads with minimum length of matced sequences and maximum distance from T-DNA borders
    filterSAMParallel(prefixlist, args.outdir, args.l1, args.d, args.l2, args.m, args.thread)
    sp.run(f"parallel -j {args.thread} \"less {args.outdir}/tmp/{{1}}.SoftClipped.tdna.filtered.txt | cut -f 1 | sort -k1,1 | uniq > \
           {args.outdir}/tmp/{{1}}.SoftClipped.tdna.filtered.readid\" ::: $(echo {prefix_join})", shell=True)    
    
    # Finding T-DNA inserted junction sites
    for s, p in zip(reps, prefixlist):
        sp.run(f"ln -s {s}.SoftClipped.rmdups.sam {p}.SoftClipped.rmdups.sam", shell=True)
    sp.run(f"parallel -j {args.thread} \"less {{1}}.SoftClipped.rmdups.sam | grep -wf {args.outdir}/tmp/{{1}}.SoftClipped.tdna.filtered.readid | \
           sort -k3,3 -k4,4n -k2,2 | awk '{{print \$3, \$4, \$4, \$1, \$2, \$6, \$10}}' OFS='\t' > {args.outdir}/tmp/{{1}}.tdnaSites.txt\" ::: $(echo {prefix_join})", shell=True)
    
    ### Map start and end positions + filtering ###
    MapEndPositionParallel(prefixlist, args.outdir, args.thread)
    empty = FilterNoise(args.outdir, prefixlist, args.paired, args.prefix)
    if empty == True:
        if args.keeptmp == False:
            Cleartmp(args.outdir, reps, args.prefix, prefixlist)
        sys.exit()
    with open(f'{args.outdir}/{args.prefix}.TDNA.txt', 'w') as o1:
        o1.write('\t'.join(['Sample', 'Position(REF)', 'Geneid', 'Readid', 'Flag(REF)', 'Cigar(REF)', 'Sequence', 'TDNA', 'Flag(TDNA)', 'Position(TDNA)', 'Cigar(TDNA)', '\n']))
    
    if args.blacklist and args.bed:
        sp.run(f"bedtools intersect -a {args.outdir}/tmp/{args.prefix}.TDNA.tmp -b {args.blacklist} -v | \
               bedtools intersect -a - -b {args.bed} -wao | \
               awk '{{ if ($9 == \"L\") {{print $15, \"TDNA-\"$1\":\"$2, $19, $4, $5, $6, $7, $8, $10, $11, $12}} \
               else {{print $15, $1\":\"$2\"-TDNA\", $19, $4, $5, $6, $7, $8, $10, $11, $12}} }}' OFS='\t' >> {args.outdir}/{args.prefix}.TDNA.txt", shell=True)
    elif args.blacklist and args.bed == False:
        sp.run(f"tail -n+2 {args.outdir}/tmp/{args.prefix}.TDNA.tmp | \
               bedtools intersect -a - -b {args.blacklist} -v | \
               awk '{{ if ($9 == \"L\") {{print $15, \"TDNA-\"$1\":\"$2, \".\", $4, $5, $6, $7, $8, $10, $11, $12}} \
               else {{print $15, $1\":\"$2\"-TDNA\", \".\", $4, $5, $6, $7, $8, $10, $11, $12}} }}' OFS='\t' >> {args.outdir}/{args.prefix}.TDNA.txt", shell=True)
    elif args.blacklist == False and args.bed == False:
        sp.run(f"tail -n+2 {args.outdir}/tmp/{args.prefix}.TDNA.tmp | \
               awk '{{ if ($9 == \"L\") {{print $15, \"TDNA-\"$1\":\"$2, \".\", $4, $5, $6, $7, $8, $10, $11, $12}} \
               else {{print $15, $1\":\"$2\"-TDNA\", \".\", $4, $5, $6, $7, $8, $10, $11, $12}} }}' OFS='\t' >> {args.outdir}/{args.prefix}.TDNA.txt", shell=True)
    elif args.blacklist == False and args.bed:
        sp.run(f"bedtools intersect -a {args.outdir}/tmp/{args.prefix}.TDNA.tmp -b {args.bed} -wao | \
               awk '{{ if ($9 == \"L\") {{print $15, \"TDNA-\"$1\":\"$2, $19, $4, $5, $6, $7, $8, $10, $11, $12}} \
               else {{print $15, $1\":\"$2\"-TDNA\", $19, $4, $5, $6, $7, $8, $10, $11, $12}} }}' OFS='\t' >> {args.outdir}/{args.prefix}.TDNA.txt", shell=True)

    if args.keeptmp == False:
        Cleartmp(args.outdir, reps, args.prefix, prefixlist)
