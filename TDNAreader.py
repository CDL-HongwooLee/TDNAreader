#!/home/hongb/anaconda3/bin/python3
"""
T-DNAreader: Identify T-DNA insertion sites (TISs) in plant genomes using RNA-seq data (or WGS).
Refactored for modularity, safety, and maintainability.
Contact: miso5103@snu.ac.kr
"""
import argparse
import logging
from pathlib import Path
from typing import Optional, List, Tuple, Dict, Any
import subprocess as sp
from multiprocessing import Pool
import pandas as pd
import numpy as np
import re
from TDNAdrawer import *#rangeBed, make_tdna_bam_bed, bam_to_bigwig_parallel, bigwig_to_bedgraph_parallel, create_track, plot_image
import sys
import os
from io import StringIO
from utils import run_cmd

_SOFT_START = re.compile(r"^(\d+)S")
_SOFT_END   = re.compile(r"(\d+)S$")
logger = logging.getLogger("TDNAreader")

def setup_logger(log_file: Path):
    logger = logging.getLogger("TDNAreader")
    logger.propagate = False
    logger.setLevel(logging.INFO)

    fh = logging.FileHandler(log_file, mode='w', encoding='utf-8')
    fh.setFormatter(logging.Formatter(
        "%(asctime)s %(levelname)s [%(name)s]: %(message)s", "%Y-%m-%d %H:%M:%S"
    ))
    logger.addHandler(fh)
    return logger

def confirm_overwrite(outdir: Path, prefix: str):
    # find anything in outdir whose name starts with the prefix
    existing = list(outdir.glob(f"{prefix}*"))
    if not existing:
        return  # nothing to overwrite

    print(f"Warning: {len(existing)} existing file(s)/dir(s) in '{outdir}' match '{prefix}*' and will be overwritten:")
    for p in existing:
        print("  ", p.name)
    ans = input("Continue and overwrite? [y/N] ").strip().lower()
    if ans not in ("y", "yes"):
        print("Aborted by user.")
        sys.exit(1)

def parse_int_tuple(s: str):
    return tuple(int(x) for x in s.split(','))

def parse_args():
    parser = argparse.ArgumentParser(description = 
    '''
    T-DNAreader is a bioinformatics tool designed to identify T-DNA insertion sites (TISs) in plant genomes using NGS data.
    Contact address: miso5103@snu.ac.kr
    ''', formatter_class = argparse.RawTextHelpFormatter)

    ### Input files 
    parser.add_argument('-fq', dest = 'fastq', required = False, default = None, help = 
    '''Input FASTQ files. Gzipped files are allowed with --gzip. 
For multiple replicates, separate files with commas (e.g. rep1.fastq,rep2.fastq).
For paired-end, separate pairs with colons (e.g. rep1_1.fastq:rep1_2.fastq,rep2_1.fastq:rep2_2.fastq).
    ''')
    parser.add_argument('-bam', dest = 'bam', required=False, default = None, help = 
    '''Input BAM or SAM files.
Separate multiple files with commas (e.g. rep1.bam,rep2.bam).
For paired-end, separate pairs with colons (e.g. rep1_1.bam:rep1_2.bam).
Cannot be used with --discordant
    ''')
    parser.add_argument('-ctl', dest = 'control', default = None, help = 
    '''(Optional) Control (WT) FASTQ or BAM files (same format as -fq/-bam).
Use commas to separate replicates and colons for paired-end files.
    ''')
    # Output and naming
    parser.add_argument('-o', dest = 'outdir', default = Path('./'), help = 
    '''Output directory where all output files will be saved.
Default = './'
    ''')
    parser.add_argument('-p', dest = 'prefix', required=True, help = 
    '''Prefix for naming all output files.
    ''')
    # Reference files
    parser.add_argument('-tdna', dest = 'tdna', required=True, help = 
    '''Path to T-DNA FASTA (LB to RB).
Requires a pre-built Bowtie2 index.
Exclude vector backbone sequences outside LB-RB.
    ''')
    parser.add_argument('-g', dest = 'genome', required = False, default = None, help = 
    '''Prefix for your genome index:
- STAR: path to genomeDir
- Bowtie2: path prefix to .bt2 files
(e.g., /path/to/genomeDir or /path/to/bowtie2/index_prefix).
    ''')
    parser.add_argument('-bed', dest='bed', default = None, help = 
    '''(Optional) BED6 file with genomic features. 
Column 4 should contain gene IDs.
Overlapping TISs will be annotated.
    ''')
    parser.add_argument('-gtf', dest='gtf', default = None, help = 
    '''(Optional) GTF file for STAR alignment and gene annotation.
    ''')
    parser.add_argument('-bl', dest = 'blacklist', default = None, help = 
    '''(Optional) BED file of regions to exclude.
If T-DNA contain endogenous genomic sequences, 
it is strongly recommended to specify the those positions to avoid misinterpretation.
    ''')
    
    ### Alignment options
    parser.add_argument('-r', dest = 'ratio', default = 0.33, help = 
    '''Minimum fraction of each read that must align to the genome (STAR: --outFilter*OverLread).
Default = 0.33 (recommended for 125-150 bp). 
Use ~0.5 for 80-125 bp, or
0.6 - 0.7 is recommend for 50-80 bp.
Only work for RNA-seq.
    ''')
    parser.add_argument('-i', dest = 'intron', default = 30000, type=int, help = 
    '''(Optional) Maximum intron size for STAR (--alignIntronMax).
Default = 30000
    ''')
    parser.add_argument('-bL', dest = 'bowtieL', default = 10, type=int, help = 
    '''(Optional) Seed length for Bowtie2 local alignment (bp).
Must be >3 and <= '-l1'. 
Default = 10
    ''')
    parser.add_argument('-mapq', dest = 'mapq', default = 10, type = int, help = 
    '''(Optional) Mapping quality threshold for WGS analysis.
Default = 30
    ''')

    ### T-DNA alignement and filtering
    parser.add_argument('-d', dest = 'd', default=200, type=int, help = 
    '''(Optional) Distance (bp) from T-DNA borders within which to apply the “-l1” match‐length threshold.
Reads mapped within this distance from T-DNA borders used the threshold specified by -l1.
Reads mapped beyond this distance use the threshold specified by l2.
Default = 200 (100<=d<=300 is recommended)
    ''')
    parser.add_argument('-l1', dest = 'l1', default = 18, type=int, help = 
    '''(Optional) Min matched sequence length (reads mapped close to T-DNA borders).
Higher increases stringency. 
Default = 18 (15<=l1<=30 is recommended)
    ''')
    parser.add_argument('-l2', dest = 'l2', default = 30, type=int, help = 
    '''(Optional) Min matched sequence length (reads mapped in internal T-DNA regions).
Higher increases stringency. 
Default = 30 (30<=l2<=50 is recommended)
    ''')
    parser.add_argument('-m', dest = 'm', default=1, type=int, help = 
    '''(Optional) Max mismatches allowed in T-DNA alignment.
Default = 1 (bp) (0<=m<=3 is recommended)
    ''')
    parser.add_argument('-w', dest = 'window', default=500, type=int, help = 
    '''(Optional) Merge distance for adjacent T-DNA insertion sites.
Default = 500 (bp)
    ''')
    parser.add_argument('-x', dest = 'weight', default=(700,300,300,1,1), type=lambda s: tuple(map(int, s.split(','))), help = 
    '''(Optional) Comma-separated integers defining the scoring weights.
Weights for group1-5 reads.
Default = 700,300,300,1,1
    ''')
    parser.add_argument('-u', dest = 'threshold', default=1000, type=int, help = 
    '''(Optional) Score threshold  to filter false positives.
Default = 1000
    ''')
    parser.add_argument('-t', dest = 'thread', default = 1, type = int, help = 
    '''Number of threads to use. 
Default = 1.
    ''')

    parser.add_argument('--paired', dest = 'paired', action='store_true', help = 
    '''Inputs are from paired-end data.
Separate pairs with ':'.
    ''')
    parser.add_argument('--gzip', dest = 'gzip', action='store_true', help = 
    '''Specify if FASTQ inputs are gzipped. 
    ''')
    parser.add_argument('--wgs', dest = 'wgs', default = False, action = 'store_true', help = 
    '''(Optional) The input is whole-genome sequencing data.
    ''')
    parser.add_argument('--discordant', dest = 'discordant', default = False, action = 'store_true', help = 
    '''(Optional) Include discordant pairs in scoring.
If specified, this increases runtime.
    ''')
    parser.add_argument('--tmp', dest = 'keeptmp', default = False, action = 'store_true', help = 
    '''(Optional) Keep intermediate temporary files''')

    return parser.parse_args()

def make_dir(outdir: Path, prefix: str):
    tmp_dir = (outdir / "tmp")
    aligned_dir = (outdir / f"{prefix}_aligned")
    tmp_dir.mkdir(parents=True, exist_ok=True)
    aligned_dir.mkdir(exist_ok=True)
    return tmp_dir, aligned_dir

def alignment_preprocessing(prefix: str, fastq: str, control: str, paired: bool):
    
    transgenic = fastq.split(',')
    controls = control.split(',') if control else []
    
    fastq_files = []
    for inp in transgenic:
        fastq_files.extend(inp.split(':') if paired else [inp])
    for inp in controls:
        fastq_files.extend(inp.split(':') if paired else [inp])

    transgenic_prefixes = [] ; control_prefixes = []
    for idx, inp in enumerate(transgenic, start=1):
        if paired and ':' in inp:
            transgenic_prefixes.extend([f"{prefix}_{idx}-F", f"{prefix}_{idx}-R"])
        else:
            transgenic_prefixes.append(f"{prefix}_{idx}")
    for idx, inp in enumerate(controls, start=1):
        if paired and ':' in inp:
            control_prefixes.extend([f"{prefix}_control_{idx}-F", f"{prefix}_control_{idx}-R"])
        else:
            control_prefixes.append(f"{prefix}_control_{idx}")
    print('\nSample information')
    for s, p in zip(fastq_files, transgenic_prefixes + control_prefixes):
        print(f'{s} : {p}')
    return fastq_files, transgenic_prefixes, control_prefixes

def star(fastq: str, prefix: str, aligned_dir: Path, genome: str, ratio: float, intron: int, t: int, gzip: bool, gtf: Optional[str] = None):
    with open(os.devnull, 'w') as devnull:
        cmd = (
            f"STAR --runMode alignReads --runThreadN {str(t)} --genomeDir {genome} --readFilesIn {fastq} "
            f"--outFileNamePrefix {aligned_dir}/{prefix} --outSAMtype BAM SortedByCoordinate --peOverlapNbasesMin 12 --peOverlapMMp 0.1 "
            f"--twopassMode Basic --outFilterMatchNminOverLread {str(ratio)} --outFilterScoreMinOverLread {str(ratio)} --outFilterMismatchNoverLmax {str(ratio)} "
            f"--outReadsUnmapped Fastx --alignIntronMax {str(intron)} --limitBAMsortRAM 20000000000 "
        )
        if gzip:
            cmd = cmd + f'--readFilesCommand zcat '
        if gtf:
            cmd = cmd + f'--sjdbGTFfile {gtf} '
        run_cmd(cmd, shell=True, silence=devnull)

def star_parallel(fastq_files: list[str], prefixes: list[str], aligned_dir: Path, genome: str, ratio: float, intron: int, t: int, gzip: bool, gtf: Optional[str] = None):
    procs = min(t, len(fastq_files))
    per_thread = max(1, t // len(fastq_files))
    with Pool(procs) as pool:
        pool.starmap(
            star,
            [(fq, pr, aligned_dir, genome, ratio, intron, per_thread, gzip, gtf)
             for fq, pr in zip(fastq_files, prefixes)]
        )

def bowtie2(fastq, prefix, aligned_dir, genome, mapq, t):
    with open(os.devnull, 'w') as devnull:
        cmd = (
            f"bowtie2 -x {genome} -U {fastq} --un {aligned_dir}/{prefix}Unmapped.out.mate1 -p {t} --local | samtools sort -O bam | "
            f"samtools view -q {mapq} -b | tee {aligned_dir}/{prefix}.sorted.bam | samtools index - {aligned_dir}/{prefix}.sorted.bam.bai"
        )
        run_cmd(cmd, shell=True, silence=devnull)
        
def bowtie2_parallel(fastq_files: list[str], prefixes: list[str], aligned_dir: Path, genome: str, mapq: int, t: int):
    procs = min(t, len(fastq_files))
    per_thread = max(1, t // len(fastq_files))

    with Pool(procs) as pool:
        pool.starmap(
            bowtie2,
            [
                (fq, pfx, aligned_dir, genome, mapq, per_thread)
                for fq, pfx in zip(fastq_files, prefixes)                
            ]
        )

def bam_preprocessing(aligned_dir: Path, prefix: str, bam: str, control: str, tag: str, paired: bool):
    
    transgenic = bam.split(',')
    controls = control.split(',') if control else []
    
    bam_files = []
    for inp in transgenic:
        bam_files.extend(inp.split(':') if paired else [inp])
    for inp in controls:
        bam_files.extend(inp.split(':') if paired else [inp])

    transgenic_prefixes = [] ; control_prefixes = []
    for idx, inp in enumerate(transgenic, start=1):
        if paired and ':' in inp:
            transgenic_prefixes.extend([f"{prefix}_{idx}-F", f"{prefix}_{idx}-R"])
        else:
            transgenic_prefixes.append(f"{prefix}_{idx}")
    for idx, inp in enumerate(controls, start=1):
        if paired and ':' in inp:
            control_prefixes.extend([f"{prefix}_control_{idx}-F", f"{prefix}_control_{idx}-R"])
        else:
            control_prefixes.append(f"{prefix}_control_{idx}")
    
    bam_links = []
    print('\nSample information')
    for b, p in zip(bam_files, transgenic_prefixes + control_prefixes):
        bam_path = Path(b).absolute()
        bam_link = aligned_dir / f'{p}{tag}'
        if bam_path == (bam_link).absolute():
            print(f'{b} : {p}')
        else:
            run_cmd(f"ln -s {bam_path} {bam_link}", shell=True)
            print(f'{b} -> {bam_link}: {p}')
            if Path(f'{bam_path}.bai').is_file():
                print(f'{bam_path}.bai')
                run_cmd(f"ln -s {bam_path}.bai {bam_link}.bai", shell=True)
            else:
                index = False
        bam_links.append(bam_link)
    try:
        index = index
    except NameError:
        index = True
    return bam_links, transgenic_prefixes, control_prefixes, index

def index_bams(prefixes: List[str], aligned_dir: Path, tag: str, threads: int):
    """Index BAM files in parallel."""
    stems = ' '.join(prefixes)
    parallel_cmd = (
        f"parallel -j {threads} "
        f"\"samtools index {aligned_dir}/{{1}}{tag}\" ::: {stems}"
    )
    run_cmd(parallel_cmd, shell=True)

def prepare_softclipped_sams(prefixes: list[str], aligned_dir: Path, tmp_dir: Path, tag: str, threads: int):
    """
    Generate SoftClipped SAMs by combining SAM header and soft-clipped lines, then mark duplicates.
    prefixes: list of sample prefixes (transgenic + control)
    """
    stems = ' '.join(prefixes)
    parallel_cmd = (
        f"parallel -j {threads} "
        f"\"samtools head {aligned_dir}/{{1}}{tag} > {tmp_dir}/{{1}}.SoftClipped.sam && "
        f"samtools view {aligned_dir}/{{1}}{tag} | awk '\$6 ~ /S/' >> {tmp_dir}/{{1}}.SoftClipped.sam && "
        f"samtools markdup -r {tmp_dir}/{{1}}.SoftClipped.sam {tmp_dir}/{{1}}.SoftClipped.rmdups.sam\" ::: $(echo {stems})"
    )
    run_cmd(parallel_cmd, shell=True)

def extract_soft_clipped_seq(sam_path: Path, fq_out: Path, min_clip: int):
    """Extract soft-clipped reads from SAM, preserve original-style header."""
    with open(sam_path) as sam, open(fq_out, 'w') as fq:
        for line in sam:
            if line.startswith('@'): continue
            cols = line.rstrip().split('\t')
            readid, flag, chrom, pos, _, cigar, _, _, _, seq_full, qual_full = (
                cols[0], cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8], cols[9], cols[10]
            )
            for match, side in (( _SOFT_START.match(cigar), 'L'), ( _SOFT_END.search(cigar), 'R')):
                if not match: continue
                length = int(match.group(1))
                if length < min_clip: continue
                seq = seq_full[:length] if side=='L' else seq_full[-length:]
                remain = seq_full[length:] if side=='L' else seq_full[:-length]
                if max(remain.count(b) for b in ("A","T","C","G")) > 0.8 * len(remain): continue
                qual= qual_full[:length] if side=='L' else qual_full[-length:]
                header = f"@{readid}::{chrom}::{pos}::{flag}::{cigar}::{seq_full}::{side}"
                fq.write(f"{header}\n{seq}\n+{header}\n{qual}\n")

def extract_soft_clipped_seq_parallel(prefixlist, tmpdir, min_clip, t):
    tasks = [(tmpdir / f'{prefix}.SoftClipped.rmdups.sam', tmpdir / f'{prefix}.SoftClipped.fastq', min_clip) for prefix in prefixlist]
    with Pool(t) as pool:    
        pool.starmap(extract_soft_clipped_seq, tasks)

def map_softclip_to_tdna_parallel(prefixes: list[str], tmp_dir: Path, tdna: str, bowtieL: int, t: int):
    cmd1 = (
        f'parallel -j {t} '
        f'"bowtie2 -x {tdna} -U {tmp_dir}/{{1}}.SoftClipped.fastq -L {bowtieL} --local | '
        f'awk \'\\$3!=\\"*\\"\' | tee {tmp_dir}/{{1}}.SoftClipped.tdna.sam | cut -f 1 | sort -u > '
        f'{tmp_dir}/{{1}}.SoftClipped.tdna.readid" ::: $(echo {" ".join(prefixes)})'
    )
    cmd2 = (
        f'parallel -j {t} '
        f'"bowtie2 -x {tdna} -U {tmp_dir}/{{1}}.SoftClipped.fastq -L {bowtieL} | '
        f'awk \'\\$3!=\\"*\\"\' | grep -wvf {tmp_dir}/{{1}}.SoftClipped.tdna.readid >> '
        f'{tmp_dir}/{{1}}.SoftClipped.tdna.sam ; samtools sort {tmp_dir}/{{1}}.SoftClipped.tdna.sam '
        f'-o {tmp_dir}/{{1}}.SoftClipped.tdna.sam" ::: $(echo {" ".join(prefixes)})'
    )
    with open(os.devnull, 'w') as devnull:
        run_cmd(cmd1, shell=True, silence=devnull)
        run_cmd(cmd2, shell=True, silence=devnull)

def get_matched_len(md: str):
    if md.startswith("MD:Z:"):
        md = md.split(':', 2)[2]
    nums = re.findall(r'\d+', md)
    return sum(int(n) for n in nums)

def get_matched_len_df(md1: Any, md2: Any):
    if isinstance(md1, str) and md1.startswith('MD:Z:'):
        return get_matched_len(md1)
    elif isinstance(md2, str) and md2.startswith('MD:Z:'):
        return get_matched_len(md2)
    else:
        print(f'Warning: wrong MD detected {md1}, {md2}')
    return 0

## Split local and global 
def filter_sam(prefix, tmpdir, l1:int, d:int, l2:int, m:int):
    infile = f'{tmpdir}/{prefix}.SoftClipped.tdna.sam'
    outfile = f'{tmpdir}/{prefix}.SoftClipped.tdna.filtered.txt'
    with open(infile) as i1, open(outfile, 'w') as o1:
        tdnaSizeDic = dict()
        for line in i1:
            if line.startswith("@SQ"):
                tdna = line.rstrip().split('\t')[1].split(':')[1]
                size = int(line.rstrip().split(':')[-1])
                tdnaSizeDic[tdna] = size
                continue
            elif line[0] == "@":
                continue

            splits = line.rstrip().split('\t')
            genomic_info, tdna_flag, tdna, tdna_pos, tdna_cigar, tdna_seq, tdna_nm, tdna_md = splits[0], splits[1], splits[2], splits[3], splits[5], splits[9], splits[-3], splits[-2]
            mismatch = int(tdna_nm.split(':')[-1])
            if mismatch > m:
                continue

            s = genomic_info.split('::')
            genomic_out = '\t'.join([s[1], s[2], s[2], s[0], s[3], s[4], s[6], s[5]])
            size = tdnaSizeDic[tdna]
            ## Position cutoff (near T-DNA border) + matched basepair cutoff ##
            near_border = int(tdna_pos) < d or int(tdna_pos) > size - d
            matched_len = get_matched_len(tdna_md)
            if near_border and matched_len >= l1:
                pos_out = f"{int(tdna_pos)}({int(tdna_pos) - size})" if int(tdna_pos) > size - d else str(tdna_pos)
            elif not near_border and matched_len >= l2:
                pos_out = str(tdna_pos)
            else:
                continue
            o1.write('\t'.join([genomic_out, tdna, tdna_flag, pos_out, tdna_cigar, str(matched_len), tdna_seq]) + '\n')

def filter_sam_parallel(prefixes: list[str], tmp_dir: str, l1: int, d: int, l2: int, m: int, t: int):
    with Pool(t) as pool:
        pool.starmap(
            filter_sam,
            [
                (pfx, tmp_dir, l1, d, l2, m)
                for pfx in prefixes
            ]
        )

def get_mapped_position(start: int, cigar: str):
    cigar_clean = re.sub(r'\d+S', '', cigar)
    cigar_clean = cigar_clean.replace('N', 'M').replace('D', 'M')
    total_M = sum(int(n) for n in re.findall(r'(\d+)M', cigar_clean))
    return start + total_M - 1

def concat_replicates(prefix: str, prefixes: List[str], tmp_dir: Path, paired: bool, control: bool):
    dflist = []
    for p in prefixes:
        df = pd.read_csv(f'{tmp_dir}/{p}.SoftClipped.tdna.filtered.txt', sep='\t', names=['chrom', 'start', 'end', 'readid', 'flag', 'cigar', 'clip_pos', 'seq', 'tdna', 'tdna_flag', 'tdna_pos', 'tdna_cigar', 'matchedL', 'tdna_seq'])
        df['prefix'] = p
        df['tdna_flag'] = np.where(df['tdna_flag'] == 0, '+', '-')
        dflist.append(df)
    mergeddf = pd.concat(dflist, ignore_index=True).sort_values(by=['chrom', 'start', 'end', 'cigar', 'readid'])
    
    if paired:
        mergeddf['sample'] = mergeddf['prefix'].str.rsplit('-', n=1).str[0]#.str[:-1].apply(lambda x : "-".join(x))
        mergeddf = (
            mergeddf
            .drop_duplicates(subset=['chrom', 'start', 'seq', 'sample'])
            .drop_duplicates(subset=['chrom', 'start', 'readid', 'clip_pos'])
            .drop_duplicates(subset=['chrom', 'end', 'readid', 'clip_pos'])
        )
    if not mergeddf.empty:
        mergeddf['end'] = mergeddf.apply(lambda row: get_mapped_position(int(row['start']), row['cigar']), axis=1)
    
    out_name = f"{prefix}_control.SoftClipped.tmp" if control else f"{prefix}.SoftClipped.tmp"
    mergeddf.iloc[:,:-1].to_csv(tmp_dir / out_name, sep='\t', index=False)
    return mergeddf.iloc[:,:-1]

def align_discordant_pairs_parallel(prefixes: List[str], aligned_dir: Path, tmp_dir: Path, tdna: str, bowtieL: int, t: int):
    cmd = (
        f'parallel -j {t} '
        f'"bowtie2 -x {tdna} -U {aligned_dir}/{{1}}Unmapped.out.mate1 -L {bowtieL} --local | '
        f'awk \'\\$3!=\\"*\\"\' | tee {tmp_dir}/{{1}}.Discordant.tdna.sam | cut -f 1 | sort -k1,1 | uniq > '
        f'{tmp_dir}/{{1}}.Discordant.tdna.readid" ::: $(echo {" ".join(prefixes)})'
    )
    run_cmd(cmd, shell=True)

def filter_discordant_pairs_parallel(discordant_prefixes: List[str], aligned_dir: Path, tmp_dir: Path, tag: str, t: int):
    cmd = (
        f'parallel -j {t} '
        f'"samtools view -h -N {tmp_dir}/{{1}}-F.Discordant.tdna.readid {aligned_dir}/{{1}}-R{tag} > {tmp_dir}/{{1}}-R.Discordant.genomic.sam && '
        f'samtools markdup -r {tmp_dir}/{{1}}-R.Discordant.genomic.sam {tmp_dir}/{{1}}-R.Discordant.genomic.rmdups.sam" ::: $(echo {" ".join(discordant_prefixes)}) &&'
        f'parallel -j {t} '
        f'"samtools view -h -N {tmp_dir}/{{1}}-R.Discordant.tdna.readid {aligned_dir}/{{1}}-F{tag} > {tmp_dir}/{{1}}-F.Discordant.genomic.sam && '
        f'samtools markdup -r {tmp_dir}/{{1}}-F.Discordant.genomic.sam {tmp_dir}/{{1}}-F.Discordant.genomic.rmdups.sam" ::: $(echo {" ".join(discordant_prefixes)})'
    )
    run_cmd(cmd, shell=True)

def filter_headerlines(file: Path):
    with open(file, 'r') as f:
        lines = [l for l in f if not l.startswith('@')]
    return StringIO(''.join(lines))

def get_discordant_pairs(tmp_dir: Path, prefix: str, discordant_prefixes: list[str], softclip_df: pd.DataFrame, control: bool):
    df_list = []
    cols_gen = [0,1,2,3,5,9]
    cols_tdna = [0,1,2,3,5,9,17,18]
    for p in discordant_prefixes:
        try:
            genomicF = pd.read_csv(filter_headerlines(f'{tmp_dir}/{p}-F.Discordant.genomic.rmdups.sam'), sep='\t', header=None, usecols=cols_gen, names=['readid', 'flag', 'chrom', 'start', 'cigar', 'seq'])
            tdnaF = pd.read_csv(filter_headerlines(f'{tmp_dir}/{p}-R.Discordant.tdna.sam'), sep='\t', header=None, usecols=cols_tdna, names=['readid', 'tdna_flag', 'tdna', 'tdna_pos', 'tdna_cigar', 'tdna_seq', 'MD1', 'MD2'])
            mergedF = genomicF.merge(tdnaF, on='readid')
            mergedF['prefix'] = f'{p}-F'
        except FileNotFoundError:
            mergedF = pd.DataFrame(columns=['readid', 'flag', 'chrom', 'start', 'cigar', 'seq', 'tdna_flag', 'tdna', 'tdna_pos', 'tdna_cigar', 'tdna_seq', 'MD1', 'MD2', 'prefix'])
        try:
            genomicR = pd.read_csv(filter_headerlines(f'{tmp_dir}/{p}-R.Discordant.genomic.rmdups.sam'), sep='\t', header=None, usecols=cols_gen, names=['readid', 'flag', 'chrom', 'start', 'cigar', 'seq'])
            tdnaR = pd.read_csv(filter_headerlines(f'{tmp_dir}/{p}-F.Discordant.tdna.sam'), sep='\t', header=None, usecols=cols_tdna, names=['readid', 'tdna_flag', 'tdna', 'tdna_pos', 'tdna_cigar', 'tdna_seq', 'MD1', 'MD2'])
            mergedR = genomicR.merge(tdnaR, on='readid')
            mergedR['prefix'] = f'{p}-R'
        except FileNotFoundError:
            mergedR = pd.DataFrame(columns=['readid', 'flag', 'chrom', 'start', 'cigar', 'seq', 'tdna_flag', 'tdna', 'tdna_pos', 'tdna_cigar', 'tdna_seq', 'MD1', 'MD2', 'prefix'])
        combined = pd.concat([mergedF, mergedR], ignore_index=True).sort_values(['chrom', 'start', 'cigar'])
        combined = combined.drop_duplicates(subset=['chrom', 'start', 'cigar', 'seq'])
        df_list.append(combined)
    discordant = pd.concat(df_list, ignore_index=True).sort_values(['chrom', 'start', 'cigar', 'prefix'])

    if not discordant.empty:
        discordant['end'] = discordant.apply(lambda row: get_mapped_position(int(row['start']), row['cigar']), axis=1)
        discordant['clip_pos'] = 'D'
        discordant['tdna_flag'] = '.'
        discordant['matchedL'] = discordant.apply(lambda row: get_matched_len_df(row['MD1'], row['MD2']), axis=1)
        uniq = discordant.merge(
            softclip_df[['chrom', 'start', 'end', 'prefix']], on=['chrom', 'start', 'end', 'prefix'], how='left', indicator=True
        )
        discordant_uniq = uniq[uniq['_merge']=='left_only'].drop(columns=['_merge'])
        discordant_uniq = discordant_uniq[['chrom', 'start', 'end', 'readid', 'flag', 'cigar', 'clip_pos', 'seq', 'tdna', 'tdna_flag', 'tdna_pos', 'tdna_cigar', 'matchedL', 'tdna_seq', 'prefix']]
    else:
        discordant_uniq = pd.DataFrame(columns = ['chrom', 'start', 'end', 'readid', 'flag', 'cigar', 'clip_pos', 'seq', 'tdna', 'tdna_flag', 'tdna_pos', 'tdna_cigar', 'matchedL', 'tdna_seq', 'prefix'])
    out_name = f"{prefix}_control.Discordant.tmp" if control else f"{prefix}.Discordant.tmp"
    discordant_uniq.to_csv(tmp_dir / out_name, sep='\t', index=False)
    return discordant_uniq

def define_tis_dict(tis_dict: Dict[str, List[str]], tdna_pos_dict: Dict[str, str], tdna_score_dict: Dict[str, np.ndarray], tis_cigar_dict: Dict[str, int], 
                    chrom: str, start: int, end: int, cigar: str, clip_pos: str, ori: str, tdna_junc: str):
    if clip_pos == "R":
        tis = f'{chrom}::{end}::{clip_pos}::{ori}'
        tis_cigar_dict[tis] = int(cigar.split('M')[-1].split('S')[0])
    else:
        tis = f'{chrom}::{start}::{clip_pos}::{ori}'
        tis_cigar_dict[tis] = int(cigar.split('S')[0])
    tis_dict[tis] = [tis]
    tdna_pos_dict[tis] = tdna_junc
    tdna_score_dict[tis] = np.zeros(5, dtype=int)
    return tis_dict, tdna_pos_dict, tdna_score_dict, tis_cigar_dict

def calculateDist(pos1:str, pos2:str):
    def parse(p: str):
        if "(" in p:
            a, b = p.split("(", 1)
            return [int(a), int(b.rstrip(")"))]
        return [int(p)]
    
    vals1, vals2 = parse(pos1), parse(pos2)
    return min(abs(a - b) for a in vals1 for b in vals2)

def soft_clip_scoring(softclip_df: pd.DataFrame, window: int):
    tis_dict = dict() ; tdna_pos_dict = dict() ; tdna_score_dict = dict() ; tis_cigar_dict = dict() ; tis_raw_dict = dict()
    for _, row in softclip_df.iterrows():
        chrom, start, end, cigar, clip_pos, tdna, ori, tdna_pos = row.iloc[[0, 1, 2, 5, 6, 8, 9, 10]]
        tdna_junc = f'{tdna}::{tdna_pos}'

        if clip_pos == "R":
            tis = f'{chrom}::{end}::{clip_pos}::{ori}'
            cigar_len = int(cigar.split('M')[-1].split('S')[0])
        else:
            tis = f'{chrom}::{start}::{clip_pos}::{ori}'
            cigar_len = int(cigar.split('S')[0])
        
        # check previously defined TISs
        matched = False
        adjacent = False
        for key, members in tis_dict.items():
            base_tdna = tdna_pos_dict[key].split('::')[0]
            base_pos = tdna_pos_dict[key].split('::')[1]
            if tdna != base_tdna:
                continue
            if calculateDist(str(tdna_pos), base_pos) >= 150: ## 150 -- less than max read length
                continue
            
            if tis in members:
                original_cigar = tis_cigar_dict[tis]
                diff = abs(cigar_len - original_cigar)
                if diff == 0:   tdna_score_dict[key][2] += 1
                elif diff == 1: tdna_score_dict[key][1] += 1
                else:           tdna_score_dict[key][0] += 1
                matched = True
                tis_raw_dict.setdefault(tis, []).append(row)
                break

            genomic_dist = abs(int(tis.split('::')[1]) - int(key.split('::')[1]))
            if chrom == key.split('::')[0] and clip_pos == key.split('::')[2] and genomic_dist < window:    
                members.append(tis)
                tis_cigar_dict[tis] = cigar_len    
                tdna_score_dict[key][3] += 1
                adjacent = True
                tis_raw_dict.setdefault(tis, []).append(row)

        if not matched and not adjacent:
            define_tis_dict(tis_dict, tdna_pos_dict, tdna_score_dict, tis_cigar_dict, chrom, start, end, cigar, clip_pos, ori, tdna_junc)   
            tis_raw_dict.setdefault(tis, []).append(row)
    return tis_dict, tdna_pos_dict, tdna_score_dict, tis_raw_dict

def discordant_scoring(discordant_df: pd.DataFrame, tdna_pos_dict: Dict[str, str], tdna_score_dict: Dict[str, np.ndarray], tis_raw_dict: Dict[str, List[pd.Series]], window: int):
    for _, row in discordant_df.iterrows():
        chrom, start, end, tdna_pos = row.iloc[[0,1,2,10]]
        for key, junc in tdna_pos_dict.items():
            base_chr, base_pos = key.split('::')[0], key.split('::')[1]
            junction_dist = calculateDist(str(tdna_pos), junc.split('::')[1])
            if chrom == base_chr and int(base_pos) >= start - window and int(base_pos) <= end + window and junction_dist < window:
                tdna_score_dict[key][4] += 1
                tis_raw_dict[key].append(row)
    return tdna_score_dict, tis_raw_dict

def DiscordantScoring(discordantDF, TDNAposDic, TDNAscoreDic, TISrawDic, window):
    for idx, row in discordantDF.iterrows():
        chrom, start, end, tdna_pos = row.iloc[[0,1,2,10]]
        for key, value in TDNAposDic.items():
            junction_dist = calculateDist(str(tdna_pos), value.split('::')[1])
            if chrom == key.split('::')[0] and int(key.split('::')[1]) >= start - window and int(key.split('::')[1]) <= end + window and junction_dist < window:
                TDNAscoreDic[key][4] += 1
                TISrawDic[key].append(row)
    return TDNAscoreDic, TISrawDic

def Scoring(out_dir: Path, prefix: str, softclip_df: pd.DataFrame, window: int, discordant: bool, discordant_df: pd.DataFrame, weight: Tuple, 
            threshold: float, blacklist: Optional[Path], bed: Optional[Path]):
    tis_dict, tdna_pos_dict, tdna_score_dict, tis_raw_dict = soft_clip_scoring(softclip_df, window)
    LRDic = {'L': 'T-DNA:REF', 'R': 'REF:T-DNA'}
    tmp_dir = out_dir / "tmp"
    score_tmp = tmp_dir / f"{prefix}.score.tmp"
    if discordant:
        tdna_score_dict, tis_raw_dict = discordant_scoring(discordant_df, tdna_pos_dict, tdna_score_dict, tis_raw_dict, window)
    with open(score_tmp, 'w') as o1:
        for key, members in tis_dict.items():
            chrom = key.split('::')[0]
            start = min([int(i.split('::')[1]) for i in members])
            end = max([int(i.split('::')[1]) for i in members])
            lr = LRDic[key.split('::')[2]]
            orientation = key.split('::')[3]
            o1.write('\t'.join([chrom, str(start), str(end), ','.join(members), tdna_pos_dict[key], lr, orientation, ','.join(map(str, tdna_score_dict[key])), str(np.dot(tdna_score_dict[key], weight)), key, prefix + '\n']))
    
    processed_tmp = tmp_dir / f"{prefix}.score.processed.tmp"
    if not blacklist and not bed:
        finaldf = pd.read_csv(score_tmp, sep='\t', header=None, names=['chrom', 'start', 'end', 'TISs', 'T-DNA', 'junction', 'orientation', 'count', 'score', 'key', 'sample'])
        finaldf.insert(loc=4, column='gene', value='.')
    else:
        if blacklist and bed:
            cmd = (
                f"bedtools intersect -a {score_tmp} -b {blacklist} -v | bedtools intersect -a - -b {bed} -wao | datamash groupby 1,2,3,4,5,6,7,8,9,10,11 collapse 15 | "
                f"awk '{{print $1, $2, $3, $4, $12, $5, $6, $7, $8, $9, $10, $11}}' OFS='\t' > {processed_tmp}"
            )
        elif bed:
            cmd = (
                f"bedtools intersect -a {score_tmp} -b {bed} -wao | datamash groupby 1,2,3,4,5,6,7,8,9,10,11 collapse 15 | "
                f"awk '{{print $1, $2, $3, $4, $12, $5, $6, $7, $8, $9, $10, $11}}' OFS='\t' > {processed_tmp}"
            )
        elif blacklist:
            cmd = (
                f"bedtools intersect -a {score_tmp} -b {blacklist} -v | awk '{{print $1, $2, $3, $4, \".\", $5, $6, $7, $8, $9, $10, $11}}' OFS='\t' > {processed_tmp}"
            )
        run_cmd(cmd, shell=True)
        finaldf = pd.read_csv(processed_tmp, sep='\t', header=None, names=['chrom', 'start', 'end', 'TISs', 'gene', 'T-DNA', 'junction', 'orientation', 'count', 'score', 'key', 'sample'])
    
    finaldf.loc[finaldf['start'] == finaldf['end'], 'end'] += 1
    read_info_list = []
    for filtered_tis in finaldf.query(f'score >= {threshold}')['key'].values:
        members = tis_dict[filtered_tis]
        for tis in members:
            for row in tis_raw_dict[tis]:
                outrow = row.copy()
                outrow['TIS'] = tis
                read_info_list.append(outrow)
    read_info_df = pd.DataFrame(read_info_list)    
    return finaldf, read_info_df

def main():
    args = parse_args()
    out_dir = Path(args.outdir)
#    confirm_overwrite(out_dir, args.prefix)
    tmp_dir, aligned_dir = make_dir(out_dir, args.prefix)
    log_path = out_dir / f"{args.prefix}.log"
    logger = setup_logger(log_path)
    tag = '.sorted.bam' if args.wgs else 'Aligned.sortedByCoord.out.bam'

    if not Path(args.tdna).is_file() and Path(f'{args.tdna}.1.bt2').is_file():
        print(f'### No T-DNA file or its bowtie2 index detected.  : {args.tdna}')
        sys.exit()
    else:
        print(f'\nT-DNA : {args.tdna}')
    print(f'Log file: {log_path}')
    
    if args.fastq:
        fastq_files, transgenic_prefixes, control_prefixes = alignment_preprocessing(args.prefix, args.fastq, args.control, args.paired)
        prefixes = transgenic_prefixes+control_prefixes if args.control else transgenic_prefixes
        print(f'Mapping reads to the reference genome...')
        if args.wgs:
            bowtie2_parallel(fastq_files, prefixes, aligned_dir, args.genome, args.mapq, args.thread)            
        else:
            star_parallel(fastq_files, prefixes, aligned_dir, args.genome, args.ratio, args.intron, args.thread, args.gzip, args.gtf)
            index_bams(prefixes, aligned_dir, tag, args.thread)
            print(f'Alignment completed. Index BAM files...')
    elif args.bam:
        print(f'\nInput bam files are detected. Skip alignment step.')
        bam_files, transgenic_prefixes, control_prefixes, index = bam_preprocessing(aligned_dir, args.prefix, args.bam, args.control, tag, args.paired)
        prefixes = transgenic_prefixes+control_prefixes if args.control else transgenic_prefixes
        if not index:
            index_bams(prefixes, aligned_dir, tag, args.thread)
            print(f'\nIndex BAM files...')
    else:
        print('## Error: No proper input files are provided')

    ### Extracting T-DNA sequences from soft clipped reads ###
    print(f'\nExtracting T-DNA sequences from soft-clipped reads...')
    prepare_softclipped_sams(prefixes, aligned_dir, tmp_dir, tag, args.thread)
    extract_soft_clipped_seq_parallel(prefixes, tmp_dir, args.l1, args.thread)

    print(f'\nMapping soft-clipped reads to T-DNA...')
    map_softclip_to_tdna_parallel(prefixes, tmp_dir, args.tdna, args.bowtieL, args.thread)

    print(f'\nFiltering the chimeric reads...')
    filter_sam_parallel(prefixes, tmp_dir, args.l1, args.d, args.l2, args.m, args.thread)
    softclip_df = concat_replicates(args.prefix, transgenic_prefixes, tmp_dir, args.paired, False)
    if args.control:
        softclip_df_ctl = concat_replicates(args.prefix, control_prefixes, tmp_dir, args.paired, True)
    
    # Discordant pairs
    if args.discordant and args.fastq:
        print(f'\n--discordant option specified. Searching for discordant pairs')
        align_discordant_pairs_parallel(prefixes, aligned_dir, tmp_dir, args.tdna, args.bowtieL, args.thread)
        discordant_prefixes = [pfx[:-2] for pfx in transgenic_prefixes[::2]]
        filter_discordant_pairs_parallel(discordant_prefixes, aligned_dir, tmp_dir, tag, args.thread)
        discordant_df = get_discordant_pairs(tmp_dir, args.prefix, discordant_prefixes, softclip_df, False)
        if args.control:
            discordant_prefixes_ctl = [discordant[:-2] for discordant in control_prefixes[::2]]
            discordant_df_ctl = get_discordant_pairs(tmp_dir, args.prefix, discordant_prefixes_ctl, softclip_df_ctl, True)
        print(f'\nEstimating the confidence score...')
        out_df, read_info_df = Scoring(out_dir, args.prefix, softclip_df, args.window, args.discordant, discordant_df, args.weight, args.threshold, args.blacklist, args.bed)
        if args.control and not softclip_df_ctl.empty:
            ctl_out_df, ctl_read_info_df = Scoring(out_dir, f'{args.prefix}_control', softclip_df_ctl, args.window, args.discordant, discordant_df_ctl, args.weight, args.threshold, args.blacklist, args.bed)
            out_df = pd.concat([out_df, ctl_out_df], ignore_index=True)
            read_info_df = pd.concat([read_info_df, ctl_read_info_df], ignore_index=True)
    else:
        if args.bam:
            print(f'--discordant option is not compatible with BAM input. The analysis will use only soft-clipped reads')
        print(f'\nEstimating the confidence score...')
        out_df, read_info_df = Scoring(out_dir, args.prefix, softclip_df, args.window, args.discordant, False, args.weight, args.threshold, args.blacklist, args.bed)
        if args.control and not softclip_df_ctl.empty:
            ctl_out_df, ctl_read_info_df = Scoring(out_dir, f'{args.prefix}_control', softclip_df_ctl, args.window, args.discordant, False, args.weight, args.threshold, args.blacklist, args.bed)
            out_df = pd.concat([out_df, ctl_out_df], ignore_index=True)
            read_info_df = pd.concat([read_info_df, ctl_read_info_df], ignore_index=True)
    out_file = out_dir / f'{args.prefix}.TDNA.txt'
    out_df.query(f'score >= {args.threshold}').to_csv(out_file, sep='\t', index=False)
    out_df.query(f'score < {args.threshold}').to_csv(f'{out_dir}/{args.prefix}.TDNA.discarded.txt', sep='\t', index=False)
    read_info = out_dir / f'{args.prefix}.read_info.txt'
    read_info_df.to_csv(read_info, sep='\t', index=False)
    
    ## Draw module
    if not out_df.query(f'score >= {args.threshold}').empty:
        print(f'\nDrawing the TIS plots...')
        track_dir, region_list, tis_dict, tis_bed = rangeBed(out_dir, out_file, args.prefix, args.bed)
        make_tdna_bam_bed(track_dir, args.prefix, tis_dict, read_info_df)        
        bam_files = [aligned_dir / f'{p}{tag}' for p in prefixes]
        bam_to_bigwig_parallel(track_dir, region_list, bam_files, prefixes, args.thread)
        if args.paired:
            prefixes = [p[:-2] for p in prefixes[::2]]
        bigwig_to_bedgraph_parallel(region_list, track_dir, prefixes, args.paired, args.thread)
        
        gtf = Path(args.gtf).absolute() if args.gtf else None
        bed = Path(args.bed).absolute() if args.bed else None
        create_track(track_dir, args.prefix, prefixes, tis_dict, Path(tis_bed).absolute(), gtf, bed)
        plot_image(track_dir, out_dir, args.prefix, region_list, tis_dict)
    else:
        print(f'\n### No T-DNA insertion sites were identified.\n')
    if args.keeptmp == False:
        run_cmd(f'rm -r {tmp_dir}', shell=True)

if __name__ == '__main__':
    main()
