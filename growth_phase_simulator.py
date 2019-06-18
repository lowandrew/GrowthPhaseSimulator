#!/usr/bin/env python

from Bio import SeqIO
import numpy as np
import subprocess
import tempfile
import argparse
import logging
import shutil
import os


def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--fasta',
                        type=str,
                        required=True,
                        help='Fasta file to simulate reads from. IMPORTANT: this should be a completed genome - using '
                             'a draft assembly will not work!')
    parser.add_argument('-d', '--depth',
                        type=float,
                        default=50.0,
                        help='Average depth to simulate reads to.')
    parser.add_argument('-r', '--ratio',
                        type=float,
                        default=2.0,
                        help='Ratio of depth at origin of replication to terminus.')
    parser.add_argument('-b', '--block_size',
                        type=int,
                        default=50000,
                        help='Number of bases to include in a block that will have the same coverage.')
    parser.add_argument('--origin_fasta',
                        default='oriC_ecoli.fasta',
                        help='Path to FASTA-formatted file with origin of replication in it.')
    parser.add_argument('-o', '--output_name',
                        help='Base name of your output FASTQ file. You\'ll end up with output_name_R1.fastq and '
                             'output_name_R2.fastq')
    parser.add_argument('-s', '--spikiness',
                        default=5,
                        help='Amount of standard deviation to add to coverage level for each block. Makes data end '
                             'up looking more spiky (to use a technical term) and hopefully therefore more realisitic.')
    return parser.parse_args()


def find_origin_of_replication(origin_fasta_file, ref_fasta_file):
    # TODO: Make sure that we aren't finding an origin of replication on a plasmid. That would be bad.
    origin_midpoint = None
    with tempfile.TemporaryDirectory() as tmpdir:
        blast_result_file = os.path.join(tmpdir, 'blast_results.tsv')
        tmpfasta = os.path.join(tmpdir, os.path.split(ref_fasta_file)[1])
        # Copy fasta file to temporary storage and make blast database, then blast
        shutil.copy(ref_fasta_file, tmpfasta)
        cmd = 'makeblastdb -dbtype nucl -in {}'.format(tmpfasta)
        subprocess.call(cmd, shell=True)
        cmd = 'blastn -query {query} -db {db} -out {out} -outfmt 6'.format(query=origin_fasta_file,
                                                                           db=tmpfasta,
                                                                           out=blast_result_file)
        subprocess.call(cmd, shell=True)
        # With BLAST done, extract the top hit.
        with open(blast_result_file) as f:
            for line in f:
                x = line.rstrip().split()
                chromosome = x[1]
                percent_id = float(x[2])
                start_position = int(x[8])
                end_position = int(x[9])
                if percent_id >= 88:  # Super arbitrary cutoff! Might actually need to look into this at some point.
                    origin_midpoint = int((start_position + end_position)/2)
                if origin_midpoint is not None:
                    break
    return origin_midpoint


if __name__ == '__main__':
    # TODO: Dependency check!
    logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    args = get_args()
    # Steps to make this work:
    # 1) Figure out where in the supplied genome the origin of replication is found.
    origin_midpoint = find_origin_of_replication(args.origin_fasta, args.fasta)
    if origin_midpoint is None:
        logging.error('ERROR: No origin of replication found. Cannot continue.')
        quit(code=1)
    logging.info('Origin of replication found at {}'.format(origin_midpoint))
    # 1.5) Create rotated version of FASTA with origin at zero coordinate.
    with tempfile.TemporaryDirectory() as tmpdir:
        rotated_fasta = os.path.join(tmpdir, 'rotated_fasta.fasta')
        sequence = SeqIO.read(args.fasta, 'fasta')
        sequence.seq = sequence.seq[:origin_midpoint] + sequence.seq[origin_midpoint:]
        sequence_length = len(sequence.seq)
        SeqIO.write(sequence, rotated_fasta, 'fasta')
        rotated_sequence = SeqIO.read(rotated_fasta, 'fasta')
        # 2) Do some math to figure out how much coverage we're doing for each block, based on the number of blocks
        # in the genome
        num_blocks = int(sequence_length/args.block_size) + 1
        mean_ratio = (args.ratio + 1)/2
        highest_level = args.depth * (args.ratio/mean_ratio)
        lowest_level = args.depth * (1/mean_ratio)
        step_size = (highest_level - lowest_level)/(num_blocks/2)
        current_coverage = highest_level
        # 3) Run ART on all the blocks
        for i in range(int(num_blocks/2) + 1):
            tmpfasta = os.path.join(tmpdir, 'tmpfasta.fasta')
            with open(tmpfasta, 'w') as f:
                f.write('>sequence\n{}\n'.format(rotated_sequence.seq[i * args.block_size: i * args.block_size + args.block_size]))
            cmd = 'art_illumina -ss MSv1 -na -i {tmpfasta} -l 250 -f {coverage} -m 400 -s 10 -o {filename}'.format(tmpfasta=tmpfasta,
                                                                                                                   coverage=np.random.normal(current_coverage, args.spikiness),
                                                                                                                   filename=os.path.join(tmpdir, 'block_' + str(i) + '_'))
            subprocess.call(cmd, shell=True)
            current_coverage -= step_size
        for j in range(int(num_blocks/2) + 1):
            tmpfasta = os.path.join(tmpdir, 'tmpfasta.fasta')
            with open(tmpfasta, 'w') as f:
                f.write('>sequence\n{}\n'.format(rotated_sequence.seq[(i + j) * args.block_size: (i + j) * args.block_size + args.block_size]))
            cmd = 'art_illumina -ss MSv1 -na -i {tmpfasta} -l 250 -f {coverage} -m 400 -s 10 -o {filename}'.format(tmpfasta=tmpfasta,
                                                                                                                   coverage=np.random.normal(current_coverage, args.spikiness),
                                                                                                                   filename=os.path.join(tmpdir, 'block_' + str(i + j) + '_'))
            subprocess.call(cmd, shell=True)
            current_coverage += step_size

        # 4) Combine all the FASTQ files created from 3 into one big FASTQ that has the coverage pattern that we want.
        cmd = 'cat {} > {}_R1.fastq'.format(os.path.join(tmpdir, '*_1.fq'), args.output_name)
        os.system(cmd)
        cmd = 'cat {} > {}_R2.fastq'.format(os.path.join(tmpdir, '*_2.fq'), args.output_name)
        os.system(cmd)
