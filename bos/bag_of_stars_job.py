#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       25/05/18
content:    Bag of stars, single bad (cluster job)
'''
import os
import sys
import argparse
import subprocess as sp
from pathlib import Path
import time


def printfl(*args):
    '''Flushed printing to stdout'''
    return print(*args, flush=True)


if __name__ == '__main__':

    pa = argparse.ArgumentParser(description='STAR + htseq-count in bags')
    pa.add_argument(
            '--samplenames',
            required=True,
            help='The samplenames')
    pa.add_argument(
            '--fastqs',
            required=True,
            help='The filenames with the fastqs (2 x number of samples)')
    pa.add_argument(
        '--dry',
        action='store_true',
        help='Dry run')
    pa.add_argument(
        '--output',
        required=True,
        help='Parent folder for the output. For each input subfolder, an output subfolder will be made')
    pa.add_argument(
            '--genomeDir',
            default='/oak/stanford/groups/quake/fzanini/postdoc/bag_of_stars/data/human_genome/STAR_DIR',
            help='Folder with the STAR genome hash')
    pa.add_argument(
            '--htseq',
            default=None,
            help='Call htseq-count on these subfolders at the end of the STAR mapping',
    pa.add_argument(
            '--annotationFile',
            default='/oak/stanford/groups/quake/fzanini/postdoc/bag_of_stars/data/human_genome/human_tiny_with_ERCC.gtf',
            help='File with the GTF feature annotations for htseq')
    args = pa.parse_args()

    args.output = args.output.rstrip('/')+'/'
    args.samplenames = args.samplenames.split(' ')
    args.fastqs = args.fastqs.split(' ')
    group = args.samplenames
    group_fastq = list(zip(args.fastqs[::2], args.fastqs[1::2]))

    printfl('Load genome into memory')
    call = [
        os.getenv('STAR', 'STAR'),
        '--genomeDir', args.genomeDir,
        '--genomeLoad', 'LoadAndExit',
        ]
    printfl(' '.join(call))
    if not args.dry:
        sp.run(call, check=True) 

    try:
        for sn, (fq1, fq2) in zip(group, group_fastq):
            printfl('Mapping sample {:}'.format(sn))
            flag_fn = args.output+sn+'/STAR.done'
            if os.path.isfile(flag_fn):
                printfl('Flag file found, skipping sample')
                continue
            call = [
                os.getenv('STAR', 'STAR'),
                '--runMode', 'alignReads',
                '--genomeDir', args.genomeDir,
                '--genomeLoad', 'LoadAndKeep',
                '--readFilesIn', fq1, fq2,
                '--readFilesCommand', 'zcat',
                '--outFilterType', 'BySJout',
                '--outFilterMultimapNmax', '20',
                '--alignSJoverhangMin', '8',
                '--alignSJDBoverhangMin', '1',
                '--outFilterMismatchNmax 999',
                '--outFilterMismatchNoverLmax 0.04',
                '--alignIntronMin', '20',
                '--alignIntronMax', '1000000',
                '--alignMatesGapMax', '1000000',
                '--outSAMstrandField', 'intronMotif',
                '--outFileNamePrefix', args.output+sn+'/',
                '--outSAMtype', 'BAM', 'Unsorted',
                '--outSAMattributes', 'NH', 'HI', 'AS', 'NM', 'MD',
                           '--outFilterMatchNminOverLread', '0.4',
                           '--outFilterScoreMinOverLread', '0.4',
                           '--clip3pAdapterSeq', 'CTGTCTCTTATACACATCT',
                           '--outReadsUnmapped', 'Fastx',
                           ]
            printfl(' '.join(call))
            if not args.dry:
                sp.run(call, check=True) 
                Path(flag_fn).touch()

    finally:
        printfl('Remove genome from memory')
        call = [
            os.getenv('STAR', 'STAR'),
            '--genomeDir', args.genomeDir,
            '--genomeLoad', 'Remove',
            ]
        printfl(' '.join(call))
        if not args.dry:
            sp.run(call, check=True) 

    if args.htseq is not None:
        samplenames = args.htseq.split()

        print('Wait for all STAR mapping from other jobs')
        flag_fns = [args.output+sn+'/STAR.done' for sn in samplenames]
        star_done = np.zeros(len(samplenames), bool)
        is_first = True
        while not star_done.all():
            for isn, flag_fn in enumerate(flag_fns):
                if os.path.isfile(flag_fn):
                    star_done[isn] = True
            
            if not is_first:
                time.sleep(60)
                is_first = False

        print('All STAR mappings done, call htseq-count')
        mapped_fns = [args.output+sn+'/Aligned.out.bam' for sn in samplenames]
        htseq_fn = args.output+'counts.tsv'
        call = [
            os.getenv('HTSEQ-COUNT', 'htseq-count'),
            '--format', 'bam',
            '--mode', 'intersection-nonempty',
            '--stranded', 'no',
            '--secondary-alignments', 'ignore',
            '--supplementary-alignments', 'ignore',
            ] + mapped_fns + [
            args.annotationFile,
            ]
        print(' '.join(call))
        if not args.dry:
            output = sp.run(call, check=True, stdout=sp.PIPE).stdout.decode()
            with open(htseq_fn, 'wt') as fout:
                fout.write('\t'.join(['feature'] + samplenames)+'\n')
                fout.write(output)
