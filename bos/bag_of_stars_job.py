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
import numpy as np
import pandas as pd
import subprocess as sp
from pathlib import Path
import time
import pysam


def printfl(*args):
    '''Flushed printing to stdout'''
    return print(*args, flush=True)


if __name__ == '__main__':

    pa = argparse.ArgumentParser(description='STAR + htseq-count in bags')
    pa.add_argument(
            '--group-name',
            required=True,
            help='Name of the group')
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
            default=None,
            help='Folder with the STAR genome hash')
    pa.add_argument(
            '--htseq',
            action='store_true',
            help='Run htseq-count at the end')
    pa.add_argument(
            '--chain-htseq',
            default=None,
            help='Chain htseq counts on these subfolders at the end of the every single htseq-count call')
    pa.add_argument(
            '--annotationFile',
            default=None,
            help='File with the GTF feature annotations for htseq')
    args = pa.parse_args()

    if args.htseq and (args.annotationFile is None):
        raise ValueError('To run htseq-count, specify --annotationFile')

    if args.chain_htseq and (args.annotationFile is None):
        raise ValueError('To run htseq-count, specify --annotationFile')

    args.output = args.output.rstrip('/')+'/'
    args.samplenames = args.samplenames.split(' ')
    args.fastqs = args.fastqs.split(' ')
    group = args.samplenames
    group_fastq = list(zip(args.fastqs[::2], args.fastqs[1::2]))
    output_tmp_star = args.output+'STAR_TMP_group_{:}/'.format(args.group_name)

    printfl('Remove previous STAR temp dir')
    if not args.dry:
        sp.run('rm -rf {:}'.format(output_tmp_star),
               check=True,
               shell=True)

    printfl('Load genome into memory')
    call = [
        os.getenv('STAR', 'STAR'),
        '--genomeDir', args.genomeDir,
        '--genomeLoad', 'LoadAndExit',
        '--outTmpDir', output_tmp_star,
        ]
    printfl(' '.join(call))
    if not args.dry:
        sp.run(call, check=True) 

    try:
        for sn, (fq1, fq2) in zip(group, group_fastq):
            printfl('Remove previous STAR temp dir')
            if not args.dry:
                sp.run('rm -rf {:}'.format(output_tmp_star),
                       check=True,
                       shell=True)

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
                '--outTmpDir', output_tmp_star,
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
        printfl('Remove previous STAR temp dir')
        if not args.dry:
            sp.run('rm -rf {:}'.format(output_tmp_star),
                   check=True,
                   shell=True)

        printfl('Remove genome from memory')
        call = [
            os.getenv('STAR', 'STAR'),
            '--genomeDir', args.genomeDir,
            '--genomeLoad', 'Remove',
            '--outTmpDir', output_tmp_star,
            ]
        printfl(' '.join(call))
        if not args.dry:
            try:
                sp.run(call, check=True)
            except sp.CalledProcessError as e:
                print(e)
                print('Moving on')

        printfl('Remove previous STAR temp dir')
        if not args.dry:
            sp.run('rm -rf {:}'.format(output_tmp_star),
                   check=True,
                   shell=True)

    if args.htseq is not None:
        print('STAR mappings done, check output BAM files')
        has_failed = False
        mapped_fns = [args.output+sn+'/Aligned.out.bam' for sn in group]
        for mapped_fn in mapped_fns:
            print(os.path.dirname(mapped_fn)+'...', end='', flush=True)
            try:
                with pysam.AlignmentFile(mapped_fn, 'rb') as bamfile:
                    for read in bamfile:
                        break
                print('OK', flush=True)
            except IOError:
                has_failed = True
                print('Failed!', flush=True)
        if has_failed:
            raise IOError('One or more BAM files failed to map')

        print('Call htseq-count')
        htseq_fn = args.output+'counts_group_{:}.tsv'.format(args.group_name)
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
                fout.write('\t'.join(['feature'] + group)+'\n')
                fout.write(output)

    if args.chain_htseq is not None:
        n_groups = int(args.chain_htseq)

        print('Wait for all STAR + htseq-count from other jobs')
        flag_fns = [args.output+'counts_group_{:}.tsv'.format(ig+1) for ig in range(n_groups)]
        jobs_done = [False for ig in range(n_groups)]
        is_first = True
        while not all(jobs_done):
            for isn, flag_fn in enumerate(flag_fns):
                if os.path.isfile(flag_fn):
                    jobs_done[isn] = True

            if not is_first:
                time.sleep(60)
                is_first = False

        print('All STAR +htseq-count, chaining counts')
        counts_fns = [args.output+'counts_group_{:}.tsv'.format(ig+1) for ig in range(n_groups)]
        htseq_fn = args.output+'counts.tsv'
        if not args.dry:
            samplenames = []
            n_samples_max = n_groups * len(group)
            first_group = True
            for ig, counts_fn in enumerate(counts_fns):
                print('Reading counts for group {:}'.format(ig+1))
                counts = pd.read_csv(counts_fn, sep='\t', index_col=0)
                if first_group:
                    first_group = False
                    n_features = counts.shape[0]
                    features = counts.index
                    counts_all = np.zeros(
                            (n_features, n_samples_max),
                            dtype=np.float32,
                            )
                print('Chaining counts for group {:}'.format(ig+1))
                counts_all[:, len(samplenames): len(samplenames) + counts.shape[1]] = counts.values.astype(np.float32)
                samplenames.extend(counts.columns.tolist())

            print('Merging all into a dataframe')
            counts_all = pd.DataFrame(
                    data=counts_all[:, :len(samplenames)],
                    index=features,
                    columns=samplenames,
                    )
            print('Writing all counts to file: {:}'.format(htseq_fn))
            counts_all.to_csv(htseq_fn, sep='\t', index=True)

