# vim: fdm=indent
'''
author:     Fabio Zanini
date:       25/05/18
content:    Call STAR on a slurm cluster without Snakemake
'''
import os
import sys
import glob
import argparse
import subprocess as sp
from pathlib import Path
import pysam



if __name__ == '__main__':

    pa = argparse.ArgumentParser(description='STAR mapping in bags')
    pa.add_argument(
        'fastq_folder',
        help='Parent folder of subfolders with 2 fastq.gz files in each.')
    pa.add_argument(
        '--dry',
        action='store_true',
        help='Dry run')
    pa.add_argument(
        '--output',
        required=True,
        help='Parent folder for the output. For each input subfolder, an output subfolder will be made')
    pa.add_argument(
            '-n',
            type=int,
            default=40,
            help='Number of samples per STAR call')
    pa.add_argument(
            '--genomeDir',
            default=None,
            help='Folder with the STAR genome hash')
    pa.add_argument(
            '--local',
            action='store_true',
            help='Do not send to cluster, do everything locally')
    pa.add_argument(
            '--cpus-per-task', type=int, default=6,
            help='Number of CPUs for each STAR call',
            )
    pa.add_argument(
            '--mem', type=int, default=32000,
            help='RAM memory in MB for each STAR call',
            )
    pa.add_argument(
            '--time', default='2-0:0:0',
            help='Time limit on each group of STAR jobs (see slurm docs for format info)',
            )
    pa.add_argument(
            '--htseq',
            action='store_true',
            help='Call htseq-count at the end of the STAR mapping',
            )
    pa.add_argument(
            '--annotationFile',
            default=None,
            help='File with the GTF feature annotations for htseq')
    pa.add_argument(
            '--delete-empty-BAM',
            action='store_true',
            help='Delete Aligned.out.bam if it contains zero reads')
    args = pa.parse_args()

    print('Find bag_of_stars_job.py')
    job_fn = os.path.dirname(os.path.abspath(__file__))+'/bag_of_stars_job.py'
    if not os.path.isfile(job_fn):
        raise IOError('bag_of_stars_job.py not found')

    print('Check input arguments')
    if args.htseq and (args.annotationFile is None):
        raise ValueError('To run htseq-count, specify --annotationFile')

    print('Make output root folder')
    args.output = args.output.rstrip('/')+'/'
    if not args.dry:
        os.makedirs(args.output, exist_ok=True)

    if not args.local:
        print('Make job log folder')
        log_fdn = args.output+'logs/'
        if not args.dry:
            os.makedirs(log_fdn, exist_ok=True)

    print('Scan input folder')
    args.fastq_folder = args.fastq_folder.rstrip('/')+'/'
    samplenames = os.listdir(args.fastq_folder)
    if len(samplenames) == 0:
        raise IOError('No input subfolders')
    fn_fastqs = []
    for sn in samplenames:
        fn_r1 = glob.glob(args.fastq_folder+sn+'/*R1_001.fastq.gz')
        fn_r2 = glob.glob(args.fastq_folder+sn+'/*R2_001.fastq.gz')
        if (len(fn_r1) != 1) or (len(fn_r2) != 1):
            raise IOError('Sample does not have two decent fastqs: {:}'.format(sn))
        fn_fastqs.append((fn_r1[0], fn_r2[0]))

    print('Make output subfolder')
    if not args.dry:
        for sn in samplenames:
            os.makedirs(args.output+sn, exist_ok=True)

    print('Split input samples into batches of {:} samples each'.format(args.n))
    n = args.n
    n_samples = len(samplenames)
    n_groups = ((n_samples - 1) // args.n) + 1
    groups = [samplenames[i * n: (i+1) * n] for i in range(n_groups)]
    group_fastqs = [fn_fastqs[i * n: (i+1) * n] for i in range(n_groups)]
    print('{:} groups'.format(n_groups))

    print('Run STAR')
    for ig, (group, group_fastq) in enumerate(zip(groups, group_fastqs)):
        groupname = 'group_{:}'.format(ig+1)
        if args.local:
            print('Group {:}, load genome into memory'.format(ig+1))
            call = [
                os.getenv('STAR', 'STAR'),
                '--genomeDir', args.genomeDir,
                '--genomeLoad', 'LoadAndExit',
                ]
            print(' '.join(call))
            if not args.dry:
                sp.run(call, check=True)

            try:
                for sn, (fq1, fq2) in zip(group, group_fastq):
                    print('Mapping sample {:}'.format(sn))
                    flag_fn = args.output+sn+'/STAR.done'
                    if os.path.isfile(flag_fn):
                        print('Flag file found, skipping sample')
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
                    print(' '.join(call))
                    if not args.dry:
                        sp.run(call, check=True)
                        Path(flag_fn).touch()

            finally:
                print('Group {:}, remove genome from memory'.format(ig+1))
                call = [
                    os.getenv('STAR', 'STAR'),
                    '--genomeDir', args.genomeDir,
                    '--genomeLoad', 'Remove',
                    ]
                print(' '.join(call))
                if not args.dry:
                    sp.run(call, check=True)

            if args.htseq:
                print('STAR mappings done, check output BAM files')
                mapped_fns = [args.output+sn+'/Aligned.out.bam' for sn in samplenames]
                has_failed = False
                for mapped_fn in mapped_fns:
                    print(os.path.dirname(mapped_fn)+'...', end='', flush=True)
                    try:
                        with pysam.AlignmentFile(mapped_fn, 'rb') as bamfile:
                            for read in bamfile:
                                break
                            else:
                                raise IOError('Zero reads in {:}'.format(mapped_fn))
                        print('OK', flush=True)
                    except IOError:
                        has_failed = True
                        print('Failed!', flush=True)
                        if args.delete_empty_BAM:
                            print('Remove file: {:}'.format(mapped_fn), flush=True)
                            os.remove(mapped_fn)
                if has_failed:
                    raise IOError('One or more BAM files failed to map')

                print('Call htseq-count')
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

        else:
            job_out_fn = log_fdn+'{:}.out'.format(groupname)
            job_err_fn = log_fdn+'{:}.err'.format(groupname)
            fastqs = []
            for gfq in group_fastq:
                fastqs.extend(list(gfq))
            call = [
                'sbatch',
                '-o', job_out_fn,
                '-e', job_err_fn,
                '-N', '1',
                '--job-name', 'bos-{:}'.format(groupname),
                '--cpus-per-task={:}'.format(args.cpus_per_task),
                '--mem={:}'.format(args.mem),
                '--time={:}'.format(args.time),
                '--partition=quake,hns,normal',
                job_fn,
                '--group-name', str(ig+1),
                '--output', args.output,
                '--genomeDir', args.genomeDir,
                '--samplenames', ' '.join(group),
                '--fastqs', ' '.join(fastqs),
                ]
            if args.dry:
                call.append('--dry')
            if args.delete_empty_BAM:
                call.append('--delete-empty-BAM')
            if args.htseq:
                call.append('--htseq')
                call.extend(['--annotationFile', args.annotationFile])
            if args.htseq and (ig == 0):
                call.extend(['--chain-htseq', str(len(groups))])
            print(' '.join(call))
            if not args.dry:
                sp.run(call, check=True)
