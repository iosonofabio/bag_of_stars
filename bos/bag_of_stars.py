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


class BagOfStars(object):
    def __init__(self, args):
        self.args = args

    def check_args(self):
        if self.args.htseq and (self.args.annotationFile is None):
            raise ValueError('To run htseq-count, specify --annotationFile')
        if self.args.fastqFolder is not None:
            self.args.fastqFolder = self.args.fastqFolder.rstrip('/')+'/'
        self.args.output = self.args.output.rstrip('/')+'/'

    def find_bos_job_script(self):
        job_fn = os.path.dirname(os.path.abspath(__file__))+'/bag_of_stars_job.py'
        if not os.path.isfile(job_fn):
            raise IOError('bag_of_stars_job.py not found')
        self.job_fn = job_fn

    def scan_input_folder(self):
        if self.args.fastqFolder is not None:
            samplenames = os.listdir(self.args.fastqFolder)
            if len(samplenames) == 0:
                raise IOError('No input subfolders')
            fn_fastqs = []
            for sn in samplenames:
                fn_r1 = glob.glob(self.args.fastqFolder+sn+'/*R1_001.fastq.gz')
                fn_r2 = glob.glob(self.args.fastqFolder+sn+'/*R2_001.fastq.gz')
                if (len(fn_r1) != 1) or (len(fn_r2) != 1):
                    raise IOError('Sample does not have two decent fastqs: {:}'.format(sn))
                fn_fastqs.append((fn_r1[0], fn_r2[0]))
            self.samplenames = samplenames
            self.fn_fastqs = fn_fastqs
        else:
            self.samplenames = os.listdir(self.args.output)
            self.fn_fastqs = None
            if len(self.samplenames) == 0:
                raise IOError('Could not find any samples in input or output folders')

    def make_output_folders(self):
        if not self.args.dry:
            os.makedirs(self.args.output, exist_ok=True)
        if not self.args.local:
            print('Make job log folder')
            self.log_fdn = self.args.output+'logs/'
            if not self.args.dry:
                os.makedirs(self.log_fdn, exist_ok=True)
        if not args.dry:
            for sn in self.samplenames:
                os.makedirs(self.args.output+sn, exist_ok=True)

    def split_into_groups(self):
        n = self.args.n
        n_samples = len(self.samplenames)
        n_groups = ((n_samples - 1) // args.n) + 1

        groups = []
        for i in range(n_groups):
            g = {'names': self.samplenames[i * n: (i+1) * n]}
            g['mapped'] = [self.args.output+sn+'/Aligned.out.bam' for sn in g['names']]
            if self.fn_fastqs is not None:
                g['fastqs'] = self.fn_fastqs[i * n: (i+1) * n]
            groups.append(g)

        self.n_groups = n_groups
        self.groups = groups

    def check_need_star_genome(self):
        for ig, g in enumerate(self.groups):
            needs_star = False
            for sn in g['names']:
                flag_fn = args.output+sn+'/STAR.done'
                if not os.path.isfile(flag_fn):
                    needs_star = True
                    break
            g['needs_star'] = needs_star
        if any(g['needs_star'] for g in self.groups) and (self.fn_fastqs is None):
            raise IOError('No fastq folder specified, but some samples need STAR')

    def star_load(self):
        call = [
            os.getenv('STAR', 'STAR'),
            '--genomeDir', self.args.genomeDir,
            '--genomeLoad', 'LoadAndExit',
            ]
        print(' '.join(call))
        if not self.args.dry:
            sp.run(call, check=True)

    def star_release(self):
        call = [
            os.getenv('STAR', 'STAR'),
            '--genomeDir', self.args.genomeDir,
            '--genomeLoad', 'Remove',
            ]
        print(' '.join(call))
        if not self.args.dry:
            sp.run(call, check=True)

    def star_map(self, group, verbose=True):
        for sn, (fq1, fq2) in zip(group['names'], group['fastqs']):
            if verbose:
                print('Mapping sample {:}'.format(sn))
            flag_fn = self.args.output+sn+'/STAR.done'
            if os.path.isfile(flag_fn):
                if verbose:
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
                '--outSAMattributes', 'NH', 'HI', 'AS', 'NM',
                '--outFilterMatchNminOverLread', '0.4',
                '--outFilterScoreMinOverLread', '0.4',
                '--clip3pAdapterSeq', 'CTGTCTCTTATACACATCT',
                '--outReadsUnmapped', 'Fastx',
                ]
            print(' '.join(call))
            if not self.args.dry:
                sp.run(call, check=True)
                Path(flag_fn).touch()

    def check_STAR_output_group(self, group):
        has_failed = False
        for mapped_fn in group['mapped']:
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
                if self.args.delete_empty_BAM:
                    print('Remove file: {:}'.format(mapped_fn), flush=True)
                    os.remove(mapped_fn)
        if has_failed:
            raise IOError('One or more BAM files failed to map')

    def htseq(self, group):
        htseq_fn = self.args.output+'counts.tsv'
        call = [
            os.getenv('HTSEQ-COUNT', 'htseq-count'),
            '--format', 'bam',
            '--mode', 'intersection-nonempty',
            '--stranded', 'no',
            '--secondary-alignments', 'ignore',
            '--supplementary-alignments', 'ignore',
            ] + group['mapped'] + [
            self.args.annotationFile,
            ]
        print(' '.join(call))
        if not self.args.dry:
            output = sp.run(call, check=True, stdout=sp.PIPE).stdout.decode()
            with open(htseq_fn, 'wt') as fout:
                fout.write('\t'.join(['feature'] + group['names'])+'\n')
                fout.write(output)

    def call_batch_job(self, group, groupname=''):
        job_out_fn = self.log_fdn+'{:}.out'.format(groupname)
        job_err_fn = self.log_fdn+'{:}.err'.format(groupname)
        fastqs = []
        for gfq in group['fastqs']:
            fastqs.extend(list(gfq))
        call = [
            'sbatch',
            '-o', job_out_fn,
            '-e', job_err_fn,
            '-N', '1',
            '--job-name', 'bos-{:}'.format(groupname),
            '--cpus-per-task={:}'.format(self.args.cpus_per_task),
            '--mem={:}'.format(self.args.mem),
            '--time={:}'.format(self.args.time),
            '--partition={:}'.format(self.args.partition),
            job_fn,
            '--group-name', str(ig+1),
            '--output', self.args.output,
            '--genomeDir', self.args.genomeDir,
            '--samplenames', "'{:}'".format(' '.join(group['names'])),
            '--fastqs', "'{:}'".format(' '.join(fastqs)),
            ]
        if self.args.dry:
            call.append('--dry')
        if self.args.delete_empty_BAM:
            call.append('--delete-empty-BAM')
        if self.args.htseq:
            call.append('--htseq')
            call.extend(['--annotationFile', self.args.annotationFile])
        if self.args.htseq and (ig == 0):
            call.extend(['--chain-htseq', str(len(self.groups))])
        print(' '.join(call))
        if not self.args.dry:
            sp.run(call, check=True)


if __name__ == '__main__':

    pa = argparse.ArgumentParser(description='STAR mapping in bags')
    pa.add_argument(
        '--fastqFolder',
        default=None,
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
    pa.add_argument(
            '--partition',
            default='quake,hns,normal',
            help='Partition to use if running on a cluster')
    args = pa.parse_args()

    bos = BagOfStars(args)

    print('Check input arguments')
    bos.check_args()

    print('Find bag_of_stars_job.py')
    bos.find_bos_job_script()
    job_fn = bos.job_fn

    print('Scan input folder')
    bos.scan_input_folder()
    samplenames, fn_fastqs = bos.samplenames, bos.fn_fastqs

    print('Make output folders')
    bos.make_output_folders()

    print('Split input samples into batches of {:} samples each'.format(args.n))
    bos.split_into_groups()
    print('{:} groups'.format(bos.n_groups))
    groups = bos.groups

    print('Check whether we still need STAR for each group')
    bos.check_need_star_genome()

    print('Run STAR')
    for ig, group in enumerate(groups):
        groupname = 'group_{:}'.format(ig+1)
        if args.local:
            if group['needs_star']:
                print('Group {:}, load genome into memory'.format(ig+1))
                bos.star_load()

                try:
                    bos.star_map(group)
                finally:
                    print('Group {:}, remove genome from memory'.format(ig+1))
                    bos.star_release()
            else:
                print('Group {:}, STAR not needed'.format(ig+1))

            if args.htseq:
                print('Check STAR output BAM files')
                bos.check_STAR_output_group(group)

                print('Call htseq-count')
                bos.htseq(group)

        else:
            bos.call_batch_job(group, groupname)
