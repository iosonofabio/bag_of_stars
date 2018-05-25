![Logo](logo.png)
# Bag of STARs
STAR mapping on slurm clusters in bags.

## Requirements
- `python3.4+`
- `STAR` in your `PATH`

## Installation
```bash
git clone https://github.com/iosonofabio/bag_of_stars.git
```

## Usage
Call bag of stars from the `bos` folder:
```bash
python bag_of_stars.py --genomeDir <your genome folder> --output <your output folder> <your fastq root folder>
```
- The fastq root folder must contain subfolders with each a read1 and read2 file.
- New subfolders with the same names will be made inside the output folder.
- The genome folder must contain STAR's hash files (e.g. `SA`, `SAindex`, `Genome`)
- `STAR` must be in your `PATH`, you can check your `.bashrc` for what folders are there.

## Help
```bash
python bag_of_stars.py --help

usage: bag_of_stars.py [-h] [--dry] --output OUTPUT [-n N]
                       [--genomeDir GENOMEDIR] [--local]
                       [--cpus-per-task CPUS_PER_TASK] [--mem MEM]
                       [--time TIME]
                       fastq_folder

STAR mapping in bags

positional arguments:
  fastq_folder          Parent folder of subfolders with 2 fastq.gz files in
                        each.

optional arguments:
  -h, --help            show this help message and exit
  --dry                 Dry run
  --output OUTPUT       Parent folder for the output. For each input
                        subfolder, an output subfolder will be made
  -n N                  Number of samples per STAR call
  --genomeDir GENOMEDIR
                        Folder with the STAR genome hash
  --local               Do not send to cluster, do everything locally
  --cpus-per-task CPUS_PER_TASK
                        Number of CPUs for each STAR call
  --mem MEM             RAM memory in MB for each STAR call
  --time TIME           Time limit on each group of STAR jobs (see slurm docs
                        for format info)
```
