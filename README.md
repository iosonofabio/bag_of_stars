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
```bash
python bag_of_stars.py --genomeDir <your STAR_DIR> --output <your output folder> <your fastq root folder>
```
- The fastq root folder must contain subfolders with each a read1 and read2 file.
- New subfolders with the same names will be made inside the output folder.

## Help
```bash
python bag_of_stars.py --help
```
