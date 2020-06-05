#!/usr/bin/env python

import sys
import os
import glob
import pandas as pd
import re
import argparse
from datetime import date
from tabulate import tabulate

def main():
    
    local_path = os.path.dirname(os.path.realpath(__file__))
    #print(local_path)
    data_path = f"{local_path}"
    canu_helper = f"{local_path}/cdhit_cluster_filter.R"
    now = date.today()

    cli = argparse.ArgumentParser()
    cli.add_argument('-i', '--InputFolder', help="Folder containing barcoded fastq", required=True,)
    cli.add_argument('-o', '--OutputFolder', help="Output Folder", required=False, default=f'~/rabies_results/output_{now}')

    cli.add_argument('--correctedErrorRate', help="Estimated Error Rate", type=float, required=False, default=0.15)
    cli.add_argument('--readCountMin', help="Minimum Read Counts in Canu to Keep", type=int, required=False, default=10)
    cli.add_argument('--genomeSize', help="Estimated Amplicon Size", required=False, type=int, default=2000)
    cli.add_argument('--verbose', help = "Keep Intermediate Files", required=False, action='store_true')
    cli.add_argument('--model', help="Basecall Model", required=False, type=str, default='r941_min_high_g303')
    cli.add_argument('--readSamplingCoverage', help="Max Sample Coverage", required=False, type=int, default=5000)
    cli.add_argument('--minReadLength', help="Minimum Read Length to Include", required=False, type=int, default=1200)
    args = cli.parse_args()

    files = sorted([f for f in glob.glob(args.InputFolder+"/**", recursive = True) if re.search(r'(.*)\.((fastq|fq)(|\.gz))$', f)])
    #InputFolder = os.path.expanduser(args.InputFolder)
    #files = sorted(glob.glob(InputFolder+"/*.fastq"))
    print(files)
    OutputFolder = os.path.expanduser(args.OutputFolder)
    os.environ["PERL5LIB"] = ""

    qc_dir= f"{OutputFolder}/qc"
    metric_dir= f"{OutputFolder}/metric"
    assembly_dir= f"{OutputFolder}/assembly"
    os.system(f"mkdir -p {qc_dir}")
    os.system(f"mkdir -p {metric_dir}")
    os.system(f"mkdir -p {assembly_dir}")
    f=open(f"{OutputFolder}/qc.log", 'w+')

    df = pd.DataFrame( columns = ['Sample', 'Reads', 'Mapped', 'Ncov', 'Gcov', 'ReadLength'])

    for i in range(0, len(files)):
        filec = files[i]

        base = os.path.splitext(os.path.basename(filec))[0]
        base = os.path.splitext(base)[0]
        #print(base)

        fastqc_cmd = f"fastqc {filec} -o {qc_dir}"
        f.write(fastqc_cmd+'\n')
        os.system(fastqc_cmd)

        minimap2_cmd = f"minimap2 -ax map-ont {data_path}/Americas2.fasta {filec} > {metric_dir}/{base}.sam "
        f.write(minimap2_cmd+'\n')
        os.system(minimap2_cmd)
        samtools_cmd1 = f"samtools view -S -b {metric_dir}/{base}.sam > {metric_dir}/{base}.bam"
        f.write(samtools_cmd1+'\n')
        os.system(samtools_cmd1)
        samtools_cmd2 = f"samtools sort {metric_dir}/{base}.bam -o {metric_dir}/{base}.sorted.bam"
        f.write(samtools_cmd2+'\n')
        os.system(samtools_cmd2)
        samtools_cmd3 = f"samtools index {metric_dir}/{base}.sorted.bam"
        f.write(samtools_cmd3+'\n')
        os.system(samtools_cmd3)


        # count total reads
        samtools_reads_cmd = f"samtools view -c {metric_dir}/{base}.sorted.bam"
        reads = os.popen(samtools_reads_cmd).read()
        # count number of mapped reads
        samtools_mapped_cmd = f"samtools view -c -F 260 {metric_dir}/{base}.sorted.bam"
        mapped = os.popen(samtools_mapped_cmd).read()
        # get average coverage for N and G genes Work in Progress
        samtools_Ndepth_cmd = f"samtools depth -b {data_path}/AmericasN.bed {metric_dir}/{base}.sorted.bam " + '| awk \'{ total += $3; count++ } END { print total/1350 }\''
        ncov = os.popen(samtools_Ndepth_cmd).read()
        samtools_Gdepth_cmd = f"samtools depth -b {data_path}/AmericasG.bed {metric_dir}/{base}.sorted.bam " + '| awk \'{ total += $3; count++ } END { print total/1575 }\''
        gcov = os.popen(samtools_Gdepth_cmd).read()
        # Get average length of mapped reads
        samtools_mappeds_cmd = f"samtools view -F 260 -b {metric_dir}/{base}.sorted.bam > {metric_dir}/{base}.mapped.bam"
        samtools_avelen_cmd = f"samtools view {metric_dir}/{base}.mapped.bam " + '| awk \'{ sum_length += length($10) } END { print sum_length/NR }\''
        os.system(samtools_mappeds_cmd)
        readlength = os.popen(samtools_avelen_cmd).read()

        df = df.append(pd.Series([base,reads.rstrip(),mapped.rstrip(),ncov.rstrip(),gcov.rstrip(),readlength.rstrip()], index=df.columns), ignore_index=True)
        
        print("progress: {}/{}".format(i+1, 2*len(files)))

    #df.to_csv(f"{OutputFolder}/covstat.txt", index = False)
    dff = open(f"{OutputFolder}/covstat.txt", 'w')
    dff.write(tabulate(df, headers="keys", showindex=False, tablefmt="psql", floatfmt=".2f"))
    dff.close()
    multiqc_cmd = f"multiqc {qc_dir} -o {OutputFolder}"
    f.write(multiqc_cmd+'\n')
    os.system(multiqc_cmd)

    f.close()

    f=open(f"{OutputFolder}/assembly.log", 'w+')

    for i in range(0, len(files)):
        filec = files[i]

        base = os.path.splitext(os.path.basename(filec))[0]
        base = os.path.splitext(base)[0]
        
        #print(base)
        canu_cmd = f"canu -correct -p asm useGrid=0 -nanopore-raw {filec} -d {assembly_dir}/canu_{base} genomeSize={args.genomeSize} minReadLength={args.minReadLength} readSamplingCoverage={args.readSamplingCoverage} correctedErrorRate={args.correctedErrorRate}"
        f.write(canu_cmd+'\n')
        cdhit_cmd = f"cd-hit-est -i {assembly_dir}/canu_{base}/asm.correctedReads.fasta.gz -o {assembly_dir}/canu_{base}/asm.correctedReads_80.fasta -c 0.80 -p 1 -d 0 -gap -2"
        f.write(cdhit_cmd+'\n')
        rbio_cmd = f"{canu_helper} {args.readCountMin} {assembly_dir}/canu_{base}/asm.correctedReads_80.fasta {assembly_dir}/canu_{base}/asm.correctedReads_80.fasta.clstr {assembly_dir}/canu_{base}/asm.correctedReads_80_filtered.fasta"
        f.write(rbio_cmd+'\n')
        medaka_cmd = f"medaka_consensus -i {filec} -d {assembly_dir}/canu_{base}/asm.correctedReads_80_filtered.fasta -o {assembly_dir}/medaka_{base} -m {args.model}"
        f.write(medaka_cmd+'\n')
        os.system(canu_cmd)
        os.system(cdhit_cmd)
        os.system(rbio_cmd)
        os.system(medaka_cmd)

        cp_cmd = f"mv {assembly_dir}/medaka_{base}/consensus.fasta {OutputFolder}/{base}_consensus.fasta"
        os.system(cp_cmd)
        f.write(cp_cmd+'\n')

        print("progress: {}/{}".format(len(files)+i+1, 2*len(files)))


    f.close()

    print(f"verbose status is {args.verbose}")
    if not args.verbose:
        os.system(f"rm -rf {qc_dir}")
        os.system(f"rm -rf {metric_dir}")
        os.system(f"rm -rf {assembly_dir}")


if __name__ == "__main__":
    sys.exit(main())
