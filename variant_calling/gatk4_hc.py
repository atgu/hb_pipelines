#!/usr/bin/env python3

__author__ = 'Lindo Nkambule'

import argparse
import os
from typing import List, Union
import hailtop.batch as hb
import hail as hl
import pandas as pd
from hailtop.batch.job import Job


utils = {"ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
         "ref_fasta_index": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
         "ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
         "calling_int_list": "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list",
         "eval_int_list": "gs://gcp-public-data--broad-references/hg38/v0/wgs_evaluation_regions.hg38.interval_list",
         "dbsnp_resource_vcf": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
         "hc_contamination": 0
         }


def get_file_size(file):
    """Get file size"""

    file_info = hl.utils.hadoop_stat(file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)

    return size_gigs


def scatter_interval_list(
        b: hb.batch.Batch,
        scatter_count: int = 50,
        break_bands_at_multiples_of: int = 1000000,
        img: str = 'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330'
) -> Job:
    """
    break the calling interval list into sub-intervals
    :param b: batch
    :param scatter_count: the number of files into which to scatter the resulting list by locus
    :param break_bands_at_multiples_of: if set to a positive value will create a new interval list with the original
    intervals broken up at integer multiples of this value. Set to 0 to NOT break up intervals
    :param img: image to use for the job
    :return:
    """

    calling_interval_list = b.read_input(utils['calling_int_list'])

    j = b.new_job(name='scatter-interval-list')

    j.image(img)
    j.command('mkdir /scatter_intervals')
    j.command(
        f"""java -Xms1g -jar /usr/gitc/picard.jar \
            IntervalListTools \
            SCATTER_COUNT={scatter_count} \
            SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
            UNIQUE=true \
            SORT=true \
            BREAK_BANDS_AT_MULTIPLES_OF={break_bands_at_multiples_of} \
            INPUT={calling_interval_list} \
            OUTPUT=/scatter_intervals"""
    )

    # file names under /scatter_intervals/temp_00_of_ are not numbered, it's just scattered.interval_list
    # use .interval_list for Picard-style interval OR .list or .intervals for GATK-style
    # https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists
    for idx in range(scatter_count):
        j.command(f"cp /scatter_intervals/temp_{str(idx + 1).zfill(4)}_of_{scatter_count}/scattered.interval_list {j[f'interval_{idx}.interval_list']}")

    return j


def haplotype_caller(
        b: hb.batch.Batch,
        input_bam: Union[str, hb.ResourceGroup],
        interval_list_file: hb.ResourceFile,
        job_name: str = None,
        img: str = 'us.gcr.io/broad-gatk/gatk:4.2.6.1',
        memory: float = 1,
        ncpu: int = 2,
        storage: int = 1
) -> Job:
    """
    Call germline SNPs and indels
    :param b: batch
    :param input_bam: input BAM file
    :param interval_list_file: interval list file with intervals to run variant calling on
    :param bam_filename_no_ext: BAM filename without extension
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :param storage: storage to use fo the job
    :return:
    """

    j = b.new_job(name=f'hc: {job_name}')

    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi'
        }
    )

    j.image(img)
    j.cpu(ncpu)
    j.memory('lowmem')
    j.storage(f'{storage}Gi')
    j.command(
        f"""gatk --java-options "-Xmx{memory}G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
            -R {utils['ref_fasta']} \
            -I {input_bam} \
            -L {interval_list_file} \
            -O {j.output_gvcf['g.vcf.gz']} \
            -contamination {utils['hc_contamination']} \
            -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
            -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            -ERC GVCF \
            --create-output-variant-index"""
    )

    return j


def merge_gvcfs(
        b: hb.batch.Batch,
        gvcf_list: List[hb.ResourceGroup] = None,
        output_gvcf_name: str = None,
        img: str = 'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330',
        memory: int = 5,
        out_dir: str = None,
        storage: int = 20
) -> Job:
    """
    Merge GVCFs from scattered HaplotypeCaller runs
    :param b: batch
    :param gvcf_list: list of GVCF files to merge
    :param output_gvcf_name: output GVCF name
    :param img: image to use for the job
    :param storage: Storage to use fo the job
    :param out_dir: output directory
    :param memory: job memory
    :return:
    """

    input_cmdl = ' '.join([f'I={v["g.vcf.gz"]}' for v in gvcf_list])

    j = b.new_job(name=f'merge: {output_gvcf_name}')

    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi'
        }
    )

    j.image(img)
    j.memory(f'{memory}Gi')
    j.storage(f'{storage}Gi')

    # merge
    j.command(
        f"""java -Xms3000m -jar /usr/gitc/picard.jar \
            MergeVcfs \
            {input_cmdl} \
            O={j.output_gvcf['g.vcf.gz']}""")

    # index
    j.command(
        f"""java -Xms2000m -jar /usr/gitc/gatk4/gatk-package-4.beta.5-local.jar \
            IndexFeatureFile \
            -F {j.output_gvcf['g.vcf.gz']} \
            -O {j.output_gvcf['g.vcf.gz.tbi']}""")

    b.write_output(j.output_gvcf, f'{out_dir}/{output_gvcf_name}')

    return j


def validate_gvcf(
        b: hb.batch.Batch,
        input_gvcf: hb.ResourceGroup,
        img: str = 'us.gcr.io/broad-gatk/gatk:4.2.0.0',
        memory: int = 7,
        storage: int = 10,
        job_name: str = None
) -> Job:
    """
    Validate the GVCF output of HaplotypeCaller
    :param b: batch
    :param input_gvcf: GVCF file to validate
    :param img: image to use for the job
    :param memory: job memory
    :param storage: storage to use fo the job
    :param output_gvcf_name: output GVCF index name
    :return:
    """

    j = b.new_job(name=f'validate: {job_name}')
    j.image(img)
    j.memory(f'{memory}Gi')
    j.storage(f'{storage}Gi')

    calling_interval_list = b.read_input(utils['calling_int_list'])

    # validate
    j.command(
        f"""gatk --java-options -Xms6000m \
            ValidateVariants \
            -V {input_gvcf['g.vcf.gz']} \
            -R {utils['ref_fasta']} \
            -L {calling_interval_list} \
            -gvcf \
            --validation-type-to-exclude ALLELES \
            --dbsnp {utils['dbsnp_resource_vcf']}""")

    return j


def run_gatk_hc(
        input_files: str = None,
        scatter_count: int = 150,
        out_bucket: str = None,
        billing_project: str = None,
):
    """
    Add jobs that perform the allele-specific VQSR variant QC
    :param input_files: text file containing paths, one per line, to files to be processed
    :param scatter_count: the number of files into which to scatter the resulting list by locus
    :param out_bucket: bucket to write output files to
    :param billing_project: billing project to use for the batch
    :return:
    """
    tmp_hc_bucket = f'{out_bucket}/gatk4_hc/'
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        remote_tmpdir=tmp_hc_bucket,
    )

    b = hb.Batch(
        'GATK-HC-Pipeline',
        backend=backend,
    )

    # get filepaths into a list
    files_df = pd.read_csv(input_files, header=None, names=['file'])
    files = files_df.file.values.tolist()

    # 1. split intervals
    intervals = scatter_interval_list(
        b=b,
        scatter_count=scatter_count
    )

    for bam in files:
        file_name = os.path.basename(bam)
        file_no_ext = os.path.splitext(file_name)[0]

        # 2. HC scattered mode
        scattered_gvcfs = [
            haplotype_caller(
                b=b,
                input_bam=bam,
                interval_list_file=intervals[f'interval_{idx}.interval_list'],
                job_name=f'{idx}_{file_no_ext}'
            ).output_gvcf
            for idx in range(scatter_count)
        ]

        # 3. Gather gvcfs and index
        merge_out_dir = out_bucket + '/gvcfs'
        merged_gvcf_j = merge_gvcfs(
            b=b,
            gvcf_list=scattered_gvcfs,
            output_gvcf_name=file_no_ext,
            out_dir=merge_out_dir
        ).output_gvcf

        # 4. Index and validate GVCF
        validate_j = validate_gvcf(
            b=b,
            input_gvcf=merged_gvcf_j,
            job_name=file_no_ext
            )

    b.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-files', type=str, required=True)
    parser.add_argument('--out-bucket', type=str, required=True)
    parser.add_argument('--billing-project', type=str, required=True)

    args = parser.parse_args()

    run_gatk_hc(input_files=args.input_files, out_bucket=args.out_bucket, billing_project=args.billing_project)


if __name__ == '__main__':
    main()
