__author__ = 'Lindo Nkambule'

import argparse
from typing import List, Union
import hailtop.batch as hb
import hailtop.fs as hfs
import pandas as pd
from hailtop.batch.job import Job


utils = {"ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
         "ref_fasta_index": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
         "ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
         "calling_int_list": "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list",
         "eval_int_list": "gs://gcp-public-data--broad-references/hg38/v0/wgs_evaluation_regions.hg38.interval_list",
         "dbsnp_resource_vcf": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
         "hc_contamination": 0,
         "haplotype_scatter_count": 150,
         "break_bands_at_multiples_of": 1000000
         }


# Perform variant calling on the sub-intervals, and then gather the results
def split_interval_list(
        b: hb.Batch = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
) -> Job:
    """
    Break the calling interval_list into sub-intervals
    :param b: Batch object to add jobs to
    :param img: image to use for the job
    :return: a Job object with a single output j.intervals of type ResourceGroup
    """
    j = b.new_job(name='scatter-interval-list')
    j.image(img)
    j.memory('standard')
    j.storage('1G')

    j.command(
        f"""set -e
        mkdir /scatter_intervals
        gatk --java-options "-Xms1g" IntervalListTools \
        --SCATTER_COUNT {utils['haplotype_scatter_count']} \
        --SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
        --UNIQUE true \
        --SORT true \
        --BREAK_BANDS_AT_MULTIPLES_OF {utils['break_bands_at_multiples_of']} \
        -I {utils['calling_int_list']} \
        -O /scatter_intervals
      """
    )

    # e.g if scatter_count = 50, /scatter_intervals will have temp_0001_of_50, temp_0002_of_50, ..., temp_0050_of_50
    # file names under /scatter_intervals/temp_00_of_ are not numbered, it's just scattered.interval_list
    # use .interval_list for Picard-style interval OR .list or .intervals for GATK-style
    # https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists
    for idx in range(utils['haplotype_scatter_count']):
        j.command(f"cp /scatter_intervals/temp_{str(idx + 1).zfill(4)}_of_{utils['haplotype_scatter_count']}/scattered.interval_list {j[f'interval_{idx}.interval_list']}")

    return j


def haplotype_caller(
        b: hb.batch.Batch,
        input_bam: Union[str, hb.ResourceGroup],
        interval_list_file: hb.ResourceFile,
        output_gvcf_name: str = None,
        # img: str = 'us.gcr.io/broad-gatk/gatk:4.2.6.1',
        img: str = 'docker.io/broadinstitute/gatk:latest',
        ncpu: int = 2,
        storage: int = 2,
        tmp_dir: str = None
) -> Job:
    """
    Call germline SNPs and indels
    """

    j = b.new_job(name=f'HaplotypeCaller: {output_gvcf_name}')

    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi'
        }
    )

    java_mem = ncpu * 4    # ‘lowmem’ ~1Gi/core, ‘standard’ ~4Gi/core, and ‘highmem’ ~7Gi/core in Hail Batch

    j.image(img)
    j.cpu(ncpu)
    j.memory('standard')
    j.storage(f'{storage}Gi')
    j.command(
        f"""gatk --java-options "-Xmx{java_mem}G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
            -R {utils['ref_fasta']} \
            -I {input_bam} \
            -L {interval_list_file} \
            -O {j.output_gvcf['g.vcf.gz']} \
            -contamination {utils['hc_contamination']} \
            -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
            -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            -ERC GVCF \
            --gcs-max-retries 100 \
            --create-output-variant-index"""
    )

    if tmp_dir:
        b.write_output(j.output_gvcf, f'{tmp_dir}/scatter_{utils["haplotype_scatter_count"]}/{output_gvcf_name}')

    return j


def merge_gvcfs(
        b: hb.batch.Batch,
        gvcf_list: List[hb.ResourceGroup] = None,
        output_gvcf_name: str = None,
        img: str = 'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330',
        memory: int = 8,
        ncpu: int = 4,
        out_dir: str = None,
        storage: int = 100
) -> Job:
    """
    Merge GVCFs from scattered HaplotypeCaller runs
    """

    input_cmdl = ' '.join([f'I={v["g.vcf.gz"]}' for v in gvcf_list])
    java_mem = ncpu * 4    # ‘lowmem’ ~1Gi/core, ‘standard’ ~4Gi/core, and ‘highmem’ ~7Gi/core in Hail Batch

    j = b.new_job(name=f'Merge gVCFs: {output_gvcf_name}')

    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi'
        }
    )

    j.image(img)
    j.memory(f'{memory}Gi')
    j.cpu(ncpu)
    j.storage(f'{storage}Gi')

    # merge
    j.command(
        f"""java -Xmx{java_mem}G -jar /usr/gitc/picard.jar \
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
        memory: int = 8,
        storage: int = 10,
        job_name: str = None
) -> Job:
    """
    Validate the gVCF output of HaplotypeCaller
    """

    j = b.new_job(name=f'Validate gVCF: {job_name}')
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
        steps: str = 'call,merge',
        out_dir: str = None,
        tmp_dir: str = None,
        billing_project: str = None,
):
    """
    Call variants from CRAM/BAM for multiple samples in parallel
    """
    steps_list = steps.split(',')
    steps_to_run = [x.lower() for x in steps_list]
    unknown_steps = [i for i in steps_to_run if i not in ['call', 'merge']]

    if len(unknown_steps) > 0:
        raise SystemExit(f'Incorrect process(es) {unknown_steps} selected. Options are [call, merge]')

    print(f'--- Steps to run: {steps_to_run} ---')

    backend = hb.ServiceBackend(
        billing_project=billing_project,
        remote_tmpdir=tmp_dir,
    )

    b = hb.Batch(
        'GATK-HC-Pipeline',
        backend=backend,
    )

    # get filepaths into a list
    bams = pd.read_csv(input_files, sep='\t', header=None, names=['id', 'files'])
    files = [(sample, bam_path) for sample, bam_path in zip(bams['id'], bams['files'])]

    # 1. split intervals
    if 'call' in steps_to_run:
        intervals = split_interval_list(
            b=b
        )

    for sample_id, bam in files:
        # 2. HC scattered mode
        if 'call' in steps_to_run:
            scattered_gvcfs = [
                haplotype_caller(
                    b=b,
                    input_bam=bam,
                    interval_list_file=intervals[f'interval_{idx}.interval_list'],
                    output_gvcf_name=f'{sample_id}_{idx}',
                    tmp_dir=tmp_dir
                ).output_gvcf
                for idx in range(utils['haplotype_scatter_count'])
            ]
        else:
            p = f'{tmp_dir}/scatter_{utils["haplotype_scatter_count"]}'
            scattered_gvcfs = [b.read_input_group(**{'g.vcf.gz': f'{p}/{sample_id}_{idx}',
                                                     'g.vcf.gz.tbi': f'{p}/{sample_id}_{idx}.bai'})
                               for idx in range(utils['haplotype_scatter_count'])]

        # 3. Gather gvcfs and index
        if 'merge' in steps_to_run:
            if not hfs.exists(f'{out_dir}/gvcfs/{sample_id}.g.vcf.gz'):
                merged_gvcf = merge_gvcfs(
                    b=b,
                    gvcf_list=scattered_gvcfs,
                    output_gvcf_name=sample_id,
                    out_dir=f'{out_dir}/gvcfs/'
                ).output_gvcf
            else:
                merged_gvcf = b.read_input_group(**{'g.vcf.gz': f'{out_dir}/gvcfs/{sample_id}.g.vcf.gz',
                                                    'g.vcf.gz.tbi': f'{out_dir}/gvcfs/{sample_id}.g.vcf.gz.tbi'})

            # 4. Index and validate GVCF
            validate_gvcf(
                b=b,
                input_gvcf=merged_gvcf,
                job_name=sample_id
                )

    b.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-files', type=str, required=True)
    parser.add_argument('--out-dir', type=str, required=True)
    parser.add_argument('--tmp-dir', type=str, required=True)
    parser.add_argument('--steps', type=str, default='call,merge')
    parser.add_argument('--billing-project', type=str, required=True)

    args = parser.parse_args()

    run_gatk_hc(input_files=args.input_files, out_dir=args.out_dir, tmp_dir=args.tmp_dir,
                billing_project=args.billing_project)


if __name__ == '__main__':
    main()
