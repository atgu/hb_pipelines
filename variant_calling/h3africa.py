__author__ = 'Sophie Parsa & Lindo Nkambule'

import argparse
import funcs as fun
import hailtop.batch as hb
import hailtop.fs as hfs
import pandas as pd

from hailtop.batch.job import Job
from wrappers import mark_duplicates_wrapper, sort_contam_bqsr_gather, more_qc_wrapper, validate_samfile_wrapper
from typing import List, Union, Tuple


# GATK Best Practices
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
# https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf-fc-hg38.wdl

# hs37d5 reference is from:
# https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/
inputs = {
    "ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
    "ref_fasta_index": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
    "ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
    "ref_alt": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt",
    "ref_sa": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa",
    "ref_amb": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb",
    "ref_bwt": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt",
    "ref_ann": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann",
    "ref_pac": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac",
    "wgs_calling_interval_list": "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list",
    "calling_int_list": "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list",
    "eval_int_list": "gs://gcp-public-data--broad-references/hg38/v0/wgs_evaluation_regions.hg38.interval_list",
    "contamination_sites_ud": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.contam.UD",
    "contamination_sites_bed": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.contam.bed",
    "contamination_sites_mu": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.contam.mu",
}


def h3a_bam_to_ubam(
        b: hb.batch.Batch,
        input_bam: hb.resource.InputResourceFile = None,
        output_bam_prefix: str = None,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: str = 'lowmem',
        ncpu: int = 4,
        tmp_dir: str = None
) -> Job:
    j = b.new_job(name=f'BamToUbam: {output_bam_prefix}')

    j.cpu(ncpu)
    j.image(img)
    j.memory(memory)
    j.storage(f'{disk_size}Gi')

    java_mem = ncpu * 1

    # filter reads using PrintReads to avoid ERROR:MATE_NOT_FOUND in some BAM files
    j.command(
        f"""
        cd /io
        mkdir tmp/
        mkdir tmp/reverted/
        
        gatk --java-options -Xmx{java_mem}G PrintReads \
            -I {input_bam} \
            -O tmp.fixed.bam \
            --read-filter ProperlyPairedReadFilter
            
        gatk --java-options -Xmx{java_mem}G RevertSam \
            -I tmp.fixed.bam  \
            -O {j.output_bam} \
            --MAX_DISCARD_FRACTION 0.005 \
            --ATTRIBUTE_TO_CLEAR XT \
            --ATTRIBUTE_TO_CLEAR XN \
            --ATTRIBUTE_TO_CLEAR AS \
            --ATTRIBUTE_TO_CLEAR OC \
            --ATTRIBUTE_TO_CLEAR OP \
            --OUTPUT_BY_READGROUP false \
            --SORT_ORDER queryname \
            --RESTORE_ORIGINAL_QUALITIES true \
            --REMOVE_DUPLICATE_INFORMATION true \
            --REMOVE_ALIGNMENT_INFORMATION true \
            --VALIDATION_STRINGENCY SILENT \
            --TMP_DIR `pwd`/tmp
        """
    )

    j.command(
        f"""
        gatk --java-options "-Xms3000m -Xmx{java_mem}G" ValidateSamFile \
            -I {j.output_bam} \
            -MODE SUMMARY
        """
    )

    b.write_output(j.output_bam, f'{tmp_dir}/{output_bam_prefix}.unmapped.bam')

    return j


# Assume out_dir and tmp_dir are different, and write every intermediate file to tmp_dir and final files to out_dir
# Wrappers
def split_by_num_reads_wrapper(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        n_reads: int = 48000000,
        additional_disk: float = 20.0,
        tmp_dir: str = None
):
    for sample_id, sample_bams in samples_and_bams:
        sample_bam_file_paths = sample_bams.split(',')
        assert len(sample_bam_file_paths) == 1
        bam_size = fun.size(sample_bam_file_paths[0])
        input_bam = b.read_input(sample_bam_file_paths[0])

        # unmap BAM to multiple uBAM files split by number of reads
        fun.split_by_num_reads(
            b=b,
            input_bam=input_bam,
            output_bam_prefix=sample_id,
            n_reads=n_reads,
            disk_size=bam_size*3 + bam_size*4 + additional_disk*3,
            tmp_dir=f'{tmp_dir}/bams_plit_by_reads'
        )

    b.run()


def bam_to_ubam_and_map_wrapper(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        bwa_ref_files: hb.ResourceGroup = None,
        ref_files_size: float = None,
        bwa_disk_multiplier: float = 2.5,
        additional_disk: float = 20.0,
        tmp_dir: str = None
):
    for sample_id, sample_bams in samples_and_bams:
        split_ls = hfs.ls(f'{tmp_dir}/bams_plit_by_reads/{sample_id}/*.bam')
        split_bam_paths = [i.path for i in split_ls]
        bam_idx = [f'shard_{i}' for i in range(len(split_bam_paths))]

        for bam, idx in zip(split_bam_paths, bam_idx):
            bam_prefix = f'{sample_id}_{idx}'
            mapped_bam_size = fun.size(bam)
            input_bam = b.read_input(bam)

            # first check if unmapped BAM exists
            ubam_exists = hfs.exists(f'{tmp_dir}/bam_to_ubam/{sample_id}/{bam_prefix}.unmapped.bam')
            if not ubam_exists:
                # unmap BAM to multiple uBAM files split by number of reads
                ubam = h3a_bam_to_ubam(
                    b=b,
                    input_bam=input_bam,
                    output_bam_prefix=bam_prefix,
                    disk_size=mapped_bam_size + mapped_bam_size*2.5,
                    tmp_dir=f'{tmp_dir}/bam_to_ubam/{sample_id}'
                ).output_bam
            else:
                ubam = b.read_input(f'{tmp_dir}/bam_to_ubam/{sample_id}/{bam_prefix}.unmapped.bam')

            # check if mapped BAM exists
            bam_exists = hfs.exists(f'{tmp_dir}/mapped_bams/{sample_id}/{bam_prefix}.aligned.unsorted.bam')
            if not bam_exists:
                fun.sam_to_fastq_and_bwa_mem_and_mba(
                    b=b,
                    input_ubam=ubam,
                    output_bam_prefix=bam_prefix,
                    bwa_ref_files=bwa_ref_files,
                    disk_size=round(mapped_bam_size + bwa_disk_multiplier*mapped_bam_size + ref_files_size + additional_disk),
                    tmp_dir=f'{tmp_dir}/mapped_bams/{sample_id}'
                )

    b.run()


def process_samples(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        fasta_ref_files: hb.ResourceGroup = None,
        bwa_ref_files: hb.ResourceGroup = None,
        haplotype_database_file: hb.ResourceGroup = None,
        contamination_sites: hb.ResourceGroup = None,
        steps: str = 'split,map,mdups,bqsr',
        ref_gb: Union[float, int] = None,
        bwa_ref_gb: Union[float, int] = None,
        bwa_disk_multiplier: float = 2.5,
        md_disk_multiplier: float = 2.25,
        sort_sam_disk_multiplier: float = 3.25,
        additional_disk: float = 20.0,
        out_dir: str = None,
        tmp_dir: str = None,
):
    steps_list = steps.split(',')
    steps_to_run = [x.lower() for x in steps_list]
    unknown_steps = [i for i in steps_to_run if i not in ['split', 'map', 'mdups', 'bqsr']]

    if len(unknown_steps) > 0:
        raise SystemExit(f'Incorrect step(s) {unknown_steps} selected. Options are [split, map, mdups, bqsr]')

    print(f'--- Steps to run: {steps_to_run} ---')

    if 'split' in steps_to_run:
        split_by_num_reads_wrapper(b=b,
                                   samples_and_bams=samples_and_bams,
                                   additional_disk=additional_disk,
                                   tmp_dir=tmp_dir)

    if 'map' in steps_to_run:
        bam_to_ubam_and_map_wrapper(b=b,
                                    samples_and_bams=samples_and_bams,
                                    bwa_ref_files=bwa_ref_files,
                                    ref_files_size=ref_gb+bwa_ref_gb,
                                    bwa_disk_multiplier=bwa_disk_multiplier,
                                    additional_disk=additional_disk,
                                    tmp_dir=tmp_dir)

    if 'mdups' in steps_to_run:
        mark_duplicates_wrapper(b=b,
                                samples_and_bams=samples_and_bams,
                                additional_disk=additional_disk,
                                md_disk_multiplier=md_disk_multiplier,
                                tmp_dir=tmp_dir)

    if 'bqsr' in steps_to_run:
        sort_contam_bqsr_gather(b=b,
                                samples_and_bams=samples_and_bams,
                                fasta_reference=fasta_ref_files,
                                haplotype_database_file=haplotype_database_file,
                                contamination_sites=contamination_sites,
                                ref_size=ref_gb,
                                additional_disk=additional_disk,
                                sort_sam_disk_multiplier=sort_sam_disk_multiplier,
                                out_dir=out_dir,
                                tmp_dir=tmp_dir)

    more_qc_wrapper(b=b,
                    samples_and_bams=samples_and_bams,
                    fasta_reference=fasta_ref_files,
                    ref_size=ref_gb,
                    additional_disk=additional_disk,
                    out_dir=out_dir,
                    tmp_dir=tmp_dir)

    validate_samfile_wrapper(b=b,
                             samples_and_bams=samples_and_bams,
                             fasta_reference=fasta_ref_files,
                             ref_size=ref_gb,
                             additional_disk=additional_disk,
                             out_dir=out_dir,
                             tmp_dir=tmp_dir)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-files', type=str, required=True)
    parser.add_argument('--out-dir', type=str, required=True)
    parser.add_argument('--steps', type=str, default='split,map,mdups,bqsr')
    parser.add_argument('--tmp-dir', type=str, required=True)
    parser.add_argument('--billing-project', type=str, required=True)

    args = parser.parse_args()

    backend = hb.ServiceBackend(
        billing_project=args.billing_project,
        regions=['us-central1'],
        remote_tmpdir=args.tmp_dir,
    )

    batch = hb.Batch(
        'H3Africa-BAM-Processing',
        backend=backend,
    )

    ref_fasta = batch.read_input_group(**{'ref_fasta': inputs['ref_fasta'],
                                          'ref_fasta_index': inputs['ref_fasta_index'],
                                          'ref_dict': inputs['ref_dict']})

    # BWA requires additional reference files
    ref_fasta_bwa = batch.read_input_group(**{'ref_fasta': inputs['ref_fasta'],
                                              'ref_fasta_index': inputs['ref_fasta_index'],
                                              'ref_dict': inputs['ref_dict'],
                                              'ref_alt': inputs['ref_alt'],
                                              'ref_sa': inputs['ref_sa'],
                                              'ref_amb': inputs['ref_amb'],
                                              'ref_bwt': inputs['ref_bwt'],
                                              'ref_ann': inputs['ref_ann'],
                                              'ref_pac': inputs['ref_pac']})

    contamination_sites = batch.read_input_group(**{'ud': inputs['contamination_sites_ud'],
                                                    'mu': inputs['contamination_sites_mu'],
                                                    'bed': inputs['contamination_sites_bed']})

    # Get the size of the standard reference files as well as the additional reference files needed for BWA
    ref_size = fun.size(inputs['ref_fasta']) + fun.size(inputs['ref_fasta_index']) + fun.size(inputs['ref_dict'])
    bwa_ref_size = ref_size + fun.size(inputs['ref_alt']) + fun.size(inputs['ref_amb']) + fun.size(inputs['ref_ann']) + \
                   fun.size(inputs['ref_bwt']) + fun.size(inputs['ref_pac']) + fun.size(inputs['ref_sa'])

    bams = pd.read_csv(args.input_files, sep='\t', header=None, names=['id', 'files'])
    files = [(sample, paths) for sample, paths in zip(bams['id'], bams['files'])]

    process_samples(
        b=batch,
        samples_and_bams=files,
        fasta_ref_files=ref_fasta,
        bwa_ref_files=ref_fasta_bwa,
        contamination_sites=contamination_sites,
        steps=args.steps,
        ref_gb=ref_size,
        bwa_ref_gb=bwa_ref_size,
        out_dir=args.out_dir,
        tmp_dir=args.tmp_dir
    )


if __name__ == '__main__':
    main()
