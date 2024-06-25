__author__ = 'Sophie Parsa & Lindo Nkambule'

import argparse
import funcs as fun
import hailtop.batch as hb
import hailtop.fs as hfs
import os
import pandas as pd

from wrappers import sort_contam_bqsr_gather, more_qc_wrapper, validate_samfile_wrapper
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
    "hs37d5_fasta": "gs://gnomaf/genome_reference/hs37d5.fa",
    "hs37d5_index": "gs://gnomaf/genome_reference/hs37d5.fa.fai",
    "hs37d5_dict": "gs://gnomaf/genome_reference/hs37d5.dict",
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


# Assume out_dir and tmp_dir are different, and write every intermediate file to tmp_dir and final files to out_dir


# Wrappers
# unamps old CRAM/BAM and maps uBAM to current reference genome
def cram_to_bam_wrapper(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        old_ref_files: hb.ResourceGroup = None,
        bwa_ref_files: hb.ResourceGroup = None,
        ref_files_size: float = None,
        bwa_disk_multiplier: float = 2.5,
        additional_disk: float = 20.0,
        tmp_dir: str = None
):
    for sample_id, sample_bams in samples_and_bams:
        sample_bam_file_paths = sample_bams.split(',')
        bam_idx = [i for i in range(len(sample_bam_file_paths))]

        for bam, idx in zip(sample_bam_file_paths, bam_idx):
            extension = os.path.splitext(os.path.basename(bam))

            bam_prefix = f'{sample_id}_{idx}'
            unmapped_bam_size = fun.size(bam)

            # check if uBAM exists
            ubam_exists = hfs.exists(f'{tmp_dir}/cram_to_ubam/{sample_id}/{bam_prefix}.unmapped.bam')
            if not ubam_exists:
                if extension == '.cram':
                    infile = b.read_input_group(**{'cram': bam,
                                                   'md5': f'{bam}.md5'})
                    ubam = fun.revert_cram_to_ubam(
                        b=b,
                        input_cram=infile,
                        output_bam_prefix=bam_prefix,
                        old_ref_files=old_ref_files,
                        disk_size=round(unmapped_bam_size + 2.0*unmapped_bam_size + ref_files_size + additional_disk),
                        tmp_dir=f'{tmp_dir}/cram_to_ubam/{sample_id}'
                    ).output_bam
                else:
                    infile = b.read_input(bam)
                    ubam = fun.revert_bam_to_ubam(
                        b=b,
                        input_bam=infile,
                        output_bam_prefix=bam_prefix,
                        disk_size=round(unmapped_bam_size + 2.0*unmapped_bam_size + additional_disk),
                        tmp_dir=f'{tmp_dir}/cram_to_ubam/{sample_id}'
                    ).output_bam
            else:
                ubam = b.read_input(f'{tmp_dir}/cram_to_ubam/{sample_id}/{bam_prefix}.unmapped.bam')

            # check if mapped BAM exists
            bam_exists = hfs.exists(f'{tmp_dir}/mapped_bams/{sample_id}/{bam_prefix}.aligned.unsorted.bam')
            if not bam_exists:
                fun.sam_to_fastq_and_bwa_mem_and_mba(
                    b=b,
                    input_ubam=ubam,
                    output_bam_prefix=bam_prefix,
                    bwa_ref_files=bwa_ref_files,
                    disk_size=round(unmapped_bam_size + bwa_disk_multiplier*unmapped_bam_size + ref_files_size + additional_disk),
                    tmp_dir=f'{tmp_dir}/mapped_bams/{sample_id}'
                )

    b.run()


def mark_duplicates_wrapper(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        additional_disk: float = 20.0,
        md_disk_multiplier: float = 2.25,
        tmp_dir: str = None
):
    for sample_id, _ in samples_and_bams:
        mapped_ls = hfs.ls(f'{tmp_dir}/mapped_bams/{sample_id}/*.bam')
        mapped_bam_paths = [i.path for i in mapped_ls]
        mapped_bam_total_size = sum([fun.size(i) for i in mapped_bam_paths])
        bams_mapped = [b.read_input(i) for i in mapped_bam_paths]

        # MarkDuplicates and SortSam currently take too long for preemptibles if the input data is too large
        use_preemptible = not mapped_bam_total_size > 110.0

        if not hfs.exists(f'{tmp_dir}/mark_duplicates/{sample_id}.aligned.unsorted.duplicates_marked.bam'):
            fun.mark_duplicates(
                b=b,
                input_bams=bams_mapped,
                output_bam_prefix=sample_id,
                use_preemptible_worker=use_preemptible,
                disk_size=(md_disk_multiplier * mapped_bam_total_size) + additional_disk,
                tmp_dir=f'{tmp_dir}/mark_duplicates'
            )

    b.run()


def process_samples(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        fasta_ref_files: hb.ResourceGroup = None,
        old_ref_files: hb.ResourceGroup = None,
        bwa_ref_files: hb.ResourceGroup = None,
        haplotype_database_file: hb.ResourceGroup = None,
        contamination_sites: hb.ResourceGroup = None,
        steps: str = 'map,mdups,bqsr',
        ref_gb: Union[float, int] = None,
        bwa_ref_gb: Union[float, int] = None,
        old_ref_gb: Union[float, int] = None,
        bwa_disk_multiplier: float = 2.5,
        md_disk_multiplier: float = 2.25,
        sort_sam_disk_multiplier: float = 3.25,
        additional_disk: float = 20.0,
        out_dir: str = None,
        tmp_dir: str = None,
):
    steps_list = steps.split(',')
    steps_to_run = [x.lower() for x in steps_list]
    unknown_steps = [i for i in steps_to_run if i not in ['map', 'mdups', 'bqsr']]

    if len(unknown_steps) > 0:
        raise SystemExit(f'Incorrect process(es) {unknown_steps} selected. Options are [map, mdups, bqsr]')

    print(f'--- Steps to run: {steps_to_run} ---')

    if 'map' in steps_to_run:
        cram_to_bam_wrapper(b=b,
                            samples_and_bams=samples_and_bams,
                            old_ref_files=old_ref_files,
                            bwa_ref_files=bwa_ref_files,
                            ref_files_size=ref_gb+bwa_ref_gb+old_ref_gb,
                            bwa_disk_multiplier=bwa_disk_multiplier,
                            additional_disk=20.0,
                            tmp_dir=tmp_dir
                            )

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
    parser.add_argument('--steps', type=str, default='map,mdups,bqsr')
    parser.add_argument('--tmp-dir', type=str, required=True)
    parser.add_argument('--billing-project', type=str, required=True)

    args = parser.parse_args()

    backend = hb.ServiceBackend(
        billing_project=args.billing_project,
        regions=['us-central1'],
        remote_tmpdir=args.tmp_dir,
    )

    batch = hb.Batch(
        'Uganda-BAM-Processing',
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

    hs37d5_ref = batch.read_input_group(**{'ref_fasta': inputs['hs37d5_fasta'],
                                           'ref_fasta_index': inputs['hs37d5_index'],
                                           'ref_dict': inputs['hs37d5_dict']})
    hs37d5_ref_size = fun.size(inputs['hs37d5_fasta']) + fun.size(inputs['hs37d5_index']) + fun.size(inputs['hs37d5_dict'])

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
        old_ref_files=hs37d5_ref,
        bwa_ref_files=ref_fasta_bwa,
        contamination_sites=contamination_sites,
        steps=args.steps,
        ref_gb=ref_size,
        bwa_ref_gb=bwa_ref_size,
        old_ref_gb=hs37d5_ref_size,
        out_dir=args.out_dir,
        tmp_dir=args.tmp_dir
    )


if __name__ == '__main__':
    main()
