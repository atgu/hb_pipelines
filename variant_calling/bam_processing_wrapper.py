__author__ = 'Sophie Parsa & Lindo Nkambule'


import hailtop.batch as hb
import hailtop.fs as hfs
from typing import Dict, List, Union, Tuple
import bam_processing as jobs


# The functions here are for processing multiple BAM/CRAM file for multiple samples in parallel

# Convert CRAM to BAM first if old CRAM version (v2.0 or older)
def cram_to_bam_wrapper(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        additional_disk: float = 20.0,
        tmp_dir: str = None
):
    for sample_id, sample_bams in samples_and_bams:
        sample_bam_file_paths = sample_bams.split(',')
        bam_idx = [i for i in range(len(sample_bam_file_paths))]

        for bam, idx in zip(sample_bam_file_paths, bam_idx):
            bam_file = b.read_input(bam)
            bam_prefix = f'{sample_id}_{idx}'
            unmapped_bam_size = jobs.size(bam)

            ref_path = 'gs://gnomaf/genome_reference/hs37d5'
            ref_tmp = b.read_input_group(**{'fasta': f'{ref_path}.fa',
                                            'idx': f'{ref_path}.fa.fai',
                                            'ref_dict': f'{ref_path}.dict'})
            jobs.cram_to_bam(
                b=b,
                input_bam=bam_file,
                output_bam_prefix=bam_prefix,
                fasta_reference=ref_tmp,
                disk_size=unmapped_bam_size + 2.0*unmapped_bam_size + additional_disk,
                tmp_dir=f'{tmp_dir}/cram_to_bam/{sample_id}'
            )

    b.run()


def get_read_groups_wrapper(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        is_cram_old: bool = False,
        tmp_dir: str = None,
):
    for sample_id, sample_bams in samples_and_bams:
        sample_bam_file_paths = sample_bams.split(',')
        bam_idx = [i for i in range(len(sample_bam_file_paths))]

        for bam, idx in zip(sample_bam_file_paths, bam_idx):
            bam_prefix = f'{sample_id}_{idx}'
            unmapped_bam_size = jobs.size(bam)

            if is_cram_old:
                input_bam = b.read_input(f'{tmp_dir}/cram_to_bam/{sample_id}/{bam_prefix}.bam')
            else:
                input_bam = b.read_input(bam)

            if not hfs.exists(f'{tmp_dir}/read_group_ids/{sample_id}/{bam_prefix}.rgs.txt'):
                jobs.get_read_groups(
                    b=b,
                    input_bam=input_bam,
                    output_bam_prefix=bam_prefix,
                    storage=unmapped_bam_size+1,
                    tmp_dir=f'{tmp_dir}/read_group_ids/{sample_id}'
                )

    b.run()


def split_by_num_reads_wrapper(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        n_reads: int = 1_000_000,
        additional_disk: float = 20.0,
        tmp_dir: str = None
):
    for sample_id, sample_bams in samples_and_bams:
        sample_bam_file_paths = sample_bams.split(',')
        bam_idx = [i for i in range(len(sample_bam_file_paths))]

        for bam, idx in zip(sample_bam_file_paths, bam_idx):
            unmapped_bam_size = jobs.size(bam)
            input_bam = b.read_input(bam)

            # unmap BAM to multiple uBAM files split by read group
            jobs.split_by_num_reads(
                b=b,
                input_bam=input_bam,
                output_bam_prefix=sample_id,
                n_reads=n_reads,
                disk_size=unmapped_bam_size + unmapped_bam_size*2.5 + additional_disk,
                tmp_dir=f'{tmp_dir}/bams_plit_by_reads'
            )

    b.run()


def bam_to_ubam_rg_wrapper(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        is_cram_old: bool = False,
        additional_disk: float = 20.0,
        tmp_dir: str = None
):
    for sample_id, sample_bams in samples_and_bams:
        sample_bam_file_paths = sample_bams.split(',')
        bam_idx = [i for i in range(len(sample_bam_file_paths))]

        for bam, idx in zip(sample_bam_file_paths, bam_idx):
            bam_prefix = f'{sample_id}_{idx}'
            unmapped_bam_size = jobs.size(bam)

            if is_cram_old:
                input_bam = b.read_input(f'{tmp_dir}/cram_to_bam/{sample_id}/{bam_prefix}.bam')
            else:
                input_bam = b.read_input(bam)

            bam_rg_ids = hfs.open(f'{tmp_dir}/read_group_ids/{sample_id}/{bam_prefix}.rgs.txt').readlines()
            rgs = [l.strip() for l in bam_rg_ids]   # rgs have a \n at the end, hence the .strip()

            # first check if unmapped BAMs exist
            ubams_exist = [hfs.exists(f'{tmp_dir}/unmapped_bams/{sample_id}/{bam_prefix}/{rg}.unmapped.bam') for rg in rgs]

            if not all(ubams_exist):
                # unmap BAM to multiple uBAM files split by read group
                jobs.bam_to_ubam_rg(
                    b=b,
                    input_bam=input_bam,
                    output_bam_prefix=bam_prefix,
                    rg_ids=rgs,
                    disk_size=unmapped_bam_size + unmapped_bam_size*2.5 + additional_disk,
                    tmp_dir=f'{tmp_dir}/unmapped_bams/{sample_id}'
                )

    b.run()


def bam_to_ubam_num_reads_wrapper(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        additional_disk: float = 20.0,
        tmp_dir: str = None
):
    for sample_id, sample_bams in samples_and_bams:
        split_ls = hfs.ls(f'{tmp_dir}/bams_plit_by_reads/{sample_id}/*.bam')
        split_bam_paths = [i.path for i in split_ls]
        bam_idx = [f'shard_{i}' for i in range(len(split_bam_paths))]

        for bam, idx in zip(split_bam_paths, bam_idx):
            bam_prefix = f'{sample_id}_{idx}'
            unmapped_bam_size = jobs.size(bam)
            input_bam = b.read_input(bam)

            # first check if unmapped BAM exists
            ubam_exist = hfs.exists(f'{tmp_dir}/unmapped_bams/{sample_id}/{bam_prefix}.unmapped.bam')
            if not ubam_exist:
                # unmap BAM to multiple uBAM files split by number of reads
                jobs.bam_to_ubam_num_reads(
                    b=b,
                    input_bam=input_bam,
                    output_bam_prefix=bam_prefix,
                    disk_size=unmapped_bam_size + unmapped_bam_size*2.5 + additional_disk,
                    tmp_dir=f'{tmp_dir}/unmapped_bams/{sample_id}'
                )

    b.run()


def sam_to_fastq_and_bwa_mem_and_mba_wrapper(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        split_by_nr: bool = False,
        bwa_reference_files: hb.ResourceGroup = None,
        bwa_ref_size: Union[float, int] = None,
        bwa_disk_multiplier: float = 2.5,
        additional_disk: float = 20.0,
        tmp_dir: str = None
):
    for sample_id, sample_bams in samples_and_bams:
        sample_bam_file_paths = sample_bams.split(',')
        bam_idx = [i for i in range(len(sample_bam_file_paths))]

        for _, idx in zip(sample_bam_file_paths, bam_idx):
            bam_prefix = f'{sample_id}_{idx}'

            # # NOTE: some files will be split by number of reads instead of read groups, incorporate that below (paths)
            if not split_by_nr:
                bam_rg_ids = hfs.open(f'{tmp_dir}/read_group_ids/{sample_id}/{bam_prefix}.rgs.txt').readlines()
                rgs = [l.strip() for l in bam_rg_ids]   # rgs have a \n at the end, hence the .strip()

                # map uBAM files in parallel
                unmapped_bams = [b.read_input(f'{tmp_dir}/unmapped_bams/{sample_id}/{bam_prefix}/{rg}.unmapped.bam') for
                                 rg in rgs]
                unmapped_bam_sizes = [jobs.size(f'{tmp_dir}/unmapped_bams/{sample_id}/{bam_prefix}/{rg}.unmapped.bam') for rg in rgs]

                bams_exist = [hfs.exists(f'{tmp_dir}/mapped_bams/{sample_id}/{bam_prefix}/{rg}.aligned.unsorted.bam') for
                              rg in rgs]

            else:
                unmapped_ls = hfs.ls(f'{tmp_dir}/unmapped_bams/{sample_id}/*.bam')
                unmapped_bam_paths = [i.path for i in unmapped_ls]
                unmapped_bams = [b.read_input(i) for i in unmapped_bam_paths]
                unmapped_bam_sizes = [jobs.size(i) for i in unmapped_bam_paths]

                # not really rgs here but shards, but we use rgs for consistency
                rgs = [f'shard_{i}' for i in range(len(unmapped_bam_paths))]
                bams_exist = [hfs.exists(f'{tmp_dir}/mapped_bams/{sample_id}/{bam_prefix}/{rg}.aligned.unsorted.bam') for
                              rg in rgs]

            if not all(bams_exist):
                # list of BAM files that do not exist
                not_mapped = [(unmapped_bams[i], rgs[i], unmapped_bam_sizes[i]) for i, val in enumerate(bams_exist) if not val]

                bams_rg_mapped = [
                    jobs.sam_to_fastq_and_bwa_mem_and_mba(
                        b=b,
                        input_bam=bam,
                        output_bam_prefix=bam_prefix,
                        rg=rg,
                        fasta_reference=bwa_reference_files,
                        disk_size=round(unmapped_bam_size+bwa_ref_size+(bwa_disk_multiplier*unmapped_bam_size)+additional_disk),
                        tmp_dir=f'{tmp_dir}/mapped_bams/{sample_id}'
                    ).output_bam
                    for bam, rg, unmapped_bam_size in not_mapped
                ]

    b.run()


def mark_duplicates_wrapper(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        additional_disk: float = 20.0,
        md_disk_multiplier: float = 2.25,
        tmp_dir: str = None
):
    for sample_id, _ in samples_and_bams:
        mapped_ls = hfs.ls(f'{tmp_dir}/mapped_bams/{sample_id}/**/*.bam')
        mapped_bam_paths = [i.path for i in mapped_ls]
        mapped_bam_total_size = sum([jobs.size(i) for i in mapped_bam_paths])
        bams_mapped = [b.read_input(i) for i in mapped_bam_paths]

        # MarkDuplicates and SortSam currently take too long for preemptibles if the input data is too large
        use_preemptible = not mapped_bam_total_size > 110.0

        if not hfs.exists(f'{tmp_dir}/mark_duplicates/{sample_id}.aligned.unsorted.duplicates_marked.bam'):
            jobs.mark_duplicates(
                b=b,
                input_bams=bams_mapped,
                output_bam_prefix=sample_id,
                use_preemptible_worker=use_preemptible,
                disk_size=(md_disk_multiplier * mapped_bam_total_size) + additional_disk,
                tmp_dir=f'{tmp_dir}/mark_duplicates'
            )

    b.run()


# sort, check contamination, BQSR (scattered) and gather reports, ApplyBQSR (scattered) and gather BAM files
def sort_contam_bqsr_gather(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        fasta_reference: hb.ResourceGroup = None,
        haplotype_database_file: hb.ResourceGroup = None,
        contamination_sites: hb.ResourceGroup = None,
        ref_size: Union[float, int] = None,
        additional_disk: float = 20.0,
        sort_sam_disk_multiplier: float = 3.25,
        out_dir: str = None,
        tmp_dir: str = None,
        utils: Dict = None
):
    for sample_id, _ in samples_and_bams:
        bam_md = b.read_input(f'{tmp_dir}/mark_duplicates/{sample_id}.aligned.unsorted.duplicates_marked.bam')
        agg_bam_size = jobs.size(f'{tmp_dir}/mark_duplicates/{sample_id}.aligned.unsorted.duplicates_marked.bam')

        sort_sam_bam = jobs.picard_sort_bam(
            b=b,
            input_bam=bam_md,
            output_bam_prefix=sample_id,
            disk_size=(sort_sam_disk_multiplier * agg_bam_size) + additional_disk
        ).output_bam

        # CrosscheckFingerprints requires a haplotype map file
        if haplotype_database_file:
            jobs.cross_check_fingerprints(b=b,
                                          input_bam=sort_sam_bam,
                                          haplotype_database_file=haplotype_database_file,
                                          output_bam_prefix=sample_id,
                                          disk_size=agg_bam_size + additional_disk
                                          )

        jobs.check_contamination(
            b=b,
            input_bam=sort_sam_bam,
            output_bam_prefix=f'{sample_id}.preBqsr',
            fasta_reference=fasta_reference,
            contamination_sites=contamination_sites,
            disk_size=agg_bam_size + ref_size + additional_disk,
            out_dir=out_dir
        )

        if not hfs.exists(f'{tmp_dir}/bqsr/{sample_id}.aligned.duplicate_marked.sorted.bqsr.bam'):
            seq_grouping = hfs.open('gs://h3africa/variant_calling_resources/hg38_sequence_grouping.txt').readlines()
            seq_grouping_subgroup = [l.split('\t') for l in seq_grouping]
            potential_bqsr_divisor = len(seq_grouping_subgroup) - 10
            bqsr_divisor = potential_bqsr_divisor if potential_bqsr_divisor > 1 else 1

            bqsr_reports = [
                jobs.base_recalibrator(
                    b=b,
                    input_bam=sort_sam_bam,
                    fasta_reference=fasta_reference,
                    output_bam_prefix=sample_id,
                    sequence_group_interval=subgroup,
                    utils=utils,
                    disk_size=agg_bam_size + ref_size + additional_disk
                ).recalibration_report
                for subgroup in seq_grouping_subgroup
            ]

            gathered_bqsr_report = jobs.gather_bqsr_reports(
                b=b,
                input_bqsr_reports=[report for report in bqsr_reports],
                output_bam_prefix=sample_id,
                disk_size=5
            ).output_bqsr_report

            seq_grouping_with_unmapped = hfs.open('gs://h3africa/variant_calling_resources/hg38_sequence_grouping_with_unmapped.txt').readlines()
            seq_grouping_with_unmapped_subgroup = [l.split('\t') for l in seq_grouping_with_unmapped]

            # Apply the recalibration model by interval
            recalibrated_bams = [
                jobs.apply_bqsr(
                    b=b,
                    input_bam=sort_sam_bam,
                    fasta_reference=fasta_reference,
                    output_bam_prefix=sample_id,
                    bsqr_report=gathered_bqsr_report,
                    sequence_group_interval=subgroup,
                    disk_size=agg_bam_size + (agg_bam_size / bqsr_divisor) + ref_size + additional_disk
                ).output_bam
                for subgroup in seq_grouping_with_unmapped_subgroup
            ]

            # Merge the recalibrated BAM files resulting from by-interval recalibration
            jobs.gather_bam_files(
                b=b,
                input_bams=recalibrated_bams,
                output_bam_prefix=sample_id,
                disk_size=(2 * agg_bam_size) + additional_disk,
                tmp_dir=f'{tmp_dir}/bqsr'
            )
            # BQSR bins the qualities which makes a significantly smaller bam. Get binned file size

    b.run()


def more_qc_wrapper(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        fasta_reference: hb.ResourceGroup = None,
        ref_size: Union[float, int] = None,
        additional_disk: float = 20.0,
        out_dir: str = None,
        tmp_dir: str = None
):
    for sample_id, _ in samples_and_bams:
        gathered_bams_path = f'{tmp_dir}/bqsr/{sample_id}.aligned.duplicate_marked.sorted.bqsr'
        gathered_bams = b.read_input_group(**{'bam': f'{gathered_bams_path}.bam',
                                              'bam.bai': f'{gathered_bams_path}.bam.bai'})
        binned_qual_bam_size = jobs.size(f'{gathered_bams_path}.bam')

        # QC the final BAM some more (no such thing as too much QC)
        agg_metrics = jobs.collect_aggregation_metrics(
            b=b,
            input_bam=gathered_bams,
            fasta_reference=fasta_reference,
            output_bam_prefix=sample_id,
            disk_size=binned_qual_bam_size + ref_size + additional_disk
        ).output

        mark_dup_metrics = b.read_input(f'{tmp_dir}/mark_duplicates/{sample_id}.marked_dup_metrics.txt')
        jobs.check_pre_validation(
            b=b,
            output_bam_prefix=sample_id,
            duplication_metrics=mark_dup_metrics,
            chimerism_metrics=agg_metrics['alignment_summary_metrics'],
            tmp_dir=f'{tmp_dir}/check_pre_validation'
        )

        # Convert the final merged recalibrated BAM file to CRAM format
        if not hfs.exists(f'{out_dir}/gatk_vc/crams/{sample_id}.cram'):
            jobs.convert_to_cram(
                b=b,
                input_bam=gathered_bams,
                fasta_reference=fasta_reference,
                output_bam_prefix=sample_id,
                disk_size=(2 * binned_qual_bam_size) + ref_size + additional_disk,
                out_dir=f'{out_dir}/crams'
            )
    b.run()


def validate_samfile_wrapper(
        b: hb.Batch,
        samples_and_bams: List[Tuple[str, str]],
        fasta_reference: hb.ResourceGroup = None,
        ref_size: Union[float, int] = None,
        additional_disk: float = 20.0,
        out_dir: str = None,
        tmp_dir: str = None
):
    for sample_id, _ in samples_and_bams:
        # Validate the CRAM file
        cram_file = b.read_input_group(**{'cram': f'{out_dir}/crams/{sample_id}.cram',
                                          'cram.crai': f'{out_dir}/crams/{sample_id}.cram.crai'})
        cram_size = jobs.size(f'{out_dir}/crams/{sample_id}.cram')

        pre_val_metric = b.read_input(f'{tmp_dir}/check_pre_validation/{sample_id}.is.outlier.txt')

        jobs.validate_samfile(
            b=b,
            input_bam=cram_file,
            fasta_reference=fasta_reference,
            output_bam_prefix=sample_id,
            is_outlier_data=pre_val_metric,
            disk_size=cram_size + ref_size + additional_disk
        )
    b.run()
