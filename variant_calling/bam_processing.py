__author__ = 'Lindo Nkambule'


from typing import List, Union
import hailtop.batch as hb
import hail as hl
from hailtop.batch.job import Job


# https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf-fc-hg38.wdl
utils = {"ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
         "ref_fasta_index": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
         "ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
         "calling_int_list": "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list",
         "eval_int_list": "gs://gcp-public-data--broad-references/hg38/v0/wgs_evaluation_regions.hg38.interval_list",
         "dbsnp_resource_vcf": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
         "gatk_img": "us.gcr.io/broad-gatk/gatk:4.2.0.0",
         "hc_contamination": 0}


def get_file_size(file):
    """Get file size"""

    file_info = hl.utils.hadoop_stat(file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)

    return size_gigs


def bam_to_cram(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup = None,
        fasta_reference: hb.ResourceGroup = None,
        cram_out_name: str = None,
        memory: int = 8,
        ncpu: int = 8,
        storage: int = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        out_dir: str = None
) -> Job:
    """
    Convert BAM to CRAM
    :param b: batch
    :param input_bam: BAM file to be converted to CRAM
    :param fasta_reference: reference files
    :param cram_out_name: output CRAM filename
    :param memory: job memory
    :param ncpu: number of CPUs
    :param storage: storage to use fo the job
    :param img: image to use for the job
    :param out_dir: output directory

    :return:
    """
    # Lindo reminder: Do the storage estimation outside of this function
    # file sizes for storage
    # bam_size = get_file_size(input_bam)
    # output_cram_size = bam_size * 0.60
    # ref_size = get_file_size(utils['ref_fasta'])
    # disk_size = round(bam_size + output_cram_size + ref_size) + 15

    j = b.new_job(name=f'BAM-TO-CRAM: {cram_out_name}')

    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        }
    )

    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}Gi')
    j.storage(f'{storage}Gi')
    j.command(
        f"""
        samtools view -@{ncpu-1} -T {fasta_reference['fasta']} -O cram -o {j.output_cram['cram']} {input_bam['bam']}
        samtools index {j.output_cram['cram']} {j.output_cram['cram.crai']}
        """
    )

    if out_dir:
        b.write_output(j.output_cram, f'{out_dir}/{cram_out_name}')

    return j


def mark_duplicates(
        b: hb.batch.Batch,
        input_cram: hb.ResourceGroup,
        fasta_reference: hb.ResourceGroup,
        bam_filename_no_ext: str = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: float = 8,
        ncpu: int = 2,
        storage: int = 200,
        output_path: str = None,
) -> Job:
    """
    Mark duplicates and convert BAM to CRAM
    :param b: batch
    :param input_cram: input CRAM file
    :param fasta_reference: reference files
    :param bam_filename_no_ext: BAM filename without extension
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :param storage: storage to use fo the job
    :param output_path: path to write output files to
    :return: Job object
    """

    j = b.new_job(name=f'MarkDups: {bam_filename_no_ext}')

    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        }
    )

    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}Gi')
    j.storage(f'{storage}Gi')
    j.command(
        f"""
        cd /io
        mkdir tmp/
        gatk --java-options -Xms4000m MarkDuplicates \
            -I {input_cram['cram']} \
            -O {j.output_cram['cram']} \
            -M {j.markdup_metrics} \
            -R {fasta_reference['fasta']} \
            --VALIDATION_STRINGENCY SILENT \
            --CLEAR_DT false \
            --ADD_PG_TAG_TO_READS false \
            --TMP_DIR `pwd`/tmp
            
        samtools index {j.output_cram['cram']} {j.output_cram['cram.crai']}
        """
    )

    if output_path:
        # only write out the metrics file
        b.write_output(j.markdup_metrics,
                       f'{output_path}/intermediate/metrics/MarkDuplicates/{bam_filename_no_ext}.txt')

    return j


def picard_sort_cram(
        b: hb.batch.Batch,
        input_cram: Union[str, hb.ResourceGroup],
        fasta_reference: hb.ResourceGroup,
        bam_filename_no_ext: str = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: float = 8,
        ncpu: int = 2,
        storage: int = 100,
) -> Job:
    """
    Use picard to sort CRAM file (slower than samtools)
    :param b: batch
    :param input_cram: input CRAM file to be sorted
    :param fasta_reference: reference files
    :param bam_filename_no_ext: BAM filename without extension
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :param storage: storage to use fo the job
    :return: Job object
    """

    j = b.new_job(name=f'Picard-Sort-CRAM: {bam_filename_no_ext}')

    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        }
    )

    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}Gi')
    j.storage(f'{storage}Gi')
    j.command(
        f"""cd /io
        mkdir tmp/
        gatk --java-options -Xms4000m SortSam \
            -I {input_cram['cram']} \
            -O {j.output_cram['cram']} \
            -R {fasta_reference['fasta']} \
            -SO coordinate \
            --CREATE_INDEX true \
            --MAX_RECORDS_IN_RAM 300000 \
            --TMP_DIR `pwd`/tmp
        """
    )

    return j


def samtools_sort_cram(
        b: hb.batch.Batch,
        input_cram: hb.ResourceGroup,
        bam_filename_no_ext: str = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: float = 8,
        ncpu: int = 8,
        storage: int = 100,
) -> Job:
    """
    Use Samtools to sort CRAM file (faster than picard)
    :param b: batch
    :param input_cram: input CRAM file to be sorted
    :param bam_filename_no_ext: BAM filename without extension
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :param storage: storage to use fo the job
    :return: Job object
    """

    j = b.new_job(name=f'Samtools-Sort-CRAM: {bam_filename_no_ext}')

    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        }
    )

    j.image(img)
    j.cpu(ncpu)
    j.memory(f'{memory}Gi')
    j.storage(f'{storage}Gi')
    j.command(
        f"""
        cd /io
        samtools sort -@{ncpu-1} -m 1G -o {j.output_cram['cram']} -O cram {input_cram['cram']}
        samtools index -@{ncpu-1} {j.output_cram['cram']} {j.output_cram['cram.crai']}
        """
    )

    # b.write_output(j.output_cram,
    #               f'{output_path}/intermediate/sorted-crams/{bam_filename_no_ext}')

    return j

