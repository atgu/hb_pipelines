__author__ = 'Sophie Parsa & Lindo Nkambule'


import argparse
import hailtop.batch as hb
import hailtop.fs as hfs
import os
import pandas as pd
from hailtop.batch.job import Job
from typing import Dict, List, Union


# GATK Best Practices
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
# https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf-fc-hg38.wdl
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
    "dbSNP_vcf": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
    "known_indels_sites_VCFs": [
        "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
    ],
    "contamination_sites_ud": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.contam.UD",
    "contamination_sites_bed": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.contam.bed",
    "contamination_sites_mu": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.contam.mu",
    }


def size(file: str):
    """
    Convert the size from bytes to GiB
    :param file: path to file, str
    :return: file size in GiB
    """

    file_info = hfs.stat(file)   # returns a named tuple
    size_bytes = file_info.size
    size_gigs = size_bytes / (1024 * 1024 * 1024)

    return size_gigs


# Assume out_dir and tmp_dir are different, and write every intermediate file to tmp_dir and final files to out_dir

# Some files are on CRAM version 2.0 is not supported which is not supported by GATK. We first convert them to
# BAM using samtools (which handles older CRAM versions) and the pass the BAM to GATK
def cram_to_bam(
        b: hb.batch.Batch,
        input_bam: hb.ResourceFile = None,
        fasta_reference: hb.ResourceGroup = None,
        output_bam_prefix: str = None,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: str = 'standard',
        ncpu: int = 8,
        tmp_dir: str = None
) -> Job:
    j = b.new_job(name=f'CramToBam: {output_bam_prefix}')

    j.cpu(ncpu)
    j.image(img)
    j.memory(memory)
    j.storage(f'{disk_size}Gi')

    j.command(f"""samtools index {input_bam}""")
    j.command(
        f"""
        samtools view -@{ncpu} -b -T {fasta_reference['fasta']} -o output_converted.bam {input_bam}
        """
    )

    j.command(f'mv output_converted.bam {j.ofile}')

    # b.write_output(j.ofile, f'{tmp_dir}/gatk_vc/cram_to_bam/{output_bam_prefix}.bam')
    b.write_output(j.ofile, f'{tmp_dir}/{output_bam_prefix}.bam')

    return j


# 1. Get read groups in a BAM file, so we can parallelize mapping
def get_read_groups(
        b: hb.batch.Batch,
        input_bam: hb.ResourceFile = None,
        output_bam_prefix: str = None,
        memory: str = 'lowmem',
        ncpu: int = 4,
        storage: Union[float, int] = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        tmp_dir: str = None
) -> Job:
    """
    Get read group IDs from a BAM file
    :param b: batch
    :param input_bam: input BAM file
    :param output_bam_prefix: output BAM filename
    :param memory: job memory
    :param ncpu: number of CPUs
    :param storage: storage to use for the job
    :param img: image to use for the job
    :param tmp_dir: output directory to write temporary file to
    :return: Job object
    """
    j = b.new_job(name=f'1.GetBamReadGroups: {output_bam_prefix}')

    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}')
    j.storage(f'{storage}Gi')
    j.command(
        f"""
        samtools view -@{ncpu} -H {input_bam} | grep ^@RG > rgs_tmp.txt
        """
    )

    # Get read group IDs only
    j.command(f"""
        awk -F'\t' '{{split($2, arr, ":"); print arr[2]}}' rgs_tmp.txt > {j.ofile}
    """)

    # write a file containing read group ID per line, so we can loop through it when mapping
    # b.write_output(j.ofile, f'{tmp_dir}/gatk_vc/read_group_ids/{output_bam_prefix}.rgs.txt')
    b.write_output(j.ofile, f'{tmp_dir}/{output_bam_prefix}.rgs.txt')

    return j


# https://gatk.broadinstitute.org/hc/en-us/articles/4403687183515--How-to-Generate-an-unmapped-BAM-from-FASTQ-or-aligned-BAM
# https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently
# 2. unmap reads
# set VALIDATION_STRINGENCY=VALIDATION_STRINGENCY to avoid issues like
# Read name HS27_10303:6:1206:18594:36926#21, Mapped mate should have mate reference name when data is paired-end but
# some reads are missing a pair
def bam_to_ubam(
        b: hb.batch.Batch,
        input_bam: hb.ResourceFile,
        output_bam_prefix: str = None,
        rg_ids: List[str] = None,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: str = 'standard',
        ncpu: int = 8,
        tmp_dir: str = None
) -> Job:
    """
    Unmap reads
    :param b: batch
    :param input_bam: input BAM file to be unmapped
    :param output_bam_prefix: BAM filename without extension
    :param rg_ids: list of unique read group IDs in the BAM file
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :param disk_size: disk size to use for the job
    :param tmp_dir: output directory to write temporary file to
    :return: Job object
    """
    j = b.new_job(name=f'2.BamToUbam: {output_bam_prefix}')

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j.storage(f'{disk_size}Gi')

    java_mem = ncpu * 4 - 10    # ‘lowmem’ ~1Gi/core, ‘standard’ ~4Gi/core, and ‘highmem’ ~7Gi/core in Hail Batch

    # set --OUTPUT_BY_READGROUP to true, so we can map each group in parallel
    j.command(
        f"""
        cd /io
        mkdir tmp/
        mkdir tmp/reverted/
        gatk --java-options -Xmx{java_mem}g RevertSam \
            -I {input_bam} \
            -O `pwd`/tmp/reverted \
            --MAX_DISCARD_FRACTION 0.005 \
            --ATTRIBUTE_TO_CLEAR XT \
            --ATTRIBUTE_TO_CLEAR XN \
            --ATTRIBUTE_TO_CLEAR AS \
            --ATTRIBUTE_TO_CLEAR OC \
            --ATTRIBUTE_TO_CLEAR OP \
            --OUTPUT_BY_READGROUP true \
            --SORT_ORDER queryname \
            --RESTORE_ORIGINAL_QUALITIES true \
            --REMOVE_DUPLICATE_INFORMATION true \
            --REMOVE_ALIGNMENT_INFORMATION true \
            --VALIDATION_STRINGENCY LENIENT \
            --TMP_DIR `pwd`/tmp
        """
    )

    for rg_id in rg_ids:
        j.command(f"cp tmp/reverted/{rg_id}.bam {j[f'{rg_id}']}")
        # b.write_output(j[f'{rg_id}'], f'{tmp_dir}/gatk_vc/unmapped_bams/{output_bam_prefix}/{rg_id}.unmapped.bam')
        b.write_output(j[f'{rg_id}'], f'{tmp_dir}/{output_bam_prefix}/{rg_id}.unmapped.bam')

    return j


# 3. Map reads to reference (still to be tested/optimized)
def sam_to_fastq_and_bwa_mem_and_mba(
        b: hb.batch.Batch,
        input_bam: hb.resource.InputResourceFile,
        output_bam_prefix: str = None,
        fasta_reference: hb.ResourceGroup = None,
        rg: str = None,
        compression_level: int = 2,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/lindonkambule/gatk-bwa:v1.0',
        memory: str = 'standard',
        ncpu: int = 8,
        tmp_dir: str = None
) -> Job:
    """
    Read unmapped BAM, convert to FASTQ and stream to BWA MEM for alignment, then stream to MergeBamAlignment
    :param b: batch
    :param input_bam: input uBAM file to be mapped
    :param output_bam_prefix: BAM filename without extension
    :param fasta_reference: reference genome files to be used in alignment
    :param compression_level: compression level
    :param rg: read group ID
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :param disk_size: disk size to use for the job
    :param tmp_dir: output directory to write temporary file to
    :return: Job object
    """
    j = b.new_job(name=f'3.SamToFastqAndBwaMemAndMba-{output_bam_prefix}: {rg}')

    j.declare_resource_group(
        output_bam={
            'bam': '{root}.bam'
        }
    )

    # input_bam = b.read_input(f'{tmp_dir}/gatk_vc/unmapped_bams/{output_bam_prefix}/{rg}.unmapped.bam')

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j.storage(f'{disk_size}Gi')
    java_mem = ncpu * 4 - 10    # ‘lowmem’ ~1Gi/core, ‘standard’ ~4Gi/core, and ‘highmem’ ~7Gi/core in Hail Batch

    j.command(f"""
        bwa_version=$(bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')
        bwa_commandline=$(echo bwa mem -K 100000000 -pt{ncpu} -v 3 -CH <(samtools view -H {input_bam}|grep ^@RG) -Y {fasta_reference['ref_fasta']})
    """)

    # https://lh3.github.io/2021/07/06/remapping-an-aligned-bam
    j.command(
        f"""
        set -o pipefail
        set -e
    
        # set the bash variable needed for the command-line
        # if ref_alt has data in it,
        if [ -s {fasta_reference['ref_alt']} ]; then
          samtools fastq -OT RG,BC {input_bam} |
          bwa mem -K 100000000 -pt{ncpu} -v 3 -CH <(samtools view -H {input_bam}|grep ^@RG) -Y {fasta_reference['ref_fasta']} - | \
          gatk --java-options "-Dsamjdk.compression_level={compression_level} -Xms5000m -Xmx{java_mem}g" MergeBamAlignment \
            --VALIDATION_STRINGENCY SILENT \
            --EXPECTED_ORIENTATIONS FR \
            --ATTRIBUTES_TO_RETAIN X0 \
            --ATTRIBUTES_TO_REMOVE NM \
            --ATTRIBUTES_TO_REMOVE MD \
            --ALIGNED_BAM /dev/stdin \
            --UNMAPPED_BAM {input_bam} \
            --OUTPUT {j.output_bam['bam']} \
            --REFERENCE_SEQUENCE {fasta_reference['ref_fasta']} \
            --PAIRED_RUN true \
            --SORT_ORDER "unsorted" \
            --IS_BISULFITE_SEQUENCE false \
            --ALIGNED_READS_ONLY false \
            --CLIP_ADAPTERS false \
            --MAX_RECORDS_IN_RAM 2000000 \
            --ADD_MATE_CIGAR true \
            --MAX_INSERTIONS_OR_DELETIONS -1 \
            --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
            --PROGRAM_RECORD_ID "bwamem" \
            --PROGRAM_GROUP_VERSION "$bwa_version" \
            --PROGRAM_GROUP_COMMAND_LINE "$bwa_commandline" \
            --PROGRAM_GROUP_NAME "bwamem" \
            --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
            --ALIGNER_PROPER_PAIR_FLAGS true \
            --UNMAP_CONTAMINANT_READS true \
            --ADD_PG_TAG_TO_READS false
    
        # else ref_alt is empty or could not be found
        else
          exit 1;
        fi
        """
    )

    # b.write_output(j.output_bam,
    #                f'{tmp_dir}/gatk_vc/mapped_bams/{output_bam_prefix}/{rg}.aligned.unsorted')
    b.write_output(j.output_bam,
                   f'{tmp_dir}/{output_bam_prefix}/{rg}.aligned.unsorted')

    return j


# 4. Mark duplicate reads to avoid counting non-independent observations
# We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
# to avoid having to spend time just merging BAM files.
def mark_duplicates(
        b: hb.batch.Batch,
        input_bams: List[hb.resource.InputResourceFile],
        output_bam_prefix: str = None,
        use_preemptible_worker: bool = None,
        compression_level: int = 2,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: float = 8,
        ncpu: int = 2,
        tmp_dir: str = None
) -> Job:
    """
    Mark duplicates
    :param b: batch
    :param input_bams: list of input BAM files. The BAMs are from the same individual, just split by read group
    :param output_bam_prefix: BAM filename without extension
    :param use_preemptible_worker: whether to use a preemptible worker
    :param compression_level: compression level
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :param disk_size: disk size to use for the job
    :param tmp_dir: output directory to write files to
    :return: Job object
    """
    j = b.new_job(name=f'4.MarkDuplicates: {output_bam_prefix}')

    j.declare_resource_group(
        output_bam={
            'bam': '{root}.bam',
        }
    )

    bams_rg = ' '.join([f"-I {f}" for f in input_bams])

    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}Gi')
    j._preemptible = use_preemptible_worker
    j.storage(f'{disk_size}Gi')
    j.command(
        f"""
        cd /io
        mkdir tmp/
        gatk --java-options "-Dsamjdk.compression_level={compression_level} -Xms4000m" MarkDuplicates \
            {bams_rg} \
            -O {j.output_bam['bam']} \
            -M {j.markdup_metrics} \
            --VALIDATION_STRINGENCY SILENT \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --ASSUME_SORT_ORDER queryname \
            --CLEAR_DT false \
            --ADD_PG_TAG_TO_READS false \
            --TMP_DIR `pwd`/tmp
        """
    )

    # b.write_output(j.output_bam,
    #                f'{tmp_dir}/gatk_vc/mark_duplicates/{output_bam_prefix}.aligned.unsorted.duplicates_marked')
    # b.write_output(j.markdup_metrics,
    #                f'{tmp_dir}/gatk_vc/mark_duplicates/{output_bam_prefix}.marked_dup_metrics.txt')
    b.write_output(j.output_bam,
                   f'{tmp_dir}/{output_bam_prefix}.aligned.unsorted.duplicates_marked')
    b.write_output(j.markdup_metrics,
                   f'{tmp_dir}/{output_bam_prefix}.marked_dup_metrics.txt')

    return j


# 5. SortSam
def picard_sort_bam(
        b: hb.batch.Batch,
        input_bam: hb.resource.InputResourceFile,
        output_bam_prefix: str = None,
        compression_level: int = 2,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: float = 5,
        ncpu: int = 1,
        disk_size: Union[float, int] = None,
        tmp_dir: str = None
) -> Job:
    """
    Sort BAM file, using Picard, by coordinate order and fix tag values for NM and UQ (slower than samtools)
    :param b: batch
    :param input_bam: input BAM file to be sorted
    :param output_bam_prefix: BAM filename without extension
    :param compression_level: compression level
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :param disk_size: disk size to use for the job
    :param tmp_dir: output directory to write files to
    :return: Job object
    """
    j = b.new_job(name=f'5.SortSam: {output_bam_prefix}')

    j.declare_resource_group(
        output_bam={
            'bam': '{root}.bam',
            'bam.bai': '{root}.bam.bai'
        }
    )

    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}Gi')
    j.storage(f'{disk_size}Gi')
    # --CREATE_INDEX was not working, hence the indexing is a separate command
    j.command(
        f"""cd /io
        mkdir tmp/
        gatk --java-options "-Dsamjdk.compression_level={compression_level} -Xms5000m" SortSam \
            -I {input_bam} \
            -O {j.output_bam['bam']} \
            -SO coordinate \
            --MAX_RECORDS_IN_RAM 300000 \
            --TMP_DIR `pwd`/tmp
        """
    )
    j.command(
        f"""cd /io
        gatk --java-options "-Xms5000m" BuildBamIndex \
            -I {j.output_bam['bam']}\
            -O {j.output_bam['bam.bai']} \
            --TMP_DIR `pwd`/tmp
        """
    )

    j.command("ls")

    if tmp_dir:
        # b.write_output(j.output_bam,
        #                f'{tmp_dir}/gatk_vc/sorted_bams/{output_bam_prefix}.aligned.duplicate_marked.sorted')
        b.write_output(j.output_bam,
                       f'{tmp_dir}/{output_bam_prefix}.aligned.duplicate_marked.sorted')

    return j


# 6. Check that the fingerprints of separate readgroups all match
def cross_check_fingerprints(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup,
        output_bam_prefix: str,
        haplotype_database_file: hb.ResourceGroup,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: str = 'standard',
        ncpu: int = 8,
        out_dir: str = None
) -> Job:
    j = b.new_job(name=f'6.CrossCheckFingerprints: {output_bam_prefix}')

    j.cpu(ncpu)
    j.image(img)
    j.memory(memory)
    j.storage(f'{disk_size}Gi')

    j.command(
        f"""
        gatk --java-options "-Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms2000m" \
            CrosscheckFingerprints \
            --OUTPUT {output_bam_prefix}.crosscheck \
            --HAPLOTYPE_MAP {haplotype_database_file} \
            --EXPECT_ALL_READ_GROUPS_TO_MATCH true \
            --INPUT {input_bam['bam']} \
            --LOD_THRESHOLD -20.0
        """
    )

    if out_dir:
        j.command(f'mv {output_bam_prefix}.crosscheck {j.output}')
        b.write_output(j.output,
                       f'{out_dir}/{output_bam_prefix}.crosscheck')

    return j


# 7. Check contamination
# We do not need to adjust the FREEMIX at this step since we'll use CHARR later
def check_contamination(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup,
        output_bam_prefix: str,
        fasta_reference: hb.ResourceGroup,
        contamination_sites: hb.ResourceGroup,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/griffan/verifybamid2:latest',
        memory: int = 8,
        ncpu: int = 8,
        out_dir: str = None
) -> Job:
    j = b.new_job(name=f'7.CheckContamination: {output_bam_prefix}')
    j.declare_resource_group(
        output={
            'selfSM': '{root}.selfSM',
            'Ancestry': '{root}.Ancestry'
        }
    )
    j.image(img)
    j.memory(f'{memory}Gi')
    j.cpu(ncpu)
    j.storage(f'{disk_size}Gi')

    j.command(
        f"""
        set -e
    
        # creates a .selfSM file, a TSV file with 2 rows, 19 columns.
        # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
        VerifyBamID \
        --NumPC 4 \
        --Output {j.output} \
        --BamFile {input_bam['bam']} \
        --Reference {fasta_reference['ref_fasta']} \
        --UDPath {contamination_sites['ud']} \
        --MeanPath {contamination_sites['mu']} \
        --BedPath {contamination_sites['bed']}
        """
    )

    j.command(
        f"""
        rows=$(tail -n +2 {j.output['selfSM']} | wc -l) # get number of rows, excluding header in the count
        rows=$((rows + 0)) # convert to int
    
        # for a single-sample BAM, number of rows should be two: header and sample stats row
        if [ "$rows" != 1 ]; then
            echo "Found $rows rows in .selfSM file. Was expecting exactly 1. This is an error" >&2
            exit 2
        fi
        
        cols=($(awk -F ' ' '{{print $7, $8, $9}}' {j.output['selfSM']} | tail -n +2))
        
        if [ "${{cols[1]}}" == 0 ] && [ "${{cols[2]}}" == 0 ]; then
            echo "Error: Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf)." >&2
            exit 1
        fi
        
        echo "FREEMIX=${{cols[0]}}"
        """
    )

    # use FREEMIX=$(tail -n +2 {j.output['selfSM']} | awk '{{print $7}}') to extract freemix estimate downstream (HC)
    # b.write_output(j.output,
    #                f'{out_dir}/gatk_vc/VerifyBamID/{output_bam_prefix}')
    b.write_output(j.output, f'{out_dir}/{output_bam_prefix}')

    return j


# Create list of sequences for scatter-gather parallelization
# Only ran this function once since it only relies on the reference dictionary file which is unchanged
def create_sequence_grouping_tsv(
        b: hb.Batch = None,
        img: str = 'docker.io/hailgenetics/python-dill:3.9-slim',
) -> Job:
    """
    Generate sets of intervals for scatter-gathering over chromosomes
    :param b: Batch object to add jobs to
    :param img: image to use for the job
    :return: a Job object
    """
    j = b.new_job("""CreateSequenceGroupingTSV""")
    ref_dict = b.read_input('gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict')
    j.image(img)
    j.memory('lowmem')
    j.storage('1G')

    j.command(f"cp {ref_dict} reference.dict")

    j.command(f'''
    python -c '
with open("reference.dict", "r") as ref_dict_file:
    sequence_tuple_list = []
    longest_sequence = 0
    for line in ref_dict_file:
        if line.startswith("@SQ"):
            line_split = line.split("""\t""")
            # (Sequence_Name, Sequence_Length)
            sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
    longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
# We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
# the last element after a :, so we add this as a sacrificial element.
hg38_protection_tag = ":1+"
# initialize the tsv string with the first sequence
tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
temp_size = sequence_tuple_list[0][1]
for sequence_tuple in sequence_tuple_list[1:]:
    if temp_size + sequence_tuple[1] <= longest_sequence:
        temp_size += sequence_tuple[1]
        tsv_string += """\t""" + sequence_tuple[0] + hg38_protection_tag
    else:
        tsv_string += """\n""" + sequence_tuple[0] + hg38_protection_tag
        temp_size = sequence_tuple[1]
# add the unmapped sequences as a separate line to ensure that they are recalibrated as well
with open("sequence_grouping.txt", "w") as tsv_file:
    tsv_file.write(tsv_string)
    tsv_file.close()

tsv_string += """\n""" + "unmapped"

with open("sequence_grouping_with_unmapped.txt", "w") as tsv_file_with_unmapped:
    tsv_file_with_unmapped.write(tsv_string)
    tsv_file_with_unmapped.close()
'
      '''
              )

    j.command(
        f"""
        cp sequence_grouping.txt {j['sequence_grouping']}
        cp sequence_grouping_with_unmapped.txt {j['sequence_grouping_with_unmapped']}
    """)

    out_dir = 'gs://h3africa/variant_calling_resources'
    b.write_output(j['sequence_grouping'], f'{out_dir}/hg38_sequence_grouping.txt')
    b.write_output(j['sequence_grouping_with_unmapped'], f'{out_dir}/hg38_sequence_grouping_with_unmapped.txt')

    return j


# 8. BaseRecalibrator
def base_recalibrator(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup,
        fasta_reference: hb.ResourceGroup,
        output_bam_prefix: str = None,
        sequence_group_interval: Union[List[bytes], str] = None,
        utils: Dict = None,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: float = 6,
        ncpu: int = 2
) -> Job:
    """
    Generate Base Quality Score Recalibration (BQSR) model
    :param b: batch
    :param input_bam: input BAM file
    :param fasta_reference: reference files
    :param output_bam_prefix: BAM filename without extension
    :param sequence_group_interval: interval containing one or more genomic intervals over which to operate
    :param utils: utils dictionary with DBSNP VCF filepath
    :param disk_size: disk size to use for the job
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :return: Job object
    """
    j = b.new_job(name=f'8.BaseRecalibrator: {output_bam_prefix}')

    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}Gi')
    j.storage(f'{disk_size}Gi')

    known_indel_sites = ' '.join([f'--known-sites {file}' for file in utils['known_indels_sites_VCFs']])
    sequence_group_interval_cmd = ' '.join([f'-L {interval}' for interval in sequence_group_interval])

    # PrintGCTimeStamps and similar flags has been removed from java
    j.command(
        f"""
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms4000m" BaseRecalibrator \
            -R {fasta_reference['ref_fasta']} \
            -I {input_bam['bam']} \
            --use-original-qualities \
            -O {j.recalibration_report}\
            --known-sites {utils['dbSNP_vcf']} \
            {known_indel_sites} \
            {sequence_group_interval_cmd}
        """
    )

    return j


# 9. GatherBQSRReports
def gather_bqsr_reports(
        b: hb.Batch,
        input_bqsr_reports: List[hb.ResourceFile],
        output_bam_prefix: str = None,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: int = 3500,
        ncpu: int = 1
) -> Job:
    """
    Combine multiple recalibration tables from scattered BaseRecalibrator runs
    :param b: Batch object to add jobs to
    :param input_bqsr_reports: list of reports from scattered BaseRecalibrator runs
    :param output_bam_prefix: BAM filename without extension
    :param disk_size: disk size to use for the job
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :return: a Job object with one output j.output_bqsr_report
    """
    j = b.new_job(f'9.GatherBqsrReports: {output_bam_prefix}')
    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}M')
    j.storage(f'{disk_size}Gi')

    inputs_cmdl = ' '.join([f'-I {report}' for report in input_bqsr_reports])
    j.command(
        f"""set -euo pipefail
        gatk --java-options "-Xms3000m" GatherBQSRReports \
          {inputs_cmdl} \
          -O {j.output_bqsr_report}
        """
    )

    return j


# 10. ApplyBQSR
def apply_bqsr(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup,
        fasta_reference: hb.ResourceGroup,
        output_bam_prefix: str = None,
        bsqr_report: hb.ResourceFile = None,
        sequence_group_interval: Union[List[bytes], str] = None,
        compression_level: int = 2,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: int = 3500,
        ncpu: int = 2
) -> Job:
    """
    Apply Base Quality Score Recalibration (BQSR) model
    :param b: batch
    :param input_bam: input BAM file
    :param fasta_reference: reference files
    :param output_bam_prefix: BAM filename without extension
    :param bsqr_report: BQSR model/report file
    :param sequence_group_interval: interval containing one or more genomic intervals over which to operate
    :param compression_level: compression level
    :param disk_size: storage to use for the job
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :return: Job object
    """

    j = b.new_job(name=f'10.ApplyBQSR: {output_bam_prefix}')

    j.declare_resource_group(
        output_bam={
            'bam': '{root}.bam',
            'md5': '{root}.bam.md5',
        }
    )

    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}M')
    j.storage(f'{disk_size}Gi')

    sequence_group_interval_cmd = ' '.join([f'-L {interval}' for interval in sequence_group_interval])

    # PrintGCTimeStamps and similar flags has been removed from java
    # "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms4000m" BaseRecalibrator
    j.command(
        f"""
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level={compression_level} -Xms3000m" \
            ApplyBQSR \
            --create-output-bam-md5 \
            --add-output-sam-program-record \
            -R {fasta_reference['ref_fasta']} \
            -I {input_bam['bam']} \
            --use-original-qualities \
            -O {j.output_bam['bam']} \
            -bqsr {bsqr_report} \
            --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
            {sequence_group_interval_cmd}

    """
    )

    return j


# 11. GatherBamFiles
def gather_bam_files(
        b: hb.batch.Batch,
        input_bams: List[hb.ResourceGroup] = None,
        output_bam_prefix: str = None,
        compression_level: int = 2,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: float = 3,
        ncpu: int = 8,
        disk_size: Union[float, int] = None,
        tmp_dir: str = None
) -> Job:
    """
    Merge the recalibrated BAM files resulting from by-interval recalibration
    :param b: batch
    :param input_bams: input BAM files from scattered ApplyBQSR runs to be merged
    :param output_bam_prefix: BAM filename without extension
    :param compression_level: compression level
    :param disk_size: storage to use for the job
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :param tmp_dir: temporary directory
    :return: Job object
    """
    input_cmdl = ' '.join([f'--INPUT {b["bam"]}' for b in input_bams])

    j = b.new_job(name=f'11.GatherBamFiles: {output_bam_prefix}')

    j.declare_resource_group(
        gathered_bams={
            'bam': '{root}.bam',
            'bam.bai': '{root}.bai',
            'bam.md5': '{root}.bam.md5'
        }
    )

    j.image(img)
    j.cpu(ncpu)
    j.memory(f'{memory}Gi')
    j.storage(f'{disk_size}Gi')
    j.command(
        f"""
        gatk --java-options "-Dsamjdk.compression_level={compression_level} -Xms2000m" GatherBamFiles \
              {input_cmdl} \
              --OUTPUT {j.gathered_bams['bam']} \
              --CREATE_INDEX true \
              --CREATE_MD5_FILE true
        """
    )

    # b.write_output(j.gathered_bams,
    #                f'{tmp_dir}/gatk_vc/bqsr/{output_bam_prefix}.aligned.duplicate_marked.sorted.bqsr')
    b.write_output(j.gathered_bams,
                   f'{tmp_dir}/{output_bam_prefix}.aligned.duplicate_marked.sorted.bqsr')

    return j


# 12. Convert BAM file to CRAM format
# Tested on a ~22GB BAM. 1 cpu = 40 minutes ($0.0536), 4 cpus = 12 mins ($0.0252)
def convert_to_cram(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup = None,
        fasta_reference: hb.ResourceGroup = None,
        output_bam_prefix: str = None,
        memory: int = 3,
        ncpu: int = 4,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        out_dir: str = None
) -> Job:
    """
    Convert BAM file to CRAM format
    :param b: batch
    :param input_bam: BAM file to be converted to CRAM
    :param fasta_reference: reference files
    :param output_bam_prefix: BAM filename without extension
    :param memory: job memory
    :param ncpu: number of CPUs
    :param disk_size: storage to use for the job
    :param img: image to use for the job
    :param out_dir: output directory
    :return: Job object
    """
    j = b.new_job(name=f'12.ConvertToCram: {output_bam_prefix}')

    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
            'md5': '{root}.cram.md5'
        }
    )

    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}Gi')
    j.storage(f'{disk_size}Gi')
    j.command(
        f"""        
        samtools view -@{ncpu} -C -T {fasta_reference['ref_fasta']} {input_bam['bam']} | \
            tee {j.output_cram['cram']} | \
            md5sum | awk '{{print $1}}' > {j.output_cram['md5']}
        
        # Create REF_CACHE. Used when indexing a CRAM
        seq_cache_populate.pl -root ./ref/cache {fasta_reference['ref_fasta']}
        export REF_PATH=:
        export REF_CACHE=./ref/cache/%2s/%2s/%s
        
        samtools index {j.output_cram['cram']} {j.output_cram['cram.crai']}
        """
    )

    # DEFINITELY write this out to a bucket
    # b.write_output(j.output_cram, f'{out_dir}/gatk_vc/crams/{output_bam_prefix}')
    b.write_output(j.output_cram, f'{out_dir}/{output_bam_prefix}')

    return j


# CheckPreValidation output is used an input in ValidateSamFile, and CheckPreValidation requires a file from
# CollectMultipleMetrics, hence the following (13 and 14) were included in the pipeline
# 13. CollectMultipleMetrics
def collect_aggregation_metrics(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup,
        fasta_reference: hb.ResourceGroup = None,
        output_bam_prefix: str = None,
        memory: int = 8,
        ncpu: int = 1,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/broadinstitute/gatk:latest'
) -> Job:

    j = b.new_job(name=f'13.CollectAggregationMetrics: {output_bam_prefix}')

    j.declare_resource_group(
        output={
            'alignment_summary_metrics': '{root}.alignment_summary_metrics' # only file required in CheckPreValidation
        }
    )

    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}Gi')
    j.storage(f'{disk_size}Gi')
    # java_mem = int((ncpu * memory) * 0.5)

    j.command(
        f"""
        gatk --java-options "-Xms5000m" CollectMultipleMetrics \
            -I {input_bam['bam']} \
            -R {fasta_reference['ref_fasta']} \
            --OUTPUT {j.output} \
            --ASSUME_SORTED true \
            --PROGRAM "null" \
            --PROGRAM "CollectAlignmentSummaryMetrics" \
            --METRIC_ACCUMULATION_LEVEL "null" \
            --METRIC_ACCUMULATION_LEVEL "SAMPLE" \
            --METRIC_ACCUMULATION_LEVEL "LIBRARY"
    """)

    return j


# 14. Check whether the data has massively high duplication or chimerism rates
def check_pre_validation(
        b: hb.Batch = None,
        output_bam_prefix: str = None,
        duplication_metrics: hb.resource.InputResourceFile = None,
        chimerism_metrics: hb.ResourceFile = None,
        max_duplication_in_reasonable_sample: float = 0.30,
        max_chimerism_in_reasonable_sample: float = 0.15,
        img: str = 'docker.io/hailgenetics/python-dill:3.9-slim',
        tmp_dir: str = None
) -> Job:
    """
    Check whether the data has massively high duplication or chimerism rates
    :param b: Batch object to add jobs to
    :param output_bam_prefix: BAM filename without extension
    :param duplication_metrics: file from 3.MarkDuplicates with duplication stats
    :param chimerism_metrics: file from 13.CollectAggregationMetrics with chimerism stats
    :param max_duplication_in_reasonable_sample: max reasonable duplication for a sample
    :param max_chimerism_in_reasonable_sample: max reasonable chimerism for a sample
    :param img: image to use for the job
    :param tmp_dir: output directory to write files to
    :return: a Job object with a single output j.intervals of type ResourceGroup
    """
    j = b.new_job(f'14.CheckPreValidation: {output_bam_prefix}')
    j.image(img)
    j.memory('lowmem')
    j.storage('1G')

    j.command(
        f'''
        grep -A 1 PERCENT_DUPLICATION {duplication_metrics} > duplication.csv
        grep -A 3 PCT_CHIMERAS {chimerism_metrics} | grep -v OF_PAIR > chimerism.csv
    ''')

    j.command(f'''
    python -c '
import csv
with open("duplication.csv") as dupfile:
    reader = csv.DictReader(dupfile, delimiter="""\t""")
    for row in reader:
        with open("duplication_value.txt","w") as file:
            file.write(row["PERCENT_DUPLICATION"])
            file.close()

with open("chimerism.csv") as chimfile:
    reader = csv.DictReader(chimfile, delimiter="""\t""")
    for row in reader:
        with open("chimerism_value.txt","w") as file:
            file.write(row["PCT_CHIMERAS"])
            file.close()
'
      '''
              )

    # we only need the boolean for validating the CRAM file
    j.command(f"""
        dup=$(awk '{{print $1}}' duplication_value.txt)
        chim=$(awk '{{print $1}}' chimerism_value.txt)
        
        # Image does not have bc, hence we're using awk
        if (( $(echo $dup {max_duplication_in_reasonable_sample} | awk '{{if ($1 > $2) print 1;}}') )) || (( $(echo $chim {max_chimerism_in_reasonable_sample}  | awk '{{if ($1 > $2) print 1;}}') )); then
            is_outlier=true
        else
            is_outlier=false
        fi
        
        echo $is_outlier > is_outlier_data.txt
    """)

    j.command(
        f"""
        cp is_outlier_data.txt {j.is_outlier_data}
    """)

    # b.write_output(j.is_outlier_data,
    #                f'{tmp_dir}/gatk_vc/check_pre_validation/{output_bam_prefix}.is.outlier.txt')
    b.write_output(j.is_outlier_data,
                   f'{tmp_dir}/{output_bam_prefix}.is.outlier.txt')

    return j


# 15. ValidateSamFile
def validate_samfile(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup,
        fasta_reference: hb.ResourceGroup,
        ignore: List[str] = ["MISSING_TAG_NM"],
        output_bam_prefix: str = None,
        max_output: int = 1000000000,
        is_outlier_data: hb.resource.InputResourceFile = None,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: int = 3500,
        ncpu: int = 2
) -> Job:
    """
    Validate the CRAM file
    :param b: batch
    :param input_bam: input CRAM file
    :param fasta_reference: reference files. dict, fast, and fasta index required
    :param ignore: list of validation error types to ignore
    :param output_bam_prefix: BAM filename without extension
    :param max_output: the maximum number of lines output in verbose mode
    :param is_outlier_data: file containing a boolean indication whether sample in an outlier or not
    :param disk_size: storage to use for the job
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :return: Job object
    """

    j = b.new_job(name=f'15.ValidateCramFile: {output_bam_prefix}')

    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}M')
    j.storage(f'{disk_size}Gi')

    ignore_cmd = ' '.join([f'--IGNORE {err}' for err in ignore]) if ignore else '''--IGNORE "null"'''

    j.command(
        f"""
        is_outlier=$(awk '{{print $1}}' {is_outlier_data})
        
        gatk --java-options "-Xms6000m" ValidateSamFile \
            -I {input_bam['cram']} \
            -O {j.report} \
            -R {fasta_reference['ref_fasta']} \
            --MAX_OUTPUT {max_output} \
            {ignore_cmd} \
            --MODE VERBOSE \
            --SKIP_MATE_VALIDATION $is_outlier \
            --IS_BISULFITE_SEQUENCED false
    """
    )

    return j


# Run a single sample (either one BAM/CRAM or multiple BAMs/CRAMs) through the entire processing workflow
def pre_process_bam(
        b: hb.Batch,
        sample: str,
        input_bams: List[str],
        is_cram_old: bool = False,
        fasta_reference: hb.ResourceGroup = None,
        bwa_reference_files: hb.ResourceGroup = None,
        haplotype_database_file: hb.ResourceGroup = None,
        contamination_sites: hb.ResourceGroup = None,
        ref_size: Union[float, int] = None,
        bwa_ref_size: Union[float, int] = None,
        bwa_disk_multiplier: float = 2.5,
        additional_disk: float = 20.0,
        md_disk_multiplier: float = 2.25,
        sort_sam_disk_multiplier: float = 3.25,
        out_dir: str = None,
        tmp_dir: str = None,
        utils: Dict = None
):
    """
    Process one BAM file
    :param b: Batch object to add jobs to
    :param sample: sample ID for the BAM file(s) to be processed
    :param input_bams: input BAM file(s) belonging to the same sample to be processed
    :param is_cram_old: whether CRAM is version 2.0. If True, a CRAM-to-BAM step is added
    :param fasta_reference: reference files. dict, fast, and fasta index required
    :param bwa_reference_files: reference files required by BWA
    :param haplotype_database_file: haplotype database file
    :param contamination_sites: ResourceGroup of the contamination sites files
    :param ref_size: size of reference files (fasta, index, dict)
    :param bwa_ref_size: size of BWA reference files
    :param bwa_disk_multiplier: full path to transmitted singletons VCF file and its index
    :param additional_disk: full path to sibling singletons VCF file and its index
    :param md_disk_multiplier: full path to transmitted singletons VCF file and its index
    :param sort_sam_disk_multiplier: full path to sibling singletons VCF file and its index
    :param out_dir: directory to write final output files to
    :param tmp_dir: directory to write intermediate files to
    :param utils: dictionary containing paths to public resources and commonly used
    :return: a final Job, and a path to the VCF with VQSR annotations
    """
    bam_idx = [i for i in range(len(input_bams))]   # indices for file output names
    mapped_bam_total_size = 0   # Sum the bam sizes to approximate the aggregated bam size
    bams_rg_mapped = []

    for bam, idx in zip(input_bams, bam_idx):
        bam_file = b.read_input(bam)
        bam_prefix = f'{sample}_{idx}'
        unmapped_bam_size = size(bam)

        # Convert CRAM to BAM first if old CRAM version
        if is_cram_old:
            ref_path = 'gs://gnomaf/genome_reference/hs37d5'
            ref_tmp = b.read_input_group(**{'fasta': f'{ref_path}.fa',
                                            'idx': f'{ref_path}.fa.fai',
                                            'ref_dict': f'{ref_path}.dict'})
            cram_to_bam(
                b=b,
                input_bam=bam_file,
                output_bam_prefix=bam_prefix,
                fasta_reference=ref_tmp,
                disk_size=unmapped_bam_size + 2.0*unmapped_bam_size + additional_disk,
                tmp_dir=f'{tmp_dir}/cram_to_bam/{sample}'
            )
            b.run()

            input_bam = b.read_input(f'{tmp_dir}/cram_to_bam/{sample}/{bam_prefix}.bam')
        else:
            input_bam = bam_file

        if not hfs.exists(f'{tmp_dir}/read_group_ids/{sample}/{bam_prefix}.rgs.txt'):
            get_read_groups(
                b=b,
                input_bam=input_bam,
                output_bam_prefix=bam_prefix,
                storage=unmapped_bam_size,
                tmp_dir=f'{tmp_dir}/read_group_ids/{sample}'
            )
            b.run()

        bam_rg_ids = hfs.open(f'{tmp_dir}/read_group_ids/{sample}/{bam_prefix}.rgs.txt').readlines()
        rgs = [l.strip() for l in bam_rg_ids]   # rgs have a \n at the end, hence the .strip()

        # first check if unmapped BAMs exist
        ubams_exist = [hfs.exists(f'{tmp_dir}/unmapped_bams/{sample}/{bam_prefix}/{rg}.unmapped.bam') for rg in rgs]

        if not all(ubams_exist):
            # unmap BAM to multiple uBAM files split by read group
            bam_to_ubam(
                b=b,
                input_bam=input_bam,
                output_bam_prefix=bam_prefix,
                rg_ids=rgs,
                disk_size=unmapped_bam_size + unmapped_bam_size*2.5 + additional_disk,
                tmp_dir=f'{tmp_dir}/unmapped_bams/{sample}'
            )
            b.run()     # call this, so we can be able to get file sizes of readgroup BAM files

        # map uBAM files in parallel
        # b.run() doesn't allow reading files from a job executed by one b.run() in a job executed by a different b.run()
        unmapped_bams = [b.read_input(f'{tmp_dir}/unmapped_bams/{sample}/{bam_prefix}/{rg}.unmapped.bam') for
                         rg in rgs]
        unmapped_bam_sizes = [size(f'{tmp_dir}/unmapped_bams/{sample}/{bam_prefix}/{rg}.unmapped.bam') for rg in rgs]

        bams_exist = [hfs.exists(f'{tmp_dir}/unmapped_bams/{sample}/{bam_prefix}/{rg}.aligned.unsorted.bam') for
                      rg in rgs]

        if not all(bams_exist):
            # list of BAM files that do not exist
            not_mapped = [(unmapped_bams[i], rgs[i], unmapped_bam_sizes[i]) for i, val in enumerate(bams_exist) if not val]

            bams_rg_mapped = [
                sam_to_fastq_and_bwa_mem_and_mba(
                    b=b,
                    input_bam=bam,
                    output_bam_prefix=bam_prefix,
                    rg=rg,
                    fasta_reference=bwa_reference_files,
                    disk_size=round(unmapped_bam_size+bwa_ref_size+(bwa_disk_multiplier*unmapped_bam_size)+additional_disk),
                    tmp_dir=f'{tmp_dir}/mapped_bams/{sample}'
                ).output_bam
                for bam, rg, unmapped_bam_size in not_mapped
            ]
            b.run()

        mapped_bam_total_size += sum([size(f'{tmp_dir}/mapped_bams/{sample}/{bam_prefix}/{rg}.aligned.unsorted.bam')
                                      for rg in rgs])
        bams_mapped = [b.read_input(f'{tmp_dir}/mapped_bams/{sample}/{bam_prefix}/{rg}.aligned.unsorted.bam') for
                       rg in rgs]
        bams_rg_mapped.append(bams_mapped)

    # MarkDuplicates and SortSam currently take too long for preemptibles if the input data is too large
    use_preemptible = not mapped_bam_total_size > 110.0

    if not hfs.exists(f'{tmp_dir}/mark_duplicates/{sample}.aligned.unsorted.duplicates_marked.bam'):
        mark_duplicates(
            b=b,
            input_bams=bams_rg_mapped,
            output_bam_prefix=sample,
            use_preemptible_worker=use_preemptible,
            disk_size=(md_disk_multiplier * mapped_bam_total_size) + additional_disk,
            tmp_dir=f'{tmp_dir}/mark_duplicates'
        )
        b.run()

    bam_md = b.read_input(f'{tmp_dir}/mark_duplicates/{sample}.aligned.unsorted.duplicates_marked.bam')
    agg_bam_size = size(f'{tmp_dir}/mark_duplicates/{sample}.aligned.unsorted.duplicates_marked.bam')

    sort_sam_bam = picard_sort_bam(
        b=b,
        input_bam=bam_md,
        output_bam_prefix=sample,
        disk_size=(sort_sam_disk_multiplier * agg_bam_size) + additional_disk
    ).output_bam

    # CrosscheckFingerprints requires a haplotype map file
    if haplotype_database_file:
        cross_check_fingerprints(b=b,
                                 input_bam=sort_sam_bam,
                                 haplotype_database_file=haplotype_database_file,
                                 output_bam_prefix=sample,
                                 disk_size=agg_bam_size + additional_disk
                                 )

    check_contamination(
        b=b,
        input_bam=sort_sam_bam,
        output_bam_prefix=f'{sample}.preBqsr',
        fasta_reference=fasta_reference,
        contamination_sites=contamination_sites,
        disk_size=agg_bam_size + ref_size + additional_disk,
        out_dir=out_dir
    )

    if not hfs.exists(f'{tmp_dir}/bqsr/{sample}.aligned.duplicate_marked.sorted.bqsr.bam'):
        seq_grouping = hfs.open('gs://h3africa/variant_calling_resources/hg38_sequence_grouping.txt').readlines()
        seq_grouping_subgroup = [l.split('\t') for l in seq_grouping]
        potential_bqsr_divisor = len(seq_grouping_subgroup) - 10
        bqsr_divisor = potential_bqsr_divisor if potential_bqsr_divisor > 1 else 1

        bqsr_reports = [
            base_recalibrator(
                b=b,
                input_bam=sort_sam_bam,
                fasta_reference=fasta_reference,
                output_bam_prefix=sample,
                sequence_group_interval=subgroup,
                utils=utils,
                disk_size=agg_bam_size + ref_size + additional_disk
            ).recalibration_report
            for subgroup in seq_grouping_subgroup
        ]

        gathered_bqsr_report = gather_bqsr_reports(
            b=b,
            input_bqsr_reports=[report for report in bqsr_reports],
            output_bam_prefix=sample,
            disk_size=5
        ).output_bqsr_report

        seq_grouping_with_unmapped = hfs.open('gs://h3africa/variant_calling_resources/hg38_sequence_grouping_with_unmapped.txt').readlines()
        seq_grouping_with_unmapped_subgroup = [l.split('\t') for l in seq_grouping_with_unmapped]

        # Apply the recalibration model by interval
        recalibrated_bams = [
            apply_bqsr(
                b=b,
                input_bam=sort_sam_bam,
                fasta_reference=fasta_reference,
                output_bam_prefix=sample,
                bsqr_report=gathered_bqsr_report,
                sequence_group_interval=subgroup,
                disk_size=agg_bam_size + (agg_bam_size / bqsr_divisor) + ref_size + additional_disk
            ).output_bam
            for subgroup in seq_grouping_with_unmapped_subgroup
        ]

        # Merge the recalibrated BAM files resulting from by-interval recalibration
        gather_bam_files(
            b=b,
            input_bams=recalibrated_bams,
            output_bam_prefix=sample,
            disk_size=(2 * agg_bam_size) + additional_disk,
            tmp_dir=f'{tmp_dir}/bqsr'
        )

        b.run()     # BQSR bins the qualities which makes a significantly smaller bam. Get binned file size

    gathered_bams_path = f'{tmp_dir}/bqsr/{sample}.aligned.duplicate_marked.sorted.bqsr'
    gathered_bams = b.read_input_group(**{'bam': f'{gathered_bams_path}.bam',
                                       'bam.bai': f'{gathered_bams_path}.bam.bai'})
    binned_qual_bam_size = size(f'{gathered_bams_path}.bam')

    # QC the final BAM some more (no such thing as too much QC)
    agg_metrics = collect_aggregation_metrics(
        b=b,
        input_bam=gathered_bams,
        fasta_reference=fasta_reference,
        output_bam_prefix=sample,
        disk_size=binned_qual_bam_size + ref_size + additional_disk
    ).output

    mark_dup_metrics = b.read_input(f'{tmp_dir}/mark_duplicates/{sample}.marked_dup_metrics.txt')
    check_pre_validation(
        b=b,
        output_bam_prefix=sample,
        duplication_metrics=mark_dup_metrics,
        chimerism_metrics=agg_metrics['alignment_summary_metrics'],
        tmp_dir=f'{tmp_dir}/check_pre_validation'
    )

    # Convert the final merged recalibrated BAM file to CRAM format
    if not hfs.exists(f'{out_dir}/gatk_vc/crams/{sample}.cram'):
        convert_to_cram(
            b=b,
            input_bam=gathered_bams,
            fasta_reference=fasta_reference,
            output_bam_prefix=sample,
            disk_size=(2 * binned_qual_bam_size) + ref_size + additional_disk,
            out_dir=f'{out_dir}/crams/'
        )
    b.run()

    # Validate the CRAM file
    cram_file = b.read_input_group(**{'cram': f'{out_dir}/crams/{sample}.cram',
                                      'cram.crai': f'{out_dir}/crams/{sample}.cram.crai'})
    cram_size = size(f'{out_dir}/crams/{sample}.cram')

    pre_val_metric = b.read_input(f'{tmp_dir}/check_pre_validation/{sample}.is.outlier.txt')

    validate_samfile(
        b=b,
        input_bam=cram_file,
        fasta_reference=fasta_reference,
        output_bam_prefix=sample,
        is_outlier_data=pre_val_metric,
        disk_size=cram_size + ref_size + additional_disk
    )
    b.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-files', type=str, required=True)
    parser.add_argument('--out-dir', type=str, required=True)
    parser.add_argument('--old-cram-version', action='store_true')
    parser.add_argument('--tmp-dir', type=str, required=True)
    parser.add_argument('--billing-project', type=str, required=True)

    args = parser.parse_args()

    backend = hb.ServiceBackend(
        billing_project=args.billing_project,
        regions=['us-central1'],
        remote_tmpdir=args.tmp_dir,
    )

    batch = hb.Batch(
        'BAM-Processing',
        backend=backend,
    )

    ref_fasta = batch.read_input_group(**{'ref_fasta': inputs['ref_fasta'],
                                          'ref_fasta_index': inputs['ref_fasta_index'],
                                          'ref_dict': inputs['ref_dict']})

    # BWA requires additional reference file
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
    ref_size = size(inputs['ref_fasta']) + size(inputs['ref_fasta_index']) + size(inputs['ref_dict'])
    bwa_ref_size = ref_size + size(inputs['ref_alt']) + size(inputs['ref_amb']) + size(inputs['ref_ann']) + size(inputs['ref_bwt']) + size(inputs['ref_pac']) + size(inputs['ref_sa'])

    bams = pd.read_csv(args.input_files, sep='\t', header=None, names=['id', 'files'])
    files = [(sample, paths) for sample, paths in zip(bams['id'], bams['files'])]

    print(files)

    for sample_id, sample_bams in files:
        sample_bam_file_paths = sample_bams.split(',')

        pre_process_bam(
            b=batch,
            sample=sample_id,
            input_bams=sample_bam_file_paths,
            is_cram_old=args.old_cram_version,
            fasta_reference=ref_fasta,
            bwa_reference_files=ref_fasta_bwa,
            contamination_sites=contamination_sites,
            ref_size=ref_size,
            bwa_ref_size=bwa_ref_size,
            out_dir=args.out_dir,
            tmp_dir=args.tmp_dir,
            utils=inputs
        )

    batch.run()


if __name__ == '__main__':
    main()
