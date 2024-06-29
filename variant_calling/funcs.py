__author__ = 'Lindo Nkambule'

import hailtop.batch as hb
import hailtop.fs as hfs
from hailtop.batch.job import Job
from typing import List, Union


# Functions
def size(file: str):
    """
    Convert the size from bytes to GiB
    :param file: path to file, str
    :return: file size in GiB
    """

    file_info = hfs.stat(file)  # returns a named tuple
    size_bytes = file_info.size
    size_gigs = size_bytes / (1024 * 1024 * 1024)

    return size_gigs


# https://github.com/gatk-workflows/five-dollar-genome-analysis-pipeline/blob/b33decd14b550ad157a948ad720166040e435890/tasks/SplitLargeReadGroup.wdl#L41
def split_by_num_reads(
        b: hb.batch.Batch,
        input_bam: hb.resource.InputResourceFile,
        output_bam_prefix: str = None,
        disk_size: Union[float, int] = None,
        n_reads: int = 48000000,
        ncpu: int = 4,
        docker: str = 'docker.io/broadinstitute/gatk:latest',
        tmp_dir: str = None
) -> Job:

    j = b.new_job(name=f'SplitSamByNumberOfReads: {output_bam_prefix}')

    j.image(docker)
    j.memory('standard')
    j.cpu(ncpu)
    j.storage(f'{disk_size}Gi')

    java_mem = ncpu * 4

    j.command(
        f"""cd /io
        mkdir tmp/
        gatk --java-options "-Xms3000m -Xmx{java_mem}g" SortSam \
            -I {input_bam} \
            -O sorted_tmp.bam \
            -SO queryname \
            --TMP_DIR `pwd`/tmp
        """
    )

    # "Splitting a coordinate sorted bam may result in invalid bams that do not always contain each read's mate in
    # the same bam" warning when BAM is coordinate sorted. Sort by queryname lexicographically
    j.command(
        f"""
        cd /io
        mkdir tmp/reverted/
                    
        total_reads=$(samtools view -c sorted_tmp.bam)
        gatk --java-options "-Xms3000m -Xmx{java_mem}g" SplitSamByNumberOfReads \
            -I sorted_tmp.bam \
            -O `pwd`/tmp/reverted \
            --SPLIT_TO_N_READS {n_reads} \
            --TOTAL_READS_IN_INPUT $total_reads \
            --TMP_DIR `pwd`/tmp

        ls `pwd`/tmp/reverted
    """
    )

    j.command(f'mv `pwd`/tmp/reverted {j.reverted}')
    b.write_output(j.reverted, f'{tmp_dir}/{output_bam_prefix}')

    return j


def revert_bam_to_ubam(
        b: hb.batch.Batch,
        input_bam: hb.resource.InputResourceFile = None,
        output_bam_prefix: str = None,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/lindonkambule/gatk-bwa:v1.0',
        memory: str = 'standard',
        ncpu: int = 4,
        tmp_dir: str = None
) -> Job:
    j = b.new_job(name=f'BamToUbam: {output_bam_prefix}')

    j.cpu(ncpu)
    j.image(img)
    j.memory(memory)
    j.storage(f'{disk_size}Gi')

    java_mem = ncpu * 4 - 10    # ‘lowmem’ ~1Gi/core, ‘standard’ ~4Gi/core, and ‘highmem’ ~7Gi/core in Hail Batch

    j.command(
        f"""
        cd /io
        mkdir tmp/
        mkdir tmp/reverted/
        gatk --java-options -Xmx{java_mem}g RevertSam \
            -I {input_bam} \
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
        gatk --java-options "-Xms3000m -Xmx3600m" ValidateSamFile \
            -I {j.output_bam} \
            -MODE SUMMARY
        """
    )

    b.write_output(j.output_bam, f'{tmp_dir}/{output_bam_prefix}.unmapped.bam')

    return j


# Some files are on CRAM version 2.0 which is not supported by GATK. We first convert them to BAM using samtools
# (which handles older CRAM versions)
def revert_cram_to_ubam(
        b: hb.batch.Batch,
        input_cram: hb.resource.InputResourceFile = None,
        old_ref_files: hb.ResourceGroup = None,
        output_bam_prefix: str = None,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/lindonkambule/gatk-bwa:v1.0',
        memory: str = 'standard',
        ncpu: int = 4,
        tmp_dir: str = None
) -> Job:
    j = b.new_job(name=f'CramToUbam: {output_bam_prefix}')

    j.cpu(ncpu)
    j.image(img)
    j.memory(memory)
    j.storage(f'{disk_size}Gi')

    java_mem = ncpu * 4 - 10    # ‘lowmem’ ~1Gi/core, ‘standard’ ~4Gi/core, and ‘highmem’ ~7Gi/core in Hail Batch

    j.command(
        f"""
        cd /io

        samtools view -@{ncpu} -b -T {old_ref_files['ref_fasta']} -o tmp.bam {input_cram}
        samtools index -@{ncpu} tmp.bam tmp.bam.bai
        """
    )

    j.command(
        f"""
        cd /io
        mkdir tmp/
        mkdir tmp/reverted/
        gatk --java-options -Xmx{java_mem}g RevertSam \
            -I tmp.bam \
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
            --VALIDATION_STRINGENCY LENIENT \
            --TMP_DIR `pwd`/tmp

        rm tmp.bam*
        """
    )

    j.command(
        f"""
        gatk --java-options "-Xms3000m -Xmx{java_mem}g" ValidateSamFile \
            -I {j.output_bam} \
            -MODE SUMMARY
        """
    )

    b.write_output(j.output_bam, f'{tmp_dir}/{output_bam_prefix}.unmapped.bam')

    return j


def sam_to_fastq_and_bwa_mem_and_mba(
        b: hb.batch.Batch,
        input_ubam: hb.resource.InputResourceFile = None,
        bwa_ref_files: hb.ResourceGroup = None,
        compression_level: int = 2,
        output_bam_prefix: str = None,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/lindonkambule/gatk-bwa:v1.0',
        memory: str = '16G',
        ncpu: int = 16,
        tmp_dir: str = None
) -> Job:
    j = b.new_job(name=f'SamToFastqAndBwaMemAndMba: {output_bam_prefix}')

    j.cpu(ncpu)
    j.image(img)
    j.memory(memory)
    j.storage(f'{disk_size}Gi')

    java_mem = ncpu * 16 - 10    # ‘lowmem’ ~1Gi/core, ‘standard’ ~4Gi/core, and ‘highmem’ ~7Gi/core in Hail Batch

    j.command(f"""
        bwa_version=$(bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')
        bwa_commandline=$(echo bwa mem -K 100000000 -pt{ncpu} -v 3 -Y {bwa_ref_files['ref_fasta']}) > {output_bam_prefix}.sam
    """)

    # Write BWA output to a SAM first since MergeBamAlignment couldn't handle file from stdin, it gave error below:
    # Error parsing text SAM file. Not enough fields; File /dev/stdin; Line ...
    j.command(
        f"""
            set -o pipefail
            set -e
        
            # set the bash variable needed for the command-line
            # if ref_alt has data in it,
            if [ -s {bwa_ref_files['ref_alt']} ]; then
              samtools fastq -OT RG,BC {input_ubam} |
              bwa mem -K 100000000 -pt{ncpu} -v 3 -Y {bwa_ref_files['ref_fasta']} /dev/stdin - > {output_bam_prefix}.sam
              gatk --java-options "-Dsamjdk.compression_level={compression_level} -Xms3000m -Xmx{java_mem}g" MergeBamAlignment \
                --VALIDATION_STRINGENCY SILENT \
                --EXPECTED_ORIENTATIONS FR \
                --ATTRIBUTES_TO_RETAIN X0 \
                --ATTRIBUTES_TO_REMOVE NM \
                --ATTRIBUTES_TO_REMOVE MD \
                --ALIGNED_BAM {output_bam_prefix}.sam \
                --UNMAPPED_BAM {input_ubam} \
                --OUTPUT {j.output_bam} \
                --REFERENCE_SEQUENCE {bwa_ref_files['ref_fasta']} \
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

    b.write_output(j.output_bam, f'{tmp_dir}/{output_bam_prefix}.aligned.unsorted.bam')

    return j


# Mark duplicate reads to avoid counting non-independent observations
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
        memory: str = '8G',
        ncpu: int = 2,
        tmp_dir: str = None
) -> Job:
    """
    Mark duplicates
    """
    j = b.new_job(name=f'MarkDuplicates: {output_bam_prefix}')

    bams_rg = ' '.join([f"-I {f}" for f in input_bams])

    j.cpu(ncpu)
    j.image(img)
    j.memory(memory)
    j._preemptible = use_preemptible_worker
    j.storage(f'{disk_size}Gi')
    j.command(
        f"""
        cd /io
        mkdir tmp/
        gatk --java-options "-Dsamjdk.compression_level={compression_level} -Xms3000m -Xmx6000m" MarkDuplicates \
            {bams_rg} \
            -O {j.output_bam} \
            -M {j.markdup_metrics} \
            --VALIDATION_STRINGENCY SILENT \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --ASSUME_SORT_ORDER queryname \
            --CLEAR_DT false \
            --ADD_PG_TAG_TO_READS false \
            --TMP_DIR `pwd`/tmp
        """
    )

    b.write_output(j.output_bam,
                   f'{tmp_dir}/{output_bam_prefix}.aligned.unsorted.duplicates_marked.bam')
    b.write_output(j.markdup_metrics,
                   f'{tmp_dir}/{output_bam_prefix}.marked_dup_metrics.txt')

    return j


# SortSam
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
    """
    j = b.new_job(name=f'SortSam: {output_bam_prefix}')

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
        gatk --java-options "-Dsamjdk.compression_level={compression_level} -Xmx5000m" SortSam \
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

    if tmp_dir:
        b.write_output(j.output_bam,
                       f'{tmp_dir}/{output_bam_prefix}.aligned.duplicate_marked.sorted')

    return j


# Check that the fingerprints of separate readgroups all match
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
    j = b.new_job(name=f'CrossCheckFingerprints: {output_bam_prefix}')

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


# Check contamination
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
    j = b.new_job(name=f'CheckContamination: {output_bam_prefix}')
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


# BaseRecalibrator
def base_recalibrator(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup,
        fasta_reference: hb.ResourceGroup,
        output_bam_prefix: str = None,
        sequence_group_interval: Union[List[bytes], str] = None,
        disk_size: Union[float, int] = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: float = 6,
        ncpu: int = 2
) -> Job:
    """
    Generate Base Quality Score Recalibration (BQSR) model
    """
    j = b.new_job(name=f'BaseRecalibrator: {output_bam_prefix}')

    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}Gi')
    j.storage(f'{disk_size}Gi')

    sites = [
        "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
        "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
    ]

    known_sites = ' '.join([f'--known-sites {file}' for file in sites])
    sequence_group_interval_cmd = ' '.join([f'-L {interval}' for interval in sequence_group_interval])

    # PrintGCTimeStamps and similar flags has been removed from java
    j.command(
        f"""
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms4000m" BaseRecalibrator \
            -R {fasta_reference['ref_fasta']} \
            -I {input_bam['bam']} \
            --use-original-qualities \
            -O {j.recalibration_report} \
            --gcs-max-retries 50 \
            {known_sites} \
            {sequence_group_interval_cmd}
        """
    )

    return j


# GatherBQSRReports
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
    """
    j = b.new_job(f'GatherBqsrReports: {output_bam_prefix}')
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


# ApplyBQSR
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
    """

    j = b.new_job(name=f'ApplyBQSR: {output_bam_prefix}')

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


# GatherBamFiles
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
    """
    input_cmdl = ' '.join([f'--INPUT {b["bam"]}' for b in input_bams])

    j = b.new_job(name=f'GatherBamFiles: {output_bam_prefix}')

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

    b.write_output(j.gathered_bams,
                   f'{tmp_dir}/{output_bam_prefix}.aligned.duplicate_marked.sorted.bqsr')

    return j


# Convert BAM file to CRAM format
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
    """
    j = b.new_job(name=f'ConvertToCram: {output_bam_prefix}')

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
    b.write_output(j.output_cram, f'{out_dir}/{output_bam_prefix}')

    return j


# CheckPreValidation output is used an input in ValidateSamFile, and CheckPreValidation requires a file from
# CollectMultipleMetrics, hence the following (13 and 14) were included in the pipeline
# CollectMultipleMetrics
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

    j = b.new_job(name=f'CollectAggregationMetrics: {output_bam_prefix}')

    j.declare_resource_group(
        output={
            'alignment_summary_metrics': '{root}.alignment_summary_metrics' # only file required in CheckPreValidation
        }
    )

    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}Gi')
    j.storage(f'{disk_size}Gi')

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
    """
    j = b.new_job(f'CheckPreValidation: {output_bam_prefix}')
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

    b.write_output(j.is_outlier_data,
                   f'{tmp_dir}/{output_bam_prefix}.is.outlier.txt')

    return j


# ValidateSamFile
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
    """

    j = b.new_job(name=f'ValidateCramFile: {output_bam_prefix}')

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
            -R {fasta_reference['ref_fasta']} \
            --MAX_OUTPUT {max_output} \
            {ignore_cmd} \
            --MODE VERBOSE \
            --SKIP_MATE_VALIDATION $is_outlier \
            --IS_BISULFITE_SEQUENCED false
    """
    )

    return j

