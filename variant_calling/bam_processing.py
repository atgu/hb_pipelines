__author__ = 'Sophie Parsa & Lindo Nkambule'


import hail as hl
import hailtop.batch as hb
import hailtop.fs as hfs
from hailtop.batch.job import Job
from typing import Dict, List, Union


# GATK Best Practices
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
# https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf-fc-hg38.wdl
inputs = {"ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
          "ref_fasta_index": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
          "ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
          "wgs_calling_interval_list": "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list",
          "calling_int_list": "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list",
          "eval_int_list": "gs://gcp-public-data--broad-references/hg38/v0/wgs_evaluation_regions.hg38.interval_list",
          "dbsnp_resource_vcf": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
          "known_indels_sites_VCFs": [
              "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
              "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
          ],
          "contamination_sites_ud": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.contam.UD",
          "contamination_sites_bed": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.contam.bed",
          "contamination_sites_mu": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.contam.mu",
          "gatk_img": "us.gcr.io/broad-gatk/gatk:4.2.0.0",
          "haplotype_scatter_count": 50,
          "break_bands_at_multiples_of": 1000000,
          "hc_contamination": 0}


def get_file_size(file):
    """Get file size"""

    file_info = hl.utils.hadoop_stat(file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)

    return size_gigs


# https://gatk.broadinstitute.org/hc/en-us/articles/4403687183515--How-to-Generate-an-unmapped-BAM-from-FASTQ-or-aligned-BAM
# 1. unmap reads
def bam_to_ubam(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup,
        output_bam_prefix: str = None,
        disk_size: int = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: str = 'standard',
        ncpu: int = 8,
        out_dir: str = None
) -> Job:
    """
    Unmap reads
    :param b: batch
    :param input_bam: input BAM file to be unmapped
    :param output_bam_prefix: BAM filename without extension
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :param disk_size: disk size to use for the job
    :param out_dir: output directory to write files to
    :return: Job object
    """
    j = b.new_job(name=f'1.BamToUbam: {output_bam_prefix}')

    j.declare_resource_group(
        output_bam={
            'bam': '{root}.bam'
        }
    )

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j.storage(f'{disk_size}Gi')

    java_mem = ncpu * 4 - 10    # ‘lowmem’ ~1Gi/core, ‘standard’ ~4Gi/core, and ‘highmem’ ~7Gi/core in Hail Batch

    j.command(
        f"""
        cd /io
        mkdir tmp/
        gatk --java-options -Xmx{java_mem}g RevertSam \
            -I {input_bam['bam']} \
            -O {j.output_bam['bam']} \
            --SANITIZE true \
            --MAX_DISCARD_FRACTION 0.005 \
            --ATTRIBUTE_TO_CLEAR XT \
            --ATTRIBUTE_TO_CLEAR XN \
            --ATTRIBUTE_TO_CLEAR AS \
            --ATTRIBUTE_TO_CLEAR OC \
            --ATTRIBUTE_TO_CLEAR OP \
            --SORT_ORDER queryname \
            --RESTORE_ORIGINAL_QUALITIES true \
            --REMOVE_DUPLICATE_INFORMATION true \
            --REMOVE_ALIGNMENT_INFORMATION true \
            --TMP_DIR `pwd`/tmp
        """
    )

    if out_dir:
        b.write_output(j.output_bam,
                       f'{out_dir}/gatk_vc/unmapped_bams/{output_bam_prefix}.unmapped')

    return j


# 2. Map reads to reference (still to be tested/optimized)
def sam_to_fastq_and_bwa_mem_and_mba(
        b: hb.batch.Batch,
        input_bam: hb.ResourceFile,
        output_bam_prefix: str = None,
        ref_genome: hb.ResourceGroup = None,
        bwa_commandline: str = 'bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta',
        compression_level: int = 2,
        disk_size: int = None,
        img: str = 'docker.io/lindonkambule/gatk-bwa:v1.0',
        memory: str = 'standard',
        ncpu: int = 16,
        out_dir: str = None
) -> Job:
    """
    Read unmapped BAM, convert to FASTQ and stream to BWA MEM for alignment, then stream to MergeBamAlignment
    :param b: batch
    :param input_bam: input uBAM file to be mapped
    :param output_bam_prefix: BAM filename without extension
    :param ref_genome: reference genome files to be used in alignment
    :param compression_level: compression level
    :param bwa_commandline: BWA command-line to be used ro map the BAM file
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :param disk_size: disk size to use for the job
    :param out_dir: output directory to write files to
    :return: Job object
    """
    j = b.new_job(name=f'2.SamToFastqAndBwaMemAndMba: {output_bam_prefix}')

    j.declare_resource_group(
        output_bam={
            'bam': '{root}.bam',
            'log': '{root}.bwa.stderr.log'
        }
    )

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j.storage(f'{disk_size}Gi')
    java_mem = ncpu * 4 - 10    # ‘lowmem’ ~1Gi/core, ‘standard’ ~4Gi/core, and ‘highmem’ ~7Gi/core in Hail Batch

    # not setting set -o pipefail here because /bwa has a rc=1 and we don't want to allow rc=1 to succeed because
    # the sed may also fail with that error and that is something we actually want to fail on.
    j.command(f"""bwa_version=$(bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')""")
    j.command(
        f"""
        set -o pipefail
        set -e
    
        # set the bash variable needed for the command-line
        bash_ref_fasta={ref_genome['ref_fasta']}
        # if ref_alt has data in it,
        if [ -s {ref_genome['ref_alt']} ]; then
          gatk --java-options "-Xms5000m -Xmx8000m" SamToFastq \
            --INPUT {input_bam} \
            --FASTQ /dev/stdout \
            --INTERLEAVE true \
            --NON_PF true | \
          {bwa_commandline} /dev/stdin - 2> >(tee {output_bam_prefix}.bwa.stderr.log >&2) | \
          gatk --java-options "-Dsamjdk.compression_level={compression_level} -Xms5000m -Xmx{java_mem}g" MergeBamAlignment \
            --VALIDATION_STRINGENCY SILENT \
            --EXPECTED_ORIENTATIONS FR \
            --ATTRIBUTES_TO_RETAIN X0 \
            --ATTRIBUTES_TO_REMOVE NM \
            --ATTRIBUTES_TO_REMOVE MD \
            --ALIGNED_BAM /dev/stdin \
            --UNMAPPED_BAM {input_bam} \
            --OUTPUT {j.output_bam['bam']} \
            --REFERENCE_SEQUENCE {ref_genome['ref_fasta']} \
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
            --PROGRAM_GROUP_COMMAND_LINE "{bwa_commandline}" \
            --PROGRAM_GROUP_NAME "bwamem" \
            --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
            --ALIGNER_PROPER_PAIR_FLAGS true \
            --UNMAP_CONTAMINANT_READS true \
            --ADD_PG_TAG_TO_READS false
    
          grep -m1 "read .* ALT contigs" {output_bam_prefix}.bwa.stderr.log | \
          grep -v "read 0 ALT contigs"
    
        # else ref_alt is empty or could not be found
        else
          exit 1;
        fi
        """
    )

    # This step is compute-heavy, we definitely want to write this file out. Make sure out_dir is specified
    if out_dir:
        b.write_output(j.output_bam,
                       f'{out_dir}/gatk_vc/mapped_bams/{output_bam_prefix}.aligned.unsorted')

    return j


# 3. Mark duplicate reads to avoid counting non-independent observations
def mark_duplicates(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup,
        output_bam_prefix: str = None,
        compression_level: int = 2,
        disk_size: int = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: float = 8,
        ncpu: int = 2,
        out_dir: str = None,
) -> Job:
    """
    Mark duplicates
    :param b: batch
    :param input_bam: input CRAM file
    :param output_bam_prefix: BAM filename without extension
    :param compression_level: compression level
    :param img: image to use for the job
    :param memory: job memory
    :param ncpu: number of CPUs
    :param disk_size: disk size to use for the job
    :param out_dir: output directory to write files to
    :return: Job object
    """
    j = b.new_job(name=f'3.MarkDuplicates: {output_bam_prefix}')

    j.declare_resource_group(
        output_bam={
            'bam': '{root}.bam',
        }
    )

    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}Gi')
    j.storage(f'{disk_size}Gi')
    j.command(
        f"""
        cd /io
        mkdir tmp/
        gatk --java-options "-Dsamjdk.compression_level={compression_level} -Xms4000m" MarkDuplicates \
            -I {input_bam['bam']} \
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

    if out_dir:
        b.write_output(j.output_bam,
                       f'{out_dir}/gatk_vc/mark_duplicates/{output_bam_prefix}.aligned.unsorted.duplicates_marked')

    return j


# 4. SortSam
def picard_sort_bam(
        b: hb.batch.Batch,
        input_bam: Union[str, hb.ResourceGroup],
        output_bam_prefix: str = None,
        compression_level: int = 2,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: float = 8,
        ncpu: int = 2,
        disk_size: int = None,
        out_dir: str = None
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
    :param out_dir: output directory to write files to
    :return: Job object
    """
    j = b.new_job(name=f'4.SortSam: {output_bam_prefix}')

    j.declare_resource_group(
        output_bam={
            'bam': '{root}.bam',
            'bam.bai': '{root}.bam.bai',
            'bam.md5': '{root}.bam.md5',
        }
    )

    j.cpu(ncpu)
    j.image(img)
    j.memory(f'{memory}Gi')
    j.storage(f'{disk_size}Gi')
    j.command(
        f"""cd /io
        mkdir tmp/
        gatk --java-options "-Dsamjdk.compression_level={compression_level} -Xms5000m" SortSam \
            -I {input_bam['bam']} \
            -O {j.output_bam['bam']} \
            -SO coordinate \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true \
            --MAX_RECORDS_IN_RAM 300000 \
            --TMP_DIR `pwd`/tmp
        """
    )

    if out_dir:
        b.write_output(j.output_bam,
                       f'{out_dir}/gatk_vc/sorted_bams/{output_bam_prefix}.aligned.duplicate_marked.sorted')

    return j


# 5. Check that the fingerprints of separate readgroups all match
def cross_check_fingerprints(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup,
        output_bam_prefix: str,
        haplotype_database_file: hb.ResourceGroup,
        disk_size: int = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: str = 'standard',
        ncpu: int = 8,
        out_dir: str = None
) -> Job:
    j = b.new_job(name=f'5.CrossCheckFingerprints: {output_bam_prefix}')

    j.cpu(ncpu)
    j.image(img)
    j.memory(memory)
    j.storage(f'{disk_size}Gi')

    j.command(
        f"""
        gatk --java-options "-Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms2000m" \
            CrosscheckReadGroupFingerprints \
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
                       f'{out_dir}/gatk_vc/cross_check_fingerprints/{output_bam_prefix}.crosscheck')

    return j


# 6. Check contamination
# We do not need to adjust the FREEMIX at this step since we'll use CHARR later
def check_contamination(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup,
        output_bam_prefix: str,
        ref_genome: hb.ResourceGroup,
        contamination_sites_ud: hb.ResourceFile,
        contamination_sites_mu: hb.ResourceFile,
        contamination_sites_bed: hb.ResourceFile,
        disk_size: float = None,
        img: str = 'docker.io/griffan/verifybamid2:latest',
        memory: int = 8,
        ncpu: int = 8
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
    j.storage(f'{disk_size + 5}Gi')

    j.command(
        f"""
        set -e
    
        # creates a .selfSM file, a TSV file with 2 rows, 19 columns.
        # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
        VerifyBamID \
        --NumPC 4 \
        --Output {j.output} \
        --BamFile {input_bam['bam']} \
        --Reference {ref_genome['ref_fasta']} \
        --UDPath {contamination_sites_ud} \
        --MeanPath {contamination_sites_mu} \
        --BedPath {contamination_sites_bed}
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

    # use FREEMIX=$(tail -n +2 {j.output['selfSM']} | awk '{{print $7}}') to extract freemix estimate downstream

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


# 7. IntervalListTools
# Perform variant calling on the sub-intervals, and then gather the results
def split_interval_list(
        b: hb.Batch = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        utils: Dict = None,
) -> Job:
    """
    Break the calling interval_list into sub-intervals
    :param b: Batch object to add jobs to
    :param img: image to use for the job
    :param utils: a dictionary containing resources (file paths and arguments) to be used to split genome
    :return: a Job object with a single output j.intervals of type ResourceGroup
    """
    j = b.new_job(f"""7.Make {utils['haplotype_scatter_count']} intervals""")
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
        -I {utils['wgs_calling_interval_list']} \
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


# 8. BaseRecalibrator
def base_recalibrator(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup,
        fasta_reference: hb.ResourceGroup,
        output_bam_prefix: str = None,
        sequence_group_interval: List[bytes] = None,
        utils: Dict = None,
        disk_size: int = None,
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

    j.command(
        f"""
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" BaseRecalibrator \
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
def gather_tranches(
        b: hb.Batch,
        input_bqsr_reports: List[hb.ResourceFile],
        output_bam_prefix: str = None,
        disk_size: int = None,
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
        sequence_group_interval: str = None,
        compression_level: int = 2,
        disk_size: int = None,
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

    j.command(
        f"""
        gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
            -XX:+PrintGCDetails -Xloggc:gc_log.log \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=${compression_level} -Xms3000m" \
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
        disk_size: int = None,
        out_dir: str = None
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
    :param out_dir: output directory
    :return: Job object
    """
    input_cmdl = ' '.join([f'--INPUT {b["bam"]}' for b in input_bams])

    j = b.new_job(name=f'11.GatherBamFiles: {output_bam_prefix}')

    j.declare_resource_group(
        gathered_bams={
            'bam': '{root}.bam',
            'index': '{root}.bai',
            'md5': '{root}.bam.md5'
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

    if out_dir:
        b.write_output(j.gathered_bams,
                       f'{out_dir}/gatk_vc/bqsr/{output_bam_prefix}.aligned.duplicate_marked.sorted.bqsr')

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
        disk_size: int = None,
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
        samtools view -@ {ncpu} -C -T {fasta_reference['ref_fasta']} {input_bam['bam']} | \
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
    if out_dir:
        b.write_output(j.output_cram, f'{out_dir}/gatk_vc/crams/{output_bam_prefix}')

    return j


# CheckPreValidation output is used an input in ValidateSamFile, and CheckPreValidation requires a file from
# CollectMultipleMetrics, hence the following (13 and 14) were included in the pipeline
# 13. CollectMultipleMetrics
def collect_aggregation_metrics(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup,
        fasta_reference: hb.ResourceGroup = None,
        output_bam_prefix: str = None,
        memory: int = 7,
        ncpu: int = 4,
        disk_size: int = None,
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
    java_mem = int((ncpu * memory) * 0.5)    # set max heap space to only up one-half of available mem

    j.command(
        f"""
        gatk --java-options "-Xms5000m -Xmx{java_mem}g" CollectMultipleMetrics \
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
        duplication_metrics: hb.ResourceFile = None,
        chimerism_metrics: hb.ResourceFile = None,
        max_duplication_in_reasonable_sample: float = 0.30,
        max_chimerism_in_reasonable_sample: float = 0.15,
        img: str = 'docker.io/hailgenetics/python-dill:3.9-slim',
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
        cp is_outlier_data.txt {j['is_outlier_data']}
    """)

    return j


# 15. ValidateSamFile
def validate_samfile(
        b: hb.batch.Batch,
        input_bam: hb.ResourceGroup,
        fasta_reference: hb.ResourceGroup,
        ignore: List[str] = ["MISSING_TAG_NM"],
        output_bam_prefix: str = None,
        max_output: int = 1000000000,
        is_outlier_data: hb.ResourceFile = None,
        disk_size: int = None,
        img: str = 'docker.io/broadinstitute/gatk:latest',
        memory: int = 3500,
        ncpu: int = 2
) -> Job:
    """
    Apply Base Quality Score Recalibration (BQSR) model
    :param b: batch
    :param input_bam: input BAM file
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
            -I {input_bam['bam']} \
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
