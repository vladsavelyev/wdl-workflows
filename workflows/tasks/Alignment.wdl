version 1.0

import "structs/DNASeqStructs.wdl" as Structs
import "Qc.wdl" as QC
import "BamProcessing.wdl" as Processing
import "Utilities.wdl" as Utils

workflow Alignment {

  input {
    Input inp
    DNASeqSingleSampleReferences references
    PapiSettings papi_settings

    Boolean check_contamination = true
    Boolean check_fingerprints = true

    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu

    String cross_check_fingerprints_by
    File haplotype_database_file
    Float lod_threshold
    Boolean hard_clip_reads = false
    
    Boolean to_cram = false
    String? subset_region
  }

  Float cutoff_for_large_rg_in_gb = 7.5

  # We don't partition input BAM for distributed alignment, to avoid extra copying
  # and extra deduplication+merge step with a large overhead.
  # However, if you want shared parallelilsm, consider adding 
  # bazam sharded parallelism like in this pipeline:
  # https://github.com/Oshlack/STRetch/blob/c5345e5dea4adfde790befb9903ec2d81ed5b2c1/pipelines/pipeline_stages.groovy#L101

  if (inp.bam_or_cram_or_fastq1 == sub(inp.bam_or_cram_or_fastq1, ".cram$", "") + ".cram" ||
    inp.bam_or_cram_or_fastq1 == sub(inp.bam_or_cram_or_fastq1, ".bam$", "") + ".bam") {
  
    call BwaFromBamOrCram {
      input:
        bam_or_cram = inp.bam_or_cram_or_fastq1,
        bai_or_crai = inp.bai_or_crai_or_fastq2,
        sample_name = inp.sample_name,
        output_bam_basename = inp.base_file_name,
        reference_fasta = references.reference_fasta,
        preemptible_tries = papi_settings.preemptible_tries,
        duplicate_metrics_fname = inp.base_file_name + ".duplicate_metrics",
        to_cram = to_cram,
        subset_region = subset_region,
    }
  }

  if (inp.bam_or_cram_or_fastq1 != sub(inp.bam_or_cram_or_fastq1, ".cram$", "") + ".cram" &&
    inp.bam_or_cram_or_fastq1 != sub(inp.bam_or_cram_or_fastq1, ".bam$", "") + ".bam") {
  
    call BwaFromFastq {
      input:
        fastq1 = inp.bam_or_cram_or_fastq1,
        fastq2 = inp.bai_or_crai_or_fastq2,
        sample_name = inp.sample_name,
        output_bam_basename = inp.base_file_name,
        reference_fasta = references.reference_fasta,
        preemptible_tries = papi_settings.preemptible_tries,
        duplicate_metrics_fname = inp.base_file_name + ".duplicate_metrics",
        to_cram = to_cram
    }
  }
  
  File mapped_file = select_first([BwaFromBamOrCram.output_file, BwaFromFastq.output_file])
  File mapped_indx = select_first([BwaFromBamOrCram.output_indx, BwaFromFastq.output_indx])
  Float mapped_file_size = size(mapped_file, "GiB")

  if (defined(haplotype_database_file) && check_fingerprints) {
    # Check identity of fingerprints across readgroups
    call QC.CrossCheckFingerprints as CrossCheckFingerprints {
      input:
        input_bams = [ mapped_file ],
        input_bam_indexes = [ mapped_indx ],
        haplotype_database_file = haplotype_database_file,
        metrics_filename = inp.base_file_name + ".crosscheck",
        total_input_size = mapped_file_size,
        lod_threshold = lod_threshold,
        cross_check_by = cross_check_fingerprints_by,
        preemptible_tries = papi_settings.agg_preemptible_tries,
        ref_dict = references.reference_fasta.ref_dict,
        ref_fasta = references.reference_fasta.ref_fasta,
        ref_fasta_index = references.reference_fasta.ref_fasta_index
    }
  }

  if (check_contamination) {
    # Estimate level of cross-sample contamination
    call Processing.CheckContamination as CheckContamination {
      input:
        input_bam = mapped_file,
        input_bam_index = mapped_indx,
        contamination_sites_ud = contamination_sites_ud,
        contamination_sites_bed = contamination_sites_bed,
        contamination_sites_mu = contamination_sites_mu,
        ref_fasta = references.reference_fasta.ref_fasta,
        ref_fasta_index = references.reference_fasta.ref_fasta_index,
        output_prefix = inp.base_file_name,
        preemptible_tries = papi_settings.agg_preemptible_tries,
        contamination_underestimation_factor = 0.75
    }
  }

  # Outputs that will be retained when execution is complete
  output {
    File? cross_check_fingerprints_metrics = CrossCheckFingerprints.cross_check_fingerprints_metrics

    File? selfSM = CheckContamination.selfSM
    Float? contamination = CheckContamination.contamination

    File duplicate_metrics = select_first([BwaFromBamOrCram.duplicate_metrics, BwaFromFastq.duplicate_metrics])

    File output_file = mapped_file
    File output_indx = mapped_indx
  }
  meta {
    allowNestedInputs: true
  }
}
  
task BwaFromFastq {
  input {
    File fastq1
    File fastq2
    String sample_name
    String output_bam_basename
    String duplicate_metrics_fname

    # reference_fasta.ref_alt is the .alt file from bwa-kit
    # (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative".
    ReferenceFasta reference_fasta

    Int preemptible_tries
    Boolean to_cram = false
  }
  
  String output_format = if to_cram then "cram" else "bam"
  
  Int bwa_cpu = 25
  Int bamsormadup_cpu = 6
  Int total_cpu = bwa_cpu + bamsormadup_cpu

  String rg_line = "@RG\\tID:~{sample_name}\\tSM:~{sample_name}"
  
  String output_file = if to_cram then "~{output_bam_basename}.cram" else "~{output_bam_basename}.bam"
  String output_indx = if to_cram then "~{output_bam_basename}.crai" else "~{output_bam_basename}.bai"
  
  # BWA command options:
  # -K     process INT input bases in each batch regardless of nThreads (for reproducibility)
  # -v3    minimum score to output [30]
  # -t16   threads
  # -Y     use soft clipping for supplementary alignments
  # -R     read group header line such as '@RG\tID:foo\tSM:bar'
  # -M     mark shorter split hits as secondary
  command <<<
    set -o pipefail
    set -ex

    (while true; do df -h; pwd; du -sh *; free -m; sleep 300; done) &
    
    bwa mem -K 100000000 -v3 -t~{bwa_cpu} -Y -R '~{rg_line}' \
      ~{reference_fasta.ref_fasta} ~{fastq1} ~{fastq2} \
      2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
    bamsormadup inputformat=sam threads=~{bamsormadup_cpu} SO=coordinate \
      M=~{duplicate_metrics_fname} \
      outputformat=sam | \
    samtools view -T ~{reference_fasta.ref_fasta} \
      -O ~{output_format} \
      -o ~{output_file}
    
    samtools index -@~{total_cpu} ~{output_file} ~{output_indx}

    df -h; pwd; du -sh *
  >>>
  
  runtime {
    # docker: "australia-southeast1-docker.pkg.dev/cpg-common/images/bazam:v2"
    # cromwell doesn't work with artifact registry:
    # java.lang.Exception: Registry australia-southeast1-docker.pkg.dev is not supported
    docker: "gcr.io/cpg-common/bwa-bazam:v1"
    preemptible: preemptible_tries
    memory: "64 GiB"
    cpu: total_cpu
    disks: "local-disk " + 300 + " HDD"
  }
  output {
    File output_file = output_file
    File output_indx = output_indx
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
    File duplicate_metrics = "~{duplicate_metrics_fname}"
  }
}

task BwaFromBamOrCram {
  input {
    File bam_or_cram
    File bai_or_crai
    String sample_name
    String output_bam_basename
    String duplicate_metrics_fname

    # reference_fasta.ref_alt is the .alt file from bwa-kit
    # (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative".
    ReferenceFasta reference_fasta

    Int preemptible_tries
    Boolean to_cram = false
    String? subset_region 
  }
  
  String output_format = if to_cram then "cram" else "bam"
  String bazam_regions = if defined(subset_region) then "--regions ~{subset_region} " else ""
  
  Int bwa_cpu = 20
  Int bazam_cpu = 5
  Int bamsormadup_cpu = 6
  Int total_cpu = bwa_cpu + bazam_cpu + bamsormadup_cpu

  String rg_line = "@RG\\tID:~{sample_name}\\tSM:~{sample_name}"
  
  String output_file = if to_cram then "~{output_bam_basename}.cram" else "~{output_bam_basename}.bam"
  String output_indx = if to_cram then "~{output_bam_basename}.crai" else "~{output_bam_basename}.bai"
  
  # BWA command options:
  # -K     process INT input bases in each batch regardless of nThreads (for reproducibility)
  # -p     smart pairing (ignoring in2.fq)
  # -v3    minimum score to output [30]
  # -t16   threads
  # -Y     use soft clipping for supplementary alignments
  # -R     read group header line such as '@RG\tID:foo\tSM:bar'
  # -M     mark shorter split hits as secondary
  command <<<
    set -o pipefail
    set -ex

    (while true; do df -h; pwd; du -sh *; free -m; sleep 300; done) &
    
    bazam -Xmx16g -Dsamjdk.reference_fasta=~{reference_fasta.ref_fasta} \
      ~{bazam_regions} -n~{bazam_cpu} -bam ~{bam_or_cram} | \
    bwa mem -K 100000000 -p -v3 -t~{bwa_cpu} -Y -R '~{rg_line}' \
      ~{reference_fasta.ref_fasta} /dev/stdin - \
      2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
    bamsormadup inputformat=sam threads=~{bamsormadup_cpu} SO=coordinate \
      M=~{duplicate_metrics_fname} \
      outputformat=sam | \
    samtools view -T ~{reference_fasta.ref_fasta} \
      -O ~{output_format} \
      -o ~{output_file}
    
    samtools index -@~{total_cpu} ~{output_file} ~{output_indx}

    df -h; pwd; du -sh *
  >>>
  
  runtime {
    # docker: "australia-southeast1-docker.pkg.dev/fewgenomes/images/bazam:v2"
    # cromwell doesn't work with artifact registry:
    # java.lang.Exception: Registry australia-southeast1-docker.pkg.dev is not supported
    docker: "gcr.io/fewgenomes/bazam:v2"
    preemptible: preemptible_tries
    memory: "64 GiB"
    cpu: total_cpu
    disks: "local-disk " + 300 + " HDD"
  }
  output {
    File output_file = output_file
    File output_indx = output_indx
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
    File duplicate_metrics = "~{duplicate_metrics_fname}"
  }
}

task SamSplitter {
  input {
    File input_bam
    Int n_reads
    Int preemptible_tries
    Int compression_level
  }

  Float unmapped_bam_size = size(input_bam, "GiB")
  # Since the output bams are less compressed than the input bam we need a disk multiplier that's larger than 2.
  Float disk_multiplier = 2.5
  Int disk_size = ceil(disk_multiplier * unmapped_bam_size + 20)

  command {
    set -e
    mkdir output_dir

    total_reads=$(samtools view -c ~{input_bam})

    java -Dsamjdk.compression_level=~{compression_level} -Xms3000m -jar /usr/gitc/picard.jar SplitSamByNumberOfReads \
      INPUT=~{input_bam} \
      OUTPUT=output_dir \
      SPLIT_TO_N_READS=~{n_reads} \
      TOTAL_READS_IN_INPUT=$total_reads
  }
  output {
    Array[File] split_bams = glob("output_dir/*.bam")
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    preemptible: preemptible_tries
    memory: "3.75 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
}
