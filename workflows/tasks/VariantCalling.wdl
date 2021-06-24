version 1.0

import "Qc.wdl" as QC
import "Utilities.wdl" as Utils
import "BamProcessing.wdl" as BamProcessing

workflow VariantCalling {

  String pipeline_version = "1.0.0"

  input {
    File calling_interval_list
    File evaluation_interval_list
    Int haplotype_scatter_count
    Int break_bands_at_multiples_of
    Float? contamination
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? dbsnp_vcf
    File? dbsnp_vcf_index
    String base_file_name
    String final_vcf_base_name
    Int agg_preemptible_tries
    Boolean make_gvcf = true
    Boolean make_bamout = false
    Boolean use_gatk3_haplotype_caller = false
    Boolean validate_gvcf = true
  }

  parameter_meta {
    make_bamout: "For CNNScoreVariants to run with a 2D model, a bamout must be created by HaplotypeCaller. The bamout is a bam containing information on how HaplotypeCaller remapped reads while it was calling variants. See https://gatkforums.broadinstitute.org/gatk/discussion/5484/howto-generate-a-bamout-file-showing-how-haplotypecaller-has-remapped-sequence-reads for more details."
  }

  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  call Utils.ScatterIntervalList as ScatterIntervalList {
    input:
      interval_list = calling_interval_list,
      scatter_count = haplotype_scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of
  }

  # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
  # If we take the number we are scattering by and reduce by 20 we will have enough disk space
  # to account for the fact that the data is quite uneven across the shards.
  Int potential_hc_divisor = ScatterIntervalList.interval_count - 20
  Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

  # Call variants in parallel over WGS calling intervals
  scatter (scattered_interval_list in ScatterIntervalList.out) {

    if (use_gatk3_haplotype_caller) {
      call HaplotypeCaller_GATK35_GVCF as HaplotypeCallerGATK3 {
        input:
          input_bam = input_bam,
          input_bam_index = input_bam_index,
          interval_list = scattered_interval_list,
          gvcf_basename = base_file_name,
          ref_dict = ref_dict,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          contamination = contamination,
          preemptible_tries = agg_preemptible_tries,
          hc_scatter = hc_divisor
      }
    }

    if (!use_gatk3_haplotype_caller) {

      # Generate GVCF by interval
      call HaplotypeCaller_GATK4_VCF as HaplotypeCallerGATK4 {
        input:
          contamination = contamination,
          input_bam = input_bam,
          input_bam_index = input_bam_index,
          interval_list = scattered_interval_list,
          vcf_basename = base_file_name,
          ref_dict = ref_dict,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          hc_scatter = hc_divisor,
          make_gvcf = make_gvcf,
          make_bamout = make_bamout,
          preemptible_tries = agg_preemptible_tries
       }

      # If bamout files were created, we need to sort and gather them into one bamout
      if (make_bamout) {
        call BamProcessing.SortSam as SortBamout {
          input:
            input_bam = select_first([HaplotypeCallerGATK4.bamout]),
            output_bam_basename = final_vcf_base_name,
            preemptible_tries = agg_preemptible_tries,
            compression_level = 2
        }
      }
    }

    File vcfs_to_merge = select_first([HaplotypeCallerGATK3.output_gvcf, HaplotypeCallerGATK4.output_vcf])
    File vcf_indices_to_merge = select_first([HaplotypeCallerGATK3.output_gvcf_index, HaplotypeCallerGATK4.output_vcf_index])
  }

  # Combine by-interval (g)VCFs into a single sample (g)VCF file
  String merge_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
  call MergeVCFs {
    input:
      input_vcfs = vcfs_to_merge,
      input_vcfs_indexes = vcf_indices_to_merge,
      output_vcf_name = final_vcf_base_name + merge_suffix,
      preemptible_tries = agg_preemptible_tries
  }

  if (make_bamout) {
    call MergeBamouts {
      input:
        bams = select_all(SortBamout.output_bam),
        output_base_name = final_vcf_base_name
    }
  }

  # Validate the (g)VCF output of HaplotypeCaller
  if (validate_gvcf && defined(dbsnp_vcf) && defined(dbsnp_vcf_index)) {
    call QC.ValidateVCF as ValidateVCF {
      input:
        input_vcf = MergeVCFs.output_vcf,
        input_vcf_index = MergeVCFs.output_vcf_index,
        dbsnp_vcf = select_first([dbsnp_vcf]),
        dbsnp_vcf_index = select_first([dbsnp_vcf_index]),
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        calling_interval_list = calling_interval_list,
        is_gvcf = make_gvcf,
        preemptible_tries = agg_preemptible_tries
    }

    # QC the (g)VCF
    call QC.CollectVariantCallingMetrics as CollectVariantCallingMetrics {
      input:
        input_vcf = MergeVCFs.output_vcf,
        input_vcf_index = MergeVCFs.output_vcf_index,
        metrics_basename = final_vcf_base_name,
        dbsnp_vcf =  select_first([dbsnp_vcf]),
        dbsnp_vcf_index =  select_first([dbsnp_vcf_index]),
        ref_dict = ref_dict,
        evaluation_interval_list = evaluation_interval_list,
        is_gvcf = make_gvcf,
        preemptible_tries = agg_preemptible_tries
    }
  }

  output {
    File? vcf_summary_metrics = CollectVariantCallingMetrics.summary_metrics
    File? vcf_detail_metrics = CollectVariantCallingMetrics.detail_metrics
    File output_vcf = MergeVCFs.output_vcf
    File output_vcf_index = MergeVCFs.output_vcf_index
    File? bamout = MergeBamouts.output_bam
    File? bamout_index = MergeBamouts.output_bam_index
  }
  meta {
    allowNestedInputs: true
  }
}

task HaplotypeCaller_GATK35_GVCF {
  input {
    File input_bam
    File input_bam_index
    File interval_list
    String gvcf_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? contamination
    Int preemptible_tries
    Int hc_scatter
  }

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(((size(input_bam, "GiB") + 30) / hc_scatter) + ref_size) + 20

  # We use interval_padding 500 below to make sure that the HaplotypeCaller has context on both sides around
  # the interval because the assembly uses them.
  #
  # Using PrintReads is a temporary solution until we update HaploypeCaller to use GATK4. Once that is done,
  # HaplotypeCaller can stream the required intervals directly from the cloud.
  command {
    /usr/gitc/gatk4/gatk --java-options "-Xms2g" \
      PrintReads \
      -I ~{input_bam} \
      --interval-padding 500 \
      -L ~{interval_list} \
      -O local.sharded.bam \
    && \
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms8000m \
      -jar /usr/gitc/GATK35.jar \
      -T HaplotypeCaller \
      -R ~{ref_fasta} \
      -o ~{gvcf_basename}.vcf.gz \
      -I local.sharded.bam \
      -L ~{interval_list} \
      -ERC GVCF \
      --max_alternate_alleles 3 \
      -variant_index_parameter 128000 \
      -variant_index_type LINEAR \
      -contamination ~{default=0 contamination} \
      --read_filter OverclippedRead
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    preemptible: preemptible_tries
    memory: "10 GiB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_gvcf = "~{gvcf_basename}.vcf.gz"
    File output_gvcf_index = "~{gvcf_basename}.vcf.gz.tbi"
  }
}

task HaplotypeCaller_GATK4_VCF {
  input {
    File input_bam
    File input_bam_index
    File interval_list
    String vcf_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? contamination
    Boolean make_gvcf
    Boolean make_bamout
    Int preemptible_tries
    Int hc_scatter
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
  }

  String output_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
  String output_file_name = vcf_basename + output_suffix

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(((size(input_bam, "GiB") + 30) / hc_scatter) + ref_size) + 20

  String bamout_arg = if make_bamout then "-bamout ~{vcf_basename}.bamout.bam" else ""

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command <<<
    set -e
    gatk --java-options "-Xms7000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
      HaplotypeCaller \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -L ~{interval_list} \
      -O ~{output_file_name} \
      -contamination ~{default=0 contamination} \
      -G StandardAnnotation -G StandardHCAnnotation ~{true="-G AS_StandardAnnotation" false="" make_gvcf} \
      -GQB 20 \
      ~{true="-ERC GVCF" false="" make_gvcf} \
      ~{bamout_arg}

    # Cromwell doesn't like optional task outputs, so we have to touch this file.
    touch ~{vcf_basename}.bamout.bam
  >>>

  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "8 GiB"
    cpu: "2"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_vcf = "~{output_file_name}"
    File output_vcf_index = "~{output_file_name}.tbi"
    File? bamout = "~{vcf_basename}.bamout.bam"
  }
}
  
# This task is here because merging bamout files using Picard produces an error.
task MergeBamouts {

  input {
    Array[File] bams
    String output_base_name
  }

  Int disk_size = ceil(size(bams, "GiB") * 2) + 10

  command {
    samtools merge ~{output_base_name}.bam ~{sep=" " bams}
    samtools index ~{output_base_name}.bam
    mv ~{output_base_name}.bam.bai ~{output_base_name}.bai
  }

  output {
    File output_bam = "~{output_base_name}.bam"
    File output_bam_index = "~{output_base_name}.bai"
  }

  runtime {
    docker: "biocontainers/samtools:1.3.1"
    memory: "4 GiB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: 3
    cpu: 1
  }
}
  
# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_name
    Int preemptible_tries
  }

  Int disk_size = ceil(size(input_vcfs, "GiB") * 2.5) + 10

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    java -Xms2000m -jar /usr/picard/picard.jar \
      MergeVcfs \
      INPUT=~{sep=' INPUT=' input_vcfs} \
      OUTPUT=~{output_vcf_name}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "3 GiB"
    disks: "local-disk ~{disk_size} HDD"
  }
  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}

task HardFilterVcf {
  input {
    File input_vcf
    File input_vcf_index
    String vcf_basename
    File interval_list
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.4.0"
  }

  Int disk_size = ceil(2 * size(input_vcf, "GiB")) + 20
  String output_vcf_name = vcf_basename + ".filtered.vcf.gz"

  command {
     gatk --java-options "-Xms3000m" \
      VariantFiltration \
      -V ~{input_vcf} \
      -L ~{interval_list} \
      --filter-expression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0" \
      --filter-name "HardFiltered" \
      -O ~{output_vcf_name}
  }
  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "3 GiB"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
  }
}
