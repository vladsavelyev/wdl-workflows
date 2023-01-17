version 1.0

## Copyright Broad Institute, 2018
##
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human whole-genome data.
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in a BAM or CRAM format
## - GVCF output names must end in ".g.vcf.gz"
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "tasks/Alignment.wdl"
import "tasks/AggregatedBamQC.wdl" as AggregatedQC
import "tasks/Qc.wdl" as QC
import "tasks/VariantCalling.wdl"
import "tasks/BamProcessing.wdl" as Processing
import "tasks/structs/DNASeqStructs.wdl"

# WORKFLOW DEFINITION
workflow SingleSample {

  String pipeline_version = "2.3.3"

  input {
    Input inp
    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    PapiSettings papi_settings
    
    Boolean realign = true
    Boolean to_cram = true
    Boolean check_contamination = true
    Boolean check_fingerprints = false

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    File? wgs_coverage_interval_list
    File? target_interval_list
    File? bait_interval_list

    Boolean output_alignment_file = true
    Boolean use_gatk3_haplotype_caller = false
    Boolean validate_gvcf = true

    String? subset_region
  }

  # Not overridable:
  Int read_length = 250
  String cross_check_fingerprints_by = "READGROUP"
  String final_gvcf_base_name = select_first([inp.final_gvcf_base_name, inp.base_file_name])

  if (realign) {
    if (defined(target_interval_list)) {
      call Processing.GenerateSubsettedContaminationResources {
        input:
          target_interval_list = select_first([target_interval_list]),
          contamination_sites_bed = references.contamination_sites_bed,
          contamination_sites_mu = references.contamination_sites_mu,
          contamination_sites_ud = references.contamination_sites_ud,
          preemptible_tries = papi_settings.preemptible_tries
      }
    }

    call Alignment.Alignment {
      input:
        inp = inp,
        references = references,
        papi_settings = papi_settings,
  
        check_contamination = check_contamination,
        check_fingerprints = check_fingerprints,
    
        contamination_sites_ud = select_first([GenerateSubsettedContaminationResources.subsetted_contamination_ud, references.contamination_sites_ud]),
        contamination_sites_bed = select_first([GenerateSubsettedContaminationResources.subsetted_contamination_bed, references.contamination_sites_bed]),
        contamination_sites_mu = select_first([GenerateSubsettedContaminationResources.subsetted_contamination_mu, references.contamination_sites_mu]),
  
        cross_check_fingerprints_by = cross_check_fingerprints_by,
        haplotype_database_file = references.haplotype_database_file,
        lod_threshold = -20.0,

        to_cram = to_cram,
        subset_region = subset_region
    }
  }

  File mapped_file = select_first([Alignment.output_file, inp.bam_or_cram_or_fastq1])
  File mapped_indx = select_first([Alignment.output_indx, inp.bai_or_crai_or_fastq2])

  call AggregatedQC.AggregatedBamQC {
    input:
      base_recalibrated_bam = mapped_file,
      base_recalibrated_bam_index = mapped_indx,
      base_name = inp.base_file_name,
      sample_name = inp.sample_name,
      haplotype_database_file = references.haplotype_database_file,
      references = references,
      fingerprint_genotypes_file = fingerprint_genotypes_file,
      fingerprint_genotypes_index = fingerprint_genotypes_index,
      papi_settings = papi_settings
  }

  if (defined(wgs_coverage_interval_list)) {
    # QC the sample WGS metrics (stringent thresholds)
    call QC.CollectWgsMetrics as CollectWgsMetrics {
      input:
        input_bam = mapped_file,
        input_bam_index = mapped_indx,
        metrics_filename = inp.base_file_name + ".wgs_metrics",
        ref_fasta = references.reference_fasta.ref_fasta,
        ref_fasta_index = references.reference_fasta.ref_fasta_index,
        wgs_coverage_interval_list = select_first([wgs_coverage_interval_list]),
        read_length = read_length,
        preemptible_tries = papi_settings.agg_preemptible_tries
    }
  
    # QC the sample raw WGS metrics (common thresholds)
    call QC.CollectRawWgsMetrics as CollectRawWgsMetrics {
      input:
        input_bam = mapped_file,
        input_bam_index = mapped_indx,
        metrics_filename = inp.base_file_name + ".raw_wgs_metrics",
        ref_fasta = references.reference_fasta.ref_fasta,
        ref_fasta_index = references.reference_fasta.ref_fasta_index,
        wgs_coverage_interval_list = select_first([wgs_coverage_interval_list]),
        read_length = read_length,
        preemptible_tries = papi_settings.agg_preemptible_tries
    }
  }

  call VariantCalling.VariantCalling {
    input:
      calling_interval_list = references.calling_interval_list,
      evaluation_interval_list = references.evaluation_interval_list,
      haplotype_scatter_count = scatter_settings.haplotype_scatter_count,
      break_bands_at_multiples_of = scatter_settings.break_bands_at_multiples_of,
      contamination = Alignment.contamination,
      input_bam = mapped_file,
      input_bam_index = mapped_indx,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      ref_dict = references.reference_fasta.ref_dict,
      dbsnp_vcf = references.dbsnp_vcf,
      dbsnp_vcf_index = references.dbsnp_vcf_index,
      base_file_name = inp.base_file_name,
      final_vcf_base_name = final_gvcf_base_name,
      agg_preemptible_tries = papi_settings.agg_preemptible_tries,
      use_gatk3_haplotype_caller = use_gatk3_haplotype_caller,
      validate_gvcf = validate_gvcf
  }
  
  if (defined(target_interval_list)) {
    call QC.CollectHsMetrics as CollectHsMetrics {
      input:
        input_bam = mapped_file,
        input_bam_index = mapped_indx,
        metrics_filename = inp.base_file_name + ".hybrid_selection_metrics",
        ref_fasta = references.reference_fasta.ref_fasta,
        ref_fasta_index = references.reference_fasta.ref_fasta_index,
        target_interval_list = select_first([target_interval_list]),
        bait_interval_list = select_first([bait_interval_list]),
        preemptible_tries = papi_settings.agg_preemptible_tries
    }
  }

  if (output_alignment_file) {
    File alignment_file_to_output = mapped_file
    File alignment_indx_to_output = mapped_indx
  }

  # Outputs that will be retained when execution is complete
  output {
    File? cross_check_fingerprints_metrics = Alignment.cross_check_fingerprints_metrics

    File? selfSM = Alignment.selfSM
    Float? contamination = Alignment.contamination

    File calculate_read_group_checksum_md5 = AggregatedBamQC.calculate_read_group_checksum_md5

    File agg_alignment_summary_metrics = AggregatedBamQC.agg_alignment_summary_metrics
    File agg_bait_bias_detail_metrics = AggregatedBamQC.agg_bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = AggregatedBamQC.agg_bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = AggregatedBamQC.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = AggregatedBamQC.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = AggregatedBamQC.agg_gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = AggregatedBamQC.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = AggregatedBamQC.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = AggregatedBamQC.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = AggregatedBamQC.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = AggregatedBamQC.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = AggregatedBamQC.agg_quality_distribution_metrics
    File agg_error_summary_metrics = AggregatedBamQC.agg_error_summary_metrics

    File? fingerprint_summary_metrics = AggregatedBamQC.fingerprint_summary_metrics
    File? fingerprint_detail_metrics = AggregatedBamQC.fingerprint_detail_metrics

    File? wgs_metrics = CollectWgsMetrics.metrics
    File? raw_wgs_metrics = CollectRawWgsMetrics.metrics
    File? hybrid_selection_metrics = CollectHsMetrics.metrics

    File? duplicate_metrics = Alignment.duplicate_metrics

    File? gvcf_summary_metrics = VariantCalling.vcf_summary_metrics
    File? gvcf_detail_metrics = VariantCalling.vcf_detail_metrics

    File? alignment_file = alignment_file_to_output
    File? alignment_indx = alignment_indx_to_output

    File output_vcf = VariantCalling.output_vcf
    File output_vcf_index = VariantCalling.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}
