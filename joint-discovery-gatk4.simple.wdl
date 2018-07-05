version 1.0
## Copyright Broad Institute, 2018
## Copyright Garvan Institute, 2018
##
## This WDL implements the joint discovery and VQSR filtering portion of the GATK
## Best Practices (June 2016) for germline SNP and Indel discovery in human
## whole-genome sequencing (WGS) and exome sequencing data.
##
## Requirements/expectations :
## - One or more GVCFs produced by HaplotypeCaller in GVCF mode
## - Bare minimum 1 WGS sample or 30 Exome samples. Gene panels are not supported.
##
## Outputs :
## - A VCF file and its index, filtered using variant quality score recalibration
##   (VQSR) with genotypes for all samples present in the input VCF. All sites that
##   are present in the input VCF are retained; filtered sites are annotated as such
##   in the FILTER field.
##
## Workflow has been simplified and will not scale like the original.
## See https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/joint-discovery-gatk4.wdl
## For a workflow intended to scale
## runtime parameters have been set from a measured single GVCF run and will likely not scale.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

workflow JointGenotyping {
  input {
    # each region in this file is scattered and run in parallel
    # EG 1 line for each chromosome + X,Y,M
    File unpadded_intervals_file

    # a prefix to name output files
    String callset_name

    # GenomicsDBImport param
    # The sample map is a tab-delimited text file with sample_name--tab--path_to_sample_vcf
    # per line. Using a sample map saves the tool from having to download the GVCF headers
    # in order to determine the sample names. Sample names in the sample name map file may
    # have non-tab whitespace, but may not begin or end with whitespace.
    File sample_name_map

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File dbsnp_vcf
    File dbsnp_vcf_index

    String gatk_path

    Array[String]+ snp_recalibration_tranche_values
    Array[String]+ snp_recalibration_annotation_values
    Array[String]+ indel_recalibration_tranche_values
    Array[String]+ indel_recalibration_annotation_values

    # These intervals are used for CollectVariantCallingMetrics
    File eval_interval_list

    File hapmap_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf
    File one_thousand_genomes_resource_vcf_index
    File mills_resource_vcf
    File mills_resource_vcf_index
    File? axiomPoly_resource_vcf
    File? axiomPoly_resource_vcf_index

    # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
    # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
    Float excess_het_threshold = 54.69
    Float snp_filter_level
    Float indel_filter_level
  }

  # possible file size issue with cromwell
  # see system.input-read-limits
  Array[String]+ unpadded_intervals = read_lines(unpadded_intervals_file)

  scatter (idx in range(length(unpadded_intervals))) {
    call GenotypeGVCFs {
      input:
        sample_name_map = sample_name_map,
        interval = unpadded_intervals[idx],
        output_vcf_filename = "output.vcf.gz",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbsnp_vcf = dbsnp_vcf,
        gatk_path = gatk_path,
        excess_het_threshold = excess_het_threshold,
        batch_size = 50
    }
  }

  call GatherVcfs {
    input:
      inputs = GenotypeGVCFs.output_vcf,
      output_vcf_filename = callset_name + ".variant_filtered.vcf.gz",
      gatk_path = gatk_path
  }

  call IndelsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf = GatherVcfs.output_vcf,
      sites_only_variant_filtered_vcf_index = GatherVcfs.output_vcf_index,
      recalibration_filename = callset_name + ".indels.recal",
      tranches_filename = callset_name + ".indels.tranches",
      recalibration_tranche_values = indel_recalibration_tranche_values,
      recalibration_annotation_values = indel_recalibration_annotation_values,
      mills_resource_vcf = mills_resource_vcf,
      mills_resource_vcf_index = mills_resource_vcf_index,
      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      gatk_path = gatk_path
  }


  call SNPsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf = GatherVcfs.output_vcf,
      sites_only_variant_filtered_vcf_index = GatherVcfs.output_vcf_index,
      recalibration_filename = callset_name + ".snps.recal",
      tranches_filename = callset_name + ".snps.tranches",
      recalibration_tranche_values = snp_recalibration_tranche_values,
      recalibration_annotation_values = snp_recalibration_annotation_values,
      hapmap_resource_vcf = hapmap_resource_vcf,
      hapmap_resource_vcf_index = hapmap_resource_vcf_index,
      omni_resource_vcf = omni_resource_vcf,
      omni_resource_vcf_index = omni_resource_vcf_index,
      one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
      one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      gatk_path = gatk_path
  }

  call ApplyRecalibration {
    input:
      output_vcf_filename = callset_name + ".vcf.gz",
      input_vcf = GatherVcfs.output_vcf,
      input_vcf_index = GatherVcfs.output_vcf_index,
      indels_recalibration = IndelsVariantRecalibrator.recalibration,
      indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
      indels_tranches = IndelsVariantRecalibrator.tranches,
      snps_recalibration = SNPsVariantRecalibrator.recalibration,
      snps_recalibration_index = SNPsVariantRecalibrator.recalibration_index,
      snps_tranches = SNPsVariantRecalibrator.tranches,
      indel_filter_level = indel_filter_level,
      snp_filter_level = snp_filter_level,
      gatk_path = gatk_path
  }

  call CollectVariantCallingMetrics {
    input:
      input_vcf = ApplyRecalibration.output_vcf,
      input_vcf_index = ApplyRecalibration.output_vcf_index,
      metrics_filename_prefix = callset_name,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      interval_list = eval_interval_list,
      ref_dict = ref_dict,
      gatk_path = gatk_path
  }

  output {
    File output_vcf = ApplyRecalibration.output_vcf
    File output_vcf_index = ApplyRecalibration.output_vcf_index
    File detail_metrics_file = CollectVariantCallingMetrics.detail_metrics_file
    File summary_metrics_file = CollectVariantCallingMetrics.summary_metrics_file
  }
}

# This is unused by might be needed
# Cromwell has file size read limits for its stdlib
task GetLines {
  input {
    File input_file
  }
  command <<<
    wc -l < ~{input_file}
  >>>
  output {
    Int sample_count = read_int(stdout())
  }
}

task GenotypeGVCFs {
  input {
    File sample_name_map
    String interval

    String output_vcf_filename

    String gatk_path
    String java_opt

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File dbsnp_vcf

    Int batch_size
    String mem_size

    Float excess_het_threshold
  }

  command <<<
    set -e

    genomicsdb="$TMPDIR"/genomicsdb
    # The memory setting here is very important and must be several GB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    "~{gatk_path}" --java-options "~{java_opt}" \
    GenomicsDBImport \
    --genomicsdb-workspace-path "$genomicsdb" \
    --batch-size "~{batch_size}" \
    -L "~{interval}" \
    --sample-name-map "~{sample_name_map}" \
    --reader-threads 5 \
    -ip 500

    tmp_vcf="$TMPDIR"/tmp.vcf.gz

    "~{gatk_path}" --java-options "~{java_opt}" \
     GenotypeGVCFs \
     -R "~{ref_fasta}" \
     -O "$tmp_vcf" \
     -D "~{dbsnp_vcf}" \
     -G StandardAnnotation \
     --only-output-calls-starting-in-intervals \
     --use-new-qual-calculator \
     -V gendb://"$genomicsdb" \
     -L "~{interval}"

    "~{gatk_path}" --java-options "~{java_opt}" \
      VariantFiltration \
      --filter-expression "~{'ExcessHet > ' + excess_het_threshold}" \
      --filter-name ExcessHet \
      -O "~{output_vcf_filename}" \
      -V "$tmp_vcf"
  >>>
  runtime {
    memory: mem_size
    cpu: "2"
  }
  output {
    File output_vcf = "~{output_vcf_filename}"
  }
}

task IndelsVariantRecalibrator {
  input {
    String recalibration_filename
    String tranches_filename

    Array[String]+ recalibration_tranche_values
    Array[String]+ recalibration_annotation_values

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

    File mills_resource_vcf
    File mills_resource_vcf_index
    File? axiomPoly_resource_vcf
    File? axiomPoly_resource_vcf_index
    File dbsnp_vcf
    File dbsnp_vcf_index

    String gatk_path
    String java_opt

    String mem_size
  }

  command <<<
    "~{gatk_path}" --java-options "~{java_opt}" \
      VariantRecalibrator \
      -V "~{sites_only_variant_filtered_vcf}" \
      -O "~{recalibration_filename}" \
      --tranches-file "~{tranches_filename}" \
      --trust-all-polymorphic \
      -tranche "~{sep='" -tranche "' recalibration_tranche_values}" \
      -an "~{sep='" -an "' recalibration_annotation_values}" \
      -mode INDEL \
      --max-gaussians 4 \
      -resource mills,known=false,training=true,truth=true,prior=12:"~{mills_resource_vcf}" \
      ~{'-resource "axiomPoly,known=false,training=true,truth=false,prior=10:' + axiomPoly_resource_vcf + '"'} \
      -resource dbsnp,known=true,training=false,truth=false,prior=2:"~{dbsnp_vcf}"
  >>>
  runtime {
    memory: mem_size
    cpu: "2"
  }
  output {
    File recalibration = "~{recalibration_filename}"
    File recalibration_index = "~{recalibration_filename}.idx"
    File tranches = "~{tranches_filename}"
  }
}

task SNPsVariantRecalibrator {
  input {
    String recalibration_filename
    String tranches_filename

    Array[String]+ recalibration_tranche_values
    Array[String]+ recalibration_annotation_values

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

    File hapmap_resource_vcf
    File omni_resource_vcf
    File one_thousand_genomes_resource_vcf
    File dbsnp_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf_index
    File dbsnp_vcf_index

    String gatk_path
    String java_opt

    String mem_size
  }

  command <<<
    "~{gatk_path}" --java-options "~{java_opt}" \
      VariantRecalibrator \
      -V "~{sites_only_variant_filtered_vcf}" \
      -O "~{recalibration_filename}" \
      --tranches-file "~{tranches_filename}" \
      --trust-all-polymorphic \
      -tranche "~{sep='" -tranche "' recalibration_tranche_values}" \
      -an "~{sep='" -an "' recalibration_annotation_values}" \
      -mode SNP \
      --max-gaussians 6 \
      -resource hapmap,known=false,training=true,truth=true,prior=15:"~{hapmap_resource_vcf}" \
      -resource omni,known=false,training=true,truth=true,prior=12:"~{omni_resource_vcf}" \
      -resource 1000G,known=false,training=true,truth=false,prior=10:"~{one_thousand_genomes_resource_vcf}" \
      -resource dbsnp,known=true,training=false,truth=false,prior=7:"~{dbsnp_vcf}"
  >>>
  runtime {
    memory: mem_size
    cpu: "2"
  }
  output {
    File recalibration = "~{recalibration_filename}"
    File recalibration_index = "~{recalibration_filename}.idx"
    File tranches = "~{tranches_filename}"
  }
}

task ApplyRecalibration {
  input {
    String output_vcf_filename
    File input_vcf
    File input_vcf_index
    File indels_recalibration
    File indels_recalibration_index
    File indels_tranches
    File snps_recalibration
    File snps_recalibration_index
    File snps_tranches

    Float indel_filter_level
    Float snp_filter_level

    String gatk_path
    String java_opt

    String mem_size
  }

  command <<<
    set -e

    "~{gatk_path}" --java-options "-Xmx5g -Xms5g" \
      ApplyVQSR \
      -O tmp.indel.recalibrated.vcf \
      -V "~{input_vcf}" \
      --recal-file "~{indels_recalibration}" \
      --tranches-file "~{indels_tranches}" \
      --truth-sensitivity-filter-level "~{indel_filter_level}" \
      --create-output-variant-index true \
      -mode INDEL

    "~{gatk_path}" --java-options "~{java_opt}" \
      ApplyVQSR \
      -O "~{output_vcf_filename}" \
      -V tmp.indel.recalibrated.vcf \
      --recal-file "~{snps_recalibration}" \
      --tranches-file "~{snps_tranches}" \
      --truth-sensitivity-filter-level "~{snp_filter_level}" \
      --create-output-variant-index true \
      -mode SNP
  >>>
  runtime {
    memory: mem_size
    cpu: "1"
  }
  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
}

task GatherVcfs {
  input {
    Array[File]+ inputs
    String output_vcf_filename
    String gatk_path
    String java_opt

    String mem_size
  }

  command <<<
    set -e

    # ignoreSafetyChecks make a big performance difference so we include it in our invocation
    "~{gatk_path}" --java-options "~{java_opt}" \
    GatherVcfsCloud \
    --input "~{sep='" --input "' inputs}" \
    --output "~{output_vcf_filename}"

    "~{gatk_path}" --java-options "~{java_opt}" \
    IndexFeatureFile \
    --feature-file "~{output_vcf_filename}"
  >>>
  runtime {
    memory: mem_size
    cpu: "1"
  }
  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
}

task CollectVariantCallingMetrics {
  input {
    File input_vcf
    File input_vcf_index

    String metrics_filename_prefix
    File dbsnp_vcf
    File dbsnp_vcf_index
    File interval_list
    File ref_dict

    String gatk_path
    String java_opt

    String mem_size
  }

  command <<<
    "~{gatk_path}" --java-options "~{java_opt}" \
      CollectVariantCallingMetrics \
      --INPUT "~{input_vcf}" \
      --DBSNP "~{dbsnp_vcf}" \
      --SEQUENCE_DICTIONARY "~{ref_dict}" \
      --OUTPUT "~{metrics_filename_prefix}" \
      --THREAD_COUNT 8 \
      --TARGET_INTERVALS "~{interval_list}"
  >>>
  output {
    File detail_metrics_file = "~{metrics_filename_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "~{metrics_filename_prefix}.variant_calling_summary_metrics"
  }
  runtime {
    memory: mem_size
    cpu: 2
  }
}
