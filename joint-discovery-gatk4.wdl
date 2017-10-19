## Copyright Broad Institute, 2017
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
## Note about VQSR wiring :
## The SNP and INDEL models are built in parallel, but then the corresponding 
## recalibrations are applied in series. Because the INDEL model is generally ready 
## first (because there are fewer indels than SNPs) we set INDEL recalibration to 
## be applied first to the input VCF, while the SNP model is still being built. By 
## the time the SNP model is available, the indel-recalibrated file is available to 
## serve as input to apply the SNP recalibration. If we did it the other way around, 
## we would have to wait until the SNP recal file was available despite the INDEL 
## recal file being there already, then apply SNP recalibration, then apply INDEL 
## recalibration. This would lead to a longer wall clock time for complete workflow 
## execution. Wiring the INDEL recalibration to be applied first solves the problem.
##
## Cromwell version support 
## - Successfully tested on v29
## - Does not work on versions < v23 due to output syntax
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

workflow JointDiscovery_GATK4 {

    File input_gvcfs_list
    File input_gvcfs_indices_list
    File sample_name_map
    String cohort_vcf_name

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File dbsnp_vcf

    Array[String] SNP_annotations
    Array[String] INDEL_annotations
    Array[Float] SNP_tranches
    Array[Float] INDEL_tranches
    Array[String] SNP_resources
    Array[String] INDEL_resources
    Array[File] resource_files
    Array[File] resource_indices

    Float SNP_filter_level
    Float INDEL_filter_level

    File calling_intervals_list

    String gatk_docker
    String picard_docker

    String gatk_launch_path
    String picard_path

    Array[File] input_gvcfs = read_lines(input_gvcfs_list)
    Array[File] input_gvcf_indices = read_lines(input_gvcfs_indices_list)
  
    Array[String] calling_intervals = read_lines(calling_intervals_list)

    # Joint-call variants in parallel over WGS calling intervals
    scatter (idx in range(length(calling_intervals))) {

        call ImportGVCFs {
            input:
                sample_name_map = sample_name_map,
                interval = calling_intervals[idx],
                workspace_dir_name = "genomicsdb",
                docker_image = gatk_docker,
                gatk_launch_path = gatk_launch_path
        }

        # Perform joint genotyping per interval
        call GenotypeGVCFs {
            input:
                workspace_tar = ImportGVCFs.output_genomicsdb,
                vcf_basename = cohort_vcf_name,
                interval = calling_intervals[idx],
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                dbsnp_vcf = dbsnp_vcf,
                docker_image = gatk_docker,
                gatk_launch_path = gatk_launch_path
        }
    }

    # Merge per-interval VCFs into a single cohort VCF file
    call MergeVCFs {
        input:
            input_vcfs = GenotypeGVCFs.genotyped_vcf,
            vcf_name = cohort_vcf_name + ".vcf.gz",
            vcf_index = cohort_vcf_name + ".vcf.gz.tbi",
            docker_image = picard_docker,
            picard_path = picard_path       
    }

    # Build SNP model 
    call BuildVQSRModel as BuildVQSRModelForSNPs {
        input:
            cohort_vcf = MergeVCFs.output_vcf,
            cohort_vcf_index = MergeVCFs.output_vcf_index,
            interval_list = calling_intervals_list,
            output_basename = cohort_vcf_name,
            annotations = SNP_annotations,
            mode = "SNP",
            tranches = SNP_tranches,
            resources = SNP_resources,
            resource_files = resource_files,
            resource_indices = resource_indices,
            docker_image = gatk_docker,
            gatk_launch_path = gatk_launch_path
    }

    # Build INDEL model 
    call BuildVQSRModel as BuildVQSRModelForINDELs {
        input:
            cohort_vcf = MergeVCFs.output_vcf,
            cohort_vcf_index = MergeVCFs.output_vcf_index,
            interval_list = calling_intervals_list,
            output_basename = cohort_vcf_name,
            annotations = INDEL_annotations,
            mode = "INDEL",
            tranches = INDEL_tranches,
            resources = INDEL_resources,
            resource_files = resource_files,
            resource_indices = resource_indices,
            docker_image = gatk_docker,
            gatk_launch_path = gatk_launch_path
    }

    # Apply INDEL filter (first because INDEL model is usually done sooner)
    call ApplyRecalibrationFilter as ApplyRecalibrationFilterForINDELs {
        input:
            cohort_vcf = MergeVCFs.output_vcf,
            cohort_vcf_index = MergeVCFs.output_vcf_index,
            interval_list = calling_intervals_list,
            output_basename = cohort_vcf_name + ".recal.INDEL",
            mode = "INDEL",
            recal_file = BuildVQSRModelForINDELs.recal_file,
            recal_file_index = BuildVQSRModelForINDELs.recal_file_index,
            tranches_file = BuildVQSRModelForINDELs.tranches_file,
            filter_level = INDEL_filter_level,
            docker_image = gatk_docker,
            gatk_launch_path = gatk_launch_path
    }

    # Apply SNP filter
    call ApplyRecalibrationFilter as ApplyRecalibrationFilterForSNPs {
        input:
            cohort_vcf = ApplyRecalibrationFilterForINDELs.recalibrated_vcf,
            cohort_vcf_index = ApplyRecalibrationFilterForINDELs.recalibrated_vcf_index,
            interval_list = calling_intervals_list,
            output_basename = cohort_vcf_name + ".recal.INDEL.SNP",
            mode = "SNP",
            recal_file = BuildVQSRModelForSNPs.recal_file,
            recal_file_index = BuildVQSRModelForSNPs.recal_file_index,
            tranches_file = BuildVQSRModelForSNPs.tranches_file,
            filter_level = SNP_filter_level,
            docker_image = gatk_docker,
            gatk_launch_path = gatk_launch_path
    }

    # Outputs that will be retained when execution is complete
    output {
        File jointcalled_vcf = MergeVCFs.output_vcf
        File jointcalled_vcf_index = MergeVCFs.output_vcf_index
        File filtered_vcf = ApplyRecalibrationFilterForSNPs.recalibrated_vcf
        File filtered_vcf_idx = ApplyRecalibrationFilterForSNPs.recalibrated_vcf_index
    }
}

# TASK DEFINITIONS

# ...
task ImportGVCFs {
  File sample_name_map
  String interval

  String workspace_dir_name

  Int disk_size
  String mem_size

  String docker_image
  String gatk_launch_path
  String java_opt

  command <<<

    rm -rf ${workspace_dir_name}

    # The memory setting here is very important and must be several GB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    ${gatk_launch_path}gatk-launch --javaOptions "${java_opt}" \
      GenomicsDBImport \
      --genomicsDBWorkspace ${workspace_dir_name} \
      --sampleNameMap ${sample_name_map} \
      --readerThreads 5 \
      -L ${interval} \
      -ip 500

    tar -cf ${workspace_dir_name}.tar ${workspace_dir_name}

  >>>
  runtime {
    docker: docker_image
    memory: mem_size
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_genomicsdb = "${workspace_dir_name}.tar"
  }
}

# Perform joint-genotyping
task GenotypeGVCFs { 

    String vcf_basename
    File workspace_tar
  	String interval

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    String dbsnp_vcf

    Int disk_size
    String mem_size

    String docker_image
    String gatk_launch_path
    String java_opt

    String workspace_name = basename(workspace_tar, ".tar")

	command <<<

        tar -xf ${workspace_tar}

		${gatk_launch_path}gatk-launch --javaOptions "${java_opt}" \
        	GenotypeGVCFs \
        	-R ${ref_fasta} \
            -V gendb://${workspace_name} \
        	-L ${interval} \
            -D ${dbsnp_vcf} \
        	-O ${vcf_basename}.vcf.gz \
            -G StandardAnnotation \
            -newQual \
            --onlyOutputCallsStartingInIntervals
	>>>

	output {
		File genotyped_vcf = "${vcf_basename}.vcf.gz"
        File genotyped_vcf_index = "${vcf_basename}.vcf.gz.tbi"
	}

	runtime {
		docker: docker_image
		memory: mem_size
    	cpu: "2"
    	disks: "local-disk " + disk_size + " HDD"
	}
}

# Combine multiple VCFs 
task MergeVCFs {
    Array [File] input_vcfs
    String vcf_name
    String vcf_index

    Int disk_size
    String mem_size

    String docker_image
    String picard_path
    String java_opt

    command {
	  java ${java_opt} -jar ${picard_path}picard.jar \
	    MergeVcfs \
	    INPUT=${sep=' INPUT=' input_vcfs} \
	    OUTPUT=${vcf_name}
    }

  	runtime {
	    docker: docker_image
	    memory: mem_size
	    disks: "local-disk " + disk_size + " HDD"
	}

    output {
    	File output_vcf = "${vcf_name}"
    	File output_vcf_index = "${vcf_index}"
    }
}

# Build VQSR model
task BuildVQSRModel {
    File cohort_vcf
    File cohort_vcf_index
    String output_basename
    File interval_list
    String mode
    Int maxGaussians
    Array[String] annotations
    Array[Float] tranches
    Array[String] resources
    Array[File] resource_files
    Array[File] resource_indices

    Int disk_size
    String mem_size

    String docker_image
    String gatk_launch_path
    String java_opt

    command {
        ${gatk_launch_path}gatk-launch --javaOptions "${java_opt}" \
            VariantRecalibrator \
            -variant ${cohort_vcf} \
            -L ${interval_list} \
            -resource ${sep=' -resource ' resources} \
            -an ${sep=' -an ' annotations} \
            -mode ${mode} \
            -tranche ${sep=' -tranche ' tranches} \
            -O ${output_basename}.${mode}.recal \
            --maxGaussians ${maxGaussians} \
            -tranchesFile ${output_basename}.${mode}.tranches \
            -rscriptFile ${output_basename}.${mode}.plots.R
    }

    runtime {
        docker: docker_image
        memory: mem_size
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File recal_file = "${output_basename}.${mode}.recal"
        File recal_file_index = "${output_basename}.${mode}.recal.idx"
        File tranches_file = "${output_basename}.${mode}.tranches"
        File rscript_file = "${output_basename}.${mode}.plots.R"
    }
}

# Apply recalibration
task ApplyRecalibrationFilter {
    File cohort_vcf
    File cohort_vcf_index
    File recal_file
    File recal_file_index
    File interval_list
    String output_basename
    String mode
    File tranches_file
    Float filter_level

    Int disk_size
    String mem_size

    String docker_image
    String gatk_launch_path
    String java_opt

    command {
        ${gatk_launch_path}gatk-launch --javaOptions "${java_opt}" \
            ApplyVQSR \
            -V ${cohort_vcf} \
            -L ${interval_list} \
            --ts_filter_level ${filter_level} \
            --createOutputVariantIndex true \
            -mode ${mode} \
            -recalFile ${recal_file} \
            -tranchesFile ${tranches_file} \
            -O ${output_basename}.vcf.gz
    }

    runtime {
        docker: docker_image
        memory: mem_size
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File recalibrated_vcf = "${output_basename}.vcf.gz"
        File recalibrated_vcf_index = "${output_basename}.vcf.gz.tbi"
    }
}

