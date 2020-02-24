# gatk4-germline-snps-indels

### Purpose : 
Workflows for [germline short variant discovery](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145) with GATK4. 

### haplotypecaller-gvcf-gatk :
The haplotypecaller-gvcf-gatk4 workflow runs the GATK4 HaplotypeCaller tool
in GVCF mode on a single sample according to GATK Best Practices. When 
executed the workflow scatters the HaplotypeCaller tool over the input bam sample 
using an interval list file. The output produced by the workflow will be a single GVCF 
file which can then be provided to the JointGenotyping workflow along with several other 
GVCF files to call for variants simultaneously, producing a multisample VCF. 
The haplotypecaller-gvcf-gatk4 workflows default GVCF mode is useful when calling variants 
for several samples efficiently. However, for instances when calling variants for one or a 
few samples it is possible to have the workflow directly call variants and output a VCF file by 
setting the `make_gvcf` input variable to `false`. 

#### Requirements/expectations
- One analysis-ready BAM file for a single sample (as identified in RG:SM)
- A file containing a set of variant calling interval list for the scatter

#### Outputs 
- One GVCF file and its index

### JointGenotyping.wdl :
This WDL implements the joint calling and VQSR filtering portion of the 
GATK Best Practices for germline SNP and Indel discovery 
in human whole-genome sequencing (WGS). The workflow requires a sample map 
file with 50 or more GVCFs and produces a multisample VCF.

*NOTE:*  
*- JointGenotyping-terra.wdl is a slightly modified version of the 
original workflow to support users interested in running the 
workflow on Terra. The changes include variables for dockers and disksize, making 
it easier to configure the workflow.*
*- Creating a sample map can be nuisance on Terra, use the [generate-sample-map](https://portal.firecloud.org/?return=terra#methods/gatk/generate-sample-map/1) to create one for you.*


#### Requirements/expectations
- One or more GVCFs produced by HaplotypeCaller in GVCF mode
- Bare minimum 50 samples. Gene panels are not supported.
- When determining disk size in the JSON, use the guideline below
  - small_disk = (num_gvcfs / 10) + 10
  - medium_disk = (num_gvcfs * 15) + 10
  - huge_disk = num_gvcfs + 10

### Outputs 
- A VCF file and its index, filtered using variant quality score recalibration  
  (VQSR) with genotypes for all samples present in the input VCF. All sites that  
  are present in the input VCF are retained; filtered sites are annotated as such  
  in the FILTER field.

### Software version requirements :
- GATK 4.1.4.0 
- Samtools 1.3.1
- Python 2.7
- Cromwell version support 
  - Successfully tested on v37
  - Does not work on versions < v23 due to output syntax

### IMPORTANT NOTE : 
- VQSR wiring. The SNP and INDEL models are built in parallel, but then the corresponding 
  recalibrations are applied in series. Because the INDEL model is generally ready 
  first (because there are fewer indels than SNPs) we set INDEL recalibration to 
  be applied first to the input VCF, while the SNP model is still being built. By 
  the time the SNP model is available, the indel-recalibrated file is available to 
  serve as input to apply the SNP recalibration. If we did it the other way around, 
  we would have to wait until the SNP recal file was available despite the INDEL 
  recal file being there already, then apply SNP recalibration, then apply INDEL 
  recalibration. This would lead to a longer wall clock time for complete workflow 
  execution. Wiring the INDEL recalibration to be applied first solves the problem.
- The current version of the posted "Generic germline short variant joint genotyping" 
  is derived from the Broad production version of the workflow, which was adapted for 
  large WGS callsets of up to 20K samples.  We believe the results of this workflow run 
  on a single WGS sample are equally accurate, but there may be some shortcomings when 
  the workflow is modified and run on small cohorts.  Specifically, modifying the SNP 
  ApplyRecalibration step for higher specificity may not be effective.  The user can verify 
  if this is an issue by consulting the gathered SNP tranches file.  If the listed 
  truthSensitivity in the rightmost column is not well matched to the targetTruthSensitivity 
  in the leftmost column, then requesting that targetTruthSensitivity from ApplyVQSR will 
  not use an accurate filtering threshold.  This workflow has not been tested on exomes.  
  The dynamic scatter interval creating was optimized for genomes.  The scattered SNP 
  VariantRecalibration may fail because of two few "bad" variants to build the negative model. 
  Also, apologies that the logging for SNP recalibration is overly verbose.
- No allele subsetting for the JointGenotyping workflow
  - for large cohorts, even exome callsets can have more than 1000 alleles at low 
    complexity/STR sites
  - for sites with more than 6 alternate alleles (by default) called genotypes will be returned, 
    but without the PLs since the PL arrays get enormous
  - allele-specific filtering could be performed if AS annotations are present, 
    but the data will still be in the VCF in one giant INFO field
- JointGenotyping output is divided into lots of shards
  - desirable for use in [Hail](https://hail.is/), which supports parallel import
  - Its possible to use [GatherVcfs](https://gatk.broadinstitute.org/hc/en-us/search?utf8=%E2%9C%93&query=GatherVcfs) to combine shards.
- Users working with large sample sets can invoke the GnarlyGenotyper task in the JointGenotyping.wdl workflow. However, the [ReblockGVCF](https://gatk.broadinstitute.org/hc/en-us/articles/360037593171-ReblockGVCF-BETA-) tool must be run for all GVCFs produced by HaplotypeCaller before they can be appropriately processed by GnarlyGenotyper. A workflow that applies the reblocking tool is provided here: [ReblockGVCF-gatk4_exomes_goodCompression](https://portal.firecloud.org/?return=terra#methods/methodsDev/ReblockGVCF-gatk4_exomes_goodCompression/4)
- GnarlyGenotyper uses a QUAL score approximation
  - dramatically improves performance compared with GenotypeGVCFs, but QUAL output (and thus
    the QD annotation) may be slightly discordant between the two tools
- The provided JSON is meant to be a ready to use example JSON template of the workflow. It is the user’s responsibility to correctly set the reference and resource input variables using the [GATK Tool and Tutorial Documentations](https://software.broadinstitute.org/gatk/documentation/).
- Relevant reference and resources bundles can be accessed in [Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle).
- Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
- For help running workflows on the Google Cloud Platform or locally please
view the following tutorial [(How to) Execute Workflows from the gatk-workflows Git Organization](https://software.broadinstitute.org/gatk/documentation/article?id=12521).
- The following material is provided by the GATK Team. Please post any questions or concerns to one of our forum sites : [GATK](https://gatkforums.broadinstitute.org/gatk/categories/ask-the-team/) , [FireCloud](https://gatkforums.broadinstitute.org/firecloud/categories/ask-the-firecloud-team) or [Terra](https://broadinstitute.zendesk.com/hc/en-us/community/topics/360000500432-General-Discussion) , [WDL/Cromwell](https://gatkforums.broadinstitute.org/wdl/categories/ask-the-wdl-team).
- Please visit the [User Guide](https://software.broadinstitute.org/gatk/documentation/) site for further documentation on our workflows and tools.

### LICENSING :
Copyright Broad Institute, 2019 | BSD-3
This script is released under the WDL open source code license (BSD-3) (full license text at https://github.com/openwdl/wdl/blob/master/LICENSE). Note however that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script.
- [GATK](https://software.broadinstitute.org/gatk/download/licensing.php)
- [Samtools](http://www.htslib.org/terms/)

