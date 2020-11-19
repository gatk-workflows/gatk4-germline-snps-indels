# gatk4-germline-snps-indels

### Purpose : 
Workflows for [germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932) with GATK4. 

### haplotypecaller-gvcf-gatk :
The haplotypecaller-gvcf-gatk4 workflow runs the GATK4 HaplotypeCaller tool
in GVCF mode on a single sample according to GATK Best Practices. When 
executed the workflow scatters the HaplotypeCaller tool over the input bam sample 
using an interval list file. The output produced by the workflow will be a single GVCF 
file which can then be provided to GenomicsDBImport along with several other 
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

### Software version requirements :
- GATK 4 
- Samtools 1.3.1
- Python 2.7
- Cromwell version support 
  - Successfully tested on v53

### IMPORTANT NOTE :
- The [JointGenotyping](https://github.com/broadinstitute/warp/tree/develop/pipelines/broad/dna_seq/germline/joint_genotyping) workflow takes the GVCF output produced by the haplotypecaller-gvcf-gatk and uses GenomicsDBImport to produce a multi-sample VCF. The JointGenotyping workflow requires GVCFs be listed in a sample map text file, this can be generated using the [generate-sample-map](https://github.com/gatk-workflows/utility-wdls/blob/main/generate-sample-map.wdl) workflow.
- The provided JSON is a ready to use example JSON template of the workflow. It is the user’s responsibility to correctly set the reference and resource input variables using the [GATK Tool and Tutorial Documentations](https://gatk.broadinstitute.org/hc/en-us/categories/360002310591).
- Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
- For help running workflows on the Google Cloud Platform or locally please
view the following tutorial [(How to) Execute Workflows from the gatk-workflows Git Organization](https://gatk.broadinstitute.org/hc/en-us/articles/360035530952).
- Please visit the [User Guide](https://gatk.broadinstitute.org/hc/en-us/categories/360002310591) site for further documentation on our workflows and tools.
- Relevant reference and resources bundles can be accessed in [Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811).

### Contact Us :
- The following material is provided by the Data Science Platforum group at the Broad Institute. Please direct any questions or concerns to one of our forum sites : [GATK](https://gatk.broadinstitute.org/hc/en-us/community/topics) or [Terra](https://support.terra.bio/hc/en-us/community/topics/360000500432).

### LICENSING :
Copyright Broad Institute, 2019 | BSD-3
This script is released under the WDL open source code license (BSD-3) (full license text at https://github.com/openwdl/wdl/blob/master/LICENSE). Note however that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script.
- [GATK](https://software.broadinstitute.org/gatk/download/licensing.php)
- [Samtools](http://www.htslib.org/terms/)

