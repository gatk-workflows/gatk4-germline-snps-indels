# gatk4-germline-snps-indels

### Purpose : 
Workflows for germline short variant discovery with GATK4. 

### haplotypecaller-gvcf-gatk :
The haplotypecaller-gvcf-gatk4 workflow runs HaplotypeCaller 
from GATK4 in GVCF mode on a single sample according to the GATK Best Practices (June 2016), 
scattered across intervals.

#### Requirements/expectations
- One analysis-ready BAM file for a single sample (as identified in RG:SM)
- Set of variant calling intervals lists for the scatter, provided in a file
#### Outputs 
- One GVCF file and its index

### joint-discovery-gatk :
The second WDL implements the joint discovery and VQSR 
filtering portion of the GATK Best Practices (June 2016) for germline SNP and Indel 
discovery in human whole-genome sequencing (WGS) and exome sequencing data.
*NOTE: joint-discovery-gatk4-fc.wdl is a slightly modified version of the original to support users interested in running the workflow on [FireCloud](https://software.broadinstitute.org/firecloud/).*

#### Requirements/expectations
- One or more GVCFs produced by HaplotypeCaller in GVCF mode
- Bare minimum 1 WGS sample or 30 Exome samples. Gene panels are not supported.
- When deteriming disk size in the json, use the guideline below
  - small_disk = (num_gvcfs / 10) + 10
  - medium_disk = (num_gvcfs * 15) + 10
  - huge_disk = num_gvcfs + 10
#### Outputs 
- A VCF file and its index, filtered using variant quality score recalibration  
  (VQSR) with genotypes for all samples present in the input VCF. All sites that  
  are present in the input VCF are retained; filtered sites are annotated as such  
  in the FILTER field.

### Software version requirements :
- GATK 4.beta.3 or later 
- Picard 2.x
- Samtools (see gotc docker)
- Python 2.7

Cromwell version support 
- Successfully tested on v29
- Does not work on versions < v23 due to output syntax

### IMPORTANT NOTE : 
- Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
- HaplotypeCaller in GATK4 is still in evaluation phase and should not
  be used in production until it has been fully vetted. In the meantime, use the GATK3 
  version for any production needs.
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
