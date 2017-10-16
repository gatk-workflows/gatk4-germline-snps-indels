# gatk4-germline-snps-indels

Workflows for germline short variant discovery with GATK4

Purpose : 
Gatk Germiline SNPs Indels is composed of two WDLs, haplotypeCaller-gvcf-gatk4 and
joint-discovery-gatk4.The haplotypecaller-gvcf-gatk4 workflow runs HaplotypeCaller 
from GATK4 in GVCF mode on a single sample according to the GATK Best Practices (June 2016), 
scattered across intervals.The second WDL implements the joint discovery and VQSR 
filtering portion of the GATK Best Practices (June 2016) for germline SNP and Indel 
discovery in human whole-genome sequencing (WGS) and exome sequencing data.

haplotypecaller-gvcf-gatk :
Requirements/expectations
- One analysis-ready BAM file for a single sample (as identified in RG:SM)
- Set of variant calling intervals lists for the scatter, provided in a file
Outputs 
- One GVCF file and its index

joint-discovery-gatk :
Requirements/expectations
- One or more GVCFs produced by HaplotypeCaller in GVCF mode
- Bare minimum 1 WGS sample or 30 Exome samples. Gene panels are not supported.
Outputs 
- A VCF file and its index, filtered using variant quality score recalibration  
  (VQSR) with genotypes for all samples present in the input VCF. All sites that  
  are present in the input VCF are retained; filtered sites are annotated as such  
  in the FILTER field.

Cromwell version support : 
- Successfully tested on v29
- Does not work on versions < v23 due to output syntax

IMPORTANT NOTE : 
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
