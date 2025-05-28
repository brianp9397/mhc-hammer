#!/usr/bin/env nextflow


//Take netMHCpan-formatted HLAs and format them for mhc_hammer
process convert_alleles_to_mhc_hammer_format {

  label 'lohhla'
  label 'lohhla_container'
  tag "${dataset}/${pat_name}/${run}"

  input:
  tuple val(pat_name), val(run_name), val(dataset), path(hla_input)
  tuple path(hla_ref), path(hla_ref_faidx)

  output:
  tuple val(pat_name), val(run_name), val(dataset), path("*formatted_hlas.txt"), emit: mhc_hammer_alleles

  """
  python /opt/hla_formatting.py -f ${hla_input} -m ${hla_ref} -o ${dataset}-${pat_name}-${run_name}.formatted_hlas.txt
  """
}

//Take Sequenza alternative_solutions file with purity/ploidy and format for mhc_hammer
process mhc_hammer_solutions {

  label 'lohhla'
  label 'lohhla_container'
  tag "${dataset}/${pat_name}/${run}"

  input:
  tuple val(pat_name), val(run), val(dataset), path(norm_bam), path(norm_bai), path(tumor_bam), path(tumor_bai), path(purity_input)

  output:
  tuple val(pat_name), val(run), val(dataset), path("*solutions.txt"), emit: formatted_solutions

  """
  f=\$(basename ${tumor_bam})
#  p=\$(find ${purity_input}/ -name '*alternative_solutions*')

  python /opt/solutions.py -f ${purity_input} -b \$f -o ${dataset}-${pat_name}-${run}.mhc_hammer_solutions.txt
  """
}


process mhc_hammer_inventory {

  label 'lohhla'
  label 'lohhla_container'
  tag "${dataset}/${pat_name}/${run}"

  input:
  tuple val(pat_name), val(run), val(dataset), path(norm_bam), path(norm_bai), path(tumor_bam), path(tumor_bai), path(rna_bam), path(rna_bai), path(formatted_hlas), path(formatted_solutions)

  output:
  tuple val(pat_name), val(run), val(dataset), path(norm_bam), path(norm_bai), path(tumor_bam), path(tumor_bai), path(rna_bam), path(rna_bai), path("*hla_alleles.csv"), path("*inventory.csv"), emit: inventory

  """
  purity=\$(awk -F',' 'NR==2{print \$3}' $formatted_solutions)
  ploidy=\$(awk -F',' 'NR==2{print \$4}' $formatted_solutions)
  norm_name=\$(basename ${norm_bam})
  norm_base=\$(basename ${norm_bam} .bam)
  norm_file=\$(readlink -f ${norm_bam})
  norm_bai_file=\$(readlink -f ${norm_bai})
  tumor_base=\$(basename ${tumor_bam} .bam)
  tumor_name=\$(basename ${tumor_bam})
  tumor_file=\$(readlink -f ${tumor_bam})
  tumor_bai_file=\$(readlink -f ${tumor_bai})
  rna_base=\$(basename ${rna_bam} .bam)
  rna_name=\$(basename ${rna_bam})
  rna_file=\$(readlink -f ${rna_bam})
  rna_bai_file=\$(readlink -f ${rna_bai})
  hla1=\$(awk 'NR==1' $formatted_hlas)
  hla2=\$(awk 'NR==2' $formatted_hlas)
  hla3=\$(awk 'NR==3' $formatted_hlas)
  hla4=\$(awk 'NR==4' $formatted_hlas)
  hla5=\$(awk 'NR==5' $formatted_hlas)
  hla6=\$(awk 'NR==6' $formatted_hlas)

  echo "A,\${hla1},\${hla2}" >> ${dataset}-${pat_name}-${run}.hla_alleles.csv
  echo "B,\${hla3},\${hla4}" >> ${dataset}-${pat_name}-${run}.hla_alleles.csv
  echo "C,\${hla5},\${hla6}" >> ${dataset}-${pat_name}-${run}.hla_alleles.csv

  wd=\$(readlink -f .)

  echo "patient,sample_name,sample_type,bam_path,bai_path,sequencing_type,purity,ploidy,normal_sample_name,hla_alleles_path" >> ${dataset}-${pat_name}-${run}.inventory.csv
  echo "${pat_name},\${tumor_base},tumour,\${tumor_file},\${tumor_bai_file},wxs,\${purity},\${ploidy},\${norm_base},\${wd}/${dataset}-${pat_name}-${run}.hla_alleles.csv" >> ${dataset}-${pat_name}-${run}.inventory.csv
  echo "${pat_name},\${norm_base},normal,\${norm_file},\${norm_bai_file},wxs,\${purity},\${ploidy},,\${wd}/${dataset}-${pat_name}-${run}.hla_alleles.csv" >> ${dataset}-${pat_name}-${run}.inventory.csv
  echo "${pat_name},\${rna_base},tumour,\${rna_file},\${rna_bai_file},rna,\${purity},\${ploidy},,\${wd}/${dataset}-${pat_name}-${run}.hla_alleles.csv" >> ${dataset}-${pat_name}-${run}.inventory.csv

  """
}

