#!/usr/bin/env nextflow

/*
  Copyright (C) 2021,  Ontario Institute for Cancer Research
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.
  
  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Authors:
    Linda Xiang
*/

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '0.1.0'  // package version

container = [
    'ghcr.io': 'ghcr.io/icgc-argo/variant-calling-tools.variant-filter'
]
default_container_registry = 'ghcr.io'
/********************************************************************/


// universal params go here
params.container_registry = ""
params.container_version = ""
params.container = ""

params.cpus = 1
params.mem = 1  // GB
params.publish_dir = ""  // set to empty string will disable publishDir


// tool specific params go here, add / change as needed
params.input_file = ""
params.regions_file = "NO_FILE_regions"
params.output_type = ""
params.apply_filters = ""
params.include = ""
params.exclude = ""

include { getSecondaryFiles as getSec } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/helper-functions@1.0.0/main.nf'

process variantFilter {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:  // input, make update as needed
    path input_file
    path input_file_idx
    path regions_file
    val apply_filters
    val include
    val exclude
    val output_type

  output:  // output, make update as needed
    path "*.filtered.vcf.gz", emit: filtered_vcf
    path "*.filtered.vcf.gz.tbi", emit: filtered_vcf_tbi

  script:
    // add and initialize variables here as needed
    arg_regions_file = regions_file.name == 'NO_FILE_regions' ? "" : " --regions-file ${regions_file}"
    arg_apply_filters = apply_filters == '' ? "" : " --apply-filters '${apply_filters}'"
    arg_include = include == '' ? "" : " --include '${include}'"
    arg_exclude = exclude == '' ? "" : " --exclude '${exclude}'"
    arg_output_type = output_type == '' ? "" : " --output-type ${output_type}"

    """
    main.py \
      -v ${input_file} \
      ${arg_output_type} \
      ${arg_regions_file} \
      ${arg_apply_filters} \
      ${arg_include} \
      ${arg_exclude} 

    """
}


// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  variantFilter(
    file(params.input_file),
    Channel.fromPath(getSec(params.input_file, ['tbi']), checkIfExists: true).collect(),
    file(params.regions_file),
    params.apply_filters,
    params.include,
    params.exclude,
    params.output_type
  )
}
