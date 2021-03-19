#!/usr/bin/env nextflow

/*
  Copyright (C) 2021,  OICR
  
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

nextflow.enable.dsl = 2
name = 'open-access-variant-filtering'
short_name = 'open-filter'
version = '0.1.0'  // package version

// universal params go here, change default value as needed
params.container = ""
params.container_registry = ""
params.container_version = ""
params.cpus = 1
params.mem = 1  // GB
params.publish_dir = ""  // set to empty string will disable publishDir

// tool specific parmas go here, add / change as needed
params.study_id = ""
params.analysis_id = ""
params.regions_file = "NO_FILE_regions"
params.output_type = ""
params.apply_filters = [
    'CaVEMan': "PASS", 
    'Pindel': "PASS", 
    'GATK-Mutect2': "PASS"
]
params.include = [
    'CaVEMan': "INFO/CLPM=0 && INFO/ASRD>=0.93", 
    'Pindel': "", 
    'GATK-Mutect2': ""
]
params.exclude = [
    'CaVEMan': "", 
    'Pindel': "", 
    'GATK-Mutect2': ""
]
params.open = true

// if provided local files will be used
params.analysis_metadata = "NO_FILE_metadata"
params.vcf_file = "NO_FILE_vcf"

params.api_token = ""
params.song_url = ""
params.score_url = ""
params.cleanup = true

params.download = [:]
params.upload = [:]
params.filter = [:]
params.payloadGenVcf = [:]

download_params = [
    'song_cpus': params.cpus,
    'song_mem': params.mem,
    'score_cpus': params.cpus,
    'score_mem': params.mem,
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    'publish_dir': params.publish_dir,
    *:(params.download ?: [:])
]

upload_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    'publish_dir': params.publish_dir,
    *:(params.upload ?: [:])
]

filter_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'publish_dir': '',
    'regions_file': params.regions_file,
    'apply_filters': params.apply_filters,
    'output_type': params.output_type,
    'include': params.include,
    'exclude': params.exclude,
    *:(params.filter ?: [:])
]

payloadGenVcf_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'publish_dir': params.publish_dir,
    *:(params.payloadGenVcf ?: [:])
]

include { songScoreDownload as dnVcf } from './song-score-utils/song-score-download' params(download_params)
include { metadataParser as mParser } from "./modules/raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/metadata-parser.0.2.0.0/tools/metadata-parser/metadata-parser.nf"
include { variantFilter as vFilter } from './wfpr_modules/github.com/icgc-argo/variant-calling-tools/variant-filter@0.1.0/main.nf' params(filter_params)
include { payloadGenVariantProcessing as pGenVar } from "./wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/payload-gen-variant-processing@0.1.0/main.nf" params(payloadGenVcf_params)
include { songScoreUpload as upVcf } from './song-score-utils/song-score-upload' params(upload_params)
include { getSecondaryFiles as getSec } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/helper-functions@1.0.0/main.nf'
include { cleanupWorkdir as cleanup } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/cleanup-workdir@1.0.0/main.nf'


// please update workflow code as needed
workflow OpenFilterWf {
  take:  // update as needed
    study_id
    analysis_id
    vcf_file
    analysis_metadata
    regions_file


  main:  // update as needed
    local_mode = false
    vcf = Channel.from()
    vcf_idx = Channel.from()

    if (analysis_id) {
      // download vcf and metadata from song/score (analysis type: variant_calling)
      dnVcf(study_id, analysis_id)
      vcf = dnVcf.out.files.flatten().first()
      vcf_idx = dnVcf.out.files.flatten().last()
      vcf_meta = dnVcf.out.song_analysis
    } else if (
      !analysis_metadata.startsWith('NO_FILE') && \
      !vcf_file.startsWith('NO_FILE')
    ) {
      if (!params.publish_dir) {
        exit 1, "When use local inputs, params.publish_dir must be specified."
      } else {
        log.info "Use local inputs, outputs will be in: ${params.publish_dir}"
      }

      local_mode = true
      vcf = file(vcf_file)
      vcf_idx = Channel.fromPath(getSec(vcf_file, ['tbi']))
      vcf_meta = file(analysis_metadata)
    } else {
      exit 1, "To download input VCF files from SONG/SCORE, please provide `params.analysis_id`.\n" +
                "Or please provide `params.analysis_metadata` and `params.vcf_file` to use local files as input."
    }

    // get analysis_tools info
    mParser(vcf_meta)

    mParser.out.analysis_tools.map{params.apply_filters[it]}.set{apply_filters}
    mParser.out.analysis_tools.map{params.include[it]}.set{include}
    mParser.out.analysis_tools.map{params.exclude[it]}.set{exclude}

    // filter variants
    vFilter(vcf, vcf_idx, file(regions_file), apply_filters, include, exclude, params.output_type)

    // genPayloadSNV
    pGenVar(vcf_meta, vFilter.out.filtered_vcf.concat(vFilter.out.filtered_vcf_tbi).collect(), name, short_name, version, params.open)

    // skip upload if in local_mode
    if (!local_mode) {
        // uploadVariant
        upVcf(study_id, pGenVar.out.payload, pGenVar.out.files_to_upload)
    }

    // cleanup, skip cleanup when running in local mode
    if (params.cleanup) {
      if (local_mode) {
        cleanup(
          vFilter.out.filtered_vcf.concat(vFilter.out.filtered_vcf_tbi, pGenVar.out).collect(), 
          true
        )
      } else {
        cleanup(
          dnVcf.out.files.concat(vFilter.out, pGenVar.out).collect(),
          upVcf.out.analysis_id
        )
      }
    }

  emit:
    payload = pGenVar.out.payload
    output_files = pGenVar.out.files_to_upload
}


// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  OpenFilterWf(
    params.study_id,
    params.analysis_id,
    params.vcf_file,
    params.analysis_metadata,
    params.regions_file
  )
}