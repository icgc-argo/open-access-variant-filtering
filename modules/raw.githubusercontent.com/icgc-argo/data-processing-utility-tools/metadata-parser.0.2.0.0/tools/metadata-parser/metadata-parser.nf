#!/usr/bin/env nextflow

/*
 * Copyright (c) 2019, Ontario Institute for Cancer Research (OICR).
 *                                                                                                               
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/*
 * Authors:
 *   Junjun Zhang <junjun.zhang@oicr.on.ca>
 *   Linda Xiang <linda.xiang@oicr.on.ca>
 */

nextflow.enable.dsl=2
version = '0.2.0.0'

params.metadata_analysis = ""
params.container_version = ""
params.cpus = 1
params.mem = 1  // GB


process metadataParser {
  container "quay.io/icgc-argo/metadata-parser:metadata-parser.${params.container_version ?: version}"
  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path metadata_analysis

  output:
    env STUDY_ID, emit: study_id
    env DONOR_ID, emit: donor_id
    env EXP, emit: experimental_strategy
    env PAIRED, emit: paired
    env ANALYSIS_TOOLS, emit: analysis_tools

  script:
    """
    set -euxo pipefail
    STUDY_ID=`cat ${metadata_analysis} | jq -er '.studyId' | tr -d '\\n'`
    DONOR_ID=`cat ${metadata_analysis} | jq -er '.samples[0].donor.donorId' | tr -d '\\n'`
    EXP=`cat ${metadata_analysis} | jq -er '.experiment | .experimental_strategy?  // .library_strategy' | tr -d '\\n'`
    VARIABLE1=`cat ${metadata_analysis} | jq -r 'if ([.read_groups[]?] | length) >0 then [.read_groups[] | .is_paired_end] | all | tostring else null end' | tr -d '\\n'`
    PAIRED=\${VARIABLE1:-'null'}
    VARIABLE2=`cat ${metadata_analysis} | jq -r '[.files[] | .info? | .analysis_tools[]?] | unique | join(",")' | tr -d '\\n'`
    ANALYSIS_TOOLS=\${VARIABLE2:-'null'}
    """
}
