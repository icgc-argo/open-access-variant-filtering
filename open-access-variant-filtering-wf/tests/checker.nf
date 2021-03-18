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

params.study_id = ""
params.analysis_id = ""
params.regions_file = "NO_FILE_regions"
params.analysis_metadata = "NO_FILE_metadata"
params.vcf_file = "NO_FILE_vcf"

include { OpenFilterWf } from '../main'

workflow {
  main:
    OpenFilterWf(
      params.study_id,
      params.analysis_id,
      params.vcf_file,
      params.analysis_metadata,
      params.regions_file
    )
}
