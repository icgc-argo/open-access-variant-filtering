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
    lindaxiang
*/

/*
 This is an auto-generated checker workflow to test the generated main template workflow, it's
 meant to illustrate how testing works. Please update to suit your own needs.
*/

nextflow.enable.dsl = 2
version = '0.1.0'  // package version

// universal params
params.publish_dir = ""
params.container = ""
params.container_registry = ""
params.container_version = ""

// tool specific parmas go here, add / change as needed
params.input_file = ""
params.expected_output = ""
params.cleanup = false

include { OpenAccessVariantFilteringWf } from '../main'
// include section starts
// include section ends

Channel
  .fromPath(params.input_file, checkIfExists: true)
  .set { input_file }


process file_smart_diff {
  input:
    path output_file
    path expected_file

  output:
    stdout()

  script:
    """
    # Note: this is only for demo purpose, please write your own 'diff' according to your own needs.
    # remove date field before comparison eg, <div id="header_filename">Tue 19 Jan 2021<br/>test_rg_3.bam</div>
    # sed -e 's#"header_filename">.*<br/>test_rg_3.bam#"header_filename"><br/>test_rg_3.bam</div>#'

    diff <( cat ${output_file} | sed -e 's#"header_filename">.*<br/>#"header_filename"><br/>#' ) \
         <( ([[ '${expected_file}' == *.gz ]] && gunzip -c ${expected_file} || cat ${expected_file}) | sed -e 's#"header_filename">.*<br/>#"header_filename"><br/>#' ) \
    && ( echo "Test PASSED" && exit 0 ) || ( echo "Test FAILED, output file mismatch." && exit 1 )
    """
}


workflow checker {
  take:
    input_file
    expected_output

  main:
    OpenAccessVariantFilteringWf(
      input_file
    )

    file_smart_diff(
      OpenAccessVariantFilteringWf.out.output_file,
      expected_output
    )
}


workflow {
  checker(
    file(params.input_file),
    file(params.expected_output)
  )
}
