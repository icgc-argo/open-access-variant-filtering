#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
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
"""

import os
import argparse
import json
import uuid
import hashlib
import copy

variant_type_to_data_type_etc = {
    'snv': ['Simple Nucleotide Variation', 'Raw SNV Calls', ['CaVEMan', 'bcftools'], ['GATK-Mutect2', 'bcftools']],   # dataCategory, dataType, analysis_tools
    'indel': ['Simple Nucleotide Variation', 'Raw InDel Calls', ['Pindel', 'bcftools'], ['GATK-Mutect2', 'bcftools']]
}

workflow_full_name = {
    'open-access-variant-filtering': 'Open Access Variant Filtering'
}

def calculate_size(file_path):
    return os.stat(file_path).st_size


def calculate_md5(file_path):
    md5 = hashlib.md5()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            md5.update(chunk)
    return md5.hexdigest()


def get_files_info(file_to_upload, args):
    fname = os.path.basename(file_to_upload).split(".")
    input_wf = fname[5]
    variant_type = fname[7]
    process_name = args.wf_short_name
    new_fname = ".".join(fname[0:8]+[process_name, 'vcf.gz']+(['tbi'] if file_to_upload.endswith('.tbi') else []))
    file_info = {
        'fileName': new_fname,
        'fileType': 'VCF' if new_fname.endswith('.vcf.gz') else new_fname.split(".")[-1].upper(),
        'fileSize': calculate_size(file_to_upload),
        'fileMd5sum': calculate_md5(file_to_upload),
        'fileAccess': 'open' if args.open else 'controlled',
        'info': {
            'data_category': variant_type_to_data_type_etc[variant_type][0]
        }
    }

    if file_to_upload.endswith('.vcf.gz'):
        file_info['dataType'] = variant_type_to_data_type_etc[variant_type][1]
    elif file_to_upload.endswith('.vcf.gz.tbi'):
        file_info['dataType'] = 'VCF Index'
    else:
        pass

    if input_wf in (['sanger-wgs', 'sanger-wxs']):
        file_info['info']['analysis_tools'] = variant_type_to_data_type_etc[variant_type][2]
    elif input_wf in (['gatk-mutect2']):
        file_info['info']['analysis_tools'] = variant_type_to_data_type_etc[variant_type][3]

    new_dir = 'out'
    try:
        os.mkdir(new_dir)
    except FileExistsError:
        pass

    dst = os.path.join(os.getcwd(), new_dir, new_fname)
    os.symlink(os.path.abspath(file_to_upload), dst)

    return file_info

def get_sample_info(sample_list):
    samples = copy.deepcopy(sample_list)
    for sample in samples:
        for item in ['info', 'sampleId', 'specimenId', 'donorId', 'studyId']:
            sample.pop(item, None)
            sample['specimen'].pop(item, None)
            sample['donor'].pop(item, None)

    return samples

def main():
    """
    Python implementation of tool: payload-gen-variant-processing
    """

    parser = argparse.ArgumentParser(description='Tool: payload-gen-variant-processing')
    parser.add_argument("-a", dest="analysis", required=True,
                        help="json file containing sequencing_alignment SONG analysis for tumour sample")
    parser.add_argument("-f", dest="files_to_upload", type=str, nargs="+", help="Files to be uploaded", required=True)
    parser.add_argument("-w", dest="wf_name", type=str, help="workflow full name", required=True)
    parser.add_argument("-s", dest="wf_short_name", type=str, help="workflow short name", required=True)
    parser.add_argument("-v", dest="wf_version", type=str, required=True, help="workflow version")
    parser.add_argument("-r", dest="wf_run", type=str, required=True, help="workflow run ID")
    parser.add_argument("-j", dest="wf_session", type=str, required=True, help="workflow session ID")
    parser.add_argument("-o", dest="open", action='store_true', help="set file to be open access")
    args = parser.parse_args()

    analysis = {}
    with open(args.analysis, 'r') as f:
        analysis = json.load(f)

    input_analysis_type = analysis.get('analysisType').get('name')
    analysis_type = "variant_processing"
    payload = {
        'analysisType': {
            'name': analysis_type
        },
        'studyId': analysis.get('studyId'),  
        'experiment': analysis.get('experiment'),
        'samples': get_sample_info(analysis.get('samples')),
        'files': [],
        'workflow': {
            'workflow_name': workflow_full_name.get(args.wf_name),
            'workflow_short_name': args.wf_short_name,
            'workflow_version': args.wf_version,
            'run_id': args.wf_run,
            'session_id': args.wf_session,
            'inputs': [
                {
                    'input_analysis_id': analysis.get('analysisId'),
                    'analysis_type': input_analysis_type
                }
            ],
            'genome_build': 'GRCh38_hla_decoy_ebv'
        },
        'variant_class': analysis.get('variant_class')
    }

    for f in args.files_to_upload:
        file_info = get_files_info(f, args)
        payload['files'].append(file_info)

    with open("%s.%s.payload.json" % (str(uuid.uuid4()), analysis_type), 'w') as f:
        f.write(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

