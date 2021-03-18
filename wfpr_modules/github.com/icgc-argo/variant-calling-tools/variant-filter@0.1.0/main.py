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
import sys
import argparse
import subprocess
import re

def run_cmd(cmd):
    try:
        p = subprocess.run([cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           shell=True, check=True)

    except subprocess.CalledProcessError as e:   # this is triggered when cmd returned non-zero code
        print(e.stdout.decode("utf-8"))
        print('Execution returned non-zero code: %s. Additional error message: %s' %
              (e.returncode, e.stderr.decode("utf-8")), file=sys.stderr)
        sys.exit(e.returncode)

    except Exception as e:  # all other errors go here, like unable to find the command
        sys.exit('Execution failed: %s' % e)

    return p  # in case the caller of this funtion needs p.stdout, p.stderr etc


def main():
    """
    Python implementation of tool: variant-filter
    """

    parser = argparse.ArgumentParser(description='Tool: variant-filter')
    parser.add_argument('-v', '--input-file', dest='input_file', type=str,
                        help='Input file', required=True)
    parser.add_argument('-R', '--regions-file', dest='regions_file', type=str,
                        help='Restrict to regions listed in a file')
    parser.add_argument('-f', '--apply-filters', dest='apply_filters', type=str,
                        help='Skip sites where FILTER column does not contain any of the strings listed in LIST.')                    
    parser.add_argument('-i', '--include', dest='include', type=str,
                        help='Select sites for which the expression is true')
    parser.add_argument('-e', '--exclude', dest='exclude', type=str,
                        help='Exclude sites for which the expression is true')                    
    parser.add_argument('-O', '--output-type ', dest='output_type', default='z', choices=['b', 'u', 'z', 'v'],
                        help='Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). ')            
    args = parser.parse_args()

    basename = re.sub(r'\.vcf\.gz$', '', os.path.basename(args.input_file))
    output = basename+'.filtered.vcf.gz'

    cmd = f'bcftools view {args.input_file} --output-file {output} --output-type {args.output_type}'

    if args.regions_file:
        cmd = cmd + f" -R {args.regions_file}"
    
    if args.apply_filters:
        cmd = cmd + f" -f '{args.apply_filters}'"

    if args.include:
        cmd = cmd + f" --include '{args.include}'"
    
    if args.exclude:
        cmd = cmd + f" --exclude '{args.exclude}'"

    run_cmd(cmd)

    cmd = f'bcftools index --tbi {output}'

    run_cmd(cmd)


if __name__ == "__main__":
    main()

