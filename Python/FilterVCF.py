#!/usr/bin/env python

from __future__ import print_function
import sys
import warnings
from pysam import VariantFile


def main(argv):
    bcf_in = VariantFile(argv[0])  # auto-detect input format
    bcf_out = VariantFile(argv[0] + '.filtered.vcf.gz', 'w', header=bcf_in.header)
    for site in bcf_in.fetch():
        keep_site = 0 # default option is to remove SNP
        for sample, rec in site.samples.items():
            if max(rec.get('GP')[1:]) > 0.9:
                keep_site = 1 # do not remove SNP if either het or non-ref homo is greater than .9 for any sample
        if keep_site:
            bcf_out.write(site)
                          

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
           
if __name__ == "__main__":
    global usage
    usage = "python FilterVCF.py VCF_file.vcf.gz (must have corresponding VCF_file.vcf.csi/tbi file)"
    warnings.formatwarning = warning_on_one_line
    main(sys.argv[1:])
