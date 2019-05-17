# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2018-11-28 17:06:41
# @Last modified by:   jsgounot
# @Last Modified time: 2018-11-28 17:07:41

import operator
from functools import reduce

from pgpy import VCFIterator as VI
from pgpy.core.processing import apply_mp_contig
from pgpy.core.aln import snps_aln

def aln_wholegenome(vcf, fasta, ncore=4) :

	aln = apply_mp_contig(vcf, snps_aln, fasta, ncore=ncore)
	aln = reduce(operator.add, aln.values())

	for record in aln :
		record.description = ""

	return aln