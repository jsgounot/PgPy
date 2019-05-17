# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2018-11-21 11:12:45
# @Last modified by:   jsgounot
# @Last Modified time: 2018-11-28 11:30:15

from collections import defaultdict

import pandas as pd

from pgpy import VCFIterator as VI
from pgpy import PyVCFError

def maf(vcf, * args, addRef=False, snps=True, ** kwargs) :

    vcf = VI(vcf)
    if snps : vcf.modifier = VI.no_indel_modifier

    data = defaultdict(int)

    # Check that ref is not a sample name
    if "Reference" in vcf.samples and addRef : 
    	raise PyVCFError("Reference already found in sample")

    for contig, position, ref, variants in vcf.fetch( * args, ** kwargs) :

    	if addRef : variants["Reference"] = (ref, )
    	acount = VI.acount(variants, freq=True)
    	data[min(acount.values())] += 1

    vcount = sum(len(alt) for alt in variants)
    
    df = pd.Series(data).to_frame().reset_index()
    df.columns = ["SampleFreq", "SiteCount"]
    df["SampleCount"] = df["SampleFreq"] * vcount
    df["SiteFreq"] = df["SiteCount"] / df["SiteCount"].sum()

    return df[sorted(df.columns)]