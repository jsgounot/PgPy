# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2018-11-21 11:12:45
# @Last modified by:   jsgounot
# @Last Modified time: 2019-05-28 11:51:39

from collections import defaultdict

import pandas as pd

from pgpy import VCFIterator as VI
from pgpy import PyVCFError
from pgpy.core.processing import apply_mp_contig

def maf_snp_modifier(alts, site) :
    return tuple(alt[0] if alt else site.ref[0] for alt in alts)

def _maf(vcf, * args, addRef=False, snps=True, ** kwargs) :

    vcf = VI(vcf) if isinstance(vcf, str) else vcf 
    if snps : vcf.modifier = maf_snp_modifier
    data = defaultdict(int)

    # Check that ref is not a sample name
    if "Reference" in vcf.samples and addRef : 
        raise PyVCFError("Reference already found in sample")

    for contig, position, ref, variants in vcf.fetch( * args, ** kwargs) :

        if addRef : variants["Reference"] = (ref, )
        acount = VI.acount(variants, freq=True)
        data[min(acount.values())] += 1

    # sometimes we can get a freq of 1 -> Case of indels and all variants are actually the reference
    data = {freq : count for freq, count in data.items() if freq <= .5}
    vcount = sum(len(alts) for alts in variants.values())
    
    df = pd.Series(data).to_frame().reset_index()
    df.columns = ["SampleFreq", "AlleleCount"]
    df["SampleCount"] = df["SampleFreq"] * vcount
    df["AlleleFreq"] = df["AlleleCount"] * 100 / df["AlleleCount"].sum()

    return df[sorted(df.columns)]

def maf(vcf, * args, ncore=1, ** kwargs) :
    if ncore == 1 : return _maf(vcf, * args, ** kwargs)

    df = apply_mp_contig(vcf, _maf, * args, ncore=ncore, ** kwargs)
    df = pd.concat(list(df.values()))
    df = df.groupby(["SampleCount", "SampleFreq"])["AlleleCount"].sum()
    df = df.rename("AlleleCount").reset_index()
    df["AlleleFreq"] = df["AlleleCount"] * 100 / df["AlleleCount"].sum()
    return df