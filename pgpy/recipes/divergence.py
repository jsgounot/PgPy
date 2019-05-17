# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2018-11-25 19:38:50
# @Last modified by:   jsgounot
# @Last Modified time: 2018-11-28 15:29:39

from itertools import combinations

import pandas as pd

from pgpy import VCFIterator as VI
from pgpy import PyVCFError
from pgpy.core.processing import apply_mp_contig

def ref_div_modifier(alts, site, snps=True, modifier=None) :

    if modifier :
        alts = modifier(alts, site)

    if snps :
        return tuple(alt[0] if alt != None else site.ref for alt in alts if alt)
    else :
        return tuple(alt if alt != None else site.ref for alt in alts if alt)

def ref_strict_div(variants, ref, data, haploid=False) :

    if not haploid : 
        variants = {sample : tuple(set(alts)) for sample, alts in variants.items()}

    for sample, alt in variants.items() :

        # we have to check cause sometime we can have
        # indels but with no snps

        if alt != (ref, ) :
            data[sample] += 1

def ref_clean_diff(variants, ref, data) :

    raise Exception("To do")

    pass

def ref_divergence(vcf, * args, snps=True, strict=True, ** kwargs) :

    # Faster than pairwise divergence

    vcf = VI(vcf) if isinstance(vcf, str) else vcf
    previous = vcf.modifier
    vcf.modifier = lambda * args : ref_div_modifier(* args, snps, previous)

    data = {sample : 0 for sample in vcf.samples}
    fun = ref_strict_div if strict else ref_clean_diff
    haploid = None

    for contig, position, ref, variants in vcf.fetch( * args, ** kwargs) :
        
        if haploid is None :
            haploid = all(len(alts) == 1 for alts in variants.values())
            fun = ref_strict_div if haploid else fun

        fun(variants, ref, data, haploid)

    data = ({"sample" : sample, "diff" : count} for sample, count in data.items())
    return pd.DataFrame(data).sort_values(by="sample")

def ref_divergence_wg(vcf, * args, ncore=4, ** kwargs) :

    results = apply_mp_contig(vcf, ref_divergence, * args, ncore=ncore, ** kwargs)

    mpres = []
    for contig, result in results.items() :
        result = result.set_index("sample")["diff"]
        result.name = contig
        mpres.append(result)

    df = pd.concat(mpres, axis=1)
    df["Total"] = df.sum(axis=1)
    return df[sorted(df.columns)]

# -------------------
# Pairwise divergence
# -------------------

def pw_div_modifier(alts, site, snps=True, modifier=None) :

    if modifier :
        alts = modifier(alts, site)

    if snps :
        return tuple(alt[0] if alt != None else site.ref for alt in alts)
    else :
        return tuple(alt if alt != None else site.ref for alt in alts)

def pairwise_strict_div(variants, combinations, data, haploid=False) :

    if not haploid : 
        variants = {sample : tuple(set(alts)) for sample, alts in variants.items()}

    for pair in combinations :
        s1, s2 = pair
        if variants[s1] != variants[s2] :
            data[pair] += 1

def pairwise_clean_diff(variants, combinations, data) :

    raise Exception("To do")

    pass

def pairwise_divergence(vcf, * args, snps=True, strict=True, addRef=True, ** kwargs) :

    vcf = VI(vcf) if isinstance(vcf, str) else vcf
    previous = vcf.modifier
    vcf.modifier = lambda * args : pw_div_modifier(* args, snps, previous)

    poss = [pair for pair in combinations(vcf.samples, 2)]
    data = {pair : 0 for pair in poss}

    fun = pairwise_strict_div if strict else pairwise_clean_diff
    haploid = None

    for contig, position, ref, variants in vcf.fetch( * args, ** kwargs) :

        if addRef :
            variants["Ref"] = (ref, )

        if haploid == None :
            haploid = all(len(alts) == 1 for alts in variants.values())
            fun = pairwise_strict_div if haploid else fun

        fun(variants, poss, data, haploid)

    data = ({"sample1" : pair[0], "sample2" : pair[1], "diff" : count} for pair, count in data.items())
    return pd.DataFrame(data).sort_values(by=["sample1", "sample2"])