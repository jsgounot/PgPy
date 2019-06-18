# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-06-18 17:58:16
# @Last modified by:   jsgounot
# @Last Modified time: 2019-06-18 18:15:05

from pgpy import VCFIterator as VI
from collections import defaultdict
import pandas as pd

def snp_modifier_ref(alts, site) :
    return tuple(alt[0] if alt else site.ref[0] for alt in alts)

def homhet_window(vcf, wsize=10000, snp_modifier=True) :
    
    vcf = VI(vcf) if isinstance(vcf, str) else vcf
    data = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(int))))
        
    if snp_modifier :
        vcf.modifier = snp_modifier_ref

    for contig, position, ref, variants in vcf.fetch() :
        window = position // wsize * wsize
            
        for sample, variant in variants.items() :
            lvariants = len(set(variant))
            if lvariants  == 1 : data[sample][contig][window]["hom"] += 1
            elif lvariants > 1 : data[sample][contig][window]["het"] += 1
            else : raise Exception("Category error : %s" %(str(variant)))
    
    data = [{"sample" : sample, "contig" : contig, "window" : window, ** values}
           for sample, contigs in data.items() for contig, windows in contigs.items()
           for window, values in windows.items()]
    
    df = pd.DataFrame(data)
    df["hom"] = df["hom"].fillna(0).astype(int)

    if "het" in df.columns : df["het"] = df["het"].fillna(0).astype(int)
    else : df["het"] = 0

    return df