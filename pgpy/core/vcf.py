# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2018-11-21 10:04:35
# @Last modified by:   jsgounot
# @Last Modified time: 2019-01-15 12:53:37

from collections import Counter
from itertools import chain, permutations

from Bio.Data.IUPACData import ambiguous_dna_values

from pysam import VariantFile

# ---------------------------------------------------

def load_iupac() :
    iupac_code = {tuple(letters) : code for code, letters in ambiguous_dna_values.items()}
    return {permutation : (code, ) for letters, code in iupac_code.items()
    for permutation in permutations(letters)}

# ---------------------------------------------------

class PyVCFError(Exception) :

    pass

class VCFIterator() :

    iupac_code = load_iupac()
    
    SHOWMODE = False
    SHOWPOSI = 100000
    show_last = None

    def __init__(self, fname, modifier=None) :

        if not isinstance(fname, str) : raise PyVCFError("fname must be a string")
        self.fname = fname
        self.vcf = VariantFile(fname)
        self.modifier = modifier

    @property
    def samples(self):
        return list(self.vcf.header.samples)

    @property
    def contigs(self):
        return list(self.vcf.header.contigs)

    @property
    def ploidies(self):
        # pmod = self.modifier
        # self.modifier = None
        ploidies = {sample : len(alts) for sample, alts
                    in self.fetch_first()[3].items()}
        # self.modifier = pmod
        return ploidies    

    def copy(self) :
        vcf = VCFIterator(self.fname)
        vcf.modifier = self.modifier
        return vcf

    def add_modifier(self) :

        self.modifier = modifier

    def search_pos(self, contig, position, alts=False) :

        result = list(self.vcf.fetch(contig=contig, start=position-1, stop=position))
        if not result : raise PyVCFError("Variant not found at : %s %i" %(contig, position))
        return VCFIterator.dict_alts(result[0], self.modifier) if alts else result[0]

    def fetch_raw(self, * args, ** kwargs) :

        for variant in self.vcf.fetch(* args, ** kwargs) :
            yield variant

    def fetch(self, * args, ** kwargs) :

        show_last_posi = None
        show_last_contig = None

        for site in self.vcf.fetch(* args, ** kwargs) :

            if VCFIterator.SHOWMODE :
                if site.contig != show_last_contig :
                    show_last_contig = site.contig
                    show_last_posi = site.pos // VCFIterator.SHOWPOSI
                    print (site.contig, site.pos)

                elif site.pos // VCFIterator.SHOWPOSI != show_last_posi :
                    show_last_posi = site.pos // VCFIterator.SHOWPOSI
                    print (site.contig, site.pos)

            data = {}
            for sample, alleles in site.samples.items() :
                
                alts = alleles.alleles
                if self.modifier : alts = self.modifier(alts, site)
                if alts : data[sample] = alts

            if not data : continue
            yield site.contig, site.pos, site.ref, data

    def fetch_first(self, * args, ** kwargs) :

        return next(self.fetch(* args, ** kwargs))

    # Custom modifier

    @staticmethod

    def fill_values(alts, site) :
        return tuple(alt or site.ref for alt in alts)

    @staticmethod
    def iupac_modifier(alts, site) :

        bref = site.ref[0]
        alts = tuple({alt[0] if alt is not None and alt != "*" else bref for alt in alts}) if len(alts) > 1 else tuple(alts)
        return VCFIterator.iupac_code.get(alts, alts)

    @staticmethod
    def no_indel_modifier(alts, site) :

        print (alts)
        exit()

        return tuple(alt[0] for alt in alts if alt)

    @staticmethod
    def only_variable_modifier(alts, site) :

        return tuple(alt for alt in alts if alt and alt != site.ref)

    @staticmethod
    def only_indel_modifier(alts, site) :

        return tuple(alt for alt in alts if alt and len(alt) != len(site.ref))

    # Usual functions

    @staticmethod
    def dict_alts(variant, modifier) :
        
        data = {}
        
        for sample, alleles in variant.samples.items() :       
            alts = alleles.alleles
            if modifier : alts = modifier(alts)
            data[sample] = alts

        return data

    @staticmethod
    def acount(variants, freq=False) :

        c = Counter(chain(variants.values()))
        if freq : c = {variant : count / sum(c.values()) for variant, count in c.items()}
        return c