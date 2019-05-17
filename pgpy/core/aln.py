# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2018-11-21 10:39:38
# @Last modified by:   jsgounot
# @Last Modified time: 2019-05-10 10:51:54

from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet.IUPAC import unambiguous_dna

from pgpy import VCFIterator as VI

DEBUG = False

class AlnException(Exception) :

    pass

def aln_maker(vcf, reference, contig, start=None, stop=None, vcf_modifier=None, add_ref=False, check=True, alphabet=None,
    gatkwc2ref=True, wosamples=[], indels=False) :

    if isinstance(vcf, str) :
        vcf = VI(vcf, vcf_modifier)
    
    if vcf_modifier :
        vcf.modifier = vcf_modifier

    with open(reference) as f :
        fdata = SeqIO.parse(f, 'fasta')
        fdata = {res.id : res.seq for res in fdata}

    samples = vcf.samples
    if add_ref and "Reference" in samples : raise AlnException("Cannot add reference : A sample has already this name")
    if add_ref : samples.append("Reference")

    if contig not in fdata :
        raise AlnException("Contig not found in the reference sequence : %s" %(contig))

    # sequence
    start = start or 0
    sequence = str(fdata[contig][start:stop])
    data = {sample : [MutableSeq(sequence) for i in range(ploidy)] for sample, ploidy in vcf.ploidies.items()}
    if add_ref : data["Reference"] = [sequence]

    if DEBUG : print (sequence)

    if indels : _infer_indels(data, vcf, sequence, contig, start, stop, check)
    else : _infer_snps(data, vcf, sequence, contig, start, stop, gatkwc2ref, check)

    alphabet = alphabet or unambiguous_dna
    data = {"%s_%i" %(sample, idx) : Seq(str(sequence), alphabet) for sample, sequences in data.items()
            for idx, sequence in enumerate(sequences, start=1) if sample not in wosamples}

    stop = stop or len(sequence)
    desc = "%s - %i - %i" %(contig, start, stop)

    return MultipleSeqAlignment(SeqRecord(data[sample], id=sample, name=sample, description=desc)
        for sample in sorted(data))

def _infer_snps(data, vcf, sequence, contig, start, stop, gatkwc2ref, check) :

    for contig, position, ref, variants in vcf.fetch(contig=contig, start=start, stop=stop) :

        nposition = position - start - 1
        if nposition < 0 : continue # indels outside bounds

        if check and ref[0] != sequence[nposition] :

            if DEBUG : print (contig, position, nposition, ref, variants)
            raise AlnException("Not the same nucleotide found : Contig : %s - Position : %i - Fasta : %s - VCF : %s" %(
                contig, position + 1, sequence[nposition], ref[0]))

        for sample, alts in variants.items() :

            for idx, alt in enumerate(alts) :
                if not alt : continue
                alt = alt[0]

                if alt == "*" and gatkwc2ref : continue # indels info from GATK
                if alt != ref : data[sample][idx][nposition] = alt

def _infer_indels(data, vcf, sequence, contig, start, stop, check) :

    def lst_get(lst, index, default) :
        # similar to dict get but for list
        try : return lst[index] or default
        except IndexError : return default

    def add_gap(string, size) :
        return string + (size - len(string)) * "-"

    offset = 0

    for contig, position, ref, variants in vcf.fetch(contig=contig, start=start, stop=stop) :

        nposition = position - start - 1 + offset

        if nposition < 0 :
            # Occurs when you have a deletion which start before
            # the given starting point (and overlapp your starting point)

            buff = abs(nposition)
            nposition = 0

            ref = ref[buff:]
            reflen = len(ref)

            variants = {sample : [add_gap(alt[buff:], reflen) for alt in alts if alt]
                for sample, alts in variants.items()}

        subrefseq = sequence[nposition : nposition + len(ref)]

        # print (ref, subrefseq, position, nposition, offset)
        # print (sequence)
        # print (variants)

        if check and ref != subrefseq :
            if DEBUG : print (contig, position, nposition, ref, variants)
            raise AlnException("Not the same nucleotide found : Contig : %s - Position : %i - Fasta : %s - VCF : %s" %(
                contig, position + 1, subrefseq, ref))

        maxlen = max(len(alt) for sample, alts in variants.items() for alt in alts if alt)
        if len(ref) > maxlen : maxlen = len(ref)

        # print (position, nposition, len(ref), maxlen)

        for sample, mseqs in data.items() :
            for idx in range(len(mseqs)) :
                
                alt = lst_get(variants.get(sample, []), idx, ref)
                if alt == "*" : alt = ref # indels info from GATK
                alt = add_gap(alt, maxlen)

                if maxlen == 1 and alt != ref : data[sample][idx][nposition] = alt 
                elif maxlen > 1 : data[sample][idx] = data[sample][idx][:nposition] + alt + data[sample][idx][nposition+len(ref):]

        refalt = add_gap(ref, maxlen)
        sequence = sequence[:nposition] + refalt + sequence[nposition+len(ref):]

        offset += maxlen - len(ref)

def snps_aln(* args, ** kwargs) :
    # for compatibility purpose
    return aln_maker(* args, indels=False, ** kwargs)