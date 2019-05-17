# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2018-11-28 13:28:46
# @Last modified by:   jsgounot
# @Last Modified time: 2018-12-17 10:43:35

import sys
from collections import namedtuple
import multiprocessing

from pgpy import PyVCFError
from pgpy import VCFIterator as VI

FKNT = namedtuple("FKwargs", ["contig", "start", "end"])

def fkwargs2nt(fkwarg) :
    # convert fetch kwargs to namedtuple
    # to be hashable

    return FKNT(fkwarg.get("contig", None), fkwarg.get("start", None), fkwarg.get("end", None))

def wrapper_process(vcf, fkwargs, fun, args, kwargs) :
    return fkwargs, fun(vcf.copy(), * args, ** kwargs)

def apply_mp(vcf, fun, fetch_kwargs, * args, ncore=4, ** kwargs) :
    
    fkwargs = fetch_kwargs
    for fkwarg in fkwargs :
        if set(fkwarg) - set(("contig", "start", "end")) :
            raise PyVCFError("fetch_kwargs can only contain contig, start and end keys")

    vcf = VI(vcf) if isinstance(vcf, str) else vcf
    args = [(vcf, fkwarg, fun, args, {** fkwarg, ** kwargs}) for fkwarg in fkwargs]

    if ncore == 1 :
        return {fkwargs2nt(fkwarg) : result for fkwarg, result in 
                (wrapper_process(* arg) for arg in args)}

    try : 
        with multiprocessing.Pool(processes=ncore) as pool:
            return {fkwargs2nt(fkwarg) : result for fkwarg, result in 
            pool.starmap(wrapper_process, args)}
    
    except KeyboardInterrupt :
        print ("Interrupt childs process")
        pool.terminate()
        sys.exit()
    
    except Exception as e :
        print ("Error : An exception was raised in a process : %s" %(e.__class__.__name__ + " : " + str(e)))
        return None

def apply_mp_contig(vcf, fun, * args, ncore=4, ** kwargs) :
    
    vcf = VI(vcf) if isinstance(vcf, str) else vcf
    fetch_kwargs = [{"contig" : contig} for contig in vcf.contigs]
    return {fkwarg.contig : result for fkwarg, result 
            in apply_mp(vcf, fun, fetch_kwargs, * args, ncore=ncore, ** kwargs).items()}