#!/usr/bin/env python

'''
    This program takes a reference DNA sequence (from a FASTA file) and
    returns modified DNA reads containing SNPs for a defined number of individuals,
    under the assumption that they belong to one population group and the distribution
    of SNPs follows the H-W equilibrium.
    
    INPUT:      Fasta file - consensus DNA sequence
    OUTPUT:     Fasta file: haploid DNA sequences peppered with SNPs
    Population SNP table: containing "true" allelic frequencies and
    individuals' genotype
    '''

import os, sys, re, argparse, subprocess
import numpy as np
from scipy.stats import nbinom, beta
import random

def WriteVcfHeader(VcfFile, opts):

    template="""##fileformat=VCFv4.1
##fileData=Veritas
##reference={fastaFile}
##phasing=known
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"""
    
    template = template + "\t".join(["".join(["INDIV_",str(i)]) for i in range(0,opts.nindiv)])+ "\n"
    
    context = {
    "fastaFile": opts.input
    }
    
    VcfFile.write(template.format(**context))


'''
    This generator parses the Fasta file line by
    line and yields one sequence entry at a time.
    
    The FASTA header entries are sequestered and
    will not be used in the simulated Fasta files.
    The original header ID will be replaced by
    a much simple header ID, e.g. ">contig_1".
    '''
def readFasta(FastaFile, HeaderOut):
    dnaSeq = ""
    i = 1
    for line in FastaFile:
        if re.match(">", line):
            HeaderOut.write("\t".join([str(i), line]))
            if dnaSeq != "":
                i = i + 1
                yield i, dnaSeq
                dnaSeq = ""
        else:
            dnaSeq = dnaSeq + line.strip()
    yield i, dnaSeq


def SNPit(id, seq, SNPsFile, VcfFile, opts):
    
    nindiv = opts.nindiv
    snpRate = opts.snprate
    #recomb = opts.recomb
    tstvRate = opts.tstv
    
    snpPos = []
    
    # randomly drawn the number of snps
    while(sum(snpPos) <= len(seq)):
        pos = np.random.geometric(p=snpRate,size=1)
        if ( pos + sum(snpPos) <= len(seq) ):
            snpPos.append(pos)
        else:
            break

    if(len(snpPos)==0):
        return 0

    nucList = list(seq)
    modelSeq = DivideRef(seq,tstvRate)
    numSNP = len(snpPos)

    output = subprocess.check_output(" ".join([opts.ms,
                                               str(nindiv*2),
                                               "1",
                                               "-s",
                                               str(numSNP),
                                               "-r",
                                               str(opts.rho),
                                               str(len(seq))]),
                                     universal_newlines=True,
                                     shell=True).split('\n')

    haplo_ready = False
    indx = 0
    haplMatrix = np.zeros(shape=(nindiv,2,numSNP), dtype=int)

    for line in output:
        if re.match("positions:", line):
            posLine = line.strip().split(" ")
            haplo_ready = True
            # Skip the 0 index because it's not a number - "Positions:"
            for i in range(1, len(posLine)):
                snpPos[i-1] = round(len(seq)*float(posLine[i]))
                if i>1 and snpPos[i-1] == snpPos[i-2]:
                    snpPos[i-1] = snpPos[i-1] +1
    
        elif haplo_ready:
            haplotype = [int(l) for l in list(line)]
            #print(int(indx/2), haplotype, indx, nindiv*2)
            haplMatrix[int(indx/2),indx%2,] = haplotype
            seq = RetrieveSeq(haplotype,snpPos,modelSeq)
            PrintIndivFasta(seq, SNPsFile[int(indx/2)], indx, id)
            indx = indx + 1
            if indx == nindiv*2:
                haplo_ready = False

    PrintRefFasta(modelSeq, SNPsFile[nindiv],id)
    PrintVcf(VcfFile, haplMatrix, snpPos, modelSeq, id)


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def RetrieveSeq(hap,snpPos,modelSeq):
    lenSeq = modelSeq.shape[1]
    adjPos = [i-1 for i in snpPos]
    
    indic = np.zeros(shape=(lenSeq), dtype=int)
    indic[adjPos] = hap
    
    return "".join(modelSeq[indic,range(lenSeq)])

def PrintIndivFasta(seq, SNPsFile, indx, id):

    template=""">contig_{id}_{hap}
{seq}\n"""
    
    context = {
    "seq": seq,
    "hap": indx%2,
    "id": id
    }

    SNPsFile.write(template.format(**context))

def PrintRefFasta(modelSeq, SNPsFile, id):
    
    template=""">contig_{id}_{hap}
{seq}\n"""
    
    context1 = {"seq": "".join(modelSeq[0]), "hap": "major", "id": id}
    context2 = {"seq": "".join(modelSeq[1]), "hap": "minor", "id": id}
    
    SNPsFile.write(template.format(**context1))
    #SNPsFile.write(template.format(**context2))

def PrintVcf(VcfFile, haplMatrix, snpPos, modelSeq, id):

    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"""
    numSNP = len(snpPos)
    hapM= haplMatrix.swapaxes(0,2).swapaxes(1,2)
    numIndiv = hapM.shape[1]
    for i in range(numSNP):
        vcfLine = "\t".join([str(id),
                             str(snpPos[i]),
                             "NA",
                             modelSeq[0,snpPos[i]-1],
                             modelSeq[1,snpPos[i]-1],
                             "100",
                             "PASS",
                             "NA",
                             "GT",
                             ""]) + "\t".join(["|".join([str(hapM[i,j,0]),
                                                         str(hapM[i,j,1])]) for j in range(numIndiv)]) + "\n"
        
        VcfFile.write(vcfLine)


# Nucleotide Dictionary :suffix 0 transition 1 - tranversion

nuclDict = {
'R' : ['G', 'A'],
'Y' : ['T',  'C'],
'K' : ['G', 'T'],
'M' : ['A',  'C'],
'S' : ['G',  'C'],
'W' : ['A' , 'T'],
'A0' : ['G'], 'G0' : ['A'],
'C0' : ['T'], 'T0' : ['C'],
'A1' : ['C', 'T'],
'C1' : ['A', 'G'],
'T1' : ['A', 'G'],
'G1' : ['C', 'T'],
'N' : ['A' ,'T', 'C', 'G']}

'''
    This function takes the sequence and transition:transversion
    ratio then returns two model sequences: one contains the dominant
    alleles while the other carries all minor alleles
    
    Current restrictions:
    Transiton vs Transversion: assume equal rate between transition or tranversion pairs
    limited to biallelic variants
    nucleotide alphabeta: only accepts: A T C G R Y K M S W N
    '''

def DivideRef(seq, tstvRate):
    
    #!!!! w/o the list(); pointer of seq is passed
    majorSeq = list(seq)
    minorSeq = list(seq)
    # 1 as Transversion
    tS = np.random.binomial(n=1, p=1.0/(1+tstvRate), size=len(seq))
    
    for i, n in enumerate(seq):
        if n == 'N':
            n = random.choice(nuclDict['N'])
            majorSeq[i] = n
        
        if re.match('[ATCG]',n):
            minorSeq[i] = random.choice(nuclDict[''.join([n,str(tS[i])])])
        
        else:
            pairs = random.sample(nuclDict[n],2)
            majorSeq[i] = pairs[0]
            minorSeq[i] = pairs[1]

    return np.array([majorSeq, minorSeq])




def makeRNFfiles(opts, id):

    path = "".join(["snake/indiv_",str(id)])
    if not os.path.exists(path):
        os.makedirs(path)
    
    SnakeFile = open(path+"/Snakefile",'w')

    template="""
import rnftools, smbl

consen_ref = "{basePath}/data/ref/ref.fasta"
indiv_ref =  "{basePath}/data/indiv_{id}.fasta"
reads = "{basePath}/sim/simu_{id}.fq"
bwa = "{basePath}/align/BWA-MEM_{id}.bam"
yara = "{basePath}/align/YARA_{id}.bam"


# READ SIMULATION

rnftools.mishmash.sample(reads[:-3],
                         reads_in_tuple=1)

rnftools.mishmash.DwgSim(fasta=indiv_ref,
                         read_length_1=100,
                         read_length_2=0,
                         number_of_read_tuples=10000, # might not be fixed
                         #haploid_mode = True, # currently the latest version contains bugs 
                         error_rate_1=0.001, #0.001-0.05; default : 0.02
                         mutation_rate =0.0001, # default : 0.001
                         indels =0.1, # default: 0.1
                         prob_indel_ext=0.1 # default: 0.3,
                         other_params="-H",
                         )
                         
# -y FLOAT      probability of a random DNA read [0.05]
# see https://github.com/nh13/DWGSIM/wiki/Simulating-Reads-with-DWGSIM for more options

# ALIGNMENTS

alignments = [
          smbl.prog.BwaMem(
                           fasta=consen_ref,
                           fastq_1=reads,
                           bam=bwa,
                           ),
          
          smbl.prog.Yara(
                         fasta=consen_ref,
                         fastq_1=reads,
                         bam=yara,
                         ),
            
          ]

# SNAKEMAKE RULES

rule basic:
    input: [aln.bam_fn() for aln in alignments]
include: rnftools.include()

"""
    context = {
    "basePath": os.path.abspath(""),
    "input": opts.input,
    "id": id
    }
    
    SnakeFile.write(template.format(**context))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='simulate SNPs on ')
    parser.add_argument('-i','--input',required=True, default=None,type=str,help='relative input FASTA file path')
    parser.add_argument('-n','--nindiv',required=False, default=4,type=int,help='number of individuals')
    parser.add_argument('-s','--snprate',required=False, default=0.01,type=float,help='expected number of SNPs per bp')
    #parser.add_argument('-r','--recomb',required=False, default=0.5,type=double,help='recombination rate (0 to 0.5) ')
    parser.add_argument('-r','--rho',required=False, default=0,type=float,help='recombination param - 4Nr')
    parser.add_argument('-t','--tstv',required=False, default=3,type=float,help='transition to transversion ratio (>0)')
    parser.add_argument('-ms','--ms',required=True, default=None, type=str, help='local ms path')
    opts = parser.parse_args()
    #/Users/work/academic/anderson/BattleOfGenotype/src/ms.folder/msdir/ms
    
    for i in ["data", "snake", "align", "data/ref", "sim"]:
        if not os.path.exists(i):
            os.makedirs(i)
    
    FastaFile = open(opts.input,'r')
    HeaderOut = open("data/fasta_header.csv", 'w')
    VcfFile = open("data/veritas.vcf", 'w')
    
    
    WriteVcfHeader(VcfFile, opts)
    
    AllFasta = []
    for i in range(0,opts.nindiv):
        AllFasta.append(open("".join(["data/indiv_",str(i),".fasta"]), 'w'))
    AllFasta.append(open("data/ref/ref.fasta", 'w'))
    
    dnaSeqs = readFasta(FastaFile, HeaderOut)

    for i, dnaSeq in dnaSeqs:
        SNPit(i, dnaSeq, AllFasta, VcfFile, opts)
    
    
    for i in range(opts.nindiv):
        makeRNFfiles(opts, i)
    
    FastaFile.close()
    HeaderOut.close()
    VcfFile.close()
    for i in range(0,opts.nindiv):
        AllFasta[i].close()
