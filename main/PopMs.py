#!/usr/bin/env python


# version required Python3
# numpy should be > v. 1.8
# python ~/src/BattleOfGenotype/main/PopMs.py -i ~/data/ref/satro_orig_amplicons.fa -ms ~/src/ms/msdir/ms -qp ~/data/ddrad/lines.fq -mR 50000
# python  ~/academic/anderson/BattleOfGenotype/src/BattleOfGenotype/main/PopMs.py -i /Users/work/academic/anderson/BattleOfGenotype/data/satro_orig_amplicons.fa -ms /Users/work/academic/anderson/BattleOfGenotype/src/ms.folder/msdir/ms -qp /Users/work/academic/anderson/BattleOfGenotype/data/satro_ddrad_S10_R1_400000_lines.fq

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

import os, sys, re, argparse, subprocess, string, gc
import numpy as np
from scipy.stats import nbinom, beta
import random
import multiprocessing as mp

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
'nA': ['C', 'G', 'T'],
'nC': ['A', 'G', 'T'],
'nG': ['A', 'C', 'T'],
'nT': ['C', 'G', 'A'],
'B': ['C', 'G', 'T'],
'D': ['A', 'G', 'T'],
'H': ['A', 'C', 'T'],
'V': ['C', 'G', 'A'],
'N' : ['A' ,'T', 'C', 'G']}


''' 
    This function counts the number of fasta entries given the Fasta file path
'''

def GetNumEntries(inFASTA):

    cmd = " ".join(["grep '>' ",
                    inFASTA,
                    "| wc -l"])
    numLine = subprocess.check_output(cmd,shell=True)
    return int(numLine)


def GetPhred(fastQCpath):
   
    FASTQC = open(fastQCpath, 'r')
    
    #get the maximum length of the read sequence
    
    i=1
    readLen = 0;
    for line in FASTQC:
        if (i%4==0):
            a = len(line.strip())
            readLen = max(a, readLen)
        i = i + 1;
        if(i>41):
            break

    qualMatrix = np.empty([readLen,45])
    qualMatrix[:] = 0

    for line in FASTQC:
        if (i%4==0):
            for pos, ascii in enumerate(list(line.strip())):
                indx = min(ord(ascii)-33,44)
                try:
                    qualMatrix[pos,indx] += 1
                except IndexError:
                    pass
        i = i + 1;

    for pos in range(readLen):
        qualMatrix[pos,] = qualMatrix[pos,]*1.0/sum(qualMatrix[pos,])

    return qualMatrix


# This function writes the beginning header for the VCF file
def WriteVcfHeader(VcfFile, opts):

    template="""##fileformat=VCFv4.1
##fileData=Veritas
##reference={fastaFile}
##phasing=known
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"""
    
    template = template + "\t".join(["".join(["INDIV_",str(i+1)]) for i in range(0,opts.nindiv)])+ "\n"
    
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
    i = 0
    for line in FastaFile:
        if re.match(">", line):
            i = i + 1
            HeaderOut.write("\t".join([str(i), line]))
            if dnaSeq != "":
                yield i-1, dnaSeq
                dnaSeq = ""
        else:
            dnaSeq = dnaSeq + line.strip()

    if dnaSeq != "":
        yield i, dnaSeq

'''
     Given a single fasta entry as a template, this module generates a fixed number of 
     haplotype sequences with derived alleles based on Coalescent Model (MS). 
     This module depends on three output functions:       
     1) PrintIndivFasta - printing individual Fasta sequence
     2) PrintRefFasta
     3) PrintVcf
'''


def SNPit(id, coverageMatrix, seq, SNPsFile, VcfFile, opts):
    
    nindiv = opts.nindiv
    snpRate = opts.snprate
    #recomb = opts.recomb
    tstvRate = opts.tstv
    indelRate = opts.ir
    extendRate = opts.er    


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
    modelSeq = DivideRef(seq, tstvRate, indelRate, extendRate)
    numSNP = len(snpPos)

    ## passing ms terminal command line via python subprocess
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

    # Storing all the haplotype information in the form a matrix 
    # The reason for storing in the format is that I need to reguritate this info
    # in vcf format
    haplMatrix = np.zeros(shape=(nindiv,2,numSNP), dtype=int)
    
    # Parsing out MS output
    for line in output:
        # I am looking out for line that contain any information about the variant positions
        if re.match("positions:", line):
            posLine = line.strip().split(" ")
            haplo_ready = True
            # Skip the 0 index because it's not a number - "Positions:"
            for i in range(1, len(posLine)):
                snpPos[i-1] = round(len(seq)*float(posLine[i]))
                if i>1 and snpPos[i-1] == snpPos[i-2]:
                    snpPos[i-1] = snpPos[i-1] +1
    
        # this section parses ms' output line on displaying individual's specificed haplotype configuration e.g. 1001101
        elif haplo_ready:
            haplotype = [int(l) for l in list(line)]
            #print(int(indx/2), haplotype, indx, nindiv*2)
            
            # storing each haplotype into a matrix
            haplMatrix[int(indx/2),indx%2,] = haplotype

            indivSeq = RetrieveSeq(haplotype,snpPos,modelSeq)
            PrintIndivFasta(id, indivSeq, indx, SNPsFile[int(indx/2)], coverageMatrix[int(indx/2),indx%2, id-1])
            indx = indx + 1

            ## stop the process when the number of individuals needed to print out is fulfilled
            if indx == nindiv*2:
                break

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

def PrintIndivFasta(id, seq, indx, SNPsFile, numRepeat):

    template=""">contig_{id}_{hap}
{seq}\n"""
    
    context = {
    "seq": seq,
    "hap": indx%2,
    "id": id
    }

    SNPsFile.write(template.format(**context)*numRepeat)

def PrintRefFasta(modelSeq, SNPsFile, id):
    
    template=""">contig_{id}
{seq}\n"""
    
    #context1 = {"seq": "".join(modelSeq[0]), "hap": "major", "id": id}
    #context2 = {"seq": "".join(modelSeq[1]), "hap": "minor", "id": id}

    context1 = {"seq": "".join(modelSeq[0]), "id": id}
    #context2 = {"seq": "".join(modelSeq[0]), "hap": 1, "id": id}
    
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



'''
    This function takes the sequence and transition:transversion
    ratio then returns two model sequences: one contains the dominant
    alleles while the other carries all minor alleles
    
    Current restrictions:
    Transiton vs Transversion: assume equal rate between transition or tranversion pairs
    limited to biallelic variants
    nucleotide alphabeta: only accepts: A T C G R Y K M S W N B D H V
    '''

def DivideRef(seq, tstvRate, indelRate, extendRate):
    
    #!!!! w/o the list(); pointer of seq is passed
    majorSeq = list(seq)
    minorSeq = list(seq)

    # 1 as Transversion
    tS = np.random.binomial(n=1, p=1.0/(1+tstvRate), size=len(seq))

    # Shall I include insertion/deletion in place of substitution? 1 - yes
    indelI = np.random.binomial(n=1, p=indelRate, size=len(seq))
    # If it's an indel, is it an extension? 1 - yes
    isInsert = np.random.binomial(n=1, p=0.5, size=len(seq))

    # How long will the base extended for an insertion?
    lenExt = np.random.geometric(p=1.0-extendRate, size=len(seq))-1        

    for i, n in enumerate(seq):
        if n == 'N':
            n = random.choice(nuclDict[n])
            majorSeq[i] = n
        
        if re.match('[ATCG]',n):
            minorSeq[i] = random.choice(nuclDict[''.join([n,str(tS[i])])])
        else:
            pairs = random.sample(nuclDict[n],2)
            majorSeq[i] = pairs[0]
            minorSeq[i] = pairs[1]

        if indelI[i]==1:
            if isInsert[i] == 1:
                minorSeq[i] = minorSeq[i] + "".join(np.random.choice(nuclDict['N'],lenExt[i]))
            else:
                minorSeq[i] = ''

    return np.array([majorSeq, minorSeq])


# I will not be using any of the MISMASH read simulator modules since none of them are
# flexible enough to handle the type of desired simulated reads

# For example: 
# None of the programs provide any option in specifying the starting position of read samples
# They are perfect for shotgun genomic simulator but not for simulating results from site-directed library assay (like ddrad)
# Art Ilumina and Mason Illumina : difficulties in generating >250 reads 
# DwgSim: Failure in generating mock read quality profiles

def makeRNFfiles(opts, id):

    path = "".join(["data/rnf/snake/indiv_",str(id)])
    if not os.path.exists(path):
        os.makedirs(path)
    
    SnakeFile = open(path+"/Snakefile",'w')

    mean =opts.meanRead
    variance =opts.varRead
    shape = mean*mean /variance
    scale = variance /mean
    nReads= np.random.gamma(shape, scale, 1)

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

#rnftools.mishmash.MasonIllumina(fasta=indiv_ref,
#rnftools.mishmash.DwgSim(fasta=indiv_ref,
rnftools.mishmash.ArtIllumina(fasta=indiv_ref,
                         read_length_1={lenR},
                         read_length_2=0,
                        #coverage=4,
			 number_of_read_tuples={numReads}, # might not be fixed
                         #haploid_mode = True, # currently the latest version contains bugs 
                         #error_rate_1=0.02, #0.001-0.05; default : 0.02
                         #mutation_rate =0.0001, # default : 0.001
                         #indels =0.1, # default: 0.1
                         #prob_indel_ext=0.1, # default: 0.3,
                         #other_params="-H",
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
    "id": id,
    "lenR": int(opts.len),
    "numReads":int(nReads)
    }
    
    SnakeFile.write(template.format(**context))


def CalculateCoverage(opts, numLoci):

    # dimension: num of Indiv X haplotype X num of Loci 
    coverageMatrix = np.zeros(shape=(opts.nindiv, 2, numLoci), dtype=int)

    mean = opts.meanRead
    variance = opts.varRead
    shape = mean*mean /variance
    scale = variance /mean
    nReads= np.random.gamma(shape, scale, opts.nindiv)

    if opts.logNormal:
        lociWeight = np.random.lognormal(opts.lm, opts.lsd, numLoci) 
    else:
        lociWeight = np.random.gamma(opts.ga, 1/opts.gb, numLoci)

    lociProb = np.random.dirichlet(alpha=lociWeight, size=1)
    for i in range(opts.nindiv):
        coverageMatrix[i,1,] = np.random.multinomial(n=nReads[i], pvals=lociProb[0], size=1)
        haploP = np.random.beta(opts.ha, opts.hb, numLoci)
        for j in range(numLoci):
            coverageMatrix[i,0,j] = np.random.binomial(n=coverageMatrix[i,1,j],p=haploP[j],size=1)
            coverageMatrix[i,1,j] = coverageMatrix[i,1,j] - coverageMatrix[i,0,j]
    return coverageMatrix

def PrintCoverageMatrix(coverageMatrix, CVFILE):

    flatIter = coverageMatrix.flat
    indxIter = np.ndindex(coverageMatrix.shape)

    CVFILE.write("#indiv,haplotype,locus,numReads\n")
    for i in flatIter:
        indx = indxIter.next()
        
        CVFILE.write(",".join([str(indx[0]+1),
                              str(indx[1]),
                              str(indx[2]+1),
                              str(i)+"\n"]))


def RoundPos(i): return int(max(i,0))

def ConvFastaToFastQc(indivID, phredMatrix, coverageMatrix, indelRate, readLen):

    random.seed()
    np.random.seed()
    print("Processing Individual ",indivID+1," \n");
    FASTA = open("".join(["data/seq/indiv_",str(indivID+1),".fasta"]), 'r')
    FASTQC = open("".join(["data/seq/indiv_",str(indivID+1),".fq"]), 'w')
    ERROR = open("".join(["data/seq/seqerr_",str(indivID+1),".txt"]), 'w')
    
    processSize = 100000#10000#100000
    totReads = sum(sum(coverageMatrix[indivID,]))
    numReads = min(processSize, totReads)

    predLen = phredMatrix.shape[0]
    adjIndx = [ max(0, round((indx+1) * predLen/readLen)-1) for indx in range(readLen)]

    indivPredMatrix = np.empty([readLen, numReads], dtype="int")
    errorIndic = np.empty([45, readLen, numReads], dtype="bool")

    roundToPos = np.vectorize(RoundPos)
    
    for i in range(readLen): 
        indivPredMatrix[i,] = np.random.choice(45, p=phredMatrix[adjIndx[i],], size=numReads)
    for j in range(45):
        errorIndic[j,] = np.random.binomial(1, pow(10, j/-10.0), readLen * numReads).reshape(readLen,numReads)

    adjPhredMatrix = roundToPos(indivPredMatrix + (np.random.sample(readLen * numReads)*3- 1.5).reshape(readLen,numReads))

    ct = 0
    for l, line in enumerate(FASTA):
    
        if(l%2==0):
            FASTQC.write(line.replace(">","@"+str(int((l/2)+1))+"_" ).strip() + "_" + str(indivID+1) +"\n")
        else:
            seq = list(line.strip()) # if i choose numpy.array beware of the dtype since some characters might get trimmed out pending on the dtype
            acceptLen = min(len(seq),readLen)
            seq = seq[:acceptLen]
            
            try:            
                phredScore = indivPredMatrix[:acceptLen,ct]
            except IndexError:
                print ("==", totReads-int(l/2), " vs ", processSize);
                numReadRevised = min(processSize, totReads-int(l/2))
       
                del indivPredMatrix
                del errorIndic
                #gc.collect()
                print("Generating ", int(l/2), "reads for individual ", indivID+1)

                if (numReadRevised != numReads):
                    numReads = numReadRevised
                
                print("ReadLen: ", readLen, " numReads: ", numReads, " for indiv ", indivID+1, "\n")
                indivPredMatrix = np.empty([readLen, numReads], dtype="int")
                errorIndic = np.empty([45, readLen, numReads], dtype="bool")

                for i in range(readLen): 
                    indivPredMatrix[i,] = np.random.choice(45, p=phredMatrix[adjIndx[i],], size=numReads)
                for j in range(45):
                    errorIndic[j,] = np.random.binomial(1, pow(10, j/-10.0), readLen * numReads).reshape(readLen,numReads)

                adjPhredMatrix = roundToPos(indivPredMatrix + (np.random.sample(readLen * numReads)*3- 1.5).reshape(readLen,numReads))
                ct = 0
                phredScore = indivPredMatrix[:acceptLen,ct]



            #phredAdj = phredScore + (np.random.sample(acceptLen)-0.5)
           
            qualSeq = [chr(i+33) for i in phredScore]
            SeqErrorI = errorIndic[adjPhredMatrix[:acceptLen,ct], np.arange(acceptLen), ct]
            
            # 1 means invite sequence error!!
            #SeqErrorI = [np.random.binomial(1, pow(10, max(phred,0)/-10.0)) for phred in phredAdj]
            SeqIndex = np.where(SeqErrorI==1)
            
            for index in SeqIndex[0]:
                ERROR.write("Entries:\t{}\tpos:\t{}\t{}\t".format(int(l/2)+1, index, seq[index]))
                if np.random.binomial(1, indelRate)==0:
                    seq[index] = random.choice(nuclDict["n"+seq[index]])
                else:
                    if np.random.binomial(1,0.5):
                        seq[index] = "" # deletion
                        qualSeq[index] = ""
                    else:
                        seq[index] = random.choice(nuclDict["n"+seq[index]]) + random.choice(nuclDict["N"])
                        qualSeq[index]=qualSeq[index] + qualSeq[index]

                ERROR.write("{}\t{}\n".format(seq[index], qualSeq[index]))
            
            template="""{seq}
+
{qual}\n"""
            context = {
            "seq":"".join(seq),
            "qual":"".join(qualSeq),
            }
            FASTQC.write(template.format(**context))
            ct = ct + 1

    del indivPredMatrix
    del errorIndic
    FASTA.close()
    FASTQC.close()
    ERROR.close()

#bwa index -p Satrovirens_amplicons_gtseq3 -a is Satrovirens_amplicons_gtseq3.fa
#http://161.55.237.25/~newmedusa/dokuwiki/doku.php?id=projects:rockfish_gtseq_run3_25may2015
#for i in {1..96}; do bwa mem -aM -v 3 -t 12 -R "@RG\tID:s${i}\tLB:amplicon\tPL:ILLUMINA\tSM:rock${i}" ./Satrovirens_amplicons_gtseq3 ../trimfilter/rock_S${i}_L001_R1_001_val_1.fq.gz ../trimfilter/rock_S${i}_L001_R2_001_val_2.fq.gz > ./satro_s${i}_aln.sam; done
#for i in {1..96}; do samtools view -bS satro_s${i}_aln.sam > satro_s${i}.bam; done

#~/src/bowtie/bowtie2-2.2.5/bowtie2-build data/veritas/ref_sequence.fasta analysis/align/ref_sequence
#~/src/bowtie/bowtie2-2.2.5/bowtie2 -x analysis/align/ref_sequence -q data/seq/indiv_1.fq -S analysis/align/indiv_1_bowtie.sam

#~/src/bwa/bwa/bwa mem -a -v 1 -t 1 -R "@RG\tID:s1\tLB:amplicon\tPL:ILLUMINA\tSM:rock1" analysis/align/ref_sequence data/seq/indiv_1.fq > analysis/align/indiv_1_aln.sam0

def MakeBWAIndex():
    cmd = "bwa index -a is -p analysis/align/ref_sequence data/veritas/ref_sequence.fasta"
    subprocess.check_output(cmd,shell=True)

def AlignIndiv(indivID):
    cmd = "".join(['bwa mem -aM -v 1 -t 1 -R "@RG\tID:s',str(indivID+1),'\tLB:amplicon\tPL:ILLUMINA\tSM:rock',str(indivID+1),'" analysis/align/ref_sequence data/seq/indiv_',str(indivID+1),'.fq > analysis/align/indiv_',str(indivID+1),'_aln.sam'])
    subprocess.check_output(cmd,shell=True)
    cmd = "".join(['samtools view -bS analysis/align/indiv_',str(indivID+1),'_aln.sam > analysis/align/indiv_',str(indivID+1),'.bam'])
    subprocess.check_output(cmd,shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='simulate population-based ddrad short read illumina sequence from consensus sequences')
    parser.add_argument('-i','--input',required=True, default=None,type=str,help='relative input FASTA file path')
    parser.add_argument('-ms','--ms',required=True, default=None, type=str, help='local ms path')
    parser.add_argument('-n','--nindiv',required=False, default=4,type=int,help='number of individuals')
    parser.add_argument('-l','--len',required=False, default=300, type=int, help='simulated read length')
    #parser.add_argument('-p','--paired',required=False, default=False, type=bool, help='paired read options')
    #parser.add_argument(-seed, '--seed', required=False, default=0, type=int, help="set random seed")

    # biological/true mutation and variation
    # assumption: limit to biallelic variants for SNPs
    
    parser.add_argument('-s','--snprate',default=0.01,type=float,help='expected number of SNPs per bp')
    #parser.add_argument('-r','--recomb',required=False, default=0.5,type=double,help='recombination rate (0 to 0.5) ')
    parser.add_argument('-r','--rho', default=0,type=float,help='recombination param - 4Nr')
    parser.add_argument('-t','--tstv', default=3,type=float,help='transition to transversion ratio (>0)')
    #parser.add_argument('-smr','--smr',required=False, default=0.0001,type=float,help='expected bp rate for somatic or gamete mutation') #ignoring this for now; since it is such an extreme infrequent case
    parser.add_argument('-ir', '--ir', default=0.05, type=float, help='fraction of mutations that are indels')
    parser.add_argument('-er', '--er', default=0.3, type=float, help='Prob that an indel is extended')    

    # Experimental and sequencing errors
    parser.add_argument('-qp', '--qp', default=None, type=str, help='path for fq file to construct quality profile ')    
    parser.add_argument('-se', '--se', default=0.01, type=float, help='fraction of sequencing errors that are indels')

    # Parameters for modulating different layers in read representations

       # inidividual coverage
    parser.add_argument('-mR','--meanRead', default=100000, type=int, help='expected number of total reads per individual')
    parser.add_argument('-vR','--varRead', default=15000, type=int, help='variance of number of total reads per individual')

       # haplotype effect
    parser.add_argument('-ha','--ha', default=2, type=float, help='alpha parameter of the binomial-beta that instantiates haplotype bias')
    parser.add_argument('-hb','--hb', default=2, type=float, help='beta parameter of the binomial-beta that instantiates haplotype bias')

       # locus effect: either described with gamma prior to support multinomial-dirchlet or log-normal (log-normal: will be used as default)
    parser.add_argument('-gm','--gModel', action='store_true', help='locus effect will be based on higher level gamma model') 
    parser.add_argument('-ga','--ga', default=620, type=float, help='alpha parameter for gamma distrib that supports each locus-based dirichlet weight')
    parser.add_argument('-gb','--gb', default=1.0/700, type=float, help='beta parameter for gamma distrib that supports each locu-based dirichlet weight ')

    parser.add_argument('-ln','--logNormal', action='store_true', help='locus effect will be based on higher level gamma model (default) ') 
    parser.add_argument('-lm','--lm', default=2, type=float, help='alpha parameter for gamma distrib that supports each locus-based dirichlet weight')
    parser.add_argument('-lsd','--lsd', default=0.1, type=float, help='beta parameter for gamma distrib that supports each locu-based dirichlet weight ')

    opts = parser.parse_args()
    #/Users/work/academic/anderson/BattleOfGenotype/src/ms.folder/msdir/ms
    if(opts.gModel + opts.logNormal == 0): 
        opts.logNormal = True 
    
    for i in ["data", "data/rnf", "data/seq", "data/veritas",
              "analysis", "analysis/align", "analysis/geno"]:
        if not os.path.exists(i):
            os.makedirs(i)
    
    FASTAFILE = open(opts.input,'r')
    HEADERFILE = open("data/veritas/header_info.csv", 'w')
    VCFFILE = open("data/veritas/genotype.vcf", 'w')
    CVFILE = open("data/veritas/coverage.csv", 'w')

    AllFasta = []
    for i in range(0,opts.nindiv):
        AllFasta.append(open("".join(["data/seq/indiv_",str(i+1),".fasta"]), 'w'))
    AllFasta.append(open("data/veritas/ref_sequence.fasta", 'w'))

    numEntries = GetNumEntries(opts.input)
    #print("Number of Entries ", numEntries)
    WriteVcfHeader(VCFFILE, opts)

    coverageMatrix = CalculateCoverage(opts, numEntries)
    PrintCoverageMatrix(coverageMatrix, CVFILE)
    np.save("data/veritas/coverage_matrix.npy", coverageMatrix)
    
    dnaSeqs = readFasta(FASTAFILE, HEADERFILE)
    for i, dnaSeq in dnaSeqs:
        #print(str(i),"--");
        SNPit(i, coverageMatrix, dnaSeq, AllFasta, VCFFILE, opts)

    # closing all file handler
    FASTAFILE.close()
    HEADERFILE.close()
    VCFFILE.close()
    CVFILE.close()
    for i in range(0,opts.nindiv+1):
        AllFasta[i].close()

    # Create quality profile for each individual FASTA file
    if opts.qp != None:
        phredMatrix = GetPhred(opts.qp)
        np.savetxt("data/veritas/quality_profile.csv", phredMatrix, delimiter=',')
        np.save("data/veritas/300bp_quality_profile.npy", phredMatrix)
    else:
        phredMatrix = np.load("data/veritas/300bp_quality_profile.npy")

    pool = mp.Pool(processes = 4)
    results = [pool.apply_async(ConvFastaToFastQc, args=(i, phredMatrix, coverageMatrix, opts.se, opts.len)) for i in range(opts.nindiv)]
    for r in results:
        r.get()

    print("Alignment Time")
    # Alignment
    MakeBWAIndex()
    results = [pool.apply_async(AlignIndiv, args=[i]) for i in range(opts.nindiv)]
    for r in results:
        r.get()



    #for i in range(opts.nindiv):
    #    makeRNFfiles(opts, i)

