#!/usr/bin/env python3

import sys, os
import argparse
import time
import logging
import re

from numpy import min_scalar_type

logging.basicConfig(level=logging.DEBUG)

# def getArgs(argList, name):
def getArgs():
    p = argparse.ArgumentParser()

    p.add_argument('-v', "--verbose", action="store_true", help="Enable verbose debugging", default=False)
    p.add_argument('--fasta',   action="store", help="fast file name", required=True)
    p.add_argument('--sam',     action="store", help="sam file name", required=True)
    p.add_argument('--minS',    action="store", type=int, help="min of gap", required=True)
    # p.add_argument('--blastDB', action="store", help="blast database", required=True)
    p.add_argument('--db2vec',  action="store_true", help="blast direction, vector to Arabidopsis")
    # p.add_argument('--blastTH', action="store", help="blast thread", required=True)
    
    return p.parse_args()

def proc_sam_files(args):
    # logging.debug(args)

    ### line 38 of tdnaFinder.py
    
    samFile        = open(args.sam)
    minS           = args.minS
    trimSFile      = open("trimS.sam.1","w")
    trimSFastaFile = open("trimSFasta.fa.1","w")
    # newSamFile     = open("new.sam","w")

    onlyNumber = ['0','1','2','3','4','5','6','7','8','9']
    for oneLine in samFile:
        # logging.debug(oneLine)
        if oneLine.startswith("@"):
            trimSFile.write(oneLine)
            continue

        # sam_fields: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
        data = oneLine.split("\t")
        # logging.debug(data)
        
        # findS = data[5].split("S")
        
        if data[5].count("S") != 1:
            continue
        # logging.debug(data[5].count("S"))

        ops = re.findall(r"\*|[0-9]+[MIDNSHPX=]", data[5])
        sequence = data[9]
        # sCount = 0
        # logging.debug(len(sequence))
        # logging.debug(ops)

        ### ref ck_filterSAM.pl: sub proc_cigar{}
        pMaxS = 0
        maxS = 0
        pos = 0
        for e in (ops):
            m = re.match("([0-9]+)(.)", e)
            l, op = int(m.group(1)), m.group(2)
            # logging.debug("%s:%s" % (l, op))
            if op == 'S':
                if l > maxS:
                    pMaxS, maxS = pos, l
            else:
                pos += l

        trimS = sequence[pMaxS:(pMaxS + maxS)]
        # logging.debug("%s->%s" % (pMaxS, pMaxS + maxS))
        # print("%s->%s" % (pMaxS, pMaxS + maxS))
        # print("%s: %s" % (pMaxS, trimS))
        
        trimSameChar = 10.0 # %
        sCount = len(trimS)
        if(sCount > minS):
            maxSameChar = len(trimS) / 100.0 * trimSameChar

            ch = ''
            cnt = 0
            save = True
            for c in trimS:
                if ch != c:
                    ch = c
                    cnt = 1
                else:
                    cnt += 1
                    if cnt >= maxSameChar:
                        save = False
                        break
            if save:
                trimSFastaFile.write(">"+data[0]+"\n"+trimS+"\n")
                trimSFile.write(oneLine)
            
    samFile.close()
    trimSFile.close()
    trimSFastaFile.close()

def proc_blast(args):
    pass
    ### line 24 of tdnaFinder.py
    
    # if not args.db2vec:
    #     choseOne = "Arabidopsis"
    #     otherOne = "vector"
    #     otherFile    = open("Vector_Blast_result.txt","w")
    #     needKeepFile = open("Arabidopsis_Blast_result.txt","w")
    # else:
    #     choseOne = "vector"
    #     otherOne = "Arabidopsis"
    #     needKeepFile = open("Vector_Blast_result.txt","w")
    #     otherFile    = open("Arabidopsis_Blast_result.txt","w")
    
    # print ("./blastn -query trimSFasta.fa -db "+blastDatabase+" -evalue 1e-10 -num_threads "+str(blastThread)+" -outfmt '6 qseqid salltitles' > Blast.out")
    # subprocess.call("./blastn -query trimSFasta.fa -db "+blastDatabase+" -evalue 1e-10 -num_threads "+str(blastThread)+" -outfmt '6 qseqid salltitles' > Blast.out",shell=True)
    # blastFile = open("Blast.out")
    # trimSFile = open("trimS.sam")

    # contigCount   = {}
    # whichINeed    = []
    # whichINotNeed = []
    # while(True):
    # 	oneLine =blastFile.readline()
    # 	if(oneLine==""):
    # 		break
    # 	data = oneLine.split("\t")
    # 	name = data[0]
    # 	blastData = data[1]
    # 	if(not contigCount.has_key(name)):
    # 		contigCount[name] = 1
    # 		if(blastData.find(choseOne) >= 0):
    # 			needKeepFile.write(oneLine)
    # 			if(name not in whichINeed):
    # 				whichINeed.append(name)
    # 		elif(blastData.find(otherOne) >= 0):
    # 			otherFile.write(oneLine)
    # 			if(name not in whichINotNeed):
    # 				whichINotNeed.append(name)

    # 	elif(contigCount[name] < topNum):
    # 		contigCount[name]+= 1
    # 		if(blastData.find(choseOne) >= 0):
    # 			needKeepFile.write(oneLine)
    # 			if(name not in whichINeed):
    # 				whichINeed.append(name)
    # 		elif(blastData.find(otherOne) >= 0):
    # 			otherFile.write(oneLine)
    # 			if(name not in whichINotNeed):
    # 				whichINotNeed.append(name)
    # 	elif(contigCount[name]==topNum): # it have Arabidopsis or vactor which I need 
    # 		continue

	
    # blastFile.close()
    # needKeepFile.close()
    # otherFile.close()

    # while(True):
    # 	oneLine = trimSFile.readline()
    # 	if(oneLine==""):
    # 		break
    # 	data = oneLine.split("\t")
    # 	if(len(data)<13):
    # 		newSamFile.write(oneLine)
    # 	name = data[0]
    # 	if(name in whichINeed):
    # 		newSamFile.write(oneLine)
    # trimSFile.close()
    # newSamFile.close()

    # subprocess.call("perl FilterSAM.pl "+fastaName+" new.sam "+str(minS),shell=True)
    # subprocess.call("rm -rf merge.ps",shell=True)
    # subprocess.call("rm -rf psFile",shell=True)
    # subprocess.call("rm -rf new.sam",shell=True)
    # subprocess.call("rm -rf trim_info.txt",shell=True)
    # subprocess.call("rm -rf trim_S.fa",shell=True)
    # subprocess.call("rm -rf trimSFasta.fa",shell=True)
    # subprocess.call("rm -rf trimS.sam",shell=True)

def main():
    args = getArgs()

    # print("Starting %s" % (sys.argv[0]))
    startTime = float(time.time())

    proc_sam_files(args)
    proc_blast(args)

    print("Finished in %0.4f seconds" % (time.time() - startTime))
    return

if __name__ == '__main__':
    #sys.argv = ["programName.py","--input","test.txt","--output","tmp/test.txt"]
    main()
    