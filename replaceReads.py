####
# Modifies a certain proportion of reads at a specified location
# Only reads with a string of 'M' (matches) running across the modification site are considered for swapping. Other reads are deleted
#quality scores are left as is
# some notes:
# 1) Reads that have the deletion too close to one of the ends are discarded
# 2) nucleotides are added to the end of the sequence. If the orignal sequence had soft clipping a the end, the added nucleotides are added as clipped
# 3) All reads are paired - otherwise they weren't included in bedtools bam2fastq
####
import argparse
import pysam
from Bio import SeqIO
import random
import re
import os
import copy
from replaceReadsUtils import replaceRead

#print("SETTING SEED")
#random.seed(123)

#cigar variables
MATCH = 0
INSERTION = 1
DELETION = 2
CLIPPING = 4

fileroot = "/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/simulations/circle/"
argparser = argparse.ArgumentParser()
argparser.add_argument("--downsampleNumber",help="Depth of downsampling",type=int,default=100)
argparser.add_argument("--swapFreq",help="Frequency of swapping reads from alternate into sample",type=float,default=.10)
argparser.add_argument("--qualAdd",help="This number is added to all quality scores",type=int,default=0)
argparser.add_argument("--swapChr",help="Swap chromosome",default="chr2")
argparser.add_argument("--swapLoc",help="Swap site",type=int,default=72933870)
argparser.add_argument("--reference",help="Reference genome location",default="/data/pinello/COMMON_DATA/REFERENCE_GENOMES/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa")
argparser.add_argument("--outfile",help="Output file",required=True)
argparser.add_argument("--unalteredBam",help="Source BWA-aligned file with unaltered reads",
		default= fileroot + "simulations/alignNA12878/1kbAroundCutsite.BWA.bam")
argparser.add_argument("--unalteredNamesortedBam",help="Source BWA-aligned file with unaltered reads, name sorted using samtools sort -n",
		default= fileroot + "simulations/alignNA12878/1kbAroundCutsite.BWA.nameSort.bam")
argparser.add_argument("--alteredBam",help="BWA-aligned file with altered reads",
		default= fileroot + "circle/SRR1046762.bam")
argparser.add_argument("--onlyIncludeAlteredWithIndel",help="Of reads from the altered file, only those with an indel are considered for swapping in alterations",action='store_true')
args = argparser.parse_args()

print(args)


#get reference sequence
refStart = args.swapLoc - 5000
refEnd = args.swapLoc + 5000
refFile = pysam.Fastafile(args.reference)
refSeq = refFile.fetch(args.swapChr,refStart,refEnd)



sourceAln = pysam.AlignmentFile(args.unalteredBam,"rb")

cigarExplodePattern = re.compile("(\d+)(\w)")
def explodeCigar(cigarString):
	explodedCigar = ""
	for (cigarCount,cigarChar) in re.findall(cigarExplodePattern,cigarString):
		explodedCigar += cigarChar * int(cigarCount)
	return explodedCigar

cigarUnexplodePattern = re.compile(r'((\w)\2{0,})')
def unexplodeCigar(explodedCigarString):
	cigar = ""
	for (cigarStr,cigarChar) in re.findall(cigarUnexplodePattern,explodedCigarString):
		cigar += str(len(cigarStr)) + cigarChar
	return cigar


unsortedOutName = ""
unsortedCtlOutName = ""
mutationsFileName = ""
#open with h == header mode
if(args.outfile.endswith(".bam")):
	unsortedOutName = args.outfile+".unsorted.bam"
	destAln = pysam.AlignmentFile(unsortedOutName,"wb",template=sourceAln)
	unsortedCtlOutName = args.outfile+".ctl.unsorted.bam"
	destAlnCtl = pysam.AlignmentFile(unsortedCtlOutName,"wb",template=sourceAln)
else:
	unsortedOutName = args.outfile+".unsorted.sam"
	destAln = pysam.AlignmentFile(unsortedOutName,"wh",template=sourceAln)
	unsortedCtlOutName = args.outfile+".ctl.unsorted.sam"
	destAlnCtl = pysam.AlignmentFile(unsortedCtlOutName,"wb",template=sourceAln)

mutatedReadsFileName = args.outfile+".mutatedReads.txt"
mutationsFileName = args.outfile+".mutations.txt"

olapBuffer = args.swapLoc + 5 #don't modify reads that don't overlap completely with this sequence

#####
#first, get reads that overlap with deletion site, make sure they are good quality, and then add them to a list
allReadAtTargetCount = 0
goodReadAtTargetCount = 0
readsAtTarget = []
readsAtTargetLookup = {} #dict for lookup
readsAtTargetToDiscard = {}
print("fetching " + args.swapChr + ":" + str(args.swapLoc))
for read in sourceAln.fetch(args.swapChr, args.swapLoc,args.swapLoc + 1):
	if read.query_name in readsAtTargetLookup:
		continue
	#include soft-clipped bases in these rStart and rEnd locations
	rStart = read.reference_start - read.query_alignment_start
	rEnd = read.reference_end + (read.query_length - read.query_alignment_end)
	rSeq = read.query_sequence
	rQual = read.query_qualities

	allReadAtTargetCount += 1
	readIsGood = 0


	if rEnd > args.swapLoc + olapBuffer or rStart < args.swapLoc - olapBuffer:
		readsAtTargetToDiscard[read.query_name + str(read.is_read1)] = 1
		continue
	readInd = 0
#	print str(read.cigar)
#	print 'read starts' + str(rStart) + ' and ends '+ str(rEnd)
	for op,count in read.cigar:
#		print "op: '"+str(op)+"' count: '"+str(count)+"'"
		#if match, add old sequence to new read
		if op == MATCH:
#			tStart = rStart + readInd
#			tEnd = rStart + readInd + count
#			print "from " + str(tStart) + " to " + str(tEnd) + "are match";

			if rStart + readInd < args.swapLoc and rStart + readInd + count > args.swapLoc :
				readsAtTarget.append(read)
				readsAtTargetLookup[read.query_name] = True
				goodReadAtTargetCount += 1
				readIsGood = 1
			readInd += count
		elif op == INSERTION:
#			tStart = rStart + readInd
#			tEnd = rStart + readInd + count
#			print "from " + str(tStart) + " to " + str(tEnd) + "are insertion"
			readInd += count
		elif op == DELETION:
#			tStart = rStart + readInd
#			tEnd = rStart + readInd + count
#			print "from " + str(tStart) + " to " + str(tEnd) + "are deletion"
			next
		elif op == CLIPPING:
#			tStart = rStart + readInd
#			tEnd = rStart + readInd + count
#			print "from " + str(tStart) + " to " + str(tEnd) + "are clipping"
			readInd += count
		else:
			raise Exception("got unrecognized op '"+str(op)+"' in cigar string " + str(read.cigar))

		if (not readIsGood):
			readsAtTargetToDiscard[read.query_name + str(read.is_read1)] = 1


print("read " + str(allReadAtTargetCount) + " reads, kept " + str(goodReadAtTargetCount) + " reads at target")
sourceAln.close()


#####
#next, randomize the list, and make edits to the first X reads
readsToReplace = {} #dict holds reads that were changed
readsAtTargetToPrint = {} #dict holds reads to print
random.shuffle(readsAtTarget)
numChanges = int(round(args.downsampleNumber * args.swapFreq))
numToPrint = args.downsampleNumber

lastCigarRE = re.compile(".*\D(\d+)([MIDNSHPX=])$")
lastCigarMatchRE = re.compile(".*\D(\d+)M$")

print("making changes in " + str(numChanges) + "/" + str(len(readsAtTarget)) + " sites")
for i in range(numChanges):
#	print 'replacement ' + str(i)
	read = readsAtTarget[i]
	readsToReplace[read.query_name+str(read.is_read1)] = read
	readsAtTargetToPrint[read.query_name] = True

for i in range(numChanges,numToPrint):
	read = readsAtTarget[i]
	readsAtTargetToPrint[read.query_name] = True



####
#next, read reads from bam that should be swapped in
alteredAln = pysam.AlignmentFile(args.alteredBam,"rb")
alteredRead1s = []
alteredRead2s = []
for read in alteredAln.fetch(args.swapChr, args.swapLoc,args.swapLoc + 1):
	if read.cigarstring is None:
		continue
	if args.onlyIncludeAlteredWithIndel is True and not ("I" in read.cigarstring or "D" in read.cigarstring):
		continue
	if read.is_read1:
		alteredRead1s.append(read)
	else:
		alteredRead2s.append(read)
random.shuffle(alteredRead1s)
random.shuffle(alteredRead2s)
alteredAln.close()


####
#finally, iterate over all the reads. Replace reads if they are in readToReplace
#reset the qualities -- when the sequence is changed, the quality is also changed.
readReads = 0
printedNotAtCutSiteReads = 0
printedAsIsReads = 0
printedChangedReads = 0
printedToControlReads = 0

#write mutations to file
mutatedReadsFile = open(mutatedReadsFileName,"w")
def addMutations(oldRead,newRead,alteredRead,mutationDict):
	oldGenomeLoc = oldRead.reference_start - oldRead.query_alignment_start

	#first, count which mutations were in the old one
	inOld = {}
	for (cigar_op,cigar_len) in oldRead.cigartuples:
		if cigar_op == MATCH:
			oldGenomeLoc += cigar_len
		elif cigar_op == INSERTION:
			mutKey = "%d\tI\t%d"%(cigar_len,oldGenomeLoc)
			if mutKey not in mutationDict:
				mutationDict[mutKey] = 0
			inOld[mutKey] += 1

		elif cigar_op == DELETION:
			mutKey = "%d\tD\t%d"%(cigar_len,oldGenomeLoc)
			if mutKey not in mutationDict:
				mutationDict[mutKey] = 0
			inOld[mutKey] += 1
			oldGenomeLoc += cigar_len
		elif cigar_op == CLIPPING:
			oldGenomeLoc += cigar_len
		else:
			print("got op: " + str(cigar_op) + " len: " + str(cigar_len) + " from " + str(oldRead.cigarstring))
			sys.exit()

	#then add mutations in altered to mutationDict if they weren't in the old
	genomeLoc = alteredRead.reference_start - alteredRead.query_alignment_start
	for (cigar_op,cigar_len) in alteredRead.cigartuples:
		if cigar_op == MATCH:
			genomeLoc += cigar_len
		elif cigar_op == INSERTION:
			mutKey = "%d\tI\t%d"%(cigar_len,genomeLoc)
			if mutKey not in inOld: #if mutation isn't in old read
				if mutKey not in mutationDict:
					mutationDict[mutKey] = 0
				mutationDict[mutKey] += 1

		elif cigar_op == DELETION:
			mutKey = "%d\tD\t%d"%(cigar_len,genomeLoc)
			if mutKey not in inOld:
				if mutKey not in mutationDict:
					mutationDict[mutKey] = 0
				mutationDict[mutKey] += 1
			genomeLoc += cigar_len
		elif cigar_op == CLIPPING:
			genomeLoc += cigar_len
		else:
			print("got op: " + str(cigar_op) + " len: " + str(cigar_len) + " from " + str(alteredRead.cigarstring))
			sys.exit()
def printReadReplacement(oldRead,newRead,alteredRead):
	minStart = oldRead.reference_start
	if newRead.reference_start < minStart:
		minStart = newRead.reference_start
	if alteredRead.reference_start < minStart:
		minStart = alteredRead.reference_start
	oldSpaces = " "*(oldRead.reference_start - minStart)
	newSpaces = " "*(newRead.reference_start - minStart)
	alteredSpaces =" "* (alteredRead.reference_start - minStart)
	thisRefSeq = refFile.fetch(args.swapChr,minStart,minStart + 300)
	mutatedReadsFile.write(
		"ref: %10s %10s %20s %s"%(minStart,minStart+300," ",thisRefSeq) +
		"\nold: %10s %10s %20s%s %s %s"%(oldRead.reference_start,oldRead.reference_end,oldRead.cigarstring, oldSpaces,oldRead.query_sequence,oldRead.query_name) +
		"\nalt: %10s %10s %20s%s %s %s"%(alteredRead.reference_start,alteredRead.reference_end,alteredRead.cigarstring, alteredSpaces,alteredRead.query_sequence,alteredRead.query_name) +
		"\nnew: %10s %10s %20s%s %s"%(newRead.reference_start,newRead.reference_end,newRead.cigarstring, newSpaces,newRead.query_sequence) +
		"\n\n")

#we will be printing read pairs, so only print downsamplePct/2
downsamplePct = (float(args.downsampleNumber)/float(allReadAtTargetCount))
sourceAlnNamesorted = pysam.AlignmentFile(args.unalteredNamesortedBam,"rb")
mutationDict = {} # holds counts of all mutations
for read1 in sourceAlnNamesorted.fetch(until_eof=True):
	read1String = read1.query_name + str(read1.is_read1)
	#read2 = sourceAlnNamesorted.next()
	read2 = next(sourceAlnNamesorted)
	read2String = read2.query_name + str(read2.is_read1)

	if (read1.query_name != read2.query_name):
		raise Exception("r1: " + read1.query_name + " is not " + read2.query_name)


	readReads += 2

	read1.query_qualities = [min(x+args.qualAdd,41) for x in read1.query_qualities]
	read2.query_qualities = [min(x+args.qualAdd,41) for x in read2.query_qualities]

	#write to control
	if(random.random() < downsamplePct):
		destAlnCtl.write(read1)
		destAlnCtl.write(read2)
		printedToControlReads += 2


	#if this read pair is at the target location,
	#read1.query_name is the same as read2.query_name
	if read1.query_name in readsAtTargetLookup:
		if read1.query_name in readsAtTargetToPrint:
			if (read1String in readsToReplace and read2String in readsToReplace):
				raise Exception("Both read1 and read2 were supposed to be swapped!")
			elif (read1String in readsToReplace):
				printedChangedReads += 1
				printedNotAtCutSiteReads += 1
				alteredRead1 = alteredRead1s.pop(0)
				origRead1 = copy.deepcopy(read1)
				newread1 = replaceRead(read1,alteredRead1,refSeq,refStart)
				newread1.query_qualities = [min(x+args.qualAdd,41) for x in newread1.query_qualities]
				printReadReplacement(origRead1,newread1,alteredRead1)
				addMutations(origRead1,newread1,alteredRead1,mutationDict)
				destAln.write(newread1)
				destAln.write(read2)
			elif (read2String in readsToReplace):
				printedChangedReads += 1
				printedNotAtCutSiteReads += 1
				alteredRead2 = alteredRead2s.pop(0)
				origRead2 = copy.deepcopy(read2)
				newread2 = replaceRead(read2,alteredRead2,refSeq,refStart)
				newread2.query_qualities = [min(x+args.qualAdd,41) for x in newread2.query_qualities]
				printReadReplacement(origRead2,newread2,alteredRead2)
				addMutations(origRead2,newread2,alteredRead2,mutationDict)
				destAln.write(newread2)
				destAln.write(read1)
			else:
				printedAsIsReads += 1
				printedNotAtCutSiteReads += 1
				destAln.write(read1)
				destAln.write(read2)
		continue

	#if the read is not at the target location, print it out with downsampling
	else:
		if(random.random() >= downsamplePct):
			continue
		else:
			destAln.write(read1)
			destAln.write(read2)
			printedNotAtCutSiteReads += 2

destAln.close()
destAlnCtl.close()
sourceAlnNamesorted.close()

mutationsFile = open(mutationsFileName,"w")
mutationsFile.write("BP\tINDEL\tLOC\tCOUNT\n")
for indelKey in mutationDict:
	mutationsFile.write(indelKey + "\t" + str(mutationDict[indelKey]) + "\n")
mutationsFile.close()

totPrinted = printedNotAtCutSiteReads + printedAsIsReads + printedChangedReads
outstring = "Finished\n";
outstring += "read "+ str(readReads) + " reads\n"
outstring += "printed " + str(printedNotAtCutSiteReads) + " reads not at the cut site (downsample pct was " + str(downsamplePct) + ")\n"
outstring += "printed " + str(printedAsIsReads) + " reads at cut site without modification\n"
outstring += "printed " + str(printedChangedReads) + " reads at cut site with modification\n"
outstring += "printed " + str(totPrinted) + " reads to the treatment bam\n"
outstring += "printed " + str(printedToControlReads) + " reads to the control bam\n"

print("sorting...")
os.system("samtools sort -o " + args.outfile + " " + unsortedOutName)
os.system("samtools index " + args.outfile)
os.system("rm " + unsortedOutName)

print("sorting control...")
os.system("samtools sort -o " + args.outfile + ".ctl.bam " + unsortedCtlOutName)
os.system("samtools index " + args.outfile + ".ctl.bam")
os.system("rm " + unsortedCtlOutName)

print (outstring)
f1 = open(args.outfile+".info","w")
f1.write(outstring)
f1.close()
