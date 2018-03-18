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
import pprint

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



def replaceRead(_oldRead,_newRead,_genome_seq,_genome_start):
	oldRead = copy.deepcopy(_oldRead)
	newRead = copy.deepcopy(_newRead)
	originalLen = oldRead.query_length

#	print("old: " + str(oldRead))
#	print("new: " + str(newRead))

	#old is NA (unaltered read)
	oldStart = oldRead.reference_start - oldRead.query_alignment_start
	oldEnd = oldRead.reference_end + (oldRead.query_length - oldRead.query_alignment_end)

	#new is EMX1 (altered read)
	newStart = newRead.reference_start - newRead.query_alignment_start
	newEnd = newRead.reference_end + (newRead.query_length - newRead.query_alignment_end)

#	print ("oldStart: " + str(oldStart) + " oldEnd " + str(oldEnd))
#	print ("newStart: " + str(newStart) + " newEnd " + str(newEnd))


	explodeOldCigar = explodeCigar(str(oldRead.cigarstring))
	explodeNewCigar = explodeCigar(str(newRead.cigarstring))

	final_read = ""
	final_cigar = ""
	#if NA read starts to the left of the EMX1 read, just take bases from the NA read

	addNewLoc = 0 #index of new read that we are adding

#	print ("old: " + oldRead.query_sequence)
#	print ("old: " + oldRead.cigarstring)
#	print ("old: " + explodeOldCigar)
#	print ("old: " + str(oldRead.query_qualities))
#	print ("new: " + newRead.query_sequence)
#	print ("new: " + newRead.cigarstring)
#	print ("new: " + explodeNewCigar)
#	print ("new: " + str(newRead.query_qualities))


	final_genomic_loc = oldStart
	oldReadLoc = 0
	oldCigarLoc = 0
	if oldStart < newStart:
		while final_genomic_loc < newStart:
			if str(explodeOldCigar[oldCigarLoc]) == "H":
				final_cigar += explodeOldCigar[oldCigarLoc]
				oldCigarLoc += 1
			elif str(explodeOldCigar[oldCigarLoc]) == "D":
				final_cigar += explodeOldCigar[oldCigarLoc]
				oldCigarLoc += 1
				final_genomic_loc += 1
			elif str(explodeOldCigar[oldCigarLoc]) == "I":
				final_cigar += explodeOldCigar[oldCigarLoc]
				final_read += oldRead.query_sequence[oldReadLoc]
				oldCigarLoc += 1
				oldReadLoc += 1
			else:
				final_cigar += explodeOldCigar[oldCigarLoc]
				final_read += oldRead.query_sequence[oldReadLoc]
				oldCigarLoc += 1
				oldReadLoc += 1
				final_genomic_loc += 1

#	print ("len is : " + str(len(oldRead.query_sequence)) + " first is: " + final_read + " oldStart " + str(oldStart) + " oldEnd " + str(oldEnd) + " new start: " + str(newStart) + " newEnd: " + str(newEnd))

	newReadLoc = 0
	newCigarLoc = 0
	#if old read starts in middle of new, catch the 'readLoc' pointer up to where we start pulling from the newRead sequence
	if oldStart > newStart:
		final_genomic_loc = newStart
		while final_genomic_loc < oldStart:
			if str(explodeNewCigar[newCigarLoc]) == "H":
				newCigarLoc += 1
			elif str(explodeNewCigar[newCigarLoc]) == "D":
				newCigarLoc += 1
				final_genomic_loc += 1
			elif str(explodeNewCigar[newCigarLoc]) == "I":
				newCigarLoc += 1
				newReadLoc += 1
			else:
				newCigarLoc += 1
				newReadLoc += 1
				final_genomic_loc += 1

#	print ("checking if first cigar is deletion (" + explodeNewCigar[newCigarLoc]+")")
#	print ("finalCigar: " + final_cigar)
#	print ("finalRead: " + final_read)

#	print ("newSeq: " + newRead.query_sequence)
#	print ("explodenewcigar: " + explodeNewCigar)
	#if first cigar is a deletion, add some bases before
	if explodeNewCigar[newCigarLoc] == "D" and len(final_read) == 0:
		#first determine how long the deletion string was
		newCigarDelStart = newCigarLoc
		newCigarDelEnd = newCigarLoc
		while newCigarDelEnd < len(explodeNewCigar):
			if explodeNewCigar[newCigarDelEnd] != "D":
				break
			newCigarDelEnd += 1

		#find start of del
		while newCigarDelStart > 0:
			if explodeNewCigar[newCigarDelStart] != "D":
				break
			newCigarDelStart -= 1

		deletionLength = (newCigarDelEnd-1) - newCigarDelStart
		tempNewCigarLoc = newCigarDelStart -1 #start right before deletion starts
		tempNewReadLoc = newReadLoc - 1 #because there weren't any bases(they were all deletions) start is the same

#		print("deletion length: " + str(deletionLength))



		#build up final_read and final_cigar from right to left, then reverse later
		addedBP = 0
		genomeShift = 0
		while tempNewCigarLoc > 0 and addedBP < deletionLength:
			if str(explodeNewCigar[tempNewCigarLoc]) == "H":
				final_cigar += explodeNewCigar[tempNewCigarLoc]
				tempNewCigarLoc -= 1
			elif str(explodeNewCigar[tempNewCigarLoc]) == "D":
				final_cigar += explodeNewCigar[tempNewCigarLoc]
				tempNewCigarLoc -= 1
				genomeShift += 1
			elif str(explodeNewCigar[tempNewCigarLoc]) == "I":
				final_cigar += explodeNewCigar[tempNewCigarLoc]
				final_read += newRead.query_sequence[tempNewReadLoc]
				tempNewCigarLoc -= 1
				tempNewReadLoc -= 1
				addedBP += 1
			else:
#				print("adding qseq[" + str(tempNewReadLoc) + "]" + newRead.query_sequence[tempNewReadLoc])
				final_cigar += explodeNewCigar[tempNewCigarLoc]
				final_read += newRead.query_sequence[tempNewReadLoc]
				tempNewCigarLoc -= 1
				tempNewReadLoc -= 1
				genomeShift += 1
				addedBP += 1

		#update start of read
		oldRead.reference_start = final_genomic_loc-genomeShift
		final_read = final_read[::-1] #reverse because we built them up from left to right
		final_cigar = final_cigar[::-1]

#	print ("finalCigar: " + final_cigar)
#	print ("finalRead: " + final_read)


#	print("newReadLoc: " + str(newReadLoc) + " newCigarLoc " + str(newCigarLoc) + " final_genomic_loc: " + str(final_genomic_loc))
#	print("at 1: newReadLoc: %s newCigarLoc: %s\n" % (newReadLoc,newCigarLoc)+\
#		"oldReadLoc: %s oldCigarLoc: %s\n" % (oldReadLoc,oldCigarLoc)+\
#		"final_genomic_loc: %s final_read: %s finalCigar: %s\n" % (final_genomic_loc,final_read,final_cigar))
	#overalpping locations
	while len(final_read) < len(oldRead.query_sequence):
		if newReadLoc >= len(newRead.query_sequence):
			break
		if str(explodeNewCigar[newCigarLoc]) == "H":
			final_cigar += explodeNewCigar[newCigarLoc]
			newCigarLoc += 1
		elif str(explodeNewCigar[newCigarLoc]) == "D":
			final_cigar += explodeNewCigar[newCigarLoc]
			newCigarLoc += 1
			final_genomic_loc += 1
		elif str(explodeNewCigar[newCigarLoc]) == "I":
			final_cigar += explodeNewCigar[newCigarLoc]
			final_read += newRead.query_sequence[newReadLoc]
			newCigarLoc += 1
			newReadLoc += 1
		else:
			final_cigar += explodeNewCigar[newCigarLoc]
			final_read += newRead.query_sequence[newReadLoc]
			newCigarLoc += 1
			newReadLoc += 1
			final_genomic_loc += 1
#	print ("len is : " + str(len(oldRead.query_sequence)) + " mid is: " + final_read + " " + final_cigar + " oldStart " + str(oldStart) + " new start: " + str(newStart))
#	print("at 2: newReadLoc: %s newCigarLoc: %s\n" % (newReadLoc,newCigarLoc)+\
#		"oldReadLoc: %s oldCigarLoc: %s\n" % (oldReadLoc,oldCigarLoc)+\
#		"final_genomic_loc: %s final_read: %s finalCigar: %s\n" % (final_genomic_loc,final_read,final_cigar))

#	pprint.pprint(str(oldRead))
#	refPos = oldRead.get_reference_positions(full_length=True)
#	print("pos: " + str(refPos))
#	pprint.pprint(str(newRead))
#	refPos = newRead.get_reference_positions(full_length=True)
#	print("pos: " + str(refPos))
#	print ("len final read: " + str(len(final_read)) + " old len: " + str(len(oldRead.query_sequence)))
	if len(final_read) < len(oldRead.query_sequence):
		#catch old pointer up
		while oldStart + oldReadLoc < oldEnd and oldStart + oldReadLoc < final_genomic_loc:
			if str(explodeOldCigar[oldCigarLoc]) == "H":
				oldCigarLoc += 1
			elif str(explodeOldCigar[oldCigarLoc]) == "D":
				oldCigarLoc += 1
			elif str(explodeOldCigar[oldCigarLoc]) == "I":
				oldCigarLoc += 1
				oldReadLoc += 1
			else:
				oldCigarLoc += 1
				oldReadLoc += 1
#		print("at 3: newReadLoc: %s newCigarLoc: %s\n" % (newReadLoc,newCigarLoc)+\
#			"oldReadLoc: %s oldCigarLoc: %s\n" % (oldReadLoc,oldCigarLoc)+\
#			"final_genomic_loc: %s final_read: %s finalCigar: %s\n" % (final_genomic_loc,final_read,final_cigar))

		#add from old
		while final_genomic_loc < oldEnd and len(final_read) < len(oldRead.query_sequence):
			if str(explodeOldCigar[oldCigarLoc]) == "H":
				final_cigar += explodeOldCigar[oldCigarLoc]
				oldCigarLoc += 1
			elif str(explodeOldCigar[oldCigarLoc]) == "D":
				final_cigar += explodeOldCigar[oldCigarLoc]
				oldCigarLoc += 1
				final_genomic_loc += 1
			elif str(explodeOldCigar[oldCigarLoc]) == "I":
				final_cigar += explodeOldCigar[oldCigarLoc]
				final_read += oldRead.query_sequence[oldReadLoc]
				oldCigarLoc += 1
				oldReadLoc += 1
			else:
				final_cigar += explodeOldCigar[oldCigarLoc]
				final_read += oldRead.query_sequence[oldReadLoc]
				oldCigarLoc += 1
				oldReadLoc += 1
				final_genomic_loc += 1

#		print("at 4: newReadLoc: %s newCigarLoc: %s\n" % (newReadLoc,newCigarLoc)+\
#			"oldReadLoc: %s oldCigarLoc: %s\n" % (oldReadLoc,oldCigarLoc)+\
#			"final_genomic_loc: %s final_read: %s finalCigar: %s\n" % (final_genomic_loc,final_read,final_cigar))

		#finally, add from genome
		while len(final_read) < len(oldRead.query_sequence):
			final_cigar += "M"
			final_read += _genome_seq[final_genomic_loc - _genome_start]
			final_genomic_loc += 1


#	print("at 5: newReadLoc: %s newCigarLoc: %s\n" % (newReadLoc,newCigarLoc)+\
#		"oldReadLoc: %s oldCigarLoc: %s\n" % (oldReadLoc,oldCigarLoc)+\
#		"final_genomic_loc: %s final_read: %s finalCigar: %s\n" % (final_genomic_loc,final_read,final_cigar))

#	if re.search("14D",newRead.cigarstring):
#		print ("oldStart: " + str(oldStart) + " oldEnd " + str(oldEnd))
#		print ("newStart: " + str(newStart) + " newEnd " + str(newEnd))
#		print ("old: " + oldRead.query_sequence)
#		print ("old: " + oldRead.cigarstring)
#		print ("old: " + explodeOldCigar)
#		print ("old: " + str(oldRead.query_qualities))
#		print ("new: " + newRead.query_sequence)
#		print ("new: " + newRead.cigarstring)
#		print ("new: " + explodeNewCigar)
#		print ("new: " + str(newRead.query_qualities))
#
#		print ("final: " + final_read)
#		print ("final: " + final_cigar)




	oldQuals = oldRead.query_qualities
	oldRead.query_sequence = final_read
	oldRead.cigarstring = unexplodeCigar(final_cigar)
	oldRead.query_qualities = oldQuals

	return (oldRead)


if __name__ == "__main__":

	genome = "G"*100
	genome_start = 3

	a = pysam.AlignedSegment()
	qseq = "A"*20
	a.query_name = "read1"
	a.query_sequence=qseq
	a.flag = 0
	a.reference_id = 0
	a.reference_start = 10
	a.mapping_quality = 20
	a.cigarstring = str(len(qseq)) + "M"
	a.query_qualities = pysam.qualitystring_to_array("<"*len(qseq))
	a.tags = (("NM", 1),("RG", "L1"))

	b = pysam.AlignedSegment()
	qseq = "T"*20
	b.query_name = "read2"
	b.query_sequence=qseq
	b.flag = 0
	b.reference_id = 0
	b.reference_start = 15
	b.mapping_quality = 20
	b.cigarstring = str(len(qseq)) + "M"
	b.query_qualities = pysam.qualitystring_to_array("<"*len(qseq))
	b.tags = (("NM", 1),("RG", "L1"))

	print("readA: " + str(a))
	print("readB: " + str(b))
	# b is old
	if(0):
		new = replaceRead(b,a)
		if new.query_sequence != "AAAAAAAAAAAAAAATTTTT":
			raise Exception("Unexpected result!")


	#    10
	# a: AAAAAAAAAAAAAAAAAAA
	# b:      TTTTTTTTTTTTTTTTTTTT
	# m: AAAAAAAAAA-----AAAAACCCCC
	# m2:       AAA---CCC
	m = pysam.AlignedSegment()
	qseq = "AAAAAAAAAAAAAAACCCCC"
	m.query_name = "mod"
	m.query_sequence=qseq
	m.flag = 0
	m.reference_id = 0
	m.reference_start = 10
	m.mapping_quality = 20
	m.cigarstring = "10M5D10M"
	m.query_qualities = pysam.qualitystring_to_array("<"*len(qseq))
	m.tags = (("NM", 1),("RG", "L1"))

	m2 = pysam.AlignedSegment()
	qseq = "AAACCC"
	m2.query_name = "mod2"
	m2.query_sequence=qseq
	m2.flag = 0
	m2.reference_id = 0
	m2.reference_start = 18
	m2.mapping_quality = 20
	m2.cigarstring = "3M3D3M"
	m2.query_qualities = pysam.qualitystring_to_array("<"*len(qseq))
	m2.tags = (("NM", 1),("RG", "L1"))

	new = replaceRead(b,m2,genome,genome_start)
	print("new:" + str(new))
	print("readA: " + str(a))
	print("readB: " + str(b))
