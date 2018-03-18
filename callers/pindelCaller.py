from caller import Caller
import numpy as np
import os
import stat
import configparser
import sys
import re

class PindelCaller(Caller):

	def __init__(self, Config):
		swapChr = Config.get('Simulation','chr') #location of simulated indels
		swapLoc = Config.getint('Simulation','loc') #location of simulated indels
		fdrTolerances = [int(i) for i in Config.get('Simulation','FDRtolerance').split(",")]
		pindelRange = max(fdrTolerances) + 100 #pindel will look for indels within this range around the swap loc

		#get mean insert size
		unalteredBam = Config.get('Simulation','unalteredBam')
		ub = open(unalteredBam+".insertSizeMetrics","r")
		line = ub.readline()
		mean_insert_size = -1
		while line:
			line = line.strip()
			lineEls = line.split()
			if line != "" and lineEls[0] == "MEDIAN_INSERT_SIZE" and lineEls[4] == "MEAN_INSERT_SIZE":
				dataline = ub.readline()
				dataEls = dataline.split()
				mean_insert_size = float(dataEls[4])
				break
			line = ub.readline()
		ub.close()

		if mean_insert_size < 0:
			raise Exception ("Could not parse mean insert size from " + unalteredBam + ".insertSizeMetrics")

		self.call_chr = swapChr
		self.call_start = swapLoc - pindelRange
		self.call_end = swapLoc + pindelRange

		self.reference = Config.get('Simulation','reference')
		self.mean_insert_size = mean_insert_size

	def get_name(self):
		return('Pindel')

	def run_caller(self,sample_bam,control_bam):
		outfileName = re.sub(r".bam$","",sample_bam)
		pindelConfigFile = outfileName + ".pindel.config"
		pFile = open(pindelConfigFile,"w")
		#pFile.write("/broad/hptmp/kendell/circle/simulations/alignNA12878/1kbAroundCutsite.BWA.bam " + str(meanInsertSize) + " NORMAL\n")
		pFile.write(control_bam + " " + str(self.mean_insert_size) + " NORMAL\n")
		pFile.write(sample_bam + " " + str(self.mean_insert_size) + " TUMOR\n")
		pFile.close()


#		r/--report_inversions           report inversions (default true)
#		-t/--report_duplications         report tandem duplications (default true)
#		-l/--report_long_insertions      report insertions of which the full sequence cannot be deduced because of their length (default true)
#		-k/--report_breakpoints          report breakpoints (default true)
#		-s/--report_close_mapped_reads   report reads of which only one end (the one closest to the mapped read of the paired-end read) could
#						    be mapped (default false)
		ignoreFlags = "-r false -t false -l false -k false -s false"
		#command = "pindel -f /data/pinello/COMMON_DATA/REFERENCE_GENOMES/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa " + \
		finished_file = outfileName + ".pindel.finished"
		command = "pindel -f " + self.reference + " " + \
				ignoreFlags + " -i " + pindelConfigFile + " -o " + outfileName + " -c " + self.call_chr + ":" + str(self.call_start) + "-" + str(self.call_end) + " && touch " + finished_file

		if os.path.isfile(finished_file):
			return ""
		else:
			return command

	def get_results(self, sample_bam,control_bam):
#From pindel documentation
#					The header line contains the following data:
#1) The index of the indel/SV (57 means that 57 insertions precede this insertion in the file)
#2) The type of indel/SV: I for insertion, D for deletion, INV for inversion, TD for tandem duplication
#3) The length of the SV
#4) "NT" (to indicate that the next number is the length of non-template sequences inserted; insertions are fully covered by the NT-fields, deletions can have NT bases if the deletion is not ‘pure’, meaning that while bases have been deleted, some bases have been inserted between the breakpoints)
#5) the length(s) of the NT fragment(s)
#6) the sequence(s) of the NT fragment(s)
#7-8) the identifier of the chromosome the read was found on
#9-10-11) BP: the start and end positions of the SV
#12-13-14) BP_range: if the exact position of the SV is unclear since bases at the edge of one read-half could equally well be appended to the other read-half. In the deletion example, ACA could be on any side of the gap, so the original deletion could have been between 1337143 and 1338815, between 1337144 and 1338816, or between 1337145 and 133817, or between 1337146 and 133818. BP-range is used to indicate this range.
#15) "Supports": announces that the total count of reads supporting the SV follow.
#16) The number of reads supporting the SV
#17) The number of unique reads supporting the SV (so not counting duplicate reads)
#18) +: supports from reads whose anchors are upstream of the SV
#19-20) total number of supporting reads and unique number of supporting reads whose anchors are upstream of the SV.
#21) -: supports from reads whose anchors are downstream of the SV
#22-23) total number of supporting reads and unique number of supporting reads whose anchors are downstream of the SV
#24-25) S1: a simple score, (“# +” + 1)* (“# -” + 1) ;
#26-27) SUM_MS: sum of mapping qualities of anchor reads, The reads with variants or unmapped are called split-read, whose mate is called anchor reads. We use anchor reads to narrow down the search space to speed up and increase sensitivity;
#28) the number of different samples scanned
#29-30-31) NumSupSamples?: the number of samples supporting the SV, as well as the number of samples having unique reads supporting the SV (in practice, these numbers are the same)
#32+) Per sample: the sample name, followed by the total number of supporting reads whose anchors are upstream, the total number of unique supporting reads whose anchors are upstream, the total number of supporting reads whose anchors are downstream, and finally the total number of unique supporting reads whose anchors are downstream.

				theseFoundMuts = {}
				outfileName = re.sub(r".bam$","",sample_bam)
				pindelConfigFile = outfileName + ".pindel.config"
				pindelFileNames = [outfileName + "_D", #deletiongs
						outfileName + "_SI", #short insertions
						outfileName + "_LI", #long insertions
						]
				for pindelFileName in pindelFileNames:
					pindelFile = open(pindelFileName,"r")
					for line in pindelFile:
						lineEls = line.split()
						if len(lineEls) > 20:
							if lineEls[25] != "SUM_MS":
								raise Exception("Can't read line " + line + " (got lineEls[25] = " + lineEls[25] + ")")

							#BP INDEL LOC => COUNT
							thisKey = lineEls[2] + " " + lineEls[1] + " " + lineEls[9]
							theseFoundMuts[thisKey] = int(lineEls[15]) #pindel Support

					pindelFile.close()
				return theseFoundMuts
