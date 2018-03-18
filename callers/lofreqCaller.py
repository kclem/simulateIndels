from caller import Caller
import configparser
import os
import gzip
import re
import io

class LofreqCaller(Caller):
	def __init__(self,Config):
		self.reference = Config.get('Simulation','reference_wg')
		self.lofreqLoc = "~/software/lofreq/lofreq_star-2.1.3.1/bin/lofreq"

	def get_name(self):
		return('Lofreq')

	def run_caller(self, sample_bam,control_bam):
		outfileName = re.sub(r".bam$","",sample_bam)
		lofreqBam = sample_bam+".lofreq.bam"
		lofreqCtl = control_bam+".lofreq.bam"
		finishedFile = outfileName + ".lofreq.finished"

		command = "rm -f " + lofreqBam + "&& " + self.lofreqLoc + " indelqual --dindel -f " + self.reference + " -o " + lofreqBam + " " + sample_bam + " && samtools index " + lofreqBam + " && "
		command += "rm -f " + lofreqCtl + " && " + self.lofreqLoc + " indelqual --dindel -f " + self.reference + " -o " + lofreqCtl + " " + control_bam + " && samtools index " + lofreqCtl + " && "
		command += "rm -f " + outfileName + "_lofreq_* && " + self.lofreqLoc + " somatic -f " + self.reference + " " + \
			"-t " + lofreqBam + " -n " + lofreqCtl + " --min-cov 0 --call-indels -o " + outfileName + "_lofreq_ " + \
		" && touch " + finishedFile

		if os.path.isfile(finishedFile):
			return ""
		else:
			return command

	def get_results(self, sample_bam,control_bam):
		lofreqFileName = re.sub(r".bam$","_lofreq_somatic_final.indels.vcf.gz",sample_bam)
		lofreqFile = io.TextIOWrapper(gzip.open(lofreqFileName,'rb'))
		theseFoundMuts = {}
		for line in lofreqFile:
			if line.startswith("#"):
				continue
			lineEls = line.split()
			numChange = len(lineEls[4]) - len(lineEls[3])
			dp = re.search(r'DP=(\d+);',lineEls[7]).group(1)
			af = re.search(r'AF=([\d\.]+);',lineEls[7]).group(1)
			count = float(dp)*float(af)
			if numChange > 0:
				#BP INDEL LOC => COUNT
				thisKey = str(numChange) + " " + "I" + " " + lineEls[1]
				theseFoundMuts[thisKey] = int(count)
			elif numChange < 0:
				#BP INDEL LOC => COUNT
				thisKey = str(-1*numChange) + " " + "D" + " " + lineEls[1]
				theseFoundMuts[thisKey] = int(count) 

		lofreqFile.close()
		return theseFoundMuts
