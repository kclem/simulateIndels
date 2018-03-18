from caller import Caller
import numpy as np
import os
import stat
import configparser
import sys
import re

class SomaticSniperCaller(Caller):

	def __init__(self, Config):

		self.reference = Config.get('Simulation','reference')

	def get_name(self):
		return('SomaticSniper')

	def run_caller(self,sample_bam,control_bam):
		outfileName = re.sub(r".bam$",".somaticSniper.vcf",sample_bam)
		finishedFile = outfileName + ".finished"

		command = "bam-somaticsniper -F vcf -f " + self.reference + " " + sample_bam + " " + control_bam + " " + outfileName + \
		" && touch " + finishedFile

		if os.path.isfile(finishedFile):
			return ""
		else:
			return command

	def get_results(self, sample_bam,control_bam):

		vcfFileName = re.sub(r".bam$",".somaticSniper.vcf",sample_bam)
		vcfFile = open(vcfFileName,'r')
		theseFoundMuts = {}
		for line in vcfFile:
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

		vcfFile.close()
		return theseFoundMuts
