from caller import Caller
import vcf
import numpy as np
import os
import stat
import configparser
import sys
import re

class VarscanCaller(Caller):

	def __init__(self, Config):

		self.reference = Config.get('Simulation','reference')

	def get_name(self):
		return('Varscan')

	def run_caller(self,sample_bam,control_bam):
		outfileRoot = re.sub(r".bam$",".varscan",sample_bam)
		finishedFile = outfileRoot + ".finished"

		command = "samtools mpileup -f " + self.reference + " " + sample_bam + " " + control_bam + " -o " + outfileRoot + ".mpileup.vcf " + \
		" && varscan somatic  " + outfileRoot + ".mpileup.vcf " + outfileRoot + " --mpileup 1 --output-vcf --min-coverage 0 --min-coverage-normal 0 --min-coverage-tumor 0 --somatic-p-value 1 " + \
		" && touch " + finishedFile

		if os.path.isfile(finishedFile):
			return ""
		else:
			return command

	def get_results(self, sample_bam,control_bam):
		outfileRoot = re.sub(r".bam$",".varscan",sample_bam)
		finishedFile = outfileRoot + ".finished"
		vcfFileName = outfileRoot + ".indel.vcf"

		vcf_reader = vcf.Reader(open(vcfFileName,'r'))
		theseFoundMuts = {}
		for record in vcf_reader:
			lenChange = len(record.ALT[0]) - len(record.REF)
			numControl = record.samples[0]['AD']
			numAlt = record.samples[1]['AD']
			
			#if tumor contains more alternat than WT, swap
			if record.samples[1]['AD'] > record.samples[1]['RD']:
				numControl = record.samples[0]['RD']
				numAlt = record.samples[1]['RD']
				lenChange = -1*lenChange


			if lenChange > 0:
				#BP INDEL LOC => COUNT
				thisKey = str(lenChange) + " " + "I" + " " + str(record.POS)
				theseFoundMuts[thisKey] = numAlt
			elif lenChange < 0:
				#BP INDEL LOC => COUNT
				thisKey = str(-1*lenChange) + " " + "D" + " " + str(record.POS)
				theseFoundMuts[thisKey] = numAlt

		return theseFoundMuts
