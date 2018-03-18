from caller import Caller

class SimCaller(Caller):

	def __init__(self):
		pass

	def run_caller(self,sample_bam,control_bam):
		return []

	def get_results(self, sample_bam,control_bam):
		mutsIn = sample_bam + ".mutations.txt"
		mutsFile = open(mutsIn,"r")
		theseSimulatedMuts = {}
		headLine = mutsFile.readline()
		#BP INDEL LOC COUNT
		for line in mutsFile:
			line.strip()
			lineEls = line.split()
			indelKey = lineEls[0] + " " + lineEls[1] + " " + lineEls[2]
			theseSimulatedMuts[indelKey] = int(lineEls[3])
		mutsFile.close()

		return theseSimulatedMuts
