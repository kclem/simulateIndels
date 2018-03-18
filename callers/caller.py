class Caller:
	def __init__(self):
		""""Instantiate a new caller object.
		Any variables unique to this caller should be set here.
		For example, if this caller requires the path to the fastq reference, that should be set here.
		"""
		pass
	def get_name(self):
		"""Get the name of the caller.
		This name will be used in output that compares this caller to other callers
		Returns: A string with a name identifying the caller
		"""
		pass
	def run_caller(self, sample_bam,control_bam):
		"""Extract indels from a pair of bams.
		sample_bam is the bam with simulated indels
		control_bam is the bam with no simulated indels
		Returns: a string with a command to run this indel for this pair of bams
		   or: '' if the calling for the pair of bams has been completed
		"""
		pass
	def get_results(self, sample_bam,control_bam):
		"""Extract indel calling results.
		sample_bam is the bam with simulated indels
		control_bam is the bam with no simulated indels
		run_caller has been run previously. This function collects the results for indel calling by this caller
		Returns: a dictionary: 'BP INDEL LOC' => count
		where INDEL is either 'I' for insertion or 'D' for deletion
		e.g. if a 25bp insertion was observed at position 100 supported by 5 reads, the key would be '25 I 100' and the value would be 5
		"""
		pass

