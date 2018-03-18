
def aggregate_indels(indels):
	"""Aggregate counts of each indel type across simulated samples"""
	aggregated_indels = {}
	for name in indels:
		for key in indels[name]:
			if key not in aggregated_indels:
				aggregated_indels[key] = 0
			aggregated_indels[key] += indels[name][key]
	return aggregated_indels

def print_indels(indels,fileName=None):
	"""Prints counts of indels for a single caller across simulated samples
	fileName: If given, results are written to the file, otherwise they are printed to stdOut"
	"""
	out_string = "Simulated Sample\tBP\tINDEL\tLOC\tCOUNT\n";
	for name in sorted(indels):
		for key in sorted(indels[name]):
			out_string += name + "\t" + key.replace(" ", "\t")

	if fileName is not None:
		fh = open(fileName,"w")
		fh.write(out_string)
		fh.close
	else:
		print(out_string)

def print_aggregate_indels(indels,fileName=None):
	"""Prints counts of aggregated indels across simulated samples
	fileName: If given, results are written to the file, otherwise they are printed to stdOut"
	"""
	out_string = "BP\tINDEL\tLOC\tCOUNT\n";
	for key in sorted(indels):
		out_string += key.replace(" ", "\t")

	if fileName is not None:
		fh = open(fileName,"w")
		fh.write(out_string)
		fh.close
	else:
		print(out_string)

def print_all_aggregate_indels(all_aggregated_indels,aggregated_simulated_indels,callers,fileName=None):
	"""For all callers, prints counts of aggregated indels
	fileName: If given, results are written to the file, otherwise they are printed to stdOut"
	"""
	all_indels = {} #keep track of all possible keys
	for indel in aggregated_simulated_indels:
		if indel not in all_indels:
			all_indels[indel] = 0

	head_string = "BP\tINDEL\tLOC\tSIMULATED";
	caller_names = []
	for caller in callers:
		caller_name = caller.get_name()
		caller_names.append(caller_name)
		head_string += "\t"+caller_name
		for indel in all_aggregated_indels[caller_name]:
			if indel not in all_indels:
				all_indels[indel] = 0

	out_string = ""
	for indel in sorted(all_indels):
		out_string += indel.replace(" ", "\t")
		sim_count = 0
		if indel in aggregated_simulated_indels:
			sim_count = aggregated_simulated_indels[indel]
		out_string += "\t" + str(sim_count)

		for caller_name in caller_names:
			caller_count = 0
			if indel in all_aggregated_indels[caller_name]:
				caller_count = all_aggregated_indels[caller_name][indel]
			out_string += "\t" + str(caller_count)
		out_string += "\n"


	out_string = head_string + "\n" + out_string
	if fileName is not None:
		fh = open(fileName,"w")
		fh.write(out_string)
		fh.close
	else:
		print(out_string)


