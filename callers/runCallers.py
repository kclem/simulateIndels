from runReplaceReads import runReplaceReads
from runCallersHelpers import aggregate_indels, print_indels, print_aggregate_indels, print_all_aggregate_indels
from simCaller import SimCaller
from pindelCaller import PindelCaller
from lofreqCaller import LofreqCaller
from somaticSniperCaller import SomaticSniperCaller
from varscanCaller import VarscanCaller
import configparser
import os
import sys
from subprocess import call

sys.path.insert(0, '/PHShome/mkc50/2017_07_DARPA_SIMULATIONS/simulations/circle/replaceReads/scripts_2')


def runCallers():
	#parse settings file
	if len(sys.argv) < 2:
		raise Exception("Couldn't find settings file")
	settingsFile = sys.argv[1]
	if settingsFile == "":
		settingsFile = "simSettings.txt"
	if not os.path.isfile(settingsFile):
		raise Exception("Couldn't find settings file")
	workFolder = os.path.dirname(os.path.abspath(settingsFile)) + "/"
	outFolder = workFolder + "out/"
	if not os.path.isdir(outFolder):
		os.mkdir(outFolder)
	Config = configparser.ConfigParser()
	Config.read(settingsFile)
	depths = [int(i) for i in Config.get('Simulation','depths').split(",")]
	pctMut = [float(i) for i in Config.get('Simulation','pctMut').split(",")]
	addQual = [int(i) for i in Config.get('Simulation','addQual').split(",")]
	maxReps = Config.getint('Simulation','reps')
	reps = range(0,maxReps)

	#first simulate samples by replacing reads
	replaceReadsScript = "/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/simulations/circle/replaceReads/scripts/replaceReads.py"
	simulateReadsCommands = runReplaceReads(Config,workFolder,outFolder,replaceReadsScript)
	if len(simulateReadsCommands) > 0:
		print ("got " + str(len(simulateReadsCommands)) + " commands to simulate samples")
		for command in simulateReadsCommands:
			call(command,shell=True)
	else:
		print('Finished simulating samples')



	#store the names and bams of each sample
	samples = []
	for d in depths:
		for p in pctMut:
			for q in addQual:
				for r in reps:
					simulated_bam = outFolder + "replace_d" + str(d) + "_p" + str(p) + "_q" + str(q)+ "_r" + str(r)+ ".bam"
					simulated_control = outFolder + "replace_d" + str(d) + "_p" + str(p) + "_q" + str(q) + "_r" + str(r) + ".bam.ctl.bam"
					name = "replace_d" + str(d) + "_p" + str(p) + "_q" + str(q)+ "_r" + str(r)
					samples.append((name,simulated_bam,simulated_control))

	#next read simulated indels
	simCaller = SimCaller()
	simulated_indels = {}
	for name, simulated_bam, simulated_control in samples:
		these_simulated_indels = simCaller.get_results(simulated_bam,simulated_control)
		simulated_indels[name] = these_simulated_indels
	aggregated_simulated_indels = aggregate_indels(simulated_indels)

	callers = []

	#create pindel caller and append it to to the list of callers
	pindelCaller = PindelCaller(Config)
	callers.append(pindelCaller)
	lofreqCaller = LofreqCaller(Config)
	callers.append(lofreqCaller)
	somaticSniperCaller = SomaticSniperCaller(Config)
	callers.append(somaticSniperCaller)
	varscanCaller = VarscanCaller(Config)
	callers.append(varscanCaller)

	for caller in callers:
		caller_commands = []
		caller_name = caller.get_name()
		for name, simulated_bam, simulated_control in samples:
			this_command = caller.run_caller(simulated_bam,simulated_control)
			if this_command != "":
				caller_commands.append(this_command)
		if len(caller_commands) > 0:
			print ("got " + str(len(caller_commands)) + " commands to run " + caller_name)
			for command in caller_commands:
				print("running command " + command)
				call(command,shell=True)
				print("Finished");
		else:
			print('Finished running ' + caller_name + ' commands')


	# Read in pindel results
	all_indels = {}
	all_aggregated_indels = {}
	for caller in callers:
		caller_name = caller.get_name()
		caller_indels = {}
		for name, simulated_bam, simulated_control in samples:
			these_simulated_indels = caller.get_results(simulated_bam,simulated_control)
			caller_indels[name] = these_simulated_indels
		all_indels[caller_name] = caller_indels
		caller_aggregated_indels = aggregate_indels(caller_indels)
		all_aggregated_indels[caller_name] = caller_aggregated_indels
		print("Finished reading results from calling indels with " + caller_name)

	print("Finished")


if __name__ == "__main__":
	runCallers()
