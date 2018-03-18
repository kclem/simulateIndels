from subprocess import call
import numpy as np
import os
import stat
import configparser
import sys

def runReplaceReads(Config,workFolder,outFolder,replaceReadsScript):

    depths = [int(i) for i in Config.get('Simulation','depths').split(",")]
    pctMut = [float(i) for i in Config.get('Simulation','pctMut').split(",")]
    addQual = [int(i) for i in Config.get('Simulation','addQual').split(",")]
    maxReps = Config.getint('Simulation','reps')
    swapChr = Config.get('Simulation','chr')
    swapLoc = Config.getint('Simulation','loc')
    useOnlyIndelFlag = Config.getboolean('Simulation','useOnlyIndels')
    reference = Config.get('Simulation','reference')
    sourceBam = Config.get('Simulation','sourceBam')
    alteredBam = Config.get('Simulation','alteredBam')
    unalteredBam = Config.get('Simulation','unalteredBam')
    unalteredNamesortedBam = Config.get('Simulation','unalteredNamesortedBam')
    simulateRange = Config.getint('Simulation','simulateRange') # bp around cut site to simulate reads for
    
    reps = range(0,maxReps)
    if useOnlyIndelFlag: useOnlyIndelString = " --onlyIncludeAlteredWithIndel "
    
    
    if not os.path.isfile(alteredBam):
        raise Exception("Cannot find altered bam '" + alteredBam + "''")
    if not os.path.isfile(reference):
        raise Exception("Cannot find altered bam '" + reference + "''")
    
    processedBamsCreated = False
    if os.path.isfile(unalteredBam) and os.path.isfile(unalteredNamesortedBam):
        processedBamsCreated = True
    
    
    if not processedBamsCreated:
            filesToDelete = []
    	#check whether source bam exists (this is really big, so sometimes it is deleted)
            if not os.path.isfile(sourceBam):
                raise Exception("Cannot find source bam '" + sourceBam + "''")
            simStart = swapLoc - simulateRange
            simEnd = swapLoc + simulateRange
            tempNovalign = workFolder + str(simulateRange*2) + "bAroundCutsite.novoalign.bam"
    	#samtools view -Sb ../../NA12878/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam chr2:72928870-72938870 > 1kbAroundCutsite.novoalign.bam
            cmd = "samtools view -Sb " + sourceBam + " " + swapChr + ":" + str(simStart) + "-" + str(simEnd) + " > " + tempNovalign
            print("Calling: "+cmd)
            call(cmd,shell=True)
    
            filesToDelete.append(tempNovalign)
    
            tempSort = tempNovalign + ".qsort.bam"
    	#samtools sort -n 1kbAroundCutsite.novoalign.bam 1kbAroundCutsite.novoalign.bam.qsort
            cmd = "samtools sort -n " + tempNovalign + " -o " +  tempSort
            print("Calling: "+cmd)
            call(cmd,shell=True)
    
            filesToDelete.append(tempSort)
    
    	#bedtools bamtofastq -i 1kbAroundCutsite.novoalign.bam.qsort.bam -fq 1kbAroundCutsite.novoalign.bam.R1.fastq -fq2 1kbAroundCutsite.novoalign.bam.R2.fastq
            cmd = "module load BEDTools_2.17 && bedtools bamtofastq -i " + tempSort + " -fq " + tempSort + ".R1.fastq -fq2 " + tempSort + ".R2.fastq"
            print("Calling: "+cmd)
            call(cmd,shell=True)
    
            filesToDelete.append(tempSort+".R1.fastq")
            filesToDelete.append(tempSort+".R2.fastq")
    
    	#bwa mem /seq/references/Homo_sapiens_assembly38/v0/Homo_sapiens_assembly38.fasta 1kbAroundCutsite.novoalign.bam.R1.fastq 1kbAroundCutsite.novoalign.bam.R2.fastq > 1kbAroundCutsite.BWA.sam
            bwaSam = workFolder + str(simulateRange*2) + "bAroundCutsite.BWA.sam"
            cmd = "bwa mem " + reference + " " + tempSort + ".R1.fastq " + tempSort + ".R2.fastq > " + bwaSam
            print("Calling: "+cmd)
            call(cmd,shell=True)
            filesToDelete.append(bwaSam)
    
    	#filter out secondary alignments (-F 2048)
    	#samtools view -F 2048 -Sb 1kbAroundCutsite.BWA.sam | samtools sort - -o 1kbAroundCutsite.BWA.bam && samtools index 1kbAroundCutsite.BWA.bam
            cmd = "samtools view -F 2048 -Sb " + bwaSam + " | samtools sort - -o " + unalteredBam + " && samtools index " + unalteredBam
            print("Calling: "+cmd)
            call(cmd,shell=True)
    
    	#create name-sorted bam for operating on read pairs
    	#samtools view -F 2048 -Sb 1kbAroundCutsite.BWA.sam | samtools sort -n -o 1kbAroundCutsite.BWA.nameSort.bam
            cmd = "samtools view -F 2048 -Sb " + bwaSam + " | samtools sort -n -o " + unalteredNamesortedBam
            print("Calling: "+cmd)
            call(cmd,shell=True)
    
    	#java -Xmx8g -jar /seq/software/picard/current.20170216/bin/picard.jar CollectInsertSizeMetrics I=1kbAroundCutsite.BWA.bam O=1kbAroundCutsite.BWA.bam.insertSizeMetrics H=1kbAroundCutsite.BWA.bam.insertSizeMetrics.pdf
            cmd = "picard CollectInsertSizeMetrics I=" + unalteredBam + " O=" + unalteredBam + ".insertSizeMetrics H=" + unalteredBam + ".insertSizeMetrics.pdf"
            print("Calling: "+cmd)
            call(cmd,shell=True)
    
#        for file in filesToDelete:
#            call(["rm",file])
    
    
    qCommandArr = []
    for d in depths:
        for p in pctMut:
            for q in addQual:
                for r in reps:
                    outfileName = outFolder + "replace_d" + str(d) + "_p" + str(p) + "_q" + str(q) + "_r" + str(r) + ".bam"
                    command = "python " + replaceReadsScript + " --downsampleNumber " + str(d) + " --swapFreq " + str(p) + " --qualAdd " + str(q) + \
                        " --swapChr " + swapChr + " --swapLoc " + str(swapLoc) + " --reference " + reference + \
                        " --unalteredBam " + unalteredBam + " --unalteredNamesortedBam " + unalteredNamesortedBam + " --alteredBam " + alteredBam + \
                        " --outfile " + outfileName + " " + useOnlyIndelString + " ; touch " + outfileName + ".finished"
                    if (not os.path.isfile(outfileName+".finished") or not os.path.isfile(outfileName)):
                        qCommandArr.append(command)
    return qCommandArr


if __name__ == "__main__":

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
    print("Read setting file " + settingsFile)

    replaceReadsScript = "/data/pinello/PROJECTS/2017_07_DARPA_SIMULATIONS/simulations/circle/replaceReads/scripts/replaceReads.py"
    commands = runReplaceReads(Config,workFolder,outFolder,replaceReadsScript)
