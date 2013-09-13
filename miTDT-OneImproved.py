import sys
import string
import gzip
import re
import time
import os
import numpy as np
from scipy import stats

#rTDT: gives out a range of chiSq TDT statistic due to admissible completions of missing trios' genotypes
#MITDT-ONE: uses marker allele frequencies over rTDT; computed as AC/AN for every ALT allele

#USAGE:  python miTDT-OneImproved.py mode=allelic model=all fm=/titan/cancerregulome14/inova_data/mergedVCF/var2vcf_anno_mie_filtered_merged/merged_sorted_filteredVCF_v4.vcf.gz ids=DF4/DF4_VCF_PretermPEDIDs #output=tdtResults log=tdtLog

#PARSE COMMAND LINE ARGS. 
#TODO: All arguments are manadatory for now, will add defaults later

while len(sys.argv) > 1 and sys.argv[1].find("=") >= 0:
        # Command line arguments containing '=' are implicitly options.
        name,value = sys.argv.pop(1).split('=')
        if __debug__: # ...just to see what's going on.
                print( "{},{}".format( name, value ) )
        name = name.lower().lstrip("-")
        # The following just allows for shorter command lines since
        # my arguments are matrix, rows, cols, and levels.
        if name.startswith("model"):  #model = a(dditive)/d(ominant)/r(ecessive)
                MODEL = value
        elif name.startswith("mode"):    #mode = allelic/genotypic
                MODE = value
        elif name.startswith("fm"):     #feature matrix path, assumes file is gzipped
                FM_FILENAME = value
        elif name.startswith("ids"):
                PEDIDS_FILENAME = value    #pedigree ids file path
	elif name.startswith("output"):       # output file path
		OUTPUT_FILENAME = value	
	elif name.startswith("log"):     #log file path
		LOG_FILENAME = value
        else:
                print "unrecognized option:", name 
                sys.exit(1)

#INPUT: vcf - only full and missing calls, no half calls;
if FM_FILENAME.endswith("gz"):
	fmFile = gzip.open(FM_FILENAME,"r")
else:
	fmFile = open(FM_FILENAME,"r")		

#INPUT: pedigree ids : only case pedigrees, TDT does not use controls
pedIDFile = open(PEDIDS_FILENAME,"r")

#OUTPUT: log file, results file 
logFile = gzip.open(LOG_FILENAME+"_Model"+MODEL.upper()+".log.gz","w")
outputFile = gzip.open(OUTPUT_FILENAME+"_Model"+MODEL.upper()+".out.gz","w")


logFile.write("Running with options:\n")
logFile.write("MODE = "+MODE+"\n")
logFile.write("MODEL = "+MODEL+"\n")
logFile.write("FEATURE MATRIX = "+os.path.abspath(FM_FILENAME)+"\n")
logFile.write("PEDIDS FILE = "+os.path.abspath(PEDIDS_FILENAME)+"\n")
logFile.write("OUTPUT FILE = "+os.path.abspath(OUTPUT_FILENAME)+"\n")
logFile.write("LOG FILE = "+os.path.abspath(LOG_FILENAME)+"\n")
logFile.write("Trio types : \n")
logFile.write("X50 = set([0,1]) -> 0\n")
logFile.write("X51 = set([0,1]) -> 1\n")
logFile.write("X40 = set([1,1]) -> 0\n")
logFile.write("X41 = set([1,1]) -> 1\n")
logFile.write("X42 = set([1,1]) -> 2\n")
logFile.write("X21 = set([1,2]) -> 1\n")
logFile.write("X22 = set([1,2]) -> 2\n")
logFile.write("Start time : "+ time.strftime("%Y-%m-%d %H:%M:%S",time.localtime())+"\n")

pedIDs = np.array(pedIDFile.read().split())
headerColumns = []
pedIndices = []

#look-up table for genotype to number encoding
genotypeCode = {'0/0':0,'0':0,'0/1':1,'1/0':1,'1':1,'1/1':2,'./.':-1,'.':-1}

#look-up table for complete trio genotypes [set('F','M'),'NB']
#Informative mating types
#Mating type 2: set(1,2)
#Mating type 4: set(1,1)
#Mating type 5: set(0,1)
#trio type notation: [X_50, X_51, X_40, X_41, X_42, X_21, X_22]
informativeTrioType = [[set([0,1]),0],
			[set([0,1]),1],
			[set([1,1]),0],
			[set([1,1]),1],
			[set([1,1]),2],
			[set([1,2]),1],
			[set([1,2]),2]]


### DOM and REC REFERENCE: Schaid, Sommer (1994) Comparison of Statistics for Candidate-Gene Association Studies Using Cases and Parents
#Note: parent set (1,2) is not counted for dominant model
#Note: parent set (0,1) is not counted for recessive model

#keep track of no. of informative trios, non informative trios, and MIEs
nInformative = 0
nNonInformative = 0
nMissing = 0

#counts of trio types
nInformativeTrio = [0]*7

#print header line 
if MODE.lower() == "allelic":
	outputColumns = ['CHR','POS','GENE','REF','ALT','TotalAlleleCount','ALTAlleleCount','nInformativeFam','nNonInformativeFam','nMissingFam','nInformativeType']
	if MODEL.lower() == "all":
		outputColumns.extend(['ChiSq_TDT_Additive','P-value_TDT_Additive','ChiSq_TDT_Dominant','P-value_TDT_Dominant','ChiSq_TDT_Recessive','P-value_TDT_Recessive'])
	elif MODEL.lower() == "a":    
		outputColumns.extend(['ChiSq_TDT_Additive','P-value_TDT_Additive'])
	elif MODEL.lower() == "d":
		outputColumns.extend(['ChiSq_TDT_Dominant','P-value_TDT_Dominant'])
	elif MODEL.lower() == "r":
		outputColumns.extend(['ChiSq_TDT_Recessive','P-value_TDT_Recessive'])
	else:
		print "Unrecognised model: "+MODEL
        	sys.exit(1)
	
outputFile.write('\t'.join(outputColumns)+'\n')

#read vcf one line at a time
for line in fmFile:
	vcfValues = line.strip().split('\t')
	if vcfValues[0].find("chr") <> -1:	## not a header row, valid marker row
		
		 
		#reset for each new marker
		nInformativeTrio = [0]*7
		nInformative = 0
		nNonInformative = 0
		nMissing = 0
		
		#loop through each pedigree from the list pedIDs
		for ped in pedIndices:
			# get numeric genotype code
			genoF = genotypeCode[vcfValues[ped[0]]]
			genoM = genotypeCode[vcfValues[ped[1]]]
			genoNB = genotypeCode[vcfValues[ped[2]]]

			#if complete trio info available
			if genoF <> -1 and genoM <> -1 and genoNB <> -1:
				try:
					nInformativeTrio[informativeTrioType.index([set([genoF,genoM]),genoNB])] += 1
					nInformative += 1
					logText = '\t'.join([vcfValues[0],vcfValues[1],'Complete and Informative',','.join([headerColumns[ped[0]],headerColumns[ped[1]],headerColumns[ped[2]]]),','.join([vcfValues[ped[0]],vcfValues[ped[1]],vcfValues[ped[2]]])])+'\n'
					logFile.write(logText)
				except (KeyError,ValueError):    # trio not informative for this model
					nNonInformative += 1
	                                logText = '\t'.join([vcfValues[0],vcfValues[1],'Complete but Non-Informative',','.join([headerColumns[ped[0]],headerColumns[ped[1]],headerColumns[ped[2]]]),','.join([vcfValues[ped[0]],vcfValues[ped[1]],vcfValues[ped[2]]])])+'\n'
					logFile.write(logText)
					pass
		
			else:     #atleast one member genotype missing; incomplete trio
				nMissing += 1
				logText = '\t'.join([vcfValues[0],vcfValues[1],'Incomplete',','.join([headerColumns[ped[0]],headerColumns[ped[1]],headerColumns[ped[2]]]),','.join([vcfValues[ped[0]],vcfValues[ped[1]],vcfValues[ped[2]]])])+'\n'
                                logFile.write(logText)

		#standard TDT
		
		if MODEL.lower() == "all" or MODEL.lower() == "a":
			try:
				#Note: nInformativeTrio[3] counts in both b and c. adding or removing it does not affect the numerator b-c, but does affect the denominatore b+c
				# so we keep it, beacuse this trio transmits 1 REF allele and 1 ALT allele
	
				b = nInformativeTrio[6] + 2*nInformativeTrio[4] + nInformativeTrio[1] + nInformativeTrio[3]
				c = nInformativeTrio[5] + 2*nInformativeTrio[2] + nInformativeTrio[0] + nInformativeTrio[3]
				chiSqTDTAdditive = (b-c)**2/float(b+c)
				pValueAdditive = 1 - stats.chi2.cdf(chiSqTDTAdditive, 1)
			except ZeroDivisionError:
	                        chiSqTDTAdditive = 'NA'
        	                pValueAdditive = 'NA'

		if MODEL.lower() == "all" or MODEL.lower() == "d":   ##equation 11 of the referenced paper
			try:
				totalN4 = nInformativeTrio[2] +  nInformativeTrio[3] +  nInformativeTrio[4]
				totalN5 = nInformativeTrio[0] +  nInformativeTrio[1]
				chiSqTDTDominant = (nInformativeTrio[1] + nInformativeTrio[3] + nInformativeTrio[4] - totalN5/float(2) - 3*totalN4/float(4))**2/(totalN5/float(4) + 3*totalN4/float(16))
				pValueDominant = 1 - stats.chi2.cdf(chiSqTDTDominant, 1)
			except ZeroDivisionError:
                        	chiSqTDTDominant = 'NA'
       	                        pValueDominant = 'NA'

			## Convolution of binomials

#			binomType4 = stats.binom.pmf(range(totalN4+1),totalN4,0.75)
#                       binomType5 = stats.binom.pmf(range(totalN5+1),totalN5,0.50)

#			if totalN4 != 0 and totalN5 == 0:
#				binomDOMPValue = sum(binomType4[nInformativeTrio[3]+nInformativeTrio[4]:len(binomType4)])
#			elif totalN5 != 0 and totalN4 == 0:
#				binomDOMPValue = sum(binomType5[nInformativeTrio[1]:len(binomType5)])
#			elif totalN4 != 0 and totalN5 != 0:
#				binomDOM = np.convolve(binomType4,binomType5)
#				print binomType5,binomType4,binomDOM
#				binomDOMPValue = sum(binomDOM[nInformativeTrio[1] + nInformativeTrio[3] + nInformativeTrio[4]:len(binomDOM)])
#			else:
#				binomDOMPValue = 'NA'
#	
#			print vcfValues[0],vcfValues[1],nInformativeTrio,pValueDominant, binomDOMPValue
		if MODEL.lower() == "all" or MODEL.lower() == "r":   ##equation 12 of the referenced paper
			try:
				totalN4 = nInformativeTrio[2] +  nInformativeTrio[3] +  nInformativeTrio[4]
                                totalN2 = nInformativeTrio[5] +  nInformativeTrio[6]
                                chiSqTDTRecessive = (nInformativeTrio[4] + nInformativeTrio[6] - totalN4/float(4) - totalN2/float(2))**2/(totalN2/float(4) + 3*totalN4/float(16))
                                pValueRecessive = 1 - stats.chi2.cdf(chiSqTDTRecessive, 1)
			except ZeroDivisionError:
				chiSqTDTRecessive = 'NA'
				pValueRecessive = 'NA'
			
			## Convolution of binomials

#                       binomType4 = stats.binom.pmf(range(totalN4+1),totalN4,0.25)
#                       binomType2 = stats.binom.pmf(range(totalN2+1),totalN2,0.50)

#                       if totalN4 != 0 and totalN2 == 0:
#                              binomRECPValue = sum(binomType4[nInformativeTrio[4]:len(binomType4)])
#                       elif totalN2 != 0 and totalN4 == 0:
#                              binomRECPValue = sum(binomType2[nInformativeTrio[6]:len(binomType2)])
#                       elif totalN4 != 0 and totalN2 != 0:
#                               binomREC = np.convolve(binomType4,binomType2)
#                               binomRECPValue = sum(binomREC[nInformativeTrio[4]+nInformativeTrio[6]:len(binomREC)])
#                       else:
#                               binomRECPValue = 'NA'

#                       print vcfValues[0],vcfValues[1],nInformativeTrio,pValueRecessive, binomRECPValue

		
		#print CHR, POS, GENE_Name,REF, ALT, AN, AC, nInformativeFam, nNonInformative, nMissing, nInformativeTypes, chiSqTDT, p-valueTDT
		#get GENE_NAMES
	        geneNameStartIndex = vcfValues[7].find('GENE_NAMES')
		if geneNameStartIndex <> -1:
			geneNameEndIndex = vcfValues[7].find(';',geneNameStartIndex)
			geneNames = vcfValues[7][geneNameStartIndex+11:geneNameEndIndex]	
		else:
			geneNames = "NA"
			
		# calculate ALT allele frequency as AC/AN; works fine IF the input data has only biallelic markers
                # INFO field index in the vcf = 8

		indexAN = vcfValues[7].find(";AN=")
                totalAlleles = vcfValues[7][indexAN+4:vcfValues[7].find(";",indexAN+1)]
                indexAC = vcfValues[7].find(";AC=")
                # THIS LOGIC NEEDS TO BE REVISITED.
                # As of 08/02/2013, if ALT allele count = 0, AC is not reported, and ALT is reported as "."
                if indexAC <> -1:
                	m = re.search(r"(\d+)",vcfValues[7][indexAC+4:indexAC+10])
			altAlleleCount = int(m.group(0))
                else:
                	altAlleleCount = 0

                alleleFreq = int(altAlleleCount)/float(totalAlleles)

		
		#concatenate output string and print to output file
		outputColumns = [vcfValues[0],vcfValues[1],geneNames,vcfValues[3],vcfValues[4],totalAlleles, str(altAlleleCount), str(nInformative),str(nNonInformative),str(nMissing),str(nInformativeTrio)]
		if MODEL.lower() == "all":
			outputColumns.extend([str(chiSqTDTAdditive),str(pValueAdditive),str(chiSqTDTDominant),str(pValueDominant),str(chiSqTDTRecessive),str(pValueRecessive)])	
		elif MODEL.lower() == "a":
			outputColumns.extend([str(chiSqTDTAdditive),str(pValueAdditive)])
		elif MODEL.lower() == "d":
			outputColumns.extend([str(chiSqTDTDominant),str(pValueDominant)])	
		elif MODEL.lower() == "r":
			outputColumns.extend([str(chiSqTDTRecessive),str(pValueRecessive)])

		outputFile.write('\t'.join(outputColumns)+'\n')

	elif vcfValues[0].find("CHROM") <> -1: #this line specifies sample IDs
		headerColumns = vcfValues 
		pedIndices = []
		for id in pedIDs:
                        memberIndices = [headerColumns.index(x) for x in headerColumns if x.find(id) <> -1]   #in order F, M, NB
			pedIndices.append(memberIndices)



logFile.write("Completion time : "+ time.strftime("%Y-%m-%d %H:%M:%S",time.localtime())+"\n")

fmFile.close()
pedIDFile.close()
logFile.close()				
outputFile.close()

