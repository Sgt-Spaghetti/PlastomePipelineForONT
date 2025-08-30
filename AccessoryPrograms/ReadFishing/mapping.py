import pandas as pd
import sys

input_file_location = sys.argv[1]

chromosomes = {
"chromosome_01" : (8225636,1),
"chromosome_02" : (8655884,2),
"chromosome_03" : (9286894,3),
"chromosome_04" : (4130073,4),
"chromosome_05" : (3682160,5),
"chromosome_06" : (8913359,6),
"chromosome_07" : (6492107,7),
"chromosome_08" : (4526983,8),
"chromosome_09" : (6807148,9),
"chromosome_10" : (6800247,10),
"chromosome_11" : (4479522,11),
"chromosome_12" : (9952739,12),
"chromosome_13" : (5281438,13),
"chromosome_14" : (4217303,14),
"chromosome_15" : (5870643,15),
"chromosome_16" : (8042475,16),
"chromosome_17" : (6954842,17)
}

step_size = int(sys.argv[2])

bins = []
FD = []

chromosome_number = 0
for c in chromosomes:
	bins.append([])
	FD.append([])
	i = 0
	while i < chromosomes[c][0]: 
		bins[chromosome_number].append(i)
		FD[chromosome_number].append(0)
		i += step_size
	chromosome_number += 1

data = pd.read_csv(input_file_location, sep=",")

chromosome_hits = data["Chromosome"].to_list()
positions = data["Location"].to_list()

for index, entry in enumerate(chromosome_hits):
	if entry in chromosomes:
		pos = positions[index]
		binned = False
		i = 0
		while binned == False:
			if i < len(bins[chromosomes[entry][1]-1]):
				if pos > bins[chromosomes[entry][1]-1][i]:
					i += 1
				else:
					binned = True
					#print(pos,bins[chromosomes[entry][1]-1][i-1], "binned", chromosomes[entry][1]-1, i-1)
					FD[chromosomes[entry][1]-1][i-1] += 1
			else:
				binned = True
				FD[chromosomes[entry][1]-1][len(bins[chromosomes[entry][1]-1])-1] += 1


max_length = 0
for i in range(len(bins)):
	max_length = max(len(bins[i][:]), max_length)	

for i in range(len(bins)):
	if len(bins[i]) < max_length:
		bins[i] += [None]*(max_length-len(bins[i]))
for j in range(len(FD)):
	if len(FD[j]) < max_length:
		FD[j] += [None]*(max_length-len(FD[j]))

out = pd.DataFrame({"chromosome_1 bins": bins[0], "chromosome_1 FD" : FD[0], "chromosome_2 bins": bins[1], "chromosome_2 FD" : FD[1], "chromosome_3 bins": bins[2], "chromosome_3 FD" : FD[2], "chromosome_4 bins": bins[3], "chromosome_4 FD" : FD[3], "chromosome_5 bins": bins[4], "chromosome_5 FD" : FD[4], "chromosome_6 bins": bins[5], "chromosome_6 FD" : FD[5], "chromosome_7 bins": bins[6], "chromosome_7 FD" : FD[6], "chromosome_8 bins": bins[7], "chromosome_8 FD" : FD[7], "chromosome_9 bins": bins[8], "chromosome_9 FD" : FD[8], "chromosome_10 bins": bins[9], "chromosome_10 FD" : FD[9], "chromosome_11 bins": bins[10], "chromosome_11 FD" : FD[10], "chromosome_12 bins": bins[11],"chromosome_12 FD": FD[11],"chromosome_13 bins": bins[12],  "chromosome_13 FD" : FD[12], "chromosome_14 bins": bins[13], "chromosome_14 FD" : FD[13], "chromosome_15 bins": bins[14], "chromosome_15 FD" : FD[14], "chromosome_16 bins": bins[15],  "chromosome_16 FD" : FD[15], "chromosome_17 bins": bins[16],  "chromosome_17 FD" : FD[16]})
print(out.to_string())
out.to_csv("chromosome_mapping_output.csv")
