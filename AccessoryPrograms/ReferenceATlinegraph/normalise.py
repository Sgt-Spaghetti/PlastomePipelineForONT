import pandas as pd
import sys

Input_Data = pd.read_csv(sys.argv[1])

print(Input_Data)

loci = Input_Data['Position'].to_list()
percentages = Input_Data['AT Percentage'].to_list()

normalised_loci = []
max_length = max(loci)
for i in loci:
	normalised_loci.append(i/max_length)

new_dataframe = pd.DataFrame({'Position':normalised_loci, 'AT Percentage':percentages})

new_dataframe.to_csv('normalised_'+sys.argv[1], index=False)
