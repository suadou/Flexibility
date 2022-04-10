import matplotlib.pyplot as plt
from pandas import *
def represent_data(csv_file, prefix):
	
	data = read_csv(csv_file, sep = '\t')
	res = []
	flex_score = data['P_1'].tolist()
	
	i = 1
	for element in flex_score:
		res.append(int(i))
		i = i+1
	
	plt.scatter(res,flex_score, color='red', marker='.')
	plt.title("MEDUSA's Flexibility score vs residue position - "+prefix)
	plt.xlabel('Residue')
	plt.ylabel('Flexibility score')
	plt.savefig(prefix+'.png')
	plt.clf()

if __name__ == "__main__":
	represent_data("CDC24_MEDUSA/prediction/S_prediction.csv", "P11433 (CDC24)")
	represent_data("EXOC1_MEDUSA/prediction/S_prediction.csv", "Q9VVG4 (EXOC1)")
	represent_data("GNAI1_MEDUSA/prediction/S_prediction.csv", "P38401 (GNAI1)")
	represent_data("MCSB_MEDUSA/prediction/S_prediction.csv", "P65206 (MCSB)")
	represent_data("PDE10_MEDUSA/prediction/S_prediction.csv", "Q9Y223 (PDE10)")
	represent_data("PRGR_MEDUSA/prediction/S_prediction.csv", "P06401 (PRGR)")
	represent_data("SEC18_MEDUSA/prediction/S_prediction.csv", "Q9P7Q4 (SEC18)")
