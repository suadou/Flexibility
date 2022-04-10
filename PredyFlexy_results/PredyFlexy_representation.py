import matplotlib.pyplot as plt

def represent_data(txt_file, prefix):
	
	res = []
	flex_score = []
	with open(txt_file) as file:
		for line in file:
		    row = line.split()
		    res.append(float(row[0]))
		    flex_score.append(float(row[9]))
			
	i = 0
	for element in flex_score:
		if element == -9.999:
			flex_score[i] = -0.75
		
		i = i+1
	
	plt.scatter(res,flex_score, color='red', marker='.')
	plt.title("PredyFlexy's flexibility score vs residue position - "+prefix)
	plt.xlabel('Residue')
	plt.ylabel('Flexibility score')
	plt.savefig(prefix+'.png')
	plt.clf()

if __name__ == "__main__":
	represent_data("predyflexy_P06401.txt", "P06401 (PRGR)")
	represent_data("predyflexy_P65206.txt", "P65206 (MCSB)")
	represent_data("predyflexy_Q9Y223.txt", "Q9Y223 (PDE10)")
	represent_data("predyflexy_P11433.txt", "P11433 (CDC24)")
	represent_data("predyflexy_Q9P7Q4.txt", "Q9P7Q4 (SEC18)")
	represent_data("predyflexy_P38401.txt", "P38401 (GNAI1)")
	represent_data("predyflexy_Q9VVG4.txt", "Q9VVG4 (EXOC1)")
