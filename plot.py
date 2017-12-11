import matplotlib.pyplot as plt

def plot(hblist,backbone_polars):

	plot_dict = {} #We will do a dict with key=residue, value = type of interaction (backbone = 1, side chain = 2)
	unique_res = [] #And we create a list to not repeat residues

	for hb in sorted (hblist,key=lambda i: i[0].get_serial_number()):

		r1 = hb[0].get_parent()
		r2 = hb[1].get_parent()


#In each iteration we will se the value, 1 for backbone interactions and 2 for side chains interactions
		if hb[1].id in backbone_polars:
			tmp_v = 1               
		else: 
			tmp_v = 2

		if r1.id[1] not in plot_dict.keys():

			plot_dict[r1.id[1]] = [tmp_v]        
		else:
			plot_dict[r1.id[1]].append(tmp_v)     

	plot_table = []


	for key in sorted (plot_dict.keys()):   #Now we look through all the observations we have in the dictionary
		if min(plot_dict[key]) == max(plot_dict[key]):
			plot_table.append([key, min(plot_dict[key])])
		else:
			plot_table.append([key, 3])  


	xlist , ylist = map(list, zip(*plot_table))  #And we split it into two lists with zip function (using map list, because we dont want a tuple)
	plt.scatter(xlist, ylist)
	plt.show()