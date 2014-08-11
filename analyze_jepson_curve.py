import numpy as np

prop = []
for i in bigv:
	cnt = []
	for j in i:
		if max(j) > 30: cnt.append(1)
		else:
			cnt.append(0)
	
	cnt = np.array(cnt)
	prop.append(sum(cnt)/len(cnt))

