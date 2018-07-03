import numpy as np
for i in range(0,100):
	p = np.random.randint(0,50)
	n = np.random.randint(4,11)
	print './flatfold_gen '+str(n)+' '+str(i)+' '+str(p)+' 1'
	print 'gnuplot fold.gnuplot'
	print 'mv fold_fig.png ./images/RadFlat_'+str(i)+'.png'
