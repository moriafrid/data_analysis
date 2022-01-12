import matplotlib.pyplot as plt
def add_figure(title,x_label,y_label):
	plt.close()
	plt.figure()
	plt.title(title,fontsize=15)
	plt.xlabel(x_label)
	plt.ylabel(y_label)
