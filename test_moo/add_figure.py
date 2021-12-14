import matplotlib.pyplot as plt
def add_figure(title,x_label,y_label,show=False):
	plt.close()
	plt.figure()
	plt.title(title)
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	if show:
		plt.show()