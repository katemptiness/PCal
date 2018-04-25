#!/usr/bin/env python

import numpy as np, numpy.linalg as linalg
import matplotlib.pyplot as plt

def diff(a, b):
	a = a.replace(',', '')
	a = a.replace('[', '')
	a = a.replace(']', '')
	a = a.split()

	b = b.replace(',', '')
	b = b.replace('[', '')
	b = b.replace(']', '')
	b = b.split()

	i = 0
	while i < len(a):
		a[i] = float(a[i])
		b[i] = float(b[i])
		i = i + 1

	diff = abs(np.asarray(a) - np.asarray(b)) * 1e6

	xlist = np.linspace(1, len(a), len(a))

	A = (np.vstack([xlist, np.ones(len(xlist))])).transpose()
	m, c = linalg.lstsq(A, diff)[0]
	trend = m * xlist + c
	
	print '\nThe slope of trend line is', "%.3f" % (m)

	plt.plot(xlist, trend)
	plt.plot(xlist, diff)
	plt.plot(xlist, diff, 'o')
	plt.grid()
	plt.xlabel('accumulation periods')
	plt.ylabel('difference between time delays, ps')
	plt.gcf().canvas.set_window_title('Difference between time delays')
	plt.show()


if __name__ == '__main__':
	print 'Hello. Please enter the 1st list of time delays:'
	a = raw_input()

	print '\nNow enter the 2nd list of time delays:'
	b = raw_input()

	if len(a.split()) == len(b.split()):
		diff(a, b)
	else:
		print 'Sorry, lists should have same length.'