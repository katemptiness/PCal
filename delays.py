#!/usr/bin/env python

import numpy as np, numpy.linalg as linalg
import matplotlib.pyplot as plt
import os
import getopt, sys

def main():
	global files, exception

	exception = -1

	try:
		opts, args = getopt.getopt(sys.argv[1:], 'hf:e:', ['files=', 'exception='])
	except getopt.GetoptError:
		print 'Looks like something is wrong...'
		usage()
		sys.exit()
	for opt, arg in opts:
		if opt in ('-h', '--help'):
			usage()
			sys.exit()
		elif opt in ('-f', '--files'):
			files = arg
		elif opt in ('-e', '--exception'):
			try:
				exception = int(arg)
				if exception <= 0:
					print 'Error: -e parameter should be more than 0.'
					sys.exit()
			except:
				print 'Error: -e parameter should be one integer value.'
				sys.exit()


def usage():
	print 'Use -f for path to your files. For example: ./delays.py -f ~/files'
	print 'You can also use -e to exclude the value. For example: ./delays.py -f ~/files -e 14'


def reading(files):
	all_files = sorted(os.listdir(files))

	spaces = []
	i = 0
	while i < (len(all_files) - 2):
		spaces.append(int((all_files[i + 2])[11 : 17]) - int((all_files[i])[11 : 17]))
		i = i + 2

	lens = []
	j = 0
	while j < len(all_files):
		lens.append(min(sum(1 for line in open(files + all_files[j])), sum(1 for line in open(files + all_files[j + 1]))))
		j = j + 2

	delay_difs = []
	k = 0
	i = 0
	while k < (len(all_files) - 1):
		file1 = open(files + all_files[k])
		file2 = open(files + all_files[k + 1])
		j = 0
		while j < lens[i]:
			delay_difs.append(abs(float(file1.readline()) - float(file2.readline())) * 1e6)
			j = j + 1
		i = i + 1
		k = k + 2

	x1 = 1
	y1 = 0
	ms = []
	cs = []
	i = 0
	while i < len(lens):
		x2 = x1 + lens[i]
		y2 = y1 + lens[i]
		xlist = np.linspace((x1 / 2), (x2 / 2), (x2 - x1))
		ylist = delay_difs[y1 : y2]
		A = (np.vstack([xlist, np.ones(len(xlist))])).transpose()
		m, c = linalg.lstsq(A, ylist, rcond = -1)[0]
		if i != exception:
			k = 0
			while k < len(delay_difs[y1 : y2]):
				ms.append(m)
				cs.append(c)
				k = k + 1
			trend = m * xlist + c
			plt.plot(xlist, ylist, 'o')
			plt.plot(xlist, trend)
		else:
			k = 0
			while k < len(delay_difs[y1 : y2]):
				ms.append(0)
				k = k + 1
		try:
			x1 = x1 + spaces[i] * 2
		except:
			x1 = x1
		y1 = y2
		i = i + 1

	j = 0
	while j < len(spaces):
		i = 0
		while i < spaces[j]:
			ms.append(0)
			i = i + 1
		j = j + 1

	m = np.mean(ms)
	c = np.mean(cs)
	super_xlist = np.linspace(0, xlist[-1], 100)
	super_trend = m * super_xlist + c
	plt.plot(super_xlist, super_trend)

	print 'The slope is', m

	plt.grid()
	plt.xlabel('time, s')
	plt.ylabel('differences between time delays, ps')
	plt.gcf().canvas.set_window_title('Difference between time delays')
	plt.show()


if __name__ == '__main__':
	main()
	reading(files)