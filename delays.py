#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np, numpy.linalg as linalg
import matplotlib.pyplot as plt
import os
import getopt, sys
import scipy.linalg

def main():
	global files, exception, trending, mode

	exception = -1
	mode = 'diff'
	trending = 'linear'

	try:
		opts, args = getopt.getopt(sys.argv[1:], 'hf:e:t:m:', ['files=', 'exception=', 'trending=', 'mode='])
	except getopt.GetoptError:
		print 'Looks like something is wrong...'
		usage()
		sys.exit()
	for opt, arg in opts:
		if opt in ('-h', '--help'):
			usage()
			sys.exit()
		elif opt in ('-e', '--exception'):
			try:
				exception = int(arg)
				if exception <= 0:
					print 'Error: -e parameter should be more than 0.'
					sys.exit()
			except:
				print 'Error: -e parameter should be one integer value.'
				sys.exit()
		elif opt in ('-f', '--files'): files = arg
		elif opt in ('-t', '--trending'): trending = arg
		elif opt in ('-m', '--mode'): mode = arg


def usage():
	print 'Use -f for path to your files. For example: ./delays.py -f ~/files'
	print 'Use -m to select display mode (diff or only 1 station). For example: ./delays.py -f ~/files -m diff'
	print 'or: ./delays.py -f ~/files -m 1'
	print 'or: ./delays.py -f ~/files -m 2'
	print 'Use -t for the trending method (linear or quadratic). For example: ./delays -f ~/files -t quadratic'
	print 'Use -e to exclude the value. For example: ./delays.py -f ~/files -e 14'


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
		if mode == 'diff':
			for j in range(lens[i]): delay_difs.append(abs(float(file1.readline()) - float(file2.readline())) * 1e6)
		elif mode == '1':
			for j in range(lens[i]): delay_difs.append(float(file1.readline()) * 1e6)
		elif mode == '2':
			for j in range(lens[i]): delay_difs.append(float(file2.readline()) * 1e6)
		i = i + 1
		k = k + 2

	if trending == 'quadratic':
		trendy = []

	x1 = 1
	y1 = 0
	ms = []
	cs = []
	for i in range(len(lens)):
		x2 = x1 + lens[i]
		y2 = y1 + lens[i]
		xlist = np.linspace((x1 / 2), (x2 / 2), (x2 - x1))
		ylist = delay_difs[y1 : y2]
		if trending == 'linear':
			A = (np.vstack([xlist, np.ones(len(xlist))])).transpose()
			m, c = linalg.lstsq(A, ylist, rcond = -1)[0]
		if i != exception:
			if trending == 'linear':
				for k in range(len(delay_difs[y1 : y2])):
					ms.append(m)
					cs.append(c)
				trend = m * xlist + c
			elif trending == 'quadratic':
				trendy.append(ylist)
			plt.plot(xlist, ylist, 'o')
			if trending == 'linear':
				plt.plot(xlist, trend)
		else:
			if trending == 'linear':
				for k in range(len(delay_difs[y1 : y2])): ms.append(0)
		try:
			x1 = x1 + spaces[i] * 2
		except:
			pass
		y1 = y2

	if trending == 'quadratic':
		new_trendy = []
		for i in range(len(trendy)): new_trendy.append(np.mean(trendy[i]))

		i_s = []

		for i in range(len(all_files) - 1):
			if (float((all_files[i + 1])[11 : 17]) - float((all_files[i])[11 : 17])) > 1e3: i_s.append((i + 1) / 2)
	
		i_s.append(int(len(all_files) / 2))

		j = 0
		mins = []

		for i in range(len(i_s)):
			k = []
			while j < i_s[i]:
				k.append(new_trendy[j])
				j = j + 1
			mins.append(np.mean(k))

		x = np.linspace(1, xlist[-1], len(mins))

		m = scipy.vstack((x ** 2, x, np.ones(len(mins)))).transpose()
		s = scipy.linalg.lstsq(m, mins)[0]
		super_xlist = np.linspace(0, xlist[-1], 100)
		
		plt.plot(super_xlist, s[0] * super_xlist ** 2 + s[1] * super_xlist + s[2])

		print 'The coefficients are:', s[0], 'and', s[1]
			
	elif trending == 'linear':
		for j in range(len(spaces)):
			for i in range(spaces[j]): ms.append(0)

		m = np.mean(ms)
		c = np.mean(cs)
		super_xlist = np.linspace(0, xlist[-1], 100)
		super_trend = m * super_xlist + c
		plt.plot(super_xlist, super_trend)

		print 'The slope is', m

	plt.grid()
	plt.xlabel(u'время, с')
	if mode == 'dif':
		plt.ylabel(u'разностная задержка, пс')
		plt.gcf().canvas.set_window_title(u'Разностная задержка')
	else:
		plt.ylabel(u'задержка, пс')
		plt.gcf().canvas.set_window_title(u'Временная задержка')
	plt.show()


if __name__ == '__main__':
	main()
	reading(files)