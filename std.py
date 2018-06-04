#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import os
import getopt, sys

def main():
	global files

	try:
		opts, args = getopt.getopt(sys.argv[1:], 'hf:', ['files='])
	except getopt.GetoptError:
		usage()
		sys.exit()
	for opt, arg in opts:
		if opt in ('-h', '--help'):
			usage()
			sys.exit()
		elif opt in ('-f', '--files'): files = arg


def usage():
	print 'Use -f for path to your files. For example: ./std.py -f ~/files'
	

def reading(files):
	all_files = sorted(os.listdir(files))

	for k in range(len(all_files)):
		ifile = open(files + all_files[k])
		std_s = []
		for i in range(512): std_s.append(float(ifile.readline()))
		plt.plot(std_s)

	plt.grid()
	plt.xlabel(u'номера тонов, ед.')
	plt.ylabel(u'СКО, градусы')
	plt.gcf().canvas.set_window_title(files)
	plt.show()


if __name__ == '__main__':
	main()
	reading(files)