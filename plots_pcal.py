#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np, numpy.linalg as linalg
import getopt, sys
import matplotlib.pyplot as plt


def main():
    global ifile, unw

    unw = 'false'

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hf:u:', ['ifile=', 'unw='])
    except getopt.GetoptError:
        usage()
        sys.exit()
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-f', '--ifile'): ifile = arg
        elif opt in ('-u', '--unw'): unw = org

    
def usage():
    print 'You should use this form for work:'
    print '-f is for path to the file;'
    print '-u is for unwraping mode (true or false);'


def unwraping(lista):
   	for k in range(len(lista)): lista[k] = ((lista[k] / 360) - ((lista[k] / 360) // 1)) * 360
   	return lista


def pcal_plot(ifile):
	ifile = open(ifile)

	q = len(ifile.readlines())
	ifile.seek(0)

	ifile_massive = []
	for k in range(q): ifile_massive.append(float((ifile.readline())[0 : -1]))
	ifile.close()

	plt.figure(1)
	plt.gcf().canvas.set_window_title(u'Зависимость разности фаз от частоты')
	plt.plot(ifile_massive, 'o')
	plt.plot(ifile_massive)
	plt.grid()
	plt.xticks(np.arange(0, (len(ifile_massive) + 32), step = 32))
	plt.yticks(np.arange(-420, 420, step = 60))
	plt.xlabel('frequency, MHz')
	plt.ylabel('phase, grad')
	plt.show(block = False)
	
	massive = []
	for k in range((len(ifile_massive) - 1)): massive.append(ifile_massive[k + 1] - ifile_massive[k])
	massive = unwraping(unwraping(massive))

	m = np.mean(massive)
	s = np.std(massive)
	print 'Mean is', m
	print 'STD is', s

	plt.figure(2)
	plt.gcf().canvas.set_window_title(u'Зависимость разности разностных фаз от частоты')
	plt.plot(massive, 'o')
	plt.plot(massive)
	plt.grid()
	plt.xticks(np.arange(0, (len(ifile_massive) + 32), step = 32))
	plt.xlabel('frequency, MHz')
	plt.ylabel('phase, grad')
	plt.show()


if __name__ == '__main__':
	main()

	pcal_plot(ifile)