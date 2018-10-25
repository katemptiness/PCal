#!/usr/bin/env python

import numpy as np
import getopt, sys
import matplotlib.pyplot as plt

def main():
    global ifile

    ifile = None

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hf:', ['ifile='])
    except getopt.GetoptError:
        usage()
        sys.exit()
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-f', '--ifile'): ifile = arg

    
def usage():
    print 'You should use this form for work:'
    print '-f is for path to the file.'


def pcal_plot(ifile):
	ifile = open(ifile)

	q = len(ifile.readlines())
	ifile.seek(0)

	ifile_massive = []
	for k in range(q): ifile_massive.append(float((ifile.readline())[0 : -1]))
	ifile.close()

	massive = []
	for k in range((len(ifile_massive) - 1)): massive.append(ifile_massive[k + 1] - ifile_massive[k])

	plt.plot(massive, 'o')
	plt.plot(massive)
	plt.grid()
	plt.xticks(np.arange(0, (len(ifile_massive) + 32), step = 32))
	plt.show()


if __name__ == '__main__':
	main()

	pcal_plot(ifile)