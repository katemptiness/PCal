#!/usr/bin/env python

import numpy as np, numpy.linalg as linalg
import matplotlib.pyplot as plt
import getopt, sys
import os

def main():
    global files, exception

    exception = -1

    try:
	opts, args = getopt.getopt(sys.argv[1:], 'hf:e:', ['files=', 'exception='])
    except getopt.GetoptError:
	print 'Something is wrong...\n'
	usage()
	sys.exit(2)
    for opt, arg in opts:
	if opt in ('-h', '--help'):
	    usage()
	    sys.exit()
	elif opt in ('-f', '--files'):
	    files = arg
	elif opt in ('-e', '--exception'):
	    exception = int(arg) * 2


def usage():
    print 'Use -f for path to your files.'
    print 'You can also use -e for number of pair of files you want to exclude from sample.'
    print
    print 'For example: ./delays.py -f ~/files -e 12'
    print 'or just : ./delays.py -f ~/files'


def reading(files):
    all_files = os.listdir(files)
    all_files = sorted(all_files)
    
    len_files = len(all_files)
    
    k = 0
    while k < len_files:
	all_files[k] = files + all_files[k]
	k = k + 1
    
    li = []
    sp = []
    x1 = 0
    ms = []
    cs = []
    
    i = 0
    k = 0
    while i < len_files:
    	if i != exception:
    	    start = ((all_files[i]).find('PCAL_') + 5) + 6
	    stop = start + 6

	    len1 = sum(1 for line in open(all_files[i]))
	    len2 = sum(1 for line in open(all_files[i + 1]))
	    lenn = min(len1, len2)

	    sp.append(float((all_files[i])[start : stop]))
	    
	    if k > 0:
	        space = (sp[k] - sp[k - 1]) * 2

	    ifile1 = open(all_files[i])
	    ifile2 = open(all_files[i + 1])
	
	    li = []
	    j = 0
	    while j < lenn:
	        ch1 = float(ifile1.readline())
	        ch2 = float(ifile2.readline())
	        ch = abs(ch1 - ch2) * 1e6

	        li.append(ch)

	        j = j + 1
	
	    x2 = x1 + lenn

	    if k > 0:
	
	        xlist = np.linspace((x1 * 0.5), ((x2 - 1) * 0.5), lenn)
	
	        A = (np.vstack([xlist, np.ones(len(xlist))])).transpose()
	        m, c = linalg.lstsq(A, li, rcond = -1)[0]
	        trend = m * xlist + c
	    
	        h = 0
	        while h < lenn:
		    ms.append(m)
		    cs.append(c)
		    h = h + 1

	        plt.plot(xlist, li, 'o')
	        plt.plot(xlist, trend)
	
	        x1 = x2 + space

	        h = 0
	        while h < space:
		    ms.append(0)
		    h = h + 1

	else:
	    start = ((all_files[i]).find('PCAL_') + 5) + 6
	    stop = start + 6

	    len1 = sum(1 for line in open(all_files[i]))
	    len2 = sum(1 for line in open(all_files[i + 1]))
	    lenn = min(len1, len2)

	    sp.append(float((all_files[i])[start : stop]))
	    
	    if k > 0:
	        space = (sp[k] - sp[k - 1]) * 2

	    x2 = x1 + lenn

	    if k > 0:
	
	        xlist = np.linspace((x1 * 0.5), ((x2 - 1) * 0.5), lenn)
	
	        x1 = x2 + space

	        h = 0
	        while h < space:
		    ms.append(0)
		    h = h + 1

	k = k + 1
	i = i + 2
    
    c = np.mean(cs)
    m = np.mean(ms)

    xlist = np.linspace(0, xlist[-1], 10)

    trend = m * xlist + c

    plt.plot(xlist, trend)

    print 'The slope is', m
    
    plt.grid()
    plt.xlabel('time, s')
    plt.ylabel('difference between time delays, ps')
    plt.gcf().canvas.set_window_title('Difference between time delays')
    plt.show()


if __name__ == '__main__':
    main()
    reading(files)