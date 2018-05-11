import numpy as np, numpy.linalg as linalg
import matplotlib.pyplot as plt
import getopt, sys
import os

def main():
    global files

    try:
	opts, args = getopt.getopt(sys.argv[1:], 'hf:', ['files='])
    except getopt.GetoptError:
	sys.exit(2)
    for opt, arg in opts:
	if opt in ('-h', '--help'):
	    sys.exit()
	elif opt in ('-f', '--files'):
	    files = arg

def reading(files):
    all_files = os.listdir(files)
    
    len_files = len(all_files)
    
    k = 0
    while k < len_files:
	all_files[k] = files + all_files[k]
	k = k + 1
    
    li = []
    
    i = 0
    while i < len_files:
	len1 = sum(1 for line in open(all_files[i]))
	len2 = sum(1 for line in open(all_files[i + 1]))
	
	lenn = min(len1, len2)
	
	ifile1 = open(all_files[i])
	ifile2 = open(all_files[i + 1])
	
	j = 0
	while j < lenn:
	    ch1 = float(ifile1.readline())
	    ch2 = float(ifile2.readline())
	    ch = abs(ch1 - ch2) * 1e6
	    li.append(ch)
	    j = j + 1
	
	i = i + 2
    
    xlist = np.linspace(1, lenn * (i / 2), lenn * (i / 2))
    A = (np.vstack([xlist, np.ones(len(xlist))])).transpose()
    m, c = linalg.lstsq(A, li, rcond = -1)[0]
    trend = m * xlist + c
    
    print '\nThe slope of trend line is', "%.3f" % (m)

    plt.plot(xlist, li, 'o')
    plt.plot(xlist, trend)
    plt.grid()
    plt.xlabel('accumulation periods')
    plt.ylabel('difference between time delays, ps')
    plt.gcf().canvas.set_window_title('Difference between time delays')
    plt.show()


if __name__ == '__main__':
    main()
    reading(files)