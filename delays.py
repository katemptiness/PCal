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
    all_files = sorted(all_files)
    
    len_files = len(all_files)
    
    k = 0
    while k < len_files:
	all_files[k] = files + all_files[k]
	k = k + 1
    
    li = []
    sp = []
    x1 = 0
    
    i = 0
    k = 0
    while i < len_files:
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
	
	    xlist = np.linspace(x1, (x2 - 1), lenn)
	
	    A = (np.vstack([xlist, np.ones(len(xlist))])).transpose()
	    m, c = linalg.lstsq(A, li, rcond = -1)[0]
	    trend = m * xlist + c
	
	    plt.plot(xlist, li, 'o')
	    plt.plot(xlist, trend)
	
	    x1 = x2 + space
	
	k = k + 1
	i = i + 2
    
    plt.grid()
    plt.xlabel('accumulation periods')
    plt.ylabel('difference between time delays, ps')
    plt.gcf().canvas.set_window_title('Difference between time delays')
    plt.show()


if __name__ == '__main__':
    main()
    reading(files)