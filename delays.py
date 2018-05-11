import numpy as np, numpy.linalg as linalg
import matplotlib.pyplot as plt
import getopt, sys

def main():
    global a, b

    try:
	opts, args = getopt.getopt(sys.argv[1:], 'ha:b:', ['aone=', 'btwo='])
    except getopt.GetoptError:
	sys.exit(2)
    for opt, arg in opts:
	if opt in ('-h', '--help'):
	    sys.exit()
	elif opt in ('-a', '--aone'):
	    a = arg
	elif opt in ('-b', '--btwo'):
	    b = arg

def reading(a, b):
    len1 = sum(1 for line in open(a))
    len2 = sum(1 for line in open(b))
    
    a = open(a)
    b = open(b)
    
    li = []
    
    lenn = min(len1, len2)
    
    i = 0
    while i < lenn:
	li.append(float(a.readline()) - float(b.readline()) * 1e6)
	i = i + 1
    
    xlist = np.linspace(1, lenn, lenn)

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
    reading(a, b)