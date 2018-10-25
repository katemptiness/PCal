#!/usr/bin/env python

import numpy as np, numpy.fft as fft, numpy.linalg as linalg
import cmath
import getopt, sys
import os
import struct

def main():
    global ntone, ifile1, ifile2, average, write

    ifile1 = None
    ifile2 = None
    ntone = '1 : 512'
    average = 'false'
    write = 'false'

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hd:r:n:a:w:', ['ifile1=', 'ifile2=', 'ntone=', 'average=', 'write='])
    except getopt.GetoptError:
        usage()
        sys.exit()
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-r', '--ifile2'): ifile2 = arg                
        elif opt in ('-d', '--ifile1'): ifile1 = arg
        elif opt in ('-n', '--ntone'): ntone = arg
    	elif opt in ('-a', '--average'):
    		if arg == 'true' or arg == 'false': average = arg
    		else:
    			usage()
    			sys.exit()
    	elif opt in ('-w', '--write'):
    		if arg == 'true' or arg == 'false': write = arg
    		else:
    			usage()
    			sys.exit()
        

def usage():
    print 'You can use this forms to work:'
    print '-d is for path to the DiFX file;'
    print '-r is for path to the RASFX file;'
    print '-n is for tone numbers;'
    print '-a is for average mode (true or false);'
    print '-w is for phases recording (true or false).'


def unwraping(lista):
    for k in range(len(lista) - 1):
        if abs(lista[k] - lista[k + 1]) > 345 or abs(lista[k + 1] - lista[k]) > 345:
            list1 = list2 = []
        
            for i in range(len(lista)):
                if lista[i] == abs(lista[i]):
                    if lista[i] > np.mean(lista): list1.append(i)
                    else: list2.append(i)
                elif lista[i] == -abs(lista[i]):
                    if lista[i] < np.mean(lista): list2.append(i)
                    else: list1.append(i)
    
            if len(list1) > len(list2): c = '+'
            elif len(list1) < len(list2): c = '-'
            elif len(list1) == len(list2):
                if lista[0] >= 0: c = '+'
                elif lista[0] < 0: c = '-'
    
            if c == '+':
                for i in range(len(lista) - 1):
                    if (lista[i + 1] - lista[i]) < -180: lista[i + 1] = lista[i + 1] + 360
                    elif (lista[i + 1] - lista[i]) > 180: lista[i] = lista[i] + 360
    
            elif c == '-':
                for i in range(len(lista) - 1):
                    if (lista[i + 1] - lista[i]) < -180: lista[i] = lista[i] - 360
                    elif (lista[i + 1] - lista[i]) > 180: lista[i + 1] = lista[i + 1] - 360

            if np.mean(lista) > 0:
                if lista[0] < 0 and lista[1] < 0:
                    lista[0] = lista[0] + 360
                    lista[1] = lista[1] + 360
            elif np.mean(lista) < 0:
                if lista[0] > 0 and lista[1] > 0:
                    lista[0] = lista[0] - 360
                    lista[1] = lista[1] - 360

    return lista


def ntone_form(ntone):
    global ntones

    ntones = (str(ntone).replace(',', '')).split()
    
    tones = []
    
    i = ntones.count(':')
    while i > 0:
        g = ntones.index(':')
        tones.append(ntones[g - 1])
        tones.append(ntones[g + 1])
        del ntones[g], ntones[g], ntones[g - 1]
        i = i - 1
    
    for i in range(len(ntones)): ntones[i] = int(ntones[i])
    
    tones_1 = []
    
    for i in range(len(tones) / 2):
        g = np.linspace(int(tones[0 + 2 * i]), int(tones[1 + 2 * i]), (int(tones[1 + 2 * i]) - int(tones[0 + 2 * i]) + 1))
        for j in range(len(g)): tones_1.append(int(g[j]))

    ntones = np.append(np.asarray(ntones), tones_1)
    ntones.sort()
    ntones = ntones - 1
    
    return ntones


def pcal_difx(ifile1, ntone):
    global table1, counter, acc_periods

    ifile = open(ifile1)

    acc_periods = len((ifile).readlines()) - 5
    ifile.seek(0)
    
    for k in range(5): ifile.readline()

    accumulation_period = 0.5
    
    ntones = ntone_form(ntone)
    counter = len(ntones)

    table1 = (np.empty((counter, 0))).tolist()
    ph = []
    
    for j in range(acc_periods):
        smh = ((ifile.readline()).split())[6:]
        for i in range(counter):
            table1[i].append(cmath.phase(complex(float(smh[(len(smh) - 1) - 1 - int(ntones[i]) * 4]),float(smh[(len(smh) - 1) - int(ntones[i]) * 4]))))
            (table1[i])[j] = ((table1[i])[j] * (180 / np.pi))
            ph.append(complex(float(smh[2 + int(ntones[i]) * 4]), float(smh[3 + int(ntones[i]) * 4])))

    for k in range(counter): table1[k] = unwraping(unwraping(table1[k]))
    
    ifile.close()


def pcal_rasfx(ifile2, ntone):
    global table2, counter, acc_periods

    r = ifile2

    ifile = open(ifile2, 'rb')
    
    pcal_version = str(ifile.read(20))
    bandwidth, = struct.unpack('i', ifile.read(4))
    bandwidth = int(bandwidth * 1e-6)
    frequency_offset, = struct.unpack('i', ifile.read(4))
    number_of_channels, = struct.unpack('i', ifile.read(4))
    accumulation_period, = struct.unpack('f', ifile.read(4))
    acc, = struct.unpack('i', ifile.read(4))

    acc_periods = acc
    
    ph = []
    
    ntones = ntone_form(ntone)
    counter = len(ntones)

    for i in range(acc_periods):
        for j in range(bandwidth):
            el1, = struct.unpack('d', ifile.read(8))
            el2, = struct.unpack('d', ifile.read(8))

            if len(ntones) == bandwidth: ph.append(complex(el1, el2))
            else:
                for k in range(len(ntones)):
                    if j == ntones[k]: ph.append(complex(el1, el2))

    ph_table = (np.empty((int(counter), 0))).tolist()
    
    for i in range(counter):
        for j in range(acc_periods): ph_table[i].append((ph[i + counter * j]))

    table2 = (np.empty((int(counter), 0))).tolist()
    
    for i in range(counter):
        for j in range(acc_periods): table2[i].append((cmath.phase(((ph_table[i])[j]))) * (180 / np.pi))

    for i in range(counter):
    	table2[i].reverse()
    	table2[i] = unwraping(unwraping(table2[i]))

    ifile.close()


def pcal_diff(table1, table2):
	global diff_massive

	diff_massive = []

	for i in range(counter):
		for j in range(acc_periods): diff_massive.append(((table1[i])[j] - (table2[i])[j]))

	new_diff_massive = []

	if average == 'true':
		k = 0
		while k < (counter * acc_periods):
			new_diff_massive.append(np.mean(diff_massive[k : (k + acc_periods)]))
			k = k + acc_periods

		diff_massive = new_diff_massive
	

if __name__ == '__main__':
	main()

	if ifile1 != None:
		pcal_difx(ifile1, ntone)

		if write == 'true':
			name = ifile1 + '_phases'
			f = open(name, 'w')
			for k in range(counter):
				for i in range(acc_periods):
					f.write(str((table1[k])[i]) + '\n')

	if ifile2 != None:
		pcal_rasfx(ifile2, ntone)

		if write == 'true':
			name = ifile2 + '_phases'
			f = open(name, 'w')
			for k in range(counter):
				for i in range(acc_periods):
					f.write(str((table2[k])[i]) + '\n')

	if ifile1 != None and ifile2 != None:
		pcal_diff(table1, table2)

		if write == 'true' and average == 'false':
			name = ifile1 + '_phases_diff'
			f = open(name, 'w')
			for k in range(counter * acc_periods):
				f.write(str(diff_massive[k]))
				f.write('\n')
		elif write == 'true' and average == 'true':
			name = ifile1 + '_phases_diff'
			f = open(name, 'w')
			for k in range(counter):
				f.write(str(diff_massive[k]) + '\n')