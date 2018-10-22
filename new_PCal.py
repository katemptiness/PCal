#!/usr/bin/env python

import numpy as np, numpy.fft as fft, numpy.linalg as linalg
import cmath
import getopt, sys
import os
import struct

def main():
    global ntone, ifile, acc_period

    ifile = None
    acc_period = None

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hf:n:a:', ['ifile=', 'ntone=', 'acc_period='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-f', '--ifile'):
            ifile = arg
            if os.path.isfile(ifile):
                ntone = '1 : 512'
        elif opt in ('-n', '--ntone'): ntone = arg
        elif opt in ('-a', '--acc_period'):
            acc_period = int(arg)
            if acc_period == 1:
                print 'Sorry, -a parameter should be more than 1.'
                sys.exit()
        

def usage():
    print 'You can use this forms to work:'
    print '-f is for path to the file or directory;'
    print '-n is for tone numbers;'
    print '-a is for accumulation periods.'


def file_read(ifile):
    where_to_go = open(ifile)
    if (where_to_go.readline())[0] == '#': 
        ifile = open(ifile)

        i = 0
        while i < 5:
            ifile.readline()
            i = i + 1

        counter = int((str(ifile.readline()).split())[5])

        return counter
    
    else:
        ifile = open(ifile, 'rb')
    
        pcal_version = str(ifile.read(20))
        bandwidth, = struct.unpack('i', ifile.read(4))
        bandwidth = int(bandwidth * 1e-6)
    
        return bandwidth


def ntone_form(ntone):
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


def pcal_difx(ifile, ntone, acc_period):
    global table

    ifile = open(ifile)

    if acc_period == None:
        acc_periods = len((ifile).readlines()) - 5
        ifile.seek(0)
    else:
        acc_periods = acc_period
    
    for k in range(5): ifile.readline()

    accumulation_period = 0.5
    
    ntones = ntone_form(ntone)
    counter = len(ntones)

    table = (np.empty((counter, 0))).tolist()
    ph = []
    
    for j in range(acc_periods):
        smh = ((ifile.readline()).split())[6:]
        for i in range(counter):
            table[i].append(cmath.phase(complex(float(smh[(len(smh) - 1) - 1 - int(ntones[i]) * 4]), float(smh[(len(smh) - 1) - int(ntones[i]) * 4]))))
            (table[i])[j] = (table[i])[j] * (180 / np.pi)
            ph.append(complex(float(smh[2 + int(ntones[i]) * 4]), float(smh[3 + int(ntones[i]) * 4])))

    ifile.close()


def pcal_rasfx(ifile, ntone, acc_period):
    global table

    r = ifile

    ifile = open(ifile, 'rb')
    
    pcal_version = str(ifile.read(20))
    bandwidth, = struct.unpack('i', ifile.read(4))
    bandwidth = int(bandwidth * 1e-6)
    frequency_offset, = struct.unpack('i', ifile.read(4))
    number_of_channels, = struct.unpack('i', ifile.read(4))
    accumulation_period, = struct.unpack('f', ifile.read(4))
    acc, = struct.unpack('i', ifile.read(4))

    if acc_period == None:
        acc_periods = acc
    else:
        acc_periods = acc_period
    
    ph = []
    
    ntones = ntone_form(ntone)
    counter = len(ntones)

    i = 0
    while i < acc_periods:
        j = 0
        while j < file_read(r):
            el1, = struct.unpack('d', ifile.read(8))
            el2, = struct.unpack('d', ifile.read(8))

            if len(ntones) == file_read(r):
                ph.append(complex(el1, el2))
            else:
                k = 0
                while k < len(ntones):
                    if j == ntones[k]:
                        ph.append(complex(el1, el2))
                    k = k + 1

            j = j + 1
        i = i + 1

    ph_table = (np.empty((int(counter), 0))).tolist()
    
    i = 0
    while i < counter:
        j = 0
        while j < acc_periods:
            ph_table[i].append((ph[i + counter * j]))
            j = j + 1
        i = i + 1

    table = (np.empty((int(counter), 0))).tolist()
    
    i = 0
    while i < counter:
        j = 0
        while j < acc_periods:
            table[i].append((cmath.phase(((ph_table[i])[j]))) * (180 / np.pi))
            j = j + 1
        i = i + 1
    table.reverse()

    ifile.close()
    
    return table


if __name__ == '__main__':
	main()

	try:
		pcal_difx(ifile, ntone, acc_period)
	except:
		pcal_rasfx(ifile, ntone, acc_period)

	name = ifile + '_phases'
	f = open(name, 'w')
	for k in range(512):
		for i in range(42):
			f.write(str((table[k])[i]))
			f.write('\n')