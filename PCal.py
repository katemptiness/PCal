#!/usr/bin/env python

import numpy as np, numpy.fft as fft, numpy.linalg as linalg
import cmath, math
import getopt, sys
import matplotlib.pyplot as plt
import os
from itertools import count, izip
import struct

def file_read(ifile):
    if ifile[0:7] == 'PCAL_58':
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

    
def main():
    global itype, dbg, ntone, ifile, acc_period

    itype = 'phase'
    dbg = 'false'
    acc_period = None

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hf:n:t:d:a:', ['ifile=', 'ntone=', 'itype=', 'dbg=', 'acc_period='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-f', '--ifile'):
            if arg[0] == ' ':
                ifile = arg[1:]
            else:
                ifile = arg
                ntone = '1 : ' + str(file_read(ifile))
        elif opt in ('-n', '--ntone'):
            if arg[0] == ' ':
                ntone = arg[1:]
            else:
                ntone = arg
        elif opt in ('-t', '--itype'):
            if arg[0] == ' ':
                itype = arg[1:]
            else:
                itype = arg
        elif opt in ('-d', '--dbg'):
            if arg[0] == ' ':
                dbg = arg[1:]
            else:
                dbg = arg
        elif opt in ('-a', '--acc_period'):
            if arg[0] == ' ':
                acc_period = int(arg[1:])
            else:
                acc_period = int(arg)
            if acc_period == 1:
                print 'Sorry, -a parameter should be more than 1.'
                sys.exit()

    
def usage():
    print 'Hello, this is the USAGE function.'
    print
    print 'You can use this forms to work:'
    print '-f is for path to the file;'
    print '-n is for tone numbers;'
    print '-a is for accumulation periods;'
    print '-t is for data type (amplitudes or phases);'
    print '-d is for graphics display mode (true or false).'
    print
    print 'For example:'
    print 'python PCal.py -f W:/Files/My_File -n "1 : 512" -a "20" -t phase - d true'
    print
    print '-n, -a, -t & -d parameters are "1 : last", "all", "phase" & "false" as default, so you can only use -f.'
    print
    print 'Warning: -a parameter (accumulation periods) should be more than 1!'


def unwraping(lista):
    k = 0
    while k < (len(lista) - 1):
        if abs(lista[k] - lista[k + 1]) > 340:

            list1 = []
            list2 = []
        
            i = 0
            while i < len(lista):
                if lista[i] == abs(lista[i]):
                    if lista[i] > np.mean(lista):
                        list1.append(i)
                    else:
                        list2.append(i)
                elif lista[i] == -abs(lista[i]):
                    if lista[i] < np.mean(lista):
                        list2.append(i)
                    else:
                        list1.append(i)
                i = i + 1
    
            if len(list1) > len(list2):
                c = '+'
            elif len(list1) < len(list2):
                c = '-'
            elif len(list1) == len(list2):
                if lista[0] >= 0:
                    c = '+'
                elif lista[0] < 0:
                    c = '-'
    
            if c == '+':
                i = 0
                while i < (len(lista) - 1):
                    if (lista[i + 1] - lista[i]) < -180:
                        lista[i + 1] = lista[i + 1] + 360
                    elif (lista[i + 1] - lista[i]) > 180:
                        lista[i] = lista[i] + 360
                    i = i + 1
    
            elif c == '-':
                i = 0
                while i < (len(lista) - 1):
                    if (lista[i + 1] - lista[i]) < -340:
                        lista[i] = lista[i] - 360
                    elif (lista[i + 1] - lista[i]) > 340:
                        lista[i + 1] = lista[i + 1] - 360
                    i = i + 1

            if np.mean(lista) > 0:
                if lista[0] < 0 and lista[1] < 0:
                    lista[0] = lista[0] + 360
                    lista[1] = lista[1] + 360
            elif np.mean(lista) < 0:
                if lista[0] > 0 and lista[1] > 0:
                    lista[0] = lista[0] - 360
                    lista[1] = lista[1] - 360

        k = k + 1

    return lista


def unwraping2(lista):
    i = 0
    while i < (len(lista) - 1):
        r = lista[i + 1] - lista[i]
        if abs(r) > 360 and lista[i + 1] < lista[i]:
            lista[i + 1] = lista[i + 1] + 360
            
            j = i
            while j < len(lista):
                lista[j] = lista[j] + 360
                j = j + 1
            
        elif abs(r) > 360 and lista[i + 1] > lista[i]:
            lista[i + 1] = lista[i + 1] - 360
            
            j = i
            while j < len(lista):
                lista[j] = lista[j] - 360
                j = j + 1
    
        i = i + 1

    return lista


def SetIndex0(i, N, iPow):
    i0 = i
    
    bb = []
    
    kk = 0
    while kk < 32:
        bb.append('')
        kk = kk + 1
    
    j = 0
    while j < iPow:
        bb.append(False)
        j = j + 1

    j = N
    k = iPow
    l = 0
    m = 0

    while i >= 1:
        bb[k] = i >= j
        if bb[k]:
            i = i - j
            l = l + m
        j = j / 2
        k = k - 1
        if m == 0:
            m = 1
        else:
            m = m * 2

    return l


def Fraq_FFT(N, Re0, Im0, Tau, bInv):
    Re1 = Im1 = 0
    iPow = 0
    
    i = 1
    while i < N:
        iPow = iPow + 1
        i = i * 2
    
    if i > N:
        return
    
    acos = []
    asin = []
    
    kk = 0
    while kk < iPow:
        acos.append('')
        asin.append('')
        kk = kk + 1
    
    Re = []
    Im = []
    
    kk = 0
    while kk < 2 * N:
        Re.append('')
        Im.append('')
        kk = kk + 1
    
    I = []
    
    kk = 0
    while kk < N:
        I.append('')
        kk = kk + 1
    
    j = 0
    while j < N:
        I[j] = SetIndex0(j, N, iPow)
        j = j + 1
        
    a = np.pi * Tau
        
    if bInv == 0:
        a = a * (-1)
        
    i = 0
    while i < iPow:
        acos[i] = np.cos(a)
        asin[i] = np.sin(a)
        a = a / 2
        i = i + 1
        
    j = 0
    while j < N:
        j1 = I[j]
        Re[j1] = Re0[j]
        Im[j1] = Im0[j]
        j = j + 1
        
    k = 2
        
    l = 0
    while l < iPow:
        j = 0
        while j < N:
            j1 = j + k / 2
            Re[j] = Re[j] + Re[j1] * acos[l] - Im[j1] * asin[l]
            Im[j] = Im[j] + Re[j1] * asin[l] + Im[j1] * acos[l]
            j = j + k
        k = k * 2
        l = l + 1
        
    Re1 = Re[0] / N
    Im1 = Im[0] / N
        
    return Im1, Re1
    
    
def pcal_read(ifile, ntone, itype, dbg, acc_period):
    global table, table2, counter, ph, ph_table, acc_periods, ntones, accumulation_period
    
    ifile = open(ifile)

    if acc_period == None:
        acc_periods = len((ifile).readlines()) - 5
        ifile.seek(0)
    else:
        acc_periods = acc_period
    
    k = 1
    while k <= 5:
        ifile.readline()
        k = k + 1

    accumulation_period = 0.5
    
    ntones = (str(ntone).replace(',', '')).split()
    
    i = ntones.count(':')
    
    tones = []
    
    while i > 0:
        g = ntones.index(':')
        tones.append(ntones[g - 1])
        tones.append(ntones[g + 1])
        del ntones[g], ntones[g], ntones[g - 1]
        i = i - 1
    
    li = []
    i = 0
    while i < len(ntones):
        li.append(int(ntones[i]))
        i = i + 1
    ntones = li
    
    tones_1 = []
    i = 0
    
    while i < len(tones) / 2:
        g = np.linspace(int(tones[0 + 2 * i]), int(tones[1 + 2 * i]), (int(tones[1 + 2 * i]) - int(tones[0 + 2 * i]) + 1))
        j = 0
        while j <= (len(g) - 1):
            tones_1.append(int(g[j]))
            j = j + 1
        i = i + 1
    
    ntones = np.append(np.asarray(ntones), tones_1)
    ntones.sort()
    ntones = ntones - 1

    counter = len(ntones)

    table = (np.empty((counter, 0))).tolist()
    table2 = (np.empty((counter, 0))).tolist()
    ph_table = (np.empty((counter, 0))).tolist()

    ph = []
    
    j = 0
    while j < acc_periods:
        smh = ((ifile.readline()).split())[6:]
        i = 0
        while i < counter:
            if itype == 'phase':
                table[i].append(cmath.phase(complex(float(smh[(len(smh) - 1) - 1 - int(ntones[i]) * 4]), float(smh[(len(smh) - 1) - int(ntones[i]) * 4]))) * (180 / np.pi))
            elif itype == 'amplitude':
                table[i].append(math.hypot(float(smh[(len(smh) - 1) - 1 - int(ntones[i] * 4)]), float(smh[(len(smh) - 1) - int(ntones[i]) * 4])) * 1000)
            elif itype == 'phase-amplitude':
                table[i].append(cmath.phase(complex(float(smh[(len(smh) - 1) - 1 - int(ntones[i]) * 4]), float(smh[(len(smh) - 1) - int(ntones[i]) * 4]))) * (180 / np.pi))
                table2[i].append(math.hypot(float(smh[(len(smh) - 1) - 1 - int(ntones[i] * 4)]), float(smh[(len(smh) - 1) - int(ntones[i]) * 4])) * 1000)
            ph.append(complex(float(smh[2 + int(ntones[i]) * 4]), float(smh[3 + int(ntones[i]) * 4])))
            #ph_table[i].append(complex(float(smh[2 + int(ntones[i]) * 4]), float(smh[3 + int(ntones[i]) * 4])))
            i = i + 1
        j = j + 1

    ifile.close()
    
    if itype == 'phase-amplitude':
        return 
    else:
        return table
    

def pcal_reading(ifile, ntone, itype, dbg, acc_period):
    global table, counter, ph, ph_table, acc_periods, ntones, accumulation_period
    
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
    
    ntones = (str(ntone).replace(',', '')).split()
    
    i = ntones.count(':')
    
    tones = []
    
    while i > 0:
        g = ntones.index(':')
        tones.append(ntones[g - 1])
        tones.append(ntones[g + 1])
        del ntones[g], ntones[g], ntones[g - 1]
        i = i - 1
    
    li = []
    i = 0
    while i < len(ntones):
        li.append(int(ntones[i]))
        i = i + 1
    ntones = li
    
    tones_1 = []
    i = 0
    
    while i < len(tones) / 2:
        g = np.linspace(int(tones[0 + 2 * i]), int(tones[1 + 2 * i]), (int(tones[1 + 2 * i]) - int(tones[0 + 2 * i]) + 1))
        j = 0
        while j <= (len(g) - 1):
            tones_1.append(int(g[j]))
            j = j + 1
        i = i + 1

    ntones = np.append(np.asarray(ntones), tones_1)
    ntones.sort()
    ntones = ntones - 1

    counter = len(ntones)

    ph = []

    i = 0
    while i < acc_periods:
        j = 0
        while j < file_read(r):
            el1, = struct.unpack('f', ifile.read(4))
            el2, = struct.unpack('f', ifile.read(4))

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
            ph_table[i].append(ph[i + counter * j])
            j = j + 1
        i = i + 1

    table = (np.empty((int(counter), 0))).tolist()
    
    i = 0
    while i < counter:
        j = 0
        while j < acc_periods:
            table[i].append(cmath.phase((((ph_table[i])[j]))) * (180 / np.pi))
            j = j + 1
        i = i + 1

    table.reverse()

    ifile.close()

    return table


def pcal_trend(ifile, ntones, itype, dbg):
    global trends, std, new_table

    if dbg == 'true' and what == '2':
        f, axar = plt.subplots(4)

    time = np.linspace(0, accumulation_period * acc_periods, acc_periods)

    trends = []

    i = 0
    while i < counter:
        A = (np.vstack([time, np.ones(len(time))])).transpose()
        m, c = linalg.lstsq(A, unwraping(table[i - 1]))[0]
        trend = m * time + c
        trends.append(trend)
                
        if dbg == 'true':
            axar[0].plot(time, unwraping(table[i - 1]), 'o')
            axar[0].plot(time, trend)
        
        i = i + 1

    AC = len(time)
    
    alphas = []

    re_trends = []
    re_table = []
    new_table = []
    
    std = []
    
    i = 0
    while i < counter:
        j = 0
        while j < acc_periods:
            re = (trends[i])[1] - (trends[i])[0]
            re_trends.append((trends[i])[j] - j * re)
            re_table.append(unwraping(table[i])[j] - j * re)
            j = j + 1
        std.append(np.std(unwraping(re_table)))
        new_table.append(re_table)
        re_table = []
        i = i + 1
            
    if dbg == 'true' and what == '2':
        axar[0].grid()
        axar[0].set_xlabel('time')
        axar[0].set_ylabel(itype)

        j = 0
        while j < counter:
            BC = (trends[j])[acc_periods - 1] - (trends[j])[0]
            AB = np.sqrt(BC * BC + AC * AC)
            alpha = np.arcsin(BC / AB) * (180 / np.pi)
            axar[1].plot((int(ntones[j]) + 1), alpha, 'o')
            j = j + 1

        axar[1].grid()
        axar[1].set_xlabel('tone numbers')
        axar[1].set_ylabel('tilt angle')

        j = 0
        while j < counter:
            axar[2].plot((ntones[j] + 1), std[j], 'o')
            j = j + 1

        axar[2].grid()
        axar[2].set_xlabel('tone numbers')
        axar[2].set_ylabel('standard deviation')
        
        axar[3].hist(std, bins = counter)
        axar[3].set_xlabel('standard deviation')
        axar[3].set_ylabel('tones')

        plt.show()


def pcal_phaseresponse(ifile, ntones, itype, dbg):
    global delay

    pcal_trend(ifile, ntones, itype, dbg)
    
    li = []
    
    std_threshold = 4 * min(std)

    good_table = []
    good_ntones = []

    j = 0
    while j < counter:
        if std[j] < std_threshold:
            good_table.append(new_table[j])
            good_ntones.append(j)
        j = j + 1

    good_ntones = np.asarray(good_ntones)

    if what == '2':
        i = 0
        while i < len(good_ntones):
            li.append(np.mean(good_table[i]))
            i = i + 1

        plt.plot(good_ntones, unwraping2(li))
        plt.plot(good_ntones, unwraping2(li), 'o')
    
        if dbg == 'true':
            plt.grid()
            plt.xlabel('frequency')
            if itype == 'phase-amplitude':
                plt.ylabel('phase')
            else:
                plt.ylabel(itype)
            plt.show()


def pcal_delay(ifile, ntones, itype, dbg):
    if itype == 'phase-amplitude':
        plt.plot(table, table2, 'o')
        plt.grid()
        plt.xlabel('phase')
        plt.ylabel('amplitude')
        plt.show()
    
    else:
        #############################################################

        #good_ntones = pcal_phaseresponse(ifile, ntones, itype, dbg)
        #ph_table_new = []

        #i = 0
        #while i < len(good_ntones):
            #ph_table_new.append(ph_table[good_ntones[i]])
            #i = i + 1

        #ph_new = []

        #j = 0
        #while j < len(good_ntones):
            #i = 0
            #while i < acc_periods:
                #ph_new.append((ph_table_new[j])[i])
                #i = i + 1
            #j = j + 1
        #############################################################

        li = []

        f, axar = plt.subplots(2)

        time = np.linspace((1 / counter), 1, counter)

        j = 0
        while j < acc_periods:
            ph1 = abs(fft.ifft(ph[(j * counter) : (j * counter + counter)]))
            axar[0].plot(time, ph1)

            number = max(izip(ph1, count()))[1]

            j0 = number

            if j == 0:
                j1 = ("%.6f" % (((j0 * 1e-6) / 512) * 1e6))
                print 'The time delay is probably', j1, 'microseconds \n'
                print 'Starting the calculation...'
            
            tau_min = j0 - 1
            tau_max = j0 + 1

            delta_tau = 0.1
            while delta_tau >= 1e-3:
                tau = tau_min
            
                tau_list = []
                cj = []
            
                while tau <= tau_max:
                    im, re = Fraq_FFT(file_read(ifile), (np.asarray(ph).real)[(j * counter) : (j * counter + counter)], (np.asarray(ph).imag)[(j * counter) : (j * counter + counter)], tau, 1)
                
                    cj.append(complex(re, im))
                    tau_list.append(tau)

                    tau = tau + delta_tau
                
                cj = abs(np.asarray(cj))

                number = max(izip(cj, count()))[1]
                tau = tau_list[number]

                tau_min = tau - delta_tau * 2
                tau_max = tau + delta_tau * 2
                delta_tau = delta_tau / 10
        
            tau = tau / 512
            tau = float("%.6f" % (tau))

            sys.stdout.write(' accumulation periods have been processed...' + ('\r%d'%(j + 1) + '/' + str(acc_periods)))
            sys.stdout.flush()

            li.append(tau)

            j = j + 1

        tau = "%.6f" % (np.mean(li))

        print '\nAnd the clarified time delay is', tau, 'microseconds'
        
        if dbg == 'true':
            axar[0].grid()
            axar[0].set_xlabel('time')
            axar[0].set_ylabel('amplitude')

            xlist = np.linspace(1, acc_periods, acc_periods)
            axar[1].cla()
            axar[1].plot(xlist, li)
            axar[1].plot(xlist, li, 'o')
            axar[1].grid()
            axar[1].set_xlabel('accumulation periods')
            axar[1].set_ylabel('time delay')

            plt.show()


if __name__ == '__main__':
    main()
    print 'Hello. Welcome to PCal interface. Please wait until your file will be ready.'
    
    where_to_go = open(ifile)
    if (where_to_go.readline())[0] == '#':
        pcal_read(ifile, ntone, itype, dbg, acc_period)
    else:
        pcal_reading(ifile, ntone, itype, dbg, acc_period)

    print '\nNow tell me what you want to do:'
    print 'press 1 if uou want to plot signal and see the time delay;'
    print 'press 2 if you want to plot tilt angle and STD graphics and see phase-frequency response.'
    
    global what
    what = raw_input()
    if what == '2':
        pcal_phaseresponse(ifile, ntones, itype, dbg)
    elif what == '1':
        pcal_delay(ifile, ntones, itype, dbg)
    else:
        print 'Error, please try again.'