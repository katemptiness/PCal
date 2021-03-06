#!/usr/bin/env python

import numpy as np, numpy.fft as fft, numpy.linalg as linalg
import cmath, math
import getopt, sys
import matplotlib.pyplot as plt
import os
from itertools import count, izip
import struct

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

    
def main():
    global itype, dbg, ntone, ifile, acc_period, write

    itype = 'phase'
    dbg = 'false'
    acc_period = None
    write = 'false'

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hf:n:t:d:a:w:', ['ifile=', 'ntone=', 'itype=', 'dbg=', 'acc_period=', 'write='])
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
                ntone = '1 : ' + str(file_read(ifile))
        elif opt in ('-n', '--ntone'):
            ntone = arg
        elif opt in ('-t', '--itype'):
            itype = arg
        elif opt in ('-d', '--dbg'):
            dbg = arg
        elif opt in ('-w', '--write'):
            write = arg
        elif opt in ('-a', '--acc_period'):
            acc_period = int(arg)
            if acc_period == 1:
                print 'Sorry, -a parameter should be more than 1.'
                sys.exit()

    
def usage():
    print 'Hello, this is the USAGE function.'
    print
    print 'You can use this forms to work:'
    print '-f is for path to the file or directory;'
    print '-n is for tone numbers;'
    print '-a is for accumulation periods;'
    print '-t is for data type (amplitudes or phases);'
    print '-d is for graphics display mode (true or false);'
    print '-w is for delays recording (true or false).'
    print
    print 'For example:'
    print 'python PCal.py -f W:/Files/My_File -n "1 : 512" -a "20" -t phase - d true -w true'
    print
    print '-n, -a, -t, -d & -w parameters are "1 : last", "all", "phase", "false" & "false" as default, so you can only use -f.'
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
    
    bb = (np.empty((32, 0))).tolist()
    
    j = 0
    while j < iPow:
        bb.append(False)
        j = j + 1

    j = N
    k = iPow
    l = m = 0

    while i >= 1:
        if i >= j:
            bb[k] = True
        if bb[k] == True:
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
    
    acos = (np.empty((iPow, 0))).tolist()
    asin = (np.empty((iPow, 0))).tolist()

    Re = (np.empty((2 * N, 0))).tolist()
    Im = (np.empty((2 * N, 0))).tolist()
    
    I = (np.empty((N, 0))).tolist()
    
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

    ph = []
    
    j = 0
    while j < acc_periods:
        smh = ((ifile.readline()).split())[6:]
        i = 0
        while i < counter:
            if itype == 'phase':
                table[i].append(cmath.phase(complex(float(smh[(len(smh) - 1) - 1 - int(ntones[i]) * 4]), float(smh[(len(smh) - 1) - int(ntones[i]) * 4]))))
                (table[i])[j] = (table[i])[j] * (180 / np.pi)
            elif itype == 'amplitude':
                table[i].append(math.hypot(float(smh[(len(smh) - 1) - 1 - int(ntones[i] * 4)]), float(smh[(len(smh) - 1) - int(ntones[i]) * 4])) * 1000)
            ph.append(complex(float(smh[2 + int(ntones[i]) * 4]), float(smh[3 + int(ntones[i]) * 4])))
            i = i + 1
        j = j + 1

    ifile.close()
    
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
            
            if itype == 'phase':
                table[i].append(cmath.phase((((ph_table[i])[j]))) * (180 / np.pi))
            elif itype == 'amplitude':
                table[i].append(math.hypot(((ph_table[i])[j]).real, ((ph_table[i])[j]).imag))

            j = j + 1
        i = i + 1

    table.reverse()

    ifile.close()

    return table


def pcal_phaseresponse(ifile, ntones, itype, dbg):
    global trends, std, new_table

    if dbg == 'true':
        f, axar = plt.subplots(4)

    time = np.linspace(0, (accumulation_period * acc_periods - accumulation_period), acc_periods)
    trends = []
    i = 0
    while i < counter:
        A = (np.vstack([time, np.ones(len(time))])).transpose()
        m, c = linalg.lstsq(A, unwraping(unwraping(table[i - 1])), rcond = -1)[0]
        trend = m * time + c
        trends.append(trend)
                
        if dbg == 'true':
            plt.figure(1)
            axar[0].plot(time, unwraping(unwraping(table[i - 1])), 'o')
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
            re_table.append(unwraping(unwraping(table[i]))[j] - j * re)
            j = j + 1
        std.append(np.std(unwraping(unwraping(re_table))))
        new_table.append(re_table)
        re_table = []
        i = i + 1
            
    if dbg == 'true':
        plt.figure(1)
        axar[0].grid()
        axar[0].set_xlabel('time, s')
        axar[0].set_ylabel('phase, grad')

        j = 0
        while j < counter:
            BC = (trends[j])[acc_periods - 1] - (trends[j])[0]
            AB = np.sqrt(BC * BC + AC * AC)
            alpha = np.arcsin(BC / AB) * (180 / np.pi)
            plt.figure(1)
            axar[1].plot((int(ntones[j]) + 1), alpha, 'o')
            j = j + 1

        plt.figure(1)
        axar[1].grid()
        axar[1].set_xlabel('tone numbers')
        axar[1].set_ylabel('tilt angle, grad')

        j = 0
        while j < counter:
            plt.figure(1)
            axar[2].plot((ntones[j] + 1), std[j], 'o')
            j = j + 1

        plt.figure(1)
        axar[2].grid()
        axar[2].set_xlabel('tone numbers')
        axar[2].set_ylabel('standard deviation, grad')
        
        axar[3].hist(std, bins = counter)
        axar[3].set_xlabel('standard deviation, grad')
        axar[3].set_ylabel('tones')

        plt.gcf().canvas.set_window_title('Phase of time graph, tilt angle & STD')

        plt.show(block = False)

    #li = []

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

    #i = 0
    #while i < len(good_ntones):
        #li.append(np.mean(good_table[i]))
        #i = i + 1
        
    #plt.figure(2)
    #plt.plot(good_ntones, unwraping2(li))
    #plt.plot(good_ntones, unwraping2(li), 'o')
    
    if dbg == 'true':
        #plt.figure(2)
        #plt.grid()
        #plt.xlabel('frequency')
        #plt.ylabel(itype)
        #plt.gcf().canvas.set_window_title('Phase-frequency responce')
        plt.show()

    new_ph = []

    j = 0
    while j < acc_periods:
        i = 0
        k = 0
        while i < counter:
            try:
                if i == good_ntones[k]:
                    new_ph.append((ph[(j * counter) : ((j + 1) * counter)])[i])
                    k = k + 1
                else:
                    new_ph.append(0)
            except:
                new_ph.append(0)
            i = i + 1
        j = j + 1

    return new_ph


def pcal_delay(ifile, ntones, itype, dbg, qwerty, write):
    new_ph = pcal_phaseresponse(ifile, ntones, itype, 'false')
    ph = new_ph

    li = []
    
    if dbg == 'true':
        f, axar = plt.subplots(2)
        axar[0].grid()

    time = np.linspace((1 / counter), 1, counter)

    ampls = []
    xlist = []

    j = 0
    while j < acc_periods:
        ph1 = abs(fft.ifft(ph[(j * counter) : (j * counter + counter)]))
            
        if dbg == 'true':
            axar[0].plot(time, ph1)
            plt.pause(0.001)
            axar[0].set_xlabel('time, us')
            axar[0].set_ylabel('amplitude, V')
                
        number = max(izip(ph1, count()))[1]

        #ampl = ph1[number]
            
        #noise = []
        #if number >= 250:
            #i = 0
            #while i < (number - 100):
                #noise.append(ph1[i])
                #i = i + 1
        #elif number < 250:
            #i = (file_read(ifile) - 1)
            #while i > (number + 100):
                #noise.append(ph1[i])
                #i = i - 1

        #m = np.mean(noise)
        #snr = ((ampl - m) / np.std(np.asarray(noise) - m))
        #sigma = (np.sqrt(12) / (2 * np.pi * file_read(ifile) * 1e6 * snr))

        j0 = number

        if j == 0:
            j1 = ("%.7f" % (((j0 * 1e-6) / 512) * 1e6))
            print '\nThe time delay is about', j1, 'microseconds \n'
            print 'Starting the calculation...'
            
        tau_min = j0 - 1
        tau_max = j0 + 1
            
        delta_tau = 0.1
        while delta_tau >= 5e-5:
            tau = tau_min
            
            tau_list = []
            cj = []

            res = (np.asarray(ph).real)[(j * counter) : (j * counter + counter)]
            ims = (np.asarray(ph).imag)[(j * counter) : (j * counter + counter)]
                
            while tau <= tau_max:
                im, re = Fraq_FFT(file_read(ifile), res, ims, tau, 1)
                
                cj.append(complex(re, im))
                tau_list.append(tau)

                tau = tau + delta_tau
            
            cj = abs(np.asarray(cj))
            
            number = max(izip(cj, count()))[1]
            tau = tau_list[number]
            tau_min = tau - delta_tau * 2
            tau_max = tau + delta_tau * 2
                
            if delta_tau == 1e-4:
                delta_tau = delta_tau / 2
            else:
                delta_tau = delta_tau / 10

        tau = tau / 512
        tau = float("%.7f" % (tau))

        sys.stdout.write(' accumulation periods have been processed...' + ('\r%d'%(j + 1) + '/' + str(acc_periods)))
        sys.stdout.flush()

        li.append(tau)

        if dbg == 'true':
            xlist.append(j * 0.5)

            axar[1].cla()
                
            plt.pause(0.001)
            
            axar[1].plot(xlist, np.asarray(li) * 1e6, 'o')
            axar[1].grid()
            axar[1].set_xlabel('time, s')
            axar[1].set_ylabel('time delay, ps')
            plt.gcf().canvas.set_window_title(ifile)
            plt.draw()
            
        j = j + 1
        
    tau = "%.7f" % (np.mean(li))

    print '\nAnd the clarified time delay is', tau, 'microseconds'

    if dbg == 'true':
        plt.gcf().canvas.set_window_title(ifile)
        if qwerty == '0':
            plt.show(block = False)
        elif qwerty == '1':
            plt.show()

    if write == 'true':
        name = ifile + '_delays.txt'
        f = open(name, 'w')
        k = 0
        while k < len(li):
            f.write(str(li[k]))
            f.write('\n')
            k = k + 1

    return li


def pcal_diff(a, b):
    if len(a) > len(b):
        d = len(a) - len(b)
        while d > 0:
            del a[-1]
            d = d - 1
    elif len(b) > len(a):
        d = len(b) - len(a)
        while d > 0:
            del b[-1]
            d = d - 1

    diff = abs(np.asarray(a) - np.asarray(b)) * 1e6

    xlist = np.linspace(1, len(a), len(a))
    xlist = np.linspace(0, (len(a) / 2 - 0.5), len(a))

    A = (np.vstack([xlist, np.ones(len(xlist))])).transpose()
    m, c = linalg.lstsq(A, diff, rcond = -1)[0]
    trend = m * xlist + c
    
    print '\nThe slope of trend line is', "%.3f" % (m)

    if dbg == 'true':
        plt.figure(3)
        plt.plot(xlist, diff, 'o')
        plt.plot(xlist, trend)
        plt.grid()
        plt.xlabel('time, s')
        plt.ylabel('difference between time delays, ps')
        plt.gcf().canvas.set_window_title('Difference between time delays')
        plt.show()


if __name__ == '__main__':
    global what, qwerty, files
    main()

    if os.path.exists(ifile):

        if os.path.isfile(ifile):

            print 'Hello. Welcome to PCal. Please wait until your file will be ready.'

            where_to_go = open(ifile)
            if (where_to_go.readline())[0] == '#':
                pcal_read(ifile, ntone, itype, dbg, acc_period)
            else:
                pcal_reading(ifile, ntone, itype, dbg, acc_period)

            print '\nNow tell me what you want to do:'
            print 'press 1 if you want to plot signal and see the time delay;'
            print 'press 2 if you want to plot tilt angle and STD graphics and see phase-frequency response;'
            print 'press 3 if you want to see the difference between time delays in 2 different files.'
            print 'press 4 if you want to leave.'
    
            what = raw_input()
            if what == '2':
                pcal_phaseresponse(ifile, ntones, itype, dbg)
            elif what == '1':
                qwerty = '1'
                pcal_delay(ifile, ntones, itype, dbg, qwerty, write)
            elif what == '3':
                qwerty = '0'
                a = pcal_delay(ifile, ntones, itype, dbg, qwerty, write)
        
                print '\nPlease enter the 2nd file...'
                new_ifile = raw_input()
        
                where_to_go = open(ifile)
                if (where_to_go.readline())[0] == '#':
                    pcal_read(new_ifile, ntone, itype, dbg, acc_period)
                else:
                    pcal_reading(new_ifile, ntone, itype, dbg, acc_period)
        
                b = pcal_delay(new_ifile, ntones, itype, dbg, qwerty, write)

                pcal_diff(a, b)

            elif what == '4':
                sys.exit()

            else:
                print 'Error, please try again.'

        elif os.path.isdir(ifile):
            print 'Hello. Welcome to PCal.'
            print 'Press 1 if you want to calculate time delays and see differences between them.'
            print 'Press 2 if you want to leave.'
            mode = raw_input()

            if mode == '1':
                qwerty = '0'
                
                all_files = os.listdir(ifile)

                files = []

                i = 0
                while i < len(all_files):
            
                    if (all_files[i])[0:5] == 'PCAL_' and (all_files[i])[-1] == 'V':
                        files.append(ifile + all_files[i])
            
                    i = i + 1

                files = sorted(files)

                try:
                    i = 0
                    while i < (len(files) - 1):
                        pcal_read(files[i], ntone, itype, dbg, acc_period)
                        a = pcal_delay(files[i], ntones, itype, dbg, qwerty, write)

                        pcal_read(files[i + 1], ntone, itype, dbg, acc_period)
                        b = pcal_delay(files[i + 1], ntones, itype, dbg, qwerty, write)

                        if (files[i])[0 : 17] == (files[i + 1])[0 : 17]:
                            pcal_diff(a, b)

                        i = i + 2
        
                except NameError:
                    ntone = '1 : 512'
            
                    i = 0
                    while i < (len(files) - 1):
                        pcal_read(files[i], ntone, itype, dbg, acc_period)
                        poiuy = '2'
                        a = pcal_delay(files[i], ntones, itype, dbg, qwerty, write)

                        pcal_read(files[i + 1], ntone, itype, dbg, acc_period)
                        poiuy = '3'
                        b = pcal_delay(files[i + 1], ntones, itype, dbg, qwerty, write)

                        if (files[i])[0 : 17] == (files[i + 1])[0 : 17]:
                            pcal_diff(a, b)

                        i = i + 2

            elif mode == '2':
                sys.exit()

    else:
        print 'Sorry, file or directory was not found. Please try again'