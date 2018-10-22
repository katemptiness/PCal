#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np, numpy.fft as fft, numpy.linalg as linalg
import cmath, math
import getopt, sys
import os
from itertools import count, izip

def main():
    global dbg, ntone, ifile, acc_period, write, mode_filtration

    ifile = None
    dbg = 'false'
    acc_period = None
    write = 'false'
    mode_filtration = None

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hf:n:d:a:w:m:', ['ifile=', 'ntone=', 'dbg=', 'acc_period=', 'write=', 'mode_filtration='])
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
                #ntone = '1 : ' + str(file_read(ifile))
                ntone = '1 : 512'
        elif opt in ('-n', '--ntone'): ntone = arg
        elif opt in ('-d', '--dbg'): dbg = arg
        elif opt in ('-a', '--acc_period'):
            acc_period = int(arg)
            if acc_period == 1:
                print 'Sorry, -a parameter should be more than 1.'
                sys.exit()
        elif opt in ('-w', '--write'): write = arg
        elif opt in ('-m', '--mode_filtration'): mode_filtration = float(arg)

    
def usage():
    print 'You can use this forms to work:'
    print '-f is for path to the file or directory;'
    print '-n is for tone numbers;'
    print '-a is for accumulation periods;'
    print '-d is for graphics display mode (true or false);'
    print '-w is for delays recording (true or false);'
    print '-m is for tones filtration (true or false).'
    print
    print 'For example:'
    print 'python PCal.py -f W:/Files/My_File -n "1 : 512" -a "20" - d true -w true -m true'
    print
    print '-n, -a, -d, -w & -m parameters are "1 : last", "all", "false", "false" & "false" as default, so you can only use -f.'
    print
    print 'Warning: -a parameter (accumulation periods) should be more than 1!'


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


def unwraping(lista):
    for k in range(len(lista) - 1):
        if abs(lista[k] - lista[k + 1]) > 330:
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
                    if (lista[i + 1] - lista[i]) < -150: lista[i + 1] = lista[i + 1] + 360
                    elif (lista[i + 1] - lista[i]) > 150: lista[i] = lista[i] + 360
    
            elif c == '-':
                for i in range(len(lista) - 1):
                    if (lista[i + 1] - lista[i]) < -330: lista[i] = lista[i] - 360
                    elif (lista[i + 1] - lista[i]) > 330: lista[i + 1] = lista[i + 1] - 360

            if np.mean(lista) > 0:
                if lista[0] < 0 and lista[1] < 0:
                    lista[0] = lista[0] + 360
                    lista[1] = lista[1] + 360
            elif np.mean(lista) < 0:
                if lista[0] > 0 and lista[1] > 0:
                    lista[0] = lista[0] - 360
                    lista[1] = lista[1] - 360

    return lista


def SetIndex0(i, N, iPow):
    i0 = i
    
    bb = (np.empty((32, 0))).tolist()

    for j in range(iPow): bb.append(False)

    j = N
    k = iPow
    l = m = 0

    while i >= 1:
        if i >= j: bb[k] = True
        if bb[k] == True:
            i = i - j
            l = l + m
        j = j / 2
        k = k - 1
        if m == 0: m = 1
        else: m = m * 2

    return l


def Fraq_FFT(N, Re0, Im0, Tau, bInv):
    Re1 = Im1 = iPow = 0
    
    i = 1
    while i < N:
        iPow = iPow + 1
        i = i * 2
    
    acos = (np.empty((iPow, 0))).tolist()
    asin = (np.empty((iPow, 0))).tolist()
    Re = (np.empty((2 * N, 0))).tolist()
    Im = (np.empty((2 * N, 0))).tolist()
    I = (np.empty((N, 0))).tolist()
    
    for j in range(N): I[j] = SetIndex0(j, N, iPow)
    
    a = np.pi * Tau
        
    if bInv == 0: a = a * (-1)
    
    for i in range(iPow):
        acos[i] = np.cos(a)
        asin[i] = np.sin(a)
        a = a / 2

    for j in range(N):
        j1 = I[j]
        Re[j1] = Re0[j]
        Im[j1] = Im0[j]
    
    k = 2
    for l in range(iPow):
        j = 0
        while j < N:
            j1 = j + k / 2
            Re[j] = Re[j] + Re[j1] * acos[l] - Im[j1] * asin[l]
            Im[j] = Im[j] + Re[j1] * asin[l] + Im[j1] * acos[l]
            j = j + k
        k = k * 2

    Re1 = Re[0] / N
    Im1 = Im[0] / N
        
    return Im1, Re1
    
    
def pcal_read(ifile, ntone, acc_period):
    global table, counter, ph, acc_periods, ntones, accumulation_period
    
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


def pcal_reading(ifile, ntone, acc_period):
    global table, counter, ph, acc_periods, ntones, accumulation_period
    
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


def pcal_phaseresponse(ifile, ntones, dbg, mode_filtration, write):
    global trends, std, new_table, table

    time = np.linspace(0, (accumulation_period * acc_periods - accumulation_period), acc_periods)

    trends = []

    if dbg == 'true':
        plt.figure(1)
        plt.gcf().canvas.set_window_title(u'Фаза от времени')
    
    A = (np.vstack([time, np.ones(len(time))])).transpose()

    for i in range(counter):
        m, c = linalg.lstsq(A, unwraping(unwraping(table[i - 1])), rcond = -1)[0]
        trend = m * time + c
        trends.append(trend)
        
        if dbg == 'true':        
            plt.plot(time, unwraping(unwraping(table[i - 1])), 'o')
            plt.plot(time, trend)

    AC = len(time)
    
    alphas = []

    re_trends = []
    re_table = []
    new_table = []
    
    std = []
            
    for i in range(counter):
        for j in range(acc_periods):
            re = (trends[i])[1] - (trends[i])[0]
            re_trends.append((trends[i])[j] - j * re)
            re_table.append(unwraping(unwraping(table[i]))[j] - j * re)
        std.append(np.std(unwraping(unwraping(re_table))))
        new_table.append(re_table)
        re_table = []

    if dbg == 'true':
        plt.grid()
        plt.xlabel(u'время, с')
        plt.ylabel(u'фаза, градусы')
        plt.show(block = False)
        
        plt.figure(2)
        plt.gcf().canvas.set_window_title(u'Угол наклона от номера тона')

        for j in range(counter):
            BC = (trends[j])[-1] - (trends[j])[0]
            AB = np.sqrt(BC * BC + AC * AC)
            alpha = np.arcsin(BC / AB) * (180 / np.pi)
            plt.plot((int(ntones[j]) + 1), alpha, 'o')

        plt.grid()
        plt.xlabel(u'номера тонов, ед.')
        plt.ylabel(u'угол наклона, градусы')
        plt.show(block = False)
        
        f, axar = plt.subplots(2)
        plt.gcf().canvas.set_window_title(u'СКО')
        
        for j in range(counter): axar[0].plot((ntones[j] + 1), std[j], 'o')
        
        axar[0].grid()
        axar[0].set_xlabel(u'номера тонов, ед.')
        axar[0].set_ylabel(u'СКО, градусы')
        axar[1].hist(std, bins = counter)
        axar[1].set_xlabel(u'СКО, градусы')
        axar[1].set_ylabel(u'число тонов, ед.')
        plt.show()

    if write == 'true':
        name = ifile + '_ph.txt'
        f = open(name, 'w')
        for k in range(len(ph)):
            f.write(str((cmath.phase(ph[k])) * (180 / np.pi)))
            f.write('\n')

    if mode_filtration != None:
        std_threshold = mode_filtration

        good_ntones = []

        for j in range(counter):
            if std[j] < std_threshold: good_ntones.append(j)

        good_ntones = np.asarray(good_ntones)

        new_ph = []

        for j in range(acc_periods):
            k = 0
            for i in range(counter):
                try:
                    if i == good_ntones[k]:
                        new_ph.append((ph[(j * counter) : ((j + 1) * counter)])[i])
                        k = k + 1
                    else:
                        new_ph.append(0)
                except:
                    new_ph.append(0)

        return new_ph


def pcal_delay(ifile, ntones, dbg, qwerty, write, mode_filtration):
    global ph

    if mode_filtration == 'true':
        new_ph = pcal_phaseresponse(ifile, ntones, 'false', mode_filtration)
        ph = new_ph

    li = []
    
    if dbg == 'true':
        f, axar = plt.subplots(2)
        axar[0].grid()

    time = np.linspace((1 / counter), 1, counter)

    ampls = []
    xlist = []

    if dbg == 'true':
        axar[0].set_xlabel(u'время, мкс')
        axar[0].set_ylabel(u'амплитуда, усл.ед.')
    
    for j in range(acc_periods):
        ph1 = abs(fft.ifft(ph[(j * counter) : (j * counter + counter)]))
        
        if dbg == 'true':
            axar[0].plot(time, ph1)
            plt.pause(0.001)
                
        number = max(izip(ph1, count()))[1]
        j0 = number
        
        res = (np.asarray(ph).real)[(j * counter) : (j * counter + counter)]
        ims = (np.asarray(ph).imag)[(j * counter) : (j * counter + counter)]

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
            axar[1].set_xlabel(u'время, с')
            axar[1].set_ylabel(u'задержка, пс')
            plt.gcf().canvas.set_window_title(ifile)
            plt.draw()
        
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
        for k in range(len(li)):
            f.write(str(li[k]))
            f.write('\n')

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

    xlist = np.linspace(0, (len(a) / 2 - 0.5), len(a))

    A = (np.vstack([xlist, np.ones(len(xlist))])).transpose()
    m, c = linalg.lstsq(A, diff, rcond = -1)[0]
    trend = m * xlist + c
    
    print '\nThe slope of trend line is', "%.3f" % (m)

    if dbg == 'true':
        plt.figure()
        plt.plot(xlist, diff, 'o')
        plt.plot(xlist, trend)
        plt.grid()
        plt.xlabel(u'время, с')
        plt.ylabel(u'разностная задержка, пс')
        plt.gcf().canvas.set_window_title(u'Разностная задержка')
        plt.show()


if __name__ == '__main__':
    global qwerty, files
    main()

    if dbg == 'true':
        import matplotlib.pyplot as plt

    if ifile == None:
        print 'You should enter the path to your file!\n'
        usage()
        sys.exit()

    if os.path.exists(ifile):

        if os.path.isfile(ifile):
            print 'Hello. Welcome to PCal. Please wait until your file will be ready...'

            try:
                pcal_read(ifile, ntone, acc_period)
            except:
                import struct
                pcal_reading(ifile, ntone, acc_period)
            
            print '\nNow tell me what you want to do:'
            print 'press 1 if you want to plot signal and see the time delay;'
            print 'press 2 if you want to plot tilt angle and STD graphics;'
            print 'press 3 if you want to see the difference between time delays in 2 different files;'
            print 'press 4 if you want to leave.'
    
            what = raw_input()
            if what == '2':
                pcal_phaseresponse(ifile, ntones, dbg, mode_filtration, write)
            elif what == '1':
                qwerty = '1'
                pcal_delay(ifile, ntones, dbg, qwerty, write, mode_filtration)
            elif what == '3':
                qwerty = '0'
                a = pcal_delay(ifile, ntones, dbg, qwerty, write, mode_filtration)
        
                print '\nPlease enter the 2nd file...'
                new_ifile = raw_input()
        
                try:
                    pcal_read(new_ifile, ntone, acc_period)
                except:
                    import struct
                    pcal_reading(new_ifile, ntone, acc_period)
        
                b = pcal_delay(new_ifile, ntones, dbg, qwerty, write, mode_filtration)

                pcal_diff(a, b)

            elif what == '4':
                sys.exit()

            else:
                print 'Error, please try again.'

        elif os.path.isdir(ifile):
            print 'Hello. Welcome to PCal.'
            print 'Press 1 if you want to calculate time delays and see differences between them.'
            print 'Press 2 if you want to analize STD graphics of all your files.'
            print 'Press 3 if you want to leave.'
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
                        pcal_read(files[i], ntone, acc_period)
                        a = pcal_delay(files[i], ntones, dbg, qwerty, write, mode_filtration)

                        pcal_read(files[i + 1], ntone, acc_period)
                        b = pcal_delay(files[i + 1], ntones, dbg, qwerty, write, mode_filtration)

                        if (files[i])[0 : 17] == (files[i + 1])[0 : 17]:
                            pcal_diff(a, b)

                        i = i + 2
        
                except NameError:
                    ntone = '1 : 512'
            
                    i = 0
                    while i < (len(files) - 1):
                        pcal_read(files[i], ntone, acc_period)
                        poiuy = '2'
                        a = pcal_delay(files[i], ntones, dbg, qwerty, write, mode_filtration)

                        pcal_read(files[i + 1], ntone, acc_period)
                        poiuy = '3'
                        b = pcal_delay(files[i + 1], ntones, dbg, qwerty, write, mode_filtration)

                        if (files[i])[0 : 17] == (files[i + 1])[0 : 17]:
                            pcal_diff(a, b)

                        i = i + 2

            elif mode == '2':
                try:
                    all_files = os.listdir(ifile)

                    files = []

                    for i in range(len(all_files)):
                        if (all_files[i])[0:5] == 'PCAL_' and (all_files[i])[-1] == 'V':
                            files.append(ifile + all_files[i])

                    files = sorted(files)

                    for i in range(len(files)):
                        pcal_read(files[i], ntone, acc_period)
                        pcal_phaseresponse(files[i], ntones, dbg, mode_filtration, write)

                        sys.stdout.write(' files have been processed...' + ('\r%d'%(i + 1) + '/' + str(len(files))))
                        sys.stdout.flush()
                   
                except NameError:
                    ntone = '1 : 512'

                    all_files = os.listdir(ifile)

                    files = []

                    for i in range(len(all_files)):
                        if (all_files[i])[0:5] == 'PCAL_' and (all_files[i])[-1] == 'V':
                            files.append(ifile + all_files[i])

                    files = sorted(files)

                    for i in range(len(files)):
                        pcal_read(files[i], ntone, acc_period)
                        pcal_phaseresponse(files[i], ntones, dbg, mode_filtration, write)

                        sys.stdout.write(' files have been processed...' + ('\r%d'%(i + 1) + '/' + str(len(files))))
                        sys.stdout.flush()

            elif mode == '3':
                sys.exit()

    else:
        print 'Sorry, file or directory was not found. Please try again'