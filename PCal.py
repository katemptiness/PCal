#!/usr/bin/env python

import numpy as np, numpy.fft as fft, numpy.linalg as linalg
import cmath, math
import getopt, sys
import matplotlib.pyplot as plt
import os
from itertools import count, izip

#ch = 0

def file_read(ifile):
    ifile = open(ifile)

    i = 0
    while i < 5:
        ifile.readline()
        i = i + 1

    counter = int((str(ifile.readline()).split())[5])

    return counter


def main():
    global itype, dbg, ntones, ifile, n

    ntones = '1 : 512'
    itype = 'phase'
    dbg = 'false'

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hf:n:t:d:', ['ifile=', 'ntones=', 'itype=', 'dbg='])
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
        elif opt in ('-n', '--ntones'):
            if arg[0] == ' ':
                ntones = arg[1:]
            else:
                ntones = arg
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

    
def usage():
    print 'Hello, this is the USAGE function.'
    print
    print 'You can use this forms to work:'
    print '-f is for path to the file;'
    print '-n is for tone numbers;'
    print '-t is for data type (amplitudes or phases);'
    print '-d is for graphics display mode (true or false).'
    print
    print 'For example:'
    print 'python PCal.py -f W:/Files/My_File -n "1 : 512" -t phase - d true'
    print
    print '-n, -t & -d parameters are "1 : 512", "phase" & "false" as default, so you can only use -f.'


def unwraping(lista):
    if np.std(lista) > 49:

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
                    lista[i + 1] = lista[i + 1] + 2 * 180
                elif (lista[i + 1] - lista[i]) > 180:
                    lista[i] = lista[i] + 2 * 180
                i = i + 1
    
        elif c == '-':
            i = 0
            while i < (len(lista) - 1):
                if (lista[i + 1] - lista[i]) < -300:
                    lista[i] = lista[i] - 2 * 180
                elif (lista[i + 1] - lista[i]) > 300:
                    lista[i + 1] = lista[i + 1] - 2 * 180
                i = i + 1

        if np.mean(lista) > 0:
            if lista[0] < 0 and lista[1] < 0:
                lista[0] = lista[0] + 2 * 180
                lista[1] = lista[1] + 2 * 180
        elif np.mean(lista) < 0:
            if lista[0] > 0 and lista[1] > 0:
                lista[0] = lista[0] - 2 * 180
                lista[1] = lista[1] - 2 * 180
    
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
    
    try:
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
        
        del Re, Im, acos, asin, I
        return Im1, Re1
    
    except:
        del Re, Im, acos, asin, I
        
        return Im1, Re1


def pcal_read(ifile, ntones, itype, dbg):
    global table, table2, acc_periods, counter, ph, ntones_full
    
    #if ifile[:5] == 'files':
        #import delays_difs
	#reload(delays_difs)
        #p = delays_difs.help_me()
        
        #a = os.listdir(p)
        
        #files = []
        #i = 0
        #while i < len(a):
            #if (a[i])[:5] == 'PCAL_':
                #files.append(a[i])
            #i = i + 1
        #if type(ifile[6]) is str:
            #k = ch
            #ch = ch + 1
            #ifile = files[int(k)]
        #else:
            #ifile = files[int(ifile[6])]
    
        #ifile = p + '/' + ifile
    
    ifile = open(ifile)
    acc_periods = len((ifile).readlines()) - 5
    ifile.seek(0)
    
    k = 1
    while k <= 5:
        ifile.readline()
        k = k + 1
    
    ntones = (str(ntones).replace(',', '')).split()
    
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
    ntones_full = ntones
    counter = len(ntones_full)
    
    k = 0
    while k < counter:
        ntones_full[k] = int(ntones_full[k]) - 1
        k = k + 1

    table = (np.empty((counter, 0))).tolist()
    ph = []

    table2 = (np.empty((counter, 0))).tolist()
    
    j = 0
    while j < acc_periods:
        smh = ((ifile.readline()).split())[6:]
        i = 0
        while i < counter:
            if itype == 'phase':
                table[i].append(cmath.phase(complex(float(smh[(len(smh) - 1) - 1 - int(ntones_full[i]) * 4]), float(smh[(len(smh) - 1) - int(ntones_full[i]) * 4]))) * (180 / np.pi))
            elif itype == 'amplitude':
                table[i].append(math.hypot(float(smh[(len(smh) - 1) - 1 - int(ntones_full[i] * 4)]), float(smh[(len(smh) - 1) - int(ntones_full[i]) * 4])) * 1000)
            elif itype == 'phase-amplitude':
                table[i].append(cmath.phase(complex(float(smh[(len(smh) - 1) - 1 - int(ntones_full[i]) * 4]), float(smh[(len(smh) - 1) - int(ntones_full[i]) * 4]))) * (180 / np.pi))
                table2[i].append(math.hypot(float(smh[(len(smh) - 1) - 1 - int(ntones_full[i] * 4)]), float(smh[(len(smh) - 1) - int(ntones_full[i]) * 4])) * 1000)
            ph.append(complex(float(smh[2 + int(ntones_full[i]) * 4]), float(smh[3 + int(ntones_full[i]) * 4])))
            i = i + 1
        j = j + 1
    
    ifile.close()

    if itype == 'phase-amplitude':
        return 
    else: 
        return table
    

def pcal_delay(ifile, ntones, itype, dbg):
    if itype == 'phase-amplitude':
        plt.plot(table, table2, 'o')
        plt.grid()
        plt.xlabel('phase')
        plt.ylabel('amplitude')
        plt.show()
    
    else:
        li = []

        f, axar = plt.subplots(2)

        #time = np.linspace(0, 1, 512)

        j = 0
        while j < acc_periods:
            ph1 = abs(fft.ifft(ph[(j * counter) : (j * counter + counter)]))
            #axar[0].plot(time, ph1)
            axar[0].plot(ph1)

            number = max(izip(ph1, count()))[1]

            j0 = number

            if j == 0:
                #print 'The time delay is probably', ((j0 * 1e-6) / 512), 'microseconds'
                print 'The time delay is probably', j0
            
            tau_min = j0 - 1
            tau_max = j0 + 1

            delta_tau = 0.1
            while delta_tau >= 1e-3:
                tau = tau_min
            
                tau_list = []
                cj = []
            
                while tau <= tau_max:
                    im, re = Fraq_FFT(file_read(ifile), (np.asarray(ph).real)[(j * counter) : (j * counter + counter)], (np.asarray(ph).imag)[(j * counter) : (j * counter + counter)], tau, 0)
                
                    cj.append(complex(re, im))
                    tau_list.append(tau)

                    tau = tau + delta_tau
                
                cj = abs(np.asarray(cj))

                number = max(izip(cj, count()))[1]
                tau = tau_list[number]

                tau_min = tau - delta_tau * 2
                tau_max = tau + delta_tau * 2
                delta_tau = delta_tau / 10
        
            #tau = tau / 512
            #tau = float("%.6f" % (tau))

            print tau

            li.append(tau)

            j = j + 1

        #print new_li[-1]

        tau = "%.6f" % (np.mean(li))

        #print 'And the clarified time delay is', tau, 'microseconds'
        print 'And the clarified time delay is', tau

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


def pcal_trend(ifile, ntones, itype, dbg):
    global trends, std, new_table

    f, axar = plt.subplots(4)

    time = np.linspace(0, 0.5 * acc_periods, acc_periods)
    
    trends = []
    
    i = 0
    while i < counter:
        axar[0].plot(time, unwraping(table[i - 1]), 'o')
        A = (np.vstack([time, np.ones(len(time))])).transpose()
        m, c = linalg.lstsq(A, unwraping(table[i - 1]))[0]
        trend = m * time + c
        trends.append(trend)
        if dbg == 'true':
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
            
    if dbg == 'true':
        axar[0].grid()
        axar[0].set_xlabel('time')
        axar[0].set_ylabel(itype)
        
        j = 0
        while j < counter:
            BC = (trends[j])[acc_periods - 1] - (trends[j])[0]
            AB = np.sqrt(BC * BC + AC * AC)
            alpha = np.arcsin(BC / AB) * (180 / np.pi)
            if dbg == 'true':
                axar[1].plot((ntones_full[j] + 1), alpha, 'o')
            j = j + 1

        axar[1].grid()
        axar[1].set_xlabel('tone numbers')
        axar[1].set_ylabel('tilt angle')

        j = 0
        while j < counter:
            axar[2].plot((ntones_full[j] + 1), std[j], 'o')
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

    i = 0
    while i < len(good_ntones):
        li.append(np.mean(good_table[i]))
        i = i + 1

    plt.plot(good_ntones, unwraping2(li))
    plt.plot(good_ntones, unwraping2(li), 'o')
    
    trends = []
    
    A = (np.vstack([good_ntones, np.ones(len(good_ntones))])).transpose()
    m, c = linalg.lstsq(A, li)[0]
    trend = m * good_ntones + c
    
    if dbg == 'true':
        plt.plot(good_ntones, trend)
        plt.grid()
        plt.xlabel('frequency')
        if itype == 'phase-amplitude':
            plt.ylabel('phase')
        else:
            plt.ylabel(itype)
        plt.show()
    
    a = abs(max(trend) - min(trend))
    b = (counter - 1) * (10 ** 6)

    delay = (a * (np.pi / 180)) / b
        
    if __name__ == '__main__':
        print 'The time delay is ', delay
    else:
        return delay


if __name__ == '__main__':
    main()
    
    pcal_read(ifile, ntones, itype, dbg)

    print 'Hello. Welcome to PCal interface.'
    print 'Now tell me what you want to do:'
    print 'press 1 if uou want to plot signal and see the time delay;'
    print 'press 2 if you want to plot tilt angle and STD graphics and see phase-frequency response.'
    
    what = raw_input()

    if what == '2':
        pcal_phaseresponse(ifile, ntones, itype, dbg)
    elif what == '1':
        pcal_delay(ifile, ntones, itype, dbg)
    else:
        print 'Error, please try again.'