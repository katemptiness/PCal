import numpy as np, numpy.fft as fft, numpy.linalg as linalg
import cmath, math
import getopt, sys
import matplotlib.pyplot as plt

def usage():
    print 'Hello, this is the USAGE function.'
    print 'Use this form to make the program work correctly: -f <the path to the ifile> -n <tone numbers> -t <phase / amplitude>.'
    print 'For example: pcal_read("-f W:/Files/My_File", "-n 1 : 20, 40, 25, 300 : 408", "-t phase")'

def pcal_read(argv):
    global table, acc_periods, counter, ph, ntones, ifile, type

    try:
        opts, args = getopt.getopt(argv, 'hf:n:t:', ['ifile=', 'ntones=', 'type='])
    except getopt.GetoptError:
        print 'Looks like something went wrong. \n'
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
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
        elif opt in ('-t', '--type'):
            if arg[0] == ' ':
                type = arg[1:]
            else:
                type = arg
    
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
    counter = len(ntones)
    
    k = 0
    while k < counter:
        ntones[k] = int(ntones[k]) - 1
        k = k + 1
    
    table = (np.empty((counter, 0))).tolist()
    ph = []
    
    j = 0
    while j < acc_periods:
        smh = ifile.readline()
        smh = (smh.split())[6:]
        i = 0
        while i < counter:
            if type == 'phase':
                table[i].append(cmath.phase(complex(float(smh[(len(smh) - 1) - 1 - int(ntones[i]) * 4]), float(smh[(len(smh) - 1) - int(ntones[i]) * 4]))) * (180 / np.pi))
            elif type == 'amplitude':
                table[i].append(math.hypot(float(smh[(len(smh) - 1) - 1 - int(ntones[i] * 4)]), float(smh[(len(smh) - 1) - int(ntones[i]) * 4])) * 1000)
            ph.append(complex(float(smh[2 + int(ntones[i]) * 4]), float(smh[3 + int(ntones[i]) * 4])))
            i = i + 1
        j = j + 1
        
    q = raw_input('Print the table (y/n)? ')
    if q == 'y':
        i = 0
        while i < counter:
            print table[i]
            i = i + 1

    ifile.close()


def pcal_plot(argv):
    pcal_read(argv)
    
    time = np.linspace(0, 0.5 * acc_periods, acc_periods)
    
    i = counter
    while i > 0:
        plt.plot(time, np.unwrap(table[i - 1]))
        i = i - 1
    
    if type == 'phase':
        plt.axis([0, acc_periods * 0.5, -200, 200])
    elif type == 'amplitude':
        plt.axis([0, acc_periods * 0.5, 0, 0.006])

    plt.grid()
    plt.xlabel('time')
    plt.ylabel(type)
    plt.show()
    
    condition = raw_input('Plot the signal (y/n)? ')
    if condition == 'y':
        time = np.linspace(0, 1e-6, counter)
        j = 0
        while j < acc_periods:
            ph1 = abs(fft.ifft(ph[(j * (counter - 1)) : (j * (counter - 1) + counter)]))
            j = j + 1
            plt.plot(time, ph1)
        plt.grid()
        plt.xlabel('time')
        plt.ylabel('amplitude')
        plt.show()


def pcal_trend(argv):
    global trends

    pcal_read(argv)
    
    time = np.linspace(0, 0.5 * acc_periods, acc_periods)
    
    trends = []
    
    i = 0
    while i < counter:
        plt.plot(time, np.unwrap(table[i - 1]), 'o')
        A = (np.vstack([time, np.ones(len(time))])).transpose()
        m, c = linalg.lstsq(A, np.unwrap(table[i - 1]))[0]
        trend = m * time + c
        trends.append(trend)
        plt.plot(time, trend)
        i = i + 1
    
    plt.grid()
    plt.xlabel('time')
    plt.ylabel(type)
    plt.show()
    
    AC = len(time)
    
    alphas = []
    
    j = 0
    while j < counter:
        BC = (trends[j])[acc_periods - 1] - (trends[j])[0]
        AB = np.sqrt(BC * BC + AC * AC)
        alpha = np.arcsin(BC / AB) * (180 / np.pi)
        plt.plot((ntones[j] + 1), alpha, 'o')
        j = j + 1
    
    plt.grid()
    plt.xlabel('tone numbers')
    plt.ylabel('tilt angle')
    plt.show()


def pcal_retrend(argv):
    global std
    
    pcal_trend(argv)
    
    re_trends = []
    re_table = []
    
    std = []
    
    i = 0
    while i < counter:
        j = 0
        while j < acc_periods:
            re = (trends[i])[1] - (trends[i])[0]
            re_trends.append((trends[i])[j] - j * re)
            re_table.append((table[i])[j] - j * re)
            j = j + 1
        std.append(np.std(np.unwrap(re_table)))
        re_table = []
        i = i + 1
            
    f, axar = plt.subplots(2)
    
    j = 0
    while j < counter:
        axar[0].plot((ntones[j] + 1), std[j], 'o')
        j = j + 1
    
    axar[0].grid()
    axar[0].set_xlabel('tone numbers')
    axar[0].set_ylabel('standard deviation')
    
    axar[1].hist(std)
    axar[1].set_xlabel('standard deviation')
    axar[1].set_ylabel('tones')
    
    plt.show()


def pcal_filt(argv):
    global std_good

    pcal_retrend(argv)

    std_good = []
    i = 0
    while i < len(std):
        if std[i] < 2:
            std_good.append(std[i])
        i = i + 1


def pcal_pfr(argv):
    pcal_read(argv)
    
    li = []
    
    i = 0
    while i < counter:
        j = 0
        sum = 0
        while j < acc_periods:
            sum = sum + (table[i])[j]
            j = j + 1
        li.append(sum / acc_periods)
        i = i + 1
    
    plt.plot(ntones, np.unwrap(li))
    
    trends = []
    
    A = (np.vstack([ntones, np.ones(len(ntones))])).transpose()
    m, c = linalg.lstsq(A, np.unwrap(li))[0]
    trend = m * ntones + c
    
    plt.plot(ntones, trend)
    
    plt.grid()
    plt.xlabel('frequency')
    plt.ylabel(type)
    plt.show()
    
    if type == 'phase':
        #AB = abs(max(trend) - min(trend))
        #BC = (counter - 1) * 10 ** 6
        #AC = math.hypot(AB, BC)
        #alpha = np.arcsin(AB / AC) * (180 / np.pi)

        #delay = np.tan(alpha)
        #print 'For this range the delay is', delay

        a = abs(max(trend) - min(trend))
        b = (counter - 1) * (10 ** 6)

        delay = (a * (np.pi / 180)) / b
        print delay