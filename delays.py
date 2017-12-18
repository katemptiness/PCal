import os, PCal
reload(PCal)

def delays_bv(a):
    global b
    b = a

    a = os.listdir(a)

    files = []

    i = 0
    while i < len(a):
        if (a[i])[:5] == 'PCAL_':
            files.append(a[i])
        i = i + 1

    delays_bv = []
    files_bv = []

    k = 0
    while k < len(files):
        if (files[k])[18:20] == 'BV':
            delays_bv.append(PCal.pcal_delay(['-f files[k]', '-n 128 : 384']))
            files_bv.append(files[k])
        elif (files[k])[18:20] == 'ZV':
            delays_zv.append(PCal.pcal_delay(['-f files[k]', '-n 128 : 384']))
            files_zv.append(files[k])
        k = k + 1

    delays_bv.sort()
    
    return delays_bv


def delays_zv(a):
    global b
    b = a
    
    a = os.listdir(a)

    files = []

    i = 0
    while i < len(a):
        if (a[i])[:5] == 'PCAL_':
            files.append(a[i])
        i = i + 1

    delays_zv = []
    files_zv = []

    k = 0
    while k < len(files):
        if (files[k])[18:20] == 'BV':
            delays_bv.append(PCal.pcal_delay(['-f files[k]', '-n 128 : 384']))
            files_bv.append(files[k])
        elif (files[k])[18:20] == 'ZV':
            delays_zv.append(PCal.pcal_delay(['-f files[k]', '-n 128 : 384']))
            files_zv.append(files[k])
        k = k + 1

    delays_zv.sort()
    
    return delays_zv


def help_me():
    return b