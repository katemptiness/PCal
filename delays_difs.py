def delays_difs(a):
    import os
    import PCal
    reload(PCal)
    
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
    delays_zv = []

    files_bv = []
    files_zv = []

    k = 0
    while k < len(files):
        if (files[k])[18:20] == 'BV':
            delays_bv.append(PCal.pcal_delay(['-f files[k]']))
            files_bv.append(files[k])
        elif (files[k])[18:20] == 'ZV':
            delays_zv.append(PCal.pcal_delay(['-f files[k]']))
            files_zv.append(files[k])
        k = k + 1

    delays_zv.append(delays_zv[0])
    del delays_zv[0]

    files_zv.append(files_zv[0])
    del files_zv[0]

    delays_diff = []
    i = 0
    while i < len(files_bv):
        if (files_bv[i])[:18] == (files_zv[i])[:18]:
            delays_diff.append(abs(delays_bv[i] - delays_zv[i]))
        i = i + 1
    
    return delays_diff

def help_me():
    return b