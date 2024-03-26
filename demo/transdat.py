def read_trans_dat(filename):
    openfile = open(filename)
    
    lines = openfile.readlines()
    openfile.close()

    lines.pop(0)    
    wavs = []
    trans = []
    for i in lines:
        i = i.split()
        
        wavs.append(float(i[0]))
        trans.append(float(i[1]))

    return wavs, trans