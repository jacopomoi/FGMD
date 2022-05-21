import os
import openmm
from pathlib import Path
import shutil
import fix_PDB
frglist = "./SSE_list/"

pdb_fold = "/media/jacopo/SAMSUNG FEDE/PDB2/"

def sort(l):
    return l[1]

def filt(fasta,index_first_res):
    src = list(fasta)
    src_filt = []

    
    count = index_first_res
    for i in src:
        if i != "-":
            src_filt.append((count,i))
            count = count +1
        else:
            src_filt.append(i)
    return src_filt

def parse(src,trg,match):
    trg_filt=filt(trg,1)
    src_filt=filt(src,1)
    
    ret = []
    for j in zip(src_filt,trg_filt,match):
        if j[0]!="-" and j[1]!="-" and j[2]==":":
            ret.append((j[0],j[1]))
    
    return ret


def Max20(text):
    high = []
    highh= []
    
    for i in text:
        try:
            sp = i.split()
            if sp[0] != "Warning!" and sp[0]!="Total":
                high.append((sp[0],sp[3]))
                
        except:
            pass
        
    high.sort(key = sort,reverse = True)
    return high[0:20]
        


#routine per estrarre frammenti dai pdb
max = []
def Takefrag(frglist,pdb_fold):
    for i in os.listdir(frglist):
        with open(frglist + i,"r") as s1:
            t = s1.readlines()[1:]
            
            max = Max20(t)
            if max!=[]:
                for i in range(0,10):
                    maxx = max[i]
                    print(maxx,i)
                    try:
                        shutil.copyfile(pdb_fold+maxx[0],Path("./sel_FRAG") / maxx[0])
                    except Exception as e:
                        print(e)
    fix_PDB.fix("./sel_FRAG/")
                 
        


