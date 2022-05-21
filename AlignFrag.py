from ost.bindings import tmtools
import ost.io as io
import ost.mol.alg
from ost import *
from ost.seq import *
from ost import bindings
from ost.bindings import tmtools as tm
import math

def dist(i,j):
    d = math.sqrt(pow((i.x-j.x),2)+pow((i.y-j.y),2)+pow((i.z-j.z),2))
    return d


def GetTmalign(ent_target,ent_subunit):
    #result = bindings.WrappedTMAlign(ent_subunit.chains[0], ent_target.chains[0])
    result = tm.TMAlign(ent_target, ent_subunit,"./TMalign")
    seq_t = result.alignment.GetSequence(1)
    seq_sub=result.alignment.GetSequence(0)

    seq_tn= seq_t.Copy()
    seq_subn=seq_sub.Copy()

    #seq_tn.Normalise()
    #seq_subn.Normalise()

    seq_tn.AttachView(ent_target)
    seq_subn.AttachView(ent_subunit)
    
    return result,seq_tn,seq_subn

def GetDistanceGlob(pdb_target,pdb_glob):
    pdb_target = io.LoadPDB(pdb_target).Select("peptide=true")
    pdb_subunit = io.LoadPDB(pdb_glob).Select("peptide=true")
    
    result,seq_tn,seq_subn = GetTmalign(pdb_target,pdb_subunit)
    
    CA_list=[]

    for i in range(0,len(result.alignment)):
        
        if(result.alignment[i] !="-"):
            if(seq_subn.GetResidue(i).IsValid() and seq_tn.GetResidue(i).IsValid()):
                
                res_target = seq_tn.GetResidue(i)
                sub_target = seq_subn.GetResidue(i)
                #print(seq_subn.GetResidue(i).index,seq_tn.GetResidue(i).index)
                CA_target = res_target.FindAtom("CA")
                #CA_subunit = sub_target.FindAtom("CA")
                #print(sub_target.GetNumber(),res_target.GetNumber())
                #CA_target_index = CA_target.index
                #CA_subunit_index = CA_subunit.index
                
                CA_target_pos = CA_target.GetPos()
                #CA_subunit_pos = CA_subunit.GetPos()
                
                #print(CA_target_index,CA_subunit_index)
                #CA_list.append((CA_subunit_index,CA_target_pos,sub_target.index))
                #print(sub_target, res_target)
                #print((sub_target.index,CA_target_pos))
                int_target = int(str(sub_target.GetNumber()))
                #print(int_target,CA_target_pos)
                CA_list.append((int_target,CA_target_pos))

    return ComputeDistanceMap(CA_list)
def GetDistanceFragNoTm(pdb_target,pdb_subunit,min_sub,max_sub):
    pdb_target_view = io.LoadPDB(pdb_target).Select('rnum<=%i and rnum>=%i'%(max_sub,min_sub))
    pdb_subunit_view = io.LoadPDB(pdb_subunit).Select('rnum<=%i and rnum>=%i'%(max_sub,min_sub))
    #pdb_ref = io.LoadPDB(pdb_reference).CreateFullView()
    
    CA_list=[]
    count_CA = min_sub-1
    for i in pdb_target_view.atoms:
        if i.name == "CA":
            count_CA=count_CA+1
            CA_target_pos = i.GetPos()
            
            CA_list.append((count_CA,CA_target_pos))
    
    return ComputeDistanceMapFilt(CA_list)
    

    
def GetDistanceFrag(pdb_target,pdb_subunit,min_sub,max_sub):
    
    pdb_target_view = io.LoadPDB(pdb_target).Select("")
    pdb_subunit_view = io.LoadPDB(pdb_subunit).Select('rnum<=%i and rnum>=%i'%(max_sub,min_sub))
    #pdb_ref = io.LoadPDB(pdb_reference).CreateFullView()

    
    result,seq_tn,seq_subn = GetTmalign(pdb_target_view,pdb_subunit_view)
    
    CA_list=[]

    for i in range(0,len(result.alignment)):
        
        if(result.alignment[i] !="-"):
            if(seq_subn.GetResidue(i).IsValid() and seq_tn.GetResidue(i).IsValid()):
                
                res_target = seq_tn.GetResidue(i)
                sub_target = seq_subn.GetResidue(i)
                #print(seq_subn.GetResidue(i).index,seq_tn.GetResidue(i).index)
                CA_target = res_target.FindAtom("CA")
                #CA_subunit = sub_target.FindAtom("CA")
                #print(sub_target.GetNumber(),res_target.GetNumber())
                #CA_target_index = CA_target.index
                #CA_subunit_index = CA_subunit.index
                
                CA_target_pos = CA_target.GetPos()
                #CA_subunit_pos = CA_subunit.GetPos()
                
                #print(CA_target_index,CA_subunit_index)
                #CA_list.append((CA_subunit_index,CA_target_pos,sub_target.index))
                #print(sub_target, res_target)
                #print((sub_target.index,CA_target_pos))
                int_target = int(str(sub_target.GetNumber()))
                #print(int_target,CA_target_pos)
                CA_list.append((int_target,CA_target_pos))

    return ComputeDistanceMapFilt(CA_list)
def GetDistanceGlob(pdb_target,pdb_subunit):
    
    pdb_target_view = io.LoadPDB(pdb_target).Select("chain=A")
    pdb_subunit_view = io.LoadPDB(pdb_subunit).Select('chain=A')
    #pdb_ref = io.LoadPDB(pdb_reference).CreateFullView()

    
    result,seq_tn,seq_subn = GetTmalign(pdb_target_view,pdb_subunit_view)
    
    CA_list=[]

    for i in range(0,len(result.alignment)):
        
        if(result.alignment[i] !="-"):
            if(seq_subn.GetResidue(i).IsValid() and seq_tn.GetResidue(i).IsValid()):
                
                res_target = seq_tn.GetResidue(i)
                sub_target = seq_subn.GetResidue(i)
                #print(seq_subn.GetResidue(i).index,seq_tn.GetResidue(i).index)
                CA_target = res_target.FindAtom("CA")
                #CA_subunit = sub_target.FindAtom("CA")
                #print(sub_target.GetNumber(),res_target.GetNumber())
                #CA_target_index = CA_target.index
                #CA_subunit_index = CA_subunit.index
                
                CA_target_pos = CA_target.GetPos()
                #CA_subunit_pos = CA_subunit.GetPos()
                
                #print(CA_target_index,CA_subunit_index)
                #CA_list.append((CA_subunit_index,CA_target_pos,sub_target.index))
                #print(sub_target, res_target)
                #print((sub_target.index,CA_target_pos))
                int_target = int(str(sub_target.GetNumber()))
                #print(int_target,CA_target_pos)
                CA_list.append((int_target,CA_target_pos))

    return ComputeDistanceMap(CA_list)

def GetListCaFrag(pdb_target,pdb_subunit,min_sub,max_sub):
    
    pdb_target_view = io.LoadPDB(pdb_target).Select("chain=A")
    pdb_subunit_view = io.LoadPDB(pdb_subunit).Select('rnum<=%i and rnum>=%i'%(max_sub,min_sub))
    #pdb_ref = io.LoadPDB(pdb_reference).CreateFullView()

    
    result,seq_tn,seq_subn = GetTmalign(pdb_target_view,pdb_subunit_view)
    
    CA_list={}

    for i in range(0,len(result.alignment)):
        
        if(result.alignment[i] !="-"):
            if(seq_subn.GetResidue(i).IsValid() and seq_tn.GetResidue(i).IsValid()):
                
                res_target = seq_tn.GetResidue(i)
                sub_target = seq_subn.GetResidue(i)
                #print(seq_subn.GetResidue(i).index,seq_tn.GetResidue(i).index)
                CA_target = res_target.FindAtom("CA")
                #CA_subunit = sub_target.FindAtom("CA")
                #print(sub_target.GetNumber(),res_target.GetNumber())
                #CA_target_index = CA_target.index
                #CA_subunit_index = CA_subunit.index
                
                CA_target_pos = CA_target.GetPos()
                #CA_subunit_pos = CA_subunit.GetPos()
                
                #print(CA_target_index,CA_subunit_index)
                #CA_list.append((CA_subunit_index,CA_target_pos,sub_target.index))
                #print(sub_target, res_target)
                #print((sub_target.index,CA_target_pos))
                int_target = int(str(sub_target.GetNumber()))
                #print(int_target,CA_target_pos)
                CA_list[int_target]=CA_target_pos

    return CA_list
        
def ComputeDistanceMap(CA_list):  
    lCA_map_temp ={}
    
    for sub_index,pos in CA_list:
        for i,j in CA_list:
            if(sub_index != i):
                #print((sub_index,i),round(dist(pos,j)/10,3))  
                lCA_map_temp.update({(sub_index,i):round(dist(pos,j)/10,3)})
    return lCA_map_temp

def ComputeDistanceMapFilt(CA_list):  
    lCA_map_temp ={}
    
    for sub_index,pos in CA_list:
        for i,j in CA_list:
            if(sub_index != i):
                
                #if round(dist(pos,j)/10,3) < 1.5:
                lCA_map_temp.update({(sub_index,i):round(dist(pos,j)/10,3)})
                    #print((sub_index,i),round(dist(pos,j)/10,3))
    return lCA_map_temp

    

