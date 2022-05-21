#!/usr/bin/env python
# coding: utf-8

# In[100]:


import csb.bio.io as io
import os 
from pdbtools import pdb_chainbows as chbw
def stride(model_file):
    
	os.system("rm pred -f")
	#os.system("pdb_chainbows "+model_file+">"+model_file+".tmp")
	#os.system("pdb_chain -A "+model_file+">"+model_file+".tmp")
	#os.system("rm "+model_file)
	#os.system("mv "+model_file+".tmp "+model_file)
	os.system("./stride "+model_file+">"+"pred")
	return

def splitSSE(model_file):
	
	stride(model_file)
	pr = io.StrideParser()
	ret = pr.parse("pred")


	main_chain = ret.get("A")
	SS = []
	for i in main_chain:
		res = main_chain.get(i).secondary_structure
		
		SS.append((i,res))



	f = SS[0][1]
	min_f = SS[0][0]
	max_f = 0

	SS_list = []

	for pos,sec in SS:
		if sec == f:
			max_f = pos
		else:
			SS_list.append([min_f,max_f,str(f)])
			f = sec
			min_f = pos
        
        



	count = 1
	same = True
	type_ss = SS[0][1]

	SS_sep = []
	SS_sep.append((SS[0],count))

	min_ss = 1
	max_ss = 1

	SS_def = []
	prev = None
	for i in SS[1:]:
		
		if i[1]==type_ss:
			prev = type_ss
			max_ss=int(i[0])
			SS_sep.append((i,count,min_ss,max_ss))
			
			
		elif i!=type_ss:
			if i[1]=="Coil" and type_ss=="Turn":
				SS_sep.append([(i.index,"Turn"),min_ss,max_ss])
			else:
				type_ss = i[1]
				SS_def.append([min_ss,max_ss,prev,])

				min_ss = max_ss
				count=count+1
				SS_sep.append((i,count,min_ss,max_ss))

	SS_def.append([min_ss,max_ss,prev])
        
	import math

	split_seq = 4
	T_count = 0
	SS_filt_coil = []
	SS_sep = []
	for i in range(0,len(SS_def)):
		cur = SS_def[i] 
		if cur[2] == "T":
			T_count=T_count+1
			SS_def[i].append(T_count % 2)

	cnt = 1
	min = SS_def[0][0]
	max = SS_def[0][1]

	for k in SS_def:
		
		if cnt==split_seq:
			SS_sep.append([min+1,max])
			min = k[0]
			
			cnt=1
		else:
			max = k[1]
			cnt = cnt +1
			
	SS_sep.append([SS_sep[len(SS_sep)-1][1],SS_def[len(SS_def)-1][1]])


	print(SS_sep)
	SS_filt_coil = SS_sep


	SS3_cont = []

	for i,j,k in zip(SS_filt_coil,SS_filt_coil[1:],SS_filt_coil[2:]):
		SS3_cont.append([int(i[0]),int(k[1])])


	print(SS3_cont)

	try:
		os.mkdir("SSE")
	except:
		pass

	count = 1
	
	for i in SS3_cont:
		os.system("pdb_selres -"+str(i[0])+":"+str(i[1])+" "+model_file+" "+">"+"./SSE/SSE_"+str(count)+"_"+str(i[0])+":"+str(i[1])+".pdb")
		count = count+1









