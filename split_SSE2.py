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

	SS_new = []
	for i,name in SS:
		if name == "T" or name == "C":
			SS_new.append((int(i),"-"))
		else:
			SS_new.append((int(i),name))


# In[32]:


	f = SS[0][1]
	min_f = SS[0][0]
	max_f = 0

	SS_list = []

	for pos,sec in SS_new:
		if sec == f:
			max_f = pos
		else:
			SS_list.append([min_f,max_f,str(f)])
			f = sec
			min_f = pos
	SS_list = SS_list[1:]


# In[39]:


	SS_only = []
	for i in SS_list:
		
		if i[2] != "T" and i[2] != "-":
			SS_only.append(i)


# In[40]:


	SS_cont = []
	for k,l in zip(SS_list,SS_list[6:]):
		SS_cont.append([k[0],l[1]])

	


	print(SS_cont)

	try:
		os.mkdir("SSE")
	except:
		pass

	count = 1
	
	for i in SS_cont:
		os.system("pdb_selres -"+str(i[0])+":"+str(i[1])+" "+model_file+" "+">"+"./SSE/SSE_"+str(count)+"_"+str(i[0])+":"+str(i[1])+".pdb")
		count = count+1









