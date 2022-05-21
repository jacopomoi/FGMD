import os

DIR = "./SSE/"
count = 0
def selfrag(PDB):
	DIR = "./SSE/"
	count=0
	for i in os.listdir(DIR):
		try:
			count=int(i.split("_")[1].split(".")[0])
			max = i.split("_")[2].split(".")[0].split(":")[1]
			min = i.split("_")[2].split(".")[0].split(":")[0]
			if ("SSE_"+str(min)+":"+str(max)) not in os.listdir("./SSE_list"):
				print(count)
				print("processing",i)
				os.system("./TMalign -dir1 "+PDB+" "+PDB+"/list -suffix .pdb "+DIR+i+" -outfmt 2 > ./SSE_list/SSE_"+min+":"+max)
			else:
				print(i+"just parsed in SSE_list folder")
		except Exception as e:
			print(e)
			pass
		
     
