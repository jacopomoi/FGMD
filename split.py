import os

for i in os.listdir("Traj_FGMD"):
	with open("Traj_FGMD/"+i,"r") as file:
		line = file.readline()
		cnt = 0
		while line:
			line=file.readline()
			cnt=cnt+1
			if "MODEL" in line.split("   "):
				num=line.split("\n")[0].split("      ")[1]
				print(num)
				if int(num)==70:
					print("ok")
					buf=[]
					while line:
					
						line = file.readline()
						buf.append(line )
						print(line)
					with open("laststep_FGMD/"+i,"w") as file2:
						for k in buf:
							file2.write(k+"\n")
				

