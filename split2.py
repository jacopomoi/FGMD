import os
cnt = 0
for i in os.listdir("Traj_FGMD"):
	print(cnt)
	os.system("cd /home/jacopo/Desktop/temp/FGMD-noTm")
	os.system("mkdir /home/jacopo/Desktop/temp/FGMD-noTm/steps/"+i)
	os.system("cd steps/"+i )
	os.system("cd steps/"+i+" |pdb_splitmodel /home/jacopo/Desktop/temp/FGMD-noTm/Traj_FGMD/"+i)
	#os.system("mkdir steps/"+i)
	#os.system("mv tt steps/"+i)
	
