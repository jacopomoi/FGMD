import os


count = 0
def selglob(PDB,model_file):

    print("processing ",model_file)
    os.system("./TMalign -dir1 "+PDB+" "+PDB+"/list -suffix .pdb "+model_file+" -outfmt 2 -fast -TMcut -1 > ./GLOB_list/list_glob")
    

