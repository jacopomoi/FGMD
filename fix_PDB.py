import os

def fix(frag_folder):
	
    print("Fixing pdb :",frag_folder)
    for i in os.listdir(frag_folder):
        #os.system("pdb_chain -A "+frag_folder+i +">"+frag_folder+i+".ok")
        os.system("pdb_chain -A "+frag_folder+i +">"+frag_folder+i+".ok")
        os.system("rm "+frag_folder+i)
        os.system("pdb_reres "+frag_folder+i+".ok" +">"+frag_folder+i+".red")
        os.system("rm "+frag_folder+i+".ok")
        os.system("pdb_keepcoord "+frag_folder+i+".red" +">"+frag_folder+i)
        os.system("rm "+frag_folder+i+".red")
        os.system("pdb_delhetatm "+frag_folder+i+">"+frag_folder+i+".red")
        os.system("rm "+frag_folder+i)
        os.system("mv "+frag_folder+i+".red"+" "+frag_folder+i)

if __name__ == "__main__":
    fix("./sel_FRAG/")
