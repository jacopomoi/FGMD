import selfrag as sel
import split_SSE2 as split
import takefrag as take
import selglob as selg
import takeglob as tkg
import FGMDveramb_grad as sim
#import takeglob as takeg
PDBDB ="/rhome/847329/PDBhres/"
PDB_FOLD = "./model/"
model = "type1.pdb"
print("Split SSE")
split.splitSSE(PDB_FOLD+model)
"""
print("TMalign PDB -> Fragment")
sel.selfrag(PDBDB)
PDBDB_take ="/rhome/847329/PDBhres/"
print("Takin PDB -> Fragment")
take.Takefrag("./SSE_list/",PDBDB_take)
print("TMalign PDB -> Global template")
selg.selglob(PDBDB,PDB_FOLD+model)
print("Takin PDB -> Global template")
tkg.Takeglob("./GLOB_list/",PDBDB_take,"./sel_GLOB/")
#sim.run()"
"""
