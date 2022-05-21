from openmm.app import *
from openmm import *
import openmm.unit as un
from openmm.openmm import XmlSerializer as xml
from sys import stdout
from openmm.openmm import CMMotionRemover, CustomBondForce, CustomNonbondedForce, State, VerletIntegrator,LangevinIntegrator,LangevinMiddleIntegrator,NoseHooverIntegrator
from openmm.openmm import LocalEnergyMinimizer
import openmmtools as aux
import itertools as iter
from openmm.openmm import AndersenThermostat
import numpy as np
import math
from openmmtools import integrators
from ost.bindings import tmtools as tm
from openmmtools.integrators import GradientDescentMinimizationIntegrator, MetropolisMonteCarloIntegrator
from openmmtools.integrators import GeodesicBAOABIntegrator
from ost import io as io
from ost import bindings as bind
import AlignFrag
import takefrag
import takeglob
import time

#INIT_PDB = "./minim_start.pdb"
INIT_PDB="./model/globf.pdb"
CA_PDB = "./model/out_caok.pdb"
OUT_NAME="MSMsimtemp.pdb"
TEMPLATE = "./model/native.pdb"
unm = []
FAIL=[]
T_MAX = 273
TEMPERATURE = 0
T_HIGH=700
TEMP_STEP=25

platform = Platform.getPlatformByName('OpenCL')
properties={}
#properties["DeviceName"] = "gfx804"
#properties['Precision']= 'mixed'

def dist(i,j):
    d = math.sqrt(pow((i.x-j.x),2)+pow((i.y-j.y),2)+pow((i.z-j.z),2))
    return d

def drmsd(list_dr):
    n = len(list_dr)
    sum = 0
    
    for i in list_dr:
        if i[2] != None and i[3]!=None:
            r = i[2]-i[3]
            
            D = pow(r,2)
            sum = sum + D
    
    tot = round(math.sqrt(sum/n),4)
    
    return tot

OUT ="fgmd-6a90-full-final"
for Target in os.listdir("./model/6a90comp/"):
#-------------------------
	print(Target)
	if Target not in os.listdir(OUT):
		try:
			INIT_PDB= "./model/6a90comp/"+ Target
			print('Loading...',Target)
			#template = PDBFile(TEMPLATE)
			init = PDBFile(INIT_PDB)
			print("Done")
			model = Modeller(init.topology, init.positions)
			  
			
			PDBFile.writeFile(model.topology,model.positions,open(OUT+"/"+Target,"w"))
			print(init.getTopology())
			#print(template.getTopology())

			print("loading Forcefield")
			#forcefield = ForceField('ff99FG.xml','tip3p.xml')
			forcefield = ForceField("./forcefield/ff99FG.xml",'./implicit/obc2.xml')
			#forcefield = ForceField("charmm36.xml")
			#forcefield = ForceField("amoeba2018.xml")
			#forcefield=ForceField("amber/amber03.xml","amber/tip3p.xml")
			#forcefield = ForceField('./forcefield/ff99FG.xml')
			#forcefield = ForceField("./forcefield/ff99FG.xml",'tip3p.xml')
			#forcefield = ForceField("./forcefield/ff99SB.xml",'tip3p.xml')
			#print("Done")
			model.addExtraParticles(forcefield)
			model.addHydrogens(forcefield)
			# In[4]:

			
			#topologie
			top_init = init.getTopology()
			#top_template = template.getTopology()

			lCA_init = []
			lCA_map_init = []
			map_frag = {}
			lCA_temp = []
			lCA_map_template = []
			lCA_map_temp={}
			frglist = "./SSE_list/"
			glob_fold = "./sel_GLOB/"
			glblist = "./GLOB_list/"

			print("Computing Distance restraint:Fragment")

			frg_fold = os.listdir(frglist)
			max = []
			proc_fail = 0
			
			for i in frg_fold:
				ss=i.split("_")[1].split(":")
				min_ss = int(ss[0])
				max_ss = int(ss[1])
				with open(frglist + i,"r") as s1:
					t = s1.readlines()[1:]
					try:
						max = takefrag.Max20(t)[0:1]
						print("Sel Frag: ",max)
					except Exception as e:
						print("NOT PROCESSED "+i)
						proc_fail=proc_fail+1
						pass
				try :
					for k in max:
						try :
							print(k)
							temp_map = AlignFrag.GetDistanceFrag("./sel_FRAG/"+k[0],INIT_PDB,min_ss,max_ss)
							#map_frag={**map_frag,**temp_map}
							#map_frag.update(temp_map)
							for i in temp_map.keys():
								if map_frag.get(i)==None:
									map_frag[i]=temp_map[i]
						except Exception as e:
							print(e)
				except Exception as e:
					print("problem with:"+i)
					print(e)

			#print(len(temp_map))
			print("Computing distance restraint: Global")


			#-----------------------DISTANCE map TEMPLATE
			lCA_map_temp = {}
			frg_fold = os.listdir(glblist)

			for i in frg_fold:
				with open(glblist + i,"r") as s1:
					t = s1.readlines()[1:]
					try:
						max = takeglob.Max20(t)[0:1]
						#print("Sel Glob: ",max)
					except Exception as e:
						print("NOT PROCESSED "+i)
						proc_fail=proc_fail+1
						pass
				try :
					for k in max:
						try :
							temp_map = AlignFrag.GetDistanceGlob("./sel_GLOB/"+k[0],INIT_PDB)
							#map_frag={**map_frag,**temp_map}
							#map_frag.update(temp_map)
							for i in temp_map.keys():
								if lCA_map_temp.get(i)==None:
									lCA_map_temp[i]=temp_map[i]
						except Exception as e:
							print(e)
				except Exception as e:
					print("problem with:"+i)
					print(e)

			#print(len(temp_map))
			

			

			print("Computing distance restraint: InitModel")
			
			top_init = model.getTopology()
			atom_init = top_init.atoms()
			pos_init = model.getPositions()

			lCA_init = []
			lCA_map_init = []

			initt= io.LoadPDB(INIT_PDB,calpha_only=True).Select("peptide=true")
			initt_CA = initt.residues
			nCA=int(str(initt_CA[0].number))
			print(nCA)
			for i,j in zip(atom_init,pos_init):
				if i.name == "CA": 
					lCA_init.append([i.index,j,nCA])
					#print(i.index)
					nCA=nCA+1



			for k in iter.combinations(lCA_init,2):
				
				d = dist(k[0][1],k[1][1])
						
						#print(k[0],l[0],d)
				if(k[0][1] != None and k[1][1] !=None):
					#if(map_frag.get((k[2],l[2]))!=None):    
						#print(k[2],l[2],map_frag.get((k[2],l[2])))
					#lCA_map_init.append([k[0][0],k[1][0],round(dist(k[0][1],k[1][1]),3),lCA_map_temp.get((k[0][2],k[1][2])),map_frag.get((k[0][2],k[1][2])),k[0][2],k[1][2]])
					lCA_map_init.append([k[0][0],k[1][0],round(dist(k[0][1],k[1][1]),3),lCA_map_temp.get((k[0][2],k[1][2])),map_frag.get((k[0][2],k[1][2])),k[0][2],k[1][2]])
			
			#model.topology.setUnitCellDimensions((30,30,30))

			#model.addSolvent(forcefield,neutralize=False)
			#model.addSolvent(forcefield,neutralize=True)
			#model.addHydrogens()
			#unm = forcefield.getUnmatchedResidues(model.topology)
			
			system = forcefield.createSystem(model.topology, nonbondedMethod=NoCutoff,nonbondedCutoff=1*un.nanometer)
			#system = forcefield.createSystem(model.topology, nonbondedMethod=NoCutoff,ignoreExternalBonds = True)
			#system = forcefield.createSystem(model.topology, nonbondedCutoff=1*un.nanometer)


			

			
			#c_forc_alfa = CustomBondForce("step(ref-r)*k1*(r-r1)^2 + step(ref-r)*k2*(r-r2)^2 + step(r3)*k3*(r-r3)^2") #tutte le forze
			#c_forc_alfa = CustomBondForce("step(ref-r)*k1*(r-r1)^2 + step(ref-r)*k2*(r-r2)^2") #no fragment
			c_forc_alfa_frag = CustomBondForce("step(r3)*k3*(r-r3)^2") #solo fragment
			c_forc_alfa = CustomBondForce("k2*step(r2)*(r-r2)^2") #diverso
			c_forc_alfa_init = CustomBondForce("k1*(r-r1)^2")
			c_forc_clash = CustomBondForce("step(d0-r)*k4*(d0-r)^2")
			#init template    
			scale=  1
			c_forc_alfa_init.addPerBondParameter("r1")
			c_forc_alfa_init.setForceGroup(1)
			c_forc_alfa_init.addGlobalParameter("k1",scale*0.5*(un.kilocalories_per_mole/(un.angstrom**2)))
			c_forc_alfa_init.addGlobalParameter("ref",1.5 )
			#c_forc_alfa_init.addGlobalParameter("k1",scale*0.5*un.kilocalorie_per_mole/(un.angstrom**2))

			#glob template
			c_forc_alfa.addPerBondParameter("r2")
			c_forc_alfa.setForceGroup(2)
			c_forc_alfa.addGlobalParameter("k2",scale*0.5*(un.kilocalories_per_mole/(un.angstrom**2)))
			c_forc_alfa.addGlobalParameter("ref",1.5 )

			#frag force
			c_forc_alfa_frag.addGlobalParameter("k3",scale*2*(un.kilocalories_per_mole/(un.angstrom**2)))
			c_forc_alfa_frag.addPerBondParameter("r3")
			c_forc_alfa_frag.setForceGroup(3)

			#Clash force
			#c_forc_clash.addGlobalParameter("k4",scale*83600) 
			c_forc_clash.addGlobalParameter("k4",scale*200*un.kilocalories_per_mole/(un.angstrom**2))
			c_forc_clash.addGlobalParameter("d0",0.36)
			c_forc_clash.setForceGroup(4)

			print("Adding Restraint Force")
			fail = 0
			frag_fail=0
			count_bond = 0
			
			for t in lCA_map_init: 
				r1 = t[2]
				r2 = t[3]
				
				
				if t[4]!= None:
					r3=t[4]
					#print("OK",t[1],t[2])
					
				else:
					r3= -1
					frag_fail=frag_fail+1
					
				if t[3]!= None:
					r2=t[3]
					#print("OK",t[1],t[2])
					
				else:
					r2=-1
					
				try:
					
					#print("ok")
					if r1<0.36:
						c_forc_clash.addBond(t[0],t[1])    
					
					if r1<1.5 and r3 !=-1 :
							c_forc_alfa_frag.addBond(t[0],t[1],[r3])
					#        print(t[len(t)-2],t[len(t)-1],r1,r2,r3)
							
					#if r1 <1.5 and r2 !=-1:
					#        c_forc_alfa.addBond(t[0],t[1],[r2]) 
					
					if(r1 < 1.5):
						
						c_forc_alfa_init.addBond(t[0],t[1],[r1]) #forcefield force
						
						
						
				except Exception as e:
					fail = fail+1
					print("oh no",t[0],t[1],[r1,r2,r3])
					
				
				is_in_dom = False
			
			
			print("DRMSD global template: ",drmsd(lCA_map_init))
			#print("DRMSD fragmental template":,drmsd(LCA_init))
			print("Restraint  fail:",fail)
			#print("Fragment Restraint fail: ",frag_fail)

			system.addForce(c_forc_alfa_init)
			system.addForce(c_forc_alfa)
			system.addForce(c_forc_alfa_frag)
			system.addForce(c_forc_clash)

			#for i in range(0,len(system.getForces())):
			#    print(system.getForce(i),i)
			sysforc = system.getForces()
			#for i in range(0,len(sysforc)):
			#    print(sysforc[i])
			hforce = system.getForce(4)

			#print(hforce)
			#for k in range(0,5):
			#    print(hforce.getPerDonorParameterName(k))
				
			"""
			scale_h=10
			hforce.setGlobalParameterDefaultValue(0,scale_h*2*un.kilocalories_per_mole/(un.angstrom**2))
			hforce.setGlobalParameterDefaultValue(1,scale_h*0.5*un.kilocalories_per_mole/(un.degrees**2))
			hforce.setGlobalParameterDefaultValue(3,scale_h*0.5*un.kilocalories_per_mole/(un.degrees**2))

			"""
			#hforce.setCutoffDistance(3*un.angstrom)
			hforce.setForceGroup(5)
			sysforc = system.getForces()

			cmm = system.getForce(len(sysforc)-5)
			
			#print(sysforc)
			for i in sysforc[0:4]:
				print(i)
				i.setForceGroup(6)
			print("Fragmental Restraints: ",sysforc[len(sysforc)-2].getNumBonds())
			print("Ca clash: ",sysforc[len(sysforc)-1].getNumBonds())
			#for i in range(0,6):
			#    print(system.getForce(i),i
			#system.removeForce(3)


			#print(system.getForces()[6].getNumBonds())
			#system.removeForce(1)
			#system.removeForce(1)
			#system.removeForce(1)
			#ystem.removeForce(len(system.getForces())-2)
			#print(system.getForces())




			TEMPERATURE = 100
			#integrator = VerletIntegrator( 1*un.femtosecond)
			#integrator = aux.integrators.VelocityVerletIntegrator(timestep = 2*un.femtosecond)  #Interesting
			#integrator = aux.integrators.AndersenVelocityVerletIntegrator(TEMPERATURE*un.kelvin, 91/un.picosecond,2*un.femtoseconds)
			#integrator = NoseHooverIntegrator(TEMPERATURE*un.kelvin, 91/un.picosecond,2*un.femtoseconds)
			#integrator = VerletIntegrator(0.001*un.picoseconds)
			#integrator = aux.integrator.GeodesicBAOABIntegrator(TEMPERATURE,1/un.picoseconds,1*un.femtoseconds)
			#integrator = aux.integrators.MetropolisMonteCarloIntegrator(temperature=TEMPERATURE)
			integrator = LangevinMiddleIntegrator(TEMPERATURE*un.kelvin, 1/un.picosecond, 2*un.femtoseconds)
			#integrator = integrators.HMCIntegrator(temperature=TEMPERATURE,nsteps=10)
			pressure = 1*un.atmospheres
			barostatInterval = 10
			#barostat = MonteCarloBarostat(1*un.atmospheres,TEMPERATURE*un.kelvin)
			#thermostat = AndersenThermostat(TEMPERATURE*un.kelvin,1/un.picoseconds)
			#system.addForce(barostat)
			#system.addForce(thermostat)
			#model_init = Modeller(init.topology, init.positions) 
			#integrator= GeodesicBAOABIntegrator(0*un.kelvin,1/un.picoseconds, 4 * un.femtosecond)
			#nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
			#nonbonded.setReciprocalSpaceForceGroup(7)
			#integrator = MTSLangevinIntegrator(273*un.kelvin, 1/un.picosecond, 4*un.femtoseconds, [(1,1),(2,1),(3,1),(4,1),(6,4),(7,4)])
			#integrator = aux.integrators.AndersenVelocityVerletIntegrator(temperature=TEMPERATURE*un.kelvin, timestep=1*un.femtosecond,collision_rate=91 / un.picoseconds)
			#integrator.setIntegrationForceGroups({6,1})
			#barostat = AndersenThermostat(0*un.kelvin, 1/un.picosecond)
			#system.addForce(barostat)
			#model_init.addSolvent(forcefield,neutralize=False)
			print("Setting simulation context..")


			#4000 step mi equilibro
			simulation = Simulation(model.topology, system, integrator,platform,properties)
			#simulation.context.computeVirtualSites()
			simulation.context.setPositions(model.positions)
			
			#print("save checkpoint..")
			#simulation.saveState("initGRAD")
			#simulation.saveCheckpoint("check1")
			#barostat = MonteCarloBarostat(pressure, 0, barostatInterval)
			#system.addForce(barostat)
			
			simulation.reporters.append(PDBReporter("Traj_FGMD/"+Target, 1000))
			simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, elapsedTime=True ))




			# In[ ]:

			#print(simulation.context.getState(getForces=True).getForces())
			#simulation.saveCheckpoint("firstMin")
			print("-----------")
			print("Name: ",OUT_NAME)
			print("Initial template: ", INIT_PDB)

			print("-----------")
				
			print("Start minimization")
			

			#state = simulation.context.getState(getPositions = True, getEnergy=True)
			#print("System En: ",state.getPotentialEnergy())

			
			
			simulation.minimizeEnergy()
			#LocalEnergyMinimizer.minimize(simulation.context, 1e-7)
			#simulation.loadState("old.xml")
			#LocalEnergyMinimizer.minimize(simulation.context, 1e-6,maxIterations=5000)
			#barostat.setDefaultTemperature(temp*un.kelvin)
			#integrator.setTemperature(TEMPERATURE*un.kelvin)
			#simulation.context.setVelocitiesToTemperature(TEMPERATURE*un.kelvin)
			simulation.context.setVelocitiesToTemperature(TEMPERATURE*un.kelvin)
			state = simulation.context.getState(getPositions = True).getPositions()
			a=Modeller(simulation.topology,state)
			a.deleteWater()
			PDBFile.writeFile(a.topology,a.positions,open("Grad_FGMD/"+Target,"w"))
			simulation.step(1000)
			
			for i in range(1,TEMPERATURE,1)[::-1]:
				integrator.setTemperature(i*un.kelvin)
				#barostat.setDefaultTemperature(i*un.kelvin)
				#thermostat.setDefaultTemperature(i*un.kelvin)
				simulation.step(100)
			
			simulation.step(1000)
			
				
			state = simulation.context.getState(getPositions = True).getPositions()
			a=Modeller(simulation.topology,state)
			a.deleteWater()
			PDBFile.writeFile(a.topology,a.positions,open(OUT+"/"+Target,"w"))
			
			for i in range(0,3):

				state_old = simulation.context.getState(getPositions = True, getEnergy=True)
				simulation.saveState("old.xml")

				LocalEnergyMinimizer.minimize(simulation.context, 1e-6,maxIterations=5000)
				state_new = simulation.context.getState(getPositions = True, getEnergy=True)
				simulation.saveState("prev.xml")

				if state_new.getPotentialEnergy() < state_old.getPotentialEnergy(): 
					
					#print("System En: - saving",state_new.getPotentialEnergy() )
					a=Modeller(simulation.topology,state_new.getPositions())
					a.deleteWater()
					PDBFile.writeFile(a.topology,a.positions,open(OUT+"/"+Target,"w"))
					#LocalEnergyMinimizer.minimize(simulation.context, 1e-4,maxIterations=1000)
					
					
				
			state = simulation.context.getState(getPositions = True, getEnergy=True)
			a=Modeller(simulation.topology,state.getPositions())
			a.deleteWater()
			PDBFile.writeFile(a.topology,a.positions,open("Gradloc_FGMD/"+Target,"w"))
		
		except Exception as e:
      
			
			with open(OUT+"/"+Target+".fail","w") as f:
				f.write(str(e))
				
			f.close()
		








