
import openmm
import openmm.app as app
from  openmm.unit import *
from pdbfixer import PDBFixer 
from openmm.app import PDBFile, PDBReporter, DCDFile
from openmm.app import modeller
import os,time

def fix_pdb(pdb_id):
    path = os.getcwd()
    if len(pdb_id) != 4:
        print("Creating PDBFixer...")
        fixer = PDBFixer(pdb_id)
        print("Finding missing residues...")
        fixer.findMissingResidues()

        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                print("ok")
                del fixer.missingResidues[key]

        print("Finding nonstandard residues...")
        fixer.findNonstandardResidues()
        print("Replacing nonstandard residues...")
        fixer.replaceNonstandardResidues()
        print("Removing heterogens...")
        fixer.removeHeterogens(keepWater=False)

        print("Finding missing atoms...")
        fixer.findMissingAtoms()
        print("Adding missing atoms...")
        fixer.addMissingAtoms()
        print("Adding missing hydrogens...")
        fixer.addMissingHydrogens(7)
        print("Writing PDB file...")

        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            open(os.path.join(path, "%s_fixed_pH_%s.pdb" % (pdb_id.split('.')[0], 7)),
                 "w"),
            keepIds=True)
        return "%s_fixed_pH_%s.pdb" % (pdb_id.split('.')[0], 7)

os.getcwd()

# pre-fix
pdbfile =sys.argv[1] #'2jof.pdb'
T = float(sys.argv[2]) #300
steps = int(sys.argv[3]) #1e6
saveframes = 50
savepdbinterval = steps // saveframes
fix_pdb(pdbfile)

# system build
pdb = app.PDBFile("%s_fixed_pH_%s.pdb" % (pdbfile.split('.')[0], 7))
ff = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
mod = modeller.Modeller(pdb.topology, pdb.positions)
mod.addSolvent(forcefield=ff, model='tip3p', padding=1.2*nanometers, 
positiveIon='Na+', negativeIon='Cl-', ionicStrength=0.1*molar, neutralize=True)
PDBFile.writeFile(mod.topology, mod.positions, open('ions.pdb', 'w'))
system = ff.createSystem(mod.topology, 
                    nonbondedMethod=app.PME, 
                    nonbondedCutoff=1*nanometer, 
                    constraints=app.HBonds)
barostat = openmm.MonteCarloBarostat(1.0*bar, 5.0*kelvin, 25)
system.addForce(barostat)
temperature = T * kelvin
friction = 1.0 / picosecond
timestep = 2.0 * femtosecond
integrator = openmm.LangevinIntegrator(temperature, friction, timestep)
platform = openmm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}
simulation = app.Simulation(mod.topology, system, integrator, platform, properties)
simulation.context.setPositions(mod.positions)
# context = openmm.Simulation(system, integrator)
# context.setPositions(mod.positions)
# openmm.LocalEnergyMinimizer.minimize(context)
# context.setVelocitiesToTemperature(temperature)
simulation.minimizeEnergy()
min_pos = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(mod.topology, min_pos, open(f'min{T}.pdb', 'w'))
simulation.context.setVelocitiesToTemperature(temperature)

simulation.reporters.append(app.StateDataReporter(
    f'md{T}.log',
    1000,
    step=True,
    time=True,
    potentialEnergy=True,
    kineticEnergy=True,
    remainingTime=True,
    speed=True,
    totalSteps=steps,# totalSteps must be set for progress estimation
))
simulation.reporters.append(app.DCDReporter(f'md{T}.dcd', savepdbinterval//10))
simulation.reporters.append(app.PDBReporter(f'md{T}.pdb', savepdbinterval))
simulation.reporters.append(app.CheckpointReporter(f'checkpoint_file{T}.chk', savepdbinterval//10))  # save progress during the simulation


start_time = time.time()
print('simulation started at:\n%s'%time.ctime())
simulation.step(steps)
end_time = time.time()
print('simulation finished at:\n%s'%time.ctime())
deltat = end_time-start_time
print('time used: %s'%(deltat))
print('Efficency: %.2f ns/day'%(steps*2*float(2e-6)/deltat*3600*24))

