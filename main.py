import sys
import hoomd
import system
import numpy as np
import json, glob, math
import itertools
import gsd.hoomd



def snapprint(snap):
    np.set_printoptions(threshold=np.inf)
    print('particle position:\n{} \
            \nparticle typeid:\n{} \
            \npartical body:\n{} \
            \nparticale diameter:\n{}' \
            .format(np.array(snap.particles.position[:]),np.array(snap.particles.typeid[:]),np.array(snap.particles.body[:]),np.array(snap.particles.diameter[:])))
    print('bond groups:\n{} \
            \nbond typeid:\n{}' \
            .format(np.array(snap.bonds.group[:]),np.array(snap.bonds.typeid[:])))
    print('angle groups:\n{} \
            \nangle typeid:\n{}' \
            .format(np.array(snap.angles.group[:]),np.array(snap.angles.typeid[:])))
    print('dihedral groups:\n{} \
            \ndihedral typeid:\n{}' \
            .format(np.array(snap.dihedrals.group[:]),np.array(snap.dihedrals.typeid[:])))

class lattice():
    def __init__(self, obj, param, shift=0.):
        self.n = param['N_'+obj.name]
        self.L = param['boxL']
        if self.n > 0:
            K = math.ceil(self.n**(1 / 3)) 
            spacing = math.floor(self.L / K)
            print(K, self.L)
            x = np.linspace(-self.L / 2, self.L / 2, K, endpoint=False) + 0.5*spacing + 5.
            y = np.linspace(-self.L / 2, self.L / 2, K, endpoint=False) + 0.5*spacing
            z = np.linspace(0., self.L / 2-5., K, endpoint=False) + 0.5*spacing - 2.
            allpos = np.array(list(itertools.product(x, y, z)))
            
            self.particles = [[x, y, z+shift] for x,y,z in allpos[allpos[:,2]>-self.L/2+5.][0: self.n]]
        else:
            self.particles = []
        self.particletypes = [obj.type] * self.n
        self.masses = [obj.mass] * self.n
        self.moment_inertias = [obj.moment_of_inertia] * self.n
        self.diameters = [obj.size] * self.n
        self.bodyids = [-1] * self.n



def run():
    gpu = hoomd.device.GPU()

    gsd_file = gsd.hoomd.open('../relaxed_'+param['rnafile']+'.gsd')
    filesnap = gsd_file[0]
    snap = hoomd.Snapshot()
    snap.particles.N = filesnap.particles.N
    snap.particles.types = filesnap.particles.types
    snap.particles.diameter[:] = filesnap.particles.diameter[:]
    snap.particles.position[:] = filesnap.particles.position[:]
    snap.particles.typeid[:] = filesnap.particles.typeid[:]
    snap.particles.body[:] = filesnap.particles.body[:]
    snap.particles.moment_inertia[:] = filesnap.particles.moment_inertia[:]
    snap.particles.mass[:] = filesnap.particles.mass[:]
    snap.particles.angmom[:] = filesnap.particles.angmom[:]
    
    snap.bonds.types = filesnap.bonds.types
    snap.bonds.N = filesnap.bonds.N
    snap.bonds.group[:] = filesnap.bonds.group[:]
    snap.bonds.typeid[:] = filesnap.bonds.typeid[:]               


    typedic = {'A':0, 'U':1, 'C':2, 'G':3, 'AT':4, 'Cap':5, 'protE':6, 'Ea':7, 'Eb':8, 'Ec':9, 'Ed':10, 'protG':11, 'Ga':12, 'Gb':13, 'Gc':14, 'M':15, 'backbone':0, 'basepair':1}


    protein = system.protein4E.CreateBody(param)
    lat = lattice(protein, param)
    snap.configuration.box = [lat.L, lat.L, lat.L, 0, 0, 0]
    system.snap.Update_Snap(snap, lat, typedic, bond = False, angle = False, dihedral = False)

    protein4G = system.protein4G.CreateBody(param)
    lat = lattice(protein4G, param, 2.)
    system.snap.Update_Snap(snap, lat, typedic, bond = False, angle = False, dihedral = False)

    
    sim = hoomd.Simulation(device = gpu, seed = param['seed'])
    sim.create_state_from_snapshot(snap)
    
    chamber = system.sys.System(sim, protein, protein4G, param)
    print('sim.particles.types:', chamber.sim.state.particle_types)

    chamber.setup_object()
    chamber.setup_potential()
    chamber.setup_potential_coef()
    chamber.setup_operation()
    chamber.setup_integrator()

    chamber.sim.run(param['runtime'], write_at_start = True)

def continue_run(file):
    gpu = hoomd.device.GPU()

    print("continue running")
    

    protein = system.protein4E.CreateBody(param)
    protein4G = system.protein4G.CreateBody(param)
    
    sim = hoomd.Simulation(device=gpu, seed=param['seed'])
    sim.create_state_from_gsd(file, frame = -1)
    chamber = system.sys.System(sim, protein, protein4G, param)

    chamber.setup_object(False)
    chamber.setup_potential()
    chamber.setup_potential_coef()
    chamber.setup_operation()
    chamber.setup_integrator()

    chamber.sim.run(param['runtime'], write_at_start= True)


param=sys.argv[1]
print(param)
with open(param) as f: param=json.load(f)

if glob.glob(param['gsdfile']+'.gsd'):
    continue_run(param['gsdfile']+'.gsd')
else:
    run()