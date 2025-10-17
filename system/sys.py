import hoomd
from hoomd import md
import numpy as np
import json
import copy


class System(object):
    def __init__(self, sim, protein, protein4G, param):
        self.sim = sim
        self.protein = protein
        self.protein4G = protein4G
        self.param = copy.copy(param)
        self.system_types = ['A', 'U', 'C', 'G', 'AT', 'Cap', 'protE', 'Ea', 'Eb', 'Ec', 'Ed', 'protG', 'Ga', 'Gb', 'Gc', 'M']
        self.particlesizes = [param['sizeC'], param['sizeC'], param['sizeC'], param['sizeC'], param['sizeC'], param['sizeCap'], param['sizeE'], param['sizeEa'], param['sizeEb'], param['sizeEc'], param['sizeEd'], param['sizeG'], param['sizeGa'], param['sizeGb'], param['sizeGc'], param['memcsize']]
        assert(len(self.system_types) == len(self.particlesizes))
        self.nl = hoomd.md.nlist.Cell(buffer = 0.4, exclusions = ('body', 'constraint'))
        
    
    def setup_object(self, create = True):
        self.rigid = hoomd.md.constrain.Rigid()
        self.rigid.body['protE'] = {
                "constituent_types": self.protein.constituent_types,
                "positions": self.protein.constituents,
                "charges": self.protein.constituent_charges,
                "orientations": self.protein.constituent_orientations,
                "diameters": self.protein.constituent_diameters
                }
        self.rigid.body['protG'] = {
                "constituent_types": self.protein4G.constituent_types,
                "positions": self.protein4G.constituents,
                "charges": self.protein4G.constituent_charges,
                "orientations": self.protein4G.constituent_orientations,
                "diameters": self.protein4G.constituent_diameters
                }
        if create:
            self.rigid.create_bodies(self.sim.state)


    def setup_potential(self):
        self.lj = hoomd.md.pair.LJ(self.nl, default_r_cut = self.param['lj_cut'])
        self.lj.mode = 'shift'
        self.yukawa = hoomd.md.pair.Yukawa(self.nl, default_r_cut = self.param['yukawa_cut'])
        self.yukawa.mode = 'xplor'
        self.bond = hoomd.md.bond.Harmonic()
        

    def setup_potential_coef(self):
        for i in range(len(self.system_types)):
            for j in range(len(self.system_types)):
                if i <= j:
                    type1, size1 = self.system_types[i], self.particlesizes[i]
                    type2, size2 = self.system_types[j], self.particlesizes[j]
                    dis = 0.5*(size1+size2)
                    lj = self.param.get('lj_'+type1+type2, 0.) + self.param.get('lj_'+type2+type1, 0.)
                    if type1 == type2:
                        lj = 0.5 * lj
                    if lj > 0.:
                        self.lj.params[(type1, type2)] = dict(sigma = dis*2**(-1/6), epsilon = lj)
                        self.lj.r_cut[(type1, type2)] = dis + self.param['lj_cut']
                        print(type1, type2, dis, 'interact', lj, ' ron', self.lj.r_on[(type1, type2)], ' rcut', self.lj.r_cut[(type1, type2)])
                    else:
                        self.lj.params[(type1, type2)] = dict(sigma = dis*2**(-1/6), epsilon = self.param['lj_repel'])
                        self.lj.r_cut[(type1, type2)] = dis
                        print(type1, type2, dis, 'exclude', lj, ' ron', self.lj.r_on[(type1, type2)], ' rcut', self.lj.r_cut[(type1, type2)])
            
        for i in range(len(self.system_types)):
            for j in range(len(self.system_types)):
                if i <= j:
                    type1, size1 = self.system_types[i], self.particlesizes[i]
                    type2, size2 = self.system_types[j], self.particlesizes[j]
                    dis = 0.5*(size1+size2)
                    kappa = self.param['yukawa_kp']
                    m_epsilon = self.param['yukawa_lb'] * self.param.get('Z_'+type1, 0.) * self.param.get('Z_'+type2, 0.) * (np.exp(kappa * 0.5*size1) / (1 + kappa * 0.5*size1)) * (np.exp(kappa * 0.5*size2) / (1 + kappa * 0.5*size2))
                    if m_epsilon != 0.:
                        print(type1, type2, dis, 'ele_interact', m_epsilon, 'E(dis)=', m_epsilon*np.exp(-self.param['yukawa_kp']*dis)/dis)
                        self.yukawa.params[(type1, type2)] = dict(epsilon = m_epsilon, kappa = self.param['yukawa_kp'])
                        self.yukawa.r_cut[(type1, type2)] = self.param['yukawa_cut']
                    else:
                        print(type1, type2, dis, 'ele_none', m_epsilon)
                        self.yukawa.params[(type1, type2)] = dict(epsilon = 0, kappa = 0)
                        self.yukawa.r_cut[(type1, type2)] = dis
     
        
        self.bond.params['backbone'] = dict(k = self.param['backbonek'], r0 = self.param['sizeC'])
        self.bond.params['basepair'] = dict(k = self.param['basepairk'], r0 = self.param['sizeC'])
    
    def setup_operation(self):
        gsd_writer = hoomd.write.GSD(filename="md.gsd", filter=hoomd.filter.All(), trigger=hoomd.trigger.Periodic(self.param['record_period']), mode='ab', dynamic=['property', 'momentum', 'attribute', 'topology'])
        self.sim.operations.writers.append(gsd_writer)

        thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(
            filter=hoomd.filter.All())
        self.sim.operations.computes.append(thermodynamic_properties)

        logger = hoomd.logging.Logger(categories=['scalar', 'string'])
        logger.add(thermodynamic_properties, quantities=['potential_energy'])
        logger.add(self.sim, quantities=['timestep', 'walltime'])
        logger.add(self.lj, quantities=['energy'])
        logger.add(self.yukawa, quantities=['energy'])
        logger.add(self.bond, quantities=['energy'])
        file = open('log.txt', mode='a', newline='\n')
        table_file = hoomd.write.Table(output=file,
            trigger=hoomd.trigger.Periodic(self.param['record_period']),
            logger=logger)
        self.sim.operations.writers.append(table_file)


    def setup_integrator(self):
        langevin = hoomd.md.methods.Langevin(filter = hoomd.filter.Type(['A', 'U', 'C', 'G', 'AT', 'Cap', 'protE', 'protG']), kT = self.param['kT'])

        self.sim.operations.integrator = hoomd.md.Integrator(dt = self.param['dt'], methods = [langevin], forces = [self.lj, self.yukawa, self.bond], integrate_rotational_dof = True, rigid = self.rigid)
