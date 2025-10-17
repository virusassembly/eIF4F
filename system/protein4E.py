from __future__ import division
import numpy as np


class CreateBody():

    def __init__(self, param):
        super(CreateBody, self).__init__()

        self.name = 'pro4E'
        self.size = param['sizeE']
        self.type = 'protE'        
        self.orientation = (1, 0, 0, 0)


        a = 0.6   
        self.ligandA_sites = [[a,0,-a], [-a,0,-a]]
        self.ligandB_sites = [[0,0,0]]
        self.ligandC_sites = [[a,a,0], [0,a,0], [-a,a,0], 
                              [a,0,0], [-a,0,0], 
                              [a,-a,0], [0,-a,0], [-a,-a,0]]
        self.ligandD_sites = [[a,a,a], [0,a,a], [-a,a,a], 
                              [a,0,a], [0,0,a], [-a,0,a], 
                              [a,-a,a], [0,-a,a], [-a,-a,a]]
        
        self.constituents = self.ligandA_sites + self.ligandB_sites + self.ligandC_sites + self.ligandD_sites
        

        self.constituent_types = ['Ea'] * len(self.ligandA_sites) + ['Eb'] * len(self.ligandB_sites) + ['Ec'] * len(self.ligandC_sites) + ['Ed'] * len(self.ligandD_sites)

        self.constituent_diameters = [param['sizeEa']] * len(self.ligandA_sites) + [param['sizeEb']] * len(self.ligandB_sites) + [param['sizeEc']] * len(self.ligandC_sites) + [param['sizeEd']] * len(self.ligandD_sites)
        

        if param['extralayer'] > 0:
            layer = [[a,a,2*a], [0,a,2*a], [-a,a,2*a], 
                    [a,0,2*a], [0,0,2*a], [-a,0,2*a], 
                    [a,-a,2*a], [0,-a,2*a], [-a,-a,2*a]]
            self.constituents = self.constituents + layer
            self.constituent_types = self.constituent_types + ['Ga'] * len(layer)
            self.constituent_diameters = self.constituent_diameters + [param['sizeGa']] * len(layer)
        if param['extralayer'] > 1:
            layer =[[a,a,3*a], [0,a,3*a], [-a,a,3*a], 
               [a,0,3*a], [0,0,3*a], [-a,0,3*a], 
               [a,-a,3*a], [0,-a,3*a], [-a,-a,3*a]]
            self.constituents = self.constituents + layer
            self.constituent_types = self.constituent_types + ['Gb'] * len(layer)
            self.constituent_diameters = self.constituent_diameters + [param['sizeGb']] * len(layer)
        if param['extralayer'] > 2:
            layer =[[a,a,4*a], [0,a,4*a], [-a,a,4*a], 
               [a,0,4*a], [0,0,4*a], [-a,0,4*a], 
               [a,-a,4*a], [0,-a,4*a], [-a,-a,4*a]]
            self.constituents = self.constituents + layer
            self.constituent_types = self.constituent_types + ['Gc'] * len(layer)
            self.constituent_diameters = self.constituent_diameters + [param['sizeGc']] * len(layer)
        
        print("Protein contains ", len(self.constituents), " particles")
        mass = 0.1
        I = np.zeros(shape=(3, 3))
        for r in self.constituents:
            I += mass * (np.dot(r, r) * np.identity(3) - np.outer(r, r))
        self.mass = mass*len(self.constituents)
        self.moment_of_inertia = [I[0,0], I[1,1], I[2,2]]
        print("moment of inertia:", I[0,0], I[1,1], I[2,2], "mass:", self.mass)

        
        self.constituent_charges = [0.] * len(self.constituents)
        self.constituent_orientations = [(1, 0, 0, 0)] * len(self.constituents)


        print("test protein:", self.constituents)
        print("test protein:", self.constituent_types)
        print("test protein:", self.constituent_diameters)
        print("test protein:", self.constituent_charges)
        print("test protein:", self.constituent_orientations)