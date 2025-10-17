from __future__ import division
import numpy as np


class CreateBody():

    def __init__(self, param):
        super(CreateBody, self).__init__()

        self.name = 'pro4G'
        self.size = param['sizeG']
        self.type = 'protG'        
        self.orientation = (1, 0, 0, 0)


        a = 0.6
        self.layera = [[a,a,-a], [0,a,-a], [-a,a,-a], 
               [a,0,-a], [0,0,-a], [-a,0,-a], 
               [a,-a,-a], [0,-a,-a], [-a,-a,-a]]
        self.layerb =[[a,a,0], [0,a,0], [-a,a,0], 
               [a,0,0], [0,0,0], [-a,0,0], 
               [a,-a,0], [0,-a,0], [-a,-a,0]]
        self.layerc = [[a,a,a], [0,a,a], [-a,a,a], 
               [a,0,a], [0,0,a], [-a,0,a], 
               [a,-a,a], [0,-a,a], [-a,-a,a]]


        self.constituents = self.layera + self.layerb + self.layerc
        print("Protein contains ", len(self.constituents), " particles")

        mass = 0.1
        I = np.zeros(shape=(3, 3))
        for r in self.constituents:
            I += mass * (np.dot(r, r) * np.identity(3) - np.outer(r, r))
        self.mass = mass*len(self.constituents)
        self.moment_of_inertia = [I[0,0], I[1,1], I[2,2]]
        print("moment of inertia:", I[0,0], I[1,1], I[2,2], "mass:", self.mass)

        self.constituent_types = ['Ga'] * len(self.layera) + ['Gb'] * len(self.layerb) + ['Gc'] * len(self.layerc)
        self.constituent_diameters = [param['sizeGa']] * len(self.layera) + [param['sizeGb']] * len(self.layerb) + [param['sizeGc']] * len(self.layerc)
        self.constituent_charges = [0.] * len(self.constituents)
        self.constituent_orientations = [(1, 0, 0, 0)] * len(self.constituents)

        print("test protein:", self.constituents)
        print("test protein:", self.constituent_types)
        print("test protein:", self.constituent_diameters)
        print("test protein:", self.constituent_charges)
        print("test protein:", self.constituent_orientations)