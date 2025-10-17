import numpy as np

def assign_snap_particles(snap, obj, typedic):
	pstart = snap.particles.N
	pend = pstart+obj.n
	snap.particles.N = pend
	print("pstart, pend:",pstart,pend)
	print(snap.particles.position[0])

	for i, (cor, d, typename, bodyid, mi, mass) in enumerate(zip(obj.particles, obj.diameters, obj.particletypes, obj.bodyids, obj.moment_inertias, obj.masses), start = pstart):
		snap.particles.diameter[i] = d
		snap.particles.position[i] = cor
		snap.particles.typeid[i] = typedic[typename]
		snap.particles.body[i] = bodyid
		snap.particles.moment_inertia[i] = mi
		snap.particles.mass[i] = mass
		snap.particles.angmom[i] = (1,0,0,0)

	return pstart, pend

def assign_snap_bonds(snap, obj, typedic, pstart):
	bstart=snap.bonds.N
	snap.bonds.N = bstart+len(obj.bonds)
	
	for i, (bond, typename) in enumerate(zip(obj.bonds, obj.bondtypes), start = bstart):
		snap.bonds.group[i] = pstart + bond
		snap.bonds.typeid[i] = typedic[typename]

def assign_snap_angles(snap, obj, typedic, pstart):
	dstart = snap.angles.N
	snap.angles.N = dstart+len(obj.angles)
	
	for i, angle in enumerate(obj.angles,start=dstart):
		snap.angles.group[i] = pstart + angle
		snap.angles.typeid[i] = typedic[obj.atypename]

def assign_snap_dihedrals(snap, obj, typedic, pstart):
	dstart = snap.dihedrals.N
	snap.dihedrals.N = dstart+len(obj.dihedrals)
	
	for i, (dihedral, typename) in enumerate(zip(obj.dihedrals, obj.dtypes),start=dstart):
		snap.dihedrals.group[i] = pstart + dihedral
		snap.dihedrals.typeid[i] = typedic[typename]


def Update_Snap(snap, obj, typedic, bond = True, angle = False, dihedral = True):
	pstart,pend = assign_snap_particles(snap,obj,typedic)
	if bond:
		assign_snap_bonds(snap,obj,typedic,pstart)
	if angle:
		assign_snap_angles(snap,obj,typedic,pstart)
	if dihedral:
		assign_snap_dihedrals(snap,obj,typedic,pstart)
	return pstart, pend