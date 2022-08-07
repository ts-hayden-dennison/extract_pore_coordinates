import mdtraj
import sys
import csv
import argparse
from shapely.geometry import Polygon, LineString
from shapely.validation import make_valid
from shapely.geometry.polygon import orient

def extract_yz(trajfile, topfile, residue, atom, chainlist):
	
	data = []
	for i in range(0, len(chainlist)):
		data.append([])
	
	chains = ' '.join([str(x) for x in chainlist])
	
	top = mdtraj.load(topfile).topology
	traj = mdtraj.load_xtc(trajfile, top)
	if type(atom) == str:
		indexes = top.select('name {} and residue {} and chainid '.format(atom, residue) + chains)
	else:
		indexes = top.select('index {} and residue {} and chainid '.format(atom, residue) + chains)
	
	for frame in traj:
		r = frame.xyz[0][indexes]
		for i, xyz in enumerate(r):
			data[i].append([xyz[1], xyz[2]])
	
	return data


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Extract specific atom coordinates from a trajectory and prepare them \
									 for analysis")
	parser.add_argument('-f', '--trajectory', type=str, help='The trajectory file to analyze', required=True)
	parser.add_argument('-s', '--topology', type=str, help='The topology file to use', required=True)
	parser.add_argument('-r', '--residue', type=int, help='The residue number to pick', required=True)
	
	
	atomgroup = parser.add_mutually_exclusive_group(required=True)
	atomgroup.add_argument('-an', '--atomn', type=str, help='The atom name on the specific residue to pick', \
							default=None)
	atomgroup.add_argument('-ai', '--atomi', type=int, help='The atom index on the specific residue to pick', \
							default=None)
	
	parser.add_argument('-c', '--chains', type=int, help='The chain indexes to pick from', nargs='+', required=True)
	parser.add_argument('-o', '--output', type=str, help='The prefix of the output files', default='f', \
						nargs='?')
	args = parser.parse_args()
	print(vars(args))
	if args.atomn != None:
		data = extract_yz(args.trajectory, args.topology, args.residue, args.atomn, args.chains)
	else:
		data = extract_yz(args.trajectory, args.topology, args.residue, args.atomi, args.chains)
	
	
	# need to convert data from being sorted by chain to being sorted by frame
	sorted_data = []
	for frame in range(0, len(data[0])):
		sorted_data.append([])
		for chaini, chain in enumerate(data):
			sorted_data[-1].append(chain[frame])
	
	newdata = []
	for frame in sorted_data:
		# perform largest pore extraction
		pore = make_valid(Polygon(frame))
		if pore.geom_type == 'MultiPolygon':
			largest = 0
			p = None
			for i, poreseg in enumerate(pore.geoms):
				if poreseg.area > largest:
					p = Polygon(poreseg)
					largest = poreseg.area
		else:
			p = pore
		newdata.append(list(p.exterior.coords)[:-1])
		
	#output as many files as there are frames
	for frame in range(0, len(newdata)):
		with open(args.output+'{}.dat'.format(frame+1), 'w') as f:
			s = ''
			for xy in newdata[frame]:
				s = s + '{} {}\n'.format(*xy)
			f.write(s)
			f.close()
	sys.exit()
