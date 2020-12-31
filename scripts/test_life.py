import glob
import numpy as np
import automata

# read in metadata file entries to a list of lists
fname = "./outfiles/metadata.txt"
file = open(fname, "r")
metadata = []
for line in file:
    line = line.split()
    metadata.append(line)

# parse metadata variables
rows = int(metadata[0][1])
cols = int(metadata[0][2])
num_procs = int(metadata[1][1])
iprocs = int(metadata[2][1])
jprocs = int(metadata[2][2])
num_gen = int(metadata[3][1])
is_periodic = int(metadata[4][1])
freq_dump = int(metadata[5][1])

paths = glob.glob('./outfiles/dump*')
paths = sorted(paths, key=lambda x: int(x.split("p")[1]))
states = []
for i, path in enumerate([paths[0], paths[-1]]):
    data = np.fromfile(path, dtype=np.int8)
    n = int(data.size / (rows*cols))
    offset = i*rows*cols*(n-1)
    state = np.reshape(data[offset:offset+rows*cols], (rows, cols))
    states.append(state)

# test final state against benchmark game of life simulator (James Percival, ACSE-1)
initial_state = states[0]
final_state = states[1]
print("Running benchmark test on initial state...")
test_state = automata.life(initial_state, num_gen, is_periodic)

if np.all(final_state == test_state):
    print("Test Passed")
else:
    print("Test Failed")
