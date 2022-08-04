import sys
import os.path as path
import os

suite = sys.argv[1]
dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns"

for file in sorted(os.listdir(path.join(dir, suite))):
    if len(file) < 4 or file[:4] != "Halo": continue
    print(file)
    os.system("python3 subhalo_positions.py %s %s %s example_plots" % (dir, suite, file))
    os.system("python3 subhalo_MAHs.py %s %s %s example_plots" % (dir, suite, file))
