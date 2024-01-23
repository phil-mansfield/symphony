import symlib
import numpy as np

base_dir = "/sdf/home/p/phil1/ZoomIns"
suite = "MWest"
n_hosts = symlib.n_hosts(suite)

n_sub_LMC = 0
n_sub_tot = 0

for i_host in range(n_hosts):
    sim_dir = symlib.get_host_directory(base_dir, suite, i_host)

    r, hist = symlib.read_rockstar(sim_dir)
    host_idx = symlib.pre_infall_host(hist)
    i_LMC = np.where(r["ok"][:,-1])[0][1]
    
    is_LMC_sub = host_idx == i_LMC

    n_sub_LMC += np.sum(is_LMC_sub)
    n_sub_tot += len(hist) - 2 # Don't count host and LMC

print("Average LMC subhalo fraction:")
print("%.3f" % (n_sub_LMC / n_sub_tot))
