import numpy as np
import symlib

base_dir = "/sdf/home/p/phil1/ZoomIns"
sim_dir = symlib.get_host_directory(base_dir, "SymphonyMilkyWay", 0)

h, hist = symlib.read_subhalos(sim_dir)
c = symlib.read_cores(sim_dir, suffix="fid2")

#idx = np.arange(len(h), dtype=int)
#print(np.where(~c["ok"][idx, hist["first_infall_snap"]])[0])
snap = np.arange(h.shape[1], dtype=int)
idxs = [250, 258, 386]
for i in idxs:
    start, end = hist["first_infall_snap"][i], np.max(snap[c["ok"][i]])
    skipped = snap[(~c["ok"][i]) & (snap >= start) & (snap <= end)]
    print("start: %d, end: %d, n_skip: %d" %
          (start, end, end-start+1 - np.sum(c["ok"][i])))
    print(skipped)

    print(i)
    for s in range(start, end):
        if s > start+5 and s < end - 6: continue
        r = np.sqrt(np.sum(c["x"][i,s]**2))
        vr = np.sum(c["x"][i,s]*c["v"][i,s]/r)
        print("%3d %6.2f %6.2f %6.2f %6.2f %8.3g %8.3g %5s %5s" %
              (s, r, vr, c["r50_bound"][i,s], c["vmax"][i,s],
               c["m_bound"][i,s], c["m_tidal"][i,s],
               str(c["ok"][i,s]), str(c["interp"][i,s])))
    print()
