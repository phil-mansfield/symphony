import numpy as np
import symlib
import matplotlib.pyplot as plt
import palette
from palette import pc
palette.configure(True)

sim_dir = symlib.get_host_directory(
    "/sdf/home/p/phil1/ZoomIns",
    "SymphonyLCluster", "Halo_047")

scale = symlib.scale_factors(sim_dir)
s, hist = symlib.read_symfind(sim_dir)
r, _ = symlib.read_rockstar(sim_dir)

fig, ax = plt.subplots()
colors = ["tab:red", "tab:orange", "tab:green",
          "tab:blue", "tab:purple"]

# First, plot the host
ok = r[0,:]["ok"]
ax.plot(scale[ok], r[0,ok]["m"], c="k")

# For now, let's only plot minor mergers which dirupt before
# the end of the simulation in the Rockstar catalogue.
is_target = (hist["merger_ratio"] < 0.1) & (~r["ok"][:,-1])
targets = np.where(is_target)[0][:5]
for i_color, i in enumerate(targets):
    print(i_color, i)
    # Plot the mass history of the rockstar subhalo
    ok = r[i,:]["ok"]
    ax.plot(scale[ok], r[i,ok]["m"], c=colors[i_color])
    # Plot the mass history of the symfind subhalo
    ok = s[i,:]["ok"]
    ax.plot(scale[ok], s[i,ok]["m_raw"], c=colors[i_color], lw=1.5)

ax.set_yscale("log")
ax.set_xlabel(r"$a(t)$")
ax.set_ylabel(r"$m(t)$")
fig.savefig("mah.png")
