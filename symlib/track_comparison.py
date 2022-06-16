import numpy as np
import matplotlib.pyplot as plt
import scipy
import lib
import palette
from palette import pc

def distance(x1, x2):
    return np.sqrt(np.sum((x1 - x2)**2, axis=1))

def main():
    import palette
    from palette import pc

    palette.configure(False)
    
    track_file = "/home/phil/code/src/github.com/phil-mansfield/symphony_pipeline/tmp_data/SymphonyMilkyWay/Halo023/halos/track_comparison.txt"
    suite = "SymphonyMilkyWay"
    sim_dir = "/home/phil/code/src/github.com/phil-mansfield/symphony_pipeline/tmp_data/SymphonyMilkyWay/Halo023"

    param = lib.parameter_table[suite]
    h, hist = lib.read_subhalos(param, sim_dir)

    mpeak = np.max(h["mvir"], axis=1)/param["h100"]
    rpeak = np.max(h["rvir"], axis=1)/param["h100"]*1e3
    
    cols = np.loadtxt(track_file).T
    idx = np.array(cols[0], dtype=int)

    rpeak, mpeak = rpeak[idx], mpeak[idx]
    x_rs, x_c1, x_c2 = cols[1:4].T, cols[4:7].T, cols[7:10].T
    r_rs, r_tidal, r_bound = cols[10:13]
    m_rs, m_tidal, m_bound = cols[13:16]

    d_r1 = distance(x_rs, x_c1)
    d_r2 = distance(x_rs, x_c2)
    d_12 = distance(x_c1, x_c2)

    ok_r1, ok_r2, ok_12 = d_r1 <= r_rs, d_r2 <= r_rs, d_12 <= r_rs
    
    print(len(idx))
    print("RS vs. 1:", idx[~ok_r1])
    print("RS vs. 2:", idx[~ok_r2])
    print(" 1 vs. 2", idx[~ok_12])

    agree = idx[d_r1 < r_rs]

    ok = ok_r2

    plt.figure()
    bins = np.linspace(0, 2, 31)
    
    plt.hist(m_tidal[ok]/m_rs[ok], bins=bins, density=True,
             label=r"$M_{\rm tidal}/M_{\rm rockstar}$",
             color=pc("r"), histtype="step", lw=2)
    plt.hist(m_bound[ok]/m_rs[ok], bins=bins, density=True,
             label=r"$M_{\rm bound}/M_{\rm rockstar}$",
             color=pc("o"), histtype="step", lw=2)
    plt.hist(m_bound[ok]/m_tidal[ok], bins=bins, density=True,
             label=r"$M_{\rm bound}/M_{\rm tidal}$",
             color=pc("b"), histtype="step", lw=2)

    ylo, yhi = plt.ylim()
    plt.legend(loc="upper left")
    plt.plot([1, 1], [ylo, yhi], "--", c="k", lw=2)
    plt.ylim(ylo, yhi)

    plt.figure()
    bins = np.linspace(0, 5, 21)
    
    plt.hist(r_tidal[ok]/r_rs[ok], bins=bins, density=True,
             label=r"$R_{\rm tidal}/R_{\rm rockstar}$",
             color=pc("r"), histtype="step", lw=2)
    plt.hist(r_bound[ok]/r_rs[ok], bins=bins, density=True,
             label=r"$R_{\rm bound,99}/R_{\rm rockstar}$",
             color=pc("o"), histtype="step", lw=2)
    plt.hist(r_bound[ok]/r_tidal[ok], bins=bins, density=True,
             label=r"$R_{\rm bound,99}/R_{\rm tidal}$",
             color=pc("b"), histtype="step", lw=2)

    ylo, yhi = plt.ylim()
    plt.legend(loc="upper right")
    plt.plot([1, 1], [ylo, yhi], "--", c="k", lw=2)
    plt.ylim(ylo, yhi)
    
    plt.show()
    
if __name__ == "__main__": main()
