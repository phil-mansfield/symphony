import numpy as np
import lib
import matplotlib.pyplot as plt
import os.path as path

def main():
    suite_name, halo_name = "SymphonyMilkyWay", "Halo023"
    sim_dir = path.join("../../symphony_pipeline/tmp_data", suite_name, halo_name)

    param = lib.parameter_table[suite_name]
    
    h, hist = lib.read_subhalos(param, sim_dir)

    print(h["x"][20,-10:])
    print(h["x_core"][20,-10:])
    
if __name__ == "__main__": main()
