import numpy as np
import sys

def main():
    targets =["trees"]
    out_file = "../download_tables.py"

    suites = [
        "SymphonyLMC",
        "SymphonyMilkyWay",
        "SymphonyGroup",
        "SymphonyLCluster",
        "SymphonyCLuster"
    ]

    lines = []

    f = open(out_file, "w+")
    
    for target in targets:
        print("%s = {" % target, file=f)
        for suite in suites:
            print("    '%s': {" % suite, file=f)

            table_name = "%s/%s_links.csv" % (target, suite)
            halo_names, gd_ids, urls = np.loadtxt(
                table_name, delimiter=",", dtype=str, skiprows=1).T

            #for halo_name, gd_id in zip(halo_names, gd_ids):
            #    url = "https://drive.google.com/uc?id=%s" % gd_id
            for halo_name, url in zip(halo_names, urls):
                if halo_name == halo_names[-1]:
                    print("        '%s': '%s'" % (halo_name, url), file=f)
                else:
                    print("        '%s': '%s'," % (halo_name, url), file=f)

            if suite == suites[-1]:
                print("    }", file=f)
            else:
                print("    },", file=f)
        print("}\n", file=f)
        

if __name__ == "__main__": main()
