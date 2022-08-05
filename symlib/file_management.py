from . import util
from os import path
import os

def pack_files(suite, halo_name, base_out_dir, target="halos",
               logging=True):
    if logging:
        print("Packing halo %s in suite %s to %s" % (
            str(halo_name), suite, base_out_dir))
    
    if suite is None:
        for suite in ["SymphonyLMC", "SymphonyMilkyWay", "SymphonyGroup",
                      "SymphonyLCluster", "SymphonyCluster"]:
            tar_files(suite, halo_name, base_out_dir, target=target)

    if halo_name is None:
        for i in range(util.n_hosts(suite)):
            tar_files(suite, i, base_out_dir, target=target)
        return

    # MUST BE RUN FROM base_dir!!!
    if isinstance(halo_name, int):
        halo_name = util.DEFAULT_HALO_NAMES[suite][halo_name]

    rel_dir = path.join(suite, halo_name)
    rel_out = path.join(base_out_dir, suite, "halos")
    os.system("mkdir -p %s" % rel_out)
    os.system("tar --exclude='%s/*core*' --exclude='%s/particles' -cf %s/%s.tar.gz %s" % (rel_dir, rel_dir, rel_out, halo_name, rel_dir))

def download_files(suite, halo_name, base_out_dir, target="halos"):
    if logging:
        print("Downloading halo %s in suite %s to %s" % (
            str(halo_name), suite, base_out_dir))
    
    if suite is None:
        for suite in ["SymphonyLMC", "SymphonyMilkyWay", "SymphonyGroup",
                      "SymphonyLCluster", "SymphonyCluster"]:
            download_files(suite, halo_name, base_out_dir, target=target)

    if halo_name is None:
        for i in range(util.n_hosts(suite)):
            download_files(suite, i, base_out_dir, target=target)
        return

    # MUST BE RUN FROM base_dir!!!
    if isinstance(halo_name, int):
        halo_name = util.DEFAULT_HALO_NAMES[suite][halo_name]
