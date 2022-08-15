from . import util
from . import lib
from . import download_tables
from os import path
import os
import gdown
import shutil

def pack_files(suite, halo_name, base_out_dir, target="halos",
               logging=True):
    if logging:
        print("Packing halo %s in suite %s to %s" % (
            str(halo_name), suite, base_out_dir))
    
    if suite is None:
        for suite in ["SymphonyLMC", "SymphonyMilkyWay", "SymphonyGroup",
                      "SymphonyLCluster", "SymphonyCluster"]:
            pack_files(suite, halo_name, base_out_dir, target=target)

    if halo_name is None:
        for i in range(util.n_hosts(suite)):
            pack_files(suite, i, base_out_dir, target=target)
        return

    # MUST BE RUN FROM base_dir!!!
    if isinstance(halo_name, int):
        halo_name = util.DEFAULT_HALO_NAMES[suite][halo_name]

    if target == "trees":
        rel_dir = path.join(suite, halo_name)
        rel_out = path.join(base_out_dir, suite, target)
        os.system("mkdir -p %s" % rel_out)
        os.system("tar --exclude='%s/halos/*core*' --exclude='%s/particles' -cf %s/%s.tar.gz %s" % (rel_dir, rel_dir, rel_out, halo_name, rel_dir))
    elif target == "halos":
        rel_dir = path.join(suite, halo_name)
        rel_out = path.join(base_out_dir, suite, target)
        os.system("mkdir -p %s" % rel_out)
        os.system("tar --exclude='%s/halos/*core*' --exclude='%s/halos/tree_[0-9].dat' --exclude='%s/particles' -cf %s/%s.tar.gz %s" % (rel_dir, rel_dir, rel_dir, rel_out, halo_name, rel_dir))

def download_packed_files(suite, halo_name, base_out_dir,
                          target="halos", logging=True):
    if logging:
        print("Downloading halo %s in suite %s to %s" % (
            str(halo_name), suite, base_out_dir))
    
    if suite is None:
        for suite in ["SymphonyLMC", "SymphonyMilkyWay", "SymphonyGroup",
                      "SymphonyLCluster", "SymphonyCluster"]:
            download_packed_files(suite, halo_name, base_out_dir,
                                  target=target)

    if halo_name is None:
        for i in range(util.n_hosts(suite)):
            download_packed_files(suite, i, base_out_dir, target=target)
        return

    if isinstance(halo_name, int):
        halo_name = util.DEFAULT_HALO_NAMES[suite][halo_name]

    if suite not in lib.parameter_table:
        raise ValueError("Unrecognized simulation suite, %s" % suite)
        
    out_dir = path.join(base_out_dir, "tar_files", suite, target)
    os.makedirs(out_dir, exist_ok=True)

    if target == "trees":
        table = download_tables.trees
    else:
        raise ValueError("Unknown download target, '%s'" % target)

    gdown.download(
        table[suite][halo_name],
        path.join(out_dir, "%s.tar" % halo_name),
        quiet=not logging,
        fuzzy=True
    )

    if not path.exists(path.join(out_dir, "%s.tar" % halo_name)):
        raise FileNotFoundError("unable to download halo %s from the suite %s" %
                                (halo_name, suite))
        
def unpack_files(suite, halo_name, base_out_dir,
                 target="halos", logging=True):
    if logging:
        print("Unpacking halo %s in suite %s to %s" % (
            str(halo_name), suite, base_out_dir))
    
    if suite is None:
        for suite in ["SymphonyLMC", "SymphonyMilkyWay", "SymphonyGroup",
                      "SymphonyLCluster", "SymphonyCluster"]:
            download_files(suite, halo_name, base_out_dir, target=target)

    if halo_name is None:
        for i in range(util.n_hosts(suite)):
            download_files(suite, i, base_out_dir, target=target)
        return

    if isinstance(halo_name, int):
        halo_name = util.DEFAULT_HALO_NAMES[suite][halo_name]

    tar_loc = path.join(base_out_dir, "tar_files", suite,
                        target, "%s.tar" % halo_name)
    shutil.unpack_archive(tar_loc, base_out_dir)
    os.remove(tar_loc)
    
    
def download_files(suite, halo_name, base_out_dir,
                   target="halos", logging=True):
    download_packed_files(suite, halo_name, base_out_dir, target, logging)
    unpack_files(suite, halo_name, base_out_dir, target, logging)
