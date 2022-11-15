from . import util
from . import lib
from os import path
import os
import shutil
import requests
import time

DOWNLOAD_TARGETS = ["halos", "trees"]
BASE_URL_FORMAT = "https://%s:%s@s3df.slac.stanford.edu/data/kipac/symphony/"
RETRY_WAIT_TIME = 5

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


        
def download_packed_files(user, password, suite, halo_name, base_out_dir,
                          target="halos", logging=True, retries=5):
    orig_halo_name = halo_name
    if logging:
        print("Downloading halo %s in suite %s to %s" % (
            str(halo_name), suite, base_out_dir))
    
    if suite is None:
        for suite in ["SymphonyLMC", "SymphonyMilkyWay", "SymphonyGroup",
                      "SymphonyLCluster", "SymphonyCluster"]:
            download_packed_files(user, password, suite, halo_name, base_out_dir,
                                  target=target, logging=logging, retries=retries)
        return
        
    if halo_name is None:
        for i in range(util.n_hosts(suite)):
            download_packed_files(user, password, suite, i, base_out_dir,
                                  target=target, logging=logging, retries=retries)
        return

    if target not in DOWNLOAD_TARGETS:
        raise ValueError("Unrecoginzed download target, %s" % target)
    
    if suite not in lib.parameter_table:
        raise ValueError("Unrecognized simulation suite, %s" % suite)
    
    if isinstance(halo_name, int):
        halo_name = util.DEFAULT_HALO_NAMES[suite][halo_name]

    if halo_name not in util.DEFAULT_HALO_NAMES[suite]:
        raise ValueError("Unrecognized halo, %s, for the suite %s" %
                         (str(halo_name), suite))
                
    out_dir = path.join(base_out_dir, "tar_files", suite, target)
    os.makedirs(out_dir, exist_ok=True)

    base_url = BASE_URL_FORMAT % (user, password)
    url = "%s%s/%s/%s.tar.gz" % (base_url, suite, target, halo_name)
    out_file = path.join(out_dir, "%s.tar" % halo_name)

    # The next three lines are adapted from
    # https://stackoverflow.com/a/39217788
    # Original author is John Zwink, edited by the users vog, Martjin Peters, and
    # Daniel F.
    with requests.get(url, stream=True) as r:
        with open(out_file, 'wb') as f:
            shutil.copyfileobj(r.raw, f)
            
    if not path.exists(out_file):
        raise FileNotFoundError("Unable to download halo %s from the suite %s" %
                                (halo_name, suite))
    
    with open(out_file, "rb") as f:
        text = f.read(6)
        if text == b"<html>":
            f.seek(0, 0)
            text = f.read()
            os.remove(out_file)

            if text == b"<html>\r\n<head><title>401 Authorization Required</title></head>\r\n<body>\r\n<center><h1>401 Authorization Required</h1></center>\r\n<hr><center>nginx</center>\r\n</body>\r\n</html>\r\n":
                raise ValueError("Unrecognized user and password combination.")
            elif text == b"<html>\r\n<head><title>502 Bad Gateway</title></head>\r\n<body>\r\n<center><h1>502 Bad Gateway</h1></center>\r\n<hr><center>nginx</center>\r\n</body>\r\n</html>\r\n":
                if retries == 0:
                    raise ValueError("Recieved network error after exhausing maximum number of retries for halo %s. You could try increasing the retries variable, but the server is probably overloaded right now. Please try again later" % orig_halo_name)
                else:
                    if logging:
                        print("Internal server error ('502 Bad Gateway') encountered while downloading halo %s: trying to download again after a short wait (tries left: %d)" % (orig_halo_name, retries))
                    time.sleep(RETRY_WAIT_TIME)
                    download_packed_files(
                        user, password, suite, halo_name, base_out_dir,
                        target=target, logging=logging, retries=retries-1)
            else:
                raise ValueError("Server sent back the following error on halo %s:\n%s" %
                                 (orig_halo_name, text))
        
            
def unpack_files(suite, halo_name, base_out_dir,
                 target="halos", logging=True):
    if logging:
        print("Unpacking halo %s in suite %s to %s" % (
            str(halo_name), suite, base_out_dir))
    
    if suite is None:
        for suite in ["SymphonyLMC", "SymphonyMilkyWay", "SymphonyGroup",
                      "SymphonyLCluster", "SymphonyCluster"]:
            unpack_files(suite, halo_name, base_out_dir, target=target)
        return
            
    if halo_name is None:
        for i in range(util.n_hosts(suite)):
            unpack_files(suite, i, base_out_dir, target=target)
        return

    if isinstance(halo_name, int):
        halo_name = util.DEFAULT_HALO_NAMES[suite][halo_name]

    tar_loc = path.join(base_out_dir, "tar_files", suite,
                        target, "%s.tar" % halo_name)
    if not path.exists(tar_loc):
        raise FileNotFoundError("The file associated with the requested target data for host %s in suite %s does not exist within the base directory %s. This is probably because it was already downloaded or already unpacked." % (halo_name, suite, base_out_dir))
    shutil.unpack_archive(tar_loc, base_out_dir)
    os.remove(tar_loc)
    
    
def download_files(user, password, suite, halo_name, base_out_dir,
                   target="halos", logging=True, retries=5):
    download_packed_files(user, password, suite, halo_name,
                          base_out_dir, target, logging, retries=retries)
    unpack_files(suite, halo_name, base_out_dir, target, logging)
