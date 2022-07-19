import matplotlib.pyplot as plt
import numpy as np
import os.path as path

def plot_circle(ax, x0, y0, r, **kwargs):
    """ plot_circle plots a circle with radius r centered on the position
    (x0, y0) onto the matplotlib axis, ax. This function accepts all the named
    arguments that matplotlib's plot function accepts.
    """
    th = np.linspace(0, 2*np.pi, 100)
    y = np.sin(th)*r + y0
    x = np.cos(th)*r + x0
    ax.plot(x, y, **kwargs)

def set_units_x(x, h, scale, param):
    """ set_units_x converts an array of positions with Gadget-2 code units
    to physical kpc, centered on the subhalo h (a single element of a
    SUBHALO_DTYPE array).
    """
    dx = np.zeros(x.shape)
    for dim in range(3): dx[:,dim] = x[:,dim] - h["x"][dim]
    dx *= scale*1e3/param["h100"]
    return dx

def set_units_v(v, h, scale, param):
    """ set_units_v converts an array of velocities with Gadget-2 code units
    to physical km/s, centered on the subhalo h (a single element of a
    SUBHALO_DTYPE array).
    """
    v *= np.sqrt(scale)
    dv = np.zeros(v.shape)
    for dim in range(3): dv[:,dim] = v[:,dim] - h["v"][dim]
    return dv

def set_units_parameters(scale, param):
    """ set_units_parameters returns particle mass, mp, and Plummer-equivalent
    force-softening scales in units of Msun and physical kpc, respectively.
    scale is the current scale factor and param is the parameter dictionary
    returned by simulation_parameters.
    """
    return param["mp"]/param["h100"], param["eps"]*scale/param["h100"]

def set_units_histories(hist, scale, param):
    """ set_units_histories takes an array of type HISTORY_DTYPE and returns a
    copy with units that have been converted from Rockstar's detail units to
    symlib's default units: Msun, physical kpc, and physical km/s.
    """
    hist = np.copy(hist)
    hist["mpeak"] /= param["h100"]
    return hist

def set_units_halos(h, scale, param):
    """ set_units_halos takes a array of type SUBHALO_DTYPE and returns a copy
    with units that have been converted from Rockstar's default units to
    symlib's default units, Msun, physical kpc centered on the host halo,
    and physical km/s centered on the host halo.
    """
    h = np.copy(h)    
    for dim in range(3):
        x0 = np.copy(h["x"][0,:,dim])
        v0 = np.copy(h["v"][0,:,dim])

        for hi in range(len(h)):
            h["x"][hi,:,dim] -= x0
            h["v"][hi,:,dim] -= v0
                
            h["x"][hi,:,dim] *= 1e3*scale/param["h100"]
            
            if dim == 0:
                h["rvir"][hi,:] *= 1e3*scale/param["h100"]
                h["mvir"][hi,:] *= 1/param["h100"]

    
    for hi in range(len(h)):
        invalid = h[hi]["rvir"] < 0
        h[hi, invalid]["rvir"] = -1
        h[hi, invalid]["x"] = -1
        h[hi, invalid]["v"] = -1
        
    return h

DEFAULT_HALO_NAMES = {
    "SymphonyLMC": sorted([
        "Halo032", "Halo097", "Halo218", "Halo374", "Halo463", "Halo567", "Halo721", "Halo853",
        "Halo059", "Halo104", "Halo296", "Halo380", "Halo4662", "Halo575", "Halo767", "Halo914",
        "Halo0662", "Halo110", "Halo301", "Halo391", "Halo511", "Halo602", "Halo802", "Halo932",
        "Halo083", "Halo202", "Halo303", "Halo405", "Halo524", "Halo697", "Halo824", "Halo933",
        "Halo088", "Halo208", "Halo340", "Halo440", "Halo539", "Halo711", "Halo850"
    ]),
    "SymphonyMilkyWay": sorted([
        "Halo023", "Halo268", "Halo364", "Halo440", "Halo558", "Halo641", "Halo797", "Halo878", "Halo939",
        "Halo088", "Halo270", "Halo374", "Halo460", "Halo567", "Halo675", "Halo800", "Halo881", "Halo967",
        "Halo119", "Halo288", "Halo414", "Halo469", "Halo570", "Halo718", "Halo825", "Halo925", "Halo9749",
        "Halo188", "Halo327", "Halo415", "Halo490", "Halo606", "Halo738", "Halo829", "Halo926", "Halo9829",
        "Halo247", "Halo349", "Halo416", "Halo530", "Halo628", "Halo749", "Halo852", "Halo937", "Halo990"
    ]),
    "SymphonyGroup": sorted([
        "Halo015", "Halo175", "Halo313", "Halo352", "Halo470", "Halo579", "Halo755", "Halo806", "Halo985",
        "Halo024", "Halo183", "Halo331", "Halo383", "Halo496", "Halo581", "Halo759", "Halo834",
        "Halo029", "Halo186", "Halo336", "Halo384", "Halo501", "Halo593", "Halo774", "Halo861",
        "Halo046", "Halo257", "Halo342", "Halo399", "Halo504", "Halo606", "Halo781", "Halo909",
        "Halo055", "Halo286", "Halo346", "Halo428", "Halo551", "Halo656", "Halo785", "Halo927",
        "Halo090", "Halo302", "Halo347", "Halo444", "Halo570", "Halo752", "Halo790", "Halo962",
    ]),
    "SymphonyLCluster": sorted([
        "Halo_000", "Halo_005", "Halo_012", "Halo_017", "Halo_023", "Halo_028", "Halo_046",
        "Halo_001", "Halo_006", "Halo_013", "Halo_019", "Halo_024", "Halo_029", "Halo_047",
        "Halo_002", "Halo_008", "Halo_014", "Halo_020", "Halo_025", "Halo_042", "Halo_050",
        "Halo_003", "Halo_009", "Halo_015", "Halo_021", "Halo_026", "Halo_043",
        "Halo_004", "Halo_010", "Halo_016", "Halo_022", "Halo_027", "Halo_044",
    ]),
    "SymphonyCluster": sorted([
        "Halo156", "Halo293", "Halo339", "Halo367", "Halo409", "Halo448", "Halo475", "Halo522", "Halo629",
        "Halo175", "Halo304", "Halo345", "Halo377", "Halo415", "Halo452", "Halo476", "Halo529", "Halo631",
        "Halo200", "Halo305", "Halo346", "Halo378", "Halo416", "Halo454", "Halo478", "Halo544", "Halo639",
        "Halo211", "Halo306", "Halo348", "Halo385", "Halo419", "Halo455", "Halo479", "Halo545", "Halo645",
        "Halo213", "Halo308", "Halo349", "Halo386", "Halo428", "Halo456", "Halo480", "Halo546", "Halo653",
        "Halo222", "Halo317", "Halo352", "Halo387", "Halo429", "Halo461", "Halo483", "Halo561", "Halo734",
        "Halo225", "Halo321", "Halo354", "Halo390", "Halo436", "Halo462", "Halo489", "Halo572",
        "Halo266", "Halo324", "Halo358", "Halo391", "Halo437", "Halo465", "Halo494", "Halo574",
        "Halo274", "Halo326", "Halo360", "Halo394", "Halo441", "Halo471", "Halo502", "Halo595",
        "Halo277", "Halo335", "Halo361", "Halo400", "Halo445", "Halo472", "Halo517", "Halo600",
        "Halo282", "Halo337", "Halo366", "Halo407", "Halo447", "Halo474", "Halo518", "Halo604"
    ]),
    "MWest": sorted([
        "Halo004", "Halo170", "Halo282", "Halo407", "Halo523", "Halo666", "Halo756", "Halo908", "Halo983",
        "Halo113", "Halo222", "Halo327", "Halo453", "Halo625", "Halo719", "Halo788", "Halo953",
        "Halo169", "Halo229", "Halo349", "Halo476", "Halo659", "Halo747", "Halo858", "Halo975"
    ]),
}

SUITE_NAMES = ["SymphonyLMC", "SymphonyMilkyWay", "SymphonyGroup", "SymphonyLCluster", "SymphonyCluster", "MWest"]

def n_hosts(suite_name):
    """ n_halos returns the number of host halos in a given suite,
    """
    if suite_name not in DEFAULT_HALO_NAMES.keys():
        raise ValueError("Unrecognized suite name, '%s'" % suite_name)
    
    return len(DEFAULT_HALO_NAMES[suite_name])

def get_host_directory(base_dir, suite_name, halo_name):
    """ get_host_directory returns the name of the directory containing a given
    halo. base_dir is the central directory, suite_name is the name of the
    suite, and halo_name is the name of the halo. halo_name can either be
    a string with the exact name of the halo or it can be an integer, indexing
    into the list of halo names in sorted order. When looping over all the
    halos in a suite, you'll want to use this second option.
    """
    
    if isinstance(halo_name, int):
        halo_name = DEFAULT_HALO_NAMES[suite_name][halo_name]

    return path.join(base_dir, suite_name, halo_name)
