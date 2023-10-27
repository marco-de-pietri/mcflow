"""
code to write the CDGS for cylindrical pipes - author j. Alguacil UNED
"""
import sys
from .endf_extract_data import ENDF_library

# Global Variables
endf = ENDF_library()

def make_pipes_cdgs(pipes,dec_lambda, dec_yield,e_bins,e_probs, full_path):
    """
    function to write the CDGS for cylindrical pipes
    """

    # Open File
    with open(full_path,"w",encoding="utf8") as out:
    # Get Isotope Features

        # Data of the header
        nmesh,act_tot = get_header(pipes)
        source_factor = act_tot * dec_yield
        line = f"num_meshes {nmesh:d} \nglobal_source {source_factor:9.7e} \n"
        out.write(line)


        # Write Each pipe CDGS
        nmesh = 0
        for pipe in pipes.values():
            if pipe.node_type == 'tank':
                print ("WARNING: tank node not included in the CDGS")
                continue

            act = pipe.tot_activity_bq
            y_pipe = act*dec_yield

            #Introduction of the common decay gamma source
            nmesh = nmesh + 1
            line  = f"mesh_id {nmesh:d}\n"
            line  = line +f"Cylinder Cell {pipe.node_id:d}\n"
            line  = line +"Cooling_time  0.0\n"
            out.write(line)

            #General Photon Info
            line = write_general_photon_info(y_pipe,e_bins)
            out.write(line)

            # Mesh Structure
            line = write_mesh_structure(pipe)
            out.write(line)

            # Source Data
            line = write_source_data(pipe,dec_yield,e_probs)
            out.write(line)

    return 0

def write_source_data(pipe,yld,e_probs):
    """
    function description
    """
    line = "source_data\n"


    # index of the mesh, i_tot (ph/s), vol(cm**3), ncel
    vol_cyl = pipe.volume_cm3

    dec_pipe       = pipe.tot_activity_bq*yld
    stat_error = pipe.mc_error

    line = line + f"1  {dec_pipe:12.8e}  {vol_cyl:12.8e} 1\n"

    # Cell activity             [celda  relative vol  ac]
    line = line + f"{pipe.mcnp_cell} 1.00000000e+00  {dec_pipe:12.8e} \n"
    gammas = []
    error  = []
    for prb in e_probs:
        gammas.append(prb*dec_pipe) # prb is already normalized to yield
        if prb == 0.0:
            error.append(0.0)
        else:
            error.append(stat_error)

    line = line + write_vec_e(gammas,6) + "\n"
    line = line + write_vec_e(error,6)  + "\n"
    line = line + "end_source_data\n"

    return line

def write_mesh_structure(pipe):
    """
    write the cylinder mesh of the pipe
    """
    radius = pipe.radius_cm

    length = pipe.length_cm

    origin = [
                pipe.origin_x_cm,
                pipe.origin_y_cm,
                pipe.origin_z_cm,
             ]

    axs =    [
                pipe.axis_x,
                pipe.axis_y,
                pipe.axis_z,
             ]
    vec_2 =  [
                pipe.axis_x_1,
                pipe.axis_y_1,
                pipe.axis_z_1,
             ]


    # Definition of the boundaries
    theta_vals = [0.0,1.0]
    z_vals     = [0.0,length]
    r_vals     = [0.0,radius]


    # Assumed cyl [r, theta,z]
    line = "mesh_type cyl\n"
    line = line + "mesh_boundaries  2 2 2\n"

    # Origin and transformations
    line = line + f"{origin[0]:12.6e}  {origin[1]:12.6e}  {origin[2]:12.6e}\n"
    line = line + f"{axs[0]:12.6e}  {axs[1]:12.6e}  {axs[2]:12.6e}\n"
    line = line + f"{vec_2[0]:12.6e}  {vec_2[1]:12.6e}  {vec_2[2]:12.6e}\n"

    # Mesh boundaries r,theta,z

    line = line + write_vec_f(r_vals,6) + "\n"
    line = line + write_vec_f(theta_vals,6) + "\n"
    line = line + write_vec_f(z_vals,6)     + "\n"

    return line

def write_general_photon_info(p_yield,e_bins):
    """
    this function writest the energy spectrum of the pipes
    """
    line = f"total_source {p_yield:12.9e} \n"
    line = line + "energy_type bins\n"    # lines are not yet supported
    line = line + f"energy_boundaries {len(e_bins):d}\n"
    line = line + write_vec_ff(e_bins,6) +"\n"
    return line

def get_bins_from_lines(energy,prob,eps=1e-5):
    """
    convert the enegy lines to bins
    """
    eb_vals = [0.0]
    for ener in energy:
        dl_val = ener*eps
        eb_vals.append(ener-dl_val)
        eb_vals.append(ener+dl_val)
    pb_vals = []
    for p_val in prob:
        pb_vals.append(0.0)
        pb_vals.append(p_val)
    return eb_vals,pb_vals



def isotope_features(isotope,kind="photon",delete=1e-6):
    """
    get the decay constant and the energy and probability of the decay
    """
    # Landa
    landa = endf.get_lambda(isotope)

    energy = []
    probs   = []
    yld    = 0.0
    if kind == "photon":
        # gamma_spectra
        spectra = endf.gamma_spectra(isotope)
    elif kind == "neutron":
        # spectra = endf.neutron_spectra(isotope) # no esta hecho
        # neutron spcetra
        if isotope.capitalize() != "N17":
            sys.exit("Only neutrons of N17 are avialable")
        en_vals = [0.3828,0.884,1.1709,1.700300]
        y_vals  = [0.348066,0.005706,0.527805,0.070374]
        ex = [0,0,0,0]
        spectra = []
        # In eV
        for i in range(4):
            spectra.append( [en_vals[i]*1e6,y_vals[i],ex[i]])

    else:
        sys.exit( f"Bad particle type {kind:s}")
    for (en_vals,y_vals,_) in spectra:
        # Filter with also Pat uses
        if (en_vals*1e-6*y_vals) < delete:
            continue
        energy.append(en_vals*1e-6)
        probs.append(y_vals)
    yld = sum(probs)
    if yld == 0.0:
        print( f"Fatal Error: No gamma from the isotope {isotope:s}")
        sys.exit()
    return landa, energy, probs, yld

def get_header(nodes):
    """
    calculate the cdgs header factor
    """
    # Get Source Factor
    nmesh   = 0
    act_tot = 0.0
    for node in nodes.values():
        if node.node_type == 'tank':
            continue
        act = node.tot_activity_bq
        act_tot = act_tot + act
        nmesh = nmesh + 1
    return nmesh,act_tot

def write_vec_ff(vec,npos):
    """
    line formatting for the energy bins
    """
    j = -1
    line = ''
    for el_val in vec:
        j = j + 1
        if j == npos:
            line = line+'\n'
            j = 0
        line = line + f"{el_val:12.7e} "
    return line

def write_vec_f(vec,npos):
    """
    function description
    """
    j = -1
    line = ''
    for el_val in vec:
        j = j + 1
        if j == npos:
            line = line+'\n'
            j = 0
        line = line + f"{el_val:12.4f} "
    return line

def write_vec_e(vec,npos):
    """
    function description
    """
    j = -1
    line = ''
    for el_val in vec:
        j = j + 1
        if j == npos:
            line = line+'\n'
            j = 0
        line = line + f"{el_val:12.4e} "
    return line
