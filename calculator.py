def get_potential(numpy, surface_dim, num_nodes_in_each_dim, z_ion_pos, ion_charge, wall_charge, f, eps):

    total_num_nodes = num_nodes_in_each_dim ** 2

    charge_per_node = wall_charge / total_num_nodes

    num_intervals = num_nodes_in_each_dim - 1
    
    lattice_param = surface_dim / num_intervals
    
    node_coord_1d = numpy.arange(num_nodes_in_each_dim) * lattice_param
    
    node_coord_1d -= surface_dim / 2
    
    ion_coord = numpy.array([0.0, 0.0, z_ion_pos])
    
    curr_node_coord = numpy.zeros(3, dtype = numpy.float)
    
    potential = 0.0
    
    charge = ion_charge * charge_per_node * f / eps
    
    for i in range(num_nodes_in_each_dim):
        for j in range(num_nodes_in_each_dim):
            curr_node_coord[0] = node_coord_1d[i]
            curr_node_coord[1] = node_coord_1d[j]
            potential += charge/numpy.linalg.norm(curr_node_coord - ion_coord)
    
    return potential


def phi(numpy, z, L):

    tmp = abs(z) / L

    tmp_s = tmp ** 2

    r1 = numpy.sqrt(0.5 + tmp_s)

    r2 = numpy.sqrt(0.25 + tmp_s)

    res = 0.5 * numpy.pi - numpy.arctan(4 * r1 * tmp)

    res = numpy.log((0.5 + r1)/r2) - tmp * res

    return 4 * L * res 