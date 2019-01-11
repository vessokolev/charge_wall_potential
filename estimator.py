#!/usr/bin/env python3

# Use Python 3.6 or later

import numpy
import calculator

surface_dim = 20.0

max_num_nodes_in_each_dim = 400

ion_charge = 0.12 # in number of electrones

wall_charge = -0.5 # in number of electrones

z_ion_pos = 0.275 # in nm

f = 138935.45838031892 # in J.nm/mol

eps = 78.6 # dimensionless, relative to the Vacuum permittivity

R = 8.3144598 # in J/mol/K
    
T = 298.15 # in K

surface_area =  surface_dim ** 2


pot_theoretical = ion_charge * wall_charge * f / eps / surface_dim / surface_dim

pot_theoretical *= calculator.phi(numpy, z_ion_pos, surface_dim)


f_obj = open('potential.txt', 'w')


for i in range(2, max_num_nodes_in_each_dim):

    p = calculator.get_potential(numpy, surface_dim, i, z_ion_pos, ion_charge, wall_charge, f, eps)

    p_diff = p - pot_theoretical

    p_diff_perc = abs(p_diff / pot_theoretical) * 100

    d = (i / surface_dim) ** 2

    str_ = '{:10d}'.format(i) + ' '
    
    str_ += '{:14.5e}'.format(p) + ' '
    
    str_ += '{:14.5e}'.format(p_diff) + ' '

    str_ += '{:14.5e}'.format(p_diff_perc) + ' '

    str_ += '{:14.5e}'.format(d) + '\n'

    f_obj.write(str_)

    print(i, p, p_diff, p_diff_perc, d)
    
f_obj.close()


pot_theoretical = ion_charge * wall_charge * f / eps / surface_dim / surface_dim

pot_theoretical *= calculator.phi(numpy, z_ion_pos, surface_dim)


print('The theoretical value of the potential: ', pot_theoretical)
