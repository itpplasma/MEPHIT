#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 16:29:35 2019

@author: patrick
"""

import sys
import os.path
import f90nml.parser

p = f90nml.parser.Parser()
nml = p.read(sys.argv[1])
config = nml['settings']

KiLCA_scale_factor = config['kilca_scale_factor']
if (KiLCA_scale_factor == 0):
    nmode = config['n']
else:
    nmode = config['n'] * KiLCA_scale_factor
ktlownfluxp1 = (2 * config['nflux'] - 1) * config['nkpol']
if (os.path.isabs(config['currn_file'])):
    currnfile = config['currn_file']
else:
    currnfile = os.path.relpath(config['currn_file'], '../FEM')
if (os.path.isabs(config['Bn_file'])):
    Bnfluxfile = config['Bn_file']
else:
    Bnfluxfile = os.path.relpath(config['Bn_file'], '../FEM')

with open('../FEM/magdif.idp', 'w') as f:
    f.write('real nmode = {};\n'.format(nmode))
    f.write('int ktlownfluxp1 = {};\n'.format(ktlownfluxp1))
    f.write('string meshfile = "../PRELOAD/inputformaxwell_ext.msh";\n')
    f.write('string currnfile = "{}";\n'.format(currnfile))
    f.write('string Bnfluxfile = "{}";\n'.format(Bnfluxfile))
    f.write('bool doplot = false;\n')
