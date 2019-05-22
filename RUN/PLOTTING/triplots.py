#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 09:15:14 2016

@author: Christopher Albert
"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import sys
try:
	from StringIO import StringIO
except ImportError:
	from io import StringIO

import glob
import os
import re
from multiprocessing import Pool

class TriPlot():
	@classmethod
	def load_mesh(cls, filename):
		with open(filename, 'r') as f:
			data = f.readlines()
		print('meshfile header: ', data[0])
		[NN, NT, NE] = np.loadtxt(StringIO(data[0]), dtype = int)
		node = np.loadtxt(StringIO(''.join(data[1:NN+1])), dtype = float)
		nlab = np.array(node[:, 2], dtype = int)
		cls.node = node[:, 0:2]
		tri = np.loadtxt(StringIO(''.join(data[NN+1:NN+NT+1])), dtype = int)
		tlab = tri[:, 3]
		cls.tri = tri[:, 0:3]
		edge = np.loadtxt(StringIO(''.join(data[NN+NT+1:])), dtype = int)
		elab = edge[:, 2]
		edge = edge[:, 0:2]

	def __init__(self, data, limits, filename):
		self.data = data
		self.limits = limits
		self.filename = filename

	def plot_and_dump(self):
		"""plots values with tripcolor and colorbar, then dumps them in PDF file"""
		plt.figure()
		plt.tripcolor(
			self.__class__.node[:, 0],
			self.__class__.node[:, 1],
			self.__class__.tri - 1,
			self.data,
			cmap = 'jet'
		)
		plt.axis('equal')
		plt.colorbar()
		plt.clim(self.limits)
		plt.savefig(self.filename, format = 'pdf')
		plt.close()

meshfile = sys.argv[1]
datadir = sys.argv[2]

TriPlot.load_mesh(meshfile)
plots = []

vector_input = glob.glob(os.path.join(datadir, 'plot_*.dat'))
scalar_input = glob.glob(os.path.join(datadir, 'presn*.dat'))
equil_input = glob.glob(os.path.join(datadir, 'fluxvar.dat'))
vector_re = re.compile(r"(?<=plot_)([^/]+)\.dat$")
scalar_re = re.compile(r"presn([^/]*)\.dat$")
equil_re = re.compile(r"fluxvar\.dat$")
# vector_infix = ['ReR', 'ImR', 'ReZ', 'ImZ', 'Rephi', 'Imphi']
vector_infix = ['Reproj', 'Improj']
scalar_infix = ['Re', 'Im']
equil_infix = ['psi', 'q', 'dens', 'temp', 'pres0']

for datafile in vector_input:
	print('reading contents of ', datafile)
	contents = np.genfromtxt(datafile)
	for column, infix in enumerate(vector_infix, 8):
		plots.append(TriPlot(
			data = contents[:, column],
			limits = [-max(abs(contents[:, column])),
					  max(abs(contents[:, column]))],
			filename = vector_re.sub(r"\1_" + infix + '.pdf', datafile)
		))
for datafile in scalar_input:
	print('reading contents of ', datafile)
	contents = np.genfromtxt(datafile)
	for column, infix in enumerate(scalar_infix):
		plots.append(TriPlot(
			data = contents[:, column],
			limits = [-max(abs(contents[:, column])),
					  max(abs(contents[:, column]))],
			filename = scalar_re.sub(r"plot_presn\1_" + infix + '.pdf', datafile)
		))
for datafile in equil_input:
	print('reading contents of ', datafile)
	contents = np.genfromtxt(datafile)
	for column, infix in enumerate(equil_infix):
		plots.append(TriPlot(
			data = contents[:, column],
			limits = [-max(abs(contents[:, column])) if infix == 'q' else 0.0,
					  0.0 if infix == 'q' else max(abs(contents[:, column]))],
			filename = equil_re.sub('plot_' + infix + '.pdf', datafile)
		))

def wrapper(index):
	plots[index].plot_and_dump()

if __name__ == '__main__':
	with Pool(max(1, os.cpu_count() - 1)) as p:
		p.map(wrapper, range(len(plots)))
