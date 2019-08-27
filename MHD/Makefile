ifndef config
config = magdif.in
endif
magdifplot := ../RUN/PLOTTING/magdifplot.py
mesh := ../PRELOAD/inputformaxwell.msh

Bn_diffs := $(sort $(wildcard Bn_diff_*.dat))

extract_psi = sed -e 's/ \+/\t/g; s/^\t//g' $(1) | cut -f 9-10 > $(2)
extract_theta = sed -e 's/ \+/\t/g; s/^\t//g' $(1) | cut -f 11-12 > $(2)
poloidal_modes = cd ../RUN && ./result_spectrum.x $(1) > $(2) && cd ../MHD

.PHONY: run plots spectra clean
.ONESHELL: convergence.dat

run: $(config) $(mesh)
	./magdif_test.sh $(config)

plots:  convergence.dat \
	$(magdifplot) $(config) $(mesh)
	python $(magdifplot) . $(config) $(mesh)

spectra: Bmn_psi.dat Bmn_vac_psi.dat currmn_001_theta.dat

clean:
	rm -f *.pdf *.out Bn*.dat presn*.dat currn*.dat eigvec*.dat plot*.dat j0phi.dat fluxvar.dat Bmn_psi.dat Bmn_vac_psi.dat currmn_001_theta.dat convergence.dat fort.333

Bmn_psi.dat: plot_Bn.dat
	$(call extract_psi, $<, Bn_psi.dat)
	$(call poloidal_modes, ../MHD/Bn_psi.dat, ../MHD/$@)
	rm Bn_psi.dat

Bmn_vac_psi.dat: plot_Bn_vac.dat
	$(call extract_psi, $<, Bn_vac_psi.dat)
	$(call poloidal_modes, ../MHD/Bn_vac_psi.dat, ../MHD/$@)
	rm Bn_vac_psi.dat

currmn_001_theta.dat: plot_currn_001.dat
	$(call extract_theta, $<, currn_001_theta.dat)
	$(call poloidal_modes, ../MHD/currn_001_theta.dat, ../MHD/$@)
	rm currn_001_theta.dat

convergence.dat: $(Bn_diffs)
	rm -f $@
	cd ../FEM
	FreeFem++-nw L2int.edp ../MHD/Bn_vac.dat ../MHD/$@ 1>/tmp/freefem.out 2>&1
	for file in $^
	do
	FreeFem++-nw L2int.edp ../MHD/$$file ../MHD/$@ 1>/tmp/freefem.out 2>&1
	done
	cd ../MHD
