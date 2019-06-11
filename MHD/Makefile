ifndef config
config = magdif.in
endif
magdifplot := ../RUN/PLOTTING/magdifplot.py
mesh := ../PRELOAD/inputformaxwell.msh

extract_proj = sed -e 's/ \+/\t/g; s/^\t//g' $(1) | cut -f 9-10 > $(2)
poloidal_modes = cd ../RUN && ./result_spectrum.x $(1) > $(2) && cd ../MHD

.PHONY: plots

plots: Bmn_psi.dat Bmn_vac_psi.dat
	python $(magdifplot) . $(config) $(mesh)

Bmn_psi.dat: plot_Bn.dat
	$(call extract_proj, $<, ../RUN/Bn_psi.dat)
	$(call poloidal_modes, Bn_psi.dat, $@)
	mv ../RUN/$@ ./

Bmn_vac_psi.dat: plot_Bn_vac.dat
	$(call extract_proj, $<, ../RUN/Bn_vac_psi.dat)
	$(call poloidal_modes, Bn_vac_psi.dat, $@)
	mv ../RUN/$@ ./
