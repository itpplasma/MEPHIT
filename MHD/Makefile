triplot := ../RUN/PLOTTING/triplot.py
mesh := ../PRELOAD/inputformaxwell.msh
plot = python $(triplot) $(mesh) $(1) $(2) $(3) 1>/dev/null

# VECTOR_PLOTS := $(shell ls plot_*.dat | sed -e 's/\([^ ]\+\)\.dat/\1_ReR.pdf \1_ImR.pdf \1_ReZ.pdf \1_ImZ.pdf \1_Rephi.pdf \1_Imphi.pdf/g')
VECTOR_PLOTS := $(shell ls plot_*.dat | sed -e 's/\([^ ]\+\)\.dat/\1_Reproj.pdf \1_Improj.pdf/g')
SCALAR_PLOTS := $(shell ls presn*.dat | sed -e 's/\([^ ]\+\)\.dat/plot_\1_Re.pdf plot_\1_Im.pdf/g')
EQUIL_PLOTS := $(addprefix plot_, psi.pdf q.pdf dens.pdf temp.pdf pres0.pdf)

.PHONY: plots

plots: $(EQUIL_PLOTS) $(SCALAR_PLOTS) $(VECTOR_PLOTS)

plot_%_ReR.pdf: plot_%.dat
	$(call plot, $<, 2, $@)
plot_%_ImR.pdf: plot_%.dat
	$(call plot, $<, 3, $@)
plot_%_ReZ.pdf: plot_%.dat
	$(call plot, $<, 4, $@)
plot_%_ImZ.pdf: plot_%.dat
	$(call plot, $<, 5, $@)
plot_%_Rephi.pdf: plot_%.dat
	$(call plot, $<, 6, $@)
plot_%_Imphi.pdf: plot_%.dat
	$(call plot, $<, 7, $@)
plot_%_Reproj.pdf: plot_%.dat
	$(call plot, $<, 8, $@)
plot_%_Improj.pdf: plot_%.dat
	$(call plot, $<, 9, $@)

plot_pres%_Re.pdf: pres%.dat
	$(call plot, $<, 0, $@)
plot_pres%_Im.pdf: pres%.dat
	$(call plot, $<, 1, $@)

plot_psi.pdf: fluxvar.dat
	$(call plot, $<, 0, $@)
plot_q.pdf: fluxvar.dat
	$(call plot, $<, 1, $@)
plot_dens.pdf: fluxvar.dat
	$(call plot, $<, 2, $@)
plot_temp.pdf: fluxvar.dat
	$(call plot, $<, 3, $@)
plot_pres0.pdf: fluxvar.dat
	$(call plot, $<, 4, $@)
