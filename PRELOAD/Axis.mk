#FC       = lf95 --ap
FC       = gfortran

#OPTS= --chk a,e,s,u,x --trace --trap -g -M OBJS
#OPTS= -O -M OBJS
OPTS= -O -J OBJS

OBJS =  OBJS/field_divB0.o \
	OBJS/bdivfree.o \
	OBJS/from_nrtype.o \
	OBJS/brent.o \
	OBJS/rtbis.o \
	OBJS/odeint.o \
	OBJS/rkck.o \
	OBJS/rkqs.o \
	OBJS/spline5_RZ.o \
	OBJS/leastsq.o \
	OBJS/lubksb.o \
	OBJS/ludcmp.o \
	OBJS/polylag.o \
	OBJS/twostreights.o \
	OBJS/mods4ax.o \
	OBJS/in_out_cw.o \
	OBJS/inside_of_wall.o \
	OBJS/axis.o 


axis.x: $(OBJS) Axis.mk
	$(FC) $(OPTS) -o axis.x $(OBJS)
OBJS/field_divB0.o: SRC/field_divB0.f90 Axis.mk
	$(FC) $(OPTS) -c SRC/field_divB0.f90
	mv field_divB0.o OBJS/
OBJS/bdivfree.o: SRC/bdivfree.f90 Axis.mk
	$(FC) $(OPTS) -c SRC/bdivfree.f90
	mv bdivfree.o OBJS/
OBJS/from_nrtype.o: SRC/from_nrtype.f90 Axis.mk
	$(FC) $(OPTS) -c SRC/from_nrtype.f90
	mv from_nrtype.o OBJS/
OBJS/spline5_RZ.o: SRC/spline5_RZ.f90 Axis.mk
	$(FC) $(OPTS) -c SRC/spline5_RZ.f90
	mv spline5_RZ.o OBJS/
OBJS/odeint.o: SRC/odeint.f Axis.mk
	$(FC) $(OPTS) -c SRC/odeint.f
	mv odeint.o OBJS/
OBJS/rkck.o: SRC/rkck.f Axis.mk
	$(FC) $(OPTS) -c SRC/rkck.f
	mv rkck.o OBJS/
OBJS/rkqs.o: SRC/rkqs.f Axis.mk
	$(FC) $(OPTS) -c SRC/rkqs.f
	mv rkqs.o OBJS/
OBJS/leastsq.o: SRC/leastsq.f90 Axis.mk
	$(FC) $(OPTS) -c SRC/leastsq.f90
	mv leastsq.o OBJS/
OBJS/lubksb.o: SRC/lubksb.f Axis.mk
	$(FC) $(OPTS) -c SRC/lubksb.f
	mv lubksb.o OBJS/
OBJS/ludcmp.o: SRC/ludcmp.f Axis.mk
	$(FC) $(OPTS) -c SRC/ludcmp.f
	mv ludcmp.o OBJS/
OBJS/polylag.o: SRC/polylag.f Axis.mk
	$(FC) $(OPTS) -c SRC/polylag.f
	mv polylag.o OBJS/
OBJS/brent.o: SRC/brent.f90 Axis.mk
	$(FC) $(OPTS) -c SRC/brent.f90
	mv brent.o OBJS/
OBJS/rtbis.o: SRC/rtbis.f90 Axis.mk
	$(FC) $(OPTS) -c SRC/rtbis.f90
	mv rtbis.o OBJS/
OBJS/twostreights.o: SRC/twostreights.f90 Axis.mk
	$(FC) $(OPTS) -c SRC/twostreights.f90
	mv twostreights.o OBJS/
OBJS/mods4ax.o: SRC/mods4ax.f90 Axis.mk
	$(FC) $(OPTS) -c SRC/mods4ax.f90
	mv mods4ax.o OBJS/
OBJS/in_out_cw.o: SRC/in_out_cw.f90 SRC/axis.f90 Axis.mk
	$(FC) $(OPTS) -c SRC/in_out_cw.f90
	mv in_out_cw.o OBJS/
OBJS/inside_of_wall.o: SRC/inside_of_wall.f90 SRC/axis.f90 Axis.mk
	$(FC) $(OPTS) -c SRC/inside_of_wall.f90
	mv inside_of_wall.o OBJS/
OBJS/axis.o: SRC/axis.f90 Axis.mk
	$(FC) $(OPTS) -c SRC/axis.f90
	mv axis.o OBJS/
