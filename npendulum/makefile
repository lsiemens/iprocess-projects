pendulum.exe :
	gfortran types.f90 pyio.f90 modules.f90 npendulum.f90 -o pendulum.exe -lblas -llapack

clear :
	-rm *.exe *.mp4

run :
	python setup_simulation.py
	./pendulum.exe
	python plot.py
	gnome-open anim.mp4
