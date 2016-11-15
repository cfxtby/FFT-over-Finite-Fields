NODES=48
MPI=/gpfs1/openmpi/bin/mpic++
MPI_RUN=/gpfs1/openmpi/bin/mpirun

serial: serialFFT.o time.o
	g++ -o FFT serialFFT.o time.o

parallel: parCut.o time.o
	${MPI} -o parFFT parCut.o time.o

parCut.o: parCut.c
	${MPI} -c parCut.c 

runPar:
	bsub -oo log.txt -K -q high_priority -n 4 -R"span[ptile=1]" ${MPI_RUN} ./parFFT 1

serialFFT.o: serialFFT.c
	g++ -c serialFFT.c

time.o: time.cpp
	g++ -c time.cpp

clean:
	rm -rf *.o *~ FFT parRep

run:
	./FFT
