all:
	make compile;make run;make raps; make clean

run:
	./ADOL-C/direct;./ADOL-C/indirect;cd ReverseAD;./reversead;cd ..;

compile:
	make -C ADOL-C; make -C ReverseAD

raps:
	make -C Rapsodia

clean:
	make -C ADOL-C clean; make -C ReverseAD clean; make -C Rapsodia clean
