ELIBS=$(EPIPHANY_HOME)/tools/host/lib
EINCS=$(EPIPHANY_HOME)/tools/host/include

ELDF=$(EPIPHANY_HOME)/bsps/current/internal.ldf
#ELDF=$(EPIPHANY_HOME)/bsps/current/fast.ldf
#ELDF=$(EPIPHANY_HOME)/bsps/current/legacy.ldf

EXES=main core.srec
OBJS=core.elf

all: $(EXES)


main: main.c common.h
	gcc -O3 -funroll-loops -ffast-math -Wall main.c -o ./bin/main -I $(EINCS) -L $(ELIBS) -le-hal -lrt
	#gcc -Wall main.c -o ./bin/main -I $(EINCS) -L $(ELIBS) -le-hal -lrt

core.elf: core.c common.h
	e-gcc -O3 -Wall -funroll-loops -ffast-math -T $(ELDF) core.c -o ./bin/core.elf -le-lib
	#e-gcc -O2 -funroll-loops -ffast-math -T $(ELDF) core.c -o ./bin/core.elf -le-lib
		

core.srec: core.elf
	e-objcopy --srec-forceS3 --output-target srec ./bin/core.elf ./bin/core.srec

clean:
	rm $(EXES) $(OBJS)
