ELIBS=$(EPIPHANY_HOME)/tools/host/lib
EINCS=$(EPIPHANY_HOME)/tools/host/include

ELDF=$(EPIPHANY_HOME)/bsps/current/internal.ldf
#ELDF=$(EPIPHANY_HOME)/bsps/current/fast.ldf
#ELDF=$(EPIPHANY_HOME)/bsps/current/legacy.ldf


EXES=main core.srec
OBJS=core.elf

all: $(EXES)


main: main.c common.h
	mpicc -O3 -funroll-loops -ffast-math -Wall main.c -o ./bin/main -I $(EINCS) -L $(ELIBS) -le-hal -lrt
	scp ./bin/main node1:SDU/cluster/bin
	scp ./bin/main node2:SDU/cluster/bin
	scp ./bin/main node3:SDU/cluster/bin
	

core.elf: core.c common.h
	e-gcc -O3 -Wall -funroll-loops -ffast-math -T $(ELDF) core.c -o ./bin/core.elf -le-lib


core.srec: core.elf
	e-objcopy --srec-forceS3 --output-target srec ./bin/core.elf ./bin/core.srec
	scp ./bin/core.srec node1:SDU/cluster/bin
	scp ./bin/core.srec node2:SDU/cluster/bin
	scp ./bin/core.srec node3:SDU/cluster/bin

clean:
	rm $(EXES) $(OBJS)
