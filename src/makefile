# just the build-commands

I =  -I../../../bin/pd/src

all: linux

linux:
	cc -DPD -O2 -DUNIX $(I) -o dtw.o -c dtw.c
	ld -shared --export-dynamic -o dtw.pd_linux dtw.o

nt:
	cl dtw.c /DNT /DPD /LD /Gd /GD /Ox /link /export:dtw_setup \system\pd\bin\pd.lib
