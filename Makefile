chromosome_maps: main.c messages.c genes.c
	gcc `pkg-config gtk+-3.0 --cflags` main.c messages.c genes.c -Wall -o chromosome_maps -lm -ggdb3 `pkg-config gtk+-3.0 --libs` -rdynamic