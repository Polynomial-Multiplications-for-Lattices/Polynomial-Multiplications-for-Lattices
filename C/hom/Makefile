
CC = gcc

CFLAGS += -O0 -Wall -Wunused -Wpedantic -Wno-pointer-arith -Wsign-compare

COMMON_PATH = ../common

CFLAGS += -I$(COMMON_PATH)

COMMON_SOURCE = $(COMMON_PATH)/tools.c $(COMMON_PATH)/naive_mult.c $(COMMON_PATH)/gen_table.c $(COMMON_PATH)/ntt_c.c

ASM_HEADERs =
ASM_SOURCEs =
C_HEADERs =
C_SOURCEs =

SOURCEs = $(ASM_SOURCEs) $(C_SOURCEs) $(COMMON_SOURCE)
HEADERs = $(ASM_HEADERs) $(C_HEADERs)

all: DWT DWT_merged_layers FNT GT Karatsuba Karatsuba-striding Karatsuba-striding_multi-layer Nussbaumer Schoenhage TC TC-striding Toeplitz-TC

DWT: DWT.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

DWT_merged_layers: DWT_merged_layers.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

FNT: FNT.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

GT: GT.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Karatsuba: Karatsuba.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Karatsuba-striding: Karatsuba-striding.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Karatsuba-striding_multi-layer: Karatsuba-striding_multi-layer.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Nussbaumer: Nussbaumer.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Schoenhage: Schoenhage.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

TC: TC.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

TC-striding: TC-striding.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Toeplitz-TC: Toeplitz-TC.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@


.PHONY: clean
clean:
	rm -f DWT
	rm -f DWT_merged_layers
	rm -f FNT
	rm -f GT
	rm -f Karatsuba
	rm -f Karatsuba-striding
	rm -f Nussbaumer
	rm -f Schoenhage
	rm -f TC
	rm -f TC-striding
	rm -f Toeplitz-TC





