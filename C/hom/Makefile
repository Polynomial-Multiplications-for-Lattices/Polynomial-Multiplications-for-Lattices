
CC = gcc

CFLAGS += -O0 -Wall -Wunused -Wpedantic -Wno-pointer-arith

COMMON_PATH = ../common

CFLAGS += -I$(COMMON_PATH)

COMMON_SOURCE = $(COMMON_PATH)/tools.c $(COMMON_PATH)/naive_mult.c $(COMMON_PATH)/gen_table.c $(COMMON_PATH)/ntt_c.c

ASM_HEADERs =
ASM_SOURCEs =
C_HEADERs =
C_SOURCEs =

SOURCEs = $(ASM_SOURCEs) $(C_SOURCEs) $(COMMON_SOURCE)
HEADERs = $(ASM_HEADERs) $(C_HEADERs)

all: Karatsuba Karatsuba-striding TC TC-striding Toeplitz-TC DWT FNT GT Nussbaumer Schoenhage

Karatsuba: Karatsuba.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Karatsuba-striding: Karatsuba-striding.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

TC: TC.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

TC-striding: TC-striding.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Toeplitz-TC: Toeplitz-TC.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

DWT: DWT.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

FNT: FNT.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

GT: GT.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Nussbaumer: Nussbaumer.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Schoenhage: Schoenhage.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@


.PHONY: clean
clean:
	rm -f Karatsuba
	rm -f Karatsuba-striding
	rm -f TC
	rm -f TC-striding
	rm -f Toeplitz-TC
	rm -f DWT
	rm -f FNT
	rm -f GT
	rm -f Nussbaumer
	rm -f Schoenhage




