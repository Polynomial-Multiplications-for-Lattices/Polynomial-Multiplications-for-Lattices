
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

all: Montgomery_acc Montgomery_sub Barrett Barrett_Montgomery_cmp

Montgomery_acc: Montgomery_acc.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Montgomery_sub: Montgomery_sub.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Barrett: Barrett.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@

Barrett_Montgomery_cmp: Barrett_Montgomery_cmp.c $(SOURCEs) $(HEADERs)
	$(CC) $(CFLAGS) $(SOURCEs) $< -o $@


.PHONY: clean
clean:
	rm -f Montgomery_acc
	rm -f Montgomery_sub
	rm -f Barrett
	rm -f Barrett_Montgomery_cmp




