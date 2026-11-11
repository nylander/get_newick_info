CC ?= cc
CFLAGS ?= -std=c99 -O2 -Wall -Wextra -pedantic
LDFLAGS ?=
LDLIBS ?= -lm

TARGET = get_newick_info

all: $(TARGET)

$(TARGET): get_newick_info.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f $(TARGET)

test: $(TARGET)
	python3 test_get_newick_info.py

.PHONY: all clean test
