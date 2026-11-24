CC ?= cc
CFLAGS ?= -std=c99 -O2 -Wall -Wextra -pedantic
LDFLAGS ?=
LDLIBS ?= -lm
INSTALL_PATH ?= /usr/local
BINDIR ?= $(INSTALL_PATH)/bin

TARGET = get_newick_info

all: $(TARGET)

$(TARGET): get_newick_info.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS)

install: $(TARGET)
	install -d $(DESTDIR)$(BINDIR)
	install -m 755 $(TARGET) $(DESTDIR)$(BINDIR)/$(TARGET)

uninstall:
	rm -f $(DESTDIR)$(BINDIR)/$(TARGET)

clean:
	rm -f $(TARGET)

.PHONY: all install uninstall clean
