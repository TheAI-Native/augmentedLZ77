CC ?= cc
CFLAGS ?= -O2 -std=c99 -Wall -Wextra -pedantic
LDLIBS ?= -lm
TARGET := lzp_codec
SOURCE := src/lzp_codec.c

.PHONY: all test clean

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

test: $(TARGET)
	sh scripts/smoke_test.sh

clean:
	rm -f $(TARGET) tests/data/sample.lzp tests/data/sample.restored.txt
