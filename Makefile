# Makefile for Forensic DNA Mixture Analysis Tool
# Compiler and flags
CC = gcc
CFLAGS = -Wall -O2 -lm
DEBUG_FLAGS = -g -Wall -lm

# Target executable
TARGET = testfunc

# Source files
SOURCES = testfunc.c json_parser.c neldermead.c
HEADERS = json_parser.h neldermead.h

# Object files
OBJECTS = $(SOURCES:.c=.o)

# Default target
all: $(TARGET)

# Link object files to create executable
$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(TARGET) $(CFLAGS)

# Compile source files to object files
%.o: %.c $(HEADERS)
	$(CC) -c $< -o $@ -Wall

# Debug build
debug: CFLAGS = $(DEBUG_FLAGS)
debug: clean $(TARGET)

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(TARGET)

# Remove backup files
cleanall: clean
	rm -f *.backup *~

# Install (optional - copies to /usr/local/bin)
install: $(TARGET)
	install -m 755 $(TARGET) /usr/local/bin/

# Uninstall
uninstall:
	rm -f /usr/local/bin/$(TARGET)

# Test with sample data (requires a.json in current directory)
test: $(TARGET)
	./$(TARGET) -i a.json -l L1

# Help
help:
	@echo "Forensic DNA Mixture Analysis - Makefile Help"
	@echo ""
	@echo "Available targets:"
	@echo "  make           - Build the program (default)"
	@echo "  make debug     - Build with debug symbols"
	@echo "  make clean     - Remove compiled objects and executable"
	@echo "  make cleanall  - Remove all build artifacts and backups"
	@echo "  make install   - Install to /usr/local/bin (requires sudo)"
	@echo "  make uninstall - Remove from /usr/local/bin (requires sudo)"
	@echo "  make test      - Build and run test with a.json"
	@echo "  make help      - Show this help message"
	@echo ""
	@echo "Example usage:"
	@echo "  make clean && make"
	@echo "  sudo make install"

.PHONY: all debug clean cleanall install uninstall test help
