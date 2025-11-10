# Compiler and flags
CXX ?= g++
CXXFLAGS ?= -O3 -std=c++20 -Iinclude -MMD -MP

# Dirs
BIN_DIR := bin
SRC_DIR := src
BUILD_DIR := build

# Executables
BINS := $(BIN_DIR)/sketch \
        $(BIN_DIR)/filter \
        $(BIN_DIR)/pwsimilarity\
		$(BIN_DIR)/expt_growth

# Sources/objects
SRCS := $(SRC_DIR)/sketch_main.cpp \
        $(SRC_DIR)/filter_main.cpp \
        $(SRC_DIR)/pwsim_main.cpp \
        $(SRC_DIR)/expt_growth.cpp

OBJS := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

# Default target
all: $(BINS)

# Ensure directories exist
$(BIN_DIR) $(BUILD_DIR):
	mkdir -p $@

# Compile objects
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link binaries
$(BIN_DIR)/sketch: $(BUILD_DIR)/sketch_main.o | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BIN_DIR)/filter: $(BUILD_DIR)/filter_main.o | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BIN_DIR)/pwsimilarity: $(BUILD_DIR)/pwsim_main.o | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BIN_DIR)/expt_growth: $(BUILD_DIR)/expt_growth.o | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

# Clean
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

# Include auto-generated dependencies
-include $(DEPS)

.PHONY: all clean
