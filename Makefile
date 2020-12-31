## project directories
SRC_DIR:=src
INCLUDE_DIR:=include
BUILD_DIR:=build
OUT_DIR=outfiles

## target executable
EXE:=lifeMPI

## set compiler, linker and preprocessor
CXX:=mpicxx
CXXFLAGS:=-Wall -Wextra
CPPFLAGS:=-MMD -MP
LDFLAGS:=
LDLIBS:=

# -----------------------------------------------------------------

SRC_FILES:=$(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES:=$(addprefix $(BUILD_DIR)/,$(addsuffix .o,$(notdir $(basename $(SRC_FILES)))))
DEP_FILES:=$(OBJ_FILES:.o=.d)

.PHONY : all
all : $(EXE)

$(EXE) : $(OBJ_FILES) | outfiles
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

.PHONY : clean
clean:
	-$(RM) $(EXE) $(OBJ_FILES) $(DEP_FILES)
	rm -rf $(BUILD_DIR) $(OUT_DIR) $(IN_DIR)	 

.PHONY : build
build:
	mkdir -p $(BUILD_DIR)

.PHONY : outfiles
outfiles:
	mkdir -p $(OUT_DIR)

.PHONY : test
test:
	python scripts/test_life.py

-include $(DEP_FILES)

$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp Makefile | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@
