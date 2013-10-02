# Ludwig Schmidt (ludwigschmidt2@gmail.com) 2013
#
# This makefile is based on http://make.paulandlesley.org/autodep.html .

CXX = g++
MEX = mex
CXXFLAGS = -Wall -Wextra -O2 -std=c++98 -ansi -fPIC
MEXCXXFLAGS = -Wall -Wextra -O2 -std=c++98 -ansi

SRCDIR = src
DEPDIR = .deps
OBJDIR = obj

SRCS = main.cc emd_flow.cc emd_flow_network.cc emd_flow_network_sap.cc

.PHONY: clean archive

clean:
	rm -rf $(OBJDIR)
	rm -rf $(DEPDIR)
	rm -f emd_flow
	rm -f emd_flow.mexa64
	rm -f emd_flow.mexmaci64
	rm -f emd_flow.tar.gz

archive:
	mkdir archive-tmp
	tar --transform='s,^\.,emd_flow,' --exclude='.git' --exclude='archive-tmp' -czf archive-tmp/emd_flow.tar.gz .
	mv archive-tmp/emd_flow.tar.gz .
	rm -rf archive-tmp

# emd_flow executable
EMD_FLOW_OBJECTS = main.o emd_flow.o emd_flow_network_factory.o emd_flow_network_sap.o

emd_flow: $(EMD_FLOW_OBJECTS:%=$(OBJDIR)/%)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_program_options


# emd_flow MEX file
MEXFILE_OBJECTS = emd_flow.o emd_flow_network_factory.o emd_flow_network_sap.o
MEXFILE_SRC = mex_wrapper.cc
MEXFILE_SRC_DEPS = $(MEXFILE_SRC) mex_helper.h emd_flow.h emd_flow_network_factory.h

mexfile: $(MEXFILE_OBJECTS:%=$(OBJDIR)/%) $(MEXFILE_SRC_DEPS:%=$(SRCDIR)/%)
	$(MEX) -v CXXFLAGS="\$$CXXFLAGS $(MEXCXXFLAGS)" -output emd_flow $(SRCDIR)/$(MEXFILE_SRC) $(MEXFILE_OBJECTS:%=$(OBJDIR)/%)


$(OBJDIR)/%.o: $(SRCDIR)/%.cc
  # Create the directory the current target lives in.
	@mkdir -p $(@D)
  # Compile and generate a dependency file.
  # See http://gcc.gnu.org/onlinedocs/gcc/Preprocessor-Options.html .
	$(CXX) $(CXXFLAGS) -MMD -MP -c -o $@ $<
  # Move dependency file to dependency file directory.
  # Create the dependency file directory if necessary.
	@mkdir -p $(DEPDIR)
	@mv $(OBJDIR)/$*.d $(DEPDIR)/$*.d

# Include the generated dependency files.
# The command replaces each file name in SRCS with its dependency file.
# See http://www.gnu.org/software/make/manual/html_node/Substitution-Refs.html#Substitution-Refs for the GNU make details.
-include $(SRCS:%.cc=$(DEPDIR)/%.d)
