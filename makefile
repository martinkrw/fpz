ASM=nasm
ASMFLAGS=-f elf64 -Wall -I include/
CXX=clang++
CXXFLAGS=-pedantic -std=c++2a -Wall -O3 -I include/
LIBS=-lm

INSTALL_PREFIX=/usr/local

TARGET=`uname`

# Compilation on linux requires libbsd
.if ${TARGET} == "Linux"
CXXFLAGS+= `pkg-config libbsd-overlay --cflags`
LIBS+= `pkg-config libbsd-overlay --libs`
.endif

ASMNAMES=add_sub mul get_ld div shift
CXXNAMES=cmp mul_cpp div_cpp debug assign io
EXNAMES=ex1 ex2 ex3 ex4
TESTNAMES=test_add_sub test_mul test_div

# Prefix names with directories and file extension (if applicable)
ASMOBJS=${ASMNAMES:C/.*/build\/&.o/}
CXXOBJS=${CXXNAMES:C/.*/build\/&.o/}
EXAMPLES=${EXNAMES:C/.*/bin\/&/}
TESTS=${TESTNAMES:C/.*/bin\/&/}

all: lib/libfpz.a ${EXAMPLES} ${TESTS}

lib/libfpz.a: ${ASMOBJS} ${CXXOBJS}
	ar rcs $@ $>

${ASMOBJS}: src/${*F}.asm include/asmmacros.inc
	${ASM} ${ASMFLAGS} -o $@ src/${*F}.asm
${CXXOBJS}: src/${*F}.cpp include/fpz_core.h
	${CXX} ${CXXFLAGS} -c -o $@ src/${*F}.cpp
${EXAMPLES}: examples/${*F}.cpp include/fpz.h lib/libfpz.a
	${CXX} ${CXXFLAGS} -o $@ examples/${*F}.cpp lib/libfpz.a ${LIBS}
${TESTS}: tests/${*F}.cpp include/fpz_core.h lib/libfpz.a include/test.h
	${CXX} ${CXXFLAGS} -o $@ tests/${*F}.cpp lib/libfpz.a ${LIBS}

.PHONY: clean all install uninstall run_tests

run_tests: ${TESTS}
	@for test in ${TESTS}; do $$test; done

clean:
	@rm -f lib/libfpz.a 
	@rm -f ${ASMOBJS} ${CXXOBJS} 
	@rm -f ${EXAMPLES} ${TESTS}

install: lib/libfpz.a
	install -C -m 644 lib/libfpz.a ${INSTALL_PREFIX}/lib/
	install -C -m 644 include/fpz.h ${INSTALL_PREFIX}/include/
	install -C -m 644 include/fpz_core.h ${INSTALL_PREFIX}/include/

uninstall:
	rm -f ${INSTALL_PREFIX}/lib/libfpz.a
	rm -f ${INSTALL_PREFIX}/include/fpz.h
	rm -f ${INSTALL_PREFIX}/include/fpz_core.h	
