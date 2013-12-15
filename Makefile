# Makefile for the wave_equation package
# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# Time-stamp: <2013-12-11 17:34:31 (jonah)>

# The default compiler is g++
CXX = g++

# Flags for the compiler. Ask for warnings. Enable the debugger.
CXXFLAGS = -Wall -g

default: travelling_wave
all: wave_equation_test_driver standing_wave travelling_wave
test: wave_equation_test_driver

travelling_wave: travelling_wave.bin
travelling_wave.bin: travelling_wave.o wave_equation.o rkf45.o
	$(CXX) $(CXXFLAGS) -o $@ $^
travelling_wave.o: wave_equation.hpp rkf45.hpp

standing_wave: standing_wave.bin
standing_wave.bin: standing_wave.o wave_equation.o rkf45.o
	$(CXX) $(CXXFLAGS) -o $@ $^
standing_wave.o: wave_equation.hpp rkf45.hpp

wave_equation_test_driver: wave_equation_test_driver.bin
wave_equation_test_driver.bin: wave_equation.o wave_equation_test_driver.o rkf45.o
	$(CXX) $(CXXFLAGS) -o $@ $^
wave_equation_test_driver.o: wave_equation.hpp rkf45.hpp

wave_equation.o: wave_equation.hpp rkf45.hpp
rkf45.o: rkf45.hpp

.PHONY: default all test wave_equation_test_driver standing_wave travelling_wave

clean:
	$(RM) rkf45.o wave_equation.o wave_equation_test_driver.o wave_equation_test_driver.bin standing_wave.bin stading_wave.o background.aux background.log background.out