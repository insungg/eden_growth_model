CXX = g++
CXXFLAGS = -std=c++11 -O2 -pthread

TARGET = ../programs/eden_simulator
SRCS = main.cc EdenSimulator.cc SimulatorHandler.cc
OBJS = $(SRCS:.cc=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	@if [ ! -d "../programs" ]; then mkdir -p ../programs; fi
	@if [ ! -d "../data" ]; then mkdir -p ../data; fi
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)
	rm -f $(OBJS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)