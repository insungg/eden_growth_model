#ifndef RUN_SIMULATION_HH
#define RUN_SIMULATION_HH

#include "EdenSimulator.hh"
#include <vector>
#include <string>

class SimulatorHandler {
public:
    SimulatorHandler(const std::string& paramFile, const std::string& outputDir);
    void loadParameters(const std::string& paramFile);
    void run(double nu, double s);
    void runParallel(double nu, double s);
    void runParameterSpace();
    void saveResults(const std::string& filenamePrefix, double nu, double s) const;
    ~SimulatorHandler(); // Destructor to free resources

private:
    int numSnapshots;
    int rows, cols, Ra, Rb, numCores;
    double phi, mu;
    bool parallelOn;
    bool hotspotsInitialized;
    std::string outputDir;
    std::string fileName; 
    std::vector<int> initialCondition;
    std::vector<std::vector<int>> initialHotspots; 
    std::vector<double> nu_values;
    std::vector<double> s_values;
    std::vector<std::vector<std::vector<int>>> snapshots;
    std::vector<std::vector<int>> hotspots;
    void initializeHotspots();
};

#endif