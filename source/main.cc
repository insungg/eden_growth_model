#include "SimulatorHandler.hh"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <parameter_file_path> <output_directory>" << std::endl;
        return 1;
    }

    std::string parameterFilePath = argv[1];
    std::string outputDir = argv[2];
    SimulatorHandler simulatorHandler(parameterFilePath, outputDir);
    simulatorHandler.runParameterSpace();

    std::cout << "Simulation successfully completed. Goodbye:)" << std::endl;
    return 0;
}