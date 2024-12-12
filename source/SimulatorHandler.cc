#include "SimulatorHandler.hh"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <thread>
#include <mutex>
#include <iostream>

SimulatorHandler::SimulatorHandler(const std::string& paramFile, const std::string& outputDir)
    : outputDir(outputDir), hotspotsInitialized(false) {
    loadParameters(paramFile);
}

SimulatorHandler::~SimulatorHandler() {
    for (auto& snapshot : snapshots) {
        snapshot.clear();
    }
    snapshots.clear();
    hotspots.clear();
    initialCondition.clear();
    initialHotspots.clear();
    nu_values.clear();
    s_values.clear();
}

void SimulatorHandler::loadParameters(const std::string& paramFile) {
    std::ifstream file(paramFile);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open parameter file");
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key;
        if (std::getline(iss, key, '=')) {
            std::string value;
            if (std::getline(iss, value)) {
                if (key == "numSnapshots") numSnapshots = std::stoi(value);
                else if (key == "rows") rows = std::stoi(value);
                else if (key == "cols") cols = std::stoi(value);
                else if (key == "Ra") Ra = std::stoi(value);
                else if (key == "Rb") Rb = std::stoi(value);
                else if (key == "phi") phi = std::stod(value);
                else if (key == "mu") mu = std::stod(value);
                else if (key == "parallelOn") parallelOn = (value == "true");
                else if (key == "initialCondition") {
                    std::istringstream icStream(value);
                    int ic;
                    while (icStream >> ic) {
                        initialCondition.push_back(ic);
                        if (icStream.peek() == ',') icStream.ignore();
                    }
                }
                else if (key == "initialHotspots") {
                    std::istringstream hsStream(value);
                    std::vector<int> row;
                    int hs;
                    while (hsStream >> hs) {
                        row.push_back(hs);
                        if (hsStream.peek() == ',') hsStream.ignore();
                        if (hsStream.peek() == ';') {
                            initialHotspots.push_back(row);
                            row.clear();
                            hsStream.ignore();
                        }
                    }
                    if (!row.empty()) {
                        initialHotspots.push_back(row);
                    }
                }
                else if (key == "numCores") numCores = std::stoi(value);
                else if (key == "nu") {
                    std::istringstream nuStream(value);
                    double nu;
                    while (nuStream >> nu) {
                        nu_values.push_back(nu);
                        if (nuStream.peek() == ',') nuStream.ignore();
                    }
                }
                else if (key == "s") {
                    std::istringstream sStream(value);
                    double s;
                    while (sStream >> s) {
                        s_values.push_back(s);
                        if (sStream.peek() == ',') sStream.ignore();
                    }
                }
                else if (key == "fileName") {
                    fileName = value;
                }
            }
        }
        if (!initialCondition.empty() && initialCondition.size() != static_cast<size_t>(cols)) {
            throw std::invalid_argument("Initial condition size does not match the number of columns");
        }
        if (!initialHotspots.empty() && initialHotspots.size() != static_cast<size_t>(rows)
            && initialHotspots[0].size() != static_cast<size_t>(cols)) {
            throw std::invalid_argument("Initial hotspots size does not match the number of rows and columns");
        }
    }
}

void SimulatorHandler::initializeHotspots() {
    EdenSimulator initialSimulation(rows, cols, Ra, Rb, phi, mu, 0, 0); // nu and s are not needed for hotspot generation
    if (!initialHotspots.empty()) {
        initialSimulation.setHotspots(initialHotspots);
    } else {
        initialSimulation.generateHotspots();
    }
    hotspots = initialSimulation.getHotspots();
    hotspotsInitialized = true;
}

void SimulatorHandler::run(double nu, double s) {
    if (!hotspotsInitialized) {
        initializeHotspots();
    }

    snapshots.clear();
    for (int i = 0; i < numSnapshots; ++i) {
        EdenSimulator simulation(rows, cols, Ra, Rb, phi, mu, nu, s);
        if (!initialCondition.empty()) {
            simulation.setInitialCondition(initialCondition);
        }
        simulation.setHotspots(hotspots);
        simulation.generateSnapshot();
        snapshots.emplace_back(simulation.getLattice());
    }
    saveResults(fileName, nu, s);
}

void SimulatorHandler::runParallel(double nu, double s) {
    if (!hotspotsInitialized) {
        initializeHotspots();
    }

    snapshots.clear();
    std::mutex mtx;
    auto worker = [&](int start, int end) {
        for (int i = start; i < end; ++i) {
            EdenSimulator simulation(rows, cols, Ra, Rb, phi, mu, nu, s);
            if (!initialCondition.empty()) {
                simulation.setInitialCondition(initialCondition);
            }
            simulation.setHotspots(hotspots);
            simulation.generateSnapshot();
            std::lock_guard<std::mutex> lock(mtx);
            snapshots.emplace_back(simulation.getLattice());
        }
    };

    std::vector<std::thread> threads;
    int snapshotsPerCore = numSnapshots / numCores;
    int remainingSnapshots = numSnapshots % numCores;

    for (int i = 0; i < numCores; ++i) {
        int start = i * snapshotsPerCore;
        int end = start + snapshotsPerCore + (i < remainingSnapshots ? 1 : 0);
        threads.emplace_back(worker, start, end);
        std::this_thread::sleep_for(std::chrono::milliseconds(449)); // 449 is a prime number
    }

    for (auto& thread : threads) {
        thread.join();
    }
    saveResults(fileName, nu, s);
}

void SimulatorHandler::runParameterSpace() {
    for (double nu : nu_values) {
        for (double s : s_values) {
            if (parallelOn) {
                runParallel(nu, s);
            } else {
                run(nu, s);
            }
        }
    }
}

void SimulatorHandler::saveResults(const std::string& filenamePrefix, double nu, double s) const {
    std::string filename = outputDir + "/" + filenamePrefix + "_N_" + std::to_string(numSnapshots) + "_rows_" + std::to_string(rows) + "_cols_" + std::to_string(cols) +
                           "_Ra_" + std::to_string(Ra) + "_Rb_" + std::to_string(Rb) + "_phi_" + std::to_string(phi) +
                           "_mu_" + std::to_string(mu) + "_nu_" + std::to_string(nu) + "_s_" + std::to_string(s) + ".bin";
    std::ofstream file(filename, std::ios::binary);
    if (file.is_open()) {
        for (const auto& snapshot : snapshots) {
               for (const auto& row : snapshot) {
                 std::vector<uint8_t> byte_row(row.begin(), row.end());
                 file.write(reinterpret_cast<const char*>(byte_row.data()), byte_row.size() * sizeof(uint8_t));
             }
        }
        // Save the hotspot configuration at the end of the file
        for (const auto& row : hotspots) {
            std::vector<uint8_t> byte_row(row.begin(), row.end());
                file.write(reinterpret_cast<const char*>(byte_row.data()), byte_row.size() * sizeof(uint8_t));
            }
          file.close();
        std::cout << "File has been successfully saved." << std::endl;
    } else {
            std::cerr << "Unable to open file: " << filename << std::endl;
    }
}