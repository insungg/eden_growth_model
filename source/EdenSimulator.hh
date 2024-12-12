#ifndef EDEN_SIMULATOR_HH
#define EDEN_SIMULATOR_HH

#include <vector>
#include <utility>
#include <random>
#include <set>

class EdenSimulator {
public:
    EdenSimulator(int L, int R, double phi, double mu, double nu, double s);
    EdenSimulator(int L, int Ra, int Rb, double phi, double mu, double nu, double s);
    EdenSimulator(int rows, int cols, int Ra, int Rb, double phi, double mu, double nu, double s);
    ~EdenSimulator(); // Destructor to free resources
    void generateSnapshot();
    void generateHotspots();

    void setHotspots(const std::vector<std::vector<int>>& hotspots);
    void setInitialCondition(const std::vector<int>& initialCondition); // New method

    const std::vector<std::vector<int>>& getLattice() const;
    const std::vector<std::vector<int>>& getHotspots() const;
    const std::set<std::pair<int, int>>& getWall() const;

private:
    int rows, cols, Ra, Rb;
    double phi, mu, nu, s;
    std::vector<std::vector<int>> lattice;
    std::vector<std::vector<int>> hotspots;
    std::set<std::pair<int, int>> wall; // occupied cells with at least one unoccupied neighbor
    std::vector<std::pair<int, int>> directions_odd;
    std::vector<std::pair<int, int>> directions_even;
    std::mt19937 gen;

    bool isInBulk(int i, int j); 
    void initializeDomainWall();
    void updateDomainWall(int i, int j, int ni, int nj);
    void updateLattice();
    std::pair<int, int> toAxial(int row, int col);
    std::pair<int, int> toReal(int q, int r);
};

#endif