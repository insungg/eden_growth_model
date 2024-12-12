#include "EdenSimulator.hh"
#include <iostream>
#include <algorithm>

EdenSimulator::EdenSimulator(int L, int R, double phi, double mu, double nu, double s)
    : EdenSimulator(L, L, R, R, phi, mu, nu, s) {}

EdenSimulator::EdenSimulator(int L, int Ra, int Rb, double phi, double mu, double nu, double s)
    : EdenSimulator(L, L, Ra, Rb, phi, mu, nu, s) {}

EdenSimulator::EdenSimulator(int rows, int cols, int Ra, int Rb, double phi, double mu, double nu, double s)
    : rows(rows), cols(cols), Ra(Ra), Rb(Rb), phi(phi), mu(mu), nu(nu), s(s), lattice(rows, std::vector<int>(cols, 0)), wall(), hotspots(),
      directions_odd({{-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, 0}, {1, 1}}),
      directions_even({{-1, -1}, {-1, 0}, {0, -1}, {0, 1}, {1, -1}, {1, 0}}),
      gen(std::random_device{}()) {
    for (int j = 0; j < cols; ++j) { 
        lattice[0][j] = (j % 2 == 0) ? 1 : 2;
    }
}

EdenSimulator::~EdenSimulator() {
    lattice.clear();
    hotspots.clear();
    wall.clear();
    directions_odd.clear();
    directions_even.clear();
}

void EdenSimulator::generateHotspots() {
    hotspots = std::vector<std::vector<int>>(rows, std::vector<int>(cols, 0));
    int num_hotspots = static_cast<int>(phi * rows * cols / (Ra * Rb) / 3);
    std::uniform_int_distribution<> dis1(0, rows - 1);
    std::uniform_int_distribution<> dis2(0, cols - 1);

    // see here for hexoganal lattice geometry: https://www.redblobgames.com/grids/hexagons/
    for (int i = 0; i < num_hotspots; ++i) {
        int cx = dis1(gen);
        int cy = dis2(gen);
        std::pair<int, int> center = toAxial(cx, cy);
        int center_q = center.first;
        int center_r = center.second;
        for (int dq = -Ra; dq <= Ra; ++dq) {
            for (int dr = -Rb; dr <= Rb; ++dr) {
                int ds = -dq - dr;
                if (std::abs(dq) <= Ra && std::abs(dr) <= Rb && std::abs(ds) <= Ra) {
                    std::pair<int, int> real = toReal(center_q + dq, center_r + dr);
                    int row = real.first;
                    int col = real.second;
                    if (row >= 0 && row < rows && col >= 0 && col < cols) {
                        hotspots[row][col] = 1;
                    }
                }
            }
        }
    }
}

bool EdenSimulator::isInBulk(int i, int j) {
    std::vector<std::pair<int, int>>& directions = (i % 2 == 0) ? directions_even : directions_odd;
    for (const auto& direction : directions) {
        int ni = i + direction.first;
        int nj = j + direction.second;
        if (ni >= 0 && ni < rows && nj >= 0 && nj < cols && lattice[ni][nj] == 0) {
            return false;
        }
    }
    return true;
}

void EdenSimulator::initializeDomainWall() {
    wall.clear();
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (lattice[i][j] != 0 && !isInBulk(i, j)) {
                wall.emplace(i, j);
            }
        }
    }
}

void EdenSimulator::updateDomainWall(int i, int j, int ni, int nj) {
    // remove (i, j) from the wall if it is in the bulk
    if (isInBulk(i, j)) {
        wall.erase(std::make_pair(i, j));
    }
    // check if (ni, nj) is in wall, along with that
    // neighbor of (ni, nj) would be no longer in the wall due to filling of (ni, nj)
    std::vector<std::pair<int, int>>& directions_n = (ni % 2 == 0) ? directions_even : directions_odd;
    for (const auto& direction : directions_n) {
        int nni = ni + direction.first;
        int nnj = nj + direction.second;
        if (nni >= 0 && nni < rows && nnj >= 0 && nnj < cols) {
            if (lattice[nni][nnj] == 0) {
                wall.emplace(ni, nj);
            } else if (isInBulk(nni, nnj)) {
                    wall.erase(std::make_pair(nni, nnj));
            }
        }
    }
}

void EdenSimulator::updateLattice() {
    std::vector<double> weights;
    for (const auto& [i, j] : wall) {
        if (lattice[i][j] == 0) {
            std::cerr << "Error: Wall contains an empty site." << std::endl;
            return;
        }
        double growth_rate = (lattice[i][j] == 1) ? 1.0 : (1.0 - s);
        growth_rate *= (hotspots[i][j] == 1) ? (1.0 + nu) : 1.0;
        weights.emplace_back(growth_rate);
    }

    // select one element from the wall following the weights
    std::discrete_distribution<> dist(weights.begin(), weights.end());
    int selected_index = dist(gen);
    auto it = std::next(wall.begin(), selected_index);
    int i = it->first;
    int j = it->second;

    // select a random empty neighbor site and fill it
    std::vector<std::pair<int, int>> empty_neighbors;
    std::vector<std::pair<int, int>>& directions = (i % 2 == 0) ? directions_even : directions_odd;
    for (const auto& direction : directions) {
        int ni = i + direction.first;
        int nj = j + direction.second;
        if (ni >= 0 && ni < rows && nj >= 0 && nj < cols && lattice[ni][nj] == 0) {
            empty_neighbors.emplace_back(ni, nj);
        }
    }
    if (!empty_neighbors.empty()) {
        std::uniform_int_distribution<> empty_dis(0, empty_neighbors.size() - 1);
        auto [new_i, new_j] = empty_neighbors[empty_dis(gen)];
        if (lattice[i][j] == 1) {
            // prevent double precision error when mu is zero
            if (mu > 0 && mu <= 1) { 
                std::bernoulli_distribution mutation(mu);
                lattice[new_i][new_j] = mutation(gen) ? 2 : 1;
            } else {
                lattice[new_i][new_j] = 1;
            }
        } else if (lattice[i][j] == 2) {
            lattice[new_i][new_j] = 2;
        }
        if (lattice[i][j] != 0  && lattice[new_i][new_j] != 0) {
            updateDomainWall(i, j, new_i, new_j); // update wall only when the move is successful
        }
    } else {
        std::cerr << "Error: No empty neighbors found." << std::endl;
        return;
    }
}

void EdenSimulator::generateSnapshot() {
    // validity checks for hotspots and initial lattice
    std::cout << "Initializing landscape generation with parameters nu: " << nu << " and s: " << s << std::endl;
    if (std::all_of(lattice.begin(), lattice.end(), [](const std::vector<int>& row) {
        return std::all_of(row.begin(), row.end(), [](int site) { return site == 0; });
    })) {
        throw std::runtime_error("Error: All lattice sites are zero. Cannot generate landscape.");
    }
    if (std::all_of(hotspots.begin(), hotspots.end(), [](const std::vector<int>& row) {
        return std::all_of(row.begin(), row.end(), [](int site) { return site == 0; });
    })) {
        throw std::runtime_error("Error: All hotspots are zero. Cannot generate landscape.");
    }
    unsigned long long iter = 0;
    initializeDomainWall();
    while (!wall.empty()) {
        updateLattice();
        iter++;
        if (iter % 10000 == 0) { // check mutant extinction
            int mutant_count = std::count_if(wall.begin(), wall.end(), [&](const std::pair<int, int>& pos) {
                return lattice[pos.first][pos.second] == 2;
            });
            if (mu <= 0 && mutant_count == 0) {
                for (auto& row : lattice) {
                    std::replace(row.begin(), row.end(), 0, 1);
                }
                std::cout << "All species on the domain wall are wild with no mutation allowed. No more updates and simulation terminated." << std::endl;
                break;
            }
        }
        if (iter >= 100 * rows*cols) {
            std::cerr << "Error: Simulation takes too long. Current wall size: " << wall.size() <<  ". Terminating the program." << std::endl;
            break;
        }
    }
    std::cout << "Finished generating landscape after " << iter << " iterations." << std::endl;
}

std::pair<int, int> EdenSimulator::toAxial(int row, int col) {
    int q = col - (row / 2);
    int r = row;
    return std::make_pair(q, r);
}

std::pair<int, int> EdenSimulator::toReal(int q, int r) {
    int col = q + (r / 2);
    int row = r;
    return std::make_pair(row, col);
}

void EdenSimulator::setHotspots(const std::vector<std::vector<int>>& newHotspots) {
    if (newHotspots.size() != static_cast<size_t>(rows) || 
        (!newHotspots.empty() && newHotspots[0].size() != static_cast<size_t>(cols))) {
        throw std::invalid_argument("New hotspots size does not match the dimensions of the lattice");
    }
    hotspots = newHotspots;
}

void EdenSimulator::setInitialCondition(const std::vector<int>& initialCondition) {
    if (initialCondition.size() != static_cast<size_t>(cols)) {
        throw std::invalid_argument("Initial condition size does not match the number of columns");
    }
    lattice[0] = initialCondition;
}

const std::vector<std::vector<int>>& EdenSimulator::getLattice() const {
    return lattice;
}

const std::vector<std::vector<int>>& EdenSimulator::getHotspots() const {
    return hotspots;
}

const std::set<std::pair<int, int>>& EdenSimulator::getWall() const {
    return wall;
}