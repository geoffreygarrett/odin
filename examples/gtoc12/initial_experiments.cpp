#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Data {
    int    ID;
    int    epoch;
    double a;
    double e;
    double i;
    double LAN;
    double argperi;
    double M;
};

std::vector<Data> importCSV(const std::string &filename) {
    std::vector<Data> data;
    std::ifstream     file(filename);

    if (!file.is_open()) {
        std::cerr << "Could not open the file - '" << filename << "'" << std::endl;
        return data;
    }

    // ignore the header
    std::string line;
    std::getline(file, line);

    // read data line by line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        Data              temp{};
        ss >> temp.ID >> temp.epoch >> temp.a >> temp.e >> temp.i >> temp.LAN >> temp.argperi >> temp.M;
        data.push_back(temp);
    }

    file.close();
    return data;
}

#include <odin/domain/astrodynamics.hpp>

int main() {

    std::vector<Data> data = importCSV("examples/gtoc12/data/GTOC12_Asteroids_Data.txt");

    for (auto &d: data) {
        std::cout << d.ID << " " << d.epoch << " " << d.a << " " << d.e << " " << d.i << " " << d.LAN << " " << d.argperi << " " << d.M << std::endl;
    }
    // output the asteroid locations as position velocity
    //    template<typename Scalar>
    //    constexpr std::tuple<Eigen::Vector3<Scalar>, Eigen::Vector3<Scalar>> coe2rv(
    //            Scalar p,    // semi-latus rectum
    //            Scalar a,    // semi-major axis
    //            Scalar e,    // eccentricity
    //            Scalar i,    // inclination
    //            Scalar Omega,// right ascension of ascending node
    //            Scalar omega,// argument of periapsis
    //            Scalar nu,   // true anomaly
    //            Scalar mu    // gravitational parameter
    //    ) {

    std::ofstream outfile("output.csv");

    // Writing the header
    outfile << "id,r_x,r_y,r_z,v_x,v_y,v_z\n";

    for (auto &d: data) {
        auto [r, v] = coe2rv<double>(d.a * (1 - d.e * d.e), d.a, d.e, d.i, d.LAN, d.argperi, d.M, 1.32712440018e11);
        outfile << d.ID << "," << r[0] << "," << r[1] << "," << r[2] << "," << v[0] << "," << v[1] << "," << v[2] << "\n";
    }

    outfile.close();
    return 0;
}