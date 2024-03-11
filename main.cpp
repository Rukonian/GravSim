#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <random>
#include <numeric>
#include <algorithm>
#include <fstream>

const double k = 1.380649E-23; // N*m/K

const double pi = 3.141592653589793238462643383279502884197; // unit-less

const double G = 6.67430E-11; // N*m^2/kg^2

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class OctNode {
public:
    std::vector<std::vector<long double>> Points = {};
    std::vector<long double> Masses = {};
    long double Mass = 0;
    std::vector<long double> Center = {0, 0, 0};
    std::vector<long double> COM = {Center[0], Center[1], Center[2]};
    long double Radius = 0;
    std::vector<int> IDs = {};
    int ID;
    std::vector<OctNode> Leaves;
    std::vector<OctNode> Children = {};
    std::vector<long double> g = {};

    OctNode(std::vector<long double> & center, long double & radius, std::vector<std::vector<long double>> & points, std::vector<long double> & masses, std::vector<int> & ids, std::vector<OctNode> & leaves) {
        //xyz coord, cube half-length in cm, coords for particles in node, masses for particles in node, node id, child node id's
        Center = center;
        Radius = radius;
        Points = points;
        Masses = masses;
        IDs = ids;
        Leaves = leaves;

        if (Points.size() == 1) { //octants with only 1 particle become nodes
            std::cout << "OctNode If-Statement start" << std::endl; /////////////////////////////////////////////////////////////
            Leaves.push_back(*this);
            std::cout << "Leaves.push_back check " << Leaves.size() << std::endl; ////////////////////////////////////////////////////////
            COM = Points[0];
            Mass = Masses[0];
            g = {0, 0, 0};
            ID = IDs[0];
            std::cout << "OctNode If-Statement end" << std::endl; /////////////////////////////////////////////////////////////
        }

        else{ //octants with 2+ particles are subdivided into further nodes
            std::cout << "Generating Children" << std::endl; ///////////////////////////////////////////////////////////////
            GenerateChildren(this->Points, this->Masses, this->IDs, this->Leaves);
            long double comTotal[3] = {0, 0, 0};
            long double mTotal = 0;
            for (const OctNode &c: Children) {
                long double m = c.Mass;
                std::vector<long double> com = c.COM;
                mTotal += m;
                comTotal[0] += m * com[0];
                comTotal[1] += m * com[1];
                comTotal[2] += m * com[2];
            }
            Mass = mTotal;
            COM[0] = comTotal[0]/Mass;
            COM[1] = comTotal[1]/Mass;
            COM[2] = comTotal[2]/Mass;
        }
    }

    void GenerateChildren(std::vector<std::vector<long double>> & genpoints, std::vector<long double> & genmasses, std::vector<int> & genids, std::vector<OctNode> & leaves)
    {
        std::vector<std::vector<int>> octantIndex = {};
        for (auto const &point: genpoints)
        {
            octantIndex.push_back({point[0]>Center[0], point[1]>Center[1], point[2]>Center[2]});
        }

        std::vector<std::vector<int>> octants = {{1,1,1}, {0,1,1}, {0,0,1}, {1,0,1}, {1,1,0}, {0,1,0}, {0,0,0}, {1,0,0}};
        for (auto const &octant: octants)
        {
            int i =0;
            std::vector<int> inOctant(octantIndex.size(), 0);
            std::vector<std::vector<long double>> genChildrenPoints = {};
            std::vector<long double> genChildrenMasses = {};
            std::vector<int> genChildrenIDS = {};
            while (i < octantIndex.size())
            {
                if (octantIndex[i] == octant)
                {
                    inOctant[i] += 1;
                    genChildrenPoints.push_back(genpoints[i]);
                    genChildrenMasses.push_back(genmasses[i]);
                    genChildrenIDS.push_back(genids[i]);
                }
                i += 1;
            }

            if (std::accumulate(inOctant.begin(), inOctant.end(), 0) != 0) {
                long double dx = Radius*(1-octant[0]);
                long double dy = Radius*(1-octant[1]);
                long double dz = Radius*(1-octant[2]);

                std::vector<long double> newCenter = {0, 0, 0};
                newCenter[0] = Center[0] + dx;
                newCenter[1] = Center[1] + dy;
                newCenter[2] = Center[2] + dz;
                long double newRadius = Radius/2;

                OctNode newNode = OctNode(newCenter, newRadius, genChildrenPoints, genChildrenMasses, genChildrenIDS, leaves);
                this->Children.push_back(newNode);
            }

            else {
                continue;
            }
        }
    }
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::tuple<std::vector<std::vector<long double>>, std::vector<std::vector<long double>>, std::vector<long double>, long double> initialConditions()
{ //Function that creates the initial volume/top-node of the simulation
    double boxSize;
    long long int n;
    double T;
    double theta;
    std::cout << "Enter box length, number of particles, temperature, theta: " << std::endl;
    std::cin >> boxSize >> n >> T >> theta;
    std::cout << "box length = " << boxSize << std::endl
              << "number of particles = " << n << std::endl
              << "temperature = " << T << std::endl
              << "theta = " << theta << std::endl;

    std::vector<std::vector<long double>> points0 = {}; //This creates a vector of points that will represent the locations of our particles
    for (int i = 0; i < n; i++) {
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(-boxSize / 2, boxSize / 2);
        std::vector<long double> tempPoint = {dis(gen), dis(gen), dis(gen)};
        points0.push_back(tempPoint);
    }

    std::vector<long double> masses0 = {};
    std::vector<std::vector<long double>> velocities0 = {};
    for (int i2 = 0; i2 < n; i2++) {
        std::random_device rd2;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen2(rd2()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis2(1.6735575E-27, 9.2732796E-26);
        double tempMass = dis2(gen2);
        masses0.push_back(tempMass);
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::normal_distribution<> dis(std::sqrt(8*k*T/(pi*tempMass)), std::sqrt(3*k*T/tempMass));
        std::vector<long double> tempVelocity = {dis(gen), dis(gen), dis(gen)};
        velocities0.push_back(tempVelocity);
    }

    return {points0, velocities0, masses0, theta};
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void TreeWalk(OctNode & node, OctNode & node0, const long double & theta) {
    std::vector<long double> dx = {node.COM[0] - node0.COM[0], node.COM[1] - node0.COM[1], node.COM[2] - node0.COM[2]};
    long double r = std::sqrt(std::pow(node.COM[0] - node0.COM[0], 2) + std::pow(node.COM[1] - node0.COM[1], 2) + std::pow(node.COM[2] - node0.COM[2], 2));
    if (r > 0)
    {
        if (node.Children.empty() or 2 * node.Radius / r < theta) {
            node0.g[0] += G * node.Mass * dx[0] / std::pow(r, 3);
            node0.g[1] += G * node.Mass * dx[1] / std::pow(r, 3);
            node0.g[2] += G * node.Mass * dx[2] / std::pow(r, 3);
        }
        else {
            for (OctNode &c: node.Children) {
                TreeWalk(c, node0, theta);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::vector<long double>> GravAccel(std::vector<std::vector<long double>> & Points, std::vector<long double> & Masses, std::vector<int> & IDs, const long double & theta) {
    std::cout << "GravAccel Started" << std::endl; //////////////////////////////////////////////////////////////////////////////////
    std::vector<long double> pointsx;
    std::vector<long double> pointsy;
    std::vector<long double> pointsz;
    for (int i = 0 ; std::vector<long double> point: Points) {
        pointsx.push_back(point[0]);
        pointsy.push_back(point[1]);
        pointsz.push_back(point[2]);
        IDs.push_back(i);
        i += 1;
    }
    std::cout << "GravAccel For-Loop Finished" << std::endl; /////////////////////////////////////////////////////////////////////////
    std::cout << "Points-IDs size check " << Points.size() << ", " << IDs.size() << std::endl; /////////////////////////////////////////////////////////////////////////

    std::vector<long double> center = {
            (*std::max_element(pointsx.begin(), pointsx.end()) + *std::min_element(pointsx.begin(), pointsx.end())) / 2,
            (*std::max_element(pointsy.begin(), pointsy.end()) + *std::min_element(pointsy.begin(), pointsy.end())) / 2,
            (*std::max_element(pointsz.begin(), pointsz.end()) + *std::min_element(pointsz.begin(), pointsz.end())) /
            2};

    long double topsize = std::max(
            {(*std::max_element(pointsx.begin(), pointsx.end()) - *std::min_element(pointsx.begin(), pointsx.end())),
             (*std::max_element(pointsy.begin(), pointsy.end()) - *std::min_element(pointsy.begin(), pointsy.end())),
             (*std::max_element(pointsz.begin(), pointsz.end()) - *std::min_element(pointsz.begin(), pointsz.end()))});
    long double radius = topsize/2;
    std::cout << "GravAccel center, topsize, and radius determined" << std::endl; //////////////////////////////////////////////////////////////////////

    std::vector<OctNode> Leaves = {};

    OctNode topNode = OctNode(center, radius, Points, Masses, IDs, Leaves);
    std::cout << "topNode created" << std::endl; /////////////////////////////////////////////////////////////////
    std::cout << "topNode Leaves size check " << Leaves.size() << std::endl; ////////////////////////////////////////////

    std::vector<std::vector<long double>> a (Leaves.size(), {0, 0, 0});
    for (OctNode &leaf: Leaves) {
        TreeWalk(topNode, leaf, theta);
        a[leaf.ID] = {leaf.g[0], leaf.g[1], leaf.g[2]};
    }

    return a;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Advance(const long double Duration, const long double dt, std::vector<std::vector<long double>> &Points, std::vector<std::vector<long double>> &Velocities, std::vector<long double> &Masses, long double &theta) {
    std::cout << "Advance Called" << std::endl; /////////////////////////////////////////////////////////////
    long double time = 0;
    while (time < Duration) {
        std::vector<int> IDs = {};
        std::vector<std::vector<long double>> a = GravAccel(Points, Masses, IDs, theta);
        std::cout << "GravAccel Finished" << std::endl; ////////////////////////////////////////////////////////////////
        std::cout << "acceleration vector check " << a.size() << std::endl; ///////////////////////////////////////////////////////////////

        std::string s = R"(C:\Users\Andym\OneDrive\Documents\GravSimData\C++\galaxy_t)" + std::to_string(time/dt) + ".csv" ;
        std::ofstream myFile;
        myFile.open(s);

        std::cout << "Points-Velocities-Acceleration size check " << Points.size() << ", " << Velocities.size() << ", " << a.size() << std::endl; ////////////////////////////////////////////////////////////////////////////
        for (int i = 0; i < Velocities.size(); i += 1) {
            Velocities[i][0] += a[i][0]*dt;
            Velocities[i][1] += a[i][1]*dt;
            Velocities[i][2] += a[i][2]*dt;

            Points[i][0] += Velocities[i][0]*dt;
            Points[i][1] += Velocities[i][1]*dt;
            Points[i][2] += Velocities[i][2]*dt;

            myFile << Points[i][0] << "," << Points[i][1] << "," << Points[i][2] << "\n";
        }
        myFile.close();
        time += dt;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
    std::vector<std::vector<long double>> Points;
    std::vector<std::vector<long double>> Velocities;
    std::vector<long double> Masses;
    long double theta;
    std::tuple<std::vector<std::vector<long double>>, std::vector<std::vector<long double>>, std::vector<long double>, long double> init = initialConditions();
    Points = get<0>(init);
    Velocities = get<1>(init);
    Masses = get<2>(init);
    theta = get<3>(init);
    Advance(10, 1, Points, Velocities, Masses, theta);
    return 0;
}
