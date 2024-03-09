#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <list>
#include <random>
#include <numeric>
#include <algorithm>
#include <fstream>

const double k = 1.380649E-23;

const double pi = 3.141592653589793238462643383279502884197;

const double G = 6.67430E-11;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class OctNode {
public:
    std::vector<std::vector<long double>> Points;
    std::vector<long double> Masses;
    std::vector<long double> Center = {0, 0, 0};
    long double Radius;
    std::vector<int> IDs;
    std::vector<OctNode> Leaves;
    std::vector<long double> g;
    std::vector<OctNode> Children = {};
    long double Mass = 0;
    std::vector<long double> COM = {0, 0 ,0};
    int ID;

    OctNode(std::vector<long double> center, long double radius, std::vector<std::vector<long double>> points, std::vector<long double> masses, std::vector<int> ids, std::vector<OctNode> leaves)
    { //xyz coord, cube half-length in cm, coords for particles in node, masses for particles in node, node id, child node id's
        Points = std::move(points);
        Masses = std::move(masses);
        Center = std::move(center);
        Radius = radius;
        IDs = std::move(ids);
        Leaves = std::move(leaves);

        if (Points.size() == 1)
        { //octants with only 1 particle become nodes
            Leaves.push_back(*this);
            COM = Points[0];
            Mass = Masses[0];
            g = {0, 0, 0};
            ID = IDs[0];
        }

        else
        { //octants with 2+ particles are subdivided into further nodes
            GenerateChildren(Points, Masses, IDs, Leaves);
            long double comTotal[3];
            long double mTotal;
            mTotal = 0;
            for (OctNode c: Children) {
                long double m = Masses[0];
                std::vector<long double> com = Points[0];
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

    void GenerateChildren(std::vector<std::vector<long double>>& genpoints, std::vector<long double>& genmasses, std::vector<int>& genids, std::vector<OctNode>& leaves)
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

            if (std::accumulate(inOctant.begin(), inOctant.end(), 0) != 0)
            {
                long double dx = Radius*(1-octant[0]);
                long double dy = Radius*(1-octant[1]);
                long double dz = Radius*(1-octant[2]);

                std::vector<long double> newCenter = {0, 0, 0};
                newCenter[0] = Center[0] + dx;
                newCenter[1] = Center[1] + dy;
                newCenter[2] = Center[2] + dz;

                OctNode newNode = OctNode(newCenter, Radius/2, genChildrenPoints, genChildrenMasses, genChildrenIDS, leaves);
                Children.push_back(newNode);
            }

            else
            {
                continue;
            }
        }
    }
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void TreeWalk(OctNode& node, OctNode& node0, double theta)
{
    std::vector<long double> dx = {node.COM[0] - node0.COM[0], node.COM[1] - node0.COM[1], node.COM[2] - node0.COM[2]};
    long double r = std::sqrt(std::pow(node.COM[0] - node0.COM[0], 2) + std::pow(node.COM[1] - node0.COM[1], 2) + std::pow(node.COM[2] - node0.COM[2], 2));
    if (r > 0)
    {
        if (node.Children.empty() or 2 * node.Radius / r < theta)
        {
            node0.g[0] += G * node.Mass * dx[0] / std::pow(r, 3);
            node0.g[1] += G * node.Mass * dx[1] / std::pow(r, 3);
            node0.g[2] += G * node.Mass * dx[2] / std::pow(r, 3);
        }
        else
        {
            for (OctNode& c: node.Children) {
                TreeWalk(c, node0, theta);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::vector<long double>> GravAccel(const std::vector<std::vector<long double>>& Points, std::vector<long double> Masses, double theta)
{
    std::vector<long double> pointsx;
    std::vector<long double> pointsy;
    std::vector<long double> pointsz;
    std::vector<int> IDS;
    for (int i = 0 ; std::vector<long double> point: Points) {
        pointsx.push_back(point[0]);
        pointsy.push_back(point[1]);
        pointsz.push_back(point[2]);
        IDS.push_back(i);
        i++;
    }

    std::vector<long double> center = {
            (*std::max_element(pointsx.begin(), pointsx.end()) + *std::min_element(pointsx.begin(), pointsx.end())) / 2,
            (*std::max_element(pointsy.begin(), pointsy.end()) + *std::min_element(pointsy.begin(), pointsy.end())) / 2,
            (*std::max_element(pointsz.begin(), pointsz.end()) + *std::min_element(pointsz.begin(), pointsz.end())) /
            2};

    long double topsize = std::max(
            {(*std::max_element(pointsx.begin(), pointsx.end()) - *std::min_element(pointsx.begin(), pointsx.end())),
             (*std::max_element(pointsy.begin(), pointsy.end()) - *std::min_element(pointsy.begin(), pointsy.end())),
             (*std::max_element(pointsz.begin(), pointsz.end()) - *std::min_element(pointsz.begin(), pointsz.end()))});

    std::vector<OctNode> Leaves = {};

    OctNode topNode = OctNode(center, topsize / 2, Points, std::move(Masses), IDS, Leaves);

    std::vector<std::vector<long double>> accel;
    for (OctNode leaf: topNode.Leaves) {
        TreeWalk(topNode, leaf, theta);
        accel[leaf.ID][0] = leaf.g[0];
        accel[leaf.ID][1] = leaf.g[1];
        accel[leaf.ID][2] = leaf.g[2];
    }

    return accel;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::tuple<std::vector<std::vector<long double>>, std::vector<std::vector<long double>>, std::vector<long double>> initialConditions()
{ //Function that creates the initial volume/top-node of the simulation
    double boxSize;
    long long int n;
    double T;
    std::cout << "Enter box length, number of particles, temperature: " << std::endl;
    std::cin >> boxSize >> n >> T;
    std::cout << "box length = " << boxSize << std::endl
              << "number of particles = " << n << std::endl
              << "temperature = " << T << std::endl;

    std::vector<std::vector<long double>> points0 = {}; //This creates a vector of points that will represent the locations of our particles
    for (int i = 0; i < n; i++)
    {
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(-boxSize / 2, boxSize / 2);
        std::vector<long double> tempPoint = {dis(gen), dis(gen), dis(gen)};
        points0.push_back(tempPoint);
    }

    std::vector<long double> masses0 = {};
    std::vector<std::vector<long double>> velocities0 = {};
    for (int i2 = 0; i2 < n; i2++)
    {
        std::random_device rd2;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen2(rd2()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis2(1.6735575E-27, 9.2732796E-26);
        double tempMass = dis2(gen2);
        masses0.push_back(tempMass);
        std::vector<long double> tempMomentum = {tempMass , 0, 0, 0};
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::normal_distribution<> dis(std::sqrt(8*k*T/(pi*tempMass)), std::sqrt(3*k*T/tempMass));
        std::vector<long double> tempVelocity = {dis(gen), dis(gen), dis(gen)};
        tempMomentum[1] = tempVelocity[0];
        tempMomentum[2] = tempVelocity[1];
        tempMomentum[3] = tempVelocity[2];
        velocities0.push_back(tempMomentum);
    }

    return {points0, velocities0, masses0};

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Advance(const double Duration, const double dt, std::vector<std::vector<long double>>& Points, std::vector<std::vector<long double>>& Velocities, std::vector<long double>& Masses)
{
    double time=0;
    while (time < Duration)
    {
        std::vector<std::vector<long double>> a = GravAccel(Points, Masses, 0.5);

        std::string s = R"(C:\Users\Andym\OneDrive\Documents\GravSimData\C++\galaxy_t)" + std::to_string(time) + ".csv" ;
        std::ofstream myFile;
        myFile.open(s);

        for (int i = 0; i < Velocities.size(); i++)
        {
            Velocities[i][1] += a[i][0]*dt;
            Velocities[i][2] += a[i][1]*dt;
            Velocities[i][3] += a[i][2]*dt;

            Points[i][0] += Velocities[i][1]*dt;
            Points[i][1] += Velocities[i][2]*dt;
            Points[i][2] += Velocities[i][3]*dt;

            myFile << Points[i][0] << "," << Points[i][1] << "," << Points[i][2] << ",\n";
        }
        myFile.close();
        time += dt;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{
    std::vector<std::vector<long double>> Points;
    std::vector<std::vector<long double>> Velocities;
    std::vector<long double> Masses;
    std::tuple<std::vector<std::vector<long double>>, std::vector<std::vector<long double>>, std::vector<long double>> init = initialConditions();
    Points = get<0>(init);
    Velocities = get<1>(init);
    Masses = get<2>(init);
    Advance(10, 1, Points, Velocities, Masses);
    return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
