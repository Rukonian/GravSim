#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <list>
#include <random>
#include <numeric>
#include <algorithm>
#include <fstream>

//const double R = 8.3145;

const double k = 1.380649E-23;

const double pi = 3.141592653589793238462643383279502884197;

const double G = 6.67430E-11;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//struct Point {
//	double X;
//	double Y;
//	double Z;
//
//	Point(double x, double y, double z)
//    {
//		X = x;
//		Y = y;
//		Z = z;
//	}
//	Point() = default;
//}; //Defines the Position vector data structure
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//struct Momentum {
//	double M;
//	double Vx;
//	double Vy;
//	double Vz;
//
//	Momentum(double m, double vx, double vy, double vz)
//    {
//		M = m;
//		Vx = vx;
//		Vy = vy;
//		Vz = vz;
//	}
//	Momentum() = default;
//}; //Defines the Momentum vector data structure
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//struct Force {
//	double M;
//	double Ax;
//	double Ay;
//	double Az;
//
//	Force(double m, double ax, double ay, double az)
//    {
//		M = m;
//		Ax = ax;
//		Ay = ay;
//		Az = az;
//	}
//	Force() = default;
//}; //Defines the Force vector data structure
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

        //std::cout << "confirm octnode creation " << 1 << std::endl; /////////////////////////////////////

        if (Points.size() == 1)
        { //octants with only 1 particle become nodes
            std::cout << "leaves push_back pre-check " << Center[0] << ", " << Center[1] << ", " << Center[2] << std::endl; /////////////////////////////////////
            Leaves.push_back(*this);
            std::cout << "leaves push_back post-check " << Leaves.back().Center[0] << ", " << Leaves.back().Center[1] << ", " << Leaves.back().Center[2] << std::endl; /////////////////////////////////////
            COM = Points[0];
            Mass = Masses[0];
            g = {0, 0, 0};
            ID = IDs[0];
        }

        else
        { //octants with 2+ particles are subdivided into further nodes
            //std::cout << "generating children" << std::endl; /////////////////////////////////////
            GenerateChildren(Points, Masses, IDs, Leaves);
            //std::cout << "confirm generatechildren worked " << 1 << std::endl; /////////////////////////////////////
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
            //std::cout << "octantindex " << (point[0]>Center[0]) << (point[1]>Center[1]) << (point[2]>Center[2]) << std::endl;
        }
        //std::cout << "confirm octantindex creation " << 1 << std::endl; /////////////////////////////////////
        //std::cout << "octindex, genpoints, genmasses, genids size: " << octantIndex.size() << ", " << genpoints.size() << ", " << genmasses.size() << ", " << genids.size() << std::endl;

        std::vector<std::vector<int>> octants = {{1,1,1}, {0,1,1}, {0,0,1}, {1,0,1}, {1,1,0}, {0,1,0}, {0,0,0}, {1,0,0}};
        for (auto const &octant: octants)
        {
            //std::cout << "octant " << octant[0] << octant[1] << octant[2] << std::endl;
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
            //std::cout << "confirm while-loop finished " << 1 << std::endl; /////////////////////////////////////

            if (std::accumulate(inOctant.begin(), inOctant.end(), 0) != 0)
            {
                //std::cout << "confirm if-statement start " << 1 << std::endl; /////////////////////////////////////
                long double dx = Radius*(1-octant[0]);
                long double dy = Radius*(1-octant[1]);
                long double dz = Radius*(1-octant[2]);
                //std::cout << "confirm differentials " << 1 << std::endl; /////////////////////////////////////

                std::vector<long double> newCenter = {0, 0, 0};
                //std::cout << "Center " << Center[0] << Center[1] << Center[2] << std::endl;
                newCenter[0] = Center[0] + dx;
                newCenter[1] = Center[1] + dy;
                newCenter[2] = Center[2] + dz;
                //std::cout << "newCenter " << newCenter[0] << ", " << newCenter[1] << ", " << newCenter[2] << std::endl;
                //std::cout << "confirm newcenter " << 1 << std::endl; /////////////////////////////////////

                OctNode newNode = OctNode(newCenter, Radius/2, genChildrenPoints, genChildrenMasses, genChildrenIDS, leaves);
                Children.push_back(newNode);
                //std::cout << "confirm child node append finished " << 1 << std::endl; /////////////////////////////////////
            }

            else
            {
                //std::cout << "confirm else-statement " << 1 << std::endl; /////////////////////////////////////
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
            std::cout << "treewalk if-statement " << std::endl;
            node0.g[0] += G * node.Mass * dx[0] / std::pow(r, 3);
            node0.g[1] += G * node.Mass * dx[1] / std::pow(r, 3);
            node0.g[2] += G * node.Mass * dx[2] / std::pow(r, 3);
            std::cout << "node0.g check " << node0.g[0] << ", " << node0.g[1] << ", " << node0.g[2] << std::endl; /////////////////////////////////////
        }
        else
        {
            std::cout << "treewalk else-statement " << std::endl;
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
    //std::cout << Points.size() << ", " << Masses.size() << ", " << std::endl; //////////////////////////////////

    std::vector<long double> pointsx;
    std::vector<long double> pointsy;
    std::vector<long double> pointsz;
    std::vector<int> IDS;
    for (int i = 0 ; std::vector<long double> point: Points) {
        pointsx.push_back(point[0]);
        pointsy.push_back(point[1]);
        pointsz.push_back(point[2]);
        IDS.push_back(i);
        std::cout << "GravAccel print " << pointsx[i] << ", " << pointsy[i] << ", " << pointsz[i] << std::endl; /////////////////////////////////////
        i++;
    }

    std::vector<long double> center = {
            (*std::max_element(pointsx.begin(), pointsx.end()) + *std::min_element(pointsx.begin(), pointsx.end())) / 2,
            (*std::max_element(pointsy.begin(), pointsy.end()) + *std::min_element(pointsy.begin(), pointsy.end())) / 2,
            (*std::max_element(pointsz.begin(), pointsz.end()) + *std::min_element(pointsz.begin(), pointsz.end())) /
            2};
    //std::cout << "center print " << center[0] << ", " << center[1] << ", " << center[2] << std::endl; /////////////////////////////////////

    long double topsize = std::max(
            {(*std::max_element(pointsx.begin(), pointsx.end()) - *std::min_element(pointsx.begin(), pointsx.end())),
             (*std::max_element(pointsy.begin(), pointsy.end()) - *std::min_element(pointsy.begin(), pointsy.end())),
             (*std::max_element(pointsz.begin(), pointsz.end()) - *std::min_element(pointsz.begin(), pointsz.end()))});
    //std::cout << "topsize print " << topsize << std::endl; /////////////////////////////////////

    std::vector<OctNode> Leaves = {};
    //std::cout << "confirm leaves creation " << 1 << std::endl; /////////////////////////////////////

    OctNode topNode = OctNode(center, topsize / 2, Points, std::move(Masses), IDS, Leaves);
    std::cout << "confirm topnode leaves size " << topNode.Leaves.size() << std::endl; /////////////////////////////////////

    //std::cout << "topNode Leaves size = " << topNode.Leaves.size() << std::endl; /////////////////////////////////////

    std::vector<std::vector<long double>> accel;
    for (OctNode leaf: topNode.Leaves) {
        TreeWalk(topNode, leaf, theta);
        std::cout << "leaf ID " << leaf.ID << std::endl; /////////////////////////////////////////
        accel[leaf.ID][0] = leaf.g[0];
        accel[leaf.ID][1] = leaf.g[1];
        accel[leaf.ID][2] = leaf.g[2];
        std::cout << leaf.g[0] << ", " << leaf.g[1] << ", " << leaf.g[2] << std::endl; /////////////////////////////////////////
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
    std::vector<int> IDs(n);
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

    std::iota (IDs.begin(), IDs.end(), 0);

    //for (Point i: points0) {                              //debug
    //    std::cout<<i.X<<", "<<i.Y<<", "<<i.Z<<std::endl;
    //}

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

    //for (double i: masses0) {                              //debug
    //    std::cout<<i<<std::endl;
    //}

    return {points0, velocities0, masses0};
    //std::cout<<topNode.Center.X<<", "<<topNode.Center.Y<<", "<<topNode.Center.Z<<std::endl;                              //debug
    //std::cout<<topNode.Radius<<std::endl;
    //std::cout<<topNode.ID<<std::endl;

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

        std::cout << "points, velocities, masses, accel sizes " << Points.size() << ", " << Velocities.size() << ", " << Masses.size() << ", " << a.size() << std::endl; /////////////////////////////////////////

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