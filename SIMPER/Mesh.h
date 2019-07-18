#ifndef __MESH__
#define __MESH__

#include <Eigen/Core>
#include <vector>

using namespace std;
using namespace Eigen;

struct Coords
{
	double x, y, z;
};

struct Node
{
	int n;
	bool IsBoundary = false;
	Coords Coordinates;
};

struct Element
{
	int e;
	vector<Node> Nodes;
	double Area;
	string SoilType;
	double SoilHydraulicConductivity;
	double SoilHeatCapacity;
	double SoilThermalConductivity;
	double SoilDensity;
	double SoilFreezingPoint;
};

class Mesh
{
public:
	vector<Element> Elements;
	vector<Node> Nodes;
	size_t NumberOfElements;
	size_t NumberOfNodes;
	int ElementNumberOfNodes;
	double MaxX;
	double MaxY;
	double MaxZ;
	double MinX;
	double MinY;
	double MinZ;

	Mesh();
	VectorXd GetNodesXCoordinates(int e, int ndoe);
	VectorXd GetNodesYCoordinates(int e, int ndoe);
	MatrixXd GetNodesXYCoordinates(int e, int ndoe);
	VectorXi GetElementDofs(int e, int ndoe);
	VectorXd GetNodalValues(VectorXd vector, VectorXi indices);
	double GetElementArea(int e, int ndoe);
};

#endif