#include "Mesh.h"

Mesh::Mesh() 
{

}

VectorXd Mesh::GetNodesXCoordinates(int e, int ndoe)
{
	VectorXd xNodes(ndoe);
	xNodes <<
		Elements[e].Nodes[0].Coordinates.x,
		Elements[e].Nodes[1].Coordinates.x,
		Elements[e].Nodes[2].Coordinates.x,
		Elements[e].Nodes[3].Coordinates.x;


	return xNodes;
}

VectorXd Mesh::GetNodesYCoordinates(int e, int ndoe)
{
	VectorXd yNodes(ndoe);
	yNodes <<
		Elements[e].Nodes[0].Coordinates.y,
		Elements[e].Nodes[1].Coordinates.y,
		Elements[e].Nodes[2].Coordinates.y,
		Elements[e].Nodes[3].Coordinates.y;

	return yNodes;
}

MatrixXd Mesh::GetNodesXYCoordinates(int e, int ndoe)
{
	MatrixXd xyNodes(ndoe, 2);
	VectorXd xNodes(ndoe);
	VectorXd yNodes(ndoe);
	xNodes <<
		Elements[e].Nodes[0].Coordinates.x,
		Elements[e].Nodes[1].Coordinates.x,
		Elements[e].Nodes[2].Coordinates.x,
		Elements[e].Nodes[3].Coordinates.x;
	yNodes <<
		Elements[e].Nodes[0].Coordinates.y,
		Elements[e].Nodes[1].Coordinates.y,
		Elements[e].Nodes[2].Coordinates.y,
		Elements[e].Nodes[3].Coordinates.y;
	xyNodes << xNodes, yNodes;

	return xyNodes;
}


VectorXi Mesh::GetElementDofs(int e, int ndoe)
{
	VectorXi elementDofs(ndoe), one(ndoe);
	elementDofs << Elements[e].Nodes[0].n, Elements[e].Nodes[1].n, Elements[e].Nodes[2].n, Elements[e].Nodes[3].n;

	return elementDofs - one.Ones(ndoe);
}

VectorXd Mesh::GetNodalValues(VectorXd vector, VectorXi indices)
{
	VectorXd nodalValues(indices.size());
	for (int n = 0; n < indices.size(); n++)
	{
		nodalValues(n) = vector(indices(n));
	}

	return nodalValues;
}

double Mesh::GetElementArea(int e, int ndoe)
{
	VectorXd xNodes(ndoe + 1);
	VectorXd yNodes(ndoe + 1);

	xNodes.head(ndoe) = GetNodesXCoordinates(e, ndoe);
	yNodes.head(ndoe) = GetNodesYCoordinates(e, ndoe);
	xNodes(ndoe) = xNodes(0);
	yNodes(ndoe) = yNodes(0);

	double area = 0.0;
	for (int i = 0; i < ndoe; i++)
	{
		area += xNodes(i)*yNodes(i + 1) - yNodes(i)*xNodes(i + 1);
	}
	
	area = 0.5 * abs(area);

	return area;
}
