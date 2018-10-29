#include <vector>

using namespace std;

class GaussPoints
{
public:
	vector<double> Points;
	vector<double> Weights;

	void GP(int nPoint)
	{
		switch (nPoint) {
		case 1:
			Points.resize(1);
			Points = { 0.0 };
			Weights.resize(1);
			Weights = { 2.0 };
			break;
		case 2:
			Points.resize(2);
			Points = { -sqrt(3.0) / 3.0 , +sqrt(3.0) / 3.0 };
			Weights.resize(2);
			Weights = { 1.0, 1.0 };
			break;
		case 3:
			Points.resize(3);
			Points = { -sqrt(15.0) / 5.0, 0.0, +sqrt(15.0) / 5.0 };
			Weights.resize(3);
			Weights = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
			break;
		case 4:
			Points.resize(4);
			Points = { -sqrt(525.0 + 70.0 * sqrt(30.0)) / 35.0, -sqrt(525.0 - 70.0 * sqrt(30.0)) / 35.0, +sqrt(525.0 + 70.0 * sqrt(30.0)) / 35.0, +sqrt(525.0 - 70.0 * sqrt(30.0)) / 35.0 };
			Weights.resize(4);
			Weights = { (18.0 - sqrt(30.0)) / 36.0 , (18.0 + sqrt(30.0)) / 36.0 , (18.0 + sqrt(30.0)) / 36.0 , (18.0 - sqrt(30.0)) / 36.0 };
			break;
		case 5:
			Points.resize(5);
			Points = { -sqrt(245.0 + 14.0 * sqrt(70.0)) / 21.0 , -sqrt(245.0 - 14.0 * sqrt(70.0)) / 21.0 , 0.0 , +sqrt(245.0 - 14.0 * sqrt(70.0)) / 21.0 , +sqrt(245.0 + 14.0 * sqrt(70.0)) / 21.0 };
			Weights.resize(5);
			Weights = { (322.0 - 13.0 * sqrt(70.0)) / 900.0, (322.0 + 13.0 * sqrt(70.0)) / 900.0 , 128.0 / 225.0, (322.0 + 13.0 * sqrt(70.0)) / 900.0 , (322.0 - 13.0 * sqrt(70.0)) / 900.0 };
			break;
		}
	}
};