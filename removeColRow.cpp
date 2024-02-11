
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#define _USE_MATH_DEFINES

using namespace std;
using namespace Eigen;

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove){
	unsigned int numRows = matrix.rows() - 1;
	unsigned int numCols = matrix.cols();

	if( rowToRemove < numRows )
		matrix.block(rowToRemove, 0, numRows-rowToRemove, numCols) = matrix.bottomRows(numRows-rowToRemove);

	matrix.conservativeResize(numRows);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove){
	unsigned int numRows = matrix.rows();
	unsigned int numCols = matrix.cols()-1;

	if( colToRemove < numCols )
		matrix.block(0, colToRemove, numRows, numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);

	matrix.conservativeResize(numRows,numCols);
}

int main()
{

MatrixXd test(4,1);
test << 1,
	5,
	9,
	13;

cout << "Before remove =\n" << test << endl;

removeRow(test, 1);

cout << "After remove =\n" << test << endl;

return 0;

}
