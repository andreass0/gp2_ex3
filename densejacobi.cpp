//
// Group Project 2/ Ex3.3
//

# include  <iostream>
# include <fstream>
# include <chrono>
# include <string.h>
# include <vector>
//# include <math.h>
# include "cmath"
# include <eigen-3.4.0/Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;


int main(){

    // Declaring and initializing inputs
    double resolution = 6;
    double iterations = 100;
    double res_plus_ghosts = resolution + 2;
//    double resol_double = resolution;
    double h = 1 / (resolution - 1);
    double pi = M_PI;
    int number_unknowns = pow(resolution, 2) - (4 * resolution - 4);

    // Creating the matrix and the needed arrays
    MatrixXd A_h(number_unknowns, number_unknowns);
    VectorXd u_h(number_unknowns);
    VectorXd b_h(number_unknowns);

    // Filling in the matrix and the arrays
    //A_h
    for (int row = 0; row < number_unknowns; row++)
    {
        for (int col = 0; col < number_unknowns; col++)
        {
            A_h(row, col) = 0;
        }
    }

    //u_h (filling in with the starting value for the iteration)
    for (int row = 0; row < number_unknowns; row++)
    {
        u_h(row) = 0;
    }

//    double x_coord =  h * (col - 1); // x-value moving along the grid in respect to the loop
//    double y_coord =  h * (row - 1); // y-value moving along the grid in respect to the loop

    //b_h

    for (int row = 0; row < number_unknowns; row++)
    {
        b_h(row) = 0;

    }

    // TODO: Finish up! Check if coordinates are alright. Fill in vector.
    // Boundaries for the corner points
    // Bottom
    // Bottom Left Corner
    double x_BL =  h;
    double y_BL =  h;
    b_h(0) = 4 * pi * pi * sin(2 * pi * x_BL) * sinh(2 * pi * y_BL);

    // Fill everything in between on the bottom
    for (int row = 1; row < 1 + resolution - 4; row++)
    {
        b_h(row) = 0;

    }

    // Bottom Right Corner
    double x_BR =  h * (number_unknowns - 1);
    double y_BR =  h;
    int row_BR = 1 + resolution - 4;
    b_h(row_BR) = 4 * pi * pi * sin(2 * pi * x_BR) * sinh(2 * pi * y_BR);

    // Fill everything for the left boundary

    //Top
    // Top Left Corner
    double x_TL =  h;
    double y_TL =  h * (number_unknowns - 1);
    int row_TL = number_unknowns - 1 - (1 + resolution - 4);
    b_h(row_TL) = 4 * pi * pi * sin(2 * pi * x_TL) * sinh(2 * pi * y_TL) +
            4 * pi * pi * sin(2 * pi * x_TL) * sinh(2 * pi);

    double x_coord_count = 2;
    // Fill everything in between on the top
    for (int row = number_unknowns - 1 - (1 + resolution - 4) + 1; row < number_unknowns - 1; row++)
    {

        double x_coord =  h * x_coord_count; // x-value moving along the grid in respect to the loop
        double y_coord =  1; // y-value moving along the grid in respect to the loop
        b_h(row) = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord) +
                   4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi);
        x_coord_count++;

    }

    // Top Right Corner
    double x_TR =  h * (number_unknowns - 1);
    double y_TR =  h * (number_unknowns - 1);
    int row_TR = number_unknowns - 1;
    b_h(row_TR) = 4 * pi * pi * sin(2 * pi * x_TR) * sinh(2 * pi * y_TR) +
            4 * pi * pi * sin(2 * pi * x_TL) * sinh(2 * pi);








    //Starting value for u_h
    // TODO: Change starting value to Zero as in handout stated -> for debugging reasons it is 1, so I can see it.
////    double u_h = 0;
//
//    // Creating the grid in matrix form based on the resolution, plus 2 points in each direction for the ghost layers
//    std::vector<std::vector<double> > grid (res_plus_ghosts, std::vector <double> (res_plus_ghosts, 0));
//
//    // Filling the boundary conditions of the grid
//    // BC (x, 0) = 0, Bottom
//    for (int col = 1; col <= res_plus_ghosts - 2; col++)
//    {
//        grid[1][col] = 0;
//    }
//    // BC (0, y) = 0, Left
//    for (int row = 1; row <= res_plus_ghosts - 2; row++)
//    {
//        grid[row][1] = 0;
//    }
//    // BC (1, y) = 0, Right
//    for (int row = 1; row <= res_plus_ghosts - 2; row++)
//    {
//        grid[row][res_plus_ghosts - 2] = 0;
//    }
//
//    // BC (x, 1) = 0, Top
//    // TODO: Check if last value is alright. Seems off.
//    for (int col = 1; col <= res_plus_ghosts - 2; col++)
//    {
//        double x_coord =  h * (col - 1); // x-value moving along the grid
//        grid[res_plus_ghosts - 2][col] = sin(2 * pi * x_coord) * sinh(2 * pi);
//    }
//
//    // Filling in the starting values for the vector entries to solve
//    for (int row = 2; row <= res_plus_ghosts - 3; row++)
//    {
//        for (int col = 2; col <= res_plus_ghosts - 3; col++)
//        {
//            grid[row][col] = u_h;
//        }
//    }
//
//    // Creating a grid where the updated values are stored
//    std::vector<std::vector<double> > grid_new = grid;
//
//    // Elements of the diagonal
//    double d_elements = 4 + 4 * pi * pi * h * h;
//
//
//
//
//
//
//
//
//
//    std::cout << h << std::endl;
//
//
//
//


    std::cout << A_h << std::endl;
    std::cout << b_h << std::endl;



    return 0;
};