//
// Group Project 2/ Ex3.3
//

# include  <iostream>
# include <fstream>
# include <chrono>
# include <string.h>
# include <vector>
# include <cmath>

int main(){

    // Declaring and initializing inputs
    double resolution = 32;
    double iterations = 100;
    double res_plus_ghosts = resolution + 2;
    double h = 1 / (resolution - 1);
    double pi = M_PI;

    //Starting value for u_h
    // TODO: Change starting value to Zero as in handout stated -> for debugging reasons it is 1, so I can see it.
    double u_h = 0;

    // Creating the grid in matrix form based on the resolution, plus 2 points in each direction for the ghost layers
    std::vector<std::vector<double> > grid (res_plus_ghosts, std::vector <double> (res_plus_ghosts, 0));

    // Filling the boundary conditions of the grid
    // BC (x, 0) = 0, Bottom
    for (int col = 1; col <= res_plus_ghosts - 2; col++)
    {
        grid[1][col] = 0;
    }
    // BC (0, y) = 0, Left
    for (int row = 1; row <= res_plus_ghosts - 2; row++)
    {
        grid[row][1] = 0;
    }
    // BC (1, y) = 0, Right
    for (int row = 1; row <= res_plus_ghosts - 2; row++)
    {
        grid[row][res_plus_ghosts - 2] = 0;
    }

    // BC (x, 1) = 0, Top
    // TODO: Check if last value is alright. Seems off.
    for (int col = 1; col <= res_plus_ghosts - 2; col++)
    {
        double x_coord =  h * (col - 1); // x-value moving along the grid
        grid[res_plus_ghosts - 2][col] = sin(2 * pi * x_coord) * sinh(2 * pi);
    }

    // Filling in the starting values for the vector entries to solve
    for (int row = 2; row <= res_plus_ghosts - 3; row++)
    {
        for (int col = 2; col <= res_plus_ghosts - 3; col++)
        {
            grid[row][col] = u_h;
        }
    }

    // Creating a grid where the updated values are stored
    std::vector<std::vector<double> > grid_new = grid;

    // Elements of the diagonal
    double d_elements = 4 + 4 * pi * pi * h * h;







    // Jacobi Method
    for (int iter = 1; iter <= iterations; iter++)
    {
        // Loop over boundary nodes

        // BC (x, 0) = 0, Bottom
        for (int col = 2; col <= res_plus_ghosts - 3; col++)
        {
            double row = 2;
            double x_coord =  h * (col - 1); // x-value moving along the grid in respect to the loop
            double y_coord =  h * (row - 1); // y-value moving along the grid in respect to the loop


            double u_s = grid[row-1][col];
            double u_c = grid[row][col];
            double b_rhs = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord) + u_s;


            // Jacobi iteration scheme
            grid_new[row][col] = (b_rhs - (( - 1 - 1 - 1)  * u_c * 1/(h * h))) * 1 / d_elements;

        }
        // BC (0, y) = 0, Left
        for (int row = 2; row <= res_plus_ghosts - 3; row++)
        {
            double col = 2;
            double x_coord =  h * (col - 1); // x-value moving along the grid in respect to the loop
            double y_coord =  h * (row - 1); // y-value moving along the grid in respect to the loop


            double u_w = grid[row][col-1];
            double u_c = grid[row][col];
            double b_rhs = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord) + u_w;


            // Jacobi iteration scheme
            grid_new[row][col] = (b_rhs - (( - 1 - 1 - 1)  * u_c * 1/(h * h))) * 1 / d_elements;

        }
        // BC (1, y) = 0, Right
        for (int row = 2; row <= res_plus_ghosts - 3; row++)
        {
            double col = 2;
            double x_coord =  h * (col - 1); // x-value moving along the grid in respect to the loop
            double y_coord =  h * (row - 1); // y-value moving along the grid in respect to the loop


            double u_e = grid[row][col+1];
            double u_c = grid[row][col];
            double b_rhs = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord) + u_e;


            // Jacobi iteration scheme
            grid_new[row][col] = (b_rhs - (( - 1 - 1 - 1)  * u_c * 1/(h * h))) * 1 / d_elements;

        }

        // BC (x, 1) = 0, Top
        for (int col = 2; col <= res_plus_ghosts - 3; col++)
        {
            double row = 2;
            double x_coord =  h * (col - 1); // x-value moving along the grid in respect to the loop
            double y_coord =  h * (row - 1); // y-value moving along the grid in respect to the loop


            double u_n = grid[1+row][col];
            double u_c = grid[row][col];
            double b_rhs = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord) + u_n;


            // Jacobi iteration scheme
            grid_new[row][col] = (b_rhs - (( - 1 - 1 - 1)  * u_c * 1/(h * h))) * 1 / d_elements;

        }

        // Loop over interior nodes
        for (int row = 3; row <= res_plus_ghosts - 4; row++)
        {
            for (int col = 3; col <= res_plus_ghosts - 4; col++)
            {
                double x_coord =  h * (col - 1); // x-value moving along the grid in respect to the loop
                double y_coord =  h * (row - 1); // y-value moving along the grid in respect to the loop


//                double u_n = grid[1+row][col];
//                double u_s = grid[row-1][col];
//                double u_e = grid[row][col+1];
//                double u_w = grid[row][col-1];
                double u_c = grid[row][col];
                double b_rhs = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord);


                // Jacobi iteration scheme
                grid_new[row][col] = (b_rhs - ((-1 - 1 - 1 - 1)  * u_c * 1/(h * h))) * 1 / d_elements;


            }
        }
        // Updating the old grid with the new values
        grid = grid_new;
    }





    std::cout << h << std::endl;








    return 0;
};