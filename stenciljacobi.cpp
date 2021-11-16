//
// Group Project 2/ Ex3.3
//

# include  <iostream>
# include <fstream>
# include <chrono>
# include <string.h>
# include <vector>
# include <cmath>
# include <sstream>
# include <assert.h>


//Reusing matrix-vector-multiplication function from task3
std::vector<double> multMatrixVector(const std::vector<std::vector<double>>& A, const std::vector<double>& v){
    int n = A.size();
    assert(n > 0);
    int m = v.size();
    assert(m > 0);
    assert(m == n);

    std::vector<double> x(n,0); //solution-vector initialization

    double aij = 0;
    double vj = 0;

    for(int i = 0; i<n; i++){ //fix row
        double sum = 0;

        for(int j = 0; j<n; j++){ //calc row


            aij = A[i][j];
            vj = v[j];
            sum = sum + aij*vj;

        }
        x[i] = sum;
    }

    return x;

}
//Decomment here to use user input
//int main(int argc, char *argv[]){

int main(){

    // Declaring and initializing inputs
    int resolution = 5;
    int iterations = 10;

//    int resolution;
//    int iterations;
//
//    //Decomment here to use user input
//    // using the stringstream class to insert a string and extract an int
//    std::stringstream ss;
//    ss << argv[1];
//    ss >> resolution;
//
//    ss << argv[2];
//    ss >> iterations;


    double res_plus_ghosts = resolution + 2;
    double h = 1. / (resolution - 1);
    double pi = M_PI;
    int number_unknowns = pow(resolution, 2) - (4 * resolution - 4);
    int number_unknowns_in_row = resolution - 4 + 2;
//    int number_interior_unknowns_in_row = resolution - 4 + 1;

    //Starting value for u_h
    // TODO: Change starting value to Zero as in handout stated -> for debugging reasons it is 1, so I can see it.
    double u_h_start = 0;

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
            grid[row][col] = u_h_start;
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

        // Corners
        int row_corners_grid = 2;
        int col_corners_grid = 2;
        double x_coord_corners = h * (col_corners_grid - 1);
        double y_coord_corners = h * (row_corners_grid - 1);
        double u_n_corners = grid[1 + row_corners_grid][col_corners_grid];
        double u_s_corners = grid[row_corners_grid - 1][col_corners_grid];
        double u_e_corners = grid[row_corners_grid][col_corners_grid + 1];
        double u_w_corners = grid[row_corners_grid][col_corners_grid - 1];
        double u_c_corners = grid[row_corners_grid][col_corners_grid];


        // Bottom left
        double b_rhs_corners = 4 * pi * pi * sin(2 * pi * x_coord_corners) * sinh(2 * pi * y_coord_corners)
                + u_s_corners + u_w_corners;
        // Jacobi iteration scheme
        grid_new[row_corners_grid][col_corners_grid] =
                (b_rhs_corners * 1/(h * h) - (( - 1 - 1)  * u_c_corners * 1/(h * h))) * 1 / d_elements;

        // Bottom Right
        row_corners_grid = 2;
        col_corners_grid = resolution-1;
        x_coord_corners = h * (col_corners_grid - 1);
        y_coord_corners = h * (row_corners_grid - 1);
        u_n_corners = grid[1 + row_corners_grid][col_corners_grid];
        u_s_corners = grid[row_corners_grid - 1][col_corners_grid];
        u_e_corners = grid[row_corners_grid][col_corners_grid + 1];
        u_w_corners = grid[row_corners_grid][col_corners_grid - 1];
        u_c_corners = grid[row_corners_grid][col_corners_grid];
        b_rhs_corners = 4 * pi * pi * sin(2 * pi * x_coord_corners) * sinh(2 * pi * y_coord_corners)
                               + u_s_corners + u_e_corners;
        // Jacobi iteration scheme
        grid_new[row_corners_grid][col_corners_grid] =
                (b_rhs_corners * 1/(h * h) - (( - 1 - 1)  * u_c_corners * 1/(h * h))) * 1 / d_elements;

        // Top Left
        row_corners_grid = resolution-1;
        col_corners_grid = 2;
        x_coord_corners = h * (col_corners_grid - 1);
        y_coord_corners = h * (row_corners_grid - 1);
        u_n_corners = grid[1 + row_corners_grid][col_corners_grid];
        u_s_corners = grid[row_corners_grid - 1][col_corners_grid];
        u_e_corners = grid[row_corners_grid][col_corners_grid + 1];
        u_w_corners = grid[row_corners_grid][col_corners_grid - 1];
        u_c_corners = grid[row_corners_grid][col_corners_grid];
        b_rhs_corners = 4 * pi * pi * sin(2 * pi * x_coord_corners) * sinh(2 * pi * y_coord_corners)
                        + u_n_corners + u_w_corners;
        // Jacobi iteration scheme
        grid_new[row_corners_grid][col_corners_grid] =
                (b_rhs_corners * 1/(h * h) - (( - 1 - 1)  * u_c_corners * 1/(h * h))) * 1 / d_elements;

        // Top Right
        row_corners_grid = resolution-1;
        col_corners_grid = resolution-1;
        x_coord_corners = h * (col_corners_grid - 1);
        y_coord_corners = h * (row_corners_grid - 1);
        u_n_corners = grid[1 + row_corners_grid][col_corners_grid];
        u_s_corners = grid[row_corners_grid - 1][col_corners_grid];
        u_e_corners = grid[row_corners_grid][col_corners_grid + 1];
        u_w_corners = grid[row_corners_grid][col_corners_grid - 1];
        u_c_corners = grid[row_corners_grid][col_corners_grid];
        b_rhs_corners = 4 * pi * pi * sin(2 * pi * x_coord_corners) * sinh(2 * pi * y_coord_corners)
                        + u_n_corners + u_e_corners;
        // Jacobi iteration scheme cv
        grid_new[row_corners_grid][col_corners_grid] =
                (b_rhs_corners * 1/(h * h) - (( - 1 - 1)  * u_c_corners * 1/(h * h))) * 1 / d_elements;

        // BC (x, 0) = 0, Bottom
         for (int col = 3; col <= res_plus_ghosts - 4; col++)
        {
            double row = 2;
            double x_coord =  h * (col - 1); // x-value moving along the grid in respect to the loop
            double y_coord =  h * (row - 1); // y-value moving along the grid in respect to the loop


            double u_s = grid[row-1][col];
            double u_c = grid[row][col];
            double b_rhs = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord) + u_s;


            // Jacobi iteration scheme
            grid_new[row][col] = (b_rhs * 1/(h * h) - (( - 1 - 1 - 1)  * u_c * 1/(h * h))) * 1 / d_elements;

        }
        // BC (0, y) = 0, Left
        for (int row = 3; row <= res_plus_ghosts - 4; row++)
        {
            double col = 2;
            double x_coord =  h * (col - 1); // x-value moving along the grid in respect to the loop
            double y_coord =  h * (row - 1); // y-value moving along the grid in respect to the loop


            double u_w = grid[row][col-1];
            double u_c = grid[row][col];
            double b_rhs = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord) + u_w;


            // Jacobi iteration scheme
            grid_new[row][col] = (b_rhs * 1/(h * h) - (( - 1 - 1 - 1)  * u_c * 1/(h * h))) * 1 / d_elements;

        }
        // BC (1, y) = 0, Right
        for (int row = 3; row <= res_plus_ghosts - 4; row++)
        {
            double col = res_plus_ghosts - 3;
            double x_coord =  h * (col - 1); // x-value moving along the grid in respect to the loop
            double y_coord =  h * (row - 1); // y-value moving along the grid in respect to the loop


            double u_e = grid[row][col+1];
            double u_c = grid[row][col];
            double b_rhs = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord) + u_e;


            // Jacobi iteration scheme
            grid_new[row][col] = (b_rhs * 1/(h * h) - (( - 1 - 1 - 1)  * u_c * 1/(h * h))) * 1 / d_elements;

        }

        // BC (x, 1) = 0, Top
        for (int col = 3; col <= res_plus_ghosts - 4; col++)
        {
            double row = res_plus_ghosts - 3;
            double x_coord =  h * (col - 1); // x-value moving along the grid in respect to the loop
            double y_coord =  h * (row - 1); // y-value moving along the grid in respect to the loop

            // TODO: DONE: Set boundary for u_n -> was already set before
            double u_n = grid[1+row][col];
            double u_c = grid[row][col];
            double b_rhs = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord) + u_n;


            // Jacobi iteration scheme
            grid_new[row][col] = (b_rhs * 1/(h * h) - (( - 1 - 1 - 1)  * u_c * 1/(h * h))) * 1 / d_elements;

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
                grid_new[row][col] = (b_rhs * 1/(h * h) - ((-1 - 1 - 1 - 1)  * u_c * 1/(h * h))) * 1 / d_elements;


            }
        }
        // Updating the old grid with the new values
        grid = grid_new;
    }


    //----------------------------------------------//
    // Calculating the residual after the iterations//
    //----------------------------------------------//
    // Putting the solution of the grid into the vector u_h
    std::vector<double> u_h(number_unknowns);
    int counter = 0;
    for (int row = 2; row <= res_plus_ghosts - 3; row++)
    {
        for (int col = 2; col <= res_plus_ghosts - 3; col++)
        {
            u_h[counter] = grid[row][col];
            counter++;
        }
    }

    // Creating the system matrix to calculate the residuals
    std::vector<std::vector<double> > A_h (number_unknowns, std::vector <double> (number_unknowns, 0));
    //----------------------------------//
    // System matrix A_h                //
    //----------------------------------//

    double D_element = 4 + 4 * pi * pi * h * h;

    // Creating the diagonal matrix and the identity matrix
    int size_mat = resolution - 2;
    std::vector<std::vector<double> > D (size_mat, std::vector <double> (size_mat, 0));

    // Filling in the diagonal
    for (int row = 0; row < size_mat; row++) {
        for (int col = 0; col < size_mat; col++) {
            if (row == col) {
                D[row][col] = D_element;
            }
        }
    }

    //Filling in the lower bound
    for (int row = 1; row < size_mat; row++) {
        for (int col = 0; col < size_mat; col++) {
            if (row - 1 == col) {
                D[row][col] = -1;
            }
        }
    }

    //Filling in the upper bound
    for (int row = 0; row < size_mat; row++) {
        for (int col = 1; col < size_mat; col++) {
            if (col - 1 == row) {
                D[row][col] = -1;
            }
        }
    }

//    std::cout << D << std::endl;


    std::vector<std::vector<double> > I (size_mat, std::vector <double> (size_mat, 0));


    // Filling in the diagonal
    for (int row = 0; row < size_mat; row++) {
        for (int col = 0; col < size_mat; col++) {
            if (row == col) {
                I[row][col] = -1;
            }
        }
    }
//    std::cout << I << std::endl;



    // Filling the system Matrix A_h

    // Filling in diagonal matrix
    int col_counter = 0;
    for (int row = 0; row < number_unknowns; row = row + size_mat)
    {
        for (int row_int = 0; row_int < size_mat; row_int++)
        {
            for (int col_int = 0; col_int < size_mat; col_int++)
                A_h[row+row_int][col_counter+col_int] = D[row_int][col_int];
        }
        col_counter = col_counter + size_mat;
    }
//    std::cout << A_h << std::endl;

    // Filling in lower bound identity matrix
    col_counter = 0;
    for (int row = size_mat; row < number_unknowns; row = row + size_mat)
    {
        for (int row_int = 0; row_int < size_mat; row_int++)
        {
            for (int col_int = 0; col_int < size_mat; col_int++)
                A_h[row+row_int][col_counter+col_int] = I[row_int][col_int];
        }
        col_counter = col_counter + size_mat;
    }

    // Filling in upper bound identity matrix
    col_counter = size_mat;
    for (int row = 0; row < number_unknowns; row = row + size_mat)
    {
        for (int row_int = 0; row_int < size_mat; row_int++)
        {
            for (int col_int = 0; col_int < size_mat; col_int++)
                A_h[row+row_int][col_counter+col_int] = I[row_int][col_int];
        }        col_counter = col_counter + size_mat;
        if (col_counter >= number_unknowns){
            break;
        }
    }

    // Creating the b_h vector
    std::vector<double> b_h(number_unknowns);

    // Boundaries for the corner points
    // Bottom
    // Bottom Left Corner
    double x_BL = h;
    double y_BL = h;
    b_h[0] = 4 * pi * pi * sin(2 * pi * x_BL) * sinh(2 * pi * y_BL);

    double x_coord_count = 2;
    // Fill everything in between on the bottom
    for (int row = 1; row < 1 + resolution - 4; row++)
    {
        double x_coord =  h * x_coord_count; // x-value moving along the grid in respect to the loop
        double y_coord =  h; // y-value moving along the grid in respect to the loop
        b_h[row] = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord);
        x_coord_count++;
    }

    // Bottom Right Corner
    double x_BR =  h * (resolution - 2);
    double y_BR =  h;
    int row_BR = 1 + resolution - 4;
    b_h[row_BR] = 4 * pi * pi * sin(2 * pi * x_BR) * sinh(2 * pi * y_BR);

    // Fill everything for the sides and interior
    x_coord_count = 1;
    double y_coord_count = 2;
    // Fill everything in between on the bottom
    for (int row = 1 + resolution - 4 + 1; row < number_unknowns - 1 - (1 + resolution - 4); row++)
    {
        double x_coord =  h * x_coord_count; // x-value moving along the grid in respect to the loop
        double y_coord =  h * y_coord_count; // y-value moving along the grid in respect to the loop
        b_h[row] = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord);

        // Check if end of grid points was reached
        if (number_unknowns_in_row - x_coord_count == 0)
        {
            //Resetting the x coordinate
            x_coord_count = 1;
            // Incrementing the y coordinate
            y_coord_count++;
        }
        else
        {
            x_coord_count++;
        }

    }

    //Top
    // Top Left Corner
    double x_TL =  h;
    double y_TL =  h * (resolution - 2);
    int row_TL = number_unknowns - 1 - (1 + resolution - 4);
    b_h[row_TL] = 4 * pi * pi * sin(2 * pi * x_TL) * sinh(2 * pi * y_TL) +
                  sin(2 * pi * x_TL) * sinh(2 * pi);

    x_coord_count = 2;
    // Fill everything in between on the top
    for (int row = number_unknowns - 1 - (1 + resolution - 4) + 1; row < number_unknowns - 1; row++)
    {

        double x_coord =  h * x_coord_count; // x-value moving along the grid in respect to the loop
        double y_coord =  h * (resolution - 2); // y-value moving along the grid in respect to the loop
        b_h[row] = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord) +
                   sin(2 * pi * x_coord) * sinh(2 * pi);
        x_coord_count++;

    }

    // Top Right Corner
    double x_TR =  h * (resolution - 2);
    double y_TR =  h * (resolution - 2);
    int row_TR = number_unknowns - 1;
    b_h[row_TR] = 4 * pi * pi * sin(2 * pi * x_TR) * sinh(2 * pi * y_TR) +
                  sin(2 * pi * x_TR) * sinh(2 * pi);

    // Multiplying with the constant
    for (int row = 0; row < number_unknowns; row++)
    {
        b_h[row] = 1/pow(h,2) * b_h[row];
    }


    // Finally calculate the residual r_h
    std::vector<double> r_h(number_unknowns);
    std::vector<double> A_h_mult_u_h(number_unknowns); // temp to save result off mvmult
    A_h_mult_u_h = multMatrixVector(A_h, u_h);
    for (int row = 0; row < number_unknowns; row++)
    {
        r_h[row] = A_h_mult_u_h[row] - b_h[row];
    }
//    std::cout << r_h << std::endl;


    //----------------------------------------------//
    // Calculating error for analytical solution    //
    //----------------------------------------------//

    // Analytical solution
    //---------------------
    std::vector<double> u_p(number_unknowns);
    x_coord_count = 1;
    y_coord_count = 1;

    for (int row = 0; row < number_unknowns; row++)
    {
        double x_coord =  h * x_coord_count; // x-value moving along the grid in respect to the loop
        double y_coord =  h * y_coord_count; // y-value moving along the grid in respect to the loop
        u_p[row] = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord);

        // Check if end of grid points was reached
        if (number_unknowns_in_row - x_coord_count == 0)
        {
            //Resetting the x coordinate
            x_coord_count = 1;
            // Incrementing the y coordinate
            y_coord_count++;
        }
        else
        {
            x_coord_count++;
        }
    }
//    std::cout << u_p[0] << std::endl;





    // Calculating the error
    //-----------------------
    std::vector<double> error_h(number_unknowns);
    for (int row = 0; row < number_unknowns; row++)
    {
        error_h[row] = u_p[row] - u_p[row];
    }


    //----------------------------------//
    // Norms                            //
    //----------------------------------//

    //Euclidean Norm
    // ----------------
    // TODO: Check how to write this into function
    // For Residual r_h
    double euclidean_norm_residual = 0;
    for (int row = 0; row < number_unknowns; row++)
    {
        euclidean_norm_residual = euclidean_norm_residual + pow(abs(r_h[row]), 2);
    }
    euclidean_norm_residual = sqrt(euclidean_norm_residual);

    // For Error error_h
    double euclidean_norm_error = 0;
    for (int row = 0; row < number_unknowns; row++)
    {
        euclidean_norm_error = euclidean_norm_error + pow(abs(r_h[row]), 2);
    }
    euclidean_norm_error = sqrt(euclidean_norm_error);

    std::cout << euclidean_norm_residual << std::endl;
    std::cout << euclidean_norm_error << std::endl;


    //Maximum Norm
    // ----------------
    // For Residual r_h
    double maximum_norm_residual = 0;
    for (int row = 0; row < number_unknowns; row++)
    {
        if(abs(r_h[row]) >= maximum_norm_residual){
            maximum_norm_residual = abs(r_h[row]);
        }
    }
    // For Error error_h
    double maximum_norm_error = 0;
    for (int row = 0; row < number_unknowns; row++)
    {
        if(abs(r_h[row]) >= maximum_norm_error){
            maximum_norm_error = abs(r_h[row]);
        }
    }
    std::cout << maximum_norm_residual << std::endl;
    std::cout << maximum_norm_error << std::endl;



    std::cout << h << std::endl;








    return 0;
};