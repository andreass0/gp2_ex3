//
// Group Project 2/ Ex3.3
//

# include  <iostream>
# include <fstream>
# include <chrono>
# include <string.h>
# include <vector>
# include "cmath"
# include <eigen-3.4.0/Eigen/Dense>




using Eigen::MatrixXd;
using Eigen::VectorXd;


int main(){

    // Declaring and initializing inputs
    int resolution = 5;
    int iterations = 100;
    double res_plus_ghosts = resolution + 2;
//    int resol_int = resolution;
    double h = 1. / (resolution - 1);
    double pi = M_PI;
    int number_unknowns = pow(resolution, 2) - (4 * resolution - 4);

    // Creating the matrix and the needed arrays
    MatrixXd A_h(number_unknowns, number_unknowns);
    VectorXd u_h(number_unknowns);
    VectorXd b_h(number_unknowns);

    // Filling in the matrix and the arrays


    //----------------------------------//
    //Unknown-vector: u_h
    //----------------------------------//

    //u_h (filling in with the starting value for the iteration)
    for (int row = 0; row < number_unknowns; row++)
    {
        u_h(row) = 0;
    }

    //----------------------------------//
    //RHS-vector: b_h
    //----------------------------------//

    // Boundaries for the corner points
    // Bottom
    // Bottom Left Corner
    double x_BL =  h;
    double y_BL =  h;
    b_h(0) = 4 * pi * pi * sin(2 * pi * x_BL) * sinh(2 * pi * y_BL);

    double x_coord_count = 2;
    // Fill everything in between on the bottom
    for (int row = 1; row < 1 + resolution - 4; row++)
    {
        double x_coord =  h * x_coord_count; // x-value moving along the grid in respect to the loop
        double y_coord =  h; // y-value moving along the grid in respect to the loop
        b_h(row) = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord);
        x_coord_count++;
    }

    // Bottom Right Corner
    double x_BR =  h * (resolution - 2);
    double y_BR =  h;
    int row_BR = 1 + resolution - 4;
    b_h(row_BR) = 4 * pi * pi * sin(2 * pi * x_BR) * sinh(2 * pi * y_BR);

    // Fill everything for the sides and interior
    int number_unknowns_in_row = resolution - 4 + 2;
    int number_interior_unknowns_in_row = resolution - 4 + 1;
    x_coord_count = 1;
    double y_coord_count = 2;
    // Fill everything in between on the bottom
    for (int row = 1 + resolution - 4 + 1; row < number_unknowns - 1 - (1 + resolution - 4); row++)
    {
        double x_coord =  h * x_coord_count; // x-value moving along the grid in respect to the loop
        double y_coord =  h * y_coord_count; // y-value moving along the grid in respect to the loop
        b_h(row) = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord);

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
    b_h(row_TL) = 4 * pi * pi * sin(2 * pi * x_TL) * sinh(2 * pi * y_TL) +
            sin(2 * pi * x_TL) * sinh(2 * pi);

    x_coord_count = 2;
    // Fill everything in between on the top
    for (int row = number_unknowns - 1 - (1 + resolution - 4) + 1; row < number_unknowns - 1; row++)
    {

        double x_coord =  h * x_coord_count; // x-value moving along the grid in respect to the loop
        double y_coord =  h * (resolution - 2); // y-value moving along the grid in respect to the loop
        b_h(row) = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord) +
                   sin(2 * pi * x_coord) * sinh(2 * pi);
        x_coord_count++;

    }

    // Top Right Corner
    double x_TR =  h * (resolution - 2);
    double y_TR =  h * (resolution - 2);
    int row_TR = number_unknowns - 1;
    b_h(row_TR) = 4 * pi * pi * sin(2 * pi * x_TR) * sinh(2 * pi * y_TR) +
            sin(2 * pi * x_TR) * sinh(2 * pi);

    // Multiplying with the constant
    b_h = 1/pow(h,2) * b_h;



    //----------------------------------//
    // System matrix A_h                //
    //----------------------------------//

    double D_element = 4 + 4 * pi * pi * h * h;

    // Creating the diagonal matrix and the identity matrix
    int size_mat = resolution - 2;
    MatrixXd D(size_mat, size_mat);
//    D <<    D_element, -1, 0,
//            -1, D_element, -1,
//            0, -1, D_element;

    // Filling in the diagonal
    for (int row = 0; row < size_mat; row++) {
        for (int col = 0; col < size_mat; col++) {
            if (row == col) {
                D(row, col) = D_element;
            }
        }
    }

    //Filling in the lower bound
    for (int row = 1; row < size_mat; row++) {
        for (int col = 0; col < size_mat; col++) {
            if (row - 1 == col) {
                D(row, col) = -1;
            }
        }
    }

    //Filling in the upper bound
    for (int row = 0; row < size_mat; row++) {
        for (int col = 1; col < size_mat; col++) {
            if (col - 1 == row) {
                D(row, col) = -1;
            }
        }
    }

//    std::cout << D << std::endl;


        MatrixXd I(size_mat, size_mat);
//    I <<    -1, 0, 0,
//            0, -1, 0,
//            0, 0, -1;

    I.setIdentity();
    I = I * -1;
//    std::cout << I << std::endl;



    // Filling the system Matrix A_h

    //Filling everything with zeros
    A_h.setZero();
//    for (int row = 0; row < number_unknowns; row++)
//    {
//        for (int col = 0; col < number_unknowns; col++)
//        {
//            A_h(row, col) = 0;
//        }
//    }

    // Filling in diagonal matrix
    int col_counter = 0;
    for (int row = 0; row < number_unknowns; row = row + size_mat)
    {
        A_h.block(row,col_counter,size_mat, size_mat) = D;
        col_counter = col_counter + size_mat;
    }
//    std::cout << A_h << std::endl;

    // Filling in lower bound identity matrix
    col_counter = 0;
    for (int row = size_mat; row < number_unknowns; row = row + size_mat)
    {
        A_h.block(row,col_counter,size_mat, size_mat) = I;
        col_counter = col_counter + size_mat;
    }

    // Filling in upper bound identity matrix
    col_counter = size_mat;
    for (int row = 0; row < number_unknowns; row = row + size_mat)
    {
        A_h.block(row,col_counter,size_mat, size_mat) = I;
        col_counter = col_counter + size_mat;
        if (col_counter >= number_unknowns){
            break;
        }
    }


    //----------------------------------//
    // Jacobi Method                    //
    //----------------------------------//

    //Declaring and initializing the residual matrix
    MatrixXd R_h(number_unknowns, number_unknowns);
    R_h.setZero();
    // Filling in lower bound identity matrix
    col_counter = 0;
    for (int row = size_mat; row < number_unknowns; row = row + size_mat)
    {
        R_h.block(row,col_counter,size_mat, size_mat) = I;
        col_counter = col_counter + size_mat;
    }

    // Filling in upper bound identity matrix
    col_counter = size_mat;
    for (int row = 0; row < number_unknowns; row = row + size_mat)
    {
        R_h.block(row,col_counter,size_mat, size_mat) = I;
        col_counter = col_counter + size_mat;
        if (col_counter >= number_unknowns)
        {
            break;
        }
    }

//    std::cout << R_h << std::endl;

    // Declaring and calculating the diagonal matrix
    MatrixXd D_h(number_unknowns, number_unknowns);
    D_h = A_h - R_h;
    MatrixXd D_h_inv(number_unknowns, number_unknowns);
    D_h_inv = D_h.inverse();
//    std::cout << D_h_inv << std::endl;

    // Creating iteration vector to update
    VectorXd u_h_new(number_unknowns);
    u_h_new = u_h;

    for (int iter = 1; iter <= iterations; iter++)
    {
        u_h_new = D_h_inv * (b_h - R_h * u_h);
        u_h = u_h_new;
//        std::cout << u_h << std::endl;
    }

    //----------------------------------------------//
    // Calculating the residual after the iterations//
    //----------------------------------------------//

    VectorXd r_h(number_unknowns);
    r_h = A_h * u_h - b_h;
//    std::cout << r_h << std::endl;

    //----------------------------------------------//
    // Calculating error for analytical solution    //
    //----------------------------------------------//

    // Analytical solution
    //---------------------
    VectorXd u_p(number_unknowns);
    x_coord_count = 1;
    y_coord_count = 1;

    for (int row = 0; row < u_p.rows(); row++)
    {
        double x_coord =  h * x_coord_count; // x-value moving along the grid in respect to the loop
        double y_coord =  h * y_coord_count; // y-value moving along the grid in respect to the loop
        u_p(row) = 4 * pi * pi * sin(2 * pi * x_coord) * sinh(2 * pi * y_coord);

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
//    std::cout << u_p << std::endl;


    // Calculating the error
    //-----------------------
    VectorXd error_h(number_unknowns);
    error_h = u_h - u_p;

    //----------------------------------//
    // Norms                            //
    //----------------------------------//

    //Euclidean Norm
    // ----------------
    // TODO: Check how to write this into function
    // For Residual r_h
    double euclidean_norm_residual = 0;
    for (int row = 0; row < r_h.rows(); row++)
    {
        euclidean_norm_residual = euclidean_norm_residual + pow(abs(r_h(row)), 2);
    }
    euclidean_norm_residual = sqrt(euclidean_norm_residual);

    // For Error error_h
    double euclidean_norm_error = 0;
    for (int row = 0; row < r_h.rows(); row++)
    {
        euclidean_norm_error = euclidean_norm_error + pow(abs(r_h(row)), 2);
    }
    euclidean_norm_error = sqrt(euclidean_norm_error);

    std::cout << euclidean_norm_residual << std::endl;
    std::cout << euclidean_norm_error << std::endl;


    //Maximum Norm
    // ----------------
    // For Residual r_h
    double maximum_norm_residual = 0;
    for (int row = 0; row < r_h.rows(); row++)
    {
        if(abs(r_h(row)) >= maximum_norm_residual){
            maximum_norm_residual = abs(r_h(row));
        }
    }
    // For Error error_h
    double maximum_norm_error = 0;
    for (int row = 0; row < r_h.rows(); row++)
    {
        if(abs(r_h(row)) >= maximum_norm_error){
            maximum_norm_error = abs(r_h(row));
        }
    }
    std::cout << maximum_norm_residual << std::endl;
    std::cout << maximum_norm_error << std::endl;








//    std::cout << b_h << std::endl;
//    std::cout << D << std::endl;
//    std::cout << I << std::endl;
//    std::cout << A_h << std::endl;
    std::cout << u_h << std::endl;





//    std::cout << sinh(2 * pi) << std::endl;
//    std::cout << sin(2 * pi) << std::endl;






    return 0;
};
