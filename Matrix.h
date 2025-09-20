// =======================================================================================================================================
// Matrix.h
// Header file for matrix operations
// Author: Runkai Xie
// Version:

// Member functions:
// 0. Matrix<data type> Target_Matrix: Define a Matrix class template.
// 1. Target_Matrix.InitialMatrix(): Initialize the matrix by user input.
// 2. Target_Matrix.DisplayMatrix(): Display the matrix.

// Global functions:
// 3.  Rank(const Matrix<ElementsType> &matrix):   Calculate the rank of the matrix.
// 4.  Transpose(const Matrix<ElementsType> &matrix):   Transpose the matrix.
// 5.  Determinant(const Matrix<ElementsType> &matrix):   Calculate the determinant of the matrix.
// 6.  Sum(const Matrix<ElementsType> &matrix1, const Matrix<ElementsType> &matrix2):   Calculate the sum of two matrices.
// 7.  Subtract(const Matrix<ElementsType> &substracter, const Matrix<ElementsType> &minuend):   Calculate the subtract of two matrices.
// 8.  Multiply(const Matrix<ElementsType> &matrix1 / Number, const Matrix<ElementsType> &matrix2):
//     Calculate the product of two matrices / number and matrix.
// 9.  Adjugate(const Matrix<ElementsType> &matrix):   Calculate the adjugate matrix of the target matrix.
// 10. Inverse(const Matrix<ElementsType> &matrix):   Calculate the inverse of the matrix.
// 11. Jacobian(std::function<std::vector<ElementsType>(const std::vector<ElementsType>&)> f, const std::vector<ElementsType>& x, ElementsType h = 1e-6):  
//     Calculate the Jacobian matrix of the target matrix.
// =======================================================================================================================================

#include <iostream>
#include <vector>

template <typename ElementsType>
class Matrix
{
public:
    // Data members:
    int rows;
    int columns;
    std::vector<std::vector<ElementsType>> elements;
    // Member function declarations:
    void InitialMatrix();
    void DisplayMatrix();
};

// =======================================================================================================================================
//                                              Function for matrix initialization.
template <typename ElementsType>
void Matrix<ElementsType>::InitialMatrix()
{
    std::cout << std::endl
              << " Please enter rows and columns of matrix: " << std::endl;
    int n, m;
    std::cout << " rows    : ";
    std::cin >> n;
    std::cout << " columns : ";
    std::cin >> m;
    this->rows = n;
    this->columns = m;
    this->elements.resize(n, std::vector<ElementsType>(m));
    std::cout << std::endl
              << " Please enter the matrix: " << std::endl;
    for (int i = 0; i < this->rows; i++)
    {
        std::cout << " ";
        for (int j = 0; j < this->columns; j++)
        {
            std::cin >> this->elements[i][j];
        }
    }
}
// =======================================================================================================================================

// =======================================================================================================================================
//                                               Function for displaying matrix.
template <typename ElementsType>
void Matrix<ElementsType>::DisplayMatrix()
{
    ElementsType zero = 0;
    for (int i = 0; i < this->rows; i++)
    {
        std::cout << " ";
        for (int j = 0; j < this->columns; j++)
        {
            if (this->elements[i][j] < 1e-8 && this->elements[i][j] > -1e-8)
            {
                std::cout << zero << " "; // Display zero for eliminating operator errors.
            }
            else
            {
                std::cout << this->elements[i][j] << " ";
            }
        }
        std::cout << std::endl;
    }
}
// =======================================================================================================================================

// =======================================================================================================================================
//                                       Function for calculating the rank of the matrix.
template <typename ElementsType>
int Rank(const Matrix<ElementsType> &matrix)
{
    int rank = 0;
    int n = matrix.rows, m = matrix.columns;
    int row = 0, col = 0;
    while (row < n && col < m)
    {
        // Find the longest non-zero row.
        int row_notzero = row;
        while (row_notzero < n && matrix.elements[row_notzero][col] == 0)
        {
            row_notzero++;
        }
        if (row_notzero == n)
        {
            col++;
            continue;
        }
        if (row_notzero != row) // Swap the current row with the longest non-zero row, transform the matrix into row echelon form.
        {
            std::swap(matrix.elements[row_notzero], matrix.elements[row]);
        }
        for (int i = row + 1; i < n; i++) // Eliminate the elements below the non-zero row.
        {
            ElementsType ratio = matrix.elements[i][col] / matrix.elements[row][col];
            for (int j = col; j < m; j++)
            {
                matrix.elements[i][j] -= ratio * matrix.elements[row][j];
            }
        }
        row++;
        col++;
        rank++;
    }
    return rank;
}
// =======================================================================================================================================

// =======================================================================================================================================
//                                              Function for transposing matrix.
template <typename ElementsType>
Matrix<ElementsType> Transpose(const Matrix<ElementsType> &matrix)
{
    int i, j;
    Matrix<ElementsType> transposed_matrix;
    transposed_matrix.rows = matrix.columns;
    transposed_matrix.columns = matrix.rows;
    transposed_matrix.elements.resize(transposed_matrix.rows, std::vector<ElementsType>(transposed_matrix.columns));
    for (j = 0; j < matrix.columns; j++)
    {
        for (i = 0; i < matrix.rows; i++)
        {
            transposed_matrix.elements[j][i] = matrix.elements[i][j];
        }
    }
    std::cout << std::endl
              << " The transposed matrix is : " << std::endl;
    transposed_matrix.DisplayMatrix();
    return transposed_matrix;
}
// =======================================================================================================================================

// =======================================================================================================================================
//                                   Function for calculating the determinant of the matrix.
template <typename ElementsType>
ElementsType Determinant(const Matrix<ElementsType> &matrix)
{
    ElementsType determinant = 0;
    if (matrix.rows != matrix.columns)
    {
        std::cout << " Error! The matrix is not a square matrix! " << std::endl;
        return determinant;
    }
    if (matrix.rows == 1)
    {
        return matrix.elements[0][0];
    }
    if (matrix.rows == 2)
    {
        return matrix.elements[0][0] * matrix.elements[1][1] - matrix.elements[0][1] * matrix.elements[1][0];
    }
    for (int i = 0; i < matrix.rows; i++)
    {
        Matrix<ElementsType> temp_matrix;
        temp_matrix.rows = matrix.rows - 1;
        temp_matrix.columns = matrix.columns - 1;
        temp_matrix.elements.resize(temp_matrix.rows, std::vector<ElementsType>(temp_matrix.columns));
        for (int j = 1; j < matrix.rows; j++)
        {
            int l = 0;
            for (int k = 0; k < matrix.columns; k++)
            {
                if (k == i)
                {
                    continue;
                }
                temp_matrix.elements[j - 1][l] = matrix.elements[j][k];
                l++;
            }
        }
        determinant += (i % 2 == 0 ? 1 : -1) * matrix.elements[0][i] * Determinant(temp_matrix);
    }
    return determinant;
}
// =======================================================================================================================================

// =======================================================================================================================================
//                                               Function for sum of two matrix.
template <typename ElementsType>
Matrix<ElementsType> Sum(const Matrix<ElementsType> &matrix1, const Matrix<ElementsType> &matrix2)
{
    Matrix<ElementsType> sum_matrix;
    if (matrix1.rows != matrix2.rows || matrix1.columns != matrix2.columns)
    {
        std::cout << " Error! The two matrix can not be added! " << std::endl;
        return sum_matrix;
    }
    sum_matrix.rows = matrix1.rows;
    sum_matrix.columns = matrix1.columns;
    sum_matrix.elements.resize(sum_matrix.rows, std::vector<ElementsType>(sum_matrix.columns));
    for (int i = 0; i < sum_matrix.rows; i++)
    {
        for (int j = 0; j < sum_matrix.columns; j++)
        {
            sum_matrix.elements[i][j] = matrix1.elements[i][j] + matrix2.elements[i][j];
        }
    }
    std::cout << std::endl
              << " The sum of two matrix is : " << std::endl;
    sum_matrix.DisplayMatrix();
    return sum_matrix;
}
// =======================================================================================================================================

// =======================================================================================================================================
//                                    Function for calculating the subtract of two matrix.
template <typename ElementsType>
Matrix<ElementsType> Subtract(const Matrix<ElementsType> &substracter, const Matrix<ElementsType> &minuend)
{
    Matrix<ElementsType> difference_matrix;
    if (substracter.rows != minuend.rows || substracter.columns != minuend.columns)
    {
        std::cout << " Error! The two matrix can not be subtracted! " << std::endl;
        return difference_matrix;
    }
    difference_matrix.rows = substracter.rows;
    difference_matrix.columns = substracter.columns;
    difference_matrix.elements.resize(difference_matrix.rows, std::vector<ElementsType>(difference_matrix.columns));
    for (int i = 0; i < difference_matrix.rows; i++)
    {
        for (int j = 0; j < difference_matrix.columns; j++)
        {
            difference_matrix.elements[i][j] = substracter.elements[i][j] - minuend.elements[i][j];
        }
    }
    std::cout << std::endl
              << " The subtract of two matrix is : " << std::endl;
    difference_matrix.DisplayMatrix();
    return difference_matrix;
}
// =======================================================================================================================================

// =======================================================================================================================================
//                                Overloaded function for calculating the product of two matrix.
template <typename ElementsType>
Matrix<ElementsType> Multiply(const Matrix<ElementsType> &matrix1, const Matrix<ElementsType> &matrix2)
{
    Matrix<ElementsType> product_matrix;
    if (matrix1.columns != matrix2.rows)
    {
        std::cout << " Error! The two matrix can not be multiplied! " << std::endl;
        return product_matrix;
    }
    product_matrix.rows = matrix1.rows;
    product_matrix.columns = matrix2.columns;
    product_matrix.elements.resize(product_matrix.rows, std::vector<ElementsType>(product_matrix.columns, 0));
    for (int i = 0; i < product_matrix.rows; i++)
    {
        for (int j = 0; j < product_matrix.columns; j++)
        {
            for (int k = 0; k < matrix1.columns; k++)
            {
                product_matrix.elements[i][j] += matrix1.elements[i][k] * matrix2.elements[k][j];
            }
        }
    }
    std::cout << std::endl
              << " The product is : " << std::endl;
    product_matrix.DisplayMatrix();
    return product_matrix;
}
//                           Overloaded function for calculating the product of a scalar and a matrix.
template <typename ElementsType, typename ScalarType>
Matrix<ElementsType> Multiply(ScalarType scalar, const Matrix<ElementsType> &matrix)
{
    ElementsType Scalar = static_cast<ElementsType>(scalar);
    Matrix<ElementsType> product_matrix;
    product_matrix.rows = matrix.rows;
    product_matrix.columns = matrix.columns;
    product_matrix.elements.resize(product_matrix.rows, std::vector<ElementsType>(product_matrix.columns, 0));
    for (int i = 0; i < product_matrix.rows; i++)
    {
        for (int j = 0; j < product_matrix.columns; j++)
        {
            product_matrix.elements[i][j] = Scalar * matrix.elements[i][j];
        }
    }
    std::cout << std::endl
              << " The product is : " << std::endl;
    product_matrix.DisplayMatrix();
    return product_matrix;
}
// =======================================================================================================================================

// =======================================================================================================================================
//                                           Function for calculating the adjugate matrix.
template <typename ElementsType>
Matrix<ElementsType> Adjugate(const Matrix<ElementsType> &matrix)
{
    Matrix<ElementsType> adjugate_matrix;
    if (matrix.rows != matrix.columns)
    {
        std::cout << " Error! The matrix is not a square matrix! " << std::endl;
        return adjugate_matrix;
    }
    int n = matrix.rows;
    adjugate_matrix.rows = n;
    adjugate_matrix.columns = n;
    adjugate_matrix.elements.resize(n, std::vector<ElementsType>(n, 0));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            Matrix<ElementsType> sub_matrix;
            sub_matrix.rows = n - 1;
            sub_matrix.columns = n - 1;
            sub_matrix.elements.resize(n - 1, std::vector<ElementsType>(n - 1));
            int sub_rows = 0;
            for (int r = 0; r < n; r++)
            {
                if (r == i) continue;
                int sub_columns = 0;
                for (int c = 0; c < n; c++)
                {
                    if (c == j) continue;
                    sub_matrix.elements[sub_rows][sub_columns] = matrix.elements[r][c];
                    sub_columns++;
                }
                sub_rows++;
            }
            adjugate_matrix.elements[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * Determinant(sub_matrix);
        }
    }
    return adjugate_matrix;
}
// =======================================================================================================================================

// =======================================================================================================================================
//                                        Function for calculating the inverse of the matrix.
template <typename ElementsType>
Matrix<ElementsType> Inverse(const Matrix<ElementsType> &matrix)
{
    Matrix<ElementsType> inverse_matrix;
    if (matrix.rows != matrix.columns)
    {
        std::cout << " Error! The matrix is not a square matrix! " << std::endl;
        return inverse_matrix;
    }
    ElementsType determinant = Determinant(matrix);
    if (determinant == 0)
    {
        std::cout << " Error! The matrix is singular and cannot be inverted! " << std::endl;
        return inverse_matrix;
    }

    // Using the adjugate matrix to calculate the inverse matrix.
    Matrix<ElementsType> adjugate_matrix = Adjugate(matrix);
    inverse_matrix.rows = matrix.rows;
    inverse_matrix.columns = matrix.columns;
    inverse_matrix.elements.resize(inverse_matrix.rows, std::vector<ElementsType>(inverse_matrix.columns, 0));
    for (int i = 0; i < inverse_matrix.rows; i++)
    {
        for (int j = 0; j < inverse_matrix.columns; j++)
        {
            inverse_matrix.elements[i][j] = adjugate_matrix.elements[i][j] / determinant;
        }
    }

    // Using Gauss-Jordan elimination to calculate the inverse matrix.
    // int n = matrix.rows;
    // Matrix<ElementsType> augmented_matrix;
    // augmented_matrix.rows = n;
    // augmented_matrix.columns = 2 * n;
    // augmented_matrix.elements.resize(augmented_matrix.rows, std::vector<ElementsType>(augmented_matrix.columns, 0));
    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++)
    //     {
    //         augmented_matrix.elements[i][j] = matrix.elements[i][j];
    //     }
    //     augmented_matrix.elements[i][i + n] = 1;
    // }
    // for (int i = 0; i < n; i++) // Normalize the pivot row.
    // {
    //     ElementsType pivot = augmented_matrix.elements[i][i];
    //     for (int j = 0; j < 2 * n; j++)
    //     {
    //         augmented_matrix.elements[i][j] /= pivot;
    //     }
    //     for (int k = 0; k < n; k++)
    //     {
    //         if (k != i)
    //         {
    //             ElementsType ratio = augmented_matrix.elements[k][i];
    //             for (int j = 0; j < 2 * n; j++) // Eliminate other elements except the pivot element in the current colum.
    //             {
    //                 augmented_matrix.elements[k][j] -= ratio * augmented_matrix.elements[i][j]; 
    //             }
    //         }
    //     }
    // }
    // inverse_matrix.rows = n;
    // inverse_matrix.columns = n;
    // inverse_matrix.elements.resize(inverse_matrix.rows, std::vector<ElementsType>(inverse_matrix.columns, 0));
    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++)
    //     {
    //         inverse_matrix.elements[i][j] = augmented_matrix.elements[i][j + n];
    //     }
    // }
    std::cout << std::endl << " The inverse matrix is : " << std::endl;
    inverse_matrix.DisplayMatrix();
    return inverse_matrix;
}
// =======================================================================================================================================

// =======================================================================================================================================
//                              Function for calculating the Jacobian matrix of the target matrix.
template <typename ElementsType>
Matrix<ElementsType> Jacobian(std::function<std::vector<ElementsType>(const std::vector<ElementsType>&)> f, const std::vector<ElementsType>& x, ElementsType h = 1e-6)
{
    int n = x.size();
    std::vector<ElementsType> fx = f(x);
    int m = fx.size();
    Matrix<ElementsType> jacobian_matrix;
    jacobian_matrix.rows = m;
    jacobian_matrix.columns = n;
    jacobian_matrix.elements.resize(m, std::vector<ElementsType>(n, 0));

    for (int j = 0; j < n; ++j)
    {
        std::vector<ElementsType> x_forward = x;
        std::vector<ElementsType> x_backward = x;
        x_forward[j] += h;
        x_backward[j] -= h;
        std::vector<ElementsType> f_forward = f(x_forward);
        std::vector<ElementsType> f_backward = f(x_backward);
        for (int i = 0; i < m; ++i)
        {
            jacobian_matrix.elements[i][j] = (f_forward[i] - f_backward[i]) / (2 * h);
        }
    }
    return jacobian_matrix;
}
// =======================================================================================================================================