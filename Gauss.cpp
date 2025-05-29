#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>

using namespace std;
using namespace chrono;

template <typename NumType>
class Vector {
    NumType* elements;
    size_t length;

public:
    Vector() : elements(nullptr), length(0) {}
    explicit Vector(size_t n) : elements(new NumType[n]()), length(n) {}
    ~Vector() { delete[] elements; }

    Vector(const Vector& other) : elements(new NumType[other.length]), length(other.length) {
        for (size_t i = 0; i < length; ++i)
            elements[i] = other.elements[i];
    }

    size_t size() const { return length; }

    void exchange(Vector& other) {
        swap(elements, other.elements);
        swap(length, other.length);
    }

    Vector& operator=(const Vector& other) {
        if (this != &other) {
            Vector temp(other);
            exchange(temp);
        }
        return *this;
    }

    NumType& check_idx(size_t idx) {
        if (idx >= length) throw out_of_range("Index out of bounds");
        return elements[idx];
    }

    const NumType& check_idx(size_t idx) const {
        if (idx >= length) throw out_of_range("Index out of bounds");
        return elements[idx];
    }

    Vector& operator/=(NumType divisor) {
        if (divisor == 0) throw invalid_argument("Division by zero");
        for (size_t i = 0; i < length; ++i)
            elements[i] /= divisor;
        return *this;
    }

    Vector& operator-=(const Vector& right) {
        if (length != right.length) throw invalid_argument("Vector sizes differ");
        for (size_t i = 0; i < length; ++i)
            elements[i] -= right.elements[i];
        return *this;
    }

    Vector& operator+=(const Vector& right) {
        if (length != right.length) throw invalid_argument("Vector sizes differ");
        for (size_t i = 0; i < length; ++i)
            elements[i] += right.elements[i];
        return *this;
    }

    Vector operator*(NumType scalar) const {
        Vector result(length);
        for (size_t i = 0; i < length; ++i)
            result.elements[i] = elements[i] * scalar;
        return result;
    }
};

template <typename NumType>
class Matrix : public Vector<Vector<NumType>> {
    size_t row_count;
    size_t col_count;

public:
    Matrix() : row_count(0), col_count(0) {}
    Matrix(size_t r, size_t c) : Vector<Vector<NumType>>(r), row_count(r), col_count(c) {
        for (size_t i = 0; i < r; ++i)
            this->check_idx(i) = Vector<NumType>(c);
    }

    size_t rows() const { return row_count; }
    size_t cols() const { return col_count; }

    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            Matrix temp(other);
            exchange(temp);
        }
        return *this;
    }

    Vector<NumType>& row(size_t r) {
        if (r >= row_count) throw out_of_range("Row index invalid");
        return Vector<Vector<NumType>>::check_idx(r);
    }

    const Vector<NumType>& row(size_t r) const {
        if (r >= row_count) throw out_of_range("Row index invalid");
        return Vector<Vector<NumType>>::check_idx(r);
    }

    void exchange(Matrix& other) {
        Vector<Vector<NumType>>::exchange(other);
        swap(row_count, other.row_count);
        swap(col_count, other.col_count);
    }

    void swap_rows(size_t r1, size_t r2) {
        if (r1 >= row_count || r2 >= row_count) throw out_of_range("Row indices invalid");
        if (r1 != r2) {
            row(r1).exchange(row(r2));
        }
    }

    Matrix operator*(NumType scalar) const {
        Matrix result(*this);
        for (size_t i = 0; i < result.rows(); ++i)
            result.row(i) *= scalar;
        return result;
    }
};

template <typename NumType>
struct Solution {
    bool is_consistent = false;
    bool is_unique = false;
    Vector<NumType> particular_solution;
    Vector<NumType>* basis_vectors = nullptr;
    size_t basis_dimension = 0;
    size_t* free_variables = nullptr;
    size_t matrix_rank = 0;

    Solution() = default;
    ~Solution() {
        delete[] basis_vectors;
        delete[] free_variables;
    }
};

template <typename NumType>
Solution<NumType> SolveLinearSystem(Matrix<NumType> coeff_matrix, Vector<NumType> const_terms) {
    const size_t equations = coeff_matrix.rows();
    const size_t variables = coeff_matrix.cols();
    const NumType tolerance = numeric_limits<NumType>::epsilon() * 1000;
    Solution<NumType> solution;
    solution.is_consistent = true;
    size_t rank = 0;
    size_t* pivot_columns = new size_t[equations];

    for (size_t i = 0; i < equations; ++i)
        pivot_columns[i] = variables;

    for (size_t col = 0; col < variables && rank < equations; ++col) {
        size_t max_row = rank;
        for (size_t i = rank; i < equations; ++i) {
            if (std::abs(coeff_matrix.row(i).check_idx(col)) > std::abs(coeff_matrix.row(max_row).check_idx(col)))
                max_row = i;
        }

        if (std::abs(coeff_matrix.row(max_row).check_idx(col)) < tolerance) continue;

            pivot_columns[rank] = col;
            coeff_matrix.swap_rows(rank, max_row);
            swap(const_terms.check_idx(rank), const_terms.check_idx(max_row));

            NumType pivot = coeff_matrix.row(rank).check_idx(col);
            coeff_matrix.row(rank) /= pivot;
        const_terms.check_idx(rank) /= pivot;

        for (size_t i = 0; i < equations; ++i) {
            if (i != rank && std::abs(coeff_matrix.row(i).check_idx(col)) > tolerance) {
                NumType factor = coeff_matrix.row(i).check_idx(col);
                coeff_matrix.row(i) -= coeff_matrix.row(rank) * factor;
                const_terms.check_idx(i) -= const_terms.check_idx(rank) * factor;
            }
        }
        ++rank;
    }

    for (size_t i = rank; i < equations; ++i) {
        if (std::abs(const_terms.check_idx(i)) > tolerance) {
            solution.is_consistent = false;
            delete[] pivot_columns;
            return solution;
        }
    }

    solution.matrix_rank = rank;
    solution.is_unique = (rank == variables);
    solution.particular_solution = Vector<NumType>(variables);

    bool* is_pivot_col = new bool[variables]();
    for (size_t i = 0; i < rank; ++i) {
        if (pivot_columns[i] < variables)
            is_pivot_col[pivot_columns[i]] = true;
    }

    for (ptrdiff_t i = rank - 1; i >= 0; --i) {
        size_t col = pivot_columns[i];
        if (col >= variables) continue;

        solution.particular_solution.check_idx(col) = const_terms.check_idx(i);
        for (size_t j = col + 1; j < variables; ++j) {
            if (!is_pivot_col[j]) {
                solution.particular_solution.check_idx(col) -= coeff_matrix.row(i).check_idx(j) * solution.particular_solution.check_idx(j);
            }
        }
    }

    solution.basis_dimension = 0;
    for (size_t j = 0; j < variables; ++j) {
        if (!is_pivot_col[j]) ++solution.basis_dimension;
    }

    if (solution.basis_dimension > 0) {
        solution.basis_vectors = new Vector<NumType>[solution.basis_dimension];
        solution.free_variables = new size_t[solution.basis_dimension];
        size_t basis_idx = 0;

        for (size_t free_col = 0; free_col < variables; ++free_col) {
            if (!is_pivot_col[free_col]) {
                solution.free_variables[basis_idx] = free_col;
                solution.basis_vectors[basis_idx] = Vector<NumType>(variables);
                solution.basis_vectors[basis_idx].check_idx(free_col) = 1;

                for (ptrdiff_t i = rank - 1; i >= 0; --i) {
                    size_t pivot_col = pivot_columns[i];
                    if (pivot_col >= variables) continue;

                    solution.basis_vectors[basis_idx].check_idx(pivot_col) = 0;
                    for (size_t j = pivot_col; j < variables; ++j) {
                        solution.basis_vectors[basis_idx].check_idx(pivot_col) -= coeff_matrix.row(i).check_idx(j) * solution.basis_vectors[basis_idx].check_idx(j);
                    }
                }
                ++basis_idx;
            }
        }
    }

    delete[] is_pivot_col;
    delete[] pivot_columns;
    return solution;
}

template <typename NumType>
void VerifySolution(const Matrix<NumType>& coeff_matrix, const Vector<NumType>& const_terms, const Solution<NumType>& solution) {
    const NumType tolerance = numeric_limits<NumType>::epsilon() * static_cast<NumType>(coeff_matrix.rows()) * 1000;
    if (!solution.is_consistent) {
        cout << "System has no solutions (inconsistent)\n";
        return;
    }
    for (size_t i = 0; i < coeff_matrix.rows(); ++i) {
        NumType sum = 0;
        for (size_t j = 0; j < coeff_matrix.cols(); ++j) {
            sum += coeff_matrix.row(i).check_idx(j) * solution.particular_solution.check_idx(j);
        }
        NumType error = std::abs(sum - const_terms.check_idx(i));
        NumType max_val = max(NumType(1), max(std::abs(sum), std::abs(const_terms.check_idx(i))));
        if (error > tolerance * max_val) {
            throw runtime_error("Solution verification failed");
        }
    }
    cout << "Solution verified successfully\n";
}

template <typename NumType>
void DisplaySolution(const Solution<NumType>& solution) {
    const NumType tolerance = numeric_limits<NumType>::epsilon() * 100;

    if (!solution.is_consistent) {
        cout << "System has no solutions\n";
        return;
    }

    cout << "Solution components:\n";

    if (solution.is_unique) {
        for (size_t i = 0; i < solution.particular_solution.size(); ++i) {
            cout << "x" << i + 1 << " = " << solution.particular_solution.check_idx(i) << "\n";
        }
    }
    else {
        for (size_t i = 0; i < solution.particular_solution.size(); ++i) {
            cout << "x" << i + 1 << " = " << solution.particular_solution.check_idx(i);

            for (size_t k = 0; k < solution.basis_dimension; ++k) {
                if (std::abs(solution.basis_vectors[k].check_idx(i)) > tolerance) {
                    cout << " + " << solution.basis_vectors[k].check_idx(i) << "*a" << k + 1;
                }
            }
            cout << "\n";
        }

        for (size_t k = 0; k < solution.basis_dimension; ++k) {
            cout << "a" << k + 1 << " is free (x" << solution.free_variables[k] + 1 << ")\n";
        }
    }
}

template<typename T>
void HandleRandomInput() {
    size_t system_size;
    cout << "Enter system size (n x n): ";
    cin >> system_size;

    if (system_size == 0) {
        throw invalid_argument("System size must be positive");
    }

    random_device rd;
    mt19937 gen(rd());
    T value_range = static_cast<T>(10000.0) / sqrt(static_cast<T>(system_size));
    uniform_real_distribution<T> dist(-value_range, value_range);

    Matrix<T> coefficients(system_size, system_size);
    Vector<T> terms(system_size);

    for (size_t i = 0; i < system_size; ++i) {
        T row_sum = 0;
        for (size_t j = 0; j < system_size; ++j) {
            if (i != j) {
                coefficients.row(i).check_idx(j) = dist(gen);
                row_sum += std::abs(coefficients.row(i).check_idx(j));
            }
        }
        coefficients.row(i).check_idx(i) = row_sum * (1 + dist(gen) / static_cast<T>(1000.0));
        terms.check_idx(i) = dist(gen);
    }

    auto start_time = steady_clock::now();
    Solution<T> solution = SolveLinearSystem(coefficients, terms);
    auto end_time = steady_clock::now();
    duration<double> elapsed = end_time - start_time;

    VerifySolution(coefficients, terms, solution);
    cout << "Computation time: " << elapsed.count() << " seconds\n";
}

template<typename T>
void HandleManualInput() {
    size_t equations, variables;
    cout << "Enter number of equations: ";
    cin >> equations;
    cout << "Enter number of variables: ";
    cin >> variables;

    if (equations == 0 || variables == 0) {
        throw invalid_argument("Dimensions must be positive");
    }

    Matrix<T> coefficients(equations, variables);
    Vector<T> terms(equations);

    cout << "Enter coefficient matrix and constant terms:\n";
    for (size_t i = 0; i < equations; ++i) {
        cout << "Equation " << i + 1 << ": ";
        for (size_t j = 0; j < variables; ++j) {
            if (!(cin >> coefficients.row(i).check_idx(j))) {
                throw invalid_argument("Invalid numeric input");
            }
        }
        if (!(cin >> terms.check_idx(i))) {
            throw invalid_argument("Invalid numeric input");
        }
    }

    Solution<T> solution = SolveLinearSystem(coefficients, terms);
    VerifySolution(coefficients, terms, solution);
    DisplaySolution(solution);
}

int main() {
    try {
        cout << "Select numeric precision:\n";
        cout << "1. float\n";
        cout << "2. double\n";
        cout << "3. exit\n";
        cout << "Enter choice: ";
        int precision_choice;
        cin >> precision_choice;

        if (precision_choice == 3) {
            cout << "Program terminated\n";
            return 0;
        }

        cout << "Select input method:\n";
        cout << "1. Generate random system\n";
        cout << "2. Enter system manually\n";
        cout << "3. Exit\n";
        cout << "Enter choice: ";
        int input_choice;
        cin >> input_choice;

        if (input_choice == 3) {
            cout << "Program terminated\n";
            return 0;
        }

        if (precision_choice == 1) {
            if (input_choice == 1) {
                HandleRandomInput<float>();
            }
            else if (input_choice == 2) {
                HandleManualInput<float>();
            }
            else {
                throw invalid_argument("Invalid input method selection");
            }
        }
        else if (precision_choice == 2) {
            if (input_choice == 1) {
                HandleRandomInput<double>();
            }
            else if (input_choice == 2) {
                HandleManualInput<double>();
            }
            else {
                throw invalid_argument("Invalid input method selection");
            }
        }
        else {
            throw invalid_argument("Invalid precision selection");
        }
    }
    catch (const exception& e) {
        cerr << "Error occurred: " << e.what() << endl;
        return 1;
    }
    _CrtDumpMemoryLeaks();
    return 0;
}
