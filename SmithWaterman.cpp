#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>


// Matrix class but implemented as linear vector due to performance considerations.
// Note that SW access things in column major order. We will use the mat[r][c] convention and internally organize the vector so that accessing the next element in the row is contiguous.
template <class T>
class Matrix {
    struct MatrixRow; // forward declare this to be used by operator[]

public:
    Matrix(size_t num_rows, size_t num_cols, const T & value) : num_rows_(num_rows), num_cols_(num_cols), vec_(std::vector<T>(num_rows * num_cols, value)) {}

    MatrixRow operator[] (size_t row_index) { return MatrixRow(row_index, *this); }

    size_t GetNumRows() { return num_rows_; }
    size_t GetNumCols() { return num_cols_; }

private:
    struct MatrixRow {
        Matrix & mat;
        size_t row_index;

        MatrixRow(size_t row_index, Matrix & mat) : mat(mat), row_index(row_index) {}

        T& operator[] (size_t col_index) {
            return mat.vec_[row_index * mat.num_cols_ + col_index];
        }
    };

    size_t num_rows_;
    size_t num_cols_;
    
    std::vector<T> vec_;
};

int main() {
    using DataType = int;

    DataType match = 5;
    DataType mismatch = -3;
    DataType gap_start_penalty = -8;
    DataType gap_extend_penalty = -1;

//    std::string seq1 = "CAGCCUCGCUUAG";
//    std::string seq2 = "AAUGCCAUUGCCGG";

    std::string seq1 = "CAGCCUCGCUUAG";
    std::string seq2 = "AAUGCCAUUGCCGG";

//    std::string seq1 = "TGTTACGG";
//    std::string seq2 = "GGTTGACTA";

    Matrix<DataType> e_mat(seq2.size() + 1, seq1.size() + 1, 0);
    Matrix<DataType> f_mat(seq2.size() + 1, seq1.size() + 1, 0);
    Matrix<DataType> h_mat(seq2.size() + 1, seq1.size() + 1, 0);
    Matrix<DataType> h_hat_mat(seq2.size() + 1, seq1.size() + 1, 0);


//// 1
//    for (int r = 1; r < e_mat.GetNumRows(); ++r) {
//        for (int c = 1; c < e_mat.GetNumCols(); ++c) {
//            auto S_ij = seq2.at(r-1) == seq1.at(c-1) ? match : mismatch;
//            auto residue_1 = std::max(0, h_mat[r-1][c-1] + S_ij);
//            auto residue_2_func = [&](){
//                auto best = 0;
//                for (int k = 0; k < c; ++k) {
//                    auto cur = h_mat[r][c-k] - gap_start_penalty - k * gap_extend_penalty;
//                    if (cur > best)
//                        best = cur;
//                }
//
//                return best;
//            };
//
//            auto residue_2 = residue_2_func();
//
//            auto residue_3_func = [&](){
//                auto best = 0;
//                for (int k = 0; k < r; ++k) {
//                    auto cur = h_mat[r-k][c] - gap_start_penalty - k * gap_extend_penalty;
//                    if (cur > best)
//                        best = cur;
//                }
//
//                return best;
//            };
//
//            auto residue_3 = residue_3_func();
//
//            h_mat[r][c] = std::max(std::max(residue_1, residue_2), residue_3);
//        }
//    }

//// 2
//    for (int r = 1; r < e_mat.GetNumRows(); ++r) {
//        for (int c = 1; c < e_mat.GetNumCols(); ++c) {
//            e_mat[r][c] = std::max(e_mat[r-1][c] - gap_extend_penalty, h_mat[r-1][c] - gap_start_penalty);
//            f_mat[r][c] = std::max(e_mat[r][c-1] - gap_extend_penalty, h_mat[r][c-1] - gap_start_penalty);
//
//            auto best_h = [&](){
//                auto S_ij = seq2.at(r-1) == seq1.at(c-1) ? match : mismatch;
//
//                auto best_a = std::max(h_mat[r-1][c-1] + S_ij, e_mat[r][c]);
//                auto best_b = std::max(f_mat[r][c], 0);
//
//                return std::max(best_a, best_b);
//            };
//
//            h_mat[r][c] = best_h();
//        }
//    }

//// 3
//    for (int r = 1; r < e_mat.GetNumRows(); ++r) {
//        for (int c = 1; c < e_mat.GetNumCols(); ++c) {
//            f_mat[r][c] = std::max(f_mat[r-1][c], h_mat[r-1][c] + gap_start_penalty) + gap_extend_penalty;
//
//            auto S_ij = seq2.at(r-1) == seq1.at(c-1) ? match : mismatch;
//            h_hat_mat[r][c] = std::max(std::max(h_mat[r-1][c-1] + S_ij, f_mat[r][c]), 0);
//
//            auto best_e = [&]() {
//                auto best = 0;
//                for (int k = 1; k < c; ++k) {
//                    auto cur = h_hat_mat[r][c-k] + k * gap_extend_penalty;
//                    if (cur > best) {
//                        best = cur;
//                    }
//                }
//
//                return best;
//            };
//
//            e_mat[r][c] = best_e();
//            h_mat[r][c] = std::max(h_hat_mat[r][c], e_mat[r][c] + gap_start_penalty);
//        }
//    }

// 4
    // Calculate padded row size, this will be useful later when calculating e_mat
    size_t padded_row_size = e_mat.GetNumCols();

    // bit twiddling to round up the row size to next power of two
    --padded_row_size;
    padded_row_size |= padded_row_size >> 1;
    padded_row_size |= padded_row_size >> 2;
    padded_row_size |= padded_row_size >> 4;
    padded_row_size |= padded_row_size >> 8;
    padded_row_size |= padded_row_size >> 16;
    padded_row_size |= padded_row_size >> 32;
    ++padded_row_size;
    std::cout << "Row size rounded up to next power of 2: " << padded_row_size << std::endl;

    for (int r = 1; r < e_mat.GetNumRows(); ++r) {
        for (int c = 1; c < e_mat.GetNumCols(); ++c) {
            f_mat[r][c] = std::max(f_mat[r-1][c], h_mat[r-1][c] + gap_start_penalty) + gap_extend_penalty;

            auto S_ij = seq2.at(r-1) == seq1.at(c-1) ? match : mismatch;
            h_hat_mat[r][c] = std::max(std::max(h_mat[r-1][c-1] + S_ij, f_mat[r][c]), 0);
        }

        // Calculate e_mat using modified parallel scan
        {
            auto pow_of_2 = [](const size_t pow)
            {
                return 1 << pow;
            };


            auto log2 = [](size_t num) { // This has no branching on GPU

                assert(num != 0);

                size_t log = 0;

                while (num != 0) {
                    num = num >> 1;
                    ++log;
                }

                return log-1;
            };

            // Pad structure to power of two
            std::vector<DataType> padded_row(padded_row_size);

            // Initialize padded_row
            for (int c = 0; c < e_mat.GetNumCols(); ++c) {
                padded_row[c] = h_hat_mat[r][c];
            }

//            std::cout << "Original padded row:\t";
//            for (size_t k = 0; k < padded_row_size; ++k) {
//                std::cout << padded_row[k] << "\t";
//            }
//            std::cout << std::endl;

//            // Upsweep (inner loop can be done in parallel)
//            for (size_t d = 0; d < log2(padded_row_size); ++d) {
//                for (size_t k = 0; k < padded_row_size; k += pow_of_2(d+1)) {
//                    //std::cout << k + pow_of_2(d) - 1 << " " << k + pow_of_2(d+1) - 1 << std::endl;
//                    int64_t left_elem = padded_row[k + pow_of_2(d) - 1];
//                    int64_t right_elem = padded_row[k + pow_of_2(d+1) - 1] - (pow_of_2(d) * gap_extend_penalty);
//                    //std::cout << left_elem << " " << right_elem << std::endl;
//                    padded_row[k + pow_of_2(d+1) - 1] = std::max(left_elem, right_elem);
//                }
//
//                std::cout << "Upswept padded row " << d << ":\t";
//                for (size_t k = 0; k < padded_row_size; ++k) {
//                    std::cout << padded_row[k] << "\t";
//                }
//                std::cout << std::endl;
//            }
//
//            // Downsweep (inner loop can be done in parallel)
//            padded_row[padded_row_size-1] = 0;
//            for (int64_t d = log2(padded_row_size) - 1; d >= 0; --d) {
//                for (size_t k = 0; k < padded_row_size; k += pow_of_2(d+1)) {
//                    //std::cout << k + pow_of_2(d) - 1 << " " << k + pow_of_2(d+1) - 1 << std::endl;
//                    DataType temp = padded_row[k + pow_of_2(d) - 1];
//                    padded_row[k + pow_of_2(d) - 1] = padded_row[k + pow_of_2(d+1) - 1];
//                    int64_t left_elem = temp;
//                    int64_t right_elem = padded_row[k + pow_of_2(d+1) - 1];
//                    padded_row[k + pow_of_2(d+1) - 1] = std::max(left_elem, right_elem) + (pow_of_2(d) * gap_extend_penalty);
//                }
//
//                std::cout << "Downswept padded row " << d << ":\t";
//                for (size_t k = 0; k < padded_row_size; ++k) {
//                    std::cout << padded_row[k] << "\t";
//                }
//                std::cout << std::endl;
//            }
//            std::cout << std::endl;

            // Upsweep (inner loop can be done in parallel)
            for (size_t d = 0; d < log2(padded_row_size); ++d) {
                for (size_t k = 0; k < padded_row_size / pow_of_2(d+1); ++k) {
                    size_t z = k * pow_of_2(d+1);
                    int64_t left_elem = padded_row[z + pow_of_2(d) - 1];
                    int64_t right_elem = padded_row[z + pow_of_2(d+1) - 1] - (pow_of_2(d) * gap_extend_penalty);
                    padded_row[z + pow_of_2(d+1) - 1] = std::max(left_elem, right_elem);
                }

//                std::cout << "Upswept padded row " << d << ":\t";
//                for (size_t k = 0; k < padded_row_size; ++k) {
//                    std::cout << padded_row[k] << "\t";
//                }
//                std::cout << std::endl;
            }

            std::cout << "r " << r << ": ";
            for (int k = 0; k < padded_row_size ; ++k) {
                std::cout << padded_row[k] << "\t";
            }
            std::cout << std::endl;

            // Downsweep (inner loop can be done in parallel)
            padded_row[padded_row_size-1] = 0;
            for (int64_t d = log2(padded_row_size) - 1; d >= 0; --d) {
                for (size_t k = 0; k < padded_row_size / pow_of_2(d+1); ++k) {
                    size_t z = k * pow_of_2(d+1);
                    DataType temp = padded_row[z + pow_of_2(d) - 1];
                    padded_row[z + pow_of_2(d) - 1] = padded_row[z + pow_of_2(d+1) - 1];
                    int64_t left_elem = temp;
                    int64_t right_elem = padded_row[z + pow_of_2(d+1) - 1];
                    padded_row[z + pow_of_2(d+1) - 1] = std::max(left_elem, right_elem) + (pow_of_2(d) * gap_extend_penalty);
                }

//                std::cout << "Downswept padded row " << d << ":\t";
//                for (size_t k = 0; k < padded_row_size; ++k) {
//                    std::cout << padded_row[k] << "\t";
//                }
//                std::cout << std::endl;
            }
            //std::cout << std::endl;

            // copy padded_row back into e_mat
            for (int c = 0; c < e_mat.GetNumCols(); ++c) {
                e_mat[r][c] = padded_row[c];
            }
        }

        // Calculate h_mat
        for (int c = 0; c < e_mat.GetNumCols(); ++c) {
            //std::cout << r << " " << c << " " << h_hat_mat[r][c] << " " << e_mat[r][c] - gap_start_penalty << std::endl;
            h_mat[r][c] = std::max(h_hat_mat[r][c], e_mat[r][c] + gap_start_penalty);
        }
    }

//    for (int r = 0; r < e_mat.GetNumRows(); ++r) {
//        for (int c = 0; c < e_mat.GetNumCols(); ++c) {
//            std::cout << h_hat_mat[r][c] << "\t";
//        }
//        std::cout << std::endl;
//    }
//    std::cout << std::endl;
//
//    for (int r = 0; r < e_mat.GetNumRows(); ++r) {
//        for (int c = 0; c < e_mat.GetNumCols(); ++c) {
//            std::cout << e_mat[r][c] << "\t";
//        }
//        std::cout << std::endl;
//    }
//    std::cout << std::endl;
//
    for (int r = 0; r < e_mat.GetNumRows(); ++r) {
        for (int c = 0; c < e_mat.GetNumCols(); ++c) {
            std::cout << h_mat[r][c] << "\t";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;

    std::cout << h_mat[h_mat.GetNumRows()-1][h_mat.GetNumCols()-1] << std::endl;

    return 0;
}
