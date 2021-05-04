#include "matrix.h"
#include "timer.h"

#include "mpi.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <unistd.h>
#include <algorithm>

#define OUT_TO_FILE 0

const int N1 = 260;
const int N2 = 168;
const int N3 = 159;

Matrix concat_matrices_by_cols(std::vector<Matrix> row_vec);

int main(int argc, char* argv[])
{
    // init static vars for output
    Matrix::out_width = 9;
    Matrix::out_precision = 1;

	int ProcNum, ProcRank, RecvRank;
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    bool isRandGen = false;

    if (ProcNum != 8) {
        if (ProcRank == 0)
            std::cout << "Need 8 proc num. terminate\n";
        MPI_Finalize();
        return -1;
    }


    const std::vector<int> hamilt_vec {0, 7, 6, 5, 2, 4, 3, 1};

    std::vector<int> submatrix_rows { 32, 32, 32, 32, 33, 33, 33, 33};
    std::vector<int> submatrix_result = submatrix_rows;
    std::vector<int> submatrix_displs_rows(ProcNum);
    std::vector<int> submatrix_displs_result(ProcNum);
    for (int i = 0; i < ProcNum; i++) {
        submatrix_rows[i] *= N2;
        submatrix_result[i] *= N3;
        submatrix_displs_rows[i] = i > 0? submatrix_displs_rows[i-1] + submatrix_rows[i-1] : 0;
        submatrix_displs_result[i] = i > 0? submatrix_displs_result[i-1] + submatrix_result[i-1] : 0;
    }

    std::vector<int> submatrix_cols { 20, 20, 20, 20, 20, 20, 20, 19};
    std::vector<int> submatrix_displs_cols(ProcNum);
    for (int i = 0; i < ProcNum; i++) {
        submatrix_cols[i] *= N2;
        submatrix_displs_cols[i] = i > 0? submatrix_displs_cols[i-1] + submatrix_cols[i-1] : 0;
    }

    custom_timer::Timer timer_seq,
                        timer_mpi;

    Matrix C_row;

    std::vector<Matrix> C_row_vec(ProcNum);
    Matrix A_temp(submatrix_rows[ProcRank] / N2, N2);
    Matrix B_temp(submatrix_cols[ProcRank] / N2, N2);


	if (ProcRank == 0 ) {

        Matrix A(N1, N2),
            B(N2, N3);

        // Variant info
        std::cout << std::left << std::setw(40) << "Hamilt loop:";
        for (auto elem : hamilt_vec) {
            std::cout << std::setw(5) << elem;
        }
        std::cout << "\n\nN1 = " << N1 << "; N2 = " << N2 <<
            std::left << std::setw(40) << "\nSubarray elem (rows)rows counts: ";
        for(auto elem : submatrix_rows) {
            std::cout << "(" << elem / N2 << ")" << std::left << std::setw(5) << elem;
        }
        std::cout << std::left << std::setw(40) << "\nSubarray data displacements: ";
        for(auto elem : submatrix_displs_rows) {
            std::cout << std::setw(9) << elem;
        }
        std::cout << "\n\nN2 = " << N2 << "; N3 = " << N3 <<
            std::left << std::setw(40) << "\nSubarray elem (cols)col counts: ";
        for(auto elem : submatrix_cols) {
            std::cout << "(" << elem / N2 << ")" << std::left << std::setw(5) << elem;
        }
        std::cout << std::left << std::setw(40) << "\nSubarray data displacements: ";
        for(auto elem : submatrix_displs_cols) {
            std::cout  << std::setw(9) << elem;
        }
        std::cout << "\n";
        // ------

        std::cout << "Generate random matrices '(y)es (n)o:  ";
        std::string ans;
        std::cin >> ans;
        std::for_each(ans.begin(), ans.end(), tolower);
        isRandGen = (ans == "yes" || ans == "y") ? true : false;

        if(isRandGen) {
            A.init_matrix(0, 20);
            B.init_matrix(0, 20);
        }
        else {
            std::cout << "A(" << A.get_row() << ", " << A.get_col() << ")=\n";
            std::cin >> A;
            std::cout << "B(" << B.get_row() << ", " << B.get_col() << ")=\n";
            std::cin >> B;
        }

        // ==================================================================================
        // Start Sequantial execution
        // ==================================================================================

        std::cout << "\n### Execute sequantial algo ###\n";
        std::cout << "- Start Seq Timer -\n";

        timer_seq.start();
        Matrix C_seq =  A * B;
        timer_seq.stop();

        std::cout << "- Exec ended. Write result ot seq_result.txt -\n";

        std::ofstream fout("input_matrices.txt");
        fout << "A=\n" << A
                << "B=\n" << B;
        fout.close();
        fout.open("seq_result.txt");
        fout << "C=\n" << C_seq;
        fout.close();

        std::cout << "### Seq Timer duration: " << std::fixed << timer_seq.getDuration() << " ms\n" << std::endl;

        // ==================================================================================
        // Start Parallel execution
        // ==================================================================================

        // used transposed matrix B for comfort use of Scatterv
        Matrix Bt = B.get_transpose();

        std::cout << "### Execute parallel algo ###\n";
        std::cout << "- Start Mpi Timer -\n";

        // Start algo after barrier
        MPI_Barrier(MPI_COMM_WORLD);
        timer_mpi.start();

        // Scatter A rows to all procs
        MPI_Scatterv(A.data(), submatrix_rows.data(), submatrix_displs_rows.data(), MPI_DOUBLE, A_temp.data(), submatrix_rows[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // Scatter B columns to all procs, use transposed B for send plain data.
        MPI_Scatterv(Bt.data(), submatrix_cols.data(), submatrix_displs_cols.data(), MPI_DOUBLE, B_temp.data(), submatrix_cols[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    }
	else {
        // wait for init input matrices
        MPI_Barrier(MPI_COMM_WORLD);

        // Recieve A(ProcRank) submatrix from Master
        MPI_Scatterv(NULL, submatrix_rows.data(), submatrix_displs_rows.data(), MPI_DOUBLE, A_temp.data(), submatrix_rows[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // Recieve B(ProcRank) submatrix from Master
        MPI_Scatterv(NULL, submatrix_cols.data(), submatrix_displs_cols.data(), MPI_DOUBLE, B_temp.data(), submatrix_cols[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Evaluate self row*col (C_row[ProcRank][ProcRank])
    C_row_vec[ProcRank] = M_x_Mt(A_temp, B_temp);
#if OUT_TO_FILE
        std::ofstream fout_first("result/" + std::to_string(ProcRank) + "_first_state.txt");
        fout_first << "A=\n" << A_temp << "B=\n" << B_temp.get_transpose();
        fout_first << "C(" << ProcRank << ")=\n" << C_row_vec[ProcRank];
        fout_first.close();
#endif

    // Find self postion in hamilt cycle and found next and previous proc regarding to hamilt
    int pos = std::find(hamilt_vec.begin(), hamilt_vec.end(), ProcRank) - hamilt_vec.begin(); // get position if current Proc in hamilt_vec
    int next_proc = hamilt_vec[(pos + ProcNum + 1) % ProcNum]; // get ID of next proc
    int prev_proc = hamilt_vec[(pos + ProcNum - 1) % ProcNum]; // get ID of prev proc

    int submatrix_index = ProcRank;
    for (int i = 1; i < ProcNum; i++) {
        int recv_submatrix_index = hamilt_vec[(pos + ProcNum - i) % ProcNum];
        Matrix B_recv(submatrix_cols[recv_submatrix_index] / N2, N2);

        // tag field keeps absolute porsition of submatrix in final row
        MPI_Sendrecv(
            B_temp.data(), submatrix_cols[submatrix_index], MPI_DOUBLE, next_proc, submatrix_index,                       // send proc's saved submatrix
            B_recv.data(), submatrix_cols[recv_submatrix_index], MPI_DOUBLE, prev_proc, recv_submatrix_index, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recv prev proc submatrix

        // Evaluate subresult
        C_row_vec[recv_submatrix_index] = M_x_Mt(A_temp, B_recv);

        B_temp = std::move(B_recv);
        submatrix_index = recv_submatrix_index;
    }

#if OUT_TO_FILE
    std::ofstream fout_temp("result/" + std::to_string(ProcRank) + "_C_row_res.txt");
    fout_temp << "A=\n" << A_temp;
    for (int i = 0; i < C_row_vec.size(); i++) {
        fout_temp << "\nC(" << i << ")=\n" << C_row_vec[i];
    }
    fout_temp.close();
#endif

    // Concat matrices
    C_row = concat_matrices_by_cols(C_row_vec);

    if (ProcRank == 0) {
        Matrix C(N1, N3);

        // Gather rows to Master
        MPI_Gatherv(C_row.data(), submatrix_result[ProcRank], MPI_DOUBLE, C.data(), submatrix_result.data(), submatrix_displs_result.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        timer_mpi.stop();
        std::cout << "- Exec ended. Write result to mpi_result.txt -\n";

        // Write result
        std::ofstream fout_seq("mpi_result.txt");
        fout_seq << "C=\n" << C;
        fout_seq.close();

        std::cout << "### Parallel Timer duration: " << std::fixed << timer_mpi.getDuration() << " ms" << std::endl;
    }
    else {
        MPI_Gatherv(C_row.data(), submatrix_result[ProcRank], MPI_DOUBLE, NULL, submatrix_result.data(), submatrix_displs_result.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

	MPI_Finalize();
	return 0;
}


Matrix concat_matrices_by_cols(std::vector<Matrix> row_vec) {
    Matrix C_row(row_vec[0].get_row(), N3);
    // iterate through all submatrices
    for (int n = 0, disp = 0; n < row_vec.size(); n++) {
        // copy each row to each position
        for (int r = 0; r < row_vec[n].get_row(); r++) // 3
            std::copy(
                row_vec[n].data() + r*row_vec[n].get_col(),
                row_vec[n].data() + r*row_vec[n].get_col() + row_vec[n].get_col(),
                C_row.data() + r * N3 + disp);
        disp += row_vec[n].get_col();
    }
    return C_row;
}