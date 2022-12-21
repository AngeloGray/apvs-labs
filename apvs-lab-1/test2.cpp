#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <vector>
using namespace std;

    //������� ��� ������ ������ � �������� � �������
    void ShowMatrix(vector<vector<float> >& a) {
    // ����� �������
        for (int i = 0; i < a.size(); i++) {
            for (int j = 0; j < a.size(); j++) {
                cout << a[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "\n";
    }
    void ShowMatrixEx(vector<vector<float> >& a) {
        for (int i = 0; i < a.size(); i++) {
            for (int j = 0; j < a.size() + 1; j++) {
                if (a.size() - j == 1) {
                    cout << a[i][j] << " | "; // ����������� ��������
                    continue;                      // � ����������� ������
                }
                cout << a[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "\n";
    }
    void ShowVector(vector<float> & a){
        //����� �������
        for (int i = 0; i < a.size(); i++) {
            cout << a[i] << "\n";
        }
        cout << "\n";
    }
    void ShowResult(vector<vector<float>>& a) {
        cout << "Linear equations system answer is:" << endl;
        for (int i = 0; i < a.size(); i++) {
            cout << "X" << i << " = " << a[i][a.size()] << "\n";
        }
        cout << "\n";
    }

    // ������� ��� ������� ������� ������-�������
    void RowElimination(vector<float>& row, vector<float>& ref_row, float element) {

        for (int i = 0; i < row.size(); i++) {
            row[i] = row[i] - ref_row[i] * element;
        }
    }
    void RowDivision(vector<float>& row, float element) {
        for (int i = 0; i < row.size(); i++) {
            row[i] = row[i] / element;
        }
    }
    void SolveMatrix(vector<vector<float> >& a, double delta) {
        for (int k = 0; k < a.size(); k++) {
            if (abs(a[k][k]) < delta) {                  // ��������.
                for (int i = k + 1; i < a.size(); i++) { // ���� �������
                    if (a[i][k] > delta) {               // ����� ����, ��
                        swap(a[i], a[k]);           // �������� �������
                        break;                           // ������ � ���������
                    }
                }
            }
            // ���� ������� �� ����� ������� - ����� ������ �� ����
            if (a[k][k] < (1 - delta) || a[k][k] > (1 + delta)) {
                RowDivision(a[k], a[k][k]);
            }
            // ������ ��� ��������� ��������� �������
            // ������ �� ���������� ���� ����� �������, ���������� �� �������

            for (int i = k + 1; i < a.size(); i++) {
                RowElimination(a[i], a[k], a[i][k]);
            }
        }
        // �������� ��� ��������� ��������� �������

        for (int k = a.size() - 1; k >= 0; k--) {
            for (int i = k - 1; i >= 0; i--) { // �������� ��� ����������� �����
                RowElimination(a[i], a[k], a[i][k]);
            }
        }
    }

    //������� ��� ��������� � �������������� ������� A
    void GenMatA(vector<vector <float>> &a_matrix){
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank); // ���� �������� ��������
        // �������� ������� �
        vector<vector<float>> a_matrix_v1(a_matrix.size());
        // ���������� ������� ���������� ���������� (�����������)
        for (int i = 0; i < a_matrix.size(); i++) {
            for (int j = 0; j < a_matrix.size(); j++) {
                a_matrix_v1[i].push_back(rand() % 10); // ��������� ������ ����������
            }                                 // a[i][j] = {0, 9}
        }
        if (rank == 0) ShowMatrix(a_matrix_v1);
        // ���������������� �������
        vector<vector<float>> a_matrix_v2(a_matrix.size(), vector <float>(a_matrix.size()));
        for (int i = 0; i < a_matrix.size(); i++) {
            for (int j = 0; j < a_matrix.size(); j++) {
                a_matrix_v2[i][j] = a_matrix_v1[j][i];
            }
        }
        if (rank == 0)ShowMatrix(a_matrix_v2);
        float result = 0;
        //Fill each cell of the matrix output.
        for (int i = 0; i < a_matrix.size(); i++) {
            for (int j = 0; j < a_matrix.size(); j++) {
                //Multiply each row of matrix 1 with each column of matrix 2.
                for (int k = 0; k < a_matrix.size(); k++) {
                    result += a_matrix_v1[i][k] * a_matrix_v2[k][j];
                }
                a_matrix[i][j] = result;
                result = 0; //Reset;
            }
        }
        if (rank == 0)cout << "Matrix A:\n";
        if (rank == 0)ShowMatrix(a_matrix);
    }

 
    //������� ��� ���������� ������� L ������� ���������
  
    void FindMatLSerial(vector<vector<float>> &l_matrix) {
        int numProc, rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank); // ���� �������� ��������
        MPI_Comm_size(MPI_COMM_WORLD, &numProc); // ���-�� ���������
        if (rank == 0) {
            // ���������� ������� L ������� ����������
            for (int j = 0; j < l_matrix.size(); j++) {
                for (int k = 0; k < j; k++) {
                    for (int i = j; i < l_matrix.size(); i++) {
                        l_matrix[i][j] = l_matrix[i][j] - l_matrix[i][k] * l_matrix[j][k];
                    }
                }
                l_matrix[j][j] = sqrt(l_matrix[j][j]);
                for (int i = j + 1; i < l_matrix.size(); i++) {
                    l_matrix[i][j] = l_matrix[i][j] / l_matrix[j][j];
                }
            }
            cout << "Serial Cholesky decomposition of matrix A: L= \n";
            ShowMatrix(l_matrix); // �������� ������ ����������� �������
        }
    }


    //������� ��� ���������� ������� L ������� ���������
   /*
    void FindMatL(vector<vector<float>>& l_matrix) {
        int numProc, rank;
        MPI_Comm_size(MPI_COMM_WORLD, &numProc);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        double start, end;
        MPI_Barrier(MPI_COMM_WORLD); 
        if (rank == 0) {
            start = MPI_Wtime();
        }

        // For each column
        int i, j, k;
        for (j = 0; j < l_matrix.size(); j++) {
  
             * Step 1:
             * Update the diagonal element


            if (j % numProc == rank) {

                for (k = 0; k < j; k++) {
                    l_matrix[j][j] = l_matrix[j][j] - l_matrix[j][k] * l_matrix[j][k];
                }

                l_matrix[j][j] = sqrt(l_matrix[j][j]);
            }


            // Broadcast row with new values to other processes
            MPI_Bcast(l_matrix[j], l_matrix.size(), MPI_FLOAT, j % numProc, MPI_COMM_WORLD);

     
             * Step 2:
             * Update the elements below the diagonal element
      

             // Divide the rest of the work
            for (i = j + 1; i < l_matrix.size(); i++) {
                if (i % numProc == rank) {
                    for (k = 0; k < j; k++) {
                        l_matrix[i][j] = l_matrix[i][j] - l_matrix[i][k] * l_matrix[j][k];
                    }

                    l_matrix[i][j] = l_matrix[i][j] / l_matrix[j][j];
                }
            }
        }


        MPI_Barrier(MPI_COMM_WORLD); 
        if (rank == 0) {
            end = MPI_Wtime();
            printf("Testing OpenMpi implementation Output: \n");
            printf("Runtime = %lf\n", end - start);
        }
    }
    */

    //������� ��� ���������� ������� L ������� ���������
    void FindMatLParallel(vector <vector <float>>& M) {
        int world_rank, world_size, n = M.size();
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        for (int k = 0; k < n; k++)//k-column index (algorithm goes below of diagonal from left to right)
        {
            if (world_rank == 0) {
                //0 procesor calculates every diagonal element
                for (int j = 0; j < k; j++)
                {
                    M[k][k] -= (M[k][j] * M[k][j]);
                }
                M[k][k] = sqrt(M[k][k]);

                for (int p = 1; p < n - k; p++) //p-processor index to which 0 processor sends a given task
                {
                    for (int c = 0; c <= k; c++)//processor 0 sends already calculated values
                    {
                        for (int r = c; r < n; r++)
                        {
                            MPI_Send(&M[r][c], 1, MPI_FLOAT, p, 0, MPI_COMM_WORLD);
                        }
                    }
                    MPI_Recv(&M[p + k][k], 1, MPI_FLOAT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                if (k == n - 1){

                    //print result to the screen
                    cout <<"Parallel Cholesky decomposition of matrix A: L= \n";
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                            cout << M[i][j];
                            if (j != n - 1)
                            {
                                cout << ", ";
                            }
                        }
                        cout << endl;
                    }
                }
            }
            else if (world_rank < n - k) {
                for (int c = 0; c <= k; c++)//receving already calculated values
                {
                    for (int r = c; r < n; r++)
                    {
                        MPI_Recv(&M[r][c], 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }

                for (int g = 0; g < k; g++)//calculating non-diagonal elements concurrently
                {
                    M[k + world_rank][k] -= M[k][g] * M[k + world_rank][g];
                }
                M[k + world_rank][k] /= M[k][k];
                MPI_Send(&M[k + world_rank][k], 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
            }
        }
    }

// �������� ������� ���������.
void doKholesky(int argc, char** argv, int a_size) {
    int numProc, rank;
    MPI_Init(&argc, &argv); // MPI-�������������
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // ���� �������� ��������
    MPI_Comm_size(MPI_COMM_WORLD, &numProc); // ���-�� ���������
    
    //�������� ������������ ������������ ����������� ������� �
    vector<vector<float>> a_matrix(a_size, vector <float>(a_size));
    GenMatA(a_matrix);

    // ����������� � ������� L ������ ����������� ������� �� �
    // �������� ���������� ��� ���������� ��������� - ����� ������� L
    vector<vector<float>> l_matrix(a_size, vector <float>(a_size));
    vector<vector<float>> l_matrix_copy(a_size, vector <float>(a_size));
        for (int i = 0; i < a_size; i++) {
            for (int j = 0; j <= i; j++) {
                l_matrix[i][j] = a_matrix[i][j];
                l_matrix_copy[i][j] = a_matrix[i][j];
            }
        }

        if (rank == 0) {
            ShowMatrix(l_matrix);
            double startTime = MPI_Wtime();
            FindMatLSerial(l_matrix);
            double endTime = MPI_Wtime();
            cout << "Serial Calculation time: " << endTime - startTime << endl;
        }
        double startTime = MPI_Wtime();
        FindMatLParallel(l_matrix_copy);
        double endTime = MPI_Wtime();
        if (rank == 0) cout << "Parallel Calculation time: " << endTime - startTime << endl;

    // �������� ������� b (����������� ������ �������)
    vector <float> b_vector(a_size);
    if (rank == 0) {
        for (int j = 0; j < a_size; j++) {
            b_vector[j] = rand() % 10; // ��������� ������ ����������
        }
        cout << "Vector [b]:\n";
        ShowVector(b_vector);
    }

    // ���������� ������� b � ����: [LY = B] (����������� ������ �������)
    vector<vector<float>> lb_matrix(a_size, vector <float>(a_size+1));
    if (rank == 0){
        for (int i = 0; i < a_size; i++) {
            for (int j = 0; j < a_size + 1; j++) {
                if (j == a_size) {
                    lb_matrix[i][j] = b_vector[i];
                    continue;
                }
                lb_matrix[i][j] = l_matrix[i][j];
            }
        }
        cout << "[LY = B] matrix\n";
        ShowMatrixEx(lb_matrix);
    }

    double delta = 1e-08;
    if (rank == 0) {
        SolveMatrix(lb_matrix, delta);
        cout << " SOLVED [LY = B] matrix\n";
        ShowMatrixEx(lb_matrix);
    }

    // ���������������� ������� L (����������� ������ �������)
    vector<vector<float>> lt_matrix(a_size, vector <float>(a_size));
    if (rank == 0) {
        for (int i = 0; i < a_size; i++) {
            for (int j = 0; j < a_size; j++) {
                lt_matrix[i][j] = l_matrix[j][i];
            }
        }
        vector<vector<float>> lx_matrix(a_size, vector <float>(a_size + 1));
        for (int i = 0; i < a_size; i++) {
            for (int j = 0; j < a_size + 1; j++) {
                if (j == a_size) {
                    lx_matrix[i][j] = lb_matrix[i][j];
                    continue;
                }
                lx_matrix[i][j] = lt_matrix[i][j];
            }
        }
        cout << "[LtX = Y] matrix\n";
        ShowMatrixEx(lx_matrix);
        SolveMatrix(lx_matrix, delta); // ������� ���� ������� ������-�������
        cout << "SOLVED [LtX = Y] matrix\n";
        ShowMatrixEx(lx_matrix); // ����� �������� ��������� �������
        ShowResult(lx_matrix); // ����� ������� ���� Ax = b
    }

    MPI_Finalize(); // ���������� MPI ���������
}




int main(int argc, char** argv)
{
    const int a_size = 3;
    doKholesky(argc, argv, a_size);
	return 0;
}

/*
#include <iostream>
#include <fstream>
#include <string>
#include <mpi.h>

using namespace std;

void mainProcWork(int, int);
void secondaryProcWork(int, int);

// ������������ ��������� ������ A � B 
int main(int argc, char* argv[])
{
    int rank, numProc;
    MPI_Init(&argc, &argv); // MPI-�������������
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // ���� �������� ��������
    MPI_Comm_size(MPI_COMM_WORLD, &numProc); // ���-�� ���������

    string arg = argv[1];
    int size, block;
    size = stoi(arg); // ������ ������
    block = size / numProc; // ������� ������� A �� ����� �����

    if (rank == 0) mainProcWork(size, block);
    else secondaryProcWork(size, block);

    MPI_Finalize(); // ���������� MPI ���������

    return 0;
}

void mainProcWork(int size, int block)
{
    int numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc); // ���-�� ���������

    int* a, * b, * c;
    a = (int*)malloc(sizeof(int) * size * size);
    b = (int*)malloc(sizeof(int) * size * size);
    c = (int*)malloc(sizeof(int) * size * size);

    // ������ ������� A �� �����
    ifstream inA("A.txt");
    for (auto i = 0; i < size; i++)
        for (auto j = 0; j < size; j++)
            inA >> a[i * size j];
    inA.close();

    // ������ ������� B �� �����
    ifstream inB("B.txt");
    for (auto i = 0; i < size; i++)
        for (auto j = 0; j < size; j++)
            inB >> b[i * size + j];
    inB.close();

    double startTime, endTime;
    startTime = MPI_Wtime();

    int* recBuf = (int*)malloc(sizeof(int) * size * block); // ����� ��� �������� ����������� ���������� ������ ���������

    // ����������� ������� B �������������� ���������
    for (auto i = 1; i < numProc; i++)
    {
        MPI_Send(b, size * size, MPI_INT, i, 0, MPI_COMM_WORLD);
    }

    // ����������� ������ ����� ������� A �������������� ���������
    for (auto j = 1; j < numProc; j++)
    {
        MPI_Send(a + (j - 1) * block * size, size * block, MPI_INT, j, 1, MPI_COMM_WORLD);
    }

    // ���������� ���������� �������� ��������� (������� � �����)
    for (auto i = (numProc - 1) * block; i < size; i++)
        for (auto j = 0; j < size; j++)
        {
            int tmp = 0;
            for (auto k = 0; k < size; k++)
                tmp += a[i * size + k] * b[k * size + j];
            c[i * size + j] = tmp;
        }

    // ��������� ����������� ���������� �� ������ ���������
    for (auto k = 1; k < numProc; k++)
    {
        MPI_Recv(recBuf, block * size, MPI_INT, k, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // ����������� �����������
        for (auto i = 0; i < block; i++)
        {
            for (auto j = 0; j < size; j++)
            {
                c[((k - 1) * block + i) * size + j] = recBuf[i * size + j];
            }

        }
    }

    // ������ ������� C � ����
    ofstream outC("C.txt");
    for (auto i = 0; i < size; i++)
    {
        for (auto j = 0; j < size; j++)
            outC << c[i * size + j] << " ";
        outC << endl;
    }
    outC.close();

    // ����� ���������� ������������ ���������
    endTime = MPI_Wtime();
    cout << "Calculation time: " << endTime - startTime << endl;

    free(a);
    free(b);
    free(c);
    free(recBuf);
}

void secondaryProcWork(int size, int block)
{
    int rank, numProc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // ���� �������� ��������
    MPI_Comm_size(MPI_COMM_WORLD, &numProc); // ���-�� ���������

    int* partA, * b;
    partA = (int*)malloc(sizeof(int) * block * size);
    b = (int*)malloc(sizeof(int) * size * size);

    // ��������� ������� B �� ��������� ��������
    MPI_Recv(b, size * size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // ��������� ����� ������� A �� ��������� ��������
    MPI_Recv(partA, size * block, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    int* sendBuf = (int*)malloc(sizeof(int) * block * size); // ����� ��� ����������� ����������� ���������� �������������� ���������

    // ���������� ���������� �������������� ���������
    for (auto i = 0; i < block; i++)
        for (auto j = 0; j < size; j++)
        {
            int tmp = 0;
            for (auto k = 0; k < size; k++)
                tmp += partA[i * size + k] * b[k * size + j];
            sendBuf[i * size + j] = tmp;
        }

    // ����������� ����������� ���������� ��������� ��������
    MPI_Send(sendBuf, block * size, MPI_INT, 0, 3, MPI_COMM_WORLD);
}

*/

/*
#include<mpi.h>
#include<iostream>
#include <cstdlib>
#include <vector>

using namespace std;

void ShowMatrix(vector<vector<float> >& a, int mode) {
    if (mode == 1) {
        for (int i = 0; i < a.size(); i++) {
            for (int j = 0; j < a[i].size(); j++) {
                if (a[i].size() - j == 2) {
                    cout << a[i][j] << " | "; // ����������� ��������
                    continue;                      // � ����������� ������
                }
                cout << a[i][j] << " ";
            }
            cout << "\n";
        }
    }
    else {
        cout << "Linear equations system answer is:" << endl;
        for (int i = 0; i < a.size(); i++) {
            cout << "X" << i << " = " << a[i][a.size()] << "\n";
        }
    }
    cout << endl;
}

int main(int argc, char* argv[]) {
    int n = 3; // ������� ������� �������
    vector<vector<float> > a(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n + 1; j++) {
            a[i].push_back(rand() % 10); // ��������� ������ ����������
        }                                 // a[i][j] = {0, 99}
    }
    cout << "hello";
    ShowMatrix(a, 1);


}
*/