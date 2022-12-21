#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <vector>
using namespace std;

const int a_size = 3;
    //������� ��� ������ ������ � �������
    void ShowMatrix(float a[a_size][a_size]) {
    // ����� �������
        for (int i = 0; i < a_size; i++) {
            for (int j = 0; j < a_size; j++) {
                cout << a[i][j] << " ";
            }
            cout << "\n";
        }
    }
    void ShowMatrixEx(float a[a_size][a_size + 1]) {
        for (int i = 0; i < a_size; i++) {
            for (int j = 0; j < a_size + 1; j++) {
                if (a_size - j == 1) {
                    cout << a[i][j] << " | "; // ����������� ��������
                    continue;                      // � ����������� ������
                }
                cout << a[i][j] << " ";
            }
            cout << "\n";
        }
    }
    void ShowVector(float a[a_size]){
        //����� ������� ������� ���������
        cout << "Linear equations system answer is:" << endl;
        for (int i = 0; i < a_size; i++) {
            cout << "X" << i << " = " << a[i] << "\n";
        }
    }

void MasterThread() {
    int numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc); // ���-�� ���������

    // �������� ������� �
    float a_matrix_v1[a_size][a_size];
    // ���������� ������� ���������� ���������� (�����������)
    for (int i = 0; i < a_size; i++) {
        for (int j = 0; j < a_size; j++) {
            a_matrix_v1[i][j] = rand() % 10; // ��������� ������ ����������
        }                                 // a[i][j] = {0, 9}
    }
    // ���������������� �������
    float a_matrix_v2[a_size][a_size];
    for (int i = 0; i < a_size; i++) {
        for (int j = 0; j < a_size; j++) {
            a_matrix_v2[i][j] = a_matrix_v1[j][i];
        }
    }


    float a_matrix[a_size][a_size], result = 0;
    //Fill each cell of the matrix output.
    for (int i = 0; i < a_size; i++) {
        for (int j = 0; j < a_size; j++) {
            //Multiply each row of matrix 1 with each column of matrix 2.
            for (int k = 0; k < a_size; k++) {
                result += a_matrix_v1[i][k] * a_matrix_v2[k][j];
            }
            a_matrix[i][j] = result;
            result = 0; //Reset;
        }
    }
    cout << "Matrix in master\n";
    ShowMatrix(a_matrix);
 

    // ���������� �������� ������� ���������� ������ (���������� ���������)
            float l_matrix[a_size][a_size];
            for (int i = 0; i < a_size; i++) {
                for (int j = 0; j < a_size; j++) {
                    l_matrix[i][j] = 0;
                }
            }
            for (int i = 0; i < a_size; i++) {
                for (int j = 0; j <= i; j++) {
                    l_matrix[i][j] = a_matrix[i][j];
                }
            }
            for (int j = 0; j < a_size; j++) {
                for (int k = 0; k < j; k++) {
                    for (int i = j; i < a_size; i++) {
                        l_matrix[i][j] = l_matrix[i][j] - l_matrix[i][k] * l_matrix[j][k];
                    }
                }
                l_matrix[j][j] = sqrt(l_matrix[j][j]);
                for (int i = j + 1; i < a_size; i++) {
                    l_matrix[i][j] = l_matrix[i][j] / l_matrix[j][j];
                }
            }
            cout << "L_matrix\n";
            ShowMatrix(l_matrix); // �������� ������ ����������� �������

    // ���������� ������� b � ���������� �������
    float b_vector[a_size];
    for (int j = 0; j < a_size; j++) {
        b_vector[j] = rand() % 10; // ��������� ������ ����������
    }
    ShowVector(b_vector);

    float lb_matrix[a_size][a_size + 1];
    for (int i = 0; i < a_size; i++) {
        for (int j = 0; j < a_size + 1; j++) {
            if (j == a_size) {
                lb_matrix[i][j] = b_vector[i];
                continue;
            }
            lb_matrix[i][j] = l_matrix[i][j];
        }
    }
    cout << "LB_matrix\n";
    ShowMatrixEx(lb_matrix);



    // ����������� ������� � �������������� ���������
    for (auto i = 1; i < numProc; i++)
    {
        MPI_Send(a_matrix, sizeof(a_matrix), MPI_FLOAT, i, 0, MPI_COMM_WORLD);
    }
    
}

void SlaveThread() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // ���� �������� ��������
    
    float recBufA[a_size][a_size];
    MPI_Recv(recBufA, sizeof(recBufA), MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //cout << "Matrix in slave numba" << rank << endl;
    //ShowMatrix(recBufA);
    

}

// ������� ������� �������

int main(int argc, char** argv)
{
    int rank, numProc;
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // ���� �������� ��������
    MPI_Comm_size(MPI_COMM_WORLD, &numProc); // ���-�� ���������

    if (rank == 0) MasterThread();// ������ �������� master ������
    else SlaveThread();// ������ ������� slave �������


	MPI_Finalize();

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