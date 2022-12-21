#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <vector>
using namespace std;

    //Функции для печати матриц и векторов в консоль
    void ShowMatrix(vector<vector<float> >& a) {
    // Вывод матрицы
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
                    cout << a[i][j] << " | "; // разделитель основной
                    continue;                      // и расширенной матриц
                }
                cout << a[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "\n";
    }
    void ShowVector(vector<float> & a){
        //Вывод вектора
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

    // Функции для решения методом Гаусса-Жордана
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
            if (abs(a[k][k]) < delta) {                  // Проверка.
                for (int i = k + 1; i < a.size(); i++) { // Если элемент
                    if (a[i][k] > delta) {               // равен нулю, то
                        swap(a[i], a[k]);           // поменять местами
                        break;                           // строки с ненулевой
                    }
                }
            }
            // Если элемент не равен единице - делим строку на него
            if (a[k][k] < (1 - delta) || a[k][k] > (1 + delta)) {
                RowDivision(a[k], a[k][k]);
            }
            // прямой ход обнуления элементов матрицы
            // отнять от оставшихся ниже строк опорную, умноженную на элемент

            for (int i = k + 1; i < a.size(); i++) {
                RowElimination(a[i], a[k], a[i][k]);
            }
        }
        // обратный ход обнуления элементов матрицы

        for (int k = a.size() - 1; k >= 0; k--) {
            for (int i = k - 1; i >= 0; i--) { // Обратный ход реализуется здесь
                RowElimination(a[i], a[k], a[i][k]);
            }
        }
    }

    //Функция для генерации и преобразования матрицы A
    void GenMatA(vector<vector <float>> &a_matrix){
        // Создание матрицы А
        vector<vector<float>> a_matrix_v1(a_matrix.size());
        // Заполнение матрицы случайными элементами (симметрично)
        for (int i = 0; i < a_matrix.size(); i++) {
            for (int j = 0; j < a_matrix.size(); j++) {
                a_matrix_v1[i].push_back(rand() % 10); // заполняем строки элементами
            }                                 // a[i][j] = {0, 9}
        }
        ShowMatrix(a_matrix_v1);
        // Транспонирование матрицы
        vector<vector<float>> a_matrix_v2(a_matrix.size(), vector <float>(a_matrix.size()));
        for (int i = 0; i < a_matrix.size(); i++) {
            for (int j = 0; j < a_matrix.size(); j++) {
                a_matrix_v2[i][j] = a_matrix_v1[j][i];
            }
        }
        ShowMatrix(a_matrix_v2);
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
        cout << "Matrix A:\n";
        ShowMatrix(a_matrix);
    }

    //Функция для нахождения матрицы L методом Холецкого
    void FindMatL(vector<vector<float>> &l_matrix) {
        int numProc, rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Ранг текущего процесса
        MPI_Comm_size(MPI_COMM_WORLD, &numProc); // Кол-во процессов
        if (rank == 0) {
            // Вычисление матрицы L методом Хоолецкого
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
            cout << "SOLVED L_Matrix\n";
            ShowMatrix(l_matrix); // получили нижнюю треугольную матрицу
        }
    }


// Основная функция программы.
void doKholesky(int argc, char** argv, int a_size) {
    int numProc, rank;
    MPI_Init(&argc, &argv); // MPI-инициализация
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Ранг текущего процесса
    MPI_Comm_size(MPI_COMM_WORLD, &numProc); // Кол-во процессов
    
    //Создание симметричной положительно определённой матрицы А
    vector<vector<float>> a_matrix(a_size, vector <float>(a_size));
    if (rank == 0) GenMatA(a_matrix);

    // Копирование в матрицу L нижней треугольной матрицы из А
    // Основные вычисления для разложения Холецкого - поиск матрицы L
    vector<vector<float>> l_matrix(a_size, vector <float>(a_size));
        for (int i = 0; i < a_size; i++) {
            for (int j = 0; j <= i; j++) {
                l_matrix[i][j] = a_matrix[i][j];
            }
        }
        if (rank == 0) ShowMatrix(l_matrix);
        FindMatL(l_matrix);

    // Создание вектора b (Выполняется главны потоком)
    vector <float> b_vector(a_size);
    if (rank == 0) {
        for (int j = 0; j < a_size; j++) {
            b_vector[j] = rand() % 10; // заполняем строки элементами
        }
        cout << "Vector [b]:\n";
        ShowVector(b_vector);
    }

    // Добавление вектора b к СЛАУ: [LY = B] (Выполняется главны потоком)
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

    // транспонирование матрицы L (Выполняется главны потоком)
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
        SolveMatrix(lx_matrix, delta); // Решение СЛАУ методом Гаусса-Жордана
        cout << "SOLVED [LtX = Y] matrix\n";
        ShowMatrixEx(lx_matrix); // Вывод конечной единичной матрицы
        ShowResult(lx_matrix); // Вывод решения СЛАУ Ax = b
    }

    MPI_Finalize(); // Завершение MPI программы
}




int main(int argc, char** argv)
{
    int a_size = 3;
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

// Параллельное умножение матриц A и B 
int main(int argc, char* argv[])
{
    int rank, numProc;
    MPI_Init(&argc, &argv); // MPI-инициализация
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Ранг текущего процесса
    MPI_Comm_size(MPI_COMM_WORLD, &numProc); // Кол-во процессов

    string arg = argv[1];
    int size, block;
    size = stoi(arg); // Размер матриц
    block = size / numProc; // Деление матрицы A на блоки строк

    if (rank == 0) mainProcWork(size, block);
    else secondaryProcWork(size, block);

    MPI_Finalize(); // Завершение MPI программы

    return 0;
}

void mainProcWork(int size, int block)
{
    int numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc); // Кол-во процессов

    int* a, * b, * c;
    a = (int*)malloc(sizeof(int) * size * size);
    b = (int*)malloc(sizeof(int) * size * size);
    c = (int*)malloc(sizeof(int) * size * size);

    // Чтение матрицы A из файла
    ifstream inA("A.txt");
    for (auto i = 0; i < size; i++)
        for (auto j = 0; j < size; j++)
            inA >> a[i * size j];
    inA.close();

    // Чтение матрицы B из файла
    ifstream inB("B.txt");
    for (auto i = 0; i < size; i++)
        for (auto j = 0; j < size; j++)
            inB >> b[i * size + j];
    inB.close();

    double startTime, endTime;
    startTime = MPI_Wtime();

    int* recBuf = (int*)malloc(sizeof(int) * size * block); // Буфер для принятия результатов вычислений других процессов

    // Отправление матрицы B второстепенным процессам
    for (auto i = 1; i < numProc; i++)
    {
        MPI_Send(b, size * size, MPI_INT, i, 0, MPI_COMM_WORLD);
    }

    // Отправление блоков строк матрицы A второстепенным процессам
    for (auto j = 1; j < numProc; j++)
    {
        MPI_Send(a + (j - 1) * block * size, size * block, MPI_INT, j, 1, MPI_COMM_WORLD);
    }

    // Проведение вычислений основным процессом (остаток с конца)
    for (auto i = (numProc - 1) * block; i < size; i++)
        for (auto j = 0; j < size; j++)
        {
            int tmp = 0;
            for (auto k = 0; k < size; k++)
                tmp += a[i * size + k] * b[k * size + j];
            c[i * size + j] = tmp;
        }

    // Получение результатов вычислений от других процессов
    for (auto k = 1; k < numProc; k++)
    {
        MPI_Recv(recBuf, block * size, MPI_INT, k, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Объединение результатов
        for (auto i = 0; i < block; i++)
        {
            for (auto j = 0; j < size; j++)
            {
                c[((k - 1) * block + i) * size + j] = recBuf[i * size + j];
            }

        }
    }

    // Запись матрицы C в файл
    ofstream outC("C.txt");
    for (auto i = 0; i < size; i++)
    {
        for (auto j = 0; j < size; j++)
            outC << c[i * size + j] << " ";
        outC << endl;
    }
    outC.close();

    // Время вычислений параллельной программы
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
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Ранг текущего процесса
    MPI_Comm_size(MPI_COMM_WORLD, &numProc); // Кол-во процессов

    int* partA, * b;
    partA = (int*)malloc(sizeof(int) * block * size);
    b = (int*)malloc(sizeof(int) * size * size);

    // Получение матрицы B от основного процесса
    MPI_Recv(b, size * size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Получение блока матрицы A от основного процесса
    MPI_Recv(partA, size * block, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    int* sendBuf = (int*)malloc(sizeof(int) * block * size); // Буфер для отправления результатов вычислений второстепенным процессом

    // Проведение вычислений второстепенным процессом
    for (auto i = 0; i < block; i++)
        for (auto j = 0; j < size; j++)
        {
            int tmp = 0;
            for (auto k = 0; k < size; k++)
                tmp += partA[i * size + k] * b[k * size + j];
            sendBuf[i * size + j] = tmp;
        }

    // Отправление результатов вычислений основному процессу
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
                    cout << a[i][j] << " | "; // разделитель основной
                    continue;                      // и расширенной матриц
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
    int n = 3; // Задание размера матрицы
    vector<vector<float> > a(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n + 1; j++) {
            a[i].push_back(rand() % 10); // заполняем строки элементами
        }                                 // a[i][j] = {0, 99}
    }
    cout << "hello";
    ShowMatrix(a, 1);


}
*/