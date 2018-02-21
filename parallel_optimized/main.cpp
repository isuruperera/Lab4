#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
#include <cmath>

int N=2000;
double min = 0.0;
double max = 1000.0;

void multiplyMatrices(double** matrix_one, double** matrix_two, double** result_matrix);
void runSimulation(double time_matrix[],int samples);
void initMatrix(double ** matrix);
double randomValue();
void printMatrix(double **matrix);
double calculateSD(double data[],int samples);
bool isEnoughSamples(int sampleCount, double stdDev, double mean);
double calculateMean(double data[],int samples);

int main()
{
    std::cout << "================******=================="<<std::endl;
    std::cout << "CS4532: Concurrent Programming "<<std::endl;
    std::cout << "Lab 3 "<<std::endl;
    std::cout << "P.D.I.T.S.K. Perera \t 140462E "<<std::endl;
    std::cout << "M.D.S.K. Ranasinghe \t 140499X "<<std::endl;
    std::cout << "================******=================="<<std::endl;

    for(int k = 200;k<=2000;k+=200)
    {
        N = k;
        int sampleCount = 10;
        bool sampleSizeFound = false;

        //Find the required sample size
        while (!sampleSizeFound)
        {
            double timeValues[sampleCount];
            runSimulation(timeValues,sampleCount);
            double stdDev = calculateSD(timeValues,sampleCount);
            double mean = calculateMean(timeValues,sampleCount);
            sampleSizeFound = isEnoughSamples(sampleCount,stdDev,mean);
            if(!sampleSizeFound)
            {
                sampleCount += 5;
            }
        }

        double finalTimeValues[sampleCount];
        runSimulation(finalTimeValues,sampleCount);
        double stdDev = calculateSD(finalTimeValues,sampleCount);
        double mean = calculateMean(finalTimeValues,sampleCount);

        //Print results
        std::cout << "======================================"<<std::endl;
        std::cout << "N \t\t\t\t: "<< N <<std::endl;
        std::cout << "Sample Size\t\t\t: "<< sampleCount <<std::endl;
        std::cout << "Std Deviation\t\t\t: "<< stdDev <<std::endl;
        std::cout << "Mean\t\t\t\t: "<< mean <<std::endl;
        std::cout << "Time Taken for executions: "<<std::endl;
        for(int i=0;i<sampleCount;i++)
        {
            std::cout << finalTimeValues[i]<<" ";
        }
        std::cout <<std::endl<< "======================================"<<std::endl;

    }
    return 0;
}

//Multiply given matrices
void multiplyMatrices(double** matrix_one, double** matrix_two, double** result_matrix)
{

    // Multiplying matrix firstMatrix and secondMatrix and storing in array result_matrix.
    //Runs as a parallel for loop using OpenMP library
    #pragma omp parallel for
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            //OPTIMIZATION
            //Sum is first calculated into a local variable and then populated into
            //the result matrix for reducing memory access operations and increase
            //cache hits
            double sum = 0.0;
            for(int k=0; k<N; k++)
            {
                sum += matrix_one[i][k] * matrix_two[k][j];
            }
            result_matrix[i][j] = sum;
        }
    }
}

//Run a single simulation with given sample size
void runSimulation(double time_matrix[],int samples)
{
    using namespace std::chrono;

    //Random seed for each simulation
    std::srand(std::time(nullptr));
    double** input_1 = new double*[N];
    double** input_2 = new double*[N];
    double** input_3 = new double*[N];
    double** output = new double*[N];

    for(int i=0;i<N;i++)
    {
        input_1[i] = new double[N];
        input_2[i] = new double[N];
        input_3[i] = new double[N];
        output[i] = new double[N];
    }



    //Multiply matrices while recording time
    for(int i=0;i<samples;i++)
    {
        initMatrix(input_1);
        initMatrix(input_2);

        //OPTIMIZATION
        //When the second matrix is accessed column wise, it would cause cache misses.
        //So, the second matrix is transposed in order to use spatial locality and
        //reduce cache misses
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                input_3[i][j] = input_2[j][i];
            }
        }

        high_resolution_clock::time_point t1 = high_resolution_clock::now();

        multiplyMatrices(input_1,input_3,output);

        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double, std::milli> time_span = t2 - t1;
        time_matrix[i] = time_span.count();

    }
}

//Generate a random value
double randomValue()
{
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

//Populate matrix with random values
void initMatrix(double **matrix)
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            matrix[i][j] = randomValue();
        }
    }
}

//Check if enough samples are generated
bool isEnoughSamples(int sampleCount, double stdDev, double mean)
{
    double n = (100 * 1.96 * stdDev)/(5*mean);
    n = n*n;
    std::cout << "Trying "<< sampleCount <<" samples, n= "<< n <<std::endl;
    if(n < sampleCount)
    {
        std::cout << "Success at: "<< sampleCount <<" samples, n= "<< n <<std::endl;
        return true;
    } else
    {
        return false;
    }
}



//Print a given matrix
void printMatrix(double **matrix)
{
    std::cout<<std::endl;
    std::cout << "Array print"<<std::endl;
    std::cout<<std::endl;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            std::cout << matrix[i][j]<<"   ";
        }
        std::cout<<std::endl;
    }
}

//Calculate standard deviation of given data set
double calculateSD(double data[],int samples)
{
    double mean, standardDeviation = 0.0;
    mean = calculateMean(data,samples);
    for(int i = 0; i < samples; ++i)
        standardDeviation += pow(data[i] - mean, 2);
    return sqrt(standardDeviation / samples);
}

//Calculate mean of given data set
double calculateMean(double data[],int samples)
{
    double mean = 0.0;
    double sum = 0.0;
    for(int i = 0; i < samples; ++i)
    {
        sum += data[i];
    }
    mean = sum/samples;
    return mean;
}


