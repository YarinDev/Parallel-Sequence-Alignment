// Yarin Getter
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>

#define ROOT 0
#define SEQ_1_MAX_SIZE 5001
#define SEQ_MAX_SIZE 3001
#define ABC_SIZE 26
#define CONSERVATIVE_SIZE 9
#define SEMI_CONSERVATIVE_SIZE 11

char **alloc_2d_char(int rows, int cols);
int **alloc_2d_int(int rows, int cols);
void strUpper(char *str);

void setScoringMatrix(int scoringMatrix[ABC_SIZE][ABC_SIZE], int *weights);
int calcSinglePair(int group[ABC_SIZE][ABC_SIZE], char c1, char c2);
int calcAlignmentScore(int scoringMatrix[ABC_SIZE][ABC_SIZE], char *seq_1, int seq_1_size, char *mutant, int mutant_size, int n);
void calcWithOpenMP(int scoringMatrix[ABC_SIZE][ABC_SIZE], char **sequence_work_arr, int **results_work, int work_size, char *seq_1, int size_of_seq_1, int my_rank);
void sequentialProgramming(int scoringMatrix[ABC_SIZE][ABC_SIZE], char *seq_1, int number_of_sequences, char **seq_arr);
void printSequences(char **seq_2_arr, int **results_arr, int size);

const char *conservative_str[ABC_SIZE] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
const char *semiConservative_str[ABC_SIZE] = {"SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"};

/*
allocates 2 dimensions array of int given by rows and columns
*/
int **alloc_2d_int(int rows, int cols)
{
    int *data = (int *)malloc(rows * cols * sizeof(int));
    int **array = (int **)malloc(rows * sizeof(int *));
    for (int i = 0; i < rows; i++)
        array[i] = &(data[cols * i]);

    return array;
}
/*
allocates 2 dimensions array of chars given by rows and columns
*/
char **alloc_2d_char(int rows, int cols)
{
    char *data = (char *)malloc(rows * cols * sizeof(char));
    char **array = (char **)malloc(rows * sizeof(char *));
    for (int i = 0; i < rows; i++)
        array[i] = &(data[cols * i]);

    return array;
}

/*
make sure that all chars are capital letters.
*/
void strUpper(char *str)
{
    int i = 0;
    while (str[i])
    {
        str[i] = toupper(str[i]);
        i++;
    }
}

/*
Calculate pairs via weights and the defaults Strings
1. by default cells == weights[3]
2. weights[2]
3. weights[1]
4. The cells are equal
the order of thr loops is importent because it creates the best possible score for each pair of letters.
*/
void setScoringMatrix(int scoringMatrix[ABC_SIZE][ABC_SIZE], int *weights)
{
    int i, x, y, size;
    // 1.
    for (i = 0; i < ABC_SIZE; i++)
        for (x = 0; x < ABC_SIZE; x++)
            scoringMatrix[i][x] = -weights[3];

    // 2.
    for (i = 0; i < SEMI_CONSERVATIVE_SIZE; i++)
    {
        size = strlen(semiConservative_str[i]);
        for (x = 0; x < size - 1; x++)
        {
            for (y = x + 1; y < size; y++)
            {
                scoringMatrix[semiConservative_str[i][x] - 'A'][semiConservative_str[i][y] - 'A'] = -weights[2];
                scoringMatrix[semiConservative_str[i][y] - 'A'][semiConservative_str[i][x] - 'A'] = -weights[2];
            }
        }
    }
    // 3.
    for (i = 0; i < CONSERVATIVE_SIZE; i++)
    {
        size = strlen(conservative_str[i]);
        for (x = 0; x < size - 1; x++)
        {
            for (y = x + 1; y < size; y++)
            {
                scoringMatrix[conservative_str[i][x] - 'A'][conservative_str[i][y] - 'A'] = -weights[1];
                scoringMatrix[conservative_str[i][y] - 'A'][conservative_str[i][x] - 'A'] = -weights[1];
            }
        }
    }

    // 4.
    for (i = 0; i < ABC_SIZE; i++)
        scoringMatrix[i][i] = weights[0];
}

/*
Calculate a single pair-> minus the A letter to get the value
*/
int calcSinglePair(int group[ABC_SIZE][ABC_SIZE], char c1, char c2)
{
    return group[c1 - 'A'][c2 - 'A'];
}

/*
Calculate the alingment score for all pairs of str2 mutant, given by offset, "n" & "k" combination.
*/
int calcAlignmentScore(int scoringMatrix[ABC_SIZE][ABC_SIZE], char *seq_1, int seq_1_size, char *mutant, int mutant_size, int n)
{
    int i, j = 0;
    int score = 0;
    if (n + mutant_size > seq_1_size)
        printf("Error: The mutant size is bigger than the Sequence #1 size");

    for (i = 0; i < mutant_size; i++, j++)
    {
        score += calcSinglePair(scoringMatrix, mutant[i], seq_1[j + n]);
    }
    return score;
}
/*
Using openmp to calculate str2 offset,n,k combination
*/
void calcWithOpenMP(int scoringMatrix[ABC_SIZE][ABC_SIZE], char **sequence_work_arr, int **results_work, int work_size, char *seq_1, int size_of_seq_1, int my_rank)
{

    int i, j, seq_2_size, n, k, offset, score;

    for (i = 0; i < work_size; i++)
    {

        results_work[i][0] = INT_MIN;

        seq_2_size = strlen(sequence_work_arr[i]);

#pragma omp parallel for private(x, y, temp_n, score)

        for (n = 0; n < seq_2_size - 1; n++)
        {
            for (k = n + 1; k <= seq_2_size - 1; k++)
            {

                char *mutant = (char *)malloc(1000000 * sizeof(char));
                int idxToDel = k;

                // copy str2 completly to mutant
                strcpy(mutant, sequence_work_arr[i]);

                // copy mutant to himself but without the first chosen index
                memmove(&mutant[n], &mutant[n + 1], strlen(mutant) - n);

                // copy mutant to himself but without the second chosen index
                memmove(&mutant[k - 1], &mutant[k], strlen(mutant));

                // offsets possible itartions for each mutant
                for (offset = 0; offset < (size_of_seq_1 - (seq_2_size - 2)); offset++)
                {
                    score = calcAlignmentScore(scoringMatrix, seq_1, size_of_seq_1, mutant, seq_2_size - 2, offset);

#pragma omp critical

                    if (score > results_work[i][0])
                    {
                        results_work[i][3] = k;
                        results_work[i][2] = n;
                        results_work[i][1] = offset;
                        results_work[i][0] = score;
                    }
                }
                free(mutant);
            }
        }
    }
}

/*
Calculate the score in the sequantial way (not efficient comparing to parallel)
*/
void sequentialProgramming(int scoringMatrix[ABC_SIZE][ABC_SIZE], char *seq_1, int number_of_sequences, char **seq_arr)
{
    int i, offset, n, k, temp_offset, temp_k, temp_n, size_of_seq_1, max_score, score, seq_2_size;
    size_of_seq_1 = strlen(seq_1);
    for (int i = 0; i < number_of_sequences; i++)
    {
        max_score = INT_MIN;
        offset = 0;
        n = 0;
        k = 0;
        int seq_2_size = strlen(seq_arr[i]);
        for (int temp_n = 0; temp_n < seq_2_size - 1; temp_n++)
        {
            for (int temp_k = temp_n + 1; temp_k <= seq_2_size - 1; temp_k++)
            {

                char *mutant = (char *)malloc(1000000 * sizeof(char));
                int idxToDel = temp_k;
                // copy str2 completly to mutant
                strcpy(mutant, seq_arr[i]);

                // copy mutant to himself but without the first chosen index
                memmove(&mutant[temp_n], &mutant[temp_n + 1], strlen(mutant) - temp_n);

                // copy mutant to himself but without the second chosen index
                memmove(&mutant[temp_k - 1], &mutant[temp_k], strlen(mutant));

                // offsets possible itartions for each mutant
                for (temp_offset = 0; temp_offset < (size_of_seq_1 - (seq_2_size - 2)); temp_offset++)
                {
                    score = calcAlignmentScore(scoringMatrix, seq_1, size_of_seq_1, mutant, seq_2_size - 2, temp_offset);

                    if (score > max_score)
                    {
                        k = temp_k;
                        n = temp_n;
                        offset = temp_offset;
                        max_score = score;
                    }
                }
                free(mutant);
            }
        }
        printf("Sequence # %d:\noffset: %d, MS(%d, %d)\n", (i + 1), offset, n, k);
        printf("The Best Alignment Score is: %d\n\n", max_score);
    }
}

void printSequences(char **seq_2_arr, int **results_arr, int size)
{
    int i;
    for (i = 0; i < size; i++)
    {
        printf("Sequence # %d:\noffset = %d, MS(%d, %d)\n", (i + 1), results_arr[i][1], results_arr[i][2], results_arr[i][3]);
        printf("The Best Alignment Score is: %d\n\n", results_arr[i][0]);
    }
}

int main(int argc, char *argv[])
{
    int my_rank, weights[4], number_of_sequences, size_of_seq_1, i, num_of_procs, remainings = 0, work_size, **results, **results_work;
    int scoringMatrix[ABC_SIZE][ABC_SIZE];
    char **seq_arr, **sequence_work_arr, **sequenceLeft_work_arr;
    char seq_1[SEQ_1_MAX_SIZE];
    double t;

    // MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_of_procs);

    if (my_rank == ROOT)
    {
        // Reading file input text
        fscanf(stdin, "%d %d %d %d", &weights[0], &weights[1], &weights[2], &weights[3]);
        fscanf(stdin, "%s", seq_1);
        strUpper(seq_1);
        fscanf(stdin, "%d", &number_of_sequences);
        size_of_seq_1 = strlen(seq_1);
        seq_arr = alloc_2d_char(number_of_sequences, SEQ_MAX_SIZE);
        results = alloc_2d_int(number_of_sequences, 4);

        for (i = 0; i < number_of_sequences; i++)
        {
            fscanf(stdin, "%s", seq_arr[i]);
            if (strlen(seq_arr[i]) > size_of_seq_1)
            {
                printf("ERROR: Sequence #1 must be bigger than Sequence #2\n");
                exit(1);
            }
            strUpper(seq_arr[i]);
        }
        // Calculating the scoringMatrix using weights and strings.
        setScoringMatrix(scoringMatrix, weights);
        printf("\n\nThe Scoring matrice: (i have print it just for better explantion while presentation)\n");
        for (int i = 0; i < ABC_SIZE; i++)
        {
            for (int j = 0; j < ABC_SIZE; j++)
            {
                if (j == 25)
                {
                    printf("%d \n", scoringMatrix[i][j]);
                }
                else
                {
                    printf("%d ", scoringMatrix[i][j]);
                }
            }
        }
        printf("\n\n******************** All sequences with their offset and MS(n,k) that producing the best Alignment scores ******************** \n");
        printf("Parallel Programming: \n");
        t = MPI_Wtime();
        remainings = number_of_sequences % num_of_procs;
    }
    /*
    Broadcasting the data to the proccess via MPI commands
    */
    MPI_Bcast(&remainings, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&size_of_seq_1, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(weights, 4, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&number_of_sequences, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(seq_1, size_of_seq_1 + 1, MPI_CHAR, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(scoringMatrix, ABC_SIZE * ABC_SIZE, MPI_INT, ROOT, MPI_COMM_WORLD);
    work_size = number_of_sequences / num_of_procs;

    // Allocating matrices for results
    sequence_work_arr = alloc_2d_char(work_size + remainings, SEQ_MAX_SIZE);
    results_work = alloc_2d_int(work_size + remainings, 4);

    if (my_rank == ROOT)
    {

        for (i = 1; i < num_of_procs; i++)
            MPI_Send(&(seq_arr[0][0]) + ((remainings + (i * work_size)) * SEQ_MAX_SIZE), work_size * SEQ_MAX_SIZE, MPI_CHAR, i, 0, MPI_COMM_WORLD);

        // taking the root results and put them in the final array
        results_work = results;
        sequence_work_arr = seq_arr;

        // openMP calculation function
        calcWithOpenMP(scoringMatrix, sequence_work_arr, results_work, work_size + remainings, seq_1, size_of_seq_1, my_rank);

        // Recieving the Results
        for (i = 1; i < num_of_procs; i++)
        {
            MPI_Recv(&(results[0][0]) + ((remainings + (i * work_size)) * 4), work_size * 4, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        printSequences(seq_arr, results, number_of_sequences);
        printf("Parallel runtime: %lf\n", MPI_Wtime() - t);
        printf("Sequential Programming: \n");
        t = MPI_Wtime();
        sequentialProgramming(scoringMatrix, seq_1, number_of_sequences, seq_arr);
        printf("Sequential runtime: %lf\n", MPI_Wtime() - t);
    }
    if (my_rank != ROOT)
    {

        // Recieving the sequences from MPI
        MPI_Recv(&(sequence_work_arr[0][0]), work_size * SEQ_MAX_SIZE, MPI_CHAR, ROOT, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Calculate scores
        calcWithOpenMP(scoringMatrix, sequence_work_arr, results_work, work_size, seq_1, size_of_seq_1, my_rank);
        MPI_Send(&(results_work[0][0]), work_size * 4, MPI_INT, ROOT, 0, MPI_COMM_WORLD);
    }

    // free memory in the results and seqs arrays
    free(results_work[0]);
    free(sequence_work_arr[0]);
    free(results_work);
    free(sequence_work_arr);

    MPI_Finalize();
    return 0;
}