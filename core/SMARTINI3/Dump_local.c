#include <stdio.h>
#include <stdlib.h>

typedef double matrix[3][3];
typedef float  fmatrix[3][3];

int main(int argc, char **argv)
{
    if (argc < 4) {
        printf("Usage: %s input_file output_file slice_size\n", argv[0]);
        return 1;
    }

    // Variables
    int i, j, k, d, dd;
    int nx, ny, nz;
    int isDoublePrecision;
    float sliceSize;
    
    float fBox[3][3];
    double box[3][3];
    
    matrix *localPressure; // 3x3 pressure tensor per grid point
    matrix *slice;         // 3x3 slice along z-axis
    fmatrix tmpMat;
    matrix average;

    FILE *fp = fopen(argv[1], "r");
    if (!fp) {
        perror("Failed to open input file");
        return 1;
    }

    sscanf(argv[3], "%f", &sliceSize);
    FILE *out = fopen(argv[2], "w");
    if (!out) {
        perror("Failed to open output file");
        fclose(fp);
        return 1;
    }

    // Read precision flag
    fread(&isDoublePrecision, sizeof(int), 1, fp);

    // Read box dimensions
    if (isDoublePrecision) {
        fread(box, sizeof(double), 9, fp);
    } else {
        fread(fBox, sizeof(float), 9, fp);
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                box[i][j] = fBox[i][j];
    }

    // Read grid dimensions
    fread(&nx, sizeof(int), 1, fp);
    fread(&ny, sizeof(int), 1, fp);
    fread(&nz, sizeof(int), 1, fp);

    // Allocate memory
    localPressure = malloc(sizeof(matrix) * nx * ny * nz);
    slice = malloc(sizeof(matrix) * nz);

    // Read local pressure data
    if (isDoublePrecision) {
        fread(localPressure, sizeof(matrix), nx * ny * nz, fp);
    } else {
        for (k = 0; k < nx * ny * nz; k++) {
            fread(tmpMat, sizeof(tmpMat), 1, fp);
            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                    localPressure[k][i][j] = tmpMat[i][j];
        }
    }

    fclose(fp);

    // Initialize average matrix
    for (d = 0; d < 3; d++)
        for (dd = 0; dd < 3; dd++)
            average[d][dd] = 0.0;

    // Process each slice along z-axis
    for (k = 0; k < nz; k++) {
        // Initialize slice
        for (d = 0; d < 3; d++)
            for (dd = 0; dd < 3; dd++)
                slice[k][d][dd] = 0.0;

        // Sum over x and y
        for (i = 0; i < nx; i++)
            for (j = 0; j < ny; j++)
                for (d = 0; d < 3; d++)
                    for (dd = 0; dd < 3; dd++)
                        slice[k][d][dd] += localPressure[i * ny * nz + j * nz + k][d][dd] / (nx * ny);

        // Update average
        for (d = 0; d < 3; d++)
            for (dd = 0; dd < 3; dd++)
                average[d][dd] += slice[k][d][dd] / nz;

        // Print slice to console
        printf("P[z%2d]: %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
               k,
               slice[k][0][0], slice[k][0][1], slice[k][0][2],
               slice[k][1][0], slice[k][1][1], slice[k][1][2],
               slice[k][2][0], slice[k][2][1], slice[k][2][2]);

        // Write to output file
        fprintf(out, "%12.6f %12.6f %12.6f %12.6f %12.6f\n",
                k * sliceSize,
                slice[k][0][0], slice[k][1][1], slice[k][2][2],
                0.5f * (slice[k][0][0] + slice[k][1][1]) - slice[k][2][2]);
    }

    // Print average pressure
    printf("Average: %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
           average[0][0], average[0][1], average[0][2],
           average[1][0], average[1][1], average[1][2],
           average[2][0], average[2][1], average[2][2]);

    fclose(out);

    // Free memory
    free(localPressure);
    free(slice);

    return 0;
}
