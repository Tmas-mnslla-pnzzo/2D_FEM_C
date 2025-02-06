#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fem_funciones.h"

#define MAX_CLAS_EXT 1
#define MAX_CLAS_INT 1

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Error: No se ingresó el nombre del archivo.\n");
        return 1;
    }
    const char *ARCHIVO_BASE = argv[1];
    printf("Iniciando el programa...\n");

    double cond_dir_ext[MAX_CLAS_EXT] = {0};
    double cond_dir_int[MAX_CLAS_INT] = {0};
    double puntos[MAX_LIST_LENGTH][2]; 
    int triangulos[MAX_LIST_LENGTH][3];  
    char clasificaciones[MAX_LIST_LENGTH][MAX_CELL_LENGTH]; 
    int n_p = 0; 
    int n_t = 0;  
    char nombre_archivo[MAX_CELL_LENGTH];
    char nombre_archivo_csv[MAX_CELL_LENGTH];

    snprintf(nombre_archivo_csv, MAX_CELL_LENGTH, "resultado_U_de_%s.csv", ARCHIVO_BASE);
    snprintf(nombre_archivo, MAX_CELL_LENGTH, "%s_puntos.csv", ARCHIVO_BASE);

    printf("Leyendo archivo de puntos...\n");
    FILE *file_puntos = fopen(nombre_archivo, "r");
    if (file_puntos == NULL) {
        perror("Error al abrir el archivo G_puntos.csv");
        return 1;
    }

    char line[MAX_LINE_LENGTH];
    while (fgets(line, sizeof(line), file_puntos)) {
        line[strcspn(line, "\n")] = 0;
        char *cell = strtok(line, ",");
        int n_celda = 0;

        while (cell != NULL && n_celda < 2) {
            puntos[n_p][n_celda] = atof(cell);
            cell = strtok(NULL, ",");
            n_celda++;
        }
        if (n_celda == 2) n_p++;
    }
    fclose(file_puntos);
    printf("Archivo de puntos leído. Total de puntos: %d\n", n_p);

    printf("Leyendo archivo de triángulos...\n");
    snprintf(nombre_archivo, MAX_CELL_LENGTH, "%s_elementos.csv", ARCHIVO_BASE);
    FILE *file_triangulos = fopen(nombre_archivo, "r");
    if (file_triangulos == NULL) {
        perror("Error al abrir el archivo G_triangulos.csv");
        return 1;
    }

    while (fgets(line, sizeof(line), file_triangulos)) {
        line[strcspn(line, "\n")] = 0;
        char *cell = strtok(line, ",");
        int n_celda = 0;

        while (cell != NULL && n_celda < 3) {
            triangulos[n_t][n_celda] = atoi(cell);
            cell = strtok(NULL, ",");
            n_celda++;
        }
        if (n_celda == 3) n_t++;
    }
    fclose(file_triangulos);
    printf("Archivo de triángulos leído. Total de triángulos: %d\n", n_t);

    printf("Leyendo archivo de clasificaciones...\n");
    snprintf(nombre_archivo, MAX_CELL_LENGTH, "%s_clas.csv", ARCHIVO_BASE);
    FILE *file_clasificaciones = fopen(nombre_archivo, "r");
    if (file_clasificaciones == NULL) {
        perror("Error al abrir el archivo G_clasificaciones.csv");
        return 1;
    }

    while (fgets(line, sizeof(line), file_clasificaciones)) {
        line[strcspn(line, "\n")] = 0;
        if (n_p < MAX_LIST_LENGTH) {
            strncpy(clasificaciones[n_p], line, MAX_CELL_LENGTH - 1);
            clasificaciones[n_p][MAX_CELL_LENGTH - 1] = '\0';
        }
    }
    fclose(file_clasificaciones);
    printf("Archivo de clasificaciones leído.\n");

    printf("Inicializando matrices globales y vector F...\n");
    double **G = malloc(n_p * sizeof(double *));
    double **M = malloc(n_p * sizeof(double *));
    double **C = malloc(n_p * sizeof(double *));
    double **B = malloc(n_p * sizeof(double *));
    double **A = malloc(n_p * sizeof(double *));
    double *F = malloc(n_p * sizeof(double));

    if (!G || !M || !C || !B || !A || !F) {
        perror("Error en la asignación de memoria");
        return 1;
    }

    for (int i = 0; i < n_p; i++) {
        G[i] = calloc(n_p, sizeof(double));
        M[i] = calloc(n_p, sizeof(double));
        C[i] = calloc(n_p, sizeof(double));
        B[i] = calloc(n_p, sizeof(double));
        A[i] = calloc(n_p, sizeof(double));
        if (!G[i] || !M[i] || !C[i] || !B[i] || !A[i]) {
            perror("Error en la asignación de memoria para las matrices");
            return 1;
        }
    }

    memset(F, 0, n_p * sizeof(double));

    double m_local[3][3] = {{2, 1, 1}, {1, 2, 1}, {1, 1, 2}};
    double a_local[3][3] = {{2, -1, -1}, {-1, 1, 0}, {-1, 0, 1}};

    double (*f_p[3])(double, double, double, double, double, double, double, double) = {f_p1, f_p2, f_p3};
    double (*c_p[3][3])(double, double, double, double, double, double, double, double) = {{c_p11, c_p12, c_p13},{c_p21, c_p22, c_p23},{c_p31, c_p32, c_p33}};
    double (*b_p[3][3])(double, double, double, double, double, double, double, double) = {{b_p11, b_p12, b_p13},{b_p21, b_p22, b_p23},{b_p31, b_p32, b_p33}};

    printf("Construyendo matrices globales...\n");
    for (int elem = 0; elem < n_t; elem++) {
        printf("Procesando triángulo %d de %d...\n", elem + 1, n_t);

        int nod1 = triangulos[elem][0];
        int nod2 = triangulos[elem][1];
        int nod3 = triangulos[elem][2];

        double x1 = puntos[nod1][0]; 
        double y1 = puntos[nod1][1];
        double x2 = puntos[nod2][0]; 
        double y2 = puntos[nod2][1];
        double x3 = puntos[nod3][0];
        double y3 = puntos[nod3][1];

        double determinante = det_jacob(x1, y1, x2, y2, x3, y3);
        double factor = determinante / 24;

        double fi[3], mi[3][3], ci[3][3], bi[3][3], ai[3][3];

        for (int j = 0; j < 3; j++) {
            FuncParams fp_f = {x1, y1, x2, y2, x3, y3, f_p[j]};
            fi[j] = determinante * gauss(x1, y1, x2, y2, x3, y3, f_p[j]);

            for (int k = 0; k < 3; k++) {
                FuncParams fp_c = {x1, y1, x2, y2, x3, y3, c_p[j][k]};
                FuncParams fp_b = {x1, y1, x2, y2, x3, y3, b_p[j][k]};
                FuncParams fp_a = {x1, y1, x2, y2, x3, y3, a};  

                ci[j][k] = determinante * gauss(x1, y1, x2, y2, x3, y3, c_p[j][k]);
                bi[j][k] = determinante * gauss(x1, y1, x2, y2, x3, y3, b_p[j][k]);
                ai[j][k] = determinante * gauss(x1, y1, x2, y2, x3, y3, a) * a_local[j][k];
                mi[j][k] = factor * m_local[j][k];
            }
        }

        for (int j = 0; j < 3; j++) {
            F[triangulos[elem][j]] += fi[j];
            for (int k = 0; k < 3; k++) {
                M[triangulos[elem][j]][triangulos[elem][k]] += mi[j][k];
                C[triangulos[elem][j]][triangulos[elem][k]] += ci[j][k];
                B[triangulos[elem][j]][triangulos[elem][k]] += bi[j][k];
                A[triangulos[elem][j]][triangulos[elem][k]] += ai[j][k];
            }
        }
    }
    
    printf("Matrices globales construidas.\n");

    for (int j = 0; j < n_p; j++) {
        for (int k = 0; k < n_p; k++) {
            G[j][k] = A[j][k] + B[j][k] + C[j][k];
        }
    }
    
    printf("Aplicando condiciones de frontera a las matrices...\n");
    aplicar_dir_M(G, clasificaciones, n_p);

    printf("Aplicando condiciones de frontera a los vectores...\n");
    aplicar_dir_ext_V(cond_dir_ext, MAX_CLAS_EXT, puntos, F, clasificaciones, n_p);
    aplicar_dir_int_V(cond_dir_int, MAX_CLAS_INT, puntos, F, clasificaciones, n_p);
    printf("Condiciones de frontera aplicadas.\n");
    
    double **invG = malloc(n_p * sizeof(double *));
    double *U = malloc(n_p * sizeof(double));

    if (!invG || !U) {
        perror("Error en la asignación de memoria para invG y U");
        return 1;
    }

    for (int i = 0; i < n_p; i++) {
        invG[i] = malloc(n_p * sizeof(double));
        if (!invG[i]) {
            perror("Error en la asignación de memoria para invG[i]");
            return 1;
        }
    }
    
    memset(U, 0, n_p * sizeof(double));

    invertMatrix(n_p,G,invG);

    for (int i = 0; i < n_p; i++) {
        free(G[i]);
        free(M[i]);
        free(C[i]);
        free(B[i]);
        free(A[i]);
    }
    free(G);
    free(M);
    free(C);
    free(B);
    free(A);
    free(F);

    mulMatVec(n_p, n_p, invG, F, U);
    guardarVectorCSV(nombre_archivo_csv, U, n_p);

    for (int i = 0; i < n_p; i++) {
        free(invG[i]);
    }
    free(invG);
    free(U);
    
    printf("Programa finalizado con éxito.\n");
    return 0;
}
