#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fem_funciones.h" 

void printMat(double **mat, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {  // Antes: j <= N, ahora j < N
            printf("%.6f\t", mat[i][j]);
        }
        printf("\n");
    }
}

void printVec(const double *vec, int N, int n) {
    if (n > N) {
        n = N;  
    }

    int start = N - n;

    for (int i = start; i < N; i++) {
        printf("%.6f\t", vec[i]);
    }
    printf("\n");
}
void guardarVectorCSV(const char *nombreArchivo, const double *vector, int n) {
    FILE *archivo = fopen(nombreArchivo, "w"); 
    if (archivo == NULL) {
        perror("Error al abrir el archivo");
        return;
    }
    for (int i = 0; i < n; i++) {
        fprintf(archivo, "%.6f\n", vector[i]);  
    }
    fflush(archivo);  // Asegura que los datos se escriban antes de cerrar
    fclose(archivo);
    printf("Vector guardado en %s\n", nombreArchivo);
}

void guardar_vectores_bin(const char *nombreArchivo, double *U, int n_p, int tiempo) {
    FILE *archivo = fopen(nombreArchivo, "ab");  // "ab" → Append en binario
    if (archivo == NULL) {
        perror("Error al abrir el archivo binario");
        return;
    }
    fwrite(U, sizeof(double), n_p, archivo);  // Guarda un vector U en cada instante
    fclose(archivo);
}

void borrar_contenido_archivo(const char *nombreArchivo) {
    FILE *archivo = fopen(nombreArchivo, "wb");  // "wb" → Write en binario (borra contenido)
    if (archivo == NULL) {
        perror("Error al abrir el archivo");
        return;
    }
    fclose(archivo);  // Cierra el archivo inmediatamente (el contenido ya se borró al abrirlo en modo "wb")
}

double f_wrapper(double u, double v, void *params) {
    FuncParams *p = (FuncParams *)params;
    return p->func(u, v, p->x1, p->y1, p->x2, p->y2, p->x3, p->y3);
}

double f_p1(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return f(u, v, x1, y1, x2, y2, x3, y3) * p1(u, v); 
}
double f_p2(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return f(u, v, x1, y1, x2, y2, x3, y3) * p2(u); 
}
double f_p3(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return f(u, v, x1, y1, x2, y2, x3, y3) * p3(v); 
}
double c_p11(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return c(u, v, x1, y1, x2, y2, x3, y3) * p1(u, v) * p1(u, v); 
}
double c_p12(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return c(u, v, x1, y1, x2, y2, x3, y3) * p1(u, v) * p2(u); 
}
double c_p13(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return c(u, v, x1, y1, x2, y2, x3, y3) * p1(u, v) * p3(v); 
}
double c_p21(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return c(u, v, x1, y1, x2, y2, x3, y3) * p2(u) * p1(u, v); 
}
double c_p22(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return c(u, v, x1, y1, x2, y2, x3, y3) * p2(u) * p2(u);
}
double c_p23(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return c(u, v, x1, y1, x2, y2, x3, y3) * p2(u) * p3(v); 
}
double c_p31(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return c(u, v, x1, y1, x2, y2, x3, y3) * p3(v) * p1(u, v); 
}
double c_p32(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return c(u, v, x1, y1, x2, y2, x3, y3) * p3(v) * p2(u); 
}
double c_p33(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return c(u, v, x1, y1, x2, y2, x3, y3) * p3(v) * p3(v); 
}
double b_p11(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return -p1(u,v) * (b1(u, v, x1, y1, x2, y2, x3, y3)+b2(u, v, x1, y1, x2, y2, x3, y3)); 
}
double b_p12(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return p1(u, v) * b1(u, v, x1, y1, x2, y2, x3, y3); 
}
double b_p13(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return p1(u, v) * b2(u, v, x1, y1, x2, y2, x3, y3); 
}
double b_p21(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return -p2(u) * (b1(u, v, x1, y1, x2, y2, x3, y3)+b2(u, v, x1, y1, x2, y2, x3, y3)); 
}
double b_p22(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return p2(u) * b1(u, v, x1, y1, x2, y2, x3, y3); 
}
double b_p23(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return p2(u) * b2(u, v, x1, y1, x2, y2, x3, y3); 
}
double b_p31(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return -p3(v) * (b1(u, v, x1, y1, x2, y2, x3, y3)+b2(u, v, x1, y1, x2, y2, x3, y3)); 
}
double b_p32(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return p3(v) * b1(u, v, x1, y1, x2, y2, x3, y3); 
}
double b_p33(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) { 
    return p3(v) * b2(u, v, x1, y1, x2, y2, x3, y3); 
}

double f_0(double x, double y){
    return 0;
}
double c_0(double x, double y){
    return 0;
}
double b1_0(double x, double y){
    return 0;
}
double b2_0(double x, double y){
    return 0;
}
double a_0(double x, double y){
    return 1;
}

double f(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) {
    double coord[2]; 
    T(x1, y1, x2, y2, x3, y3, u, v, coord);
    return f_0(coord[0], coord[1]);  
}
double c(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) {
    double coord[2]; 
    T(x1, y1, x2, y2, x3, y3, u, v, coord);
    return c_0(coord[0], coord[1]);  
}

double b1(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) {
    double coord[2]; 
    T(x1, y1, x2, y2, x3, y3, u, v, coord);
    return b1_0(coord[0], coord[1]);  
}

double b2(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) {
    double coord[2]; 
    T(x1, y1, x2, y2, x3, y3, u, v, coord);
    return b2_0(coord[0], coord[1]);  
}

double a(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3) {
    double coord[2]; 
    T(x1, y1, x2, y2, x3, y3, u, v, coord);
    return a_0(coord[0], coord[1]);  
}

double p1(double u, double v) {
    return 1 - u - v;
}

double p2(double u) {
    return u;
}

double p3(double v) {
    return v;
}
double det_jacob(double x1, double y1, double x2, double y2, double x3, double y3){
    double D;
    D=(x1-x3)*(y3-y1)-(x3-x1)*(y1-y2);
    return fabs(D);
}
void T(double x1, double y1, double x2, double y2, double x3, double y3, double u, double v, double coord[2]){
    double x;
    double y;
    x=(x1-x2)*u+(x3-x1)*v+x1;
    y=(y1-y2)*u+(y3-y1)*v+y1;
    coord[0]=x; 
    coord[1]=y;
}

void productoPunto(int n, double vector1[n], double vector2[n], double *resultado) {
    *resultado = 0.0; 
    for (int j = 0; j < n; j++) {
        *resultado += vector1[j] * vector2[j];  
    }
}

void mulMatVec(int m, int n, double **matriz, double *vector, double *resultado) {
    for (int i = 0; i < m; i++) {
        resultado[i] = 0.0; 
        for (int j = 0; j < n; j++) {
            resultado[i] += matriz[i][j] * vector[j];  
        }
    }
}

int invertMatrix(int n, double **mat, double **inv) {
    double **tempMat = malloc(n * sizeof(double *));
    if (!tempMat) {
        perror("Error de asignación de memoria");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n; i++) {
        tempMat[i] = malloc(n * sizeof(double));
        if (!tempMat[i]) {
            perror("Error de asignación de memoria");
            exit(EXIT_FAILURE);
        }
        memcpy(tempMat[i], mat[i], n * sizeof(double));
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int i = 0; i < n; i++) {
        if (tempMat[i][i] == 0.0) {
            int swapRow = -1;
            for (int k = i + 1; k < n; k++) {
                if (tempMat[k][i] != 0.0) {
                    swapRow = k;
                    break;
                }
            }
            if (swapRow == -1) {
                printf("La matriz no es invertible (determinante cero).\n");
                exit(EXIT_FAILURE);
            }
            for (int j = 0; j < n; j++) {
                double temp = tempMat[i][j];
                tempMat[i][j] = tempMat[swapRow][j];
                tempMat[swapRow][j] = temp;

                temp = inv[i][j];
                inv[i][j] = inv[swapRow][j];
                inv[swapRow][j] = temp;
            }
        }
        double pivot = tempMat[i][i];
        for (int j = 0; j < n; j++) {
            tempMat[i][j] /= pivot;
            inv[i][j] /= pivot;
        }
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = tempMat[k][i];
                for (int j = 0; j < n; j++) {
                    tempMat[k][j] -= factor * tempMat[i][j];
                    inv[k][j] -= factor * inv[i][j];
                }
            }
        }
    }

    for (int i = 0; i < n; i++) {
        free(tempMat[i]);
    }
    free(tempMat);

    return 1;
}

double gauss(double x1, double y1, double x2, double y2, double x3, double y3, double (*func)(double, double, double, double, double, double, double, double)) {
    int g_p=5;
    double gauss_points[7] = {
    -0.9491079123427585,
    -0.7415311855993945,
    -0.4058451513773972,
    0.0,
    0.4058451513773972,
    0.7415311855993945,
    0.9491079123427585
    };

    double gauss_weights[7] = {
    0.1294849661688697,
    0.2797053914892766,
    0.3818300505051189,
    0.4179591836734694,
    0.3818300505051189,
    0.2797053914892766,
    0.1294849661688697
    };

    double integral = 0.0;
    for (int i = 0; i < g_p; i++) {
        for (int j = 0; j < g_p; j++) {
            double xi = gauss_points[i];
            double eta = gauss_points[j];

            double u = (1.0 + xi) / 2.0;
            double v = (1.0 + eta) * (1.0 - u) / 2.0;

            double jacobian = (1.0 - xi) / 8.0;

            double weight = gauss_weights[i] * gauss_weights[j] * jacobian;

            integral += weight * func(u, v, x1, y1, x2, y2, x3, y3);
        }
    }
    return integral;
}

void aplicar_init(double *Tbl, int nTbl, double P[][2], double *U, int clas[][1], int nP, const char *prefix, int y) {
    if (nTbl == 0) {
        return;  
    }
    for (int i = 0; i < nP; i++) {  
        for (int j = 0; j < nTbl; j++) { 
            if (clas[i][0] == y*(j+1) && Tbl[j] != -1) {
                U[i] = Tbl[j];  
            }
        }
    }
}

void aplicar_dir_V(double *F, int clas[][1], int n, double *Tbl, int nTbl, const char *prefix, int y) {
    if (nTbl == 0) {
        return;  
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < nTbl; j++) {
            if (clas[i][0] == y*(j+1) && !isnan(Tbl[j])) {
                F[i] = Tbl[j];  
            }
        }
    }
}

void aplicar_dir_M(double **G, int clas[][1], int n, double *Tbl, int nTbl, const char *prefix, int y) {
    if (nTbl == 0) {
        return;  
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < nTbl; j++) {
            if (clas[i][0] == y*(j+1) && !isnan(Tbl[j])) {
                for (int k = 0; k < n; k++) {
                    G[i][k] = 0.0;
                }
                G[i][i] = 1.0;
            }
        }
    }
}
