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

void printVec(const double *vec, int N) {
    for (int i = 0; i < N; i++) {
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
    return x*y*(1-x)*(1-y);
}
double c_0(double x, double y){
    return x*y;
}
double b1_0(double x, double y){
    return x;
}
double b2_0(double x, double y){
    return y;
}
double a_0(double x, double y){
    return x+y;
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

    // Inicializar la identidad en inv
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Algoritmo de Gauss-Jordan para invertir la matriz
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

    // Liberar memoria de la copia
    for (int i = 0; i < n; i++) {
        free(tempMat[i]);
    }
    free(tempMat);

    return 1;
}


double trapecio(double (*func)(double, double, void *), double a, double b, double y, int n, void *params) {
    double h = (b - a) / n;  
    double sum = 0.5 * (func(a, y, params) + func(b, y, params));  

    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += func(x, y, params);
    }

    return sum * h;
}

double gauss(double x1, double y1, double x2, double y2, double x3, double y3, double (*func)(double, double, double, double, double, double, double, double)) {
    double gauss_points[5] = {-0.90618, -0.53847, 0.0, 0.53847, 0.90618};
    double gauss_weights[5] = {0.23693, 0.47863, 0.56888, 0.47863, 0.23693};

    double integral = 0.0;
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            double u = gauss_points[i];
            double v = gauss_points[j];
            double weight = gauss_weights[i] * gauss_weights[j];
            integral += weight * func(u, v, x1, y1, x2, y2, x3, y3);
        }
    }
    return integral;
}

double integrar(double (*func)(double, double, void *), void *params) {
    double a = 0.0, b = 1.0;
    int n = 1000;
    double h = (b - a) / n;
    double sum = 0.0;

    for (int i = 0; i < n; i++) {
        double s = a + i * h;
        double integral_y = 0.0;

        for (int j = 0; j < n; j++) {
            double r = a + j * h;
            integral_y += func(r, s, params);
        }
        sum += integral_y * h;  
    }
    return sum * h; 
}

void aplicar_dir_ext_V(double *Tb, int nTb, double P[][2], double *U, char cl[][MAX_CELL_LENGTH], int nP) {
    for (int i = 0; i < nP; i++) { 
        for (int j = 0; j < nTb; j++) { 
            char buffer[10];
            snprintf(buffer, sizeof(buffer), "ext%d", j); 
            if (strcmp(cl[i], buffer) == 0 && Tb[j] != -1) {  
                U[i] = Tb[j];  
            }
        }
    }
}

void aplicar_dir_int_V(double *Tbl, int nTbl, double P[][2], double *U, char cl[][MAX_CELL_LENGTH], int nP) {
    if (nTbl == 0) return;  

    for (int i = 0; i < nP; i++) {  
        for (int j = 0; j < nTbl; j++) { 
            char buffer[10];
            snprintf(buffer, sizeof(buffer), "int%d", j);  
            if (strcmp(cl[i], buffer) == 0 && Tbl[j] != -1) {
                U[i] = Tbl[j];  
            }
        }
    }
}

void aplicar_dir_M(double **G, char cl[][MAX_CELL_LENGTH], int n) {
    for (int j = 0; j < n; j++) {
        if (strcmp(cl[j], "None") != 0) {
            for (int k = 0; k < n; k++) {
                G[j][k] = 0;
            }
            G[j][j] = 1;
        }
    }
}
