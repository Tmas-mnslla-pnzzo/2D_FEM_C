#ifndef FEMFUNCIONES_H
#define FEMFUNCIONES_H
#define MAX_LINE_LENGTH 1024
#define MAX_CELL_LENGTH 400
#define MAX_LIST_LENGTH 4000

typedef struct {
    double x1, y1, x2, y2, x3, y3;
    double (*func)(double, double, double, double, double, double, double, double);
} FuncParams;

void printMat(double **mat, int N);
void printVec(const double *vec, int N, int n);
void guardarVectorCSV(const char *nombreArchivo, const double *vector, int n);

void guardar_vectores_bin(const char *nombreArchivo, double *U, int n_p, int tiempo);
void borrar_contenido_archivo(const char *nombreArchivo);

double f_wrapper(double u, double v, void *params);

double f_p1(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double f_p2(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double f_p3(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);

double c_p11(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double c_p12(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double c_p13(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double c_p21(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double c_p22(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double c_p23(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double c_p31(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double c_p32(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double c_p33(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);

double b_p11(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double b_p12(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double b_p13(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double b_p21(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double b_p22(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double b_p23(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double b_p31(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double b_p32(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double b_p33(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);

double f_0(double x, double y);
double c_0(double x, double y);
double b1_0(double x, double y);
double b2_0(double x, double y);
double a_0(double x, double y);

double f(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double c(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double b1(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double b2(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);
double a(double u, double v, double x1, double y1, double x2, double y2, double x3, double y3);

double p1(double u, double v);
double p2(double u);
double p3(double v);

double det_jacob(double x1, double y1, double x2, double y2, double x3, double y3);
void T(double x1, double y1, double x2, double y2, double x3, double y3, double r, double s, double coord[2]);

void productoPunto(int n, double vector1[n], double vector2[n], double *resultado);
void mulMatVec(int m, int n, double **matriz, double *vector, double *resultado);
int invertMatrix(int n, double **mat, double **inv);

double gauss(double x1, double y1, double x2, double y2, double x3, double y3, double (*func)(double, double, double, double, double, double, double, double));

void aplicar_init(double *Tbl, int nTbl, double P[][2], double *U, int clas[][1], int nP, const char *prefix, int y);
void aplicar_dir_V(double *F, int clas[][1], int n, double *Tbl, int nTbl, const char *prefix, int y);
void aplicar_dir_M(double **G, int clas[][1], int n, double *Tbl, int nTbl, const char *prefix, int y);

#endif
