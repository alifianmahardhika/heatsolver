//Function
//Function to calculate determinant
double area(double x11, double x12, double x21, double x22, double x31, double x32){
  //double det = (((x21*x32) - (x31*x22)) - (x11*(x32-x22)) + (x12*(x31-x21)));
  double a = x21-x11;
  double b = x22-x12;
  double c = x31-x11;
  double d = x32-x12;
  return (a*d-b*c)/2.0;
}
//fx in RHS
double bfunction(double ax, double bx, double nu, double t){
  double fx = (1.0+(nu*2.0*M_PI*M_PI*t))*sin(M_PI*ax)*sin(M_PI*bx);
  return fx;
}
//exact function
double exactfunction(double ay, double by, double t){
  double yx = t*sin(M_PI*ay)*sin(M_PI*by);
  return yx;
}
//print dmatrix
void printDoubleMatrix(int nrow, int ncol, double mat[nrow][ncol]) {
  for (int i = 1; i <= nrow; i++) {
    for (int j = 1; j <= ncol; j++) {
      printf("%.2f ", mat[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}
//print dvector
void printDoubleVector(int ncol, double mat[ncol]) {
  for (int i = 1; i <= ncol; i++) {
    printf("%f\n", mat[i]);
  }
  printf("\n");
}
//find maximum value
double maxValueVector(int ncol, double mat[ncol]) {
  double max = 0.0;
  for (int i = 1; i <= ncol; i++) {
    if (mat[i] > max) {
      max = mat[i];
    }
  }
  return max;
}
//find Sum of vector
double sumVector(int ncol, double mat[ncol]) {
  double sum = 0.0;
  for (int i = 1; i <= ncol; i++) {
    sum += mat[i];
  }
  return sum;
}
