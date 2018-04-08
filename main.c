//Initialize Library needed for calculation
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"./inc/functions.c"

//Function Prototype
double area();
double bfunction();
double exactfunction();
double maxValueVector();
double sumVector();

//Void Prototype
void printDoubleMatrix();
void printDoubleVector();

//Main Program
int main(int argc, char *argv[]){
  //initialize matrix
  double **npxy,**localA,**globalA,**cmatrix,**xmatrix;
  double **m;

  //initialize vector
  double *areamatrix,*localB,*globalB;
  double *v,*uh,*ue,*mv,*vhilbert;
  double *residu,*pk,*adotx,*adotpk,*residuold;

  //initialize variable
  double nablaphi = 0.0, meas = 0.0;
  double residunorm,residunormnew,residunormold;
  double alpha,alpha1,alpha2,beta;
  double epsilon,deltanew,deltaold,deltanol,errcheck;
  double x11,x12,x21,x22,x31,x32;
  int DIM=2, **elnp, **bdnp, label, i, j, k, l, np, ne, nb, iter, itermax;

  //addition for solving heat equation
  //initialize matrix
  double **locm,**globalM,**heatLHS;
  //initialize vector
  double *locn,*globalN,*heatRHS;
  //initialize variable
  int t;
  FILE *plotsolution;
  char fname[100],pngname[100],pictloc[100];
  FILE *fp;
  //plot log error
  FILE *plotlogerror; //begin plot log error
  plotlogerror = fopen("plotlogerror.txt","w");
  if (plotlogerror==NULL) {
    printf("Error Creating File plotlogerror.txt\n");
    exit(0);
  }
  //open or read mesh file
  fp=fopen(argv[1],"r");
  //error handling
  if (fp==NULL) {
    printf("Cannot open file: %s.\n", argv[1]);
    exit(1);
  }
  //read the first line of file
  fscanf(fp,"%d %d %d", &np, &ne, &nb);

  // allocation for matrix
  npxy = dmatrix(1,np,1,DIM);
  elnp = imatrix(1,ne,1,DIM+1);
  bdnp = imatrix(1,nb,1,DIM);
  localA = dmatrix(1,3,1,3);
  localB = dvector(1,3);
  globalA = dmatrix(1,np,1,np);
  globalB = dvector(1,np);
  areamatrix = dvector(1,ne);
  xmatrix = dmatrix(1,2,1,3);
  cmatrix = dmatrix(1,2,1,3);
  residu = dvector(1,np);
  pk = dvector(1,np);
  adotpk = dvector(1,np);
  residuold = dvector(1,np);
  adotx = dvector(1,np);
  v = dvector(1,np);
  m = dmatrix(1,3,1,3);
  uh = dvector(1,np);
  ue = dvector(1,np);
  mv = dvector(1,3);
  vhilbert = dvector(1,2);
  //addition for heat equation
  locm = dmatrix(1,3,1,3);
  locn = dvector(1,3);
  globalM = dmatrix(1,np,1,np);
  globalN = dvector(1,np);
  heatLHS = dmatrix(1,np,1,np);
  heatRHS = dvector(1,np);

  //print out number of nodal point, element of nodal point, and border line of nodal point
  printf("Number of Nodal Point, Element, and Border\n np = %d\n ne = %d\n nb = %d\n", np, ne, nb);
  //read the value of nodal point coordinate
  for(i=1; i<=np; i++){
    fscanf(fp, "%lf %lf %d", &npxy[i][1], &npxy[i][2], &label);
  }
  //read the value of elemental nodal Point
  for(i=1; i<=ne; i++){
    fscanf(fp, "%d %d %d %d", &elnp[i][1], &elnp[i][2], &elnp[i][3], &label);
  }
  //read the value of boundary nodal point
  for(i=1; i<=nb; i++){
    fscanf(fp, "%d %d %d", &bdnp[i][1], &bdnp[i][2], &label);
  }
  //Close reading mesh
  fclose(fp);
  //Calculating Area of the Element
for (i=1; i<=ne; i++) {
  x11 = npxy[elnp[i][1]][1];
  x12 = npxy[elnp[i][1]][2];
  x21 = npxy[elnp[i][2]][1];
  x22 = npxy[elnp[i][2]][2];
  x31 = npxy[elnp[i][3]][1];
  x32 = npxy[elnp[i][3]][2];
  areamatrix[i] = area(x11,x12,x21,x22,x31,x32);
}
double xdt;
double dt;
dt = 0.0625;
//delta time
double tmax, nu;
int nt;
tmax = 1;
nt = tmax/dt;
printf("NT = %d\n", nt);
nu = 1;
//do { //Begin of do while
  //addition for heat problem
  //given u0
  for (i = 1; i <= np; i++) {
    uh[i] = 0.0;
  }

  //calculate local and global LHS
  for (k = 1; k <= ne; k++) { //Begin Loop Each Element
    //calculating c1 and c2
    for (i = 1; i <= 2; i++) {
      for (j = 1; j <= 3; j++) {
        xmatrix[i][j] = npxy[elnp[k][j]][i];
      }
    }
    meas = areamatrix[k];
    cmatrix[1][1] = (xmatrix[2][2]-xmatrix[2][3])/(2.0*meas);
    cmatrix[1][2] = (xmatrix[2][3]-xmatrix[2][1])/(2.0*meas);
    cmatrix[1][3] = (xmatrix[2][1]-xmatrix[2][2])/(2.0*meas);
    cmatrix[2][1] = (xmatrix[1][3]-xmatrix[1][2])/(2.0*meas);
    cmatrix[2][2] = (xmatrix[1][1]-xmatrix[1][3])/(2.0*meas);
    cmatrix[2][3] = (xmatrix[1][2]-xmatrix[1][1])/(2.0*meas);
    //calculating nablaphi, local LHS
    for (i = 1; i <= 3; i++) {
      for (j = 1; j <= 3; j++) {
        nablaphi = (cmatrix[1][i]*cmatrix[1][j])+(cmatrix[2][i]*cmatrix[2][j]);
        localA[i][j] = areamatrix[k]*nablaphi*nu;
        globalA[elnp[k][i]][elnp[k][j]] += localA[i][j];
        if (i==j) {
          locm[i][j] = areamatrix[k]/(6.0*dt);
        } else {
          locm[i][j] = areamatrix[k]/(12.0*dt);
        }
        globalM[elnp[k][i]][elnp[k][j]] += locm[i][j];
      }
    }
  }//End Loop Each Element
  //modifying LHS for heat equation
  for (i = 1; i <= np; i++) {
    for (j = 1; j <= np; j++) {
      heatLHS[i][j] = globalA[i][j] + globalM[i][j];
    }
  }
  //imposing boundary
  for (i = 1; i <= nb; i++) {
    for (j = 1; j <= np; j++) {
      heatLHS[bdnp[i][1]][j] = 0.0;
      heatLHS[j][bdnp[i][1]] = 0.0;
    }
    heatLHS[bdnp[i][1]][bdnp[i][1]] = 1.0;
  } //End Loop Each Element

  //beginning solve heat equation
  double time;
  // //begin to plot png
  // char plottxtname[100];
  // FILE *pngplotter;
  // sprintf(plottxtname, "plotpng.gpl");
  // pngplotter = fopen(plottxtname,"w");
  // fprintf(pngplotter, "set terminal png\n");

  for (t = 1; t <=nt; t++) {//Begin solve heat on t step
    time = t*dt;
    //calculating exact Solution
    for (i = 1; i <= np; i++) {
      ue[i] = exactfunction(npxy[i][1],npxy[i][2],time);
    }
    //creating stiffness matrix M for heat equation
    for (k = 1; k <= ne; k++) {
      //calculating phi, local RHS
      for (i = 1; i <= 3; i++) {
        for (j = 1; j <= 3; j++) {
          if (i==j) {
            localB[i] = (bfunction(npxy[elnp[k][i]][1],npxy[elnp[k][i]][2],nu,time)*areamatrix[k])/6.0;
            locn[i] = (uh[elnp[k][i]]*areamatrix[k])/(6.0*dt);
          } else {
            localB[i] = (bfunction(npxy[elnp[k][i]][1],npxy[elnp[k][i]][2],nu,time)*areamatrix[k])/12.0;
            locn[i] = (uh[elnp[k][i]]*areamatrix[k])/(12.0*dt);
          }
          //store local b to global b
          globalB[elnp[k][i]] += localB[i];
          globalN[elnp[k][i]] += locn[i];
        }
      }
    }
    //modifying RHS for heat equation
    for (i = 1; i <= np; i++) {
      heatRHS[i] = globalB[i] + globalN[i];
    }
    //imposing boundary
    for (i = 1; i <= nb; i++) {
      heatRHS[bdnp[i][1]] = 0.0;
    } //End Loop Each Element

    //Conjugate Gradient Method
    epsilon = 0.00000001;
    iter = 0;
    itermax = 300;
    deltanew = 0.0;
    //initialize adotx
    for (i = 1; i <= np; i++) {
      adotx[i] = 0.0;
    }
    //calculate A.x
    for (i = 1; i <= np; i++) {
      for (j = 1; j <= np; j++) {
        adotx[i] += heatLHS[i][j]*uh[j];
      }
    }
    //initialize residu and pk
    for (i = 1; i <= np; i++) {
      residu[i] = heatRHS[i]-adotx[i];
      pk[i] = residu[i];
    }
    for (i = 1; i <= np; i++) {
      deltanew += residu[i]*residu[i];
    }
    deltanol = deltanew;
    errcheck = epsilon*epsilon*deltanol;
    //starting conjugate gradient process
    while (iter < itermax && deltanew > errcheck) { //Begin CG Iteration
      //initialize adotpk
      for (i = 1; i <= np; i++) {
        adotpk[i] = 0.0;
      }
      //calculate A . pk
      for (i = 1; i <= np; i++) {
        for (j = 1; j <= np; j++) {
          adotpk[i] += heatLHS[i][j]*pk[j];
        }
      }
      //alpha2 initialize
      alpha2 = 0.0;
      for (i = 1; i <= np; i++) {
        alpha2 += pk[i]*adotpk[i];
      }
      alpha = (deltanew/alpha2);
      //updating solution value
      for (i = 1; i <= np; i++) {
        uh[i] = uh[i]+(alpha*pk[i]);
      }
      //initialize adotx
      for (i = 1; i <= np; i++) {
        adotx[i] = 0.0;
      }
      if (iter % 50 == 0) {
        //calculate A.x
        for (i = 1; i <= np; i++) {
          for (j = 1; j <= np; j++) {
            adotx[i] += heatLHS[i][j]*uh[j];
          }
        }
        //calculate new residu
        for (i = 1; i <= np; i++) {
          residu[i] = heatRHS[i] - adotx[i];
        }
      } else {
        //calculating new residu
        for (i = 1; i <= np; i++) {
          residu[i] = residu[i]-(alpha*adotpk[i]);
        }
      }
      deltaold = deltanew;
      deltanew = 0.0;
      for (i = 1; i <= np; i++) {
        deltanew += residu[i]*residu[i];
      }
      //calculating beta
      beta = (deltanew/deltaold);
      //updating pk
      for (i = 1; i <= np; i++) {
        pk[i] = residu[i]+(beta*pk[i]);
      }
      iter = iter + 1;
    } //End CG Iteration

    // //plotting 3D solution
    // //write solution to a text file
    // sprintf(fname, "./data/plotsolution_%d.txt",t);
    // sprintf(pictloc, "./animation/");
    // // puts(fname);
    // // printf("\n");
    // plotsolution = fopen(fname,"w");
    // if (plotsolution==NULL) {
    //   printf("Error Creating File %s\n",fname);
    //   exit(1);
    // }
    // //print and write the coordinate of the nodal point and solution to plotsolution.txt file
    // for(i=1; i<=ne; i++){
    //   j = elnp[i][1];
    //   k = elnp[i][2];
    //   l = elnp[i][3];
    //   fprintf(plotsolution, "%f\t %f\t %f\t %f\n", npxy[j][1], npxy[j][2], uh[j], ue[j]);
    //   fprintf(plotsolution, "%f\t %f\t %f\t %f\n", npxy[k][1], npxy[k][2], uh[k], ue[k]);
    //   fprintf(plotsolution, "%f\t %f\t %f\t %f\n", npxy[l][1], npxy[l][2], uh[l], ue[l]);
    //   fprintf(plotsolution, "%f\t %f\t %f\t %f\n\n\n", npxy[j][1], npxy[j][2], uh[j], ue[j]);
    // }
    // //3D plotting png using gpl file
    // sprintf(pngname,"%sheatplot_frame-%03d.png",pictloc,t);
    // fprintf(pngplotter, "set output \'%s\'\n",pngname);
    // fprintf(pngplotter, "set title '3D Heat Equation Overtime'\n");
    // fprintf(pngplotter, "set xlabel 'X'\n");
    // fprintf(pngplotter, "set ylabel 'Y'\n");
    // fprintf(pngplotter, "set xrange [0:1]\n");
    // fprintf(pngplotter, "set yrange [0:1]\n");
    // fprintf(pngplotter, "set zrange [0:1]\n");
    // fprintf(pngplotter, "set hidden3d\n");
    // fprintf(pngplotter, "unset key\n");
    // fprintf(pngplotter, "set key\n");
    // fprintf(pngplotter, "unset label\n");
    // fprintf(pngplotter, "set label \"t = %.1fs\" at 0.0, 0.0, 24.5\n",time);
    // fprintf(pngplotter, "set palette defined (0 'blue', 20 'red')\n");
    // fprintf(pngplotter, "splot \'%s\' u 1:2:3 w l lw 1 lc 1, \'%s\' u 1:2:4 w l lw 1 lc 3\n\n",fname,fname);
    // //close file plotsolution.txt
    // fclose(plotsolution);

    //initialize global RHS
    for (i = 1; i <= np; i++) {
      heatRHS[i] = 0.0;
      globalB[i] = 0.0;
      globalN[i] = 0.0;
    }
    // printf("Maximum Value Numeric = %f\n", maxValueVector(np,uh));
    // printf("Maximum Value Exact = %f\n", maxValueVector(np,ue));
  }//End of solving Heat equation
  // fclose(pngplotter);

  //preparing variable to calculate L2_norm and H10_norm
  //calculating exact Solution
  for (i = 1; i <= np; i++) {
    ue[i] = exactfunction(npxy[i][1],npxy[i][2],time);
  }
  //calculating v
  for (i = 1; i <= np; i++) {
    v[i] = uh[i] - ue[i];
  }
  //initialize M matriks
  for (i = 1; i <= 3; i++) {
    for (j = 1; j <= 3; j++) {
      if (i == j) {
        m[i][j] = 2.0;
      } else {
        m[i][j] = 1.0;
      }
    }
  }

  //calculating L2_norm
  //preparing variable for L2_norm
  double vmv,sigmavmv,vh;
  sigmavmv = 0.0;
  vh = 0.0;
  for (k = 1; k <= ne; k++) {
    //initialize mv
    for (i = 1; i <= 3; i++) {
      mv[i] = 0.0;
    }
    //resetting vmv
    vmv = 0.0;
    //calculating m.v
    for (i = 1; i <= 3; i++) {
      for (j = 1; j <= 3; j++) {
        mv[i] += m[i][j]*v[elnp[k][i]];
      }
      vmv += v[elnp[k][i]]*mv[i];
    }
    //calculating square root (vMv/12area)
    sigmavmv += (vmv*areamatrix[k])/12.0;
    vh = sqrt(sigmavmv);
  }
  printf("L2_norm of %d Division Number: %f\n", nb/4,vh);

  //calculating H10_norm
  double sumvhil,sumvh10,vh10;
  for (k = 1; k <= ne; k++) {
    //calculating c1 and c2
    for (i = 1; i <= 2; i++) {
      for (j = 1; j <= 3; j++) {
        xmatrix[i][j] = npxy[elnp[k][j]][i];
      }
    }
    meas = areamatrix[k];
    cmatrix[1][1] = (xmatrix[2][2]-xmatrix[2][3])/(2.0*meas);
    cmatrix[1][2] = (xmatrix[2][3]-xmatrix[2][1])/(2.0*meas);
    cmatrix[1][3] = (xmatrix[2][1]-xmatrix[2][2])/(2.0*meas);
    cmatrix[2][1] = (xmatrix[1][3]-xmatrix[1][2])/(2.0*meas);
    cmatrix[2][2] = (xmatrix[1][1]-xmatrix[1][3])/(2.0*meas);
    cmatrix[2][3] = (xmatrix[1][2]-xmatrix[1][1])/(2.0*meas);
    //initialize v.c
    for (i = 1; i <= 2; i++) {
      vhilbert[i] = 0.0;
    }
    //calculating v.c
    for (i = 1; i <= 2; i++) {
      for (j = 1; j <= 3; j++) {
        vhilbert[i] += v[elnp[k][j]]*cmatrix[i][j];
      }
    }
    //calculate norm vhilbert
    sumvhil = 0.0;
    for (i = 1; i <= 2; i++) {
      sumvhil += vhilbert[i]*vhilbert[i];
    }
    sumvh10 += sumvhil*areamatrix[k];
    vh10 = sqrt(sumvh10);
  }
  printf("H10_norm of %d Division Number: %f\n", nb/4,vh10);
  double decN = 1.0/(nb/4);
  printf("1/N = %f\n", decN);

  // //gnuplotting
  // char command1[100],encoder[100],command2[100];
  // sprintf(command1, "gnuplot %s",plottxtname);
  // system(command1);
  // sprintf(command2, "mencoder mf://animation/*.png -mf fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o movie.avi");
  // system(command2);
  //plotting log error
  fprintf(plotlogerror, "%f\t %f\t %f\n", decN, vh, vh10);
// fclose(plotlogerror);

//free the memory
free_imatrix(elnp,1,ne,1,DIM+1);
free_imatrix(bdnp,1,nb,1,DIM);
free_dmatrix(npxy,1,np,1,DIM);
free_dmatrix(globalA,1,np,1,np);
free_dvector(globalB,1);
free_dmatrix(localA,1,3,1,3);
free_dvector(localB,1);
free_dvector(areamatrix,1);
free_dmatrix(cmatrix,1,2,1,3);
free_dmatrix(xmatrix,1,2,1,3);
free_dvector(uh,1);
free_dvector(residu,1);
free_dvector(pk,1);
free_dvector(adotpk,1);
free_dvector(residuold,1);
free_dvector(adotx,1);
free_dvector(vhilbert,1);
free_dmatrix(globalM,1,np,1,np);
free_dvector(globalN,1);
free_dmatrix(m,1,3,1,3);
free_dvector(locn,1);
free_dmatrix(locm,1,3,1,3);
free_dmatrix(heatLHS,1,np,1,np);
free_dvector(heatRHS,1);
return 0;
}

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
