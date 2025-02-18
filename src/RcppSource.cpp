#include <bits/stdc++.h>
using namespace std;

void merge(double arr[], double aux[], int arr_indx[], int aux_indx[],
           int low, int mid, int high, int tl_count[]){
  int k=low;
  int i=low;
  int j=mid+1;
  int count=0;
  while(i<=mid && j<=high)
  {
    if(arr[i]>arr[j])
    {
      tl_count[arr_indx[i]]+=count;
      aux_indx[k]=arr_indx[i];
      aux[k++]=arr[i++];
    }
    else{
      aux_indx[k]=arr_indx[j];
      aux[k++]=arr[j++];
      count++;
    }
  }
  while(i<=mid)
  {
    tl_count[arr_indx[i]]+=count;
    aux_indx[k]=arr_indx[i];
    aux[k++]=arr[i++];
  }
  for(int i=low; i<=high; i++){
    arr_indx[i]=aux_indx[i];
    arr[i]=aux[i];
  }
}

void merge_sort(double arr[], double aux[], int aux_indx[], int arr_indx[],
                int low, int high, int tl_count[]){
  if(high<=low){
    return;
  }
  int mid=(low+((high-low)>>1));
  merge_sort(arr, aux, aux_indx, arr_indx, low, mid, tl_count);
  merge_sort(arr, aux, aux_indx, arr_indx, mid+1, high, tl_count);
  merge(arr, aux, aux_indx, arr_indx, low, mid, high, tl_count);
}

void mergeC(double arr[], double aux[], int arr_indx[], int aux_indx[],
            int low, int mid, int high){
  int k=low;
  int i=low;
  int j=mid+1;
  while(i<=mid && j<=high)
  {
    if(arr[i]>arr[j])
    {
      aux_indx[k]=arr_indx[i];
      aux[k++]=arr[i++];
    }
    else{
      aux_indx[k]=arr_indx[j];
      aux[k++]=arr[j++];
    }
  }
  while(i<=mid)
  {
    aux_indx[k]=arr_indx[i];
    aux[k++]=arr[i++];
  }
  for(int i=low; i<=high; i++){
    arr_indx[i]=aux_indx[i];
    arr[i]=aux[i];
  }
}

void merge_sortC(double arr[], double aux[], int aux_indx[], int arr_indx[],
                 int low, int high){
  if(high<=low){
    return;
  }
  int mid=(low+((high-low)>>1));
  merge_sortC(arr, aux, aux_indx, arr_indx, low, mid);
  merge_sortC(arr, aux, aux_indx, arr_indx, mid+1, high);
  mergeC(arr, aux, aux_indx, arr_indx, low, mid, high);
}

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector trailcount(NumericVector &x){
  int n=x.length();
  double aux[n];
  double arr[n];
  int aux_indx[n];
  int arr_indx[n];
  int tl_count[n];
  for (int i=0; i<n; i++){
    aux[i]=-x[n-i-1];
    arr[i]=-x[n-i-1];
    aux_indx[i]=i;
    arr_indx[i]=i;
    tl_count[i]=0;
  }
  merge_sort(arr, aux, aux_indx, arr_indx, 0, n-1, tl_count);
  IntegerVector o(n);
  for (int i=0; i<n; i++){
    o[i]=1+tl_count[n-i-1];
  }
  return o;
}

// [[Rcpp::export]]
IntegerVector orderC(NumericVector &x){
  int n=x.length();
  double aux[n];
  double arr[n];
  int aux_indx[n];
  int arr_indx[n];
  for (int i=0; i<n; i++){
    aux[i]=x[i];
    arr[i]=x[i];
    aux_indx[i]=i;
    arr_indx[i]=i;
  }
  merge_sortC(arr, aux, aux_indx, arr_indx, 0, n-1);
  IntegerVector o(n);
  for (int i=0; i<n; i++){
    o[n-1-i]=arr_indx[i];
  }
  return o;
}

// [[Rcpp::export]]
IntegerVector trailCount_QxRx_ordered_abs_Y(LogicalVector &Qx_TF_ordered, NumericVector &ordered_abs_Y, IntegerVector &trailcount_ordered_abs_Y){
  int n=Qx_TF_ordered.length();
  NumericVector Qx_ordered_abs_Y=ordered_abs_Y[Qx_TF_ordered];
  IntegerVector trailcount_Qx_ordered_abs_Y=trailcount(Qx_ordered_abs_Y);
  NumericVector Rx_ordered_abs_Y=ordered_abs_Y[!Qx_TF_ordered];
  IntegerVector trailcount_Rx_ordered_abs_Y=trailcount(Rx_ordered_abs_Y);
  IntegerVector trailcount_QxRx_ordered_abs_Y(n);
  int l=0;
  int m=0;
  for(int k=0; k<n; k++){
    if(Qx_TF_ordered[k]==true){
      trailcount_QxRx_ordered_abs_Y[k]=trailcount_Qx_ordered_abs_Y[l++];
    }else{
      trailcount_QxRx_ordered_abs_Y[k]=trailcount_ordered_abs_Y[k]-trailcount_Rx_ordered_abs_Y[m++];
    }
  }
  return(trailcount_QxRx_ordered_abs_Y);
}


// [[Rcpp::export]]
double phi(double a, double b, double c, double d){
  double N=a*d-b*c;
  double D=(a+b)*(a+c)*(b+d)*(c+d);
  return sqrt(N*N/D);
}

// [[Rcpp::export]]
NumericVector Tk_phiC(NumericVector &x, NumericVector &y, int n, NumericMatrix &DMatX){
  NumericMatrix DMatY(n,n);
  NumericMatrix DMat(n,n);
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      double dx=DMatX(i,j);
      double dy=y[i]-y[j];
      DMatY(i,j)=dy;
      DMat(i,j)=dx*dx+dy*dy;
    }
  }
  IntegerMatrix NNMat(n,n-1);
  for(int i=0; i<n; i++){
    NumericVector dveci=DMat(i,_);
    IntegerVector oveci=orderC(dveci);
    oveci.erase(0);
    NNMat(i,_)=oveci;
  }
  NumericMatrix TMat(n,n-1);
  if(n<=60){
    for(int i=0; i<n; i++){
      for(int j0=0; j0<n-1; j0++){
        int j=NNMat(i,j0);
        double a=0;
        double b=0;
        double c=0;
        double d=0;
        for(int k=0; k<n; k++){
          if(x[k]!=x[i] && y[k]!=y[i]){
            double Dxki=DMatX(k,i);
            double Dyki=DMatY(k,i);
            if(abs(Dxki)<=abs(DMatX(j,i)) && abs(Dyki)<=abs(DMatY(j,i)) ){
              if(Dxki<0 && Dyki>0){
                a=a+1;
              }else if(Dxki>0 && Dyki>0){
                b=b+1;
              }else if(Dxki<0 && Dyki<0){
                c=c+1;
              }else if(Dxki>0 && Dyki<0){
                d=d+1;
              }
            }
          }
        }
        if(a+b && c+d && a+c && b+d){
          TMat(i,j0)=phi(a,b,c,d);
        }
      }
    }
  }else{
    for(int i=0; i<n; i++){
      NumericVector DXi=DMatX(_,i);
      NumericVector DYi=DMatY(_,i);
      NumericVector abs_X=abs(DXi);
      NumericVector abs_Y=abs(DYi);
      for(int k=0; k<n; k++){
        if(abs_X[k]*abs_Y[k]==0){
          abs_Y[k]=DBL_MAX;
        }
      }
      IntegerVector indx_abs_X=orderC(abs_X);
      NumericVector ordered_abs_Y=abs_Y[indx_abs_X];
      IntegerVector trailcount_ordered_abs_Y=trailcount(ordered_abs_Y);
      NumericVector DXi_ordered=DXi[indx_abs_X];
      NumericVector DYi_ordered=DYi[indx_abs_X];
      LogicalVector Q1_TF_ordered=(DXi_ordered>0) & (DYi_ordered>0);
      LogicalVector Q2_TF_ordered=(DXi_ordered<0) & (DYi_ordered>0);
      LogicalVector Q3_TF_ordered=(DXi_ordered<0) & (DYi_ordered<0);
      IntegerVector avec_ordered=trailCount_QxRx_ordered_abs_Y(Q2_TF_ordered,ordered_abs_Y,trailcount_ordered_abs_Y);
      IntegerVector bvec_ordered=trailCount_QxRx_ordered_abs_Y(Q1_TF_ordered,ordered_abs_Y,trailcount_ordered_abs_Y);
      IntegerVector cvec_ordered=trailCount_QxRx_ordered_abs_Y(Q3_TF_ordered,ordered_abs_Y,trailcount_ordered_abs_Y);
      IntegerVector dvec_ordered=trailcount_ordered_abs_Y-avec_ordered-bvec_ordered-cvec_ordered;
      NumericVector fvec_ordered(n);
      for(int k=0; k<n; k++){
        int a=avec_ordered[k];
        int b=bvec_ordered[k];
        int c=cvec_ordered[k];
        int d=dvec_ordered[k];
        if(ordered_abs_Y[k]==DBL_MAX){
          fvec_ordered[k]=0;
        }else if(a+b && c+d && a+c && b+d){
          fvec_ordered[k]=phi(a,b,c,d);
        }else{
          fvec_ordered[k]=0;
        }
      }
      IntegerVector inv_indx_abs_X(n);
      for(int k=0; k<n; k++){
        inv_indx_abs_X[indx_abs_X[k]]=k;
      }
      NumericVector fvec=fvec_ordered[inv_indx_abs_X];
      for(int k=0; k<n-1; k++){
        int j=NNMat(i,k);
        TMat(i,k)=fvec[j];
      }
    }
  }
  return colMeans(TMat);
}

// [[Rcpp::export]]
double T_phiC(NumericVector x, NumericVector y, int ITR){
  int n=x.length();
  NumericMatrix DMatX(n,n);
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      DMatX(i,j)=x[i]-x[j];
    }
  }
  NumericVector Tkvec=Tk_phiC(x,y,n,DMatX);
  NumericMatrix Tkmat(n-1,ITR);
  Function sampleC("sample.int");
  for(int itr=0;itr<ITR;itr++){
    IntegerVector indperm=sampleC(n);
    NumericVector yperm=y[indperm-1];
    NumericVector Tkvec_perm=Tk_phiC(x,yperm,n,DMatX);
    Tkmat(_,itr)=Tkvec_perm;
  }
  NumericVector Tkvec_mean=rowMeans(Tkmat);
  NumericVector Tkvec_sd(n-1);
  for(int i=0; i<n-1; i++){
    double sum=0;
    for(int j=0; j<ITR; j++){
      double df=Tkmat(i, j)-Tkvec_mean[i];
      sum+=df*df;
    }
    Tkvec_sd[i]=sqrt(sum/ITR);
  }
  bool all_sd_zero=std::all_of(Tkvec_sd.begin(), Tkvec_sd.end(), [](double x){ return x==0;});
  if(all_sd_zero){
    return 1.0;
  }else{
    NumericVector Tkvec_std(n-1);
    NumericMatrix Tkmat_std(n-1, ITR);
    for(int i=0; i<n-1; i++){
      if(Tkvec_sd[i]>0){
        Tkvec_std[i]=(Tkvec[i]-Tkvec_mean[i])/Tkvec_sd[i];
        for (int j=0; j<ITR; j++) {
          Tkmat_std(i, j)=(Tkmat(i, j)-Tkvec_mean[i])/Tkvec_sd[i];
        }
      }
    }
    double Tn=0;
    for(int i=0; i<n-1; i++){
      double tx=Tkvec_std[i];
      double tm=std::max(tx,0.0);
      Tn+=tm*tm;
    }
    NumericVector Tn_perms(ITR);
    for (int j=0; j<ITR; j++){
      double sum=0;
      for (int i=0; i<n-1; i++){
        double tx=Tkmat_std(i, j);
        double tm=std::max(tx,0.0);
        sum+=tm*tm;
      }
      Tn_perms[j]=sum;
    }
    double pvalue=std::count_if(Tn_perms.begin(), Tn_perms.end(), [Tn](double u) {
      return Tn<=u; })/static_cast<double>(ITR);
    return pvalue;
  }
}

// [[Rcpp::export]]
double corC(NumericVector x, NumericVector y) {
  double min=pow(2,-50);
  int n=x.size();
  double sum_x=0, sum_y=0, sum_xy=0;
  double sum_x_sq=0, sum_y_sq=0;
  for(int i=0; i<n; i++){
    sum_x+=x[i];
    sum_y+=y[i];
    sum_xy+=x[i]*y[i];
    sum_x_sq+=x[i]*x[i];
    sum_y_sq+=y[i]*y[i];
  }
  double numerator=n*sum_xy-sum_x*sum_y;
  double denominator_x=n*sum_x_sq-sum_x*sum_x;
  if(denominator_x<=min){
    return 0;
  }
  denominator_x=std::sqrt(denominator_x);
  double denominator_y=n*sum_y_sq-sum_y*sum_y;
  if(denominator_y<=min){
    return 0;
  }
  denominator_y=std::sqrt(denominator_y);
  double cor=numerator/denominator_x/denominator_y;
  return abs(cor);
}

// [[Rcpp::export]]
NumericVector Tk_corC(NumericVector &x, NumericVector &y, int n, NumericMatrix &DMatX){
  NumericMatrix DMatY(n,n);
  NumericMatrix DMat(n,n);
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      double dx=DMatX(i,j);
      double dy=abs(y[i]-y[j]);
      DMatY(i,j)=dy;
      DMat(i,j)=dx*dx+dy*dy;
    }
  }
  NumericMatrix TMat(n,n-1);
  for(int i=0; i<n; i++){
    NumericVector dveci=DMat(i,_);
    IntegerVector oveci=orderC(dveci);
    for(int j0=0; j0<n-1; j0++){
      int j=oveci[j0+1];
      if(x[i]!=x[j] && y[i]!=y[j]){
        LogicalVector points_in_rectangle(n);
        for(int k=0; k<n; k++){
          if(DMatX(k,i)<=DMatX(j,i) && DMatY(k,i)<=DMatY(j,i)){
            points_in_rectangle[k]=true;
          }
        }
        TMat(i,j0)=corC(x[points_in_rectangle], y[points_in_rectangle]);
      }
    }
  }
  return colMeans(TMat);
}

// [[Rcpp::export]]
double T_corC(NumericVector x, NumericVector y, int ITR){
  int n=x.length();
  NumericMatrix DMatX(n,n);
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      DMatX(i,j)=abs(x[i]-x[j]);
    }
  }
  NumericVector Tkvec=Tk_corC(x,y,n,DMatX);
  NumericMatrix Tkmat(n-1,ITR);
  Function sampleC("sample.int");
  for(int itr=0;itr<ITR;itr++){
    IntegerVector indperm=sampleC(n);
    NumericVector yperm=y[indperm-1];
    NumericVector Tkvec_perm=Tk_corC(x,yperm,n,DMatX);
    Tkmat(_,itr)=Tkvec_perm;
  }
  NumericVector Tkvec_mean=rowMeans(Tkmat);
  NumericVector Tkvec_sd(n-1);
  for(int i=0; i<n-1; i++){
    double sum=0;
    for(int j=0; j<ITR; j++){
      double df=Tkmat(i, j)-Tkvec_mean[i];
      sum+=df*df;
    }
    Tkvec_sd[i]=sqrt(sum/ITR);
  }
  bool all_sd_zero=std::all_of(Tkvec_sd.begin(), Tkvec_sd.end(), [](double x){ return x==0;});
  if(all_sd_zero){
    return 1.0;
  }else{
    NumericVector Tkvec_std(n-1);
    NumericMatrix Tkmat_std(n-1, ITR);
    for(int i=0; i<n-1; i++){
      if(Tkvec_sd[i]>0){
        Tkvec_std[i]=(Tkvec[i]-Tkvec_mean[i])/Tkvec_sd[i];
        for (int j=0; j<ITR; j++) {
          Tkmat_std(i, j)=(Tkmat(i, j)-Tkvec_mean[i])/Tkvec_sd[i];
        }
      }
    }
    double Tn=0;
    for(int i=0; i<n-1; i++){
      double tx=Tkvec_std[i];
      double tm=std::max(tx,0.0);
      Tn+=tm*tm;
    }
    NumericVector Tn_perms(ITR);
    for (int j=0; j<ITR; j++){
      double sum=0;
      for (int i=0; i<n-1; i++){
        double tx=Tkmat_std(i, j);
        double tm=std::max(tx,0.0);
        sum+=tm*tm;
      }
      Tn_perms[j]=sum;
    }
    double pvalue=std::count_if(Tn_perms.begin(), Tn_perms.end(), [Tn](double u) {
      return Tn<=u; })/static_cast<double>(ITR);
    return pvalue;
  }
}

// [[Rcpp::export]]
NumericVector Tk_dcorC(NumericVector &x, NumericVector &y, int n, NumericMatrix &DMatX){
  double min=pow(2,-50);
  NumericMatrix DMatY(n,n);
  NumericMatrix DMat(n,n);
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      double dx=DMatX(i,j);
      double dy=abs(y[i]-y[j]);
      DMatY(i,j)=dy;
      DMat(i,j)=dx*dx+dy*dy;
    }
  }
  NumericMatrix TMat(n,n-1);
  for(int i=0; i<n; i++){
    NumericVector dveci=DMat(i,_);
    IntegerVector oveci=orderC(dveci);
    for(int j0=0; j0<n-1; j0++){
      int j=oveci[j0+1];
      if(x[i]!=x[j] && y[i]!=y[j]){
        IntegerVector points_in_rectangle(0);
        for(int k=0; k<n; k++){
          if(DMatX(k,i)<=DMatX(j,i) && DMatY(k,i)<=DMatY(j,i)){
            points_in_rectangle.push_back(k);
          }
        }
        int m=points_in_rectangle.length();
        if(m>1){
          NumericMatrix Auv(m,m);
          NumericMatrix Buv(m,m);
          for(int u=0; u<m; u++){
            for(int v=0; v<m; v++){
              Auv(u,v)=DMatX(points_in_rectangle[u],points_in_rectangle[v]);
              Buv(u,v)=DMatY(points_in_rectangle[u],points_in_rectangle[v]);
            }
          }
          NumericVector Au_=colMeans(Auv);
          NumericVector Bu_=colMeans(Buv);
          double A__=mean(Au_);
          double B__=mean(Bu_);
          NumericMatrix A(m,m);
          NumericMatrix B(m,m);
          for(int u=0; u<m-1; u++){
            for(int v=u+1; v<m; v++){
              A(u,v)=Auv(u,v)-Au_[u]-Au_[v]+A__;
              B(u,v)=Buv(u,v)-Bu_[u]-Bu_[v]+B__;
            }
          }
          for(int u=0; u<m; u++){
            A(u,u)=-2*Au_[u]+A__;
            B(u,u)=-2*Bu_[u]+B__;
          }
          double dcov2xy=0;
          double dcov2xx=0;
          double dcov2yy=0;
          for(int u=0; u<m-1; u++){
            for(int v=u+1; v<m; v++){
              dcov2xy=dcov2xy+A(u,v)*B(u,v);
              dcov2xx=dcov2xx+A(u,v)*A(u,v);
              dcov2yy=dcov2yy+B(u,v)*B(u,v);
            }
          }
          dcov2xy=2*dcov2xy;
          dcov2xx=2*dcov2xx;
          dcov2yy=2*dcov2yy;
          for(int u=0; u<m; u++){
            dcov2xy=dcov2xy+A(u,u)*B(u,u);
            dcov2xx=dcov2xx+A(u,u)*A(u,u);
            dcov2yy=dcov2yy+B(u,u)*B(u,u);
          }
          if(dcov2xy>min && dcov2xx>min && dcov2yy>min){
            TMat(i,j0)=sqrt(dcov2xy/sqrt(dcov2xx*dcov2yy));
          }
        }
      }
    }
  }
  return colMeans(TMat);
}

// [[Rcpp::export]]
double T_dcorC(NumericVector x, NumericVector y, int ITR){
  int n=x.length();
  NumericMatrix DMatX(n,n);
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      DMatX(i,j)=abs(x[i]-x[j]);
    }
  }
  NumericVector Tkvec=Tk_dcorC(x,y,n,DMatX);
  NumericMatrix Tkmat(n-1,ITR);
  Function sampleC("sample.int");
  for(int itr=0;itr<ITR;itr++){
    IntegerVector indperm=sampleC(n);
    NumericVector yperm=y[indperm-1];
    NumericVector Tkvec_perm=Tk_dcorC(x,yperm,n,DMatX);
    Tkmat(_,itr)=Tkvec_perm;
  }
  NumericVector Tkvec_mean=rowMeans(Tkmat);
  NumericVector Tkvec_sd(n-1);
  for(int i=0; i<n-1; i++){
    double sum=0;
    for(int j=0; j<ITR; j++){
      double df=Tkmat(i, j)-Tkvec_mean[i];
      sum+=df*df;
    }
    Tkvec_sd[i]=sqrt(sum/ITR);
  }
  bool all_sd_zero=std::all_of(Tkvec_sd.begin(), Tkvec_sd.end(), [](double x){ return x==0;});
  if(all_sd_zero){
    return 1.0;
  }else{
    NumericVector Tkvec_std(n-1);
    NumericMatrix Tkmat_std(n-1, ITR);
    for(int i=0; i<n-1; i++){
      if(Tkvec_sd[i]>0){
        Tkvec_std[i]=(Tkvec[i]-Tkvec_mean[i])/Tkvec_sd[i];
        for (int j=0; j<ITR; j++) {
          Tkmat_std(i, j)=(Tkmat(i, j)-Tkvec_mean[i])/Tkvec_sd[i];
        }
      }
    }
    double Tn=0;
    for(int i=0; i<n-1; i++){
      double tx=Tkvec_std[i];
      double tm=std::max(tx,0.0);
      Tn+=tm*tm;
    }
    NumericVector Tn_perms(ITR);
    for (int j=0; j<ITR; j++){
      double sum=0;
      for (int i=0; i<n-1; i++){
        double tx=Tkmat_std(i, j);
        double tm=std::max(tx,0.0);
        sum+=tm*tm;
      }
      Tn_perms[j]=sum;
    }
    double pvalue=std::count_if(Tn_perms.begin(), Tn_perms.end(), [Tn](double u) {
      return Tn<=u; })/static_cast<double>(ITR);
    return pvalue;
  }
}

// [[Rcpp::export]]
NumericVector Tk_phi(NumericVector x, NumericVector y){
  int n=x.length();
  NumericMatrix DMatX(n,n);
  NumericMatrix DMatY(n,n);
  NumericMatrix DMat(n,n);
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      double dx=x[i]-x[j];
      double dy=y[i]-y[j];
      DMatX(i,j)=dx;
      DMatY(i,j)=dy;
      DMat(i,j)=dx*dx+dy*dy;
    }
  }
  IntegerMatrix NNMat(n,n-1);
  for(int i=0; i<n; i++){
    NumericVector dveci=DMat(i,_);
    IntegerVector oveci=orderC(dveci);
    oveci.erase(0);
    NNMat(i,_)=oveci;
  }
  NumericMatrix TMat(n,n-1);
  if(n<=60){
    for(int i=0; i<n; i++){
      for(int j0=0; j0<n-1; j0++){
        int j=NNMat(i,j0);
        double a=0;
        double b=0;
        double c=0;
        double d=0;
        for(int k=0; k<n; k++){
          if(x[k]!=x[i] && y[k]!=y[i]){
            double Dxki=DMatX(k,i);
            double Dyki=DMatY(k,i);
            if(abs(Dxki)<=abs(DMatX(j,i)) && abs(Dyki)<=abs(DMatY(j,i)) ){
              if(Dxki<0 && Dyki>0){
                a=a+1;
              }else if(Dxki>0 && Dyki>0){
                b=b+1;
              }else if(Dxki<0 && Dyki<0){
                c=c+1;
              }else if(Dxki>0 && Dyki<0){
                d=d+1;
              }
            }
          }
        }
        if(a+b && c+d && a+c && b+d){
          TMat(i,j0)=phi(a,b,c,d);
        }
      }
    }
  }
  else{
    for(int i=0; i<n; i++){
      NumericVector DXi=DMatX(_,i);
      NumericVector DYi=DMatY(_,i);
      NumericVector abs_X=abs(DXi);
      NumericVector abs_Y=abs(DYi);
      for(int k=0; k<n; k++){
        if(abs_X[k]*abs_Y[k]==0){
          abs_Y[k]=DBL_MAX;
        }
      }
      IntegerVector indx_abs_X=orderC(abs_X);
      NumericVector ordered_abs_Y=abs_Y[indx_abs_X];
      IntegerVector trailcount_ordered_abs_Y=trailcount(ordered_abs_Y);
      NumericVector DXi_ordered=DXi[indx_abs_X];
      NumericVector DYi_ordered=DYi[indx_abs_X];
      LogicalVector Q1_TF_ordered=(DXi_ordered>0) & (DYi_ordered>0);
      LogicalVector Q2_TF_ordered=(DXi_ordered<0) & (DYi_ordered>0);
      LogicalVector Q3_TF_ordered=(DXi_ordered<0) & (DYi_ordered<0);
      IntegerVector avec_ordered=trailCount_QxRx_ordered_abs_Y(Q2_TF_ordered,ordered_abs_Y,trailcount_ordered_abs_Y);
      IntegerVector bvec_ordered=trailCount_QxRx_ordered_abs_Y(Q1_TF_ordered,ordered_abs_Y,trailcount_ordered_abs_Y);
      IntegerVector cvec_ordered=trailCount_QxRx_ordered_abs_Y(Q3_TF_ordered,ordered_abs_Y,trailcount_ordered_abs_Y);
      IntegerVector dvec_ordered=trailcount_ordered_abs_Y-avec_ordered-bvec_ordered-cvec_ordered;
      NumericVector fvec_ordered(n);
      for(int k=0; k<n; k++){
        int a=avec_ordered[k];
        int b=bvec_ordered[k];
        int c=cvec_ordered[k];
        int d=dvec_ordered[k];
        if(ordered_abs_Y[k]==DBL_MAX){
          fvec_ordered[k]=0;
        }else if(a+b && c+d && a+c && b+d){
          fvec_ordered[k]=phi(a,b,c,d);
        }else{
          fvec_ordered[k]=0;
        }
      }
      IntegerVector inv_indx_abs_X(n);
      for(int k=0; k<n; k++){
        inv_indx_abs_X[indx_abs_X[k]]=k;
      }
      NumericVector fvec=fvec_ordered[inv_indx_abs_X];
      for(int k=0; k<n-1; k++){
        int j=NNMat(i,k);
        TMat(i,k)=fvec[j];
      }
    }
  }
  return colMeans(TMat);
}

// [[Rcpp::export]]
NumericVector Tk_cor(NumericVector x, NumericVector y){
  int n=x.length();
  NumericMatrix DMatX(n,n);
  NumericMatrix DMatY(n,n);
  NumericMatrix DMat(n,n);
  for(int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      DMatX(i,j)=abs(x[i]-x[j]);
      DMatY(i,j)=abs(y[i]-y[j]);
      DMat(i,j)=DMatX(i,j)*DMatX(i,j)+DMatY(i,j)*DMatY(i,j);
      DMatX(j,i)=DMatX(i,j);
      DMatY(j,i)=DMatY(i,j);
      DMat(j,i)=DMat(i,j);
    }
  }
  NumericMatrix TMat(n,n-1);
  for(int i=0; i<n; i++){
    NumericVector dveci=DMat(i,_);
    IntegerVector oveci=orderC(dveci);
    for(int j0=0; j0<n-1; j0++){
      int j=oveci[j0+1];
      if(x[i]!=x[j] && y[i]!=y[j]){
        LogicalVector points_in_rectangle(n);
        for(int k=0; k<n; k++){
          if(DMatX(k,i)<=DMatX(j,i) && DMatY(k,i)<=DMatY(j,i)){
            points_in_rectangle[k]=true;
          }
        }
        TMat(i,j0)=corC(x[points_in_rectangle], y[points_in_rectangle]);
      }
    }
  }
  return colMeans(TMat);
}

// [[Rcpp::export]]
NumericVector Tk_dcor(NumericVector x, NumericVector y){
  double min=pow(2,-50);
  int n=x.length();
  NumericMatrix DMatX(n,n);
  NumericMatrix DMatY(n,n);
  NumericMatrix DMat(n,n);
  for(int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      DMatX(i,j)=abs(x[i]-x[j]);
      DMatY(i,j)=abs(y[i]-y[j]);
      DMat(i,j)=DMatX(i,j)*DMatX(i,j)+DMatY(i,j)*DMatY(i,j);
      DMatX(j,i)=DMatX(i,j);
      DMatY(j,i)=DMatY(i,j);
      DMat(j,i)=DMat(i,j);
    }
  }
  NumericMatrix TMat(n,n-1);
  for(int i=0; i<n; i++){
    NumericVector dveci=DMat(i,_);
    IntegerVector oveci=orderC(dveci);
    for(int j0=0; j0<n-1; j0++){
      int j=oveci[j0+1];
      if(x[i]!=x[j] && y[i]!=y[j]){
        IntegerVector points_in_rectangle(0);
        for(int k=0; k<n; k++){
          if(DMatX(k,i)<=DMatX(j,i) && DMatY(k,i)<=DMatY(j,i)){
            points_in_rectangle.push_back(k);
          }
        }
        int m=points_in_rectangle.length();
        if(m>1){
          NumericMatrix Auv(m,m);
          NumericMatrix Buv(m,m);
          for(int u=0; u<m; u++){
            for(int v=0; v<m; v++){
              Auv(u,v)=DMatX(points_in_rectangle[u],points_in_rectangle[v]);
              Buv(u,v)=DMatY(points_in_rectangle[u],points_in_rectangle[v]);
            }
          }
          NumericVector Au_=colMeans(Auv);
          NumericVector Bu_=colMeans(Buv);
          double A__=mean(Au_);
          double B__=mean(Bu_);
          NumericMatrix A(m,m);
          NumericMatrix B(m,m);
          for(int u=0; u<m-1; u++){
            for(int v=u+1; v<m; v++){
              A(u,v)=Auv(u,v)-Au_[u]-Au_[v]+A__;
              B(u,v)=Buv(u,v)-Bu_[u]-Bu_[v]+B__;
            }
          }
          for(int u=0; u<m; u++){
            A(u,u)=-2*Au_[u]+A__;
            B(u,u)=-2*Bu_[u]+B__;
          }
          double dcov2xy=0;
          double dcov2xx=0;
          double dcov2yy=0;
          for(int u=0; u<m-1; u++){
            for(int v=u+1; v<m; v++){
              dcov2xy=dcov2xy+A(u,v)*B(u,v);
              dcov2xx=dcov2xx+A(u,v)*A(u,v);
              dcov2yy=dcov2yy+B(u,v)*B(u,v);
            }
          }
          dcov2xy=2*dcov2xy;
          dcov2xx=2*dcov2xx;
          dcov2yy=2*dcov2yy;
          for(int u=0; u<m; u++){
            dcov2xy=dcov2xy+A(u,u)*B(u,u);
            dcov2xx=dcov2xx+A(u,u)*A(u,u);
            dcov2yy=dcov2yy+B(u,u)*B(u,u);
          }
          if(dcov2xy>min && dcov2xx>min && dcov2yy>min){
            TMat(i,j0)=sqrt(dcov2xy/sqrt(dcov2xx*dcov2yy));
          }
        }
      }
    }
  }
  return colMeans(TMat);
}
