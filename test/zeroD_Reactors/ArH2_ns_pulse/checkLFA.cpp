#include<iostream>
#define eV 11604.25
#define LFA_NDATA 20
#define twothird 0.666666666
#define NUM_SPECIES 15
const double redf_ee_array[LFA_NDATA][2]=
{
    {2.4848,1.1964},
    {3.4393,1.4586},
    {4.6837,1.9177},
    {6.2755,2.5694},
    {8.8284,3.4426},
    {12.0226,4.1577},
    {16.3726,4.6562},
    {22.2964,5.0213},
    {30.8611,5.2639},
    {42.7158,5.6234},
    {59.1241,5.8397},
    {81.8354,6.1799},
    {111.4445,6.6645},
    {146.9125,7.1197},
    {203.3458,7.678},
    {281.4568,8.2023},
    {389.5724,9.1858},
    {505.2769,10.1905},
    {746.3477,12.1917},
    {968.0156,15.0048}
};

// compute electron temperature from EbyN using local field approximation
double etempLFA(const double Tg, double C[], const double EN) 
{
  int i=0;
  int loc=0;
  int found=0;
  double EE_interp=0.5*eV;

  for(i=0;i<LFA_NDATA;i++)
  {
    if(redf_ee_array[i][0]>EN)
    {
        found=1;
        break;
    }
  }

  loc=std::min(LFA_NDATA-1,i);
  std::cout<<"loc:"<<loc<<"\n";

  if(loc==0)
  {
      EE_interp=redf_ee_array[loc][1]
      +(redf_ee_array[loc+1][1]-redf_ee_array[loc][1])/(redf_ee_array[loc+1][0]-redf_ee_array[loc][0])*
      (EN-redf_ee_array[loc][0]);
  }
  else
  {
      EE_interp=redf_ee_array[loc][1]
      +(redf_ee_array[loc][1]-redf_ee_array[loc-1][1])/(redf_ee_array[loc][0]-redf_ee_array[loc-1][0])*
      (EN-redf_ee_array[loc][0]);
  }

  return(EE_interp*twothird*eV);
}
int main()
{
    double C[NUM_SPECIES]={0.0};
    double Tg=300.0;
    double EbyN;
    double etemp;

    std::cout<<"enter EbyN in Td:";
    std::cin>>EbyN;

    etemp=etempLFA(Tg,C,EbyN);
    std::cout<<"electron temp K,eV:"<<etemp<<"\t"<<etemp/eV<<"\n";
    std::cout<<"electron energy eV:"<<1.5*etemp/eV<<"\n";

    return(0);
}
