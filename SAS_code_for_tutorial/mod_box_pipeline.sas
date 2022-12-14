/* Trying to import dataset from url */

filename cardiac url "https://raw.githubusercontent.com/DylanDijk/Simon-Skene-Mixed-Models-Tutorial/main/Data_Images_Figures/The%20Cardiac%20Enzyme%20Data%20-%20Reduced%20Data%20set.csv";


proc import file=cardiac out=cardiac_data dbms=csv;
run;

proc print data=cardiac_data;
run;

/* Alternative importing from local */

proc import datafile="S:\Shared_Projects\RO06_Surrey_CTU_Statistics\Dylan\Small samples\R code\Mixed-Models-Small-Sample-Tutorial-bs4\Data_Images_Figures\The Cardiac Enzyme Data - Reduced Data set.csv"
        out=cardiac_data
        dbms=csv
        replace;
    
run;


proc mixed data=cardiac_data;
class dog trt time;
model atp=trt time trt*time/ ddfm=kr s;
repeated time/type=un subject=dog r=4;
run;




title2 "Inference Using Box Corrections";

proc iml;

m=12; /* number of subjects */
p=9;  /* number of time points */


use d var{atp};
read all;

y=atp;

use rmat var{col1 col2 col3 col4 col5 col6 col7 col8 col9};
read all;
rmat=col1||col2||col3||col4||col5||col6||col7||col8||col9;

SIGMA=rmat;


ic=j(m*p,1,1);

trt1=j(m*p/2,1,0)//j(m*p/2,1,1);
trt2=j(m*p/2,1,1)//j(m*p/2,1,0);
trt=trt2;

time=i(p);
do i=2 to m by 1;
	time=time//i(p);
    end;
time=time[,2:9];

int=j(m*p,(p-1),.);
do i=1 to (p-1);
   int[,i]=time[,i]#trt;
   end;

X=ic||trt||time||int;
Xr=ic||trt||time;

parm=ncol(X);c=ncol(X)-ncol(Xr);


tot=m*p;
ind=j(tot,1,1);
do i=1 to tot by 1;
	if  y[i]=. then ind[i]=0;
end;
nobs=sum(ind);

mobs=j(m,1,.);
do i=1 to m by 1;
	mobs[i]=sum(ind[(i-1)*p+1:(i*p)]);
end;

cnt=j(m,1,.);
do i=1 to m by 1;
	cnt[i]=sum(mobs[1:i]);
end;

yrem=y[loc(ind=1)];
Xrem=X[loc(ind=1),];
Xr_rem=Xr[loc(ind=1),];


start block_V(Mat) global(m);
	ans=i(m)@Mat;
	return(ans);
finish block_V;

start block_Vrem(Mat) global(m,p,cnt,mobs,nobs);
	ans=j(nobs,nobs,0);
	do i=1 to m by 1;
		if mobs[i]=p then Mati=Mat;
		else Mati=Mat[1:mobs[i],1:mobs[i]];
		if i=1 then ans[1:cnt[1],1:cnt[1]]=Mati;
		else ans[cnt[i-1]+1:cnt[i],cnt[i-1]+1:cnt[i]]=Mati;
	end;
	return(ans);
finish block_Vrem;

Vrem=block_Vrem(SIGMA);

/* ANOVA F-test*/

A=i(nobs)-Xrem*inv(t(Xrem)*Xrem)*t(Xrem);
B=Xrem*inv(t(Xrem)*Xrem)*t(Xrem)-Xr_rem*inv(t(Xr_rem)*Xr_rem)*t(Xr_rem);
F=(nobs-parm)#(t(yrem)*B*yrem)/(c#t(yrem)*A*yrem);

numdf=c;dendf=nobs-parm;

/* Box Correction */

psi=(nobs-parm)#trace(B*Vrem)/(c#trace(A*Vrem));
v1_BOX=((trace(B*Vrem))##2)/trace(B*Vrem*B*Vrem);
v2_BOX=((trace(A*Vrem))##2)/trace(A*Vrem*A*Vrem);

F_BOX=F/psi;

/* Modified Box Correction */

E=trace(B*Vrem)/trace(A*Vrem);
V=(trace(B*Vrem*B*Vrem)/((trace(B*Vrem))##2))+(trace(A*Vrem*A*Vrem)/((trace(A*Vrem))##2));
v1_MOD=c;
v2_MOD=(c#(4#V+1)-2)/(c#V-1);
lambda=((nobs-parm)/c)#((v2_MOD-2)/v2_MOD)#E;
*lambda=((nobs-parm)/c)#((v2-2)/v2)#E;

F_MOD=F/lambda;

prob_F=1-cdf("F",F,numdf,dendf);
prob_F_BOX=1-cdf("F",F_BOX,v1_BOX,v2_BOX);
prob_F_MOD=1-cdf("F",F_MOD,v1_MOD,v2_MOD);

/* Printing of results */

print / "ANOVA F Statistic";
print F numdf dendf prob_F;

print "Box Correction";
print psi F_BOX v1_BOX v2_BOX prob_F_BOX;

print "Modified Box Correction";
print lambda F_MOD v1_mod v2_mod prob_F_MOD;






proc import datafile="S:\Shared_Projects\RO06_Surrey_CTU_Statistics\Dylan\Small samples\R code\Mixed-Models-Small-Sample-Tutorial-bs4\Data_Images_Figures\The Cardiac Enzyme Data - Reduced Data set - Artificial missing.csv"
        out=cardiac_data
        dbms=csv
        replace;
    
run;

