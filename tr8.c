#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void matrixmultiplication(float av[3],float bv[3][3],float outputc[3]);

void crossproduct(float axv,float ayv,float azv,float bxv,float byv,float bzv,float outputc[3]);

void dotproduct(float axv,float ayv,float azv,float bxv,float byv,float bzv,float *outputc);

int main(int argc,char *argv[])
{
	FILE *fpt1,*fpt2;
	int i=0,j=0,ww=0;
	fpt1=fopen(argv[1],"rb");
	fpt2=fopen("output.ppm","wb");
	fprintf(fpt2,"P5 256 256 255\n");
	unsigned char chrc[15];
	int vertex,face,xyzz;
	float anglex,angley,anglez;

	anglex=atoi(argv[2]);
	angley=atoi(argv[3]);
	anglez=atoi(argv[4]);


	while(strcmp(chrc,"list")!=0)
	{
		if(strcmp(chrc,"vertex")==0)
			{
			fscanf(fpt1,"%d",&vertex);
			}
		else if(strcmp(chrc,"face")==0)
			{
			fscanf(fpt1,"%d",&face);
			}
		fscanf(fpt1,"%s",chrc);
		
	}
	
	while(strcmp(chrc,"end_header")!=0)
	{
		fscanf(fpt1,"%s",chrc);
	}

	float xv[vertex-1],yv[vertex-1],zv[vertex-1];
	int xf[face-1],yf[face-1],zf[face-1];

	for(i=0;i<vertex;i++)
	{
		fscanf(fpt1,"%f",&xv[i]);
		fscanf(fpt1,"%f",&yv[i]);
		fscanf(fpt1,"%f",&zv[i]);
	}
	
	for(i=0;i<face;i++)
	{
		fscanf(fpt1,"%d",&xyzz);
		fscanf(fpt1,"%d",&xf[i]);
		fscanf(fpt1,"%d",&yf[i]);
		fscanf(fpt1,"%d",&zf[i]);
	}
	
	float minx=0,miny=0,minz=0,maxx=0,maxy=0,maxz=0;

	minx=xv[0];
	miny=yv[0];
	minz=zv[0];
	maxx=xv[0];
	maxy=yv[0];
	maxz=zv[0];


	for(i=0;i<vertex;i++)
	{
		if(xv[i]<minx)
		{
			minx=xv[i];
		}
	}
	for(i=0;i<vertex;i++)
	{
		if(yv[i]<miny)
		{
			miny=yv[i];
		}
	}
	for(i=0;i<vertex;i++)
	{
		if(zv[i]<minz)
		{
			minz=zv[i];
		}
	}
	for(i=0;i<vertex;i++)
	{
		if(xv[i]>maxx)
		{
			maxx=xv[i];
		}
	}
	for(i=0;i<vertex;i++)
	{
		if(yv[i]>maxy)
		{
			maxy=yv[i];
		}
	}
	for(i=0;i<vertex;i++)
	{
		if(zv[i]>maxz)
		{
			maxz=zv[i];
		}
	}
	
	float centerx,centery,centerz;

	centerx=(minx+maxx)/2;
	centery=(miny+maxy)/2;
	centerz=(minz+maxz)/2;

	printf("centerx%f centery%f centerz%f\n",centerx,centery,centerz);

	float extentx,extenty,extentz,maxextentxyz;

	extentx=maxx-minx;
	extenty=maxy-miny;
	extentz=maxz-minz;

	if(extentx>extenty && extentx>extentz)
	{
		maxextentxyz=extentx;
	}
	else if(extenty>extentx && extenty>extentz)
	{
		maxextentxyz=extenty;
	}
	else if(extentz>extenty && extentz>extentx)
	{
		maxextentxyz=extentz;
	}

	printf("maxextent=%f\n",maxextentxyz);

	
	float rotatex[3][3],rotatey[3][3],rotatez[3][3];

	rotatex[0][0]=1;
	rotatex[0][1]=0;
	rotatex[0][2]=0;
	rotatex[1][0]=0;
	rotatex[1][1]=cos((3.14*anglex)/180.0);
	rotatex[1][2]=-sin((3.14*anglex)/180.0);
	rotatex[2][0]=0;
	rotatex[2][1]=sin((3.14*anglex)/180.0);
	rotatex[2][2]=cos((3.14*anglex)/180.0);

	
	rotatey[0][0]=cos((3.14*angley)/180.0);
	rotatey[0][1]=0;
	rotatey[0][2]=sin((3.14*angley)/180.0);
	rotatey[1][0]=0;
	rotatey[1][1]=1;
	rotatey[1][2]=0;
	rotatey[2][0]=-sin((3.14*angley)/180.0);
	rotatey[2][1]=0;
	rotatey[2][2]=cos((3.14*angley)/180.0);

	
	rotatez[0][0]=cos((3.14*anglez)/180.0);
	rotatez[0][1]=-sin((3.14*anglez)/180.0);
	rotatez[0][2]=0;
	rotatez[1][0]=sin((3.14*anglez)/180.0);
	rotatez[1][1]=cos((3.14*anglez)/180.0);
	rotatez[1][2]=0;
	rotatez[2][0]=0;
	rotatez[2][1]=0;
	rotatez[2][2]=1;

	float camerax=1,cameray=0,cameraz=0;
	float upx=0,upy=0,upz=1;
	float a;
	int rr,cc;
	int COLS=256,ROWS=256;
	int rrcc=ROWS*COLS;
	float outimage[rrcc][3];
	unsigned char pixels[rrcc];
	float A,B,C,D;
	float eqn[3];
	float n,d;
	float intersect[3];
	float dot1,dot2,dot3;
	float qwerty1[3],qwerty2[3];
	float camera[3];
	float up[3];
	float tempcam1[3],tempcam2[3];
	float left[3],right[3],top[3],bottom[3],topleft[3];

	camera[0]=camerax;
	camera[1]=cameray;
	camera[2]=cameraz;

	matrixmultiplication(camera,rotatex,tempcam1);
	matrixmultiplication(tempcam1,rotatey,tempcam2);
	matrixmultiplication(tempcam2,rotatez,camera);

	up[0]=upx;
	up[1]=upy;
	up[2]=upz;
	
	camerax=camera[0];
	cameray=camera[1];
	cameraz=camera[2];

	matrixmultiplication(up,rotatex,tempcam1);
	matrixmultiplication(tempcam1,rotatey,tempcam2);
	matrixmultiplication(tempcam2,rotatez,up);

	upx=up[0];
	upy=up[1];
	upz=up[2];

	camerax=((1.5*maxextentxyz*camerax)+centerx);
	cameray=((1.5*maxextentxyz*cameray)+centery);
	cameraz=((1.5*maxextentxyz*cameraz)+centerz);

	crossproduct(upx,upy,upz,(centerx-camerax),(centery-cameray),(centerz-cameraz),left);

	a=sqrt((left[0]*left[0])+(left[1]*left[1])+(left[2]*left[2]));

	left[0]=(((maxextentxyz/(2*a))*left[0])+centerx);
	left[1]=(((maxextentxyz/(2*a))*left[1])+centery);
	left[2]=(((maxextentxyz/(2*a))*left[2])+centerz);

	crossproduct((centerx-camerax),(centery-cameray),(centerz-cameraz),upx,upy,upz,right);
	
	right[0]=(((maxextentxyz/(2*a))*right[0])+centerx);
	right[1]=(((maxextentxyz/(2*a))*right[1])+centery);
	right[2]=(((maxextentxyz/(2*a))*right[2])+centerz);

	top[0]=(((maxextentxyz/(2.0))*upx)+centerx);
	top[1]=(((maxextentxyz/(2.0))*upy)+centery);
	top[2]=(((maxextentxyz/(2.0))*upz)+centerz);

	bottom[0]=(((-maxextentxyz/(2.0))*upx)+centerx);
	bottom[1]=(((-maxextentxyz/(2.0))*upy)+centery);
	bottom[2]=(((-maxextentxyz/(2.0))*upz)+centerz);

	topleft[0]=(((maxextentxyz/(2.0))*upx)+left[0]);
	topleft[1]=(((maxextentxyz/(2.0))*upy)+left[1]);
	topleft[2]=(((maxextentxyz/(2.0))*upz)+left[2]);

	printf("camx=%f camy=%f camz=%f\n",camerax,cameray,cameraz);
	printf("left0=%f left1=%f left2=%f\n",left[0],left[1],left[2]);
	printf("right0=%f right1=%f right2=%f\n",right[0],right[1],right[2]);
	printf("top0=%f top1=%f top2=%f\n",top[0],top[1],top[2]);
	printf("bottom0=%f bottom1=%f bottom2=%f\n",bottom[0],bottom[1],bottom[2]);
	printf("tl0=%f tl1=%f tl2=%f\n",topleft[0],topleft[1],topleft[2]);
	
	for(i=0;i<rrcc;i++)
	{
		pixels[i]=0;
	}

	for(rr=0;rr<(ROWS);rr++)
	{
		for(cc=0;cc<(COLS);cc++)
		{
			outimage[(int)((float)(rr*COLS)+(float)cc)][0]=topleft[0]+((float)cc/(float)(COLS-1))*(right[0]-left[0])+((float)rr/((float)(ROWS-1)))*(bottom[0]-top[0]);
			outimage[(int)((float)(rr*COLS)+(float)cc)][1]=topleft[1]+((float)cc/(float)(COLS-1))*(right[1]-left[1])+((float)rr/((float)(ROWS-1)))*(bottom[1]-top[1]);
			outimage[(int)((float)(rr*COLS)+(float)cc)][2]=topleft[2]+((float)cc/(float)(COLS-1))*(right[2]-left[2])+((float)rr/((float)(ROWS-1)))*(bottom[2]-top[2]);
		}
	}

	float aaa,bbb,ccc,ddd,eee,fff;
	float absol;
	float zbuffer=999999;
	float buffer=0.0;
	
	for(ww=0;ww<face;ww++)
	{
			aaa=(xv[yf[ww]]-xv[xf[ww]]);
			bbb=(yv[yf[ww]]-yv[xf[ww]]);
			ccc=(zv[yf[ww]]-zv[xf[ww]]);
			ddd=(xv[zf[ww]]-xv[xf[ww]]);
			eee=(yv[zf[ww]]-yv[xf[ww]]);
			fff=(zv[zf[ww]]-zv[xf[ww]]);

		crossproduct(aaa , bbb , ccc , ddd , eee , fff , eqn);
		A=eqn[0];
		B=eqn[1];
		C=eqn[2];

		dotproduct(-A,-B,-C,xv[xf[ww]],yv[xf[ww]],zv[xf[ww]],&D);
		
		dotproduct(-A,-B,-C,(camerax),(cameray),(cameraz),&n);

		n=n-D;
		buffer=n/D;


		for(i=0;i<rrcc;i++)
		{
				
			dotproduct(A,B,C,(outimage[i][0]-camerax),(outimage[i][1]-cameray),(outimage[i][2]-cameraz),&d);

			absol=abs(d);
			
			if(absol<0.000001)
				{
				continue;
				}

			intersect[0]=camerax+(n/d)*(outimage[i][0]-camerax);
			intersect[1]=cameray+(n/d)*(outimage[i][1]-cameray);
			intersect[2]=cameraz+(n/d)*(outimage[i][2]-cameraz);

				aaa=(xv[zf[ww]]-xv[xf[ww]]);
				bbb=(yv[zf[ww]]-yv[xf[ww]]);
				ccc=(zv[zf[ww]]-zv[xf[ww]]);
				ddd=(xv[yf[ww]]-xv[xf[ww]]);
				eee=(yv[yf[ww]]-yv[xf[ww]]);
				fff=(zv[yf[ww]]-zv[xf[ww]]);

			crossproduct(aaa , bbb , ccc , ddd , eee , fff , qwerty1);

				aaa=(intersect[0]-xv[xf[ww]]);
				bbb=(intersect[1]-yv[xf[ww]]);
				ccc=(intersect[2]-zv[xf[ww]]);
				ddd=(xv[yf[ww]]-xv[xf[ww]]);
				eee=(yv[yf[ww]]-yv[xf[ww]]);
				fff=(zv[yf[ww]]-zv[xf[ww]]);

			crossproduct(aaa,bbb,ccc,ddd,eee,fff,qwerty2);
			dotproduct(qwerty1[0] , qwerty1[1] , qwerty1[2] , qwerty2[0] , qwerty2[1] , qwerty2[2] ,  &dot1);

				aaa=(xv[xf[ww]]-xv[yf[ww]]);
				bbb=(yv[xf[ww]]-yv[yf[ww]]);
				ccc=(zv[xf[ww]]-zv[yf[ww]]);
				ddd=(xv[zf[ww]]-xv[yf[ww]]);
				eee=(yv[zf[ww]]-yv[yf[ww]]);
				fff=(zv[zf[ww]]-zv[yf[ww]]);

			crossproduct(aaa , bbb , ccc , ddd , eee , fff , qwerty1);

				aaa=(intersect[0]-xv[yf[ww]]);
				bbb=(intersect[1]-yv[yf[ww]]);
				ccc=(intersect[2]-zv[yf[ww]]);
				ddd=(xv[zf[ww]]-xv[yf[ww]]);
				eee=(yv[zf[ww]]-yv[yf[ww]]);
				fff=(zv[zf[ww]]-zv[yf[ww]]);

			crossproduct(aaa , bbb , ccc , ddd , eee , fff , qwerty2);
			dotproduct(qwerty1[0] , qwerty1[1] , qwerty1[2] , qwerty2[0] , qwerty2[1] , qwerty2[2] , &dot2);

				aaa=(xv[yf[ww]]-xv[zf[ww]]);
				bbb=(yv[yf[ww]]-yv[zf[ww]]);
				ccc=(zv[yf[ww]]-zv[zf[ww]]);
				ddd=(xv[xf[ww]]-xv[zf[ww]]);
				eee=(yv[xf[ww]]-yv[zf[ww]]);
				fff=(zv[xf[ww]]-zv[zf[ww]]);

			crossproduct(aaa , bbb , ccc , ddd , eee , fff , qwerty1);

				aaa=(intersect[0]-xv[zf[ww]]);
				bbb=(intersect[1]-yv[zf[ww]]);
				ccc=(intersect[2]-zv[zf[ww]]);
				ddd=(xv[xf[ww]]-xv[zf[ww]]);
				eee=(yv[xf[ww]]-yv[zf[ww]]);
				fff=(zv[xf[ww]]-zv[zf[ww]]);

			crossproduct(aaa , bbb , ccc , ddd , eee , fff , qwerty2);
			dotproduct(qwerty1[0] , qwerty1[1] , qwerty1[2] , qwerty2[0] , qwerty2[1] , qwerty2[2] , &dot3);

			if((dot1>0) && (dot2>0) && (dot3>0))
			{
				pixels[i]=(float)(155 + (ww % 100));
			}
			else
			{
				continue;
			}

			if(buffer>zbuffer)
			{
				continue;
			}
	
		}
		
	}
	
	fwrite(pixels,1,rrcc,fpt2);
	

	fclose(fpt1);
	fclose(fpt2);
}

void matrixmultiplication(float av[3],float bv[3][3],float outputc[3])
{

        outputc[0]=(av[0]*bv[0][0])+(av[1]*bv[1][0])+(av[2]*bv[2][0]);
        outputc[1]=(av[0]*bv[0][1])+(av[1]*bv[1][1])+(av[2]*bv[2][1]);
        outputc[2]=(av[0]*bv[0][2])+(av[1]*bv[1][2])+(av[2]*bv[2][2]);

}


void crossproduct(float axv,float ayv,float azv,float bxv,float byv,float bzv,float outputc[3])
{
        outputc[0] = ((ayv*bzv)-(azv*byv));
        outputc[1] = ((azv*bxv)-(axv*bzv));
        outputc[2] = ((axv*byv)-(ayv*bxv));

}


void dotproduct(float axv,float ayv,float azv,float bxv,float byv,float bzv,float *outputc)
{
        *outputc=((axv*bxv)+(ayv*byv)+(azv*bzv));
}

