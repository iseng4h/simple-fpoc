#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#define IMAGE_WIDTH   	256
#define IMAGE_HEIGHT  	256
#define BAND_WIDTH    	128
#define BAND_HEIGHT   	128
#define NUMANGLE	2000

int main(int argc,char **argv) {

  int i, w, x, y, z;
  int zs, ze;
  int bv, bx, by, bz;
  int sx, sy;
  unsigned char *inbuf;
  FILE *fp;
  char fn[100];
  fftw_plan plan;
  fftw_complex *im_i, *phase;
  double *inphase, *refphase, difphase;
  double smax, smax1, smax2, stmp[NUMANGLE];

  //angle
  zs=0;
  ze=NUMANGLE;
  inbuf = malloc(IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(unsigned char));

  //read image
  if((fp=fopen(argv[1],"rb"))==NULL){
    printf("file not found");
    exit(1);
  }
  fread(inbuf,sizeof(unsigned char),IMAGE_WIDTH*IMAGE_HEIGHT,fp);
  fclose(fp);

  for (i=0;i<IMAGE_WIDTH*IMAGE_HEIGHT;i++){
    if(inbuf[i]>9)
      inbuf[i]=255;
  }
  
  //Generate Input Image
  im_i=fftw_malloc(IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(fftw_complex));
  for(y=0;y<IMAGE_HEIGHT;y++){
    for(x=0;x<IMAGE_WIDTH;x++){
      im_i[y*IMAGE_WIDTH+x][0]=(double)inbuf[y*IMAGE_WIDTH+x];
      im_i[y*IMAGE_WIDTH+x][1]=0.0;
    }
  }
  free(inbuf);

  //FFT 2D sensor image
  plan = fftw_plan_dft_2d(IMAGE_WIDTH, IMAGE_HEIGHT, im_i, im_i, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);

  //Calculate phase of sensor image
  inphase=malloc(IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(double));
  for(y=0;y<IMAGE_HEIGHT;y++){
    for(x=0;x<IMAGE_WIDTH;x++){
      inphase[y*IMAGE_WIDTH+x]=atan2(im_i[y*IMAGE_WIDTH+x][1], im_i[y*IMAGE_WIDTH+x][0]);
    }
  }
  fftw_free(im_i);

  w=0;
  refphase=malloc(IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(double));
  phase=fftw_malloc(BAND_WIDTH*BAND_HEIGHT*sizeof(fftw_complex));

  for(z=zs;z<ze;z++) {

    /******* Read database *******/
    //sprintf(fn, "after.raw", z);

    if((fp=fopen(argv[1], "rb"))==NULL){
      printf("No such file (%s)\n", fn);
      exit(1);
    }
    fread(refphase, sizeof(double), IMAGE_WIDTH*IMAGE_HEIGHT, fp);
    fclose(fp);

    /******* Calculate differential of phase *******/
    /* Bottom left */
    by=(int)(BAND_HEIGHT/2);
    bx=(int)(BAND_WIDTH/2);
    for(y=0;y<=by;y++){
      for(x=0;x<=bx;x++){
        difphase=refphase[y*IMAGE_WIDTH+x]-inphase[y*IMAGE_WIDTH+x];
        phase[y*BAND_WIDTH+x][0]=cos(difphase);
        phase[y*BAND_WIDTH+x][1]=sin(difphase);
      }
    }

    /* Bottom right */
    by=(int)(BAND_HEIGHT/2);
    bx=IMAGE_WIDTH-(int)(BAND_WIDTH/2);
    bv=(int)(BAND_WIDTH/2)+1;
    for(y=0;y<=by;y++){
      for(x=bx;x<IMAGE_WIDTH;x++){
        difphase=refphase[y*IMAGE_WIDTH+x]-inphase[y*IMAGE_WIDTH+x];
        phase[y*BAND_WIDTH+(x-bx+bv)][0]=cos(difphase);
        phase[y*BAND_WIDTH+(x-bx+bv)][1]=sin(difphase);
      }
    }
							    
    /* Top left */
    by=IMAGE_HEIGHT-(int)(BAND_HEIGHT/2);
    bz=(int)(BAND_HEIGHT/2)+1;
    bx=(int)(BAND_WIDTH/2);
    for(y=by;y<IMAGE_HEIGHT;y++){
      for(x=0;x<=bx;x++){
        difphase=refphase[y*IMAGE_WIDTH+x]-inphase[y*IMAGE_WIDTH+x];
        phase[(y-by+bz)*BAND_WIDTH+x][0]=cos(difphase);
        phase[(y-by+bz)*BAND_WIDTH+x][1]=sin(difphase);
      }
    }
							    
    /* Top right */
    by=IMAGE_HEIGHT-(int)(BAND_HEIGHT/2);
    bz=(int)(BAND_HEIGHT/2)+1;
    bx=IMAGE_WIDTH-(int)(BAND_WIDTH/2);
    bv=(int)(BAND_WIDTH/2)+1;
    for(y=by;y<IMAGE_HEIGHT;y++){
      for(x=bx;x<IMAGE_WIDTH;x++){
        difphase=refphase[y*IMAGE_WIDTH+x]-inphase[y*IMAGE_WIDTH+x];
        phase[(y-by+bz)*BAND_WIDTH+(x-bx+bv)][0]=cos(difphase);
        phase[(y-by+bz)*BAND_WIDTH+(x-bx+bv)][1]=sin(difphase);
      }
    }

    /******* IFFT the phase (matching result) *******/
    plan = fftw_plan_dft_2d(BAND_WIDTH, BAND_HEIGHT, phase, phase, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    /******* Calculate score *******/
    smax1=smax2=0.0;
    for(sx=0;sx<BAND_WIDTH;sx++){
      for(sy=0;sy<BAND_HEIGHT;sy++){
        if(phase[sy+BAND_WIDTH*sx][0]>smax2){
          if(phase[sy+BAND_WIDTH*sx][0]>smax1){
            smax1=phase[sy+BAND_WIDTH*sx][0];
          }
          else{
            smax2=phase[sy+BAND_WIDTH*sx][0];
          }
        }
      }
    }
    stmp[w]=(smax1+smax2)/(double)(BAND_HEIGHT*BAND_WIDTH);
    w++;

  } /******* end loop of angle *******/
  
  /******* Free memories *******/
  free(inphase);
  free(refphase);
  fftw_free(phase);
  fftw_destroy_plan(plan);

  /******* Finding local maximum score *******/
  smax=stmp[0];
  for(sx=1; sx<w; sx++){
    if(smax<stmp[sx]){
      smax=stmp[sx];
    }
  }
		    
  printf("Maximum score = %5.5f\n",smax);


  return 0;

}
