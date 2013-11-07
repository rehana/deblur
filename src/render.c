/* GIMP Plug-in Template
 * Copyright (C) 2000  Michael Natterer <mitch@gimp.org> (the "Author").
 * All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHOR BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * Except as contained in this notice, the name of the Author of the
 * Software shall not be used in advertising or otherwise to promote the
 * sale, use or other dealings in this Software without prior written
 * authorization from the Author.
 */

#include "config.h"

#include <gtk/gtk.h>

#include <libgimp/gimp.h>

#include <complex.h>
#include <fftw3.h>

#include "main.h"
#include "render.h"

#include "plugin-intl.h"

#include <string.h>
#include <stdlib.h>
#include <math.h>

/* We don't know the width at compile time, but most of the arrays we use
   will be the same size as the image. */
#define PIXEL(arr, x, y) (arr)[(y)*width + (x)]

double *partial_x = 0;
double *partial_y = 0;

/*  Public functions  */

/* Reduces noise in the image by bilateral filtering: Take a weighted average 
   of nearby pixels, giving them less weight if they're farther away in space
   _or_ in color. This helps preserve edges between differently-colored regions
   while smoothing similarly-colored regions.
*/
void noise_filter(double *image, double *newimage, int width, int height){
  int i,j,x,y;
  for(y = 0; y < height; y++){
    for(x = 0; x < width; x++){
      double acc, totalweight;
      acc = 0;
      totalweight = 0;
      double pixel = PIXEL(image, x, y);
      for(i=-2; i<=2; i++){
	for(j=-2; j<=2; j++){
	  double weight, e;
	  if(y+i < 0 || x+j < 0 || y+i >= height || x+j >= width)
	    continue;
	  e = (i*i + j*j)/25; // account for distance...
	  e += pow(PIXEL(image, x+j, y+i)-pixel, 2)/2500; // and color
	  weight = exp(-e);
	  totalweight += weight;
	  acc += (PIXEL(image, x+j, y+i)) * weight;
	}
      }
      PIXEL(newimage, x, y) = acc/totalweight;
    }
  }
}

/* Calculates the partial derivatives w.r.t. x and y at every point (which
   can be combined to get the gradient)
   Assumes globals partial_x and partial_y arrays are allocated
*/
void calculate_gradient(double *image, 
			double *xresult, double *yresult,
			int width, int height){
  int x, y;
  double px, py;

  for(y = 0; y < height; y++){
    for(x = 0; x < width; x++){
      if(x == width){
	/* At the edges, just use difference to adjacent pixel */
	px = PIXEL(image, x, y) - PIXEL(image, x-1, y);
      }
      else if (x == 0){
	px = PIXEL(image, x+1, y) - PIXEL(image, x, y);
      }
      else{
	/* Everywhere else, average difference to both adjacent pixels.
	   The center pixel ends up not affecting this calculation. */
	px = (PIXEL(image, x+1, y) - PIXEL(image, x-1, y))/2;
      }
      /* ...then do the same in the y direction. */
      if(y == height){
	py = PIXEL(image, x, y) - PIXEL(image, x, y-1);
      }
      else if (y == 0){
	py = PIXEL(image, x, y+1) - PIXEL(image, x, y);
      }
      else{
	py = (PIXEL(image, x, y+1) - PIXEL(image, x, y+1))/2;
      }
      PIXEL(xresult, x, y) = px;
      PIXEL(yresult, x, y) = py;
    }
  }
}

/* This is a filter that accentuates edges. It doesn't produce a very realistic
   image, so it's only used as an intermediate step.*/
void shock_filter(double *image, double *newimage, int width, int height){
  int x, y;

  calculate_gradient(image, partial_x, partial_y, width, height);

  for(y = 0; y < height; y++){
    for(x = 0; x < width; x++){
      double px; 
      double py;
      double L;

      px = PIXEL(partial_x, x, y);
      py = PIXEL(partial_y, x, y);
      PIXEL(newimage, x, y) = sqrt(px*px + py*py)/2;
      if(x > 1 && y > 1){
	// Compute Laplacian
	L = PIXEL(partial_x, x, y) - PIXEL(partial_x, x-1, y) +
	  PIXEL(partial_y, x, y) - PIXEL(partial_y, x, y-2);
	if(L > 0){
	  PIXEL(newimage, x, y) *= -1;
	}
	PIXEL(newimage, x, y) += PIXEL(image, x, y);
	if(PIXEL(newimage, x, y) > 255)
	  PIXEL(newimage, x, y) = 255;
	if(PIXEL(newimage, x, y) < 0)
	  PIXEL(newimage, x, y) = 0;
      }
    }
  }
  for(y = 0; y < height; y++){
    for(x = 0; x < 2; x++){
      PIXEL(newimage,x,y) = PIXEL(image,x,y);
    }
  }
  for(x = 0; x < width; x++){
    for(y = 0; y < 2; y++){
      PIXEL(newimage,x,y) = PIXEL(image,x,y);
    }
  }
}

/* This is the step the paper calls "truncating gradients". The idea
   is that you only want to pay attention to the strongest edges, 
   so you zero out the rest.
*/
void flatten(double *image, double *newimage, int width, int height){
  int histogram_x[256];
  int histogram_y[256];
  int histogram_xy[256];
  int histogram_x_minus_y[256];
  int x, y;
  calculate_gradient(image, partial_x, partial_y, width, height);
  for(x=0; x < width; x++){
    for(y=0; y < height; y++){
      int px = PIXEL(partial_x, x, y);
      int py = PIXEL(partial_y, x, y);
      if(abs(px) > abs(py)/2 && abs(py) > abs(px)/2){
	if(px * py > 0){
	  histogram_xy[(px + py)/2]++;
	}
	else{
	  histogram_x_minus_y[(py+py)/2]++;
	}
      }
      else if(abs(px) > abs(py)){
	histogram_x[px]++;
      }
      else{
	histogram_y[py]++;
      }
    }
  }
  // TODO: set gradients below a threshold to 0
}

/*
void find_kernel(double *image, double *predicted, double *kernel,
		 double *fft_in, double *fft_out,
		 fftw_plan plan, fftw_plan reverse_plan;
		 int width, int height)
{
  //TODO: compute various partials

  // TODO: precompute (complex conjugate of fft(predicted)) * fft(predicted)
  // TODO: look up weights for various partials.

  //TODO: find optimal K using conjugate gradient method


}
*/
void
render (gint32              image_ID,
	GimpDrawable       *drawable,
	PlugInVals         *vals,
	PlugInImageVals    *image_vals,
	PlugInDrawableVals *drawable_vals)
{
  double *scratch_image_1, *scratch_image_2;
  fftw_complex *fft_out, *fft_out2;
  fftw_plan plan, reverse_plan;
  int x1, y1, x2, y2, width, height;
  int i, j, k, channel;
  int num_channels, kernel_channels;
  guchar *rect;
  gint32 file;
  GimpDrawable *drawable_k;

  GimpPixelRgn rgn_in, rgn_out;

  /* Get selection boundary */
  gimp_drawable_mask_bounds (drawable->drawable_id,
			     &x1, &y1,
			     &x2, &y2);

  width = x2 - x1;
  height = y2 - y1;

  num_channels = gimp_drawable_bpp(drawable->drawable_id);

  rect = g_malloc(sizeof(guchar)* width * height * num_channels);

  scratch_image_1 = fftw_malloc(sizeof(double) * width * height);
  scratch_image_2 = fftw_malloc(sizeof(double) * width * height);
  fft_out = fftw_malloc(sizeof(fftw_complex) * width * (height/2 + 1));
  fft_out2 = fftw_malloc(sizeof(fftw_complex) * width * (height/2 + 1));

  plan = fftw_plan_dft_r2c_2d(height, width, scratch_image_1, fft_out, FFTW_ESTIMATE);
  reverse_plan = fftw_plan_dft_c2r_2d(height, width, fft_out, scratch_image_2, FFTW_ESTIMATE);
  partial_x = malloc(sizeof(double) * width * height);
  partial_y = malloc(sizeof(double) * width * height);

  file = gimp_file_load(GIMP_RUN_NONINTERACTIVE, 
			"/home/rehana/projects/deblur/kernel.xcf",
			"/home/rehana/projects/deblur/kernel.xcf");

  drawable_k = gimp_drawable_get(gimp_image_get_active_drawable(file));
  gimp_pixel_rgn_init(&rgn_in, drawable_k,
		      0, 0, 0, 0, FALSE, FALSE);

  kernel_channels = gimp_drawable_bpp(drawable_k->drawable_id);
  gimp_pixel_rgn_get_rect(&rgn_in, rect, 0, 0, 22, 16);
  memset(scratch_image_1, 0, sizeof(double) * width * height);

  /*  k = 0;
  for(i = 0; i < 22*16; i++){
    k += rect[i];
  }
  g_message("%d", k);
  for(i = 0; i < 16; i++){
    for(j=0; j<22; j++){
      scratch_image_1[i * width + j] = ((double)rect[kernel_channels*(i * 22 + j)])/8000;
    }
  }

    fftw_execute(plan);

  memcpy(fft_out2, fft_out, sizeof(fftw_complex) * width * (height/2 + 1));
  */
  gimp_pixel_rgn_init(&rgn_in, drawable,
		      x1, y1,
		      width, height, 
		      FALSE, FALSE);
  gimp_pixel_rgn_init(&rgn_out, drawable, 
		      x1, y1, 
		      width, height,
		      TRUE, TRUE);

  // gimp_tile_cache_ntiles(width / 64);
 
  gimp_pixel_rgn_get_rect(&rgn_in, rect, 
			  x1, y1, 
			  width, height);

  for(channel = 0; channel < num_channels; channel++){
    for(k = 0; k < width * height; k++){
      scratch_image_1[k] = (double)rect[k*num_channels+channel];
    }
    /* TODO: noise filter will be better using combined color distance, not
     separate channels */
    noise_filter(scratch_image_1, scratch_image_2, width, height);
    for(i = 0; i < 5; i++){
      memcpy(scratch_image_1, scratch_image_2, sizeof(double) * width * height);
      shock_filter(scratch_image_1, scratch_image_2, width, height);
    }
    for(k=0; k < width * height; k++){
      rect[k*num_channels+channel] = (guchar)(scratch_image_2[k]);
    }
  } 

  /* Deconvolution
  for(channel = 0; channel < num_channels; channel++){
    for(k = 0; k < width * height; k++){
      scratch_image_1[k] = (double)rect[k*num_channels+channel];
    }
    
    fftw_execute(plan);
    for(k = 0; k< width * (height/2+1); k++){
      fft_out[k] = fft_out[k] / fft_out2[k];
    }
    fftw_execute(reverse_plan);
    
    for(k=0; k < width * height; k++){
      rect[k*num_channels+channel] = (guchar)(scratch_image_2[k]/(width * height));
    }
    }*/
 

  /* (Debug)
  for(k = 0; k < width * height; k++){
    for(i = 0; i < num_channels; i++){
      rect[k* num_channels+i] = scratch_image_1[k]*8000;
      if(scratch_image_1[k] > 0)
	g_message("%d", k);
    }
    }*/

  gimp_pixel_rgn_set_rect(&rgn_out, rect, 
			  x1, y1, 
			  width, height);

  /*  Update the modified region */
  gimp_drawable_flush (drawable);
  gimp_drawable_merge_shadow (drawable->drawable_id, TRUE);
  gimp_drawable_update (drawable->drawable_id,
			x1, y1,
			width, height);

  g_free(rect);
  g_message("Done");

  fftw_destroy_plan(plan);
  fftw_destroy_plan(reverse_plan);
  fftw_free(scratch_image_1);
  fftw_free(scratch_image_2);
  fftw_free(fft_out);
  fftw_free(fft_out2);

  free(partial_x);
  free(partial_y);
}
