#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "utils.h"
#include "gnuplotScript.h"

/* writes a gnuplot script for the file based on command line options */
void genPlotScript(char *fileName, optStruct *opts, plotStruct *gp)
{
  int i, j, k;
  char fileHead[100], outputFile[100], gnuPlot[200], *head;
  
  // create file name to generate gnuplot script to plot data in this file
  strcpy(fileHead, fileName);   // store the filename in another file for safety
  head = strtok(fileHead, "."); // get the core name of the file
  createOutputFileName("plot", head, "gp", gnuPlot);

  // open file
  FILE *gpFile = fopen(gnuPlot, "w");

  // for 1D flamelet
  if (gp->plotDim == 1)
    {
      // plot all the species
      fprintf(gpFile, "set terminal %s 1%s\n",opts->termType,opts->postTerm);
      if (opts->outTerm)
        {
          createOutputFileName("Y", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }
      fprintf(gpFile, "set title \"Mass Fraction Profiles\"\n");
      fprintf(gpFile, "set ylabel \"Y\"\n");
      fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
      fprintf(gpFile, "set label 1 \"%s\" at screen 0.05,0.96\n",gp->timeLabel);
      if (opts->Yrange)
        fprintf(gpFile, "set yrange [%f:%f]\n", opts->Ymin[0],opts->Ymax[0]);
      else
        fprintf(gpFile, "set yrange [0:*]\n");
      fprintf(gpFile, "set key out vert center right\n");
      fprintf(gpFile, "\n");
      
      fprintf(gpFile, "plot '%s' using %i:%i %s title \"%s\"", 
              fileName, opts->xIndex, gp->iY+gp->plotIndex[0], opts->rcIndex,
              gp->plotSpecies[0]);
      if (opts->nPlotSpecies > 1)
        fprintf(gpFile, ", \\\n");
      for (i=1; i < opts->nPlotSpecies; i++)
        {
          fprintf(gpFile, "      '' using %i:%i %s title \"%s\"", 
                  opts->xIndex, gp->iY+gp->plotIndex[i], opts->rcIndex,
                  gp->plotSpecies[i]);
          if (i == opts->nPlotSpecies - 1)
            fprintf(gpFile, "\n");
          else
            fprintf(gpFile, ", \\\n");
        }
      fprintf(gpFile, "\n");
      
      // plot scalar dissipation rate
      if (opts->plotChi && gp->chiPlot[1])
        {
          fprintf(gpFile, "\n");
          fprintf(gpFile, "set terminal %s 2%s\n", opts->termType,
                  opts->postTerm);
          if (opts->outTerm)
            {
              createOutputFileName("chi", head, opts->termType, outputFile);
              fprintf(gpFile, "set output \"%s\"\n", outputFile);
            }
          fprintf(gpFile, "set title \"Scalar Dissipation\"\n");
          fprintf(gpFile, "set ylabel \"chi\"\n");
          if (opts->chiRange[1])
            fprintf(gpFile, "set yrange [%f:%f]\n",opts->chiMin[1],opts->chiMax[1]);
          else
            fprintf(gpFile, "set autoscale\n");
          fprintf(gpFile, "plot '%s' using %i:%i %s notitle\n", fileName,
                  opts->xIndex, gp->iChi[1], opts->rcIndex);
        }
      fprintf(gpFile, "\n");

      // plot temperature (and enthalpy)
      fprintf(gpFile, "set terminal %s 3%s\n", opts->termType, opts->postTerm);
      if (opts->outTerm)
        {
          createOutputFileName("T", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }
      fprintf(gpFile, "set title \"Temperature and Enthalpy Profiles\"\n");
      fprintf(gpFile, "set ylabel \"T [K]\"\n");
      fprintf(gpFile, "set autoscale\n");
      if (opts->Trange)
        fprintf(gpFile, "set yrange [%f:%f]\n", opts->Tmin, opts->Tmax);
      if (gp->HPlot == 1 && opts->plotH == 1)
        {
          fprintf(gpFile, "set y2label \"H [kJ/kg]\"\n");
          fprintf(gpFile, "set ytics nomirror\n");
          fprintf(gpFile, "set y2tics\n");
          fprintf(gpFile, "set key inside vert bottom center\n");
        }
      fprintf(gpFile, "\n");
      fprintf(gpFile,"plot '%s' using %i:%i %s", fileName, opts->xIndex,
              gp->iT, opts->rcIndex);
      if (gp->HPlot == 1 && opts->plotH == 1)
        {
          fprintf(gpFile, " title \"Temperature\" axis x1y1, \\\n");
          fprintf(gpFile, " '' using %i:($%i*1e-3) %s title \"Total enthalpy\" axis x1y2\n", opts->xIndex, gp->iHTOT, opts->rcIndex);
          fprintf(gpFile, "\nreset y2\n");
        }
      else
        fprintf(gpFile, " notitle\n");
      fprintf(gpFile, "\n");
      
      // plot all miscellaneous variables
      for (i=0; i < gp->nMisc; i++)
        {
          fprintf(gpFile, "\n");
          fprintf(gpFile, "set terminal %s %i%s\n", opts->termType,
                  4+i,opts->postTerm);
          if (opts->outTerm)
            {
              createOutputFileName(gp->miscNames[i],head,opts->termType,
                                   outputFile);
              fprintf(gpFile, "set output \"%s\"\n", outputFile);
            }
          fprintf(gpFile, "\n");
          fprintf(gpFile, "set title \"%s vs %s\"\n", gp->miscNames[i],
                  gp->xLabel);
          fprintf(gpFile, "set ylabel \"%s\"\n", gp->miscNames[i]);
          fprintf(gpFile, "set autoscale\n");
          fprintf(gpFile, "set nokey\n");
          fprintf(gpFile, "plot '%s' using 1:%i %s\n", fileName, gp->iMisc+i,
                  opts->rcIndex);
        }
    } // end 1D plot script
  else if (gp->plotDim == 2)
    {
      fprintf(gpFile, "set style line 1 lt 7 lc -1\n");
      fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
      fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
      fprintf(gpFile, "set xrange [%f:%f]\n", opts->x1min, opts->x1max);
      fprintf(gpFile, "set yrange [%f:%f]\n", opts->x2min, opts->x2max);
     
      fprintf(gpFile, "set label 1 \"%s\" at screen 0.05,0.96\n",gp->timeLabel);

      if ( (opts->plotCol < -1) && (opts->plotRow < -1) )
        {
          fprintf(gpFile, "set pm3d map");
          if (opts->smooth > 0)
            fprintf(gpFile, " interpolate %i,%i", opts->smooth, opts->smooth);
          fprintf(gpFile, "\n");
          fprintf(gpFile, "set palette rgbformulae 22,13,-31\n");
        }
      fprintf(gpFile, "\n");

      if (opts->outTerm)
        fprintf(gpFile,"set terminal %s %s\n",opts->termType,opts->postTerm);

      // write plot scripts for all species      
      for (i=0; i < opts->nPlotSpecies; i++)
        {
          if (opts->outTerm)
            {
              createOutputFileName(gp->plotSpecies[i], head, opts->termType,
                                   outputFile);
              fprintf(gpFile, "set output \"%s\"\n", outputFile);
            }
          else
            fprintf(gpFile, "set terminal %s %i%s\n", opts->termType, i+1,
                    opts->postTerm);

          fprintf(gpFile, "\n");
          fprintf(gpFile,"set title \"Mass fraction %s\"\n",gp->plotSpecies[i]);
          if (opts->Yrange)
            fprintf(gpFile,"set cbrange [%f:%f]\n",opts->Ymin[i],opts->Ymax[i]);
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", 
                  fileName, gp->iY+gp->plotIndex[i]);
          fprintf(gpFile, "\n");
        }
      fprintf(gpFile, "\n");

      // write plot script for temperature
      if (opts->outTerm)
        {
          createOutputFileName("T", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }
      else
        fprintf(gpFile, "set terminal %s %i%s\n", opts->termType, 
                opts->nPlotSpecies+1, opts->postTerm);

      fprintf(gpFile, "\n");
      fprintf(gpFile, "set title \"Temperature Field\"\n");
      fprintf(gpFile, "set zlabel \"T [K]\"\n");
      if (opts->Trange)
        fprintf(gpFile, "set cbrange [%f:%f]\n", opts->Tmin,opts->Tmax);
      else
        fprintf(gpFile, "set autoscale\n");

      fprintf(gpFile, "set xrange [%f:%f]\n", opts->x1min, opts->x1max);
      fprintf(gpFile, "set yrange [%f:%f]\n", opts->x2min, opts->x2max);
      fprintf(gpFile, "splot '%s' using 1:2:%i notitle", fileName, gp->iT);
      fprintf(gpFile, "\n");
      
      if (gp->HPlot && opts->plotH)
        {
          if (opts->outTerm)
            {
              createOutputFileName("H", head, opts->termType, outputFile);
              fprintf(gpFile, "set output \"%s\"\n", outputFile);
            }
          else
            fprintf(gpFile, "set terminal %s %i%s\n", opts->termType, 
                    opts->nPlotSpecies+2, opts->postTerm);

          fprintf(gpFile, "set title \"Total Enthalpy Field\"\n");
          fprintf(gpFile, "set autoscale\n");
          fprintf(gpFile, "splot '%s' using 1:2:($%i*1e-3) notitle\n", fileName,
                  gp->iHTOT);
        }
      fprintf(gpFile, "\n");
      
      // write plot script for scalar dissipation rates
      // chi1
      if (gp->chiPlot[1] && opts->plotChi)
        {
          if (opts->outTerm)
            {
              createOutputFileName("chi1", head, opts->termType, outputFile);
              fprintf(gpFile, "set output \"%s\"\n", outputFile);
            }
          else
            fprintf(gpFile, "set terminal %s %i%s\n", opts->termType, 
                    opts->nPlotSpecies+3, opts->postTerm);

          fprintf(gpFile, "set title \"Scalar dissipation rate for Z1\"\n");
          if (opts->chiRange[1])
            fprintf(gpFile,"set cbrange [%f:%f]\n",opts->chiMin[1],opts->chiMax[1]);
          else
            fprintf(gpFile, "set autoscale\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName,
                  gp->iChi[1]);
        }
      fprintf(gpFile, "\n");

      // chi12
      if (gp->chiPlot[3] && opts->plotChi)
        {
          if (opts->outTerm)
            {
              createOutputFileName("chi12", head, opts->termType, outputFile);
              fprintf(gpFile, "set output \"%s\"\n", outputFile);
            }
          else
            fprintf(gpFile, "set terminal %s %i%s\n", opts->termType, 
                    opts->nPlotSpecies+4, opts->postTerm);

          fprintf(gpFile, "set title \"Scalar dissipation rate for Z1\"\n");
          if (opts->chiRange[3])
            fprintf(gpFile,"set cbrange [%f:%f]\n",opts->chiMin[3],opts->chiMax[3]);
          else
            fprintf(gpFile, "set autoscale\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName,
                  gp->iChi[3]);
        }
      fprintf(gpFile, "\n");
      
      // chi2
      if (gp->chiPlot[2] && opts->plotChi)
        {
          if (opts->outTerm)
            {
              createOutputFileName("chi2", head, opts->termType, outputFile);
              fprintf(gpFile, "set output \"%s\"\n", outputFile);
            }
          else
            fprintf(gpFile, "set terminal %s %i%s\n", opts->termType, 
                    opts->nPlotSpecies+5, opts->postTerm);

          fprintf(gpFile, "set title \"Scalar dissipation rate for Z2\"\n");
          if (opts->chiRange[2])
            fprintf(gpFile,"set cbrange [%f:%f]\n",opts->chiMin[2],opts->chiMax[2]);
          else
            fprintf(gpFile, "set autoscale\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName,
                  gp->iChi[2]);
        }
      fprintf(gpFile, "\n");
    } // end 2D plot script
  
  fclose(gpFile);
}
/* -------------------------------------------------------------------------- */

void genMultiPlotScript(char *fileName, char *cFile, optStruct *opts,
                        plotStruct *gp)
{
  int i, j, k, iTmp;
  int nYWindows, curWindow, Yremain, curY;
  double rangeFactor;
  char fileHead[100], outputFile[100], gnuPlot[200], speciesFile[50], *head;
  
  // create file name to generate gnuplot script to plot data in this file
  strcpy(fileHead, fileName);   // store the filename in another file for safety
  head = strtok(fileHead, "."); // get the core name of the file
  createOutputFileName("plot", head, "gp", gnuPlot);
  
  FILE *gpFile;
  gpFile = fopen(gnuPlot, "w");

  // set margins and size
  fprintf(gpFile, "set rmargin 0\n");
  fprintf(gpFile, "set lmargin 0\n");
  fprintf(gpFile, "set tmargin 0\n");
  fprintf(gpFile, "set bmargin 0\n");

  // for 1D flamelet
  if (gp->plotDim == 1)
    {
      if (opts->outTerm)
        {
          fprintf(gpFile, "set terminal %s %s\n",opts->termType,
                  opts->postTerm);
          createOutputFileName("multi1D", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }
      else
        fprintf(gpFile, "set terminal %s 1%s\n",opts->termType,
                opts->postTerm);

      fprintf(gpFile, "set size 0.37, 0.37\n");
      fprintf(gpFile, "set multiplot\n");

      // set labels for time
      fprintf(gpFile, "set label 1 \"%s\" at screen 0.2,0.95\n",gp->timeLabel);

      // plot temperature (and enthalpy) in bottom left corner
      fprintf(gpFile, "set origin 0.1,0.1\n");
      fprintf(gpFile, "set xlabel \"%s\" offset 0,0.5\n", gp->xLabel); 
      fprintf(gpFile, "set ylabel \"T [K]\" offset 0.5,0\n");
      if (opts->Trange)
        fprintf(gpFile, "set yrange [%f:%f]\n", opts->Tmin, opts->Tmax);
      else
        fprintf(gpFile, "set yrange [0:*]\n");

      if (gp->HPlot == 1 && opts->plotH == 1)
        {
          fprintf(gpFile, "set y2label \"H [kJ/kg]\"\n");
          fprintf(gpFile, "set ytics nomirror\n");
          fprintf(gpFile, "set y2tics\n");
          //fprintf(gpFile, "set key inside vert bottom center\n");
        }
      fprintf(gpFile, "\n");
      fprintf(gpFile,"plot '%s' using %i:%i %s %s", fileName, opts->xIndex,
              gp->iT, opts->rcIndex, opts->wLine);
      if (gp->HPlot == 1 && opts->plotH == 1)
        {
          fprintf(gpFile, " title \"Temperature\" axis x1y1, \\\n");
          fprintf(gpFile, " '' using %i:($%i*1e-3) %s title \"Total enthalpy\" axis x1y2\n", opts->xIndex, gp->iHTOT, opts->rcIndex);
          fprintf(gpFile, "\nunset y2tics\nunset y2label");
        }
      else
        fprintf(gpFile, " notitle");
      if (opts->compFile)
        {
          fprintf(gpFile, ", \\\n");
          fprintf(gpFile, " '%s' using %i:%i %s notitle\n", cFile,
                  opts->xIndex, gp->iT, opts->rcIndex);
        }
      else
        fprintf(gpFile, "\n");
      fprintf(gpFile, "\n");

      // plot scalar dissipation rate in bottom left corner
      fprintf(gpFile, "set origin 0.6,0.1\n");
      if (gp->chiPlot[2])
        {
          fprintf(gpFile, "set xlabel \"%s\" offset 0,0.5\n", gp->xLabel);
          fprintf(gpFile, "set ylabel \"chi2\" offset 0.5,0\n");
          if (opts->chiRange[2])
            fprintf(gpFile, "set yrange [%f:%f]\n",opts->chiMin[2],opts->chiMax[2]);
          else
            fprintf(gpFile, "set autoscale\n");
          fprintf(gpFile, "plot '%s' using %i:%i %s %s notitle\n", fileName,
                  opts->xIndex, gp->iChi[2], opts->rcIndex, opts->wLine);
        }

      // plot all the species in the top left corner
      fprintf(gpFile, "set origin 0.1,0.55\n");
      fprintf(gpFile, "unset xlabel\n");
      fprintf(gpFile, "set ylabel \"Y\" offset 0.5,0\n");
      if (opts->Yrange)
        fprintf(gpFile, "set yrange [%f:%f]\n", opts->Ymin[0],opts->Ymax[0]);
      else
        fprintf(gpFile, "set autoscale\n");
      fprintf(gpFile, "set key out vert center right\n");
      fprintf(gpFile, "\n");
      
      fprintf(gpFile, "plot '%s' using %i:%i %s %s title \"%s\"", 
              fileName, opts->xIndex, gp->iY+gp->plotIndex[0], opts->rcIndex,
              opts->wLine, gp->plotSpecies[0]);
      if (opts->nPlotSpecies > 1)
        fprintf(gpFile, ", \\\n");
      for (i=1; i < opts->nPlotSpecies; i++)
        {
          fprintf(gpFile, "      '' using %i:%i %s %s title \"%s\"", 
                  opts->xIndex, gp->iY+gp->plotIndex[i], opts->rcIndex,
                  opts->wLine, gp->plotSpecies[i]);
          if (i < opts->nPlotSpecies - 1)
            fprintf(gpFile, ", \\\n");
          /* else */
/*             fprintf(gpFile, ", \\\n"); */
        }
      if (opts->compFile)
        {
          fprintf(gpFile, ", \\\n");
          fprintf(gpFile, " '%s' using %i:%i %s notitle", 
                  cFile, opts->xIndex, gp->iY+gp->plotIndex[0], opts->rcIndex);
          if (opts->nPlotSpecies > 1)
            fprintf(gpFile, ", \\\n");
          for (i=1; i < opts->nPlotSpecies; i++)
            {
              fprintf(gpFile, "      '' using %i:%i %s title \"%s\"", 
                      opts->xIndex, gp->iY+gp->plotIndex[i], opts->rcIndex,
                      gp->plotSpecies[i]);
              if (i == opts->nPlotSpecies - 1)
                fprintf(gpFile, "\n");
              else
                fprintf(gpFile, ", \\\n");
            }
        }
      else
        fprintf(gpFile, "\n");
      fprintf(gpFile, "\n");

      // plot scalar dissipation rate in top right corner
      fprintf(gpFile, "set origin 0.6,0.55\n");
      if (gp->chiPlot[1])
        {
          fprintf(gpFile, "set xlabel \"%s\" offset 0,0.5\n", gp->xLabel);
          fprintf(gpFile, "set ylabel \"chi\" offset 0.5,0\n");
          if (opts->chiRange[1])
            fprintf(gpFile, "set yrange [%f:%f]\n",opts->chiMin[1],opts->chiMax[1]);
          else
            fprintf(gpFile, "set autoscale\n");
          fprintf(gpFile, "plot '%s' using %i:%i %s %s notitle", fileName,
                  opts->xIndex, gp->iChi[1], opts->rcIndex, opts->wLine);
        }
      if (opts->compFile)
        {
          fprintf(gpFile, ", \\\n");
          fprintf(gpFile, " '%s' using %i:%i %s w l notitle\n", cFile,
                  opts->xIndex, gp->iChi[1], opts->rcIndex);
        }
      else
        fprintf(gpFile, "\n");
      fprintf(gpFile, "unset multiplot\n");

      /* ----- plot any miscellaneous variables in new window ----- */
      if (opts->plotMisc)
        {
          fprintf(gpFile, "set autoscale\n");
          if (gp->nMisc%4 == 0)
            nYWindows = gp->nMisc/4;
          else
            nYWindows = gp->nMisc/4 + 1;

          Yremain = gp->nMisc; // number of total species yet to be plotted
          curY=0;                       // index of current species
          // plot 4 species per window, starting from bottom left
          for (k=0; k < nYWindows; k++)
            {
              // set output parameters
              if (opts->outTerm)
                {
                  fprintf(gpFile, "set terminal %s %s\n",
                          opts->termType,opts->postTerm);
                  sprintf(speciesFile, "multi1D_Misc%i", k+1);
                  createOutputFileName(speciesFile, head, opts->termType,
                                       outputFile);
                  fprintf(gpFile, "set output \"%s\"\n", outputFile);
                }
              else
                fprintf(gpFile,"set terminal %s %i%s\n",opts->termType,k+2,
                        opts->postTerm);
              fprintf(gpFile, "set multiplot\n");

              for (j=0; j<4; j++, curY++, Yremain--)
                { 
                  // if no species left for this window, leave
                  if (Yremain == 0)
                    break;

                  if (j == 0) // bottom left
                    {
                      fprintf(gpFile, "set origin 0.1,0.1\n");
                      fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);

                    }
                  else if (j == 1) // bottom right
                    {
                      fprintf(gpFile, "set origin 0.6,0.1\n");
                      fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
                      fprintf(gpFile, "unset ylabel\n");
                    }
                  else if (j == 2) // top left
                    {
                      fprintf(gpFile, "set origin 0.1,0.55\n");
                      fprintf(gpFile, "unset xlabel\n");
                      fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
                    }
                  else if (j == 3) // top right
                    {
                      fprintf(gpFile, "set origin 0.6,0.55\n");
                      fprintf(gpFile, "unset xlabel\n");
                      fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
                    }
                  fprintf(gpFile, "\n");

                  fprintf(gpFile, "set ylabel \"%s\"\n", gp->miscNames[curY]);
                  fprintf(gpFile, "plot '%s' using %i:%i %s notitle\n",
                          fileName,opts->xIndex,gp->iMisc+curY,opts->rcIndex);
                  fprintf(gpFile, "\n");
                } // end loop over species for this window

              fprintf(gpFile, "unset multiplot\n");
            } // end loop over windows
          fprintf(gpFile, "\n");
        }
      
      fprintf(gpFile, "unset multiplot\n");
    } // end 1D plot script
  else if (gp->plotDim == 2)
    {
      // set output parameters
      if (opts->outTerm)
        {
          fprintf(gpFile, "set terminal %s %s\n",opts->termType,
                  opts->postTerm);
          createOutputFileName("multi2D", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }
      else
        fprintf(gpFile,"set terminal %s 1%s\n",opts->termType,opts->postTerm);

      fprintf(gpFile, "set size 0.48, 0.48\n");
      fprintf(gpFile, "set multiplot\n");

      // set labels for time
      fprintf(gpFile, "set label 1\"%s\" at screen 0.8,0.97\n",gp->timeLabel);

      fprintf(gpFile, "set style line 1 lt 7 lc -1\n");
      fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
      fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
      fprintf(gpFile, "set xrange [%f:%f]\n",opts->x1min, opts->x1max);
      fprintf(gpFile, "set yrange [%f:%f]\n",opts->x2min, opts->x2max);
      fprintf(gpFile, "set xtics 0.1\n");
      fprintf(gpFile, "set ytics 0.1\n");
      
      if ( (opts->plotCol < -1) && (opts->plotRow < -1) )
        {
          fprintf(gpFile, "set pm3d map");
          if (opts->smooth > 0)
            fprintf(gpFile, " interpolate %i,%i", opts->smooth, opts->smooth);
          fprintf(gpFile, "\n");
          fprintf(gpFile, "set palette rgbformulae 22,13,-31\n");
        }
      fprintf(gpFile, "\n");

      /* ----- use standard output for first window ----- */
      // write plot script for temperature
      fprintf(gpFile, "set origin 0.02,0.05\n");
      fprintf(gpFile, "\n");
      fprintf(gpFile, "set title \"Temperature Field\"\n");
      if (opts->Trange)
        fprintf(gpFile, "set cbrange [%f:%f]\n", opts->Tmin,opts->Tmax);
      else
        fprintf(gpFile, "set cbrange [850:*]\n");

      fprintf(gpFile, "splot '%s' using 1:2:%i notitle, 'stoichContour.dat' u 1:2:3 with lines ls 1 notitle\n", 
              fileName, gp->iT);
      fprintf(gpFile, "\n");

      // write plot script for enthalpy
      fprintf(gpFile, "set origin 0.5,0.05\n");
      fprintf(gpFile, "\n");
      if (opts->plotH)
        {
          fprintf(gpFile, "set title \"Total Enthalpy\"\n");
          fprintf(gpFile, "unset ylabel\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle, 'stoichContour.dat' u 1:2:3 with lines ls 1 notitle\n",
                  fileName, gp->iHTOT);
          fprintf(gpFile, "\n");
        }

      // write plot script for progress variable
      fprintf(gpFile, "set origin 0.02,0.5\n");
      fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
      fprintf(gpFile, "unset xlabel\n");
      fprintf(gpFile, "\n");
      if (opts->plotProg)
        {
          fprintf(gpFile, "set title \"Progress Variable\"\n");
          fprintf(gpFile, "set cbrange [0:1]\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle, 'stoichContour.dat' u 1:2:3 with lines ls 1 notitle\n", fileName, gp->iProg);
        }
      fprintf(gpFile, "\n");
      fprintf(gpFile, "unset multiplot\n");

      /* ----- plot scalar dissipation rates ----- */
      // set output parameters
      if (opts->outTerm)
        {
          fprintf(gpFile, "set terminal %s %s\n",opts->termType,
                  opts->postTerm);
          createOutputFileName("chi2D_", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }
      else
        fprintf(gpFile,"set terminal %s 2%s\n",opts->termType,opts->postTerm);

      fprintf(gpFile, "set multiplot\n");
      // write scripts for scalar dissipation rates along the bottom
      fprintf(gpFile, "set origin 0.02,0.05\n");
      fprintf(gpFile, "set title \"chiA\"\n");
      if (opts->chiRange[1])
        fprintf(gpFile, "set cbrange [%f:%f]\n", opts->chiMin[1],opts->chiMax[1]);
      else
        fprintf(gpFile, "set cbrange [0:*]\n");
      fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName, gp->iChi[1]);
      fprintf(gpFile, "\n");

      // chi12
      fprintf(gpFile, "set origin 0.5,0.05\n");
      fprintf(gpFile, "set title \"chiB\"\n");
      fprintf(gpFile, "unset ylabel\n");
      if (opts->chiRange[3])
        fprintf(gpFile, "set cbrange [%f:%f]\n", opts->chiMin[3],opts->chiMax[3]);
      else
        fprintf(gpFile, "set cbrange [*:*]\n");
      fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName, gp->iChi[3]);
      fprintf(gpFile, "\n");

      // chi2
      fprintf(gpFile, "set origin 0.02,0.5\n");
      fprintf(gpFile, "set title \"chiC\"\n");
      fprintf(gpFile, "unset ylabel\n");
      if (opts->chiRange[2])
        fprintf(gpFile, "set cbrange [%f:%f]\n", opts->chiMin[2],opts->chiMax[2]);
      else
        fprintf(gpFile, "set cbrange [0:*]\n");
      fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName, gp->iChi[2]);
      fprintf(gpFile, "\n");

      fprintf(gpFile, "unset multiplot\n");

      /* ----- if species were specified, plot 4 in a window ----- */
      if (opts->nPlotSpecies%4 == 0)
        nYWindows = opts->nPlotSpecies/4;
      else
        nYWindows = opts->nPlotSpecies/4 + 1;

      Yremain = opts->nPlotSpecies; // number of total species yet to be plotted
      curY=0;                       // index of current species
      // plot 4 species per window, starting from bottom left
      for (k=0; k < nYWindows; k++)
        {
          // set output parameters
          if (opts->outTerm)
            {
              fprintf(gpFile, "set terminal %s %s\n",
                      opts->termType,opts->postTerm);
              sprintf(speciesFile, "multi2D_Y%i", k+1);
              createOutputFileName(speciesFile, head, opts->termType,
                                   outputFile);
              fprintf(gpFile, "set output \"%s\"\n", outputFile);
            }
          else
            fprintf(gpFile,"set terminal %s %i%s\n",opts->termType,3+k,
                    opts->postTerm);
          fprintf(gpFile, "set multiplot\n");

          for (j=0; j<4; j++, curY++, Yremain--)
            { 
              // if no species left for this window, leave
              if (Yremain == 0)
                break;

              if (j == 0) // bottom left
                {
                  fprintf(gpFile, "set origin 0.02,0.05\n");
                  fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
                  fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
                }
              else if (j == 1) // bottom right
                {
                  fprintf(gpFile, "set origin 0.5,0.05\n");
                  fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
                  fprintf(gpFile, "unset ylabel\n");
                }
              else if (j == 2) // top left
                {
                  fprintf(gpFile, "set origin 0.02,0.5\n");
                  fprintf(gpFile, "unset xlabel\n");
                  fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
                }
              else if (j == 3) // top right
                {
                  fprintf(gpFile, "set origin 0.5,0.5\n");
                  fprintf(gpFile, "unset xlabel\n");
                  fprintf(gpFile, "unset ylabel\n");
                }
              fprintf(gpFile, "\n");

              fprintf(gpFile, "set title \"Mass fraction %s\"\n",
                      gp->plotSpecies[curY]);
              if (opts->Yrange)
                fprintf(gpFile, "set cbrange [%f:%f]\n",
                        opts->Ymin[curY],opts->Ymax[curY]);
              else
                fprintf(gpFile, "set cbrange [0:*]\n");

              fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n",
                      fileName, gp->iY+gp->plotIndex[curY]);
              fprintf(gpFile, "\n");
            } // end loop over species for this window

          fprintf(gpFile, "unset multiplot\n");
        } // end loop over windows
      fprintf(gpFile, "\n");

      /* ----- plot miscellaneous variables ----- */
      if (gp->nMisc%4 == 0)
        nYWindows = gp->nMisc/4;
      else
        nYWindows = gp->nMisc/4 + 1;

      Yremain = gp->nMisc; // number of variables yet to be plotted
      curY=0;                // index of current variable
      // plot 4 species per window, starting from bottom left
      for (k=0; k < nYWindows; k++)
        {
          // set output parameters
          if (opts->outTerm)
            {
              fprintf(gpFile, "set terminal %s %s\n",
                      opts->termType,opts->postTerm);
              sprintf(speciesFile, "multi2D_misc%i", k+1);
              createOutputFileName(speciesFile, head, opts->termType,
                                   outputFile);
              fprintf(gpFile, "set output \"%s\"\n", outputFile);
            }
          else
            fprintf(gpFile,"set terminal %s %i%s\n",opts->termType,
                    3+opts->nPlotSpecies/4+1+k, opts->postTerm);
          fprintf(gpFile, "set multiplot\n");

          for (j=0; j<4; j++, curY++, Yremain--)
            { 
              // if no species left for this window, leave
              if (Yremain == 0)
                break;

              if (j == 0) // bottom left
                {
                  fprintf(gpFile, "set origin 0.02,0.05\n");
                  fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
                  fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
                }
              else if (j == 1) // bottom right
                {
                  fprintf(gpFile, "set origin 0.5,0.05\n");
                  fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
                  fprintf(gpFile, "unset ylabel\n");
                }
              else if (j == 2) // top left
                {
                  fprintf(gpFile, "set origin 0.02,0.5\n");
                  fprintf(gpFile, "unset xlabel\n");
                  fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
                }
              else if (j == 3) // top right
                {
                  fprintf(gpFile, "set origin 0.5,0.5\n");
                  fprintf(gpFile, "unset xlabel\n");
                  fprintf(gpFile, "unset ylabel\n");
                }
              fprintf(gpFile, "\n");

              fprintf(gpFile, "set title \"%s\"\n", gp->miscNames[curY]);
              if (opts->Yrange)
                fprintf(gpFile, "set cbrange [%f:%f]\n",
                        opts->miscMin[curY],opts->miscMax[curY]);
              else
                fprintf(gpFile, "set cbrange [0:*]\n");

              fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n",
                      fileName, gp->iMisc+curY);
              fprintf(gpFile, "\n");
            } // end loop over species for this window

          fprintf(gpFile, "unset multiplot\n");
        } // end loop over windows
      fprintf(gpFile, "\n");
       
    } // end 2D plot script

  if (opts->nPlotSpecies == 0)
    fprintf(gpFile, "unset multiplot\n");
  fclose(gpFile);
}
/* -------------------------------------------------------------------------- */

void genSAE_Script2D(char *fileName, char *cFile, optStruct *opts,
                     plotStruct *gp)
{
  int i, j, k, iTmp;
  int nYWindows, curWindow, Yremain, curY;
  double rangeFactor;
  char fileHead[100], outputFile[100], gnuPlot[200], speciesFile[50], *head;
  
  // create file name to generate gnuplot script to plot data in this file
  strcpy(fileHead, fileName);   // store the filename in another file for safety
  head = strtok(fileHead, "."); // get the core name of the file
  createOutputFileName("plot", head, "gp", gnuPlot);
  
  FILE *gpFile;
  gpFile = fopen(gnuPlot, "w");

  // set margins and size
  fprintf(gpFile, "set rmargin 0\n");
  fprintf(gpFile, "set lmargin 0\n");
  fprintf(gpFile, "set tmargin 0\n");
  fprintf(gpFile, "set bmargin 0\n");

  
  // set output parameters
  if (opts->outTerm)
    {
      fprintf(gpFile, "set terminal %s %s\n",opts->termType,
              opts->postTerm);
      createOutputFileName("multi2D", head, opts->termType, outputFile);
      fprintf(gpFile, "set output \"%s\"\n", outputFile);
    }
  else
    fprintf(gpFile,"set terminal %s 1%s\n",opts->termType,opts->postTerm);

  fprintf(gpFile, "set size 0.48, 1.1\n");
  fprintf(gpFile, "set multiplot\n");

  // set labels for time
  fprintf(gpFile, "set label 1\"%s\" at screen 0.45,0.97\n",gp->timeLabel);

  fprintf(gpFile, "set xlabel \"%s\"\n", "{/Symbol x} = Z_1");
  fprintf(gpFile, "set ylabel \"%s\"\n", "{/Symbol h} = Z_2/(1-Z_1)");
  fprintf(gpFile, "set xrange [%f:%f]\n",0.0, 0.2);
  fprintf(gpFile, "set yrange [%f:%f]\n",0.0, 0.4);
  fprintf(gpFile, "set xtics 0.1\n");
  fprintf(gpFile, "set ytics 0.1\n");
      
  if ( (opts->plotCol < -1) && (opts->plotRow < -1) )
    {
      fprintf(gpFile, "set pm3d map");
      if (opts->smooth > 0)
        fprintf(gpFile, " interpolate %i,%i", opts->smooth, opts->smooth);
      fprintf(gpFile, "\n");
      fprintf(gpFile, "set palette rgbformulae 22,13,-31\n");
    }
  fprintf(gpFile, "\n");

  /* ----- use standard output for first window ----- */
  // write plot script for temperature
  fprintf(gpFile, "set origin 0.02,0.0\n");
  fprintf(gpFile, "\n");
  fprintf(gpFile, "set title \"Temperature Field\"\n");
  /* if (opts->Trange) */
/*     fprintf(gpFile, "set cbrange [%f:%f]\n", opts->Tmin,opts->Tmax); */
/*   else */
    fprintf(gpFile, "set cbrange [850:2600]\n");

  fprintf(gpFile, "splot '%s' using 1:2:%i notitle, 'stoichContour.dat' u 1:2:3 with lines lt -1 notitle\n", 
          fileName, gp->iT);
  fprintf(gpFile, "\n");

  // write plot script for a species
  fprintf(gpFile, "set origin 0.5,0.0\n");
  fprintf(gpFile, "unset ylabel\n");
  fprintf(gpFile, "\n");

  fprintf(gpFile, "set title \"Mass fraction %s\"\n",
          gp->plotSpecies[0]);
/*   if (opts->Yrange) */
/*     fprintf(gpFile, "set cbrange [%f:%f]\n", */
/*             opts->Ymin[0],opts->Ymax[0]); */
/*   else */
    fprintf(gpFile, "set cbrange [0:0.2]\n");

  fprintf(gpFile, "splot '%s' using 1:2:%i notitle, 'stoichContour.dat' u 1:2:3 with lines lt -1 notitle\n",
          fileName, gp->iY+gp->plotIndex[0]);
  fprintf(gpFile, "\n");

  fprintf(gpFile, "unset multiplot\n");
  fclose(gpFile);
}
/* -------------------------------------------------------------------------- */

void genFlowScript(char *fileName, optStruct *opts, plotStruct *gp)
{
  int i, j, k, iTmp;
  int nYWindows, curWindow, Yremain, curY;
  double rangeFactor;
  char fileHead[100], outputFile[100], gnuPlot[200], speciesFile[50], *head;
  char margFile[100], buf[100];
  
  // create file name to generate gnuplot script to plot data in this file
  strcpy(fileHead, fileName);   // store the filename in another file for safety
  head = strtok(fileHead, "."); // get the core name of the file
  createOutputFileName("plot", head, "gp", gnuPlot);
  strcpy(margFile, head);
  strcat(margFile, ".marg");
  
  FILE *gpFile;
  gpFile = fopen(gnuPlot, "w");

  // set margins and size
  fprintf(gpFile, "set rmargin 0\n");
  fprintf(gpFile, "set lmargin 0\n");
  fprintf(gpFile, "set tmargin 0\n");
  fprintf(gpFile, "set bmargin 0\n");

  // for 1D flamelet
  if (gp->plotDim == 1)
    {
      if (opts->outTerm)
        {
          fprintf(gpFile, "set terminal %s %s\n",opts->termType,
                  opts->postTerm);
          createOutputFileName("flow1D", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }
      else
        fprintf(gpFile, "set terminal %s 1%s\n",opts->termType,
                opts->postTerm);

      fprintf(gpFile, "set size 0.8, 0.37\n");
      fprintf(gpFile, "set multiplot\n");

      // set labels for time
      fprintf(gpFile, "set label 1 \"%s\" at screen 0.2,0.95\n",gp->timeLabel);

      // plot scalar global pdf at the bottom
      if (gp->pdfPlot)
        {
          fprintf(gpFile, "set origin 0.1,0.1\n");
          fprintf(gpFile, "set xlabel \"%s\" offset 0,0.5\n", gp->xLabel); 
          fprintf(gpFile, "set ylabel \"pdf\" offset 0.5,0\n");
          if (opts->logScale) fprintf(gpFile, "set log y\n");
          if (opts->pdfRange)
            fprintf(gpFile, "set yrange [%f:%f]\n", opts->pdfMin[0],
                    opts->pdfMax[1]);
          else
            fprintf(gpFile, "set yrange [0:*]\n");
          fprintf(gpFile, "plot '%s' using %i:%i %s notitle\n", fileName,
                  opts->xIndex, gp->iPDF, opts->rcIndex);
          if (opts->logScale) fprintf(gpFile, "unset log y\n");
        }
      fprintf(gpFile, "\n");

      // plot scalar dissipation rate on top
      fprintf(gpFile, "set origin 0.1,0.55\n");
      if (gp->chiPlot[1])
        {
          fprintf(gpFile, "unset xlabel\n");
          fprintf(gpFile, "set ylabel \"chi\" offset 0.5,0\n");
          if (opts->chiRange[1])
            fprintf(gpFile, "set yrange [%f:%f]\n",
                    opts->chiMin[1], opts->chiMax[1]);
          else
            fprintf(gpFile, "set autoscale\n");
          fprintf(gpFile, "plot '%s' using %i:%i %s notitle\n", fileName,
                  opts->xIndex, gp->iChi[1], opts->rcIndex);
        }
      fprintf(gpFile, "\n");

      fprintf(gpFile, "unset multiplot\n");
    } // end 1D plot script
  else if (gp->plotDim == 2)
    {
      fprintf(gpFile, "reset\n");
      // set output parameters
      if (opts->outTerm)
        {
          fprintf(gpFile,"set terminal %s %s\n",opts->termType,opts->postTerm);
          createOutputFileName("pdf", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }
      else
        fprintf(gpFile,"set terminal %s 1%s\n",opts->termType,opts->postTerm);

      // set labels
      fprintf(gpFile, "set multiplot\n");
      fprintf(gpFile, "set label 1 \"%s\" at screen 0.8,0.97\n",gp->timeLabel);
      fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
      fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
      
      /* ----- plot the pdf with its marginals  ----- */
      // marginals along the bottom
      fprintf(gpFile, "set size 0.45, 0.45\n");
      fprintf(gpFile, "set origin 0.02,0.02\n");
      fprintf(gpFile, "set xrange [%f:%f]\n", opts->x1min, opts->x1max);
      fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
      fprintf(gpFile, "set ylabel \"P|[%s]\"\n", gp->xLabel);
      if (opts->logScale) fprintf(gpFile, "set log y\n");
      if (opts->pdfRange)
        fprintf(gpFile, "set yrange [%f:%f]\n",opts->pdfMin[1],opts->pdfMax[1]);
      else
        fprintf(gpFile, "set yrange [0:*]\n");
      fprintf(gpFile, "plot '%s' using 1:2 index 0 notitle\n", margFile);
      fprintf(gpFile, "\n");

      fprintf(gpFile, "set origin 0.5,0.02\n");
      fprintf(gpFile, "set xlabel \"%s\"\n", gp->yLabel);
      fprintf(gpFile, "set ylabel \"P|[%s]\"\n", gp->yLabel);
      fprintf(gpFile, "set xrange [%f:%f]\n", opts->x2min, opts->x2max);
      if (opts->pdfRange)
        fprintf(gpFile, "set yrange [%f:%f]\n",opts->pdfMin[2],opts->pdfMax[2]);
      else
        fprintf(gpFile, "set yrange [0:*]\n");
      fprintf(gpFile, "plot '%s' using 1:2 index 1 notitle\n", margFile);
      fprintf(gpFile, "\n");
      if (opts->logScale) fprintf(gpFile, "unset log y\n");

      if ( opts->plotCol < -1 && opts->plotRow < -1 )
        {
          fprintf(gpFile, "set pm3d map");
          if (opts->smooth > 0)
            fprintf(gpFile, " interpolate %i,%i", opts->smooth, opts->smooth);
          fprintf(gpFile, "\n");
          fprintf(gpFile, "set palette rgbformulae 22,13,-31\n");
          fprintf(gpFile, "\n");
        }

      // joint-pdf in top left
      fprintf(gpFile, "set size 0.6, 0.6\n");
      fprintf(gpFile, "set origin 0.2,0.45\n");
      fprintf(gpFile, "set xrange [%f:%f]\n", opts->x1min, opts->x1max);
      fprintf(gpFile, "set yrange [%f:%f]\n", opts->x2min, opts->x2max);
      fprintf(gpFile, "set title \"P[%s,%s]\"\n", gp->xLabel, gp->yLabel);
      fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
      fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
      if (opts->logScale) fprintf(gpFile, "set log cb\n");
      if (opts->pdfRange || opts->logScale)
        fprintf(gpFile,"set cbrange [%e:%e]\n",opts->pdfMin[0],opts->pdfMax[0]);
      else
        fprintf(gpFile, "set cbrange [0:*]\n");
      fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName, gp->iPDF);
      fprintf(gpFile, "\n");

      if (opts->logScale) fprintf(gpFile, "unset log cb\n");
      fprintf(gpFile, "unset multiplot\n");

      /* ----- plot scalar dissipation rates ----- */
      if (opts->outTerm)
        {
          fprintf(gpFile, "set terminal %s %s\n",opts->termType,opts->postTerm);
          createOutputFileName("chi", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }
      else
        fprintf(gpFile,"set terminal %s 2%s\n", opts->termType, opts->postTerm);

      fprintf(gpFile, "set multiplot\n");
      fprintf(gpFile, "set size 0.48, 0.48\n");
      fprintf(gpFile, "set origin 0.02,0.05\n");
      if (gp->chiPlot[1])
        {
          fprintf(gpFile, "set title \"chi1\"\n");
          if (opts->chiRange[1])
            fprintf(gpFile, "set cbrange [%f:%f]\n",
                    opts->chiMin[1], opts->chiMax[1]);
          else
            fprintf(gpFile, "set cbrange [0:*]\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName,
                  gp->iChi[1]);
          fprintf(gpFile, "\n");
        }

      if (gp->chiPlot[3])
        {
          fprintf(gpFile, "set origin 0.5,0.05\n");
          fprintf(gpFile, "set title \"chi12\"\n");
          fprintf(gpFile, "unset ylabel\n");
          if (opts->chiRange[1])
            fprintf(gpFile, "set cbrange [%f:%f]\n", opts->chiMin[3],
                    opts->chiMax[3]);
          else
            fprintf(gpFile, "set cbrange [*:*]\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName,
                  gp->iChi[3]);
          fprintf(gpFile, "\n");
        }

      // chi12
      if (gp->chiPlot[2])
        {
          fprintf(gpFile, "set origin 0.02,0.5\n");
          fprintf(gpFile, "set title \"chi2\"\n");
          fprintf(gpFile, "unset xlabel\n");
          fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
          if (opts->chiRange[2])
            fprintf(gpFile, "set cbrange [%f:%f]\n",
                    opts->chiMin[2],opts->chiMax[2]);
          else
            fprintf(gpFile, "set cbrange [0:*]\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName,
                  gp->iChi[2]);
          fprintf(gpFile, "\n");
        }

      // chi0
      if (gp->chiPlot[1] && gp->chiPlot[2] && gp->chiPlot[3])
        {
          fprintf(gpFile, "set origin 0.5,0.5\n");
          fprintf(gpFile, "set title \"chi0\"\n");
          fprintf(gpFile, "unset ylabel\n");
          if (opts->chiRange[2])
            fprintf(gpFile, "set cbrange [%f:%f]\n",
                    fmin(opts->chiMin[1], opts->chiMin[2]),
                    fmax(opts->chiMax[1], opts->chiMax[2]));
          else
            fprintf(gpFile, "set cbrange [0:*]\n");
          sprintf(buf, "($%i+2*$%i+$%i)",gp->iChi[1],gp->iChi[3],gp->iChi[2]);
          fprintf(gpFile, "splot '%s' using 1:2:%s notitle\n", fileName, buf);
          fprintf(gpFile, "\n");
        }
      fprintf(gpFile, "unset multiplot\n");
    } // end 2D plot script
}
/* -------------------------------------------------------------------------- */

void genPDFscript(char *fileName, optStruct *opts, plotStruct *gp)
{
  int i, j, k, curY, Yremain, nYWindows, nMarg;
  char fileHead[100], outputFile[100], gnuPlot[200], header[100], *head;
  char buffer[100], margFile[100];
    
  // create file name to generate gnuplot script to plot data in this file
  strcpy(fileHead, fileName);
  head = strtok(fileHead, ".");
  createOutputFileName("plot", head, "gp", gnuPlot);
  strcpy(margFile, head);
  strcat(margFile, ".marg");

  FILE *gpFile = fopen(gnuPlot, "w");

  if (gp->plotDim == 1)
    {
      fprintf(gpFile, "set terminal %s 1%s\n",opts->termType,opts->postTerm);
      if (opts->outTerm)
        {
          createOutputFileName("PDF", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }

      fprintf(gpFile, "set title \"PDF\"\n");
      fprintf(gpFile, "set ylabel \"P\"\n");
      fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
      fprintf(gpFile, "set xrange [%f:%f]\n", opts->x1min, opts->x1max);
      if (opts->pdfRange)
        fprintf(gpFile, "set yrange [%f:%f]\n",opts->pdfMin[0],opts->pdfMax[0]);
      else
        fprintf(gpFile, "set yrange [0:*]\n");

      // print all pdfs on same plot
      fprintf(gpFile, "plot '%s' u %i:%i %s t \"%s\"", fileName, 
              opts->xIndex, gp->iMisc, opts->rcIndex, gp->miscNames[0]);

      if (gp->nMisc > 1)
        fprintf(gpFile, ", \\\n");
      for (i=1; i < gp->nMisc; i++)
        {
          fprintf(gpFile, "      '' u %i:%i %s t \"%s\"", 
                  opts->xIndex, gp->iMisc+i, opts->rcIndex, gp->miscNames[i]);
          if (i == gp->nMisc - 1)
            fprintf(gpFile, "\n");
          else
            fprintf(gpFile, ", \\\n");
        }
      fprintf(gpFile, "\n");
    }
  else if (gp->plotDim == 2)
    {
      // set margins and size
      fprintf(gpFile, "set rmargin 0\n");
      fprintf(gpFile, "set lmargin 0\n");
      fprintf(gpFile, "set tmargin 0\n");
      fprintf(gpFile, "set bmargin 0\n");

      fprintf(gpFile, "reset\n");
      fprintf(gpFile, "set size 0.48, 0.48\n");
      fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
      fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
      fprintf(gpFile, "set xrange [%f:%f]\n", opts->x1min, opts->x1max);
      fprintf(gpFile, "set yrange [%f:%f]\n", opts->x2min, opts->x2max);
      fprintf(gpFile, "set label \"%s\" at screen 0.45,0.97\n",gp->timeLabel);

      if ( (opts->plotCol < -1) && (opts->plotRow < -1) )
        {
          fprintf(gpFile, "set pm3d map");
          if (opts->smooth > 0)
            fprintf(gpFile, " interpolate %i,%i", opts->smooth, opts->smooth);
          fprintf(gpFile, "\n");
          fprintf(gpFile, "set palette rgbformulae 22,13,-31\n");
        }
      fprintf(gpFile, "\n");

      /* ----- plot miscellaneous variables ----- */
      if (gp->nMisc%4 == 0)
        nYWindows = gp->nMisc/4;
      else
        nYWindows = gp->nMisc/4 + 1;

      Yremain = gp->nMisc; // number of variables yet to be plotted
      curY=0;                // index of current variable
      // plot 4 species per window, starting from bottom left
      for (k=0; k < nYWindows; k++)
        {
          // set output parameters
          if (opts->outTerm)
            {
              fprintf(gpFile, "set terminal %s %s\n",
                      opts->termType,opts->postTerm);
              sprintf(buffer, "pdf2D_%i", k+1);
              createOutputFileName(buffer, head, opts->termType,
                                   outputFile);
              fprintf(gpFile, "set output \"%s\"\n", outputFile);
            }
          else
            fprintf(gpFile,"set terminal %s %i%s\n",opts->termType,
                    3+opts->nPlotSpecies/4+1+k, opts->postTerm);
          fprintf(gpFile, "set multiplot\n");

          if (opts->logScale) fprintf(gpFile, "set log cb\n");
          for (j=0; j<4; j++, curY++, Yremain--)
            { 
              // if no species left for this window, leave
              if (Yremain == 0)
                break;

              if (j == 0) // bottom left
                {
                  fprintf(gpFile, "set origin 0.02,0.05\n");
                  fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
                  fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
                }
              else if (j == 1) // bottom right
                {
                  fprintf(gpFile, "set origin 0.5,0.05\n");
                  fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
                  fprintf(gpFile, "unset ylabel\n");
                }
              else if (j == 2) // top left
                {
                  fprintf(gpFile, "set origin 0.02,0.5\n");
                  fprintf(gpFile, "unset xlabel\n");
                  fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
                }
              else if (j == 3) // top right
                {
                  fprintf(gpFile, "set origin 0.5,0.5\n");
                  fprintf(gpFile, "unset xlabel\n");
                  fprintf(gpFile, "unset ylabel\n");
                }
              fprintf(gpFile, "\n");

              fprintf(gpFile, "set title \"%s\"\n", gp->miscNames[curY]);

              if (opts->pdfRange)
                fprintf(gpFile, "set cbrange [%f:%f]\n",
                        opts->pdfMin[0],opts->pdfMax[0]);
              else if (opts->Yrange || opts->logScale)
                fprintf(gpFile, "set cbrange [%f:%f]\n",
                        opts->miscMin[curY],opts->miscMax[curY]);
              else
                fprintf(gpFile, "set cbrange [0:*]\n");

              fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n",
                      fileName, gp->iMisc+curY);
              fprintf(gpFile, "\n");
            } // end loop over species for this window

          fprintf(gpFile, "unset multiplot\n");
        } // end loop over windows
      if (opts->logScale) fprintf(gpFile, "unset log cb\n");

      /* plot marginals */
      fprintf(gpFile, "set terminal %s 2%s\n",opts->termType,opts->postTerm);
      if (opts->outTerm)
        {
          createOutputFileName("pdf_marg", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }

      // plot the marginals
      fprintf(gpFile, "unset pm3d\n");
      fprintf(gpFile, "unset title\n");
      fprintf(gpFile, "set multiplot\n");
      fprintf(gpFile, "set size 0.85, 0.35\n");

      // marginal 1
      fprintf(gpFile, "set origin 0.1,0.1\n");
      fprintf(gpFile, "set xrange [%f:%f]\n",opts->x1min, opts->x1max);
      fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
      fprintf(gpFile, "set ylabel \"P|[%s]\"\n", gp->xLabel);
      if (opts->pdfRange)
        fprintf(gpFile, "set yrange [%f:%f]\n",opts->pdfMin[1],opts->pdfMax[1]);
      else
        fprintf(gpFile, "set yrange [0:*]\n");
      fprintf(gpFile, "plot '%s' using 1:2 index 0 w l notitle", margFile);

      if (gp->nMisc > 1 || opts->margBeta)
        fprintf(gpFile, ", \\\n");
      for (i=1; i < gp->nMisc; i++)
        {
          if (!strcmp(gp->miscNames[i], "interp") ||
              !strcmp(gp->miscNames[i], "dns"))
            fprintf(gpFile, "      '' u 1:%i index 0 notitle", i+2);
          else
            fprintf(gpFile, "      '' u 1:%i index 0 w l notitle", i+2);

          if (i == gp->nMisc - 1 && !opts->margBeta)
            fprintf(gpFile, "\n");
          else
            fprintf(gpFile, ", \\\n");
        }
      if (opts->margBeta)
        fprintf(gpFile, "      '' u 1:%i index 0 w l t \"beta\"\n", gp->nMisc+2);
      fprintf(gpFile, "\n");

      // marginal 2
      fprintf(gpFile, "set origin 0.1,0.6\n");
      fprintf(gpFile, "set xlabel \"%s\"\n", gp->yLabel);
      fprintf(gpFile, "set ylabel \"P|[%s]\"\n", gp->yLabel);
      fprintf(gpFile, "set xrange [%f:%f]\n",opts->x2min, opts->x2max);
      if (opts->pdfRange)
        fprintf(gpFile, "set yrange [%f:%f]\n",opts->pdfMin[2],opts->pdfMax[2]);
      else
        fprintf(gpFile, "set yrange [0:*]\n");
      fprintf(gpFile, "plot '%s' using 1:2 index 1 w l t \"%s\"", margFile,
              gp->miscNames[0]);

      if (gp->nMisc > 1 || opts->margBeta)
        fprintf(gpFile, ", \\\n");
      for (i=1; i < gp->nMisc; i++)
        {
          if (!strcmp(gp->miscNames[i], "interp") ||
              !strcmp(gp->miscNames[i], "dns"))
            fprintf(gpFile, "      '' u 1:%i index 1 t \"%s\"", i+2,
                    gp->miscNames[i]);
          else
            fprintf(gpFile, "      '' u 1:%i index 1 w l t \"%s\"", i+2,
                    gp->miscNames[i]);

          if (i == gp->nMisc - 1 && !opts->margBeta)
            fprintf(gpFile, "\n");
          else
            fprintf(gpFile, ", \\\n");
        }
      if (opts->margBeta)
        fprintf(gpFile, "      '' u 1:%i index 1 w l t \"beta\"\n", gp->nMisc+2);
      fprintf(gpFile, "\n");

      fprintf(gpFile, "unset multiplot\n");
    }

  fclose(gpFile);
}
/* -------------------------------------------------------------------------- */

void genChiScript(char *fileName, optStruct *opts, plotStruct *gp)
{
  int i, j, k, curY, Yremain, nYWindows, nMarg;
  char fileHead[100], outputFile[100], gnuPlot[200], header[100], *head;
  char buffer[100], margFile[100];
    
  // create file name to generate gnuplot script to plot data in this file
  strcpy(fileHead, fileName);
  head = strtok(fileHead, ".");
  createOutputFileName("plot", head, "gp", gnuPlot);
  strcpy(margFile, head);
  strcat(margFile, ".chi_marg");

  FILE *gpFile = fopen(gnuPlot, "w");

  if (gp->plotDim == 1)
    {
      fprintf(stderr, "chi script not implemented for 1D\n");
    }
  else if (gp->plotDim == 2)
    {
      // set margins and size
      fprintf(gpFile, "set rmargin 0\n");
      fprintf(gpFile, "set lmargin 0\n");
      fprintf(gpFile, "set tmargin 0\n");
      fprintf(gpFile, "set bmargin 0\n");


      /* if (1 < 0) */ {
      fprintf(gpFile, "reset\n");
      fprintf(gpFile, "set size 0.48, 0.48\n");
      fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
      fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
      fprintf(gpFile, "set xrange [%f:%f]\n", opts->x1min, opts->x1max);
      fprintf(gpFile, "set yrange [%f:%f]\n", opts->x2min, opts->x2max);
      fprintf(gpFile, "set label \"%s\" at screen 0.45,0.97\n",gp->timeLabel);

      if ( (opts->plotCol < -1) && (opts->plotRow < -1) )
        {
          fprintf(gpFile, "set pm3d map");
          if (opts->smooth > 0)
            fprintf(gpFile, " interpolate %i,%i", opts->smooth, opts->smooth);
          fprintf(gpFile, "\n");
          fprintf(gpFile, "set palette rgbformulae 22,13,-31\n");
        }
      fprintf(gpFile, "\n");

      /* ----- plot scalar dissipation rates ----- */
      if (opts->outTerm)
        {
          fprintf(gpFile, "set terminal %s %s\n",opts->termType,opts->postTerm);
          createOutputFileName("chi", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }
      else
        fprintf(gpFile,"set terminal %s 1%s\n", opts->termType, opts->postTerm);

      fprintf(gpFile, "set multiplot\n");
      fprintf(gpFile, "set size 0.48, 0.48\n");
      fprintf(gpFile, "set origin 0.02,0.05\n");
      if (gp->chiPlot[1] && gp->chiPlot[2] && gp->chiPlot[3])   // chi0
        {

          fprintf(gpFile, "set title \"chi0\"\n");
          if (opts->chiRange[2])
            fprintf(gpFile, "set cbrange [%f:%f]\n",
                    fmin(opts->chiMin[1], opts->chiMin[2]),
                    fmax(opts->chiMax[1], opts->chiMax[2]));
          else
            fprintf(gpFile, "set cbrange [0:*]\n");
          sprintf(buffer,"($%i+2*$%i+$%i)",gp->iChi[1],gp->iChi[3],gp->iChi[2]);
          fprintf(gpFile, "splot '%s' using 1:2:%s notitle\n",fileName,buffer);
          fprintf(gpFile, "\n");
        }
      if (gp->chiPlot[1])
        {
          fprintf(gpFile, "set origin 0.5,0.05\n");
          fprintf(gpFile, "set title \"chi1\"\n");
          fprintf(gpFile, "unset ylabel\n");
          if (opts->chiRange[1])
            fprintf(gpFile, "set cbrange [%f:%f]\n",
                    opts->chiMin[1], opts->chiMax[1]);
          else
            fprintf(gpFile, "set cbrange [0:*]\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName,
                  gp->iChi[1]);
          fprintf(gpFile, "\n");
        }
      if (gp->chiPlot[2])
        {
          fprintf(gpFile, "set origin 0.02,0.5\n");
          fprintf(gpFile, "set title \"chi2\"\n");
          fprintf(gpFile, "unset xlabel\n");
          fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
          if (opts->chiRange[2])
            fprintf(gpFile, "set cbrange [%f:%f]\n",
                    opts->chiMin[2],opts->chiMax[2]);
          else
            fprintf(gpFile, "set cbrange [0:*]\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName,
                  gp->iChi[2]);
          fprintf(gpFile, "\n");
        }
      if (gp->chiPlot[3]) // chi12
        {

          fprintf(gpFile, "set origin 0.5,0.5\n");
          fprintf(gpFile, "set title \"chi12\"\n");
          fprintf(gpFile, "unset ylabel\n");
          if (opts->chiRange[1])
            fprintf(gpFile, "set cbrange [%f:%f]\n", opts->chiMin[3],
                    opts->chiMax[3]);
          else
            fprintf(gpFile, "set cbrange [*:*]\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName,
                  gp->iChi[3]);
          fprintf(gpFile, "\n");
        }
      fprintf(gpFile, "unset multiplot\n\n");

      /* ----- plot log chi ----- */
      if (opts->outTerm)
        {
          fprintf(gpFile, "set terminal %s %s\n",opts->termType,opts->postTerm);
          createOutputFileName("chi_log", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }
      else
        fprintf(gpFile,"set terminal %s 2%s\n", opts->termType, opts->postTerm);

      fprintf(gpFile, "set multiplot\n");
      fprintf(gpFile, "set size 0.48, 0.48\n");
      fprintf(gpFile, "set origin 0.02,0.05\n");
      if (gp->chiPlot[1] && gp->chiPlot[2] && gp->chiPlot[3])   // chi0
        {

          fprintf(gpFile, "set title \"chi0-log\"\n");
          if (opts->chiRange[2])
            fprintf(gpFile, "set cbrange [%f:%f]\n",
                    fmin(opts->miscMin[1], opts->miscMin[3]),
                    fmax(opts->miscMax[1], opts->miscMax[3]));
          else
            fprintf(gpFile, "set cbrange [0:*]\n");
          sprintf(buffer,"($%i+2*$%i+$%i)",gp->iMisc+1,gp->iMisc+2,gp->iMisc+3);
          fprintf(gpFile, "splot '%s' using 1:2:%s notitle\n",fileName,buffer);
          fprintf(gpFile, "\n");
        }
      if (gp->chiPlot[1])
        {
          fprintf(gpFile, "set origin 0.5,0.05\n");
          fprintf(gpFile, "set title \"chi1-log\"\n");
          fprintf(gpFile, "unset ylabel\n");
          if (opts->chiRange[1])
            fprintf(gpFile, "set cbrange [%f:%f]\n",
                    opts->miscMin[1], opts->miscMax[1]);
          else
            fprintf(gpFile, "set cbrange [0:*]\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName,
                  gp->iMisc+1);
          fprintf(gpFile, "\n");
        }
      if (gp->chiPlot[2])
        {
          fprintf(gpFile, "set origin 0.02,0.5\n");
          fprintf(gpFile, "set title \"chi2-log\"\n");
          fprintf(gpFile, "unset xlabel\n");
          fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
          if (opts->chiRange[2])
            fprintf(gpFile, "set cbrange [%f:%f]\n",
                    opts->miscMin[3],opts->miscMax[3]);
          else
            fprintf(gpFile, "set cbrange [0:*]\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName,
                  gp->iMisc+3);
          fprintf(gpFile, "\n");
        }
      if (gp->chiPlot[3]) // chi12
        {

          fprintf(gpFile, "set origin 0.5,0.5\n");
          fprintf(gpFile, "set title \"chi12-log\"\n");
          fprintf(gpFile, "unset ylabel\n");
          if (opts->chiRange[1])
            fprintf(gpFile, "set cbrange [%f:%f]\n", opts->miscMin[2],
                    opts->miscMax[2]);
          else
            fprintf(gpFile, "set cbrange [*:*]\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName,
                  gp->iMisc+2);
          fprintf(gpFile, "\n");
        }
      fprintf(gpFile, "unset multiplot\n\n");

      /* ----- plot log erfc ----- */
      if (opts->outTerm)
        {
          fprintf(gpFile, "set terminal %s %s\n",opts->termType,opts->postTerm);
          createOutputFileName("chi_erfc", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }
      else
        fprintf(gpFile,"set terminal %s 3%s\n", opts->termType, opts->postTerm);

      fprintf(gpFile, "set multiplot\n");
      fprintf(gpFile, "set size 0.48, 0.48\n");
      fprintf(gpFile, "set origin 0.02,0.05\n");
      if (gp->chiPlot[1] && gp->chiPlot[2] && gp->chiPlot[3])   // chi0
        {

          fprintf(gpFile, "set title \"chi0-erfc\"\n");
          if (opts->chiRange[2])
            fprintf(gpFile, "set cbrange [%f:%f]\n",
                    fmin(opts->miscMin[5], opts->miscMin[7]),
                    fmax(opts->miscMax[5], opts->miscMax[7]));
          else
            fprintf(gpFile, "set cbrange [0:*]\n");
          sprintf(buffer,"($%i+2*$%i+$%i)",gp->iMisc+5,gp->iMisc+6,gp->iMisc+7);
          fprintf(gpFile, "splot '%s' using 1:2:%s notitle\n",fileName,buffer);
          fprintf(gpFile, "\n");
        }
      if (gp->chiPlot[1])
        {
          fprintf(gpFile, "set origin 0.5,0.05\n");
          fprintf(gpFile, "set title \"chi1-erfc\"\n");
          fprintf(gpFile, "unset ylabel\n");
          if (opts->chiRange[1])
            fprintf(gpFile, "set cbrange [%f:%f]\n",
                    opts->miscMin[5], opts->miscMax[5]);
          else
            fprintf(gpFile, "set cbrange [0:*]\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName,
                  gp->iMisc+5);
          fprintf(gpFile, "\n");
        }
      if (gp->chiPlot[2])
        {
          fprintf(gpFile, "set origin 0.02,0.5\n");
          fprintf(gpFile, "set title \"chi2-erfc\"\n");
          fprintf(gpFile, "unset xlabel\n");
          fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
          if (opts->chiRange[2])
            fprintf(gpFile, "set cbrange [%f:%f]\n",
                    opts->miscMin[7],opts->miscMax[7]);
          else
            fprintf(gpFile, "set cbrange [0:*]\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName,
                  gp->iMisc+7);
          fprintf(gpFile, "\n");
        }
      if (gp->chiPlot[3]) // chi12
        {

          fprintf(gpFile, "set origin 0.5,0.5\n");
          fprintf(gpFile, "set title \"chi12-erfc\"\n");
          fprintf(gpFile, "unset ylabel\n");
          if (opts->chiRange[1])
            fprintf(gpFile, "set cbrange [%f:%f]\n", opts->miscMin[6],
                    opts->miscMax[6]);
          else
            fprintf(gpFile, "set cbrange [*:*]\n");
          fprintf(gpFile, "splot '%s' using 1:2:%i notitle\n", fileName,
                  gp->iMisc+6);
          fprintf(gpFile, "\n");
        }
      fprintf(gpFile, "unset multiplot\n\n");

      /* ----- plot log thetas ----- */
      if (opts->outTerm)
        {
          fprintf(gpFile, "set terminal %s %s\n",opts->termType,opts->postTerm);
          createOutputFileName("chi_theta", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }
      else
        fprintf(gpFile,"set terminal %s 4%s\n", opts->termType, opts->postTerm);

      fprintf(gpFile, "set multiplot\n");
      fprintf(gpFile, "set size 0.48, 0.48\n");
      fprintf(gpFile, "set origin 0.02,0.05\n");

      /* theta dns*/
      fprintf(gpFile, "set title \"theta dns\"\n");
      if (opts->chiRange[2])
        fprintf(gpFile, "set cbrange [%f:%f]\n", opts->miscMin[0],
                opts->miscMax[0]);
      else
        fprintf(gpFile, "set cbrange [0:*]\n");
      sprintf(buffer,"($%i*180/pi)",gp->iMisc);
      fprintf(gpFile, "splot '%s' using 1:2:%s notitle\n",fileName,buffer);
      fprintf(gpFile, "\n");

      /* theta log */
      fprintf(gpFile, "set origin 0.5,0.05\n");
      fprintf(gpFile, "set title \"theta log\"\n");
      fprintf(gpFile, "unset ylabel\n");
      if (opts->chiRange[1])
        fprintf(gpFile, "set cbrange [%f:%f]\n",
                opts->miscMin[4], opts->miscMax[4]);
      else
        fprintf(gpFile, "set cbrange [0:*]\n");
      sprintf(buffer,"($%i*180/pi)",gp->iMisc+4);
      fprintf(gpFile, "splot '%s' using 1:2:%s notitle\n",fileName,buffer);
      fprintf(gpFile, "\n");

      /* theta erfc */
      fprintf(gpFile, "set origin 0.02,0.5\n");
      fprintf(gpFile, "set title \"theta erfc\"\n");
      fprintf(gpFile, "unset xlabel\n");
      fprintf(gpFile, "set ylabel \"%s\"\n", gp->yLabel);
      if (opts->chiRange[2])
        fprintf(gpFile, "set cbrange [%f:%f]\n",
                opts->miscMin[8],opts->miscMax[8]);
      else
        fprintf(gpFile, "set cbrange [0:*]\n");
      sprintf(buffer,"($%i*180/pi)",gp->iMisc+8);
      fprintf(gpFile, "splot '%s' using 1:2:%s notitle\n",fileName,buffer);
      fprintf(gpFile, "\n");
      fprintf(gpFile, "unset multiplot\n\n");
      }

      /* ----- plot the marginals ----- */
      if (opts->outTerm)
        {
          fprintf(gpFile, "set terminal %s %s\n",opts->termType,opts->postTerm);
          createOutputFileName("chi_marg", head, opts->termType, outputFile);
          fprintf(gpFile, "set output \"%s\"\n", outputFile);
        }
      else
        fprintf(gpFile,"set terminal %s 5%s\n", opts->termType, opts->postTerm);

      fprintf(gpFile, "unset pm3d\n");
      fprintf(gpFile, "unset title\n");

      fprintf(gpFile, "set size 0.85, 0.35\n");

/*       // marginal 1 */
/*       fprintf(gpFile, "set origin 0.1,0.1\n"); */
/*       fprintf(gpFile, "set xrange [%f:%f]\n",opts->x1min, opts->x1max); */
/*       fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel); */
/*       fprintf(gpFile, "set ylabel \"chi1|[%s]\"\n", gp->xLabel); */
/*       fprintf(gpFile, "set yrange [0:*]\n"); */
/*       fprintf(gpFile, "plot '%s' u 1:3 index 0 t \"dns\", \\\n",margFile); */
/*       fprintf(gpFile, " '' u 1:4 index 0 w l lc -1 t \"log\", \\\n"); */
/*       fprintf(gpFile, " '' u 1:5 index 0 w l t \"erfc\"\n"); */
/*       fprintf(gpFile, "\n"); */

/*       // marginal 2 */
/*       fprintf(gpFile, "set origin 0.1,0.6\n"); */
/*       fprintf(gpFile, "set xlabel \"%s\"\n", gp->yLabel); */
/*       fprintf(gpFile, "set ylabel \"chi2|[%s]\"\n", gp->yLabel); */
/*       fprintf(gpFile, "set xrange [%f:%f]\n",opts->x2min, opts->x2max); */
/*       fprintf(gpFile, "set yrange [0:*]\n"); */
/*       fprintf(gpFile, "plot '%s' u 1:3 index 1 t \"dns\", \\\n",margFile); */
/*       fprintf(gpFile, " '' u 1:4 index 1 w l lc -1 t \"log\", \\\n"); */
/*       fprintf(gpFile, " '' u 1:5 index 1 w l t \"erfc\"\n"); */
/*       fprintf(gpFile, "unset multiplot\n"); */

      /* ----- plot the marginal chi distributions ----- */
      // marginal 1
      fprintf(gpFile, "set multiplot\n");
      fprintf(gpFile, "set pointsize 0.5\n");
      fprintf(gpFile, "set origin 0.1,0.1\n");
      fprintf(gpFile, "set xrange [0:1]\n");
      fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
      fprintf(gpFile, "set ylabel \"chi1|[%s]\"\n", gp->xLabel);
      fprintf(gpFile, "set yrange [0:0.05]\n");
      fprintf(gpFile, "plot '%s' u 2:3 lc -1 t \"dns\", \\\n",fileName);
      fprintf(gpFile, " '%s' u 1:8 index 1 w l lc 1 lw 2 t \"dns-mean\", \\\n",
              margFile);
      fprintf(gpFile, " '' u 1:10 index 1 w l t \"log\", \\\n");
      fprintf(gpFile, " '' u 1:12 index 1 w l t \"erfc\", \\\n");
      fprintf(gpFile, " '' u 1:14 index 1 w l t \"poly\"\n");
      fprintf(gpFile, "\n");

      // marginal 2
      fprintf(gpFile, "set origin 0.1,0.6\n");
      fprintf(gpFile, "set xlabel \"%s\"\n", gp->yLabel);
      fprintf(gpFile, "set ylabel \"chi2|[%s]\"\n", gp->yLabel);
      fprintf(gpFile, "set xrange [0:1]\n");
      fprintf(gpFile, "plot '%s' u 1:3 lc -1 notitle, \\\n",fileName);
      fprintf(gpFile, " '%s' u 1:8 index 0 w l lc 1 lw 2 notitle, \\\n",
              margFile);
      fprintf(gpFile, " '' u 1:10 index 0 w l notitle, \\\n");
      fprintf(gpFile, " '' u 1:12 index 0 w l notitle, \\\n");
      fprintf(gpFile, " '' u 1:14 index 0 w l notitle\n");
     
      fprintf(gpFile, "unset multiplot\n\n");

    }

  fclose(gpFile);
}
/* -------------------------------------------------------------------------- */

void genChiMargScript(char *fileName, optStruct *opts, plotStruct *gp)
{
  int i, j, k, curY, Yremain, nYWindows, nMarg;
  char fileHead[100], outputFile[100], gnuPlot[200], header[100], *head;
  char buffer[100], margFile[100];
    
  // create file name to generate gnuplot script to plot data in this file
  strcpy(fileHead, fileName);
  head = strtok(fileHead, ".");
  createOutputFileName("plot", head, "gp", gnuPlot);

  FILE *gpFile = fopen(gnuPlot, "w");

  if (opts->outTerm)
    {
      fprintf(gpFile, "set terminal %s %s\n",opts->termType,opts->postTerm);
      createOutputFileName("chi_marg", head, opts->termType, outputFile);
      fprintf(gpFile, "set output \"%s\"\n", outputFile);
    }
  else
    fprintf(gpFile,"set terminal %s 1%s\n", opts->termType, opts->postTerm);

  fprintf(gpFile, "set multiplot\n");
  fprintf(gpFile, "set size 0.85, 0.35\n");

  // marginal 1
  fprintf(gpFile, "set origin 0.1,0.1\n");
  fprintf(gpFile, "set xrange [%f:%f]\n",opts->x1min, opts->x1max);
  fprintf(gpFile, "set xlabel \"%s\"\n", gp->xLabel);
  fprintf(gpFile, "set ylabel \"chi1|[%s]\"\n", gp->xLabel);
  fprintf(gpFile, "set yrange [0:*]\n");
  fprintf(gpFile, "plot '%s' u 1:3 index 0 t \"dns\", \\\n",margFile);
  fprintf(gpFile, " '' u 1:4 index 0 w l lc -1 t \"log\", \\\n");
  fprintf(gpFile, " '' u 1:5 index 0 w l t \"erfc\"\n");
  fprintf(gpFile, "\n");

  // marginal 2
  fprintf(gpFile, "set origin 0.1,0.6\n");
  fprintf(gpFile, "set xlabel \"%s\"\n", gp->yLabel);
  fprintf(gpFile, "set ylabel \"chi2|[%s]\"\n", gp->yLabel);
  fprintf(gpFile, "set xrange [%f:%f]\n",opts->x2min, opts->x2max);
  fprintf(gpFile, "set yrange [0:*]\n");
  fprintf(gpFile, "plot '%s' u 1:3 index 1 t \"dns\", \\\n",margFile);
  fprintf(gpFile, " '' u 1:4 index 1 w l lc -1 t \"log\", \\\n");
  fprintf(gpFile, " '' u 1:5 index 1 w l t \"erfc\"\n");
  fprintf(gpFile, "unset multiplot\n");


  fclose(gpFile);
}
/* -------------------------------------------------------------------------- */
