/*----------------------------------------------------------------------------
! This file is adapted from psOpen, modified by Tobias Falkenstein
! Modificatinons only tested on INTEL platforms:
!
! Makefile.in:
!               CFLAGS   = -DINTEL
! 
!
! Version 0.9
!
! Copyright (C) 2013 Jens Henrik Goebbert <jens.henrik.goebbert()rwth-aachen.de>
!
!    psOpen is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!    Please consider mentioning the copyright owner in any publication
!    if results of psOpen are used.
!----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#ifdef IBM

#include <spi/include/kernel/memory.h>

extern
void print_memusage_c()
{
  uint64_t shared, persist, heapavail, stackavail, stack, heap, guard, mmap;

  Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &guard);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &persist);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &stackavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_MMAP, &mmap);

  printf("   current MEMSIZE heap  : %.2f/%.2f stack: %.2f/%.2f mmap: %.2f mbyte\n", (double)heap/(1024*1024), (double)heapavail/(1024*1024),
                                                                              (double)stack/(1024*1024), (double)stackavail/(1024*1024),
                                                                              (double)mmap/(1024*1024));
  printf("   current MEMSIZE shared: %.2f persist: %.2f guard: %.2f mbyte\n", (double)shared/(1024*1024),
                                                                       (double)persist/(1024*1024),
                                                                        (double)guard/(1024*1024));

}

extern
void get_memusage_c(double *rss_out, double *vsize_out)
{
  if (rss_out != NULL) *rss_out = 0.0;
  if (vsize_out != NULL) *vsize_out = 0.0;
}
#elif INTEL

#include <unistd.h> // for sysconf(_SC_PAGESIZE)
#include <linux/limits.h> // for PATH_MAX

void readone(FILE *f, long long int *x)            { fscanf(f, "%lld ", x); }
void readunsigned(FILE *f, unsigned long long *x) { fscanf(f, "%llu ", x); }
void readstr(FILE *f, char *x)                      { fscanf(f, "%s ", x);}
void readchar(FILE *f, char *x)                     { fscanf(f, "%c ", x);}

extern
void get_memusage_c(double *rss_out, double *vsize_out)
{

  long page_size = sysconf(_SC_PAGESIZE);
  long s = -1;
  FILE *f = fopen("/proc/self/stat", "r");
  if (!f) return ;

  // example: 24773 (cat) R
  long long int pid;    readone(f,&pid);    // process id
  char tcomm[PATH_MAX];  readstr(f,tcomm);   // process name
  char state;            readchar(f,&state); // process status (R==running, ...)

  // example: 7627 24773 7627 34827 24773
  long long int ppid;     readone(f,&ppid);
  long long int pgid;     readone(f,&pgid);
  long long int sid;      readone(f,&sid);
  long long int tty_nr;   readone(f,&tty_nr);
  long long int tty_pgrp; readone(f,&tty_pgrp);

  // example: 4202496 220 0 0 0 0 0
  long long int flags;    readone(f,&flags);
  long long int min_flt;  readone(f,&min_flt);
  long long int cmin_flt; readone(f,&cmin_flt);
  long long int maj_flt;  readone(f,&maj_flt);
  long long int cmaj_flt; readone(f,&cmaj_flt);
  long long int utime;    readone(f,&utime);
  long long int stimev;   readone(f,&stimev);

  // example: 0 0 20 0 1 0
  long long int cutime;        readone(f,&cutime);
  long long int cstime;        readone(f,&cstime);
  long long int priority;      readone(f,&priority);    // process priority
  long long int nicev;         readone(f,&nicev);       // process nice value
  long long int num_threads;   readone(f,&num_threads); // no. process threads
  long long int it_real_value; readone(f,&it_real_value);

  // example: 29350452
  unsigned long long start_time; readunsigned(f,&start_time); // start time as UNIX time

  // example: 4210688 123 18446744073709551615 4194304 4235780 140734041535920 140734041533112 215179179504
  long long int vsize;      readone(f,&vsize);                                // virtual memory size (in bytes)
  long long int rss;        readone(f,&rss);        rss   = rss   *page_size; // resident set size (rss) is the portion of a process's memory that is held in RAM.
  long long int rsslim;     readone(f,&rsslim);     rsslim= rsslim*page_size; // limit for resident set size (rss) (in pages)
  long long int start_code; readone(f,&start_code);
  long long int end_code;   readone(f,&end_code);
  long long int start_stack;readone(f,&start_stack);
  long long int esp;        readone(f,&esp);
  long long int eip;        readone(f,&eip);

  // example: 0 0 0 0 0 0 0 17 3 0 0
  long long int pending;     readone(f,&pending);
  long long int blocked;     readone(f,&blocked);
  long long int sigign;      readone(f,&sigign);
  long long int sigcatch;    readone(f,&sigcatch);
  long long int wchan;       readone(f,&wchan);
  long long int zero1;       readone(f,&zero1);
  long long int zero2;       readone(f,&zero2);
  long long int exit_signal; readone(f,&exit_signal);
  long long int cpu;         readone(f,&cpu);
  long long int rt_priority; readone(f,&rt_priority);
  long long int policy;      readone(f,&policy);

  // example: 0 0 0
  // unknown

  fclose (f);

  if (rss_out != NULL) *rss_out = (double)rss/(1024*1024);
  if (vsize_out != NULL) *vsize_out = (double)vsize/(1024*1024);

}

extern
void print_memusage_c()
{
  double rss = 0.0;
  double vsize = 0.0;

  get_memusage_c(&rss,&vsize);
  printf("   current MEMSIZE RSS  : %.2f mbyte\n", rss);
  printf("   current MEMSIZE VSIZE: %.2f mbyte\n", vsize);
}

#else
extern
void get_memusage_c(double *rss_out, double *vsize_out)
{
  if (rss_out != NULL) *rss_out = 0.0;
  if (vsize_out != NULL) *vsize_out = 0.0;
}

extern
void print_memusage_c()
{
  double rss = 0.0;
  double vsize = 0.0;

  get_memusage_c(&rss,&vsize);
  printf("   current MEMSIZE RSS  : %.2f mbyte\n", rss);
  printf("   current MEMSIZE VSIZE: %.2f mbyte\n", vsize);
}
#endif

// fortran interfaces

void print_memusage_f_() { print_memusage_c(); }
void print_memusage_f__() { print_memusage_f_(); }
void print_memusage_f() { print_memusage_f_(); }

// NOTE: currently only available for intel platforms
void get_memusage_f_(double *rss, double *vsize) { get_memusage_c(rss,vsize); }
void get_memusage_f__(double *rss, double *vsize) { get_memusage_f_(rss,vsize); }
void get_memusage_f(double *rss, double *vsize) { get_memusage_f_(rss,vsize); }
