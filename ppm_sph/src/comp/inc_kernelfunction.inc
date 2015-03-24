#ifndef PARAM_KERNEL_TYPE_QUINTICSPLINE
#define PARAM_KERNEL_TYPE_QUINTICSPLINE 1
#endif
#ifndef PARAM_KERNEL_TYPE_BSPLINE
#define PARAM_KERNEL_TYPE_BSPLINE 2
#endif
#ifndef PARAM_KERNEL_TYPE_CUBIC
#define PARAM_KERNEL_TYPE_CUBIC 3
#endif
#ifndef PARAM_KERNEL_TYPE_WENDLANDQUARTIC
#define PARAM_KERNEL_TYPE_WENDLANDQUARTIC 4
#endif
#if __KERNEL==PARAM_KERNEL_TYPE_QUINTICSPLINE
#if __DIM==2
#include "inc_kernelfunction_quinticspline_2D.inc"
#elif __DIM==3
#include "inc_kernelfunction_quinticspline_3D.inc"
#endif

#elif __KERNEL==PARAM_KERNEL_TYPE_CUBIC
#if __DIM==2
#include "inc_kernelfunction_cubicspline_2D.inc"
#elif __DIM==3
#include "inc_kernelfunction_cubicspline_3D.inc"
#endif

#elif __KERNEL==PARAM_KERNEL_TYPE_BSPLINE

#elif __KERNEL==PARAM_KERNEL_TYPE_WENDLANDQUARTIC
#if __DIM==2
#include "inc_kernelfunction_wendlandquartic_2D.inc"
#endif

#endif
