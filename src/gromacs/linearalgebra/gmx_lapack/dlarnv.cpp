#include <cmath>

#include "../gmx_lapack.h"

void F77_FUNC(dlarnv, DLARNV)(int* idist, int* iseed, int* n, double* x)
{
    int i__1, i__2, i__3;

    int    i__;
    double u[128];
    int    il, iv, il2;

    --x;
    --iseed;

    i__1 = *n;
    for (iv = 1; iv <= i__1; iv += 64)
    {
        i__2 = 64, i__3 = *n - iv + 1;
        il = (i__2 < i__3) ? i__2 : i__3;
        if (*idist == 3)
        {
            il2 = il << 1;
        }
        else
        {
            il2 = il;
        }

        F77_FUNC(dlaruv, DLARUV)(&iseed[1], &il2, u);

        if (*idist == 1)
        {

            i__2 = il;
            for (i__ = 1; i__ <= i__2; ++i__)
            {
                x[iv + i__ - 1] = u[i__ - 1];
            }
        }
        else if (*idist == 2)
        {

            i__2 = il;
            for (i__ = 1; i__ <= i__2; ++i__)
            {
                x[iv + i__ - 1] = u[i__ - 1] * 2. - 1.;
            }
        }
        else if (*idist == 3)
        {

            i__2 = il;
            for (i__ = 1; i__ <= i__2; ++i__)
            {
                x[iv + i__ - 1] = std::sqrt(std::log(u[(i__ << 1) - 2]) * -2.)
                                  * std::cos(u[(i__ << 1) - 1] * (double)6.2831853071795864769252867663);
            }
        }
    }
    return;
}
