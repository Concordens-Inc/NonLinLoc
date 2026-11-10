/*
 * Copyright (C) 2026 Concordens Inc.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser Public License for more details.

 * You should have received a copy of the GNU Lesser Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */


/*   MeshTimeRead.c

        Read travel times from NonLinLoc time grids at given points, through
        the reader NLLoc uses (OpenGrid3dFile, ReadAbsInterpGrid3d), for
        regular (TIME) and mesh (TIME_MESH) grids alike.

        usage: MeshTimeRead <grid root> [<grid root> ...] < points
          <grid root>  e.g. times/NET.P.STA.time (without .hdr / .buf)
          points       one "x y z" (km) per line on standard input
        output: one line per point, the time (s) from each grid, "nan" where
        the grid has no time (outside a mesh).
 */

#define PNAME  "MeshTimeRead"

#include "GridLib.h"
#include "MeshTimeLib.h"

int main(int argc, char *argv[]) {

    if (argc < 2) {
        fprintf(stderr, "usage: %s <grid root> [<grid root> ...] < points\n", PNAME);
        return (EXIT_ERROR_USAGE);
    }

    message_flag = 0;
    SetConstants();

    int ngrid = argc - 1;
    GridDesc *grids = (GridDesc *) calloc((size_t) ngrid, sizeof (GridDesc));
    FILE **fp_grid = (FILE **) calloc((size_t) ngrid, sizeof (FILE *));
    FILE **fp_hdr = (FILE **) calloc((size_t) ngrid, sizeof (FILE *));
    SourceDesc srce;
    int n;
    for (n = 0; n < ngrid; n++) {
        if (OpenGrid3dFile(argv[n + 1], fp_grid + n, fp_hdr + n, grids + n, "time", &srce, 0) < 0) {
            fprintf(stderr, "%s: cannot open grid %s\n", PNAME, argv[n + 1]);
            return (EXIT_ERROR_FILEIO);
        }
    }

    char line[MAXLINE_LONG];
    double x, y, z;
    while (fgets(line, MAXLINE_LONG, stdin) != NULL) {
        if (sscanf(line, "%lf %lf %lf", &x, &y, &z) != 3)
            continue;
        for (n = 0; n < ngrid; n++) {
            GRID_FLOAT_TYPE t = ReadAbsInterpGrid3d(fp_grid[n], grids + n, x, y, z, 0);
            if (t < -LARGE_FLOAT)
                printf("%snan", n ? " " : "");
            else
                printf("%s%.9g", n ? " " : "", (double) t);
        }
        printf("\n");
    }

    for (n = 0; n < ngrid; n++)
        CloseGrid3dFile(grids + n, fp_grid + n, fp_hdr + n);

    return (0);
}
