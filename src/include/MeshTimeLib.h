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


/*   MeshTimeLib.h

        Travel-time grids on a tetrahedral mesh (grid type TIME_MESH).

        A TIME_MESH grid holds one travel time per node of a tetrahedral mesh
        shared by every station and phase. Times between nodes are read
        linearly in the containing tetrahedron, after removing the
        straight-ray time s0 |x - x_s| from the source (station) when the
        header gives the slowness at the source (factored read). The mesh,
        its point-location index and every time buffer are memory-mapped
        read-only, so the processes of one machine share one physical copy
        through the page cache.

        Files of a mesh, base path M (no extension):
          M.mesh.hdr        TET_MESH <num_nodes> <num_tets>
                            BINS <nbx> <nby> <nbz> <x0> <y0> <z0> <bdx> <bdy> <bdz> <num_entries>
          M.nodes.buf       float64 (num_nodes, 3): x, y, z (km, the grids' frame)
          M.tets.buf        int32 (num_tets, 4): node indices
          M.bin_start.buf   int64 (nbx * nby * nbz + 1): CSR offsets, bin (ix, iy, iz)
                            at index (ix * nby + iy) * nbz + iz
          M.bin_tets.buf    int32 (num_entries): tetrahedra overlapping each bin

        Files of a time grid, root R (e.g. <root>.P.STA.time):
          R.hdr             2 2 2  xmin ymin zmin  dx dy dz  TIME_MESH FLOAT
                            <station label> <x> <y> <z>
                            TRANSFORM ...
                            MESH <M, absolute or relative to the header's directory>
                            SOURCE_SLOWNESS <s0 (s/km); 0 reads the times linearly>
                            SOURCE_REGION <n>   (optional)
          R.buf             GRID_FLOAT_TYPE (num_nodes): time (s) per node
          R.fixed.buf       int32 (n), sorted: the source-region nodes (optional). In a
                            tetrahedron whose four nodes are all source-region nodes the
                            time is |x - x_s| times the linear read of T_i / |x_i - x_s|.

        The first header line's box is the mesh's bounding box, so the
        regular-grid checks (search grid inside the time grid, station
        distance) apply unchanged. A point inside the box but outside the
        mesh (above the free surface) has no time: the read returns
        -VERY_LARGE_FLOAT, and the oct-tree gives the sample zero likelihood.
 */

#ifndef MESH_TIME_LIB_H
#define MESH_TIME_LIB_H

#include <stdint.h>
#include <stdio.h>

/* GridLib.h (GridDesc, GRID_FLOAT_TYPE) has no include guard: include it first. */

typedef struct {
    char path[FILENAME_MAX]; /* base path of the mesh files (no extension) */
    int num_nodes;
    int num_tets;
    const double *nodes; /* (num_nodes, 3) */
    const int32_t *tets; /* (num_tets, 4) */
    int nbx, nby, nbz; /* point-location bins */
    double bx0, by0, bz0; /* origin of bin (0, 0, 0) */
    double bdx, bdy, bdz; /* bin size */
    int64_t num_entries;
    const int64_t *bin_start; /* (nbx * nby * nbz + 1) */
    const int32_t *bin_tets; /* (num_entries) */
    /* last located point; NonLinLoc reads every arrival at one sample in turn */
    int last_valid;
    double last_x, last_y, last_z;
    int last_tet;
    double last_w[4];
    /* statistics */
    long num_locate;
    long num_locate_cached;
} MeshDesc;

/* Map the mesh and the time buffer of the TIME_MESH grid whose header file is
 * fn_hdr (already open as fp_hdr); fills the grid's mesh fields. 0 on success. */
int MeshTime_Open(GridDesc *pgrid, const char *fn_hdr, FILE *fp_hdr);

/* Travel time at (x, y, z); -VERY_LARGE_FLOAT outside the mesh. */
GRID_FLOAT_TYPE MeshTime_Read(GridDesc *pgrid, double x, double y, double z);

/* Containing tetrahedron and barycentric weights of (x, y, z); 1 if found. */
int MeshTime_Locate(MeshDesc *mesh, double x, double y, double z, int *ptet, double w[4]);

/* Load (once per process) the mesh with base path `path`; NULL on error. */
MeshDesc *MeshTime_LoadMesh(const char *path);

#endif /* MESH_TIME_LIB_H */
