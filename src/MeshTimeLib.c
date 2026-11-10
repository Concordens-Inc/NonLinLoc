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


/*   MeshTimeLib.c

        Travel-time grids on a tetrahedral mesh (grid type TIME_MESH); the
        file formats are described in include/MeshTimeLib.h.

        Every buffer is mapped read-only (PROT_READ, MAP_SHARED) and stays
        mapped for the life of the process: the meshes and the time buffers
        are registered by path, so a grid re-opened for the next event, or
        for the final hypocentre's predicted times, reuses its mapping.

        Not thread-safe: the last located point is cached in the mesh
        (NonLinLoc locates one sample at a time).
 */

#include <errno.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "GridLib.h"
#include "MeshTimeLib.h"

#define MESH_BARY_TOL 1.0e-9

/* ------------------------------------------------------------------------- */
/* read-only file mappings, registered by path */

typedef struct {
    char path[FILENAME_MAX];
    const void *ptr;
    size_t bytes;
} MappedFile;

static MappedFile *mapped_files = NULL;
static int num_mapped_files = 0;
static int cap_mapped_files = 0;

static MeshDesc **meshes = NULL;
static int num_meshes = 0;
static int cap_meshes = 0;

/* Map `path` read-only; its size must be `bytes` (> 0). NULL on error. */
static const void *map_file(const char *path, size_t bytes) {

    int n;
    for (n = 0; n < num_mapped_files; n++) {
        if (strcmp(mapped_files[n].path, path) == 0) {
            if (mapped_files[n].bytes != bytes) {
                sprintf(MsgStr, "ERROR: mesh grid file %s: %zu bytes mapped, %zu expected",
                        path, mapped_files[n].bytes, bytes);
                nll_puterr(MsgStr);
                return (NULL);
            }
            return (mapped_files[n].ptr);
        }
    }

    if (bytes == 0) {
        nll_puterr2("ERROR: empty mesh grid file:", path);
        return (NULL);
    }
    int fd = open(path, O_RDONLY);
    if (fd < 0) {
        nll_puterr2("ERROR: opening mesh grid file:", path);
        return (NULL);
    }
    struct stat st;
    if (fstat(fd, &st) != 0 || (size_t) st.st_size != bytes) {
        sprintf(MsgStr, "ERROR: mesh grid file %s: %lld bytes, %zu expected",
                path, (long long) st.st_size, bytes);
        nll_puterr(MsgStr);
        close(fd);
        return (NULL);
    }
    void *ptr = mmap(NULL, bytes, PROT_READ, MAP_SHARED, fd, 0);
    close(fd);
    if (ptr == MAP_FAILED) {
        sprintf(MsgStr, "ERROR: mapping mesh grid file %s: %s", path, strerror(errno));
        nll_puterr(MsgStr);
        return (NULL);
    }

    if (num_mapped_files == cap_mapped_files) {
        int cap = cap_mapped_files > 0 ? 2 * cap_mapped_files : 256;
        MappedFile *grown = (MappedFile *) realloc(mapped_files, (size_t) cap * sizeof (MappedFile));
        if (grown == NULL) {
            munmap(ptr, bytes);
            nll_puterr("ERROR: allocating the mesh grid file register.");
            return (NULL);
        }
        mapped_files = grown;
        cap_mapped_files = cap;
    }
    strncpy(mapped_files[num_mapped_files].path, path, FILENAME_MAX - 1);
    mapped_files[num_mapped_files].path[FILENAME_MAX - 1] = '\0';
    mapped_files[num_mapped_files].ptr = ptr;
    mapped_files[num_mapped_files].bytes = bytes;
    num_mapped_files++;

    return (ptr);
}

/* ------------------------------------------------------------------------- */
/* meshes */

MeshDesc *MeshTime_LoadMesh(const char *path) {

    int n;
    for (n = 0; n < num_meshes; n++)
        if (strcmp(meshes[n]->path, path) == 0)
            return (meshes[n]);

    char fname[FILENAME_MAX];
    snprintf(fname, sizeof (fname), "%s.mesh.hdr", path);
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
        nll_puterr2("ERROR: opening mesh header file:", fname);
        return (NULL);
    }

    MeshDesc *mesh = (MeshDesc *) calloc(1, sizeof (MeshDesc));
    if (mesh == NULL) {
        fclose(fp);
        nll_puterr("ERROR: allocating a mesh description.");
        return (NULL);
    }
    strncpy(mesh->path, path, FILENAME_MAX - 1);

    char tag[MAXLINE];
    long long num_entries = 0;
    int istat = fscanf(fp, "%s %d %d", tag, &(mesh->num_nodes), &(mesh->num_tets));
    if (istat != 3 || strcmp(tag, "TET_MESH") != 0 || mesh->num_nodes < 4 || mesh->num_tets < 1) {
        nll_puterr2("ERROR: mesh header file: expected 'TET_MESH <num_nodes> <num_tets>':", fname);
        goto error;
    }
    istat = fscanf(fp, "%s %d %d %d %lf %lf %lf %lf %lf %lf %lld", tag,
            &(mesh->nbx), &(mesh->nby), &(mesh->nbz),
            &(mesh->bx0), &(mesh->by0), &(mesh->bz0),
            &(mesh->bdx), &(mesh->bdy), &(mesh->bdz), &num_entries);
    if (istat != 11 || strcmp(tag, "BINS") != 0 || mesh->nbx < 1 || mesh->nby < 1 || mesh->nbz < 1
            || mesh->bdx <= 0.0 || mesh->bdy <= 0.0 || mesh->bdz <= 0.0 || num_entries < 1) {
        nll_puterr2("ERROR: mesh header file: expected 'BINS <nbx> <nby> <nbz> <x0> <y0> <z0> <bdx> <bdy> <bdz> <num_entries>':", fname);
        goto error;
    }
    fclose(fp);
    fp = NULL;
    mesh->num_entries = (int64_t) num_entries;
    size_t num_bins = (size_t) mesh->nbx * (size_t) mesh->nby * (size_t) mesh->nbz;

    snprintf(fname, sizeof (fname), "%s.nodes.buf", path);
    if ((mesh->nodes = (const double *) map_file(fname, (size_t) mesh->num_nodes * 3 * sizeof (double))) == NULL)
        goto error;
    snprintf(fname, sizeof (fname), "%s.tets.buf", path);
    if ((mesh->tets = (const int32_t *) map_file(fname, (size_t) mesh->num_tets * 4 * sizeof (int32_t))) == NULL)
        goto error;
    snprintf(fname, sizeof (fname), "%s.bin_start.buf", path);
    if ((mesh->bin_start = (const int64_t *) map_file(fname, (num_bins + 1) * sizeof (int64_t))) == NULL)
        goto error;
    snprintf(fname, sizeof (fname), "%s.bin_tets.buf", path);
    if ((mesh->bin_tets = (const int32_t *) map_file(fname, (size_t) mesh->num_entries * sizeof (int32_t))) == NULL)
        goto error;
    if (mesh->bin_start[0] != 0 || mesh->bin_start[num_bins] != mesh->num_entries) {
        nll_puterr2("ERROR: mesh bin index is inconsistent with its header:", path);
        goto error;
    }

    if (num_meshes == cap_meshes) {
        int cap = cap_meshes > 0 ? 2 * cap_meshes : 4;
        MeshDesc **grown = (MeshDesc **) realloc(meshes, (size_t) cap * sizeof (MeshDesc *));
        if (grown == NULL) {
            nll_puterr("ERROR: allocating the mesh register.");
            goto error;
        }
        meshes = grown;
        cap_meshes = cap;
    }
    meshes[num_meshes++] = mesh;

    if (message_flag >= 2) {
        sprintf(MsgStr, "INFO: mesh %s: %d nodes, %d tetrahedra, %d x %d x %d bins",
                path, mesh->num_nodes, mesh->num_tets, mesh->nbx, mesh->nby, mesh->nbz);
        nll_putmsg(2, MsgStr);
    }

    return (mesh);

error:
    if (fp != NULL)
        fclose(fp);
    free(mesh); /* mappings stay registered: the files may be read again */
    return (NULL);
}

/* Barycentric coordinates of (x, y, z) in tetrahedron `tet`; 0 if degenerate. */
static int barycentric(const MeshDesc *mesh, int tet, double x, double y, double z, double lam[4]) {

    const int32_t *nd = mesh->tets + 4 * (size_t) tet;
    const double *v0 = mesh->nodes + 3 * (size_t) nd[0];
    const double *v1 = mesh->nodes + 3 * (size_t) nd[1];
    const double *v2 = mesh->nodes + 3 * (size_t) nd[2];
    const double *v3 = mesh->nodes + 3 * (size_t) nd[3];

    double a[3] = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
    double b[3] = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};
    double c[3] = {v3[0] - v0[0], v3[1] - v0[1], v3[2] - v0[2]};
    double p[3] = {x - v0[0], y - v0[1], z - v0[2]};

    /* Cramer's rule for [a b c] l = p */
    double bxc[3] = {b[1] * c[2] - b[2] * c[1], b[2] * c[0] - b[0] * c[2], b[0] * c[1] - b[1] * c[0]};
    double det = a[0] * bxc[0] + a[1] * bxc[1] + a[2] * bxc[2];
    if (det == 0.0 || !isfinite(det))
        return (0);
    double pxc[3] = {p[1] * c[2] - p[2] * c[1], p[2] * c[0] - p[0] * c[2], p[0] * c[1] - p[1] * c[0]};
    double bxp[3] = {b[1] * p[2] - b[2] * p[1], b[2] * p[0] - b[0] * p[2], b[0] * p[1] - b[1] * p[0]};
    double l1 = (p[0] * bxc[0] + p[1] * bxc[1] + p[2] * bxc[2]) / det;
    double l2 = (a[0] * pxc[0] + a[1] * pxc[1] + a[2] * pxc[2]) / det;
    double l3 = (a[0] * bxp[0] + a[1] * bxp[1] + a[2] * bxp[2]) / det;
    lam[0] = 1.0 - l1 - l2 - l3;
    lam[1] = l1;
    lam[2] = l2;
    lam[3] = l3;

    return (1);
}

int MeshTime_Locate(MeshDesc *mesh, double x, double y, double z, int *ptet, double w[4]) {

    int n;

    mesh->num_locate++;
    if (mesh->last_valid && x == mesh->last_x && y == mesh->last_y && z == mesh->last_z) {
        mesh->num_locate_cached++;
        *ptet = mesh->last_tet;
        for (n = 0; n < 4; n++)
            w[n] = mesh->last_w[n];
        return (mesh->last_tet >= 0);
    }

    int found = -1;
    double lam[4];
    double fx = (x - mesh->bx0) / mesh->bdx;
    double fy = (y - mesh->by0) / mesh->bdy;
    double fz = (z - mesh->bz0) / mesh->bdz;
    /* a point on the far face of the bin box belongs to the last bin */
    if (fx >= 0.0 && fy >= 0.0 && fz >= 0.0
            && fx <= (double) mesh->nbx && fy <= (double) mesh->nby && fz <= (double) mesh->nbz) {
        int ix = (int) fx, iy = (int) fy, iz = (int) fz;
        if (ix == mesh->nbx)
            ix--;
        if (iy == mesh->nby)
            iy--;
        if (iz == mesh->nbz)
            iz--;
        size_t bin = ((size_t) ix * (size_t) mesh->nby + (size_t) iy) * (size_t) mesh->nbz + (size_t) iz;
        int64_t k;
        for (k = mesh->bin_start[bin]; k < mesh->bin_start[bin + 1]; k++) {
            int tet = mesh->bin_tets[k];
            if (!barycentric(mesh, tet, x, y, z, lam))
                continue;
            if (lam[0] >= -MESH_BARY_TOL && lam[1] >= -MESH_BARY_TOL
                    && lam[2] >= -MESH_BARY_TOL && lam[3] >= -MESH_BARY_TOL) {
                /* clip and renormalise, as eikonal's SimplexMesh.locate */
                double sum = 0.0;
                for (n = 0; n < 4; n++) {
                    if (lam[n] < 0.0)
                        lam[n] = 0.0;
                    sum += lam[n];
                }
                for (n = 0; n < 4; n++)
                    lam[n] /= sum;
                found = tet;
                break;
            }
        }
    }

    mesh->last_valid = 1;
    mesh->last_x = x;
    mesh->last_y = y;
    mesh->last_z = z;
    mesh->last_tet = found;
    for (n = 0; n < 4; n++)
        mesh->last_w[n] = found >= 0 ? lam[n] : 0.0;
    *ptet = found;
    for (n = 0; n < 4; n++)
        w[n] = mesh->last_w[n];

    return (found >= 0);
}

/* ------------------------------------------------------------------------- */
/* time grids */

static int is_fixed(const GridDesc *pgrid, int node) {

    int lo = 0, hi = pgrid->mesh_num_fixed - 1;
    while (lo <= hi) {
        int mid = (lo + hi) / 2;
        int v = pgrid->mesh_fixed[mid];
        if (v == node)
            return (1);
        if (v < node)
            lo = mid + 1;
        else
            hi = mid - 1;
    }
    return (0);
}

int MeshTime_Open(GridDesc *pgrid, const char *fn_hdr, FILE *fp_hdr) {

    char line[MAXLINE_LONG];
    char tag[MAXLINE_LONG];
    char value[MAXLINE_LONG];
    char mesh_path[FILENAME_MAX] = "";
    double s0 = 0.0;
    int num_fixed = 0;
    int nline = 0;
    char label[MAXLINE_LONG];
    double sx = 0.0, sy = 0.0, sz = 0.0;
    int have_source = 0;

    pgrid->mesh = NULL;
    pgrid->mesh_values = NULL;
    pgrid->mesh_fixed = NULL;
    pgrid->mesh_num_fixed = 0;
    pgrid->mesh_s0 = 0.0;

    rewind(fp_hdr);
    while (fgets(line, MAXLINE_LONG, fp_hdr) != NULL) {
        nline++;
        if (nline == 2) {
            if (sscanf(line, "%s %lf %lf %lf", label, &sx, &sy, &sz) == 4)
                have_source = 1;
            continue;
        }
        if (sscanf(line, "%s %s", tag, value) != 2)
            continue;
        if (strcmp(tag, "MESH") == 0) {
            if (value[0] == '/') {
                strncpy(mesh_path, value, FILENAME_MAX - 1);
            } else {
                /* relative to the header's directory */
                char dir[FILENAME_MAX];
                strncpy(dir, fn_hdr, FILENAME_MAX - 1);
                dir[FILENAME_MAX - 1] = '\0';
                char *slash = strrchr(dir, '/');
                if (slash != NULL) {
                    *slash = '\0';
                    snprintf(mesh_path, sizeof (mesh_path), "%s/%s", dir, value);
                } else {
                    strncpy(mesh_path, value, FILENAME_MAX - 1);
                }
            }
        } else if (strcmp(tag, "SOURCE_SLOWNESS") == 0) {
            s0 = atof(value);
        } else if (strcmp(tag, "SOURCE_REGION") == 0) {
            num_fixed = atoi(value);
        }
    }
    rewind(fp_hdr);

    if (!have_source) {
        nll_puterr2("ERROR: mesh time grid header: no source (station) line:", fn_hdr);
        return (-1);
    }
    if (mesh_path[0] == '\0') {
        nll_puterr2("ERROR: mesh time grid header: no MESH line:", fn_hdr);
        return (-1);
    }
    if (s0 < 0.0 || num_fixed < 0) {
        nll_puterr2("ERROR: mesh time grid header: negative SOURCE_SLOWNESS or SOURCE_REGION:", fn_hdr);
        return (-1);
    }

    MeshDesc *mesh = MeshTime_LoadMesh(mesh_path);
    if (mesh == NULL)
        return (-1);

    /* <root>.hdr -> <root>.buf, <root>.fixed.buf */
    char root[FILENAME_MAX];
    strncpy(root, fn_hdr, FILENAME_MAX - 1);
    root[FILENAME_MAX - 1] = '\0';
    size_t len = strlen(root);
    if (len > 4 && strcmp(root + len - 4, ".hdr") == 0)
        root[len - 4] = '\0';
    char fname[FILENAME_MAX];
    snprintf(fname, sizeof (fname), "%s.buf", root);
    const GRID_FLOAT_TYPE *values =
            (const GRID_FLOAT_TYPE *) map_file(fname, (size_t) mesh->num_nodes * sizeof (GRID_FLOAT_TYPE));
    if (values == NULL)
        return (-1);
    const int *fixed = NULL;
    if (num_fixed > 0) {
        snprintf(fname, sizeof (fname), "%s.fixed.buf", root);
        if ((fixed = (const int *) map_file(fname, (size_t) num_fixed * sizeof (int32_t))) == NULL)
            return (-1);
    }

    pgrid->mesh = mesh;
    pgrid->mesh_values = values;
    pgrid->mesh_fixed = fixed;
    pgrid->mesh_num_fixed = num_fixed;
    pgrid->mesh_s0 = s0;
    pgrid->mesh_src_x = sx;
    pgrid->mesh_src_y = sy;
    pgrid->mesh_src_z = sz;

    return (0);
}

GRID_FLOAT_TYPE MeshTime_Read(GridDesc *pgrid, double x, double y, double z) {

    MeshDesc *mesh = (MeshDesc *) pgrid->mesh;
    if (mesh == NULL || pgrid->mesh_values == NULL)
        return (-VERY_LARGE_FLOAT);

    int tet;
    double w[4];
    if (!MeshTime_Locate(mesh, x, y, z, &tet, w))
        return (-VERY_LARGE_FLOAT);

    const int32_t *nd = mesh->tets + 4 * (size_t) tet;
    const GRID_FLOAT_TYPE *t = pgrid->mesh_values;
    int n;

    double s0 = pgrid->mesh_s0;
    if (s0 <= 0.0) {
        double time = 0.0;
        for (n = 0; n < 4; n++)
            time += w[n] * (double) t[nd[n]];
        return ((GRID_FLOAT_TYPE) (time > 0.0 ? time : 0.0));
    }

    /* factored read: s0 |x - x_s| plus the linear read of u = T - s0 |x_i - x_s| */
    double e[4];
    for (n = 0; n < 4; n++) {
        const double *v = mesh->nodes + 3 * (size_t) nd[n];
        double dx = v[0] - pgrid->mesh_src_x, dy = v[1] - pgrid->mesh_src_y, dz = v[2] - pgrid->mesh_src_z;
        e[n] = sqrt(dx * dx + dy * dy + dz * dz);
    }
    double dx = x - pgrid->mesh_src_x, dy = y - pgrid->mesh_src_y, dz = z - pgrid->mesh_src_z;
    double rho = sqrt(dx * dx + dy * dy + dz * dz);

    if (pgrid->mesh_num_fixed > 0 && is_fixed(pgrid, nd[0]) && is_fixed(pgrid, nd[1])
            && is_fixed(pgrid, nd[2]) && is_fixed(pgrid, nd[3])) {
        /* source region: |x - x_s| times the linear read of T_i / |x_i - x_s| */
        double ratio = 0.0;
        for (n = 0; n < 4; n++)
            ratio += w[n] * (e[n] > 0.0 ? (double) t[nd[n]] / e[n] : s0);
        double time = rho * ratio;
        return ((GRID_FLOAT_TYPE) (time > 0.0 ? time : 0.0));
    }

    double u = 0.0;
    for (n = 0; n < 4; n++)
        u += w[n] * ((double) t[nd[n]] - s0 * e[n]);
    /* a time is never negative: float rounding of u near the source can be (NLLoc rejects < 0) */
    double time = u + s0 * rho;
    return ((GRID_FLOAT_TYPE) (time > 0.0 ? time : 0.0));
}
