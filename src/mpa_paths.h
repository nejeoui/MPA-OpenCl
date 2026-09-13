#ifndef MPA_PATHS_H
#define MPA_PATHS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *mpaKernelDir(void)
{
    static char dir[512];
    static int  resolved = 0;
    static const char *cand[] = { ".", "src", "../src", "../../src" };
    const char *env;

    if (resolved) return dir;
    resolved = 1;

    env = getenv("MPA_KERNEL_DIR");
    if (env && *env) { snprintf(dir, sizeof dir, "%s", env); return dir; }

    for (size_t i = 0; i < sizeof cand / sizeof *cand; i++) {
        char  probe[512];
        FILE *f;
        snprintf(probe, sizeof probe, "%s/mpaKernel_32bits_opt.cl", cand[i]);
        f = fopen(probe, "rb");
        if (f) { fclose(f); snprintf(dir, sizeof dir, "%s", cand[i]); return dir; }
    }

    snprintf(dir, sizeof dir, "%s", ".");
    return dir;
}

static const char *mpaKernelPath(const char *name)
{
    static char buf[4][512];
    static int  slot = 0;
    char *p = buf[slot];
    slot = (slot + 1) & 3;

    if (strchr(name, '/')) { snprintf(p, 512, "%s", name); return p; }
    snprintf(p, 512, "%s/%s", mpaKernelDir(), name);
    return p;
}

#endif
