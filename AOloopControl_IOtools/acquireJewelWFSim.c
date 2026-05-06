/**
 * @file    acquireJewelWFSim.c
 * @brief   acquire and preprocess WFS image for Jewel/NRM data. 
 *
 */


#include <math.h>
#include <fftw3.h>

#include "CommandLineInterface/CLIcore.h"

#define FFTWOPTMODE FFTW_ESTIMATE

// Local variables pointers

static char *insname;
long fpi_insname;

static uint32_t *AOloopindex;
static long      fpi_AOloopindex;

static uint32_t *semindex;
static long      fpi_semindex;

static char *cropname;
long fpi_cropname;

static char *bandPSFlocs; //rock on
long fpi_bandPSFlocs;

static char *mfname; 
long fpi_mfname;

static CLICMDARGDEF farg[] =
{
    {
        CLIARG_STREAM,
        ".insname",
        "input stream name",
        "inV",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &insname,
        &fpi_insname
    },
    {
        CLIARG_UINT32,
        ".AOloopindex",
        "loop index",
        "0",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &AOloopindex,
        &fpi_AOloopindex
    },
    {
        CLIARG_UINT32,
        ".semindex",
        "input semaphore index",
        "1",
        CLIARG_HIDDEN_DEFAULT,
        (void **) &semindex,
        &fpi_semindex
    },
    {
        CLIARG_IMG,
        ".cropname",
        "Shm name of crop info to read from",
        "jewel_frame_szs",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &cropname,
        &fpi_cropname
    },
    {
        CLIARG_IMG,
        ".bandPSFlocs",
        "Sub-psf locations for this MBI band",
        "NAME_psfLocs",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &bandPSFlocs,
        &fpi_bandPSFlocs
    },
    {
        CLIARG_IMG,
        ".mfname",
        "Shm name of the matched filter to apply",
        "NAME_mf",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &mfname,
        &fpi_mfname
    },
};



// Optional custom configuration setup.
// Runs once at conf startup
//
static errno_t customCONFsetup()
{
    // if(data.fpsptr != NULL)
    // {
    //     data.fpsptr->parray[fpi_insname].fpflag |=
    //         FPFLAG_STREAM_RUN_REQUIRED | FPFLAG_CHECKSTREAM;


    //     data.fpsptr->parray[fpi_WFStaveragegain].fpflag  |= FPFLAG_WRITERUN;
    //     data.fpsptr->parray[fpi_WFStaveragemult].fpflag  |= FPFLAG_WRITERUN;
    //     data.fpsptr->parray[fpi_WFSnormfloor].fpflag     |= FPFLAG_WRITERUN;
    //     data.fpsptr->parray[fpi_compWFSsubdark].fpflag   |= FPFLAG_WRITERUN;
    //     data.fpsptr->parray[fpi_compWFSnormalize].fpflag |= FPFLAG_WRITERUN;
    //     data.fpsptr->parray[fpi_compWFSmask].fpflag      |= FPFLAG_WRITERUN;
    //     data.fpsptr->parray[fpi_compWFSrefsub].fpflag    |= FPFLAG_WRITERUN;
    //     data.fpsptr->parray[fpi_compWFSsigav].fpflag     |= FPFLAG_WRITERUN;
    //     data.fpsptr->parray[fpi_compWFSrefc].fpflag      |= FPFLAG_WRITERUN;
    //     data.fpsptr->parray[fpi_resetWFSrefc].fpflag     |= FPFLAG_WRITERUN;
    //     data.fpsptr->parray[fpi_WFSrefcgain].fpflag      |= FPFLAG_WRITERUN;
    //     data.fpsptr->parray[fpi_WFSrefcmult].fpflag      |= FPFLAG_WRITERUN;

    //     // reset WFS ave at startup
    //     data.fpsptr->parray[fpi_resetWFSrefc].fpflag |= FPFLAG_ONOFF;
    // }

    return RETURN_SUCCESS;
}

// Optional custom configuration checks.
// Runs at every configuration check loop iteration
//
static errno_t customCONFcheck()
{
    return RETURN_SUCCESS;
}

static CLICMDDATA CLIcmddata =
{
    "acquireJewelWFS", "acquire Jewel WFS signal", CLICMD_FIELDS_DEFAULTS
};

// detailed help
static errno_t help_function()
{
    return RETURN_SUCCESS;
}

/**
 * @brief  Pre-compute clamped crop bounds for an array of PSF centers.
 *
 * @param bounds     Pre-allocated output array of length n_psfs
 * @param locs       Flat UI64 array of PSF centers [cx0,cy0, cx1,cy1, ...]
 * @param n_psfs     Number of PSF centers
 * @param imcropsz   Requested square crop side length
 * @param sizexWFS   Source image width
 * @param sizeyWFS   Source image height
 */
 typedef struct
{
    long src_x0, src_y0;
    long crop_w, crop_h;
} CropBounds;

static void compute_crop_bounds(CropBounds      *bounds,
                                const uint64_t  *locs,
                                long             n_psfs,
                                long             imcropsz,
                                long             sizexWFS,
                                long             sizeyWFS)
{
    long half = imcropsz / 2;

    for(long p = 0; p < n_psfs; p++)
    {
        long cy = (long)locs[2*p]; // [row, col] = [y, x]
        long cx = (long)locs[2*p + 1];

        long x0 = cx - half;
        long y0 = cy - half;
        long x1 = x0 + imcropsz;
        long y1 = y0 + imcropsz;

        bounds[p].src_x0 = x0 < 0 ? 0 : x0;
        bounds[p].src_y0 = y0 < 0 ? 0 : y0;
        long src_x1 = x1 > sizexWFS ? sizexWFS : x1;
        long src_y1 = y1 > sizeyWFS ? sizeyWFS : y1;
        bounds[p].crop_w = src_x1 - bounds[p].src_x0;
        bounds[p].crop_h = src_y1 - bounds[p].src_y0;
    }
}

static errno_t compute_function()
{
    DEBUG_TRACE_FSTART();

    // read in crop and mf params from shm
    IMGID cropDim = mkIMGID_from_name(cropname); 
    resolveIMGID(&cropDim, ERRMODE_ABORT);
    long imcropsz = (long)cropDim.im->array.UI64[0];

    IMGID bandLocs = mkIMGID_from_name(bandPSFlocs);
    resolveIMGID(&bandLocs, ERRMODE_ABORT);
    long n_psfs = (long)bandLocs.md->size[1];  // px center per PSF

    IMGID mf = mkIMGID_from_name(mfname);
    resolveIMGID(&mf, ERRMODE_ABORT);
    long n_bl = (long)mf.md->size[0]; // number of baselines to samp from

    // connect to WFS image
    IMGID imgwfsim = stream_connect(insname);
    if(imgwfsim.ID == -1)
    {
        printf("ERROR: no WFS input\n");
        return RETURN_FAILURE;
    }

    uint32_t sizexWFS = imgwfsim.md->size[0];
    uint32_t sizeyWFS = imgwfsim.md->size[1];
    uint64_t sizeWFS  = sizexWFS * sizeyWFS;
    uint8_t  WFSatype = imgwfsim.md->datatype;

    // pre-calc crop coordinates
    long slice_npix  = imcropsz * imcropsz;  // pixels per crop slice

    CropBounds *bounds = malloc(sizeof(CropBounds) * n_psfs);
    if(bounds == NULL)
    {
        PRINT_ERROR("malloc() error");
        return RETURN_FAILURE;
    }
    compute_crop_bounds(bounds,
                        bandLocs.im->array.UI64,
                        n_psfs,
                        imcropsz,
                        (long)sizexWFS,
                        (long)sizeyWFS);


    // create/read images
    IMGID imgimWFS0;
    IMGID imgimWFS1;
    IMGID imgimWFS2;
    IMGID imgimWFS3;
    {
        char name[STRINGMAXLEN_IMGNAME];

        WRITE_IMAGENAME(name, "aol%u_imgimWFS0", *AOloopindex);
        imgimWFS0 = stream_connect_create_3Df32(name,
                                                imcropsz,
                                                imcropsz,
                                                n_psfs);

        WRITE_IMAGENAME(name, "aol%u_imgimWFS1", *AOloopindex);
        imgimWFS1 = stream_connect_create_3Df32(name,
                                                imcropsz,
                                                imcropsz,
                                                n_psfs);

        WRITE_IMAGENAME(name, "aol%u_imgimWFS2", *AOloopindex);
        imgimWFS2 = stream_connect_create_3Df32(name,
                                                n_bl*n_psfs,
                                                1,
                                                1);
        WRITE_IMAGENAME(name, "aol%u_imgimWFS3", *AOloopindex);
        imgimWFS3 = stream_connect_create_3Df32(name,
                                                imcropsz,
                                                imcropsz,
                                                n_psfs);
    }

    // Set-up FT plan and memory allocs
    fftwf_complex *fft_in  = fftwf_malloc(sizeof(fftwf_complex) * slice_npix * n_psfs);
    fftwf_complex *fft_out = fftwf_malloc(sizeof(fftwf_complex) * slice_npix * n_psfs);
    if(fft_in == NULL)
    {
        PRINT_ERROR("fftwf_malloc() error");
        return RETURN_FAILURE;
    }
    if(fft_out == NULL)
    {
        PRINT_ERROR("fftwf_malloc() error");
        return RETURN_FAILURE;
    }
    fftwf_plan plan = fftwf_plan_dft_2d(
                            imcropsz,                   // rows (slow axis)
                            imcropsz,                   // cols (fast axis)
                            fft_in,                     // dummy inptr
                            fft_out,                    // dummy outptr
                            FFTW_FORWARD,
                            FFTW_ESTIMATE 
                            );
    if(plan == NULL)
    {
        PRINT_ERROR("fftwf plan failed");
        return RETURN_FAILURE;
    }

    list_image_ID();

    int wfsim_semwaitindex =
        ImageStreamIO_getsemwaitindex(imgwfsim.im, *semindex);
    if(wfsim_semwaitindex > -1)
    {
        *semindex = wfsim_semwaitindex;
    }

    // initialize camera averaging arrays if not already done
    void *__restrict array_tmp;
    array_tmp = malloc(sizeof(float) * sizeWFS);
    if(array_tmp == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }
    float *__restrict arrayftmp = (float *) array_tmp;

    struct timespec time1, time2;
    long n_print_timings = 5000;

    INSERT_STD_PROCINFO_COMPUTEFUNC_START
    {
        // ================================================
        // COPY FRAME TO LOCAL MEMORY BUFFER + CAST TO FLOAT
        // ================================================
        int slice = 0;


        DEBUG_TRACEPOINT(" ");

        if(processinfo->loopcnt % n_print_timings == 0)
        {
            clock_gettime(CLOCK_MILK, &time1);
        }

        switch(WFSatype)
        {
            case _DATATYPE_FLOAT:
            {
                memcpy(arrayftmp, imgwfsim.im->array.F, sizeof(float) * sizeWFS);
            }
            break;

            case _DATATYPE_UINT16:
            {
                for(uint64_t ii = 0; ii < sizeWFS; ii++)
                    arrayftmp[ii] = (float)imgwfsim.im->array.UI16[ii];
            }
            break;

            case _DATATYPE_INT16:
            {
                for(uint64_t ii = 0; ii < sizeWFS; ii++)
                    arrayftmp[ii] = (float)imgwfsim.im->array.SI16[ii];
            }
            break;

            default:
                PRINT_ERROR("DATA TYPE NOT SUPPORTED");
                abort();
        }

        if(processinfo->loopcnt % n_print_timings == 0)
        {
            clock_gettime(CLOCK_MILK, &time2);
            printf("Pre-copy time: %f us\n", timespec_diff_double(time1, time2) * 1e6);
        }

        // ===================================================
        // wfsim -> imgimWFS0 (cropped interferogram images)
        // ===================================================
        DEBUG_TRACEPOINT(" ");

        if(processinfo->loopcnt % n_print_timings == 0)
        {
            clock_gettime(CLOCK_MILK, &time1);
        }

        imgimWFS0.md->write = 1;

        for(long p = 0; p < n_psfs; p++)
        {
            long src_x0    = bounds[p].src_x0;
            long src_y0    = bounds[p].src_y0;
            long crop_w    = bounds[p].crop_w;
            long crop_h    = bounds[p].crop_h;
            long slice_off = p * slice_npix;

            for(long row = 0; row < crop_h; row++)
            {
                float *src_row = arrayftmp + (src_y0 + row) * sizexWFS + src_x0;
                float *dst_row = imgimWFS0.im->array.F + slice_off + row * imcropsz;
                memcpy(dst_row, src_row, crop_w * sizeof(float));
            }
        }

        processinfo_update_output_stream(processinfo, imgimWFS0.ID);
        if(processinfo->loopcnt % n_print_timings == 0)
        {
            clock_gettime(CLOCK_MILK, &time2);
            printf("Test apply power 0.2 to imWFS0: %f us\n", timespec_diff_double(time1, time2) * 1e6);
        }


        DEBUG_TRACEPOINT(" ");

        // ===================================================
        // imgimWFS0 -> imgimWFS1 (phase)
        // ===================================================
        DEBUG_TRACEPOINT(" ");
        if(processinfo->loopcnt % n_print_timings == 0)
        {
            clock_gettime(CLOCK_MILK, &time1);
        }

        imgimWFS1.md->write = 1;

        // fft
        for(long p = 0; p < n_psfs; p++)
        {
            float *in = imgimWFS0.im->array.F + p * slice_npix;
            fftwf_complex *fin  = fft_in  + p * slice_npix;
            fftwf_complex *out = fft_out + p * slice_npix;

            // fftshift 
            long half = imcropsz / 2;

            for(long row = 0; row < imcropsz; row++)
            {
                long src_row = (row + half) % imcropsz;
                for(long col = 0; col < imcropsz; col++)
                {
                    long src_col = (col + half) % imcropsz;
                    fin[row * imcropsz + col][0] = in[src_row * imcropsz + src_col];
                    fin[row * imcropsz + col][1] = 0.0f;
                }
            }

            fftwf_execute_dft(plan, fin, out);

            // get phase and fftshift again
            float *phase_slice = imgimWFS1.im->array.F + p * slice_npix;

            for(long jj = 0; jj < imcropsz; jj++)
            {
                long src_row = (jj + half) % imcropsz;
                for(long ii = 0; ii < imcropsz; ii++)
                {
                    long src_col = (ii + half) % imcropsz;
                    float re = out[src_row * imcropsz + src_col][0];
                    float im = out[src_row * imcropsz + src_col][1];
                    phase_slice[jj * imcropsz + ii] = atan2f(im, re);
                }
            }


        }

        processinfo_update_output_stream(processinfo, imgimWFS1.ID);
        if(processinfo->loopcnt % n_print_timings == 0)
        {
            clock_gettime(CLOCK_MILK, &time2);
            printf("Test apply power 0.2 to imWFS1: %f us\n", timespec_diff_double(time1, time2) * 1e6);
        }

        DEBUG_TRACEPOINT(" ");

        // ===================================================
        // imgimWFS1 -> imgimWFS2 (sample, flatten, stack)
        // ===================================================
        DEBUG_TRACEPOINT(" ");
        if(processinfo->loopcnt % n_print_timings == 0)
        {
            clock_gettime(CLOCK_MILK, &time1);
        }

        imgimWFS2.md->write = 1;

        float         *dst      = imgimWFS2.im->array.F;
        const uint64_t *mf_coords = mf.im->array.UI64;

        for(long bl = 0; bl < n_bl; bl++)
        {
            for(long p = 0; p < n_psfs; p++)
            {
                long row = (long)mf_coords[0 * n_psfs * n_bl + p * n_bl + bl];
                long col = (long)mf_coords[1 * n_psfs * n_bl + p * n_bl + bl];;

                float sample = imgimWFS1.im->array.F[p * slice_npix
                                                     + row * imcropsz
                                                     + col];

                dst[bl * n_psfs + p] = sample;
            }
        }

        processinfo_update_output_stream(processinfo, imgimWFS2.ID);
        if(processinfo->loopcnt % n_print_timings == 0)
        {
            clock_gettime(CLOCK_MILK, &time2);
            printf("Test apply power 0.2 to imWFS2: %f us\n", timespec_diff_double(time1, time2) * 1e6);
        }

        DEBUG_TRACEPOINT(" ");

        // ===================================================
        // imgimWFS3 to check sampling 
        // ===================================================
        DEBUG_TRACEPOINT(" ");
        if(processinfo->loopcnt % n_print_timings == 0)
        {
            clock_gettime(CLOCK_MILK, &time1);
        }

        imgimWFS3.md->write = 1;

        // copy phase slices as base
        memcpy(imgimWFS3.im->array.F,
               imgimWFS1.im->array.F,
               sizeof(float) * slice_npix * n_psfs);

        for(long bl = 0; bl < n_bl; bl++)
        {
            for(long p = 0; p < n_psfs; p++)
            {
                long row = (long)mf_coords[0 * n_psfs * n_bl + p * n_bl + bl];
                long col = (long)mf_coords[1 * n_psfs * n_bl + p * n_bl + bl];

                imgimWFS3.im->array.F[p * slice_npix
                                      + row * imcropsz
                                      + col] = 4.0f;
            }
        }

        processinfo_update_output_stream(processinfo, imgimWFS3.ID);
        if(processinfo->loopcnt % n_print_timings == 0)
        {
            clock_gettime(CLOCK_MILK, &time2);
            printf("Test apply power 0.2 to imWFS3: %f us\n", timespec_diff_double(time1, time2) * 1e6);
        }

        DEBUG_TRACEPOINT(" ");
    }
    INSERT_STD_PROCINFO_COMPUTEFUNC_END

    // cleanup
    fftwf_destroy_plan(plan);
    fftwf_free(fft_in);
    fftwf_free(fft_out);
    free(bounds);
    free(array_tmp);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}



INSERT_STD_FPSCLIfunctions




// Register function in CLI
errno_t
CLIADDCMD_AOloopControl_IOtools__acquireJewelWFSim()
{

    CLIcmddata.FPS_customCONFsetup = customCONFsetup;
    CLIcmddata.FPS_customCONFcheck = customCONFcheck;
    INSERT_STD_CLIREGISTERFUNC

    return RETURN_SUCCESS;
}
