/**
 * @file    acquireJewelWFSim.c
 * @brief   acquire and preprocess WFS image for Jewel/NRM data. 
 *
 */


#include <math.h>

#include "CommandLineInterface/CLIcore.h"

// Local variables pointers

static char *insname;
long fpi_insname;

static uint32_t *AOloopindex;
static long      fpi_AOloopindex;

static uint32_t *semindex;
static long      fpi_semindex;

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


static errno_t compute_function()
{
    DEBUG_TRACE_FSTART();


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


    // create/read images
    IMGID imgimWFS0;
    {
        char name[STRINGMAXLEN_IMGNAME];

        WRITE_IMAGENAME(name, "aol%u_JewelTest", *AOloopindex);
        imgimWFS0 = stream_connect_create_2Df32(name, sizexWFS, sizeyWFS);
    }

    list_image_ID();

    int wfsim_semwaitindex =
        ImageStreamIO_getsemwaitindex(imgwfsim.im, *semindex);
    if(wfsim_semwaitindex > -1)
    {
        *semindex = wfsim_semwaitindex;
    }

    struct timespec time1, time2;
    long n_print_timings = 5000;

    INSERT_STD_PROCINFO_COMPUTEFUNC_START
    {
        // ===========================================
        // COPY FRAME TO LOCAL MEMORY BUFFER
        // ===========================================
        // void *__restrict array_tmp;
        // array_tmp = malloc(sizeof(float) * sizeWFS);
        // if(array_tmp == NULL)
        // {
        //     PRINT_ERROR("malloc returns NULL pointer");
        //     abort();
        // }
        // float *__restrict arrayftmp = (float *) array_tmp;
        // uint16_t *__restrict arrayutmp = (uint16_t *) array_tmp;
        // int16_t *__restrict arraystmp = (int16_t *) array_tmp;

        // int slice = 0;


        // DEBUG_TRACEPOINT(" ");

        // if(processinfo->loopcnt % n_print_timings == 0)
        // {
        //     clock_gettime(CLOCK_MILK, &time1);
        // }

        // void *ptrv = NULL;
        // switch(WFSatype)
        // {
        // case _DATATYPE_FLOAT:
        // case _DATATYPE_UINT16:
        // case _DATATYPE_INT16:
        // {
        //     int ts = ImageStreamIO_typesize(imgwfsim.md->datatype);
        //     ptrv = imgwfsim.im->array.raw + ts * slice * sizeWFS;
        //     memcpy(array_tmp, ptrv, ts * sizeWFS);
        // }
        // break;

        // default:
        //     PRINT_ERROR("DATA TYPE NOT SUPPORTED");
        //     abort();
        //     break;
        // }

        // if(processinfo->loopcnt % n_print_timings == 0)
        // {
        //     clock_gettime(CLOCK_MILK, &time2);
        //     printf("Pre-copy time: %f us\n", timespec_diff_double(time1, time2) * 1e6);
        // }

        // ===================================================
        // TEST FUNCTION
        // ===================================================
        DEBUG_TRACEPOINT(" ");

        if(processinfo->loopcnt % n_print_timings == 0)
        {
            clock_gettime(CLOCK_MILK, &time1);
        }

        imgimWFS0.md->write = 1;

        // dummy power law apply
        for(uint_fast64_t ii = 0; ii < sizeWFS; ii++)
        {
            imgimWFS0.im->array.F[ii] = 15;// powf(imgimWFS0.im->array.F[ii], 0.2f);
        }

        processinfo_update_output_stream(processinfo, imgimWFS0.ID);
        if(processinfo->loopcnt % n_print_timings == 0)
        {
            clock_gettime(CLOCK_MILK, &time2);
            printf("Test apply power 0.2 to imWFS0: %f us\n", timespec_diff_double(time1, time2) * 1e6);
        }


        DEBUG_TRACEPOINT(" ");

        // processinfo_WriteMessage_fmt(
        //     processinfo, "d%d n%d s%d a%d c%d",
        //     status_darksub,
        //     status_normalize,
        //     status_refsub,
        //     status_ave,
        //     status_wfsrefc
        // );
    }
    INSERT_STD_PROCINFO_COMPUTEFUNC_END


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
