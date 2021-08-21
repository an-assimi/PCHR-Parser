
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <winsock.h>
#include <zlib.h>

#include "chr_config.h"
#include "chr_parser.h"
#include "xml.h"
#include "sbbt.h"
#include "chr_export.h"
#include "unzip.h"
#include "hmap_cell.h"


#define UMTS  1
#define GSM   2

/// Global variables

int log_domain=0;
int log_version=0;
int log_sub_version=0;
tree_t *cfg_fh, *cfg_lh, *cfg_chr[3]={NULL};
tree_t *log_fh, *log_lh, *log_chr;
char cfg_path[1024]="../config/";
char log_path[1024]="../data/";
char out_path[1024]="../output/";
char cfg_base_path[1024]={0};
expt_var_t *expt_var[3]={NULL};
hmaps_t *hmap_names;
hmap_t **table;
int num;
sbbt_node_t *sbbt=NULL;

hmapc_t *hmapc;
///==========================================================================
///==========================================================================
int load_cfg(void)
{
    char filename[256]= {0};
    char config_name[128]= {0};
    int code;

    if(log_domain==UMTS)
    {
        strcpy(filename,"Log_FileHead.txt");
        cfg_fh=read_cfg_file(cfg_path, filename, "Log_FileHead", 0);


        strcpy(filename,"Log_LogHead.txt");
        cfg_lh=read_cfg_file(cfg_path, filename, "Log_LogHead", 0);

        strcpy(filename,"Log_Performance User Log_6.txt");
        sscanf(filename, "Log_%[^_]_%d.txt", config_name, &code);
        cfg_chr[0]=read_cfg_file(cfg_path,filename, config_name, code);

        strcpy(filename,"Log_Performance Cell Log_7.txt");
        sscanf(filename, "Log_%[^_]_%d.txt", config_name, &code);
        cfg_chr[1]=read_cfg_file(cfg_path,filename, config_name, code);

        strcpy(filename,"Log_Performance Spec Log_41.txt");
        sscanf(filename, "Log_%[^_]_%d.txt", config_name, &code);
        cfg_chr[2]=read_cfg_file(cfg_path,filename, config_name, code);
    }
    else if(log_domain==GSM)
    {
        strcpy(filename,"Log_FileHead.txt");
        cfg_fh=read_cfg_file(cfg_path, filename, "Log_FileHead", 0);

        strcpy(filename,"Log_LogHead.txt");
        cfg_lh=read_cfg_file(cfg_path, filename, "Log_LogHead", 0);

        strcpy(filename,"Log_GSM CS XPU CHR LOG_16.txt");
        sscanf(filename, "Log_%[^_]_%d.txt", config_name, &code);
        cfg_chr[0]=read_cfg_file(cfg_path,filename, config_name, code);

        strcpy(filename,"Log_GSM DSPTC GCHR LOG_18.txt");
        sscanf(filename, "Log_%[^_]_%d.txt", config_name, &code);
        cfg_chr[1]=read_cfg_file(cfg_path,filename, config_name, code);
    }
    else
    {
        printf("Invalid log domain or non supported version..\n");
        return(1);
    }

    return(0);
}

int get_log_index(int ldoamin, int ltype)
{ int log_index=-1;
        if (ldoamin==UMTS)
        {
            switch(ltype)
            {
            case 6:
                log_index=0;
                break;
            case 7:
                log_index=1;
                break;
            case 41:
                log_index=2;
                break;
            default:
                log_index=-1;
                break;
            }
        }
        else if (ldoamin==GSM)
        {
            switch(ltype)
            {
            case 16:
                log_index=0;
                break;
            case 18:
                log_index=1;
                break;
            default:
                log_index=-1;
                break;
            }
        }
  return log_index;
}

///==========================================================================
///==========================================================================
int main(int argc, char *argv[])
{
    char chrfile[1024];
    char tmp_file[1024];
  //  FILE *fptr;
    uint8_t *buffer=NULL;
    int i;
    double offset, off;
    uint64_t log_id;
    FILE *res_file;//, *exp_file;


    if(argc<4)
    {
        printf("Incorrect number of parameters..\n");
        printf("COMMAND [UMTS|GSM] [V100|V900] [R19|R20|R21] [C00|C10] <chr_file.log>\n");
        return(1);
    }
    if(strcmp(argv[1],"UMTS")==0)
    {
        log_domain=1;
        strcat(cfg_path,"umts_chr/");
    }
    else if (strcmp(argv[1],"GSM")==0)
    {
        log_domain=2;
        strcat(cfg_path,"gsm_chr/");
    }
    strcpy(cfg_base_path, cfg_path);

    if(strcmp(argv[2],"V100")==0)
    {
        log_version=19;
        strcat(cfg_path,"V100");
    }
    else if (strcmp(argv[2],"V900")==0)
    {
        log_version=20;
        strcat(cfg_path,"V900");
    }

    if(strcmp(argv[3],"R19")==0)
    {
        log_version=19;
        strcat(cfg_path,"R019");
    }
    else if (strcmp(argv[3],"R20")==0)
    {
        log_version=20;
        strcat(cfg_path,"R020");
    }
    else if (strcmp(argv[3],"R21")==0)
    {
        log_version=21;
        strcat(cfg_path,"R021");
    }

    if(strcmp(argv[4],"C00")==0)
    {
        log_sub_version=0;
        strcat(cfg_path,"C00/");
    }
    else if (strcmp(argv[4],"C10")==0)
    {
        log_sub_version=10;
        strcat(cfg_path,"C10/");
    }

    strcpy(chrfile, log_path);
    strcat(chrfile,argv[5]);
    unsigned long fsz;
    /*
    fptr = fopen(chrfile,"rb");
    if(fptr == NULL)
    {
        printf("Error opening file %s\n",chrfile);
        return(1);
    }
    fseek(fptr, 0L, SEEK_END);
    fsz = ftell(fptr);
    fseek(fptr, 0L, SEEK_SET);
    uint8_t *buffer = (uint8_t*) malloc(fsz);
    long count =fread(buffer, 1, fsz, fptr);
    if(fsz!=count)
    {
        free(buffer0);
        printf("Problem reading chr_file: file_size=%ld  read=%ld\n",fsz, count);
        exit(0);
    }
    fclose(fptr);
    */
    //unsigned long  bufsize;
    buffer=decompress(chrfile, &fsz);



    printf("Loading configuration files....\n");
    cfg_fh=NULL; cfg_lh=NULL;
    cfg_chr[0]=cfg_chr[1]=cfg_chr[2]=NULL;
    load_cfg();

    printf("Loading enrichment database ...\n");
    read_db_file(cfg_base_path, "Logdb.xml");

    printf("Loading exported cdrs....\n");
    if(log_domain==UMTS)
    {
      expt_var[0]=(expt_var_t *) malloc(sizeof(expt_var_t));
      read_exp_file(cfg_base_path, "export_pchr_6.txt", expt_var[0]);
      //read_exp_file(cfg_base_path, "export_pchr_sho.txt", expt_var[0]);
    }
    else if (log_domain==GSM)
    {
      expt_var[0]=(expt_var_t *) malloc(sizeof(expt_var_t));
      read_exp_file(cfg_base_path, "export_gchr_16.txt", expt_var[0]);
    }

    strcpy(tmp_file,out_path);
    strcat(tmp_file, argv[5]);
    strcat(tmp_file, ".xml");
    res_file=fopen(tmp_file,"w");

    //strcpy(tmp_file,out_path);
    //strcat(tmp_file, argv[5]);
    //strcat(tmp_file, ".csv");
    //exp_file=fopen(tmp_file,"w");

    hmapc=(hmapc_t *) hmapc_new(256*256, 256);

    printf("Start parsing file %s\nPlease wait...\n", chrfile);
    log_fh=log_lh=log_chr=NULL;
    offset=0.0;

    /// Parse file header
    log_fh=tree_new(cfg_fh->name, cfg_fh->code, res_file, log_free, log_print2);
    pars_tree(buffer, cfg_fh->root, log_fh, log_fh->root, offset, &off);
    offset+=off;
    //printf("file header offset=%f   off=%f\n", offset, off);
    print_tree(log_fh,1);


    ///Parse logs
    int req_id=100000;
    log_id=0;

    //export_log_hdr(exp_file, expt_var[0]->names, expt_var[0]->num);

    while ((offset<fsz) &&(log_id<req_id))
    {
        log_id++;
        if(log_id %1000==0)
            printf("log=%lu\n",(unsigned long) log_id);

        ///Parse current log header
        log_lh=tree_new(cfg_lh->name, cfg_lh->code, res_file, log_free, log_print2);
        pars_tree(buffer, cfg_lh->root, log_lh, log_lh->root, offset, &off);
        offset+=off;

        ///Parse current log content
        node_t *n=tree_find(log_lh->root, "logtype");
        chr_log_t * vinfo=(chr_log_t *)n->info;
        int logtype=*(long *)vinfo->val;
        n=tree_find(log_lh->root, "loglength");
        vinfo=(chr_log_t *)n->info;
        int loglen=*(long *)vinfo->val;

        int log_index=get_log_index(log_domain,logtype);

        if (log_index>=0)
        {
            log_chr=tree_new(cfg_chr[log_index]->name, cfg_chr[log_index]->code, res_file, log_free, log_print2);
            pars_tree(buffer, cfg_chr[log_index]->root->childs, log_chr, log_chr->root, offset, &off);
            print_tree(log_chr,0);

            tree_free(log_lh);
            tree_free(log_chr);
        }
        if (off!=loglen)
            printf("error at %lu  logtype=%d len=%d   off=%0.1f\n",(unsigned long) log_id,logtype, loglen, off );

        offset+=loglen;

    }



    /// Free section
    tree_free(log_fh);

    free(buffer);

    tree_free(cfg_fh);
    tree_free(cfg_lh);

    hmapc_free(hmapc);

    int cfg_count=0;
    if (log_domain==UMTS)
    { cfg_count=3;
      for(i=0; i<cfg_count; i++) tree_free(cfg_chr[i]);
    }
    else if(log_domain==UMTS)
    {
      cfg_count=2;
      for(i=0; i<cfg_count; i++) tree_free(cfg_chr[i]);
    }


    hmap_str_free(hmap_names);
    for(int i=0; i<num; i++)
        hmap_free(table[i]);

    sbbt_free(sbbt);

    for(int i=0; i<3; i++)
    {
      if(expt_var[i])
      {
        if(expt_var[i]->vars) free(expt_var[i]->vars);
        if(expt_var[i]->names) free(expt_var[i]->names);
        free(expt_var[i]);
      }
    }

    fclose(res_file);
   // fclose(exp_file);


     printf("Program finished successfully.\n");
    return 0;

}
