

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <limits>
#include <algorithm>
#include "VcfFileReader.h"
#include "VcfFileWriter.h"
#include "m3vcfHeader.h"
#include "m3vcfBlockHeader.h"
#include "m3vcfFileWriter.h"
#include "Unique.h"
#include "m3vcfBlock.h"
#define VCF 2
#define ZIPVCF 1
#define M3VCF 4
#define ZIPM3VCF 3
using namespace std;


struct convert_args_t
{

    //VCF File Variables
    IFILE vcfFileStream;
    VcfHeader myVcfHeader;
    
    // M3VCF File Variables
    IFILE m3vcfFileStream;
    m3vcfHeader out_hdr;
    m3vcfBlock myM3vcfBlock;
    bool VariantListOnly;

    std::vector<bool> sample_mask;

    //Common File and Argument Variables
    char **argv;
    const char *output_fname, *fname;
    int argc, output_type, record_cmd_line;
    char *sample_include_list, *sample_exlcude_list;
    int region_beg;
    int region_end;

};

static void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}
 

static std::string createCommandLine(convert_args_t &args, const char *optionName)
{
    string str = "##m3vcftools_Command=" + (string)optionName;
    for (int i=1; i<args.argc; i++){str+=" "; str += +args.argv[i];}
    time_t tt;
    str+="; Date=";
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    str+= asctime (timeinfo)  ;
    return  str.substr(0,str.size()-1);
}


static void CopyandPrintHeader(convert_args_t &args)
{
    args.myVcfHeader.reset();

    // Start from 1 since we don't need to convert the M3VCF format version
    for(int i=1; i<args.out_hdr.getNumMetaLines(); i++)
    {
        args.myVcfHeader.appendMetaLine(args.out_hdr.getMetaLine(i));
    }
    
    if(args.record_cmd_line==1)
        args.myVcfHeader.appendMetaLine(createCommandLine(args,"convert").c_str());

    //----------------------------------------------------------------//
    // Subset samples.
    std::string header_line;
    header_line.reserve(args.out_hdr.getHeaderLineSize());
    const char* s = args.out_hdr.getHeaderLine();
    const char*const e = s + args.out_hdr.getHeaderLineSize();

    std::size_t col = 0;
    const char* d = nullptr;
    while ((d = std::find(s, e, '\t')) != e)
    {
      assert(d);
      assert(col < args.sample_mask.size() + 9);
      if (col < 9 || args.sample_mask[col - 9])
      {
        if (col)
          header_line += '\t';
        header_line.append(s, d);
      }
      s = d + 1;
      ++col;
    }

    assert(col < args.sample_mask.size() + 9);
    if (col < 9 || args.sample_mask[col - 9])
    {
      if (col)
        header_line += '\t';
      header_line.append(s, e);
    }
    //----------------------------------------------------------------//

    args.myVcfHeader.addHeaderLine(header_line.c_str());
    args.myVcfHeader.write(args.vcfFileStream);

}


static void InitializeHaplotypeData(convert_args_t &args)
{ 
    args.vcfFileStream = ifopen(args.output_fname, "w", args.output_type==2? InputFile::UNCOMPRESSED : InputFile::GZIP);
    if ( !args.vcfFileStream ) error("[ERROR:] Failed to write in: %s\n", args.output_fname);
    
    args.m3vcfFileStream = ifopen(args.fname, "r"); 
    if ( !args.m3vcfFileStream ) error("[ERROR:] Failed to open: %s\n", args.fname);
    
    args.out_hdr.read(args.m3vcfFileStream);

    args.sample_mask.resize(args.out_hdr.getNumSamples(), true);
    if (args.sample_include_list)
    {
      std::unordered_set<std::string> include_set;
      std::ifstream ifs(args.sample_include_list);
      std::string sample_id;
      while (std::getline(ifs, sample_id))
        include_set.insert(sample_id);

      for (std::size_t i = 0; i < args.out_hdr.getNumSamples(); ++i)
      {
        if (include_set.find(args.out_hdr.getSampleName(i)) == include_set.end())
          args.sample_mask[i] = false;
      }
    }
}
           

static void Analyse(convert_args_t &args)
{
    int FirstBlock=0;
    while(args.myM3vcfBlock.read(args.m3vcfFileStream, args.out_hdr, args.VariantListOnly))
    {
        if (args.region_end < args.myM3vcfBlock.getStartBasePosition())
          break;

        if (args.region_beg <= args.myM3vcfBlock.getEndBasePosition())
        {
            for(int i=FirstBlock++==0?0:1; i<args.myM3vcfBlock.getNumMarkers(); i++)
            {
                int pos = args.myM3vcfBlock.getM3vcfRecord(i)->getBasePosition();
                if (pos >= args.region_beg && pos <= args.region_end)
                    args.myM3vcfBlock.writeVcfRecordGenotypes(args.vcfFileStream, args.sample_mask, i, args.VariantListOnly );
            }
        }
    }
       
    ifclose(args.m3vcfFileStream);
    ifclose(args.vcfFileStream);
}

static void convert_Data(convert_args_t &args)
{
    // Initialize and open file streams
    InitializeHaplotypeData(args);
    
    // Copy over header information
    CopyandPrintHeader(args);
   
    // Convert Main Data
    Analyse(args);
   
}


static void usage(convert_args_t &args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, " About:   Convert M3VCF file to VCF file  \n");
    fprintf(stderr, " Usage:   m3vcftools convert [options] <A.m3vcf.gz> \n");
    fprintf(stderr, "\n");
    fprintf(stderr, " Options:\n");
    fprintf(stderr, "       --no-version               do not append version and command line to the header\n");
    fprintf(stderr, "   -G, --drop-genotypes <file>            Write output to a file [standard output]\n");
    fprintf(stderr, "   -o, --output <file>            Write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type <z|v>        z: compressed VCF, v: uncompressed VCF [z] \n");
    fprintf(stderr, "   -S, --samples-file             file of samples to include \n"); // (or exclude with \"^\" prefix) \n");
    fprintf(stderr, "\n");
    exit(1);
}



int main_m3vcfconvert(int argc, char *argv[])
{
    int c;
    convert_args_t args;
    args.argc    = argc; args.argv = argv;
    args.output_fname = "-";
    args.output_type = ZIPVCF;
    args.record_cmd_line = 1;
    args.sample_include_list = NULL;
    args.sample_exlcude_list = NULL;
    args.VariantListOnly=false;
    args.region_beg = 0;
    args.region_end = std::numeric_limits<int>::max();
    
    static struct option loptions[] =
    {
        {"drop-genotypes",no_argument,NULL,'G'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"no-version",no_argument,NULL,8},
        {"region-beg",required_argument,NULL,'b'},
        {"region-end",required_argument,NULL,'e'},
        {"samples-file",required_argument,NULL,'S'},
        {NULL,0,NULL,0}
    };

    while ((c = getopt_long(argc, argv, "b:e:o:O:GS:",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'b': args.region_beg = atoi(optarg); break;
            case 'e': args.region_end = atoi(optarg); break;
            case 'o': args.output_fname = optarg; break;
            case 'S': args.sample_include_list = optarg; break;
            case 'G': args.VariantListOnly=true; break;
            case 'O':
                switch (optarg[0]) {
                    case 'z': args.output_type = ZIPVCF; break;
                    case 'v': args.output_type = VCF; break;
                    default: error("[ERROR:] The output type \"%s\" not recognised\n", optarg);
                };
                break;
            case  8 : args.record_cmd_line = 0; break;
            case '?': usage(args); break;
            default: error("[ERROR:] Unknown argument: %s\n", optarg);
        }
    }


    args.fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args.fname = "-";  // reading from stdin
        else usage(args);
    }
    else args.fname = argv[optind];
 
    if(args.sample_include_list!=NULL)
    {
        if(args.sample_include_list[0]=='^')
        {
            args.sample_exlcude_list=&args.sample_include_list[1];
            args.sample_include_list=NULL;
        }
    }
    convert_Data(args);
    return 0;
}

