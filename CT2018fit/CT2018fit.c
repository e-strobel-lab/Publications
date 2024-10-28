//
//  CT2018fit.c
//  
//
//  Created by Eric Strobel on 10/18/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <getopt.h>
#include <math.h>

#define MAX_LINE 4096   //maximum line length
#define MAX_FIELD 128   //maximum length of a field in input file

#define MAX_IPT 32      //maximum number of data columns in input file
#define COLUMN_HEADERS 4     //number of lines in input file

#define NAME_COL 0
#define YMIN_COL 1
#define YMAX_COL 2
#define EC50_COL 3

//parameters: struct to store sample parameters
typedef struct parameters {
    char name[MAX_FIELD];
    double ymin;
    double ymax;
    double EC50;
    double FL;
} parameters;

//structure declarations
int get_file(FILE **ifp, char *ipt);
int get_line(char *line, FILE *ifp);
int parse_input(FILE * ifp, parameters * p);
int compare_float(double a, double b, double precision);

int main (int argc, char *argv[])
{
    int i = 0;             //general purpose index
    int j = 0;             //general puprose index
    int ipt_provided = 0;  //tracks number of input files
    
    FILE * ifp; //input file pointer
    FILE * ofp; //output file pointer
    
    double min_conc =  0.000001;      //minimum ligand concentration
    double max_conc =  0.01;          //maximum ligand concentration
    double conc_step = min_conc/100;  //initial concentration step to use when calculating curve
    
    parameters p[MAX_IPT] = {{0}};  //array to store sample parameters
    
    int samples = 0; //number of samples
    
    char out_nm[MAX_LINE] = {0}; //output file name
    int out_nm_provided = 0;     //flag that output file name was provided
    
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"input",    required_argument, 0, 'i'},  //input values file
            {"out-name", required_argument, 0, 'o'},  //output file name
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "i:o:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
                
            case 'i': //get values file
                if (!ipt_provided) {                    //check that input file was not yet provided
                    get_file(&ifp, argv[optind-1]);     //set file pointer to input file
                    samples = parse_input(ifp, &p[0]);  //parse input and set sample count
                    ipt_provided++;                     //count input files provided
                } else {
                    printf("error - more than one input file was provided. aborting...\n");
                    abort();
                }
                break;

            case 'o': //output file name
                if (!out_nm_provided) {             //checkthat output file name was not yet provided
                    strcpy(out_nm, argv[optind-1]); //store output file name
                    out_nm_provided++;              //count output file names provided
                } else {
                    printf("error - more than one output file name was provided. aborting...\n");
                    abort();
                }
                break;
                
            default:
                printf("error: unrecognized option. Aborting program...\n");
                abort();
        }
    }
    
    if (optind < argc) {
        printf("\nnon-option ARGV-elements:\n");
        while (optind < argc)
            printf("\n%s \n", argv[optind++]);
        putchar('\n');
        printf("Aborting program.\n\n");
        abort();
    }
    /*********** end of option parsing ***********/
    
    if (!out_nm_provided) {        //if no output file name was provided
        strcpy(out_nm, "out.txt"); //set output file nameto default
    } else {
        strcat(out_nm, ".txt");    //otherwise, append suffix
    }
    
    //open output file
    if ((ofp = fopen(out_nm, "w")) == NULL) {
        printf("error - could not open output file. Aborting program...\n");
        abort();
    }
    
    //print header ine
    fprintf(ofp, "conc");                  //print concentration column header
    for (j = 0; j < samples; j++) {  //print sample name headers
        fprintf(ofp, "\t%s", p[j].name);
    }
    fprintf(ofp, "\n");
    
    double crnt_conc = min_conc; //set current concentration to the minimum concentration
    
    for (i = 0; !compare_float(crnt_conc, max_conc, 0.000000001); i++) { //until crnt_conc == max_conc
        
        fprintf(ofp, "%f", log10(crnt_conc));  //print current concentration
        for (j = 0; j < samples; j++) {  //calculate and print percent full length for every sample
            p[j].FL = ((p[j].ymax - p[j].ymin) * (crnt_conc/(p[j].EC50+crnt_conc))) + p[j].ymin;
            fprintf(ofp, "\t%f", p[j].FL);
        }
        fprintf(ofp, "\n");
        
        //if at next order of magnitude, increase conc_step ten-fold
        if (compare_float(crnt_conc, (conc_step * 100), 0.000000001)) {
            conc_step *= 10;
        }
        
        crnt_conc += conc_step; //increment crnt_conc by conc_step
    }
    
    /* close output file */
    if (fclose(ofp) == EOF) {
        printf("error - error occurred when closing output file. Aborting program...\n");
        abort();
    }
}


/* get_file: open file input and assign to file pointer */
int get_file(FILE **ifp, char *ipt)
{
    if ((*ifp = fopen(ipt, "r")) == NULL) {
        printf("get_file: error - could not open %s. Aborting program...\n", ipt);
        abort();
    }
    return 1;
}


/* get_line: get line from file, place into array, remove trailing newline, and return
 line length if successful */
int get_line(char *line, FILE *ifp)
{
    /* function gets line and replaces terminal newline with null character.
     there is no need for buffering in this case because lines that exceed
     MAX_LINE should not exist and if they do, they are an error.
     the only acceptable mode of failure is to reach the end of the file
     without getting any preceeding characters */
    
    int i = 0;
    char c = 0;
    
    for (i = 0; (c = fgetc(ifp)) != '\n' && c != EOF && c &&  i < MAX_LINE; i++) {
        line[i] = c;
    }
    if (c == '\n' && i != MAX_LINE) {
        line[i] = '\0';            //remove trailing newline
        return i;                //success
    } else if (c == EOF) {        //reached end of file
        if (i == 0) {            //EOF is expected at the start of a line
            return 0;
        } else {                //unexpected EOF
            printf("get_line: error -last line in file not terminated by newline\n");
            abort();
        }
    }else if (!c) {                //unexpected null character
        printf("get_line: error - unanticipated null character\n");
        abort();
    } else if (i == MAX_LINE) {    //unexpected long line
        printf("get_line: error - unanticipated long line\n");
        abort();
    } else {
        printf("get_line: error - reached unreachable code. figure out why. aborting...\n");
        abort();
    }
}

/* compare_float: test equality for floating point values a and b by
 testing whether the the absolute value of a-b is less than a user-provided
 precision setting. */
int compare_float(double a, double b, double precision)
{
    return ((fabs(a-b) < precision) ? 1 : 0);
}


/* parse_input: parse input values file and store sample parameters*/
int parse_input(FILE * ifp, parameters * p)
{
    int i = 0;  //general purpose index
    int j = 0;  //general purpose index
    int k = 0;  //general purpose index
    
    char line[MAX_LINE] = {0};                      //array to store file lines
    char tmp_str[COLUMN_HEADERS][MAX_FIELD] = {0};  //array to store strings temporarily
    char col_hdr[COLUMN_HEADERS][MAX_FIELD] = {0};  //array to store column headers
    
    int col = 0;           //column index
    int col_cnt = 0;       //column count
    int fnd_last_col = 0;  //flag that last column was found
    
    char xpctd_hdrs[4][5] = {"name", "ymin", "ymax", "EC50"}; //expected column headers
    
    get_line(line, ifp);  //get input file header line
    
    //get column headers and compare to expected column headers
    for (i = 0, j = 0; i < COLUMN_HEADERS; i++) {
        
        if (i > 0) { //if not first iteration,
            j++;     //increment j to bypass field terminator
        }
        
        //get the column header
        for (k = 0; line[j] != '\t' && line[j] != '\n' && line[j] != '\r' && line[j] && j < MAX_LINE && k < MAX_FIELD; j++) {
            col_hdr[i][k++] = line[j];
        }
        col_hdr[i][k] = '\0';
        
        if (line[j] != '\t' && line[j] != '\n' && line[j] != '\r' && line[j]) {
            printf("error - column field did not end on a field terminator. aborting...");
            abort();
        }
        
        if (strcmp(col_hdr[i], xpctd_hdrs[i])) {
            printf("error - col_id %s doesn't match expected header %s. aborting...\n", col_hdr[i], xpctd_hdrs[i]);
            abort();
        }
    }
    
    if (line[j] != '\n' && line[j] != '\r' && !line[j]) {
        printf("error - column header line has more fields than expected. the expected fields are: name, ymin, ymax, and EC50. aborting...");
        abort();
    }
    
    //get data line values and store in parameters struct
    for (i = 0; get_line(line, ifp); i++) {
        for (col = 0, j = 0; col < COLUMN_HEADERS; col++) {
            
            if (col > 0) { //if not first iteration of the current line,
                j++;       //increment j to bypass field terminator
            }
            
            //get current field string
            for (k = 0; line[j] != '\t' && line[j] != '\n' && line[j] != '\r' && line[j] && j < MAX_LINE && k < MAX_FIELD; j++) {
                if (col != NAME_COL   &&
                    (!isdigit(line[j]) &&
                     line[j] != '.'    &&
                     line[j] != '-'    &&
                     line[j] != 'e'    &&
                     line[j] != 'E')) {
                    printf("error - numerical data column value contains invalid character %c. aborting...\n", line[j]);
                    abort();
                } else {
                    tmp_str[col][k++] = line[j];
                }
            }
            tmp_str[col][k] = '\0';
                        
            if (line[j] != '\t' && line[j] != '\n' && line[j] != '\r' && line[j]) {
                printf("error - column field did not end on a field terminator. aborting...");
                abort();
            }
        }
        
        if (line[j] != '\n' && line[j] != '\r' && line[j]) {
            printf("error - data line has more fields than expected. the expected fields are: name, ymin, ymax, and EC50. aborting...");
            abort();
        }
        
        //store field values in parameters struct
        strcpy(p[i].name, tmp_str[NAME_COL]);
        p[i].ymin = strtod(tmp_str[YMIN_COL], NULL);
        p[i].ymax = strtod(tmp_str[YMAX_COL], NULL);
        p[i].EC50 = strtod(tmp_str[EC50_COL], NULL);
    }
    
    return i; //return data line count
}
