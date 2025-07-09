#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/trim.hpp"
#include "flexaid.h"
#include "boinc.h"

void read_emat(FA_Global* FA, char* emat_file)
{
    FILE* infile_ptr;
    char buffer[100];
    
    infile_ptr=NULL;
    if (!OpenFile_B(emat_file,"r",&infile_ptr)){
        fprintf(stderr,"ERROR: Could not read input file: %s\n", emat_file);
        Terminate(8);
    }
    
    std::vector<std::string> lines;
    std::string pairwiseline = "";
    
    // builds a string from a buffered line
    while(fgets(buffer,sizeof(buffer),infile_ptr) != NULL){
        pairwiseline += std::string(buffer);
        if(pairwiseline.find("\n") != std::string::npos){
            lines.push_back(pairwiseline);
            pairwiseline = "";
        }
    }
    
    printf("read %d lines in <%s>\n", (int)lines.size(), emat_file);
    
    float z = zero(1.0/2.0f, 1.0/2.0f, (float)(-(int)lines.size()));
    if(fabs(z - (float)((int)z)) > 0.001) {
        fprintf(stderr,"ERROR: Number of lines read in energy matrix file <%s> is incorrect (%d)\n", 
            emat_file, (int)lines.size());
        Terminate(12);
    }
    
    FA->ntypes = (int)(z + 0.001);
    printf("number of atom types: %d\n", FA->ntypes);
    
    FA->energy_matrix = (struct energy_matrix*)malloc(FA->ntypes*FA->ntypes*sizeof(struct energy_matrix));
    if(!FA->energy_matrix){
        fprintf(stderr,"ERROR: could not allocate memory for energy_matrix\n");
        Terminate(2);
    }
    
    for(int i=0;i<FA->ntypes;i++){
        for(int j=i;j<FA->ntypes;j++){
            std::string ori_line, line = *lines.begin();
            ori_line = line;
            lines.erase(lines.begin());
            
            if(line.find("=") != std::string::npos){
                line = line.substr(line.find("=")+1);
            }
            
            std::vector<std::string> values;
            boost::trim_if( line, boost::is_any_of(" \t\n") );
            boost::split( values, line, boost::is_any_of(" \t\n"), boost::token_compress_on );
            
            FA->energy_matrix[i*FA->ntypes+j].type1 = i+1;
            FA->energy_matrix[i*FA->ntypes+j].type2 = j+1;
            FA->energy_matrix[j*FA->ntypes+i].type1 = j+1;
            FA->energy_matrix[j*FA->ntypes+i].type2 = i+1;
            
            if(values.size() == 1){
                FA->energy_matrix[i*FA->ntypes+j].weight = 1;
                FA->energy_matrix[j*FA->ntypes+i].weight = 1;

                struct energy_values* weightval = (struct energy_values*)malloc(sizeof(struct energy_values));
                if(!weightval){
                    fprintf(stderr,"ERROR: could not allocate memory for weightval\n");
                    Terminate(2);
                }

                weightval->x = -1;
                weightval->y = atof((*values.begin()).c_str());
                weightval->next_value = NULL;
                
                FA->energy_matrix[i*FA->ntypes+j].energy_values = weightval;
                FA->energy_matrix[j*FA->ntypes+i].energy_values = weightval;
                
            }else if(values.size() % 2 == 0){
                FA->energy_matrix[i*FA->ntypes+j].weight = 0;
                FA->energy_matrix[j*FA->ntypes+i].weight = 0;
                
                struct energy_values* xyval_prev = NULL;
                
                for(std::vector<std::string>::iterator it=values.begin(); it!=values.end(); it+=2){
                    std::vector<std::string>::iterator xit = it;
                    std::vector<std::string>::iterator yit = it+1;
                    
                    struct energy_values* xyval = (struct energy_values*)malloc(sizeof(struct energy_values));
                    if(!xyval){
                        fprintf(stderr,"ERROR: could not allocate memory for xyval\n");
                        Terminate(2);
                    }
                    
                    xyval->x = atof((*xit).c_str());
                    xyval->y = atof((*yit).c_str());
                    xyval->next_value = NULL;
                    
					// multiply by a factor otherwise only conformer-driven
					//xyval->y *= 10.0;

					/*
					  printf("new entry[%d][%d]: x=%.3f y=%.3f\n", i+1, j+1,
					         xyval->x, xyval->y);
					*/

					// second or more xy values
                    if(xyval_prev != NULL)
                        xyval_prev->next_value = xyval;
                    else { 
                        FA->energy_matrix[i*FA->ntypes+j].energy_values = xyval;
                        FA->energy_matrix[j*FA->ntypes+i].energy_values = xyval;
                    }
                    
                    xyval_prev = xyval;
                }
            }else{
                fprintf(stderr,"ERROR: invalid number of xy-values for atom pairwise %d-%d\n", i, j);
                Terminate(12);
            }
        }
    }
    
    CloseFile_B(&infile_ptr,"r");
}