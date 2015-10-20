#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include<iostream>
#include<vector>
#include<string>
#include<algorithm>
#include<iomanip>

#include<thrust/device_vector.h>
#include<thrust/host_vector.h>

const int TGM_VALUE_BASE = 5;
const int TGM_VALUE_CB = TGM_VALUE_BASE * TGM_VALUE_BASE * TGM_VALUE_BASE;
const int NEIGHBORHOOD_RADIUS = 128;
const int MINIMUM_WEIGHT = 10;
const int BIG_NUMBER = 2000000000;
const int STRAND_LENGTH = 60;
const int GAP_PENALTY = -2;
const int BLOCKS_PER_INVOCATION = 16;


/*******************
 *  HOST FUNCTIONS *
 *******************/


int nucleotide_value(char x);

class Strand{
private:
    bool complementary;
    long int id;
    std::string original_form;

    std::vector<int> key;
    static int trigram_counts[TGM_VALUE_CB];

    static bool trigram_comparator(const int& a, const int& b) {
        if (trigram_counts[a] == trigram_counts[b]){
            return a < b;
        }
        else{
            return trigram_counts[a] > trigram_counts[b];
        }
    }


public:
    Strand() :
        id(0),
        original_form(""),
        complementary(false)
    {
    }


    Strand(int id, const std::string& original_form, bool complementary) :
        id(id),
        original_form(original_form),
        complementary(complementary)
    {
        // create the key by parsing and sorting nucleotides

        memset(trigram_counts, 0, sizeof(trigram_counts));
        const int tgm_sq = TGM_VALUE_BASE * TGM_VALUE_BASE;
        int current_trigram_value = nucleotide_value(original_form[0]) * TGM_VALUE_BASE +
            nucleotide_value(original_form[1]);
        for (unsigned int i = 2; i < original_form.length(); ++i){
            current_trigram_value %= tgm_sq;
            current_trigram_value *= TGM_VALUE_BASE;
            current_trigram_value += nucleotide_value(original_form[i]);
            trigram_counts[current_trigram_value]++;
            key.push_back(current_trigram_value);
        }
        std::sort(key.begin(), key.end(), trigram_comparator);
    }


    bool operator< (const Strand& other) const{
        for (unsigned int i = 0; i < key.size(); ++i){
            if (i == other.key.size()){
                // other's key is a substring of our key
                return false;
            }
            if (key[i] != other.key[i]){
                return key[i] < other.key[i];
            }
        }
        if (key.size() == other.key.size()){
            return false;
        }
        // end of key; our key is a substring of other's key
        return true;
    }

    std::string& get_original_form(){
        return original_form;
    }

    int get_id() const{
        return id;
    }

    bool is_complementary() const{
        return complementary;
    }
};

int Strand::trigram_counts[TGM_VALUE_CB];

struct ResultEdge{
    int weight;
};


int nucleotide_value(char x){
    switch (x){
    case 'a':
        return 0;
    case 'c':
        return 1;
    case 'g':
        return 2;
    case 't':
        return 3;
    default:
        return 4;
    }
}


std::string complementary(std::string strand){
    std::string result = "";
    for (unsigned int i = 0; i < strand.length(); ++i){
        switch (strand[i]){
        case 'a':
            result += 't';
            break;
        case 'c':
            result += 'g';
            break;
        case 'g':
            result += 'c';
            break;
        case 't':
            result += 'a';
            break;
        default:
            result += strand[i];
            break;
        }
    }
    // return reversed result
    return std::string(result.rbegin(), result.rend());
}


// get a vector of DNA fragments from standard input
void parse_input(std::vector<Strand>& strands){
    const char INFO_LINE_SYMBOL = '>';
    std::string line;
    std::string id_string;

    long int id = 0;
    std::string strand = "";

    while (std::getline(std::cin, line)){
        if (line[0] == INFO_LINE_SYMBOL){
            // previous strand finished; push it into the vector
            if (strand != ""){
                strands.push_back(Strand(id, strand, false));
                strands.push_back(Strand(id, complementary(strand), true));
                if (id % 100 == 0){
                    std::cerr << "Parsed: " << std::setfill('0') << std::setw(7) << id << std::endl;
                }
            }
            strand = "";

            // get id - it starts at line[3]
            id_string = "";
            for (unsigned int i = 3; i < line.length(); ++i){
                id_string += line[i];
            }
            id = std::atoi(id_string.c_str());
        }
        else{
            // this line contains fragment of a strand.
            strand += line;
        }
    }
    // last strand is finished - we should push it into the vector, too.
    strands.push_back(Strand(id, strand, false));
    strands.push_back(Strand(id, complementary(strand), true));
    std::cerr << "Parsed: " << std::setfill('0') << std::setw(7) << id << std::endl;
}

/********************
*  DEVICE FUNCTIONS *
*********************/

struct SimpleStrand{
    char nucleotides[STRAND_LENGTH];
};

__device__ int device_max(int a, int b){
    if (a > b){
        return a;
    }
    return b;
}

__global__ void calculate_values(int starting_n, SimpleStrand *sstrands, int n_strands, ResultEdge *resultEdges){

    int first_strand_id = blockIdx.x + starting_n;
    if (first_strand_id >= n_strands){
        resultEdges[blockIdx.x * NEIGHBORHOOD_RADIUS * 2 + threadIdx.x].weight = -1;
        return;
    }
    int start = first_strand_id - NEIGHBORHOOD_RADIUS;
    start -= device_max(0, NEIGHBORHOOD_RADIUS - n_strands + first_strand_id + 1);
    start = device_max(0, start);
    int second_strand_id = start + threadIdx.x;

    bool should_calculate = false;
    resultEdges[blockIdx.x * NEIGHBORHOOD_RADIUS * 2 + threadIdx.x].weight = -1;

    if (second_strand_id < first_strand_id){
        // we shouldn't compare two same strands twice
        int second_start = device_max(0, second_strand_id - NEIGHBORHOOD_RADIUS);
        int second_end = second_start + 2 * NEIGHBORHOOD_RADIUS;
        should_calculate = (second_end < first_strand_id);
    }
    else{
        second_strand_id += 1;
        should_calculate = (second_strand_id < n_strands);
    }

    if (should_calculate){
        SimpleStrand *first = sstrands + first_strand_id;
        SimpleStrand *second = sstrands + second_strand_id;

        // printf("%d %d\n-1- %c %c --\n-2- %c %c --\n", first_strand_id, second_strand_id, first->nucleotides[0], first->nucleotides[1], second->nucleotides[0], second->nucleotides[1]);

        int best_value = -BIG_NUMBER;
        int matrix_row[2][STRAND_LENGTH];
        int current_row = 0;
        int prev_row = 1;
        // fill first row
        for (int i = 0; i < STRAND_LENGTH; ++i){
            matrix_row[current_row][i] = 0;
        }
        for (int i = 1; i <= STRAND_LENGTH; ++i){
            prev_row = current_row;
            current_row += 1;
            current_row %= 2;
            matrix_row[current_row][0] = 0;
            for (int j = 1; j <= STRAND_LENGTH; ++j){
                if (first->nucleotides[i - 1] == second->nucleotides[j - 1]){
                    matrix_row[current_row][j] = 1;
                }
                else{
                    matrix_row[current_row][j] = -2;
                }
                int m = matrix_row[prev_row][j - 1] + matrix_row[current_row][j];
                m = device_max(m, matrix_row[prev_row][j] + GAP_PENALTY);
                m = device_max(m, matrix_row[current_row][j - 1] + GAP_PENALTY);
                matrix_row[current_row][j] = m;
                if (j == STRAND_LENGTH || i == STRAND_LENGTH){
                    if (m > best_value){
                        best_value = m;
                    }
                }
            }
        }
        resultEdges[blockIdx.x * NEIGHBORHOOD_RADIUS * 2 + threadIdx.x].weight = best_value;
    }
}


/****************
 *     MAIN     *
 ****************/


int main(){
    std::ios_base::sync_with_stdio(false);
    int counter = 0;
    std::cerr << "start" << std::endl;

    std::vector<Strand> strands;
    parse_input(strands);
    std::cerr << "parsed" << std::endl;

    std::sort(strands.begin(), strands.end());
    std::cerr << "sorted" << std::endl;

    int n_strands = strands.size();

    // create array of SimpleStrands on device
    thrust::host_vector<SimpleStrand> host_strands;
    thrust::device_vector<SimpleStrand> device_strands;
    SimpleStrand temp;
    for (unsigned int i = 0; i < n_strands; ++i){
        for (unsigned int j = 0; j < STRAND_LENGTH; ++j){
            temp.nucleotides[j] = strands[i].get_original_form()[j];
        }
        host_strands.push_back(temp);
    }
    device_strands = host_strands;
    SimpleStrand *simple_strands = thrust::raw_pointer_cast(device_strands.data());

    ResultEdge *dev_result_edges, *host_result_edges;
    host_result_edges = (ResultEdge*) malloc(BLOCKS_PER_INVOCATION * sizeof(ResultEdge) * 2 * NEIGHBORHOOD_RADIUS);
    cudaMalloc(&dev_result_edges, BLOCKS_PER_INVOCATION * sizeof(ResultEdge) * 2 * NEIGHBORHOOD_RADIUS);

    int last_i = 0;

    std::cerr << "Memory ready" << std::endl;

    // invoke kernel
    calculate_values <<< BLOCKS_PER_INVOCATION, 2 * NEIGHBORHOOD_RADIUS >>>(0, simple_strands, n_strands, dev_result_edges);
    cudaMemcpy(host_result_edges, dev_result_edges, BLOCKS_PER_INVOCATION * sizeof(ResultEdge) * 2 * NEIGHBORHOOD_RADIUS, cudaMemcpyDeviceToHost);
    for (int i = BLOCKS_PER_INVOCATION; i < n_strands; i += BLOCKS_PER_INVOCATION){
        std::cerr << "next loop  (" << std::setfill('0') << std::setw(7) << i << ")" << std::endl;
        // invoke next kernel
        calculate_values <<< BLOCKS_PER_INVOCATION, 2 * NEIGHBORHOOD_RADIUS >>>(i, simple_strands, n_strands, dev_result_edges);

        // output previous data
        for (int j = 0; j < BLOCKS_PER_INVOCATION; ++j){
            int index = j + last_i;
            int base = j * 2 * NEIGHBORHOOD_RADIUS; // where in the host_result_edges buffer does the data begin
            int start = index - NEIGHBORHOOD_RADIUS;
            start -= std::max(0, NEIGHBORHOOD_RADIUS - n_strands + index + 1);
            start = std::max(0, start);
            for (int k = 0; k < NEIGHBORHOOD_RADIUS * 2; ++k){
                int current = start + k;
                if (current >= index){
                    current += 1;
                }
                if (current >= n_strands){
                    break;
                }
                if (host_result_edges[base + k].weight < MINIMUM_WEIGHT){
                    continue;
                }
                counter++;
                std::cout << "FG " << strands[index].get_id() << "; FG " << strands[current].get_id() << "; ";
                std::cout << host_result_edges[base + k].weight << "; " << strands[index].is_complementary() << "; ";
                std::cout << strands[current].is_complementary() << std::endl;
            }

        }

        // sync with device, get next data
        last_i = i;
        cudaMemcpy(host_result_edges, dev_result_edges, BLOCKS_PER_INVOCATION * sizeof(ResultEdge) * 2 * NEIGHBORHOOD_RADIUS, cudaMemcpyDeviceToHost);
    }
    // output last batch of data
    for (int j = 0; j < BLOCKS_PER_INVOCATION; ++j){
        int index = j + last_i - BLOCKS_PER_INVOCATION;
        int base = j * 2 * NEIGHBORHOOD_RADIUS;
        int start = index - NEIGHBORHOOD_RADIUS;
        start -= std::max(0, NEIGHBORHOOD_RADIUS - n_strands + index + 1);
        start = std::max(0, start);
        for (int k = 0; k < NEIGHBORHOOD_RADIUS * 2; ++k){
            int current = start + k;
            if (current >= index){
                current += 1;
            }
            if (current >= n_strands){
                break;
            }
            if (host_result_edges[base + k].weight < MINIMUM_WEIGHT){
                continue;
            }
            counter++;
            std::cout << "FG " << strands[index].get_id() << "; FG " << strands[current].get_id() << "; ";
            std::cout << host_result_edges[base + k].weight << "; " << strands[index].is_complementary() << "; ";
            std::cout << strands[current].is_complementary() << std::endl;
        }

    }
    std::cerr << counter << std::endl;
    return 0;
}