#include<iostream>
#include<vector>
#include<string>
#include<algorithm>
#include<functional>

const int TGM_VALUE_BASE = 5;
const int NEIGHBORHOOD_RADIUS = 128;
const int MINIMUM_WEIGHT = 10;
const int BIG_NUMBER = 2000000000;


int nucleotide_value(char x);

class Strand{
private:
    bool complementary;
    long int id;
    std::string original_form;

    int trigram_counts[125];
    std::vector<int> key;

    bool trigram_comparator(const int& a, const int& b) const{
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


    Strand(int id, std::string& original_form, bool complementary) :
        id(id),
        original_form(original_form),
        complementary(complementary)
    {
        // create the key by parsing and sorting nucleotides
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
        auto cmp = std::bind(&Strand::trigram_comparator, this, std::placeholders::_1, std::placeholders::_2);
        std::sort(key.begin(), key.end(), cmp);
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
        // end of key; our key is a substring of other's key
        return true;
    }

    std::string get_original_form() const{
        return original_form;
    }

    int get_id() const{
        return id;
    }

    bool is_complementary() const{
        return complementary;
    }
};


struct ResultEdge{
    long int first_id;
    long int second_id;
    int weight;
    bool first_complementary;
    bool second_complementary;
    std::vector< std::pair<int, int> > alignment;
    friend std::ostream& operator<<(std::ostream& os, const ResultEdge& re);
};


std::ostream& operator<<(std::ostream& os, const ResultEdge& re){
    os << "FG " << re.first_id << "; FG " << re.second_id << "; " << re.weight << "; ";
    os << re.first_complementary << "; " << re.second_complementary << "; {";
    for (unsigned int i = 0; i < re.alignment.size(); ++i){
        if (re.alignment[i].first == -1){
            os << "_";
        }
        else{
            os << re.alignment[i].first;
        }
        os << "-";
        if (re.alignment[i].second == -1){
            os << "_";
        }
        else{
            os << re.alignment[i].second;
        }
        os << ";";
    }
    os << "}";
    return os;
}



int nucleotide_value(char x){
    switch (x){
    case 'A':
        return 0;
        break;
    case 'C':
        return 1;
        break;
    case 'G':
        return 2;
        break;
    case 'N':
        return 3;
        break;
    case 'T':
        return 4;
        break;
    default:
        return 3;
        break;
    }
}


std::string complementary(std::string strand){
    std::string result = "";
    for (unsigned int i = 0; i < strand.length(); ++i){
        switch (strand[i]){
        case 'A':
            result += 'T';
            break;
        case 'C':
            result += 'G';
            break;
        case 'G':
            result += 'C';
            break;
        case 'T':
            result += 'A';
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
            }
            strand = "";

            // get id - it starts at line[3]
            bool getting_number = false;
            id_string = "";
            for (unsigned int i = 3; i < line.length(); ++i){
                    id_string += line[i];
            }
            id = std::stoi(id_string);
        }
        else{
            // this line contains fragment of a strand.
            strand += line;
        }
    }
    // last strand is finished - we should push it into the vector, too.
    strands.push_back(Strand(id, strand, false));
    strands.push_back(Strand(id, complementary(strand), true));
}


// calculate 'weight' - distance between strands a and b;
// insert the edge to the given vector if it meets the criteria
void calculate_weight(const Strand& a, const Strand& b, std::vector<ResultEdge>& result_vector){
    const int GAP_PENALTY = -2;
    int l_a = a.get_original_form().length();
    int l_b = b.get_original_form().length();

    int** matrix = new int*[l_a + 1];
    int** aux_matrix = new int*[l_a + 1];
    for (int i = 0; i <= l_a; ++i){
        matrix[i] = new int[l_b + 1];
        aux_matrix[i] = new int[l_b + 1];
    }

    for (int i = 0; i <= l_a; ++i){
        for (int j = 0; j <= l_b; ++j){
            aux_matrix[i][j] = 0;
            if (i == 0 || j == 0){
                matrix[i][j] = 0;
            }
            else{
                if (a.get_original_form()[i] == b.get_original_form()[j]){
                    matrix[i][j] = 1;
                }
                else{
                    matrix[i][j] = -2;
                }
            }
        }
    }

    for (int i = 1; i <= l_a; ++i){
        for (int j = 1; j <= l_b; ++j){
            int way[4];
            way[1] = matrix[i - 1][j - 1] + matrix[i][j];
            way[2] = matrix[i][j - 1] + GAP_PENALTY;
            way[3] = matrix[i - 1][j] + GAP_PENALTY;
            int which_way = 1;
            for (int w = 2; w <= 3; ++w){
                if (way[w] > way[which_way]){
                    which_way = w;
                }
            }
            matrix[i][j] = way[which_way];
            aux_matrix[i][j] = which_way;
        }
    }

    int best_value = -BIG_NUMBER;
    int best_i = 0;
    int best_j = 0;
    for (int i = 1; i <= l_a; ++i){
        if (matrix[i][l_b] > best_value){
            best_value = matrix[i][l_b];
            best_i = i;
            best_j = l_b;
        }
    }
    for (int j = 1; j <= l_b; ++j){
        if (matrix[l_a][j] > best_value){
            best_value = matrix[l_a][j];
            best_i = l_a;
            best_j = j;
        }
    }

    if (best_value >= MINIMUM_WEIGHT){
        ResultEdge re;
        
        int current_i = best_i;
        int current_j = best_j;

        while (current_i != 0 && current_j != 0){
            switch (aux_matrix[current_i][current_j]){
            case 1:
                re.alignment.push_back(std::pair<int, int>(current_i, current_j));
                current_i--;
                current_j--;
                break;
            case 2:
                re.alignment.push_back(std::pair<int, int>(-1, current_j));
                current_j--;
                break;
            case 3:
                re.alignment.push_back(std::pair<int, int>(current_i, -1));
                current_i--;
                break;
            }
        }

        std::reverse(re.alignment.begin(), re.alignment.end());
        re.first_id = a.get_id();
        re.first_complementary = a.is_complementary();
        re.second_id = b.get_id();
        re.second_complementary = b.is_complementary();
        re.weight = best_value;

        // add the edge to the result vector
        result_vector.push_back(re);
    }

    // cleanup
    for (int i = 0; i <= l_a; ++i){
        delete[] matrix[i];
        delete[] aux_matrix[i];
    }
    delete[] matrix;
    delete[] aux_matrix;
}



std::vector<Strand> strands;
std::vector<ResultEdge> result_edges;

int main(){
    parse_input(strands);
    std::sort(strands.begin(), strands.end());

    for (int i = 0; i < strands.size(); ++i){
        int start = std::max(0, i - NEIGHBORHOOD_RADIUS);
        for (unsigned int j = start; (j <= start + 2 * NEIGHBORHOOD_RADIUS) && (j < strands.size()); ++j){
            if (j == i){
                continue;
            }
            calculate_weight(strands[i], strands[j], result_edges);
        }
    }

    for (unsigned int i = 0; i < result_edges.size(); ++i){
        std::cout << result_edges[i] << std::endl;
    }

    return 0;
}