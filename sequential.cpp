#include<iostream>
#include<vector>
#include<string>
#include<algorithm>
#include<functional>

const int TGM_VALUE_BASE = 5;
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
        id(0)
    {
    }

    // create the key by parsing and sorting nucleotides
    Strand(int id, std::string& original_form, bool complementary) :
        id(id),
        original_form(original_form),
        complementary(complementary)
    {
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
};


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



/*
std::string parse_input(){
    std::string x;
    lol:
    std::getline(std::cin, x);
    if (!(x[0] == 'A' || x[0] == 'T' || x[0] == 'C' || x[0] == 'G' || x[0] == 'N')){
        goto lol;
    }
    return x;
}
*/

std::string complementary(std::string strand){
    std::string result = "";
    for (int i = 0; i < strand.length(); ++i){
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


// add a strand and a complementary strand to a vector
void add_strands(long int id, std::string& strand, std::vector<Strand>& strands){
    strands.push_back(Strand(id, strand, false));
    strands.push_back(Strand(id, strand, true));
}


// get a vector of DNA fragments from standard input
void parse_input(std::vector<Strand>& strands){
    const char INFO_LINE_SYMBOL = '>';
    std::string line;
    std::string id_string;

    long int id = 0;
    std::string strand = "";

    while (std::getline(std::cin, line)){
        std::cout << "Read a line" << std::endl;
        if (line[0] == INFO_LINE_SYMBOL){
            // previous strand finished; push it into the vector
            if (strand != ""){
                strands.push_back(Strand(id, strand, false));
                strands.push_back(Strand(id, complementary(strand), true));
            }
            strand = "";

            std::cout << "Pushed strand into the vector" << std::endl;
            // get id - it will be located between first two '|' characters.
            bool getting_number = false;
            id_string = "";
            for (int i = 0; i < line.length(); ++i){
                if (line[i] == '|'){
                    if (getting_number){
                        break;
                    }
                    else{
                        getting_number = true;
                        continue;
                    }
                }
                if (getting_number){
                    id_string += line[i];
                }
            }
            std::cout << id_string << std::endl;
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


std::vector<Strand> strands;


int main(){
    parse_input(strands);

    std::sort(strands.begin(), strands.end());
    std::cout << std::endl;
    for (int i = 0; i < 100; ++i){
        std::cout << strands[i].get_original_form() << std::endl;
    }
    return 0;
}