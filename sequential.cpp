#include<iostream>
#include<vector>
#include<string>
#include<algorithm>
#include<functional>

const int TGM_VALUE_BASE = 5;

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

std::string parse_input(){
    std::string x;
    lol:
    std::getline(std::cin, x);
    if (!(x[0] == 'A' || x[0] == 'T' || x[0] == 'C' || x[0] == 'G' || x[0] == 'N')){
        goto lol;
    }
    return x;
}

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
    return result;
}

class Strand{
private:
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
    id(0){
    }

    // create the key by parsing and sorting nucleotides
    Strand(int id, std::string& original_form) :
    id(id),
    original_form(original_form)
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

std::vector<Strand> strands;


int main(){
    for (int i = 0; i < 50; ++i){
        std::string nct = parse_input();
        nct.pop_back();
        strands.push_back(Strand(0, nct));
        strands.push_back(Strand(0, complementary(nct)));
    }

    std::sort(strands.begin(), strands.end());
    std::cout << std::endl;
    for (int i = 0; i < 50; ++i){
        std::cout << strands[i].get_original_form() << std::endl;
    }
    return 0;
}