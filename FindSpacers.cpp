#include "impl.h"
using namespace std;

void findSpacers(const string &PromoterSequence, vector<string> &spacers, int &count, string a, string b)
{

    size_t l = PromoterSequence.find(a);
    if (l == string::npos)
    { // cout << "ACGT Absent" << endl;
        return;
    }

    // cout << "ACGT found" << endl;

    size_t h = PromoterSequence.find(b, l + 4);
    if (h == string::npos)
    {
        // cout << "AAAG Absent" << endl;
        return;
    }

    // cout << "AAAG found" << endl;

    size_t SpaceLen = h - 1 - (l + 4);
    // cout << "Spacer len : " << SpaceLen << endl;

    if (SpaceLen < 21)
    {
        count++;
        // cout << "Spacer sequence found" << endl;
        if (SpaceLen == 0)
        {
          //  outputFile << " " << endl;
            spacers.push_back("");
        }
        else
           // outputFile << PromoterSequence.substr(l + 4, SpaceLen) << endl;
           spacers.push_back(PromoterSequence.substr(l + 4, SpaceLen));
    }

    if (h + 1 < PromoterSequence.size())
    {
        findSpacers(PromoterSequence.substr(h + 1), spacers, count,a,b);
    }
}
