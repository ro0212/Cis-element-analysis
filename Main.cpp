#include "impl.h"
using namespace std;

int main()
{
    string A, B;
    cout << "Enter first Cis-element: ";
    cin >> A;
    cout << "Enter second Cis-element: ";
    cin >> B;
    // A = "ACGT";
    // B= "AAAG";
    cout << "Enter Number of chromosomes in the genome: ";
    int chromosome_num;
    cin >> chromosome_num;

     cout << "Enter name of Folder containing chromosome data: ";
     string foldername;
     cin >> foldername;
    //foldername = "SweetPotato";

    ofstream outputFile("promoter.txt");
    vector<string> promoter;                                      // TO STORE ALL THE PROMOTERS FOUND
    vector<string> genomeSpacers, genomeSpacersReverse;           // Normal and Reverse Orientation
    vector<string> SpacersPerChromosome, SpacersPerChromosomeRev; // To store Spacers per chromosome;

    // -----------------------------
    vector<int> TotalSpacerFreqInGenome(21, 0);
    vector<int> TotalSpacerFreqInPromoter(21, 0);
    // ---------------------------------
    string genomeSequence;
    long int TotalSpacersInGnme = 0, TotalSpacersInPrmtr = 0;
    long long TotalNumPromoters = 0;

    int Start = 1;
    printf("Enter the starting Chromosome:\t");
   cin >> Start;

    for (int i = Start; i <= chromosome_num; ++i)
    {

        string filename = "/Users/rohitpagar/Desktop/SOP_BIO/Data/" + foldername + "/chr" + to_string(i) + ".fasta";
        ifstream inputFile(filename);
        if (!inputFile)
        {
            cerr << "Failed to open file: " << filename << endl;
            continue;
        }

        string ChromosomeSeq, line, partialSequence;

        int CountForPartialSeq = 0;
        int CountGenomeSpacers = 0;
        while (getline(inputFile, line))
        {
            if (!line.empty() && line[0] != '>')
            {
                ChromosomeSeq += line;
                partialSequence += line;
                CountForPartialSeq++;
            }
            if (CountForPartialSeq == 10000) // Partial Sequence is Now 10000 lines
            {
                findSpacers(partialSequence, SpacersPerChromosome, CountGenomeSpacers, A, B);
                CountForPartialSeq = 0;
                partialSequence.clear();
            }
        }
        // TO TAKE CARE OF REMAINING LINES
        if (CountForPartialSeq > 0)
        {
            findSpacers(partialSequence, SpacersPerChromosome, CountGenomeSpacers, A, B);
            CountForPartialSeq = 0;
            partialSequence.clear();
        }
        // ---- I HAVE SPACERS FOR iTH CHROMOSOME ------
        // --- MAINTAINING A FREQUENCY OF SPACERS FOR N .. 0-20
        vector<int> ChromosomeSpacerFreq(21, 0);
        for (const auto &spacer : SpacersPerChromosome)
        {
            ChromosomeSpacerFreq[spacer.size()]++;
        }

        for (int i = 0; i <= 20; ++i)
        {
            TotalSpacerFreqInGenome[i] += ChromosomeSpacerFreq[i];
        }

        // --------------------------
        // To input fre per chromosome in the file
        string S = "ChromosomeSpacerFreq" + to_string(i) + ".txt";
        ofstream ChromosomeSpacerFrqFile(S);
        for (int k = 0; k <= 20; ++k)
        {
            ChromosomeSpacerFrqFile << ChromosomeSpacerFreq[k] << endl;
        }

        // ---------------------------

        // ---------------

        vector<string> PromoterPerChr; // TO STORE PROMTOTERS FOUND

        // --------- FIND PROMOTERS PER CHROMOSOME -----
        findPromoter(ChromosomeSeq, PromoterPerChr);
        TotalNumPromoters += PromoterPerChr.size();
        printf("Size of one Promoter Seq : %zu\n", PromoterPerChr[1].size());

        // -------- FIND SPACERS FOR ALL PROMOTERS------
        int Count = 0;
        vector<string> SpacersInPromoter;
        for (const auto &it : PromoterPerChr)
        {
            // cout << it << "\n";
            findSpacers(it, SpacersInPromoter, Count, A, B);
        }

        // -----------  GET THE FREQ OF THE SPACERS ---
        vector<int> ChromosomeSpacerFreqInPromoter(21, 0);
        for (const auto &spacer : SpacersInPromoter)
        {
            // cout << spacer << "\n";
            ChromosomeSpacerFreqInPromoter[spacer.size()]++;
        }

        for (int i = 0; i <= 20; ++i)
        {
            TotalSpacerFreqInPromoter[i] += ChromosomeSpacerFreqInPromoter[i];
        }

        // ----------- INPUT THE FREQ IN TO FILE ----
        string St = "ChromosomeSpacerFreqPromoter" + to_string(i) + ".txt";
        ofstream ChromosomeSpacerFrqInPromoterFile(St);
        for (int k = 0; k <= 20; ++k)
        {
            ChromosomeSpacerFrqInPromoterFile << ChromosomeSpacerFreqInPromoter[k] << endl;
        }
        // ---------------------------

        // FREE MEM FOR SPACERS AND PROMOTERS
        PromoterPerChr.clear();
        SpacersInPromoter.clear();
        SpacersPerChromosome.clear();
        ChromosomeSeq.clear();

        // CLOSE THE DATA FILE
        inputFile.close();

        // -------------------------------------------------------------------

        // printf("Number of Genome Spacers in chromosome %d: %ld\n", i, genomeSpacers.size());
        // printf("Number of Genome Spacers in chromosome %d in Reverse-Orientation : %ld\n", i, genomeSpacersReverse.size());
        // cout << "Number of bases in chromosome " << i << ": " << genomeSequence.size() << " Nucleotides" << endl;
    }
    for (auto &count : TotalSpacerFreqInGenome)
        TotalSpacersInGnme += count;
    for (auto &count : TotalSpacerFreqInPromoter)
        TotalSpacersInPrmtr += count;

    printf("Total Number of Spacers in Entire Genome: %ld\n", TotalSpacersInGnme);
    printf("Total Number of Spacers in  Promoters: %ld\n", TotalSpacersInPrmtr);

    string Str = "TotalSpacerFreqinGenome.txt";
    ofstream TotalSpacerFreqinGenomeFile(Str);
    for (int i = 0; i <= 20; ++i)
    {
        TotalSpacerFreqinGenomeFile<< TotalSpacerFreqInGenome[i] << endl;
    }

    string PromSpacers = "TotalSpacersInPromoters.txt";
    ofstream SpacerfrqInPromters(PromSpacers);
    for (int i = 0; i <= 20; ++i)
    {
        SpacerfrqInPromters<< TotalSpacerFreqInPromoter[i] << endl;
    }

    printf("TOTAL NUMBER OF PROMOTERS FOUND: %lld\n", TotalNumPromoters);
   


    return 0;
}
