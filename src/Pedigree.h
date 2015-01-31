
#ifndef __PEDIGREE_H
#define __PEDIGREE_H


#include <string>
#include <vector>
#include <sstream>
#include <map>
using namespace std;


class Phenotype
{
    public:

    int ageOfOnset;
    bool censor;
    bool status;
    Phenotype()
    {
        ageOfOnset=0;
        censor=false;
        status=false;
    };
    Phenotype(int x,bool y)
    {
        ageOfOnset=x;
        censor=y;
    };
    Phenotype(bool x)
    {
        status=x;
    };


};

class parents_type
{
    string father;
    string mother;
};

class data
{
    public:
    string self;
    string father;
    string mother;
    vector<string> spouse;
    string sex;
    bool founder;
    Phenotype pheno;
    string geno;
    map<string,vector<string> > child;
    friend class Pedigree;
    friend class Node;
};

class Node
{
    public:
    string ID;
    Node* father;
    Node* mother;
    vector<Node*> spouse;
    bool founder;
    string sex;
    map<Node*,vector<Node*> > child;
    Phenotype pheno;
    string geno;



    void addNodeChild(Node* current, Node* thisSpouse,int child_no, map<string,data> &info);
    void addNodeSpouseWithChild(Node* current, int spouse_no, map<string,data> &info,string leavechild_index);
    void addNodeFather(Node* current, map<string,data> &info);
    Node(string id, map<string,data> &info);    //constructor

    ~Node();  //destructor
    Node* insertNode(Node* current, map<string,data> &info,string direction,string from);

    Node* motherOfChild(Node* father, Node *child, map<string,data> &info);
    bool search();
    friend class Pedigree; // allow full access to myTree class
};




class Pedigree
{
    protected:


    public:
        map<string,data> info;

        Node *first_founder;
        int DataModel;
        // 0 = for case-control study, No Age Strat, No Sex Strat
        // 1 = for survival analysis study, No Age Strat, No Sex Strat
        // 2 = for case-control study, No Age Strat, ONLY Sex Strat
        // 3 = for survival analysis study, No Age Strat, ONLY Sex Strat
        // 5 = for survival analysis study, ONLY Age Strat, No Sex Strat
        // 0 = for case-control study, No Age Strat, No Sex Strat
        // 0 = for case-control study, No Age Strat, No Sex Strat

        string Allele1,Allele2;
        int censorIndicator;
        double phenoMin,phenoMax;


        //map<string,Phenotype> PhenoData;
        string findRootIndex(vector<vector<string> > &ped);
        void creatInfo(vector<vector<string> > &ped);

        void addspouse(string s1,string s2);
        void countSamples(Node* head,string from,string direction,int &count);
        void importPheno(int pos,vector<string> &data);

        void importGeno(int pos,vector<string> &data);
        int NoSample();
        void printInfo(vector<vector<string> > &ped);
        void printPedigree(Node* head,string from,string direction);
        void printPedigree();
        //Pedigree(vector<vector<string> > &ped,int PhenoNo,int GenoNo);// constructor
        Pedigree(vector<vector<string> > &ped,int Model,string A1,string A2,double af1,int censor,double pMin,double pMax);

        ~Pedigree(); // destrcutor
        void insert();
        bool search();
        //void printPedigree(Node* head);
    friend class Node;
};






#endif
