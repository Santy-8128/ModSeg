#include<algorithm>
#include <iostream>
#include "Pedigree.h"
#include "string.h"
using namespace std;


void Pedigree::importGeno(int pos,vector<string> &data)
{

    string thisGeno=data[5];
    if(thisGeno.compare("P/N")==0)
        info[data[0]].geno="1";

    if(thisGeno.compare("P/P")==0)
        info[data[0]].geno="0";

    if(thisGeno.compare("N/N")==0)
        info[data[0]].geno="2";

    else
        info[data[0]].geno="";




}



void Pedigree::importPheno(int pos,vector<string> &data)
{


if(DataModel==1 || DataModel==3 || DataModel==5 || DataModel==7)
        {
            int temp=-1,temp2;


            bool temp3;
            istringstream buffer1(data[6]);
            if(data[6].compare("")!=0)
            {

                buffer1>>temp;
                if(temp>phenoMax)
                    temp=phenoMax;
                if(temp<phenoMin)
                    temp=phenoMin;
            }

            //cout<<temp<<" ";
            istringstream buffer2(data[7]);
            buffer2>>temp2;
            if(temp2==0)
            temp3=false;
            else
            temp3=true;
            //cout<<temp<<"\t"<<temp3<<endl;
            Phenotype tempPheno(temp,temp3);
            info[data[0]].pheno=tempPheno;

        }







}



void Pedigree::countSamples(Node* head,string from,string direction,int &count)
{
    count++;

    //cout<<head->ID<<"\t"<<from<<" "<<direction<<" "<<count<<endl;

    if(direction.compare("child")==0)
        for(int i=0;i<(int)head->spouse.size();i++)
        {
            countSamples(head->spouse[i],head->ID,"spouse",count);
            for(int j=0;j<(int)head->child[head->spouse[i]].size();j++)
                countSamples(head->child[head->spouse[i]][j],head->ID,"child",count);


        }



    if(direction.compare("spouse")==0)
    {
        for(int i=0;i<(int)head->spouse.size();i++)
            if(head->spouse[i]->ID.compare(from)!=0)
            {
                countSamples(head->spouse[i],head->ID,"spouse",count);
                for(int j=0;j<(int)head->child[head->spouse[i]].size();j++)
                    countSamples(head->child[head->spouse[i]][j],head->ID,"child",count);
            }
        if(head->founder==false)
                countSamples(head->father,head->ID,"parent",count);
    }


    if(direction.compare("parent")==0)
    {
        for(int i=0;i<(int)head->spouse.size();i++)
        {
            countSamples(head->spouse[i],head->ID,"spouse",count);
            for(int j=0;j<(int)head->child[head->spouse[i]].size();j++)
            {
                if(head->child[head->spouse[i]][j]->ID.compare(from)!=0)
                    countSamples(head->child[head->spouse[i]][j],head->ID,"child",count);
            }
        }

        if(head->founder==false)
            countSamples(head->father,head->ID,"parent",count);
    }

}



int Pedigree::NoSample()
{

    int count=1;

    if(first_founder==NULL)
        return 0;

    else
    {
        for(int i=0;i<(int)info[first_founder->ID].spouse.size();i++)
        {

            countSamples(first_founder->spouse[i],first_founder->ID,"spouse",count);
            for(int j=0;j<(int)first_founder->child[first_founder->spouse[i]].size();j++)
            {
                countSamples(first_founder->child[first_founder->spouse[i]][j],first_founder->ID,"child",count);
            }
        }
    }


    return count;

}
void Pedigree::addspouse(string s1,string s2)
{

    if(find( info[s1].spouse.begin(),  info[s1].spouse.end(), s2)==info[s1].spouse.end())
    {
        info[s1].spouse.push_back(s2);
        info[s2].spouse.push_back(s1);
    }

}



void Pedigree::creatInfo(vector<vector<string> > &ped)
{


    for(int i=0;i<(int)ped.size();i++)
    {

        info[ped[i][0]].self=ped[i][0];

        info[ped[i][0]].sex=ped[i][3];
        info[ped[i][0]].father=ped[i][1];
        info[ped[i][0]].mother=ped[i][2];
        info[ped[i][0]].founder=false;
        if(ped[i][1].compare("0")==0 && ped[i][2].compare("0")==0)
            info[ped[i][0]].founder=true;
        else
            addspouse(ped[i][1],ped[i][2]);

        info[ped[i][1]].child[ped[i][2]].push_back(ped[i][0]);
        info[ped[i][2]].child[ped[i][1]].push_back(ped[i][0]);


        importPheno(i,ped[i]);
        importGeno(i,ped[i]);
    }

}

void Pedigree::printInfo(vector<vector<string> > &ped)
{


    cout<<" \n ---- INFORMATION ABOUT PEDIGREE ----\n";
    cout<<"\nID\tFATHER\tMOTHER\tSEX\tFOUNDER\n";

    for(int i=0;i<(int)info.size();i++)
    {
        cout<<endl<<info[ped[i][0]].self<<"\t";

        cout<<info[ped[i][0]].father<<"\t";
        cout<<info[ped[i][0]].mother<<"\t";
        cout<<info[ped[i][0]].sex<<"\t";
        cout<<info[ped[i][0]].founder<<"\t";

        if(info[ped[i][0]].spouse.size()>0)
            cout<<"SPOUSE : \tCHILDREN :";

        for(int j=0;j<(int)info[ped[i][0]].spouse.size();j++)
        {
            cout<<"\n\t\t\t\t\t";
            cout<<info[ped[i][0]].spouse[j]<<"\t\t";
            for(int k=0;k<(int)info[ped[i][0]].child[info[ped[i][0]].spouse[j]].size();k++)
            {
                cout<<info[ped[i][0]].child[info[ped[i][0]].spouse[j]][k]<<"\t";
            }
        }
    }


    cout<<" \n ---- END OF PEDIGREE INFORMATION ---- \n\n";
}


string Pedigree::findRootIndex(vector<vector<string> > &ped)
{


    string root_index;
    for(int i=0;i<(int)ped.size();i++)
    {
        string id=ped[i][0];
        if(info[ped[i][0]].founder==true)
        {
            int j=0;
            while(j!=(int)info[id].spouse.size())
            {
                if(info[info[id].spouse[j]].founder==true)
                    {
                        root_index=ped[i][0];
                        first_founder=new Node(ped[i][0],info);
                        return root_index;


                    }
                j++;

            }
        }
    }

    return root_index;
}


Pedigree::Pedigree(vector<vector<string> > &ped,int Model,string A1,string A2,double af1,int censor,double pMin,double pMax)
{



    first_founder=NULL;


    DataModel=Model;
    Allele1=A1;
    Allele2=A2;
    censorIndicator=censor;
    phenoMin=pMin;
    phenoMax=pMax;


    creatInfo(ped);
    string root_index=findRootIndex(ped);
    //printInfo(ped);
    first_founder->pheno=info[first_founder->ID].pheno;
    first_founder->geno=info[first_founder->ID].geno;
    for(int i=0;i<(int)info[first_founder->ID].spouse.size();i++)
    {
        Node* thisSpouse=new Node(info[root_index].spouse[i],info);
        thisSpouse->insertNode(thisSpouse,info,"spouse",first_founder->ID);
        first_founder->spouse.push_back(thisSpouse);
        thisSpouse->spouse.push_back(first_founder);

        for(int j=0;j<(int)info[root_index].child[thisSpouse->ID].size();j++)
        {
            Node* thisChild=new Node(info[root_index].child[thisSpouse->ID][j],info);
            first_founder->child[thisSpouse].push_back(thisChild);
            thisChild->insertNode(thisChild,info,"child",first_founder->ID);
            if(first_founder->sex.compare("M")==0)
            {
                thisChild->father=first_founder;
                thisChild->mother=thisSpouse;
            }
            else
            {
                thisChild->father=thisSpouse;
                thisChild->mother=first_founder;
            }
            thisSpouse->child[first_founder].push_back(thisChild);

        }
    }

    //printPedigree();

};


Pedigree::~Pedigree()
{

};





Node::Node(string id, map<string,data> &info)
{


    ID=id;
    father=NULL;
    mother=NULL;
    sex=info[id].sex;
    founder=info[id].founder;
    pheno.ageOfOnset=-1.0;
};


Node::~Node()
{

};


void Pedigree::printPedigree()
{


    cout<<" \n ---- INFORMATION ABOUT PEDIGREE FROM STRUCTURE ----\n";

    if(first_founder!=NULL)
    {

        Node* head=first_founder;
        cout<<"\nID\t FATHER\t MOTHER\t SEX\t FOUNDER\n";

        cout<<"\n"<<head->ID<<"\t ";
        cout<<(head->father==NULL? "0": head->father->ID)<<"\t ";
        cout<<(head->mother==NULL? "0": head->mother->ID)<<"\t ";
        cout<<head->sex<<"\t "<<head->founder<< "\t\t";

        cout<<"SPOUSE : \tCHILDREN :";


        for(int i=0;i<(int)head->spouse.size();i++)
        {

            cout<<"\n\t\t\t\t\t\t";

            cout<<head->spouse[i]->ID<<"\t\t";
            for(int j=0;j<(int)head->child[head->spouse[i]].size();j++)
            {
                cout<<head->child[head->spouse[i]][j]->ID<<"\t";
            }


        }
        for(int i=0;i<(int)info[head->ID].spouse.size();i++)
        {

            printPedigree(head->spouse[i],head->ID,"spouse");
            for(int j=0;j<(int)head->child[head->spouse[i]].size();j++)
                {
                    printPedigree(head->child[head->spouse[i]][j],head->ID,"child");
                }

        }
    }


    cout<<" \n ---- END OF PEDIGREE INFORMATION FROM STRUCTURE ---- \n\n";


}
void Pedigree::printPedigree(Node* head,string from,string direction)
{


    cout<<"\n"<<head->ID<<"\t ";
    cout<<(head->father==NULL? "0": head->father->ID)<<"\t ";
    cout<<(head->mother==NULL? "0": head->mother->ID)<<"\t ";
    cout<<head->sex<<"\t "<<head->founder<< "\t\t";


    cout<<"\t"<<head->pheno.ageOfOnset<<" \t"<<head->pheno.censor;


    if(head->spouse.size()>0)
        cout<<"SPOUSE : \tCHILDREN :";


    for(int i=0;i<(int)head->spouse.size();i++)
    {
        cout<<"\n\t\t\t\t\t\t";

        cout<<head->spouse[i]->ID<<"\t\t";
        for(int j=0;j<(int)head->child[head->spouse[i]].size();j++)
        {
            cout<<head->child[head->spouse[i]][j]->ID<<"\t";
        }


    }

    if(direction.compare("child")==0)
        for(int i=0;i<(int)head->spouse.size();i++)
            {
                printPedigree(head->spouse[i],head->ID,"spouse");
                for(int j=0;j<(int)head->child[head->spouse[i]].size();j++)
                    printPedigree(head->child[head->spouse[i]][j],head->ID,"child");


            }



    if(direction.compare("spouse")==0)
    {

        for(int i=0;i<(int)head->spouse.size();i++)
        {
            if(head->spouse[i]->ID.compare(from)!=0)
            {
                printPedigree(head->spouse[i],head->ID,"spouse");
                for(int j=0;j<(int)head->child[head->spouse[i]].size();j++)
                    printPedigree(head->child[head->spouse[i]][j],head->ID,"child");
            }

            if(head->founder==false)
                printPedigree(head->father,head->ID,"parent");
        }
    }




    if(direction.compare("parent")==0)
    {
        for(int i=0;i<(int)head->spouse.size();i++)
        {
            printPedigree(head->spouse[i],head->ID,"spouse");
            for(int j=0;j<(int)head->child[head->spouse[i]].size();j++)
            {
                if(head->child[head->spouse[i]][j]->ID.compare(from)!=0)
                    printPedigree(head->child[head->spouse[i]][j],head->ID,"child");
            }
        }

        if(head->founder==false)
            printPedigree(head->father,head->ID,"parent");
    }


}


Node* Node::motherOfChild(Node* father, Node *child, map<string,data> &info)
{

    for(int i=0;i<(int)info[father->ID].spouse.size();i++)
        for(int j=0;j<(int)info[father->ID].child[info[father->ID].spouse[i]].size();j++)
            if(info[father->ID].child[info[father->ID].spouse[i]][j].compare(child->ID)==0)
                return father->spouse[i];

    return father->spouse[0];

}

void Node::addNodeChild(Node* current, Node* thisSpouse,int child_no, map<string,data> &info)
{

    Node* thisChild=new Node(info[current->ID].child[thisSpouse->ID][child_no],info);
    current->child[thisSpouse].push_back(thisChild);
    thisChild->insertNode(thisChild,info,"child",current->ID);
    if(current->sex.compare("M")==0)
    {
        thisChild->father=current;
        thisChild->mother=thisSpouse;
    }
    else
    {
        thisChild->father=thisSpouse;
        thisChild->mother=current;
    }

    thisSpouse->child[current].push_back(thisChild);


}


void Node::addNodeSpouseWithChild(Node* current, int spouse_no, map<string,data> &info,string leavechild_index)
{

    Node* thisSpouse=new Node(info[current->ID].spouse[spouse_no],info);
    thisSpouse=insertNode(thisSpouse,info,"spouse",current->ID );
    current->spouse.push_back(thisSpouse);
    thisSpouse->spouse.push_back(current);

    if(leavechild_index.compare("-1")==0)
        for(int j=0;j<(int)info[current->ID].child[info[current->ID].spouse[spouse_no]].size();j++)
            addNodeChild(current,thisSpouse,j,info);
    else
        for(int j=0;j<(int)info[current->ID].child[info[current->ID].spouse[spouse_no]].size();j++)
            if(info[current->ID].child[thisSpouse->ID][j].compare(leavechild_index)!=0)
                addNodeChild(current,thisSpouse,j,info);
}



void Node::addNodeFather(Node* current, map<string,data> &info)
{

    Node* thisFather=new Node(info[current->ID].father,info);
    thisFather=insertNode(thisFather,info,"parent",current->ID);
    Node* thisMother=thisFather->motherOfChild(thisFather,current,info);
    current->father=thisFather;
    current->mother=thisMother;
    thisFather->child[thisMother].push_back(current);
    thisMother->child[thisFather].push_back(current);

}

Node* Node::insertNode(Node* current, map<string,data> &info,string direction,string from)
{

    current->geno=info[current->ID].geno;

    current->pheno=info[current->ID].pheno;
    if(direction.compare("spouse")==0)
    {
        for(int i=0;i<(int)info[current->ID].spouse.size();i++)
            if(info[current->ID].spouse[i].compare(from)!=0)
                addNodeSpouseWithChild(current,i,info,"-1");
        if(info[current->ID].founder==false)
                addNodeFather(current,info);
    }

    if(direction.compare("child")==0)
        for(int i=0;i<(int)info[current->ID].spouse.size();i++)
            addNodeSpouseWithChild(current,i,info,"-1");

    if(direction.compare("parent")==0)
    {
        for(int i=0;i<(int)info[current->ID].spouse.size();i++)
            addNodeSpouseWithChild(current,i,info,from);
        if(info[current->ID].founder==false)
            addNodeFather(current,info);
    }

    return current;
}



