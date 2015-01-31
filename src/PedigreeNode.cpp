

#include "PedigreeNode.h"



Node::Node(char id)
{


    ID=id;
    father=NULL;
    mother=NULL;
    spouse=NULL;
    child=NULL;
    sex=info[id].sex;


};


Node::~Node()
{

};


Node* Node::insert(char id, map<char,data> info)
{

    ID=id;
    sex=info[id].sex;

    if(info[id].spouse=="")
        spouse=NULL;
    else
        {
            spouse=new Node(info[id].spouse);



            for(int i=0;i<(int)info[id].child.size();i++)
            {


                child.push_back(new Node(info[id].child[i]);
                child[i].insert(info[id].child[i]));
                if(sex=="M")
                {
                    child[i].father=this;
                    child[i].mother=spouse;
                }
                else
                {
                    child[i].mother=this;
                    child[i].father=spouse;
                }
            }


            spouse.child=child;
        }





    return this;


};



bool Node::search()
{
    return true;
};

