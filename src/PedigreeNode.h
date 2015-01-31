#ifndef __PEDIGREE__NODE_H
#define __PEDIGREE__NODE_H

#include <map>
#include <vector>
#include "Pedigree.h"
#include "Individual.h"
using namespace std;


class Node
{
    char ID;
    Node* father;
    Node* mother;
    Node* spouse;
    char sex;
    vector<Node*> child;

    Node(char id);    //constructor

    ~Node();  //destructor

    Node* insert(char id, map<char,data> info);
    bool search();
    friend class Pedigree; // allow full access to myTree class
};


#endif




