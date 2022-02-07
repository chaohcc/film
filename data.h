//
// Created by chaochao on 2021/12/17.
//

#ifndef BPLUSTREE_DATA_H
#define BPLUSTREE_DATA_H

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <fstream>
#include<ctime>
#include<iostream>
#include <string>
#include<fstream>
#include <vector>



using namespace std;
template<typename type_key>

/// load the 'dset' dataset, the number of loaded keys is 'numkeys'

vector<type_key> loaddata(string dset, size_t numkeys) {
    string filepath;
    std::vector<type_key> keys (numkeys);
    if (!dset.compare("books")) {
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/dataset/books.txt";

    } else if (!dset.compare("random0.5")) {
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/dataset/random0.5_83886080.txt";
    } else if (!dset.compare("wiki_ts")) {
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/dataset/wiki_ts_100000000.txt";

    } else if (!dset.compare("astro_ra")) {
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/dataset/astro_ra_18_21_89523861.txt";

    } else {
        printf("data name is wrong");
        exit(0);
    }

    ifstream fin(filepath);
    string s;
    double a;
    for (unsigned int i = 0; i < numkeys; i++) {
        getline(fin, s);
        istringstream istr1(s);
        istr1 >> a;
        keys[i] = a;
    }

    return keys;
}

std::vector<std::string> split(std::string str, std::string pattern)
{
    std::string::size_type pos;
    std::vector<std::string> result;
    str += pattern;//扩展字符串以方便操作
    int size = str.size();
    for (int i = 0; i < size; i++)
    {
        pos = str.find(pattern, i);
        if (pos < size)
        {
            std::string s = str.substr(i, pos - i);
            result.push_back(s);
            i = pos + pattern.size() - 1;
        }
    }
    return result;
}

template<typename key_type>
vector<key_type> loadpquery(string queryname, size_t queryn, int datanum, double a,string querytype) {
    string filepath;
    std::vector<key_type> queries(queryn);

    std::ostringstream ossb,ossc;
    ossb << queryn;
    ossc << datanum;
    string str_a;
    if (a==0)
        str_a = "0.0";
    else if (a==0.25)
        str_a = "0.25";
    else if (a==0.5)
        str_a = "0.5";
    else if (a==0.75)
        str_a = "0.75";
    else if (a==1.0)
        str_a = "1.0";
    else if (a==1.25)
        str_a = "1.25";
    else if (a==1.5)
        str_a = "1.5";
    else
        str_a = "1.25";


    filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/txtPQuery/workload" + queryname + "point_Zipf_" + str_a +
               "_" + ossb.str() + "_" + ossc.str() + ".txt";
    ifstream fin(filepath);
    std::string s;
    double tmp;
    for (unsigned int i = 0; i < queryn; i++) {
        getline(fin, s);
        istringstream istr1(s);
        istr1 >> tmp;
        queries[i] = tmp;
    }
    return queries;
}

template<typename key_type, typename singlkey_type>
vector<key_type> loadrquery(string queryname, size_t queryn, int datanum, double a,string querytype) {
    string filepath;
    std::vector<key_type> queries(queryn);

    std::ostringstream ossb,ossc;
    ossb << queryn;
    ossc << datanum;
    string str_a;
    if (a==0)
        str_a = "0.0";
    else if (a==0.25)
        str_a = "0.25";
    else if (a==0.5)
        str_a = "0.5";
    else if (a==0.75)
        str_a = "0.75";
    else if (a==1.0)
        str_a = "1.0";
    else if (a==1.25)
        str_a = "1.25";
    else if (a==1.5)
        str_a = "1.5";
    else
        str_a = "1.25";

    filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/txtRQuery/workload"+queryname+"range_Zipf_"+str_a+"_"+ossb.str()+"_"+ossc.str()+".txt";
    ifstream fin(filepath);
    std::string s;
    double tmp;
    vector<singlkey_type> tmps(2);
    for (unsigned int i = 0; i < queryn; i++) {
        getline(fin, s);

        vector<string> result = split(s," ");
        for (int k =0; k< result.size();k++){
            istringstream istr1(result[k]);

            istr1 >> tmp;
            tmps[k] = tmp;
            queries[i] = tmps;
        }
    }

    return queries;
}


#endif //BPLUSTREE_DATA_H
