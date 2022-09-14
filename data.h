//
// Created by chaochao on 2021/12/17.
//

#ifndef DATA_H
#define DATA_H

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
#include <random>
#include "zipf.h"


template <class T>
bool load_binary_data(T data[], int length, const std::string& file_path) {
    std::ifstream is(file_path.c_str(), std::ios::binary | std::ios::in);
    if (!is.is_open()) {
        return false;
    }
    is.read(reinterpret_cast<char*>(data), std::streamsize(length * sizeof(T)));
    is.close();
    return true;
}

template <class T>
bool load_text_data(T array[], int length, const std::string& file_path) {
    std::ifstream is(file_path.c_str());
    if (!is.is_open()) {
        return false;
    }
    int i = 0;
    std::string str;
    double a;
    while (std::getline(is, str) && i < length) {
        std::istringstream ss(str);
        ss >> a;
        array[i] = a;
        ss >> array[i];
        i++;
    }
    is.close();

    return true;

}



/// load the 'dset' dataset, the number of loaded keys is 'numkeys'
template<typename type_key>
bool loaddata(type_key data[], std::string dset, size_t numkeys) {
    std::string filepath;
    //std::vector<type_key> keys (numkeys);

    if (!dset.compare("books")) {
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/dataset/books.txt";
        load_text_data(data, numkeys, filepath);
    } else if (!dset.compare("random0.5")) {
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/dataset/random0.5_83886080.txt";
        load_text_data(data, numkeys, filepath);
    } else if (!dset.compare("wiki_ts")) {
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/dataset/wiki_ts_100000000.txt";
        load_text_data(data, numkeys, filepath);
    } else if (!dset.compare("astro_ra")) {
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/dataset/astro_ra_18_21_89523861.txt";
        load_text_data(data, numkeys, filepath);
    }
    else if (!dset.compare("lognormal")){
        filepath = "/home/wamdm/chaohong/dataAchive/ALEX_1D/lognormal-190M.bin.data";
        load_binary_data(data, numkeys, filepath);
    }
    else if (!dset.compare("longlat")){
        filepath = "/home/wamdm/chaohong/dataAchive/ALEX_1D/longlat-200M.bin.data";
        load_binary_data(data, numkeys, filepath);
    }
    else if (!dset.compare("longitudes")){
        filepath = "/home/wamdm/chaohong/dataAchive/ALEX_1D/longitudes-200M.bin.data";
        load_binary_data(data, numkeys, filepath);
    }
    else if (!dset.compare("ycsb")){
        filepath = "/home/wamdm/chaohong/dataAchive/ALEX_1D/ycsb-200M.bin.data";
        load_binary_data(data, numkeys, filepath);
    }
    else {
        printf("data name is wrong");
        exit(0);
    }
    return true;
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
std::vector<key_type> loadpquery(std::string queryname, size_t queryn, int datanum, double a,std::string workload) {
    std::string filepath;
    std::vector<key_type> queries(queryn);

    std::ostringstream ossb,ossc;
    ossb << queryn;
    ossc << datanum;
    std::string str_a;
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


//    filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/txtPQuery/workload" + queryname + "point_Zipf_" + str_a +
//               "_" + ossb.str() + "_" + ossc.str() + ".txt";
    if (workload == "zipf"){
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/txtRQuery/workload"+queryname+"point_Zipf_"+str_a+"_"+ossb.str()+"_"+ossc.str()+".txt";
    }
    else if (workload == "random"){
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/txtPQuery/workload"+queryname+"point_Zipf_0.0"+"_"+ossb.str()+"_"+ossc.str()+".txt";
    }
    else if (workload == "zipfrandom"){
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/txtPQuery/shuffle"+queryname+"point_Zipf_0.75"+"_"+ossb.str()+"_"+ossc.str()+".txt";
    }
    else {  // workload == "hotspot"
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/txtPQuery/workload"+queryname+"point_hotspot_0.15"+"_"+ossb.str()+"_"+ossc.str()+".txt";
    }
    std::ifstream fin(filepath);
    std::string s;
    double tmp;
    for (unsigned int i = 0; i < queryn; i++) {
        getline(fin, s);
        std::istringstream istr1(s);
        istr1 >> tmp;
        queries[i] = tmp;
    }
    return queries;
}

template<typename key_type, typename singlkey_type>
std::vector<key_type> loadrquery(std::string queryname, size_t queryn, int datanum, double a,std::string workload) {
    std::string filepath;
    std::vector<key_type> queries(queryn);

    std::ostringstream ossb,ossc;
    ossb << queryn;
    ossc << datanum;
    std::string str_a;
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
    if (workload == "zipf"){
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/txtRQuery/workload"+queryname+"range_Zipf_"+str_a+"_"+ossb.str()+"_"+ossc.str()+".txt";
    }
    else if (workload == "random"){
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/txtRQuery/workload"+queryname+"range_Zipf_0.0"+"_"+ossb.str()+"_"+ossc.str()+".txt";
    }
    else if (workload == "zipfrandom"){
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/txtRQuery/shuffle"+queryname+"range_Zipf_0.75"+"_"+ossb.str()+"_"+ossc.str()+".txt";
    }
    else {  // workload == "hotspot"
        filepath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/txtRQuery/workload"+queryname+"range_hotspot_0.15"+"_"+ossb.str()+"_"+ossc.str()+".txt";
    }

    std::ifstream fin(filepath);
    std::string s;
    double tmp;
    std::vector<singlkey_type> tmps(2);
    for (unsigned int i = 0; i < queryn; i++) {
        getline(fin, s);

        std::vector<std::string> result = split(s," ");
        for (int k =0; k< result.size();k++){
            std::istringstream istr1(result[k]);

            istr1 >> tmp;
            tmps[k] = tmp;
            queries[i] = tmps;
        }
    }

    return queries;
}

template <class T>
T* get_search_keys(std::vector<T> array, int num_keys, int num_searches) {
    std::mt19937_64 gen(std::random_device{}());
    std::uniform_int_distribution<int> dis(0, num_keys - 1);
    auto* keys = new T[num_searches];
    for (int i = 0; i < num_searches; i++) {
        int pos = dis(gen);
        keys[i] = array[pos];
    }
    return keys;
}

template <class T>
T* get_search_keys_scrambledzipf(std::vector<T> array, unsigned long int num_keys, unsigned int num_searches) {
    auto* keys = new T[num_searches];
    ScrambledZipfianGenerator zipf_gen(num_keys);
    for (unsigned int i = 0; i < num_searches; i++) {
        unsigned int pos = zipf_gen.nextValue();
        keys[i] = array[pos];
    }
    return keys;
}

template <class T>
T* get_search_keys_zipf(std::vector<T> array, unsigned long int num_keys, unsigned long int num_searches, double zipf_factor) {

    auto* keys = new T[num_searches];
    zipf_distribution<long int, double > zipf_gen(num_keys,zipf_factor);
    std::mt19937 gen_;
    std::vector<T> zipf_res;
    for (unsigned int i = 0; i< num_searches; i ++){
        unsigned int pos = num_keys - zipf_gen.operator()(gen_);
        keys[i] = array[pos];
    }
    return keys;
}


template <class T>
T* get_search_keys_hotspot(std::vector<T> array, unsigned long int num_keys, unsigned long int num_searches, double hotratio=0.2,double accessratio = 0.9) {
    unsigned int hotspotlen = num_keys * hotratio;
    unsigned int hotqueryn = num_searches*accessratio;
    unsigned int randomqueryn = num_searches - hotqueryn;
    auto* keys = new T[num_searches];
    zipf_distribution<long int, double > zipf_gen(num_keys,0.75);
    std::mt19937_64 gen_random(std::random_device{}());
    std::mt19937 gen_;
    std::vector<T> hospot_res;
    unsigned int hot_start = num_keys - zipf_gen.operator()(gen_);
    while (num_keys-hot_start < hotspotlen)
        hot_start = hot_start-hotspotlen;
    std::uniform_int_distribution<int> dis1(hot_start, hot_start+hotspotlen);
    for (unsigned int i = 0; i< hotqueryn; i ++){
        unsigned int pos  = dis1(gen_random);
        keys[i] = array[pos];
    }

    std::uniform_int_distribution<int> dis2(0, num_keys - 1);
    for (unsigned int i = 0;i<randomqueryn;i++){
        unsigned int pos = dis2(gen_random);
        keys[hotqueryn+i] = array[pos];
    }
    return keys;
}



// generate the range queries
template <class T>
T** get_search_ranges(std::vector<T> array, int num_keys, int num_searches,int minlen = 0,int maxlen =100) {
    std::mt19937_64 gen(std::random_device{}());
    std::uniform_int_distribution<int> dis(0, num_keys - maxlen);
    std::mt19937_64 gen_random(std::random_device{}());
    std::uniform_int_distribution<int> disrange(minlen, maxlen);
    auto* ranges = new T*[num_searches];
    for (int i = 0; i < num_searches; i++) {
        unsigned int pos = dis(gen);
        ranges[i] = new T[2];
        ranges[i][0] = array[pos];
        ranges[i][1] = array[pos+disrange(gen_random)];
    }
    return ranges;
}


template <class T>
T** get_search_ranges_zipf(std::vector<T> array, unsigned long int num_keys, unsigned long int num_searches, double zipf_factor,unsigned int minlen = 0,unsigned int maxlen=100) {

    auto* ranges = new T*[num_searches];
    long int upper = num_keys-maxlen;
    zipf_distribution<long int, double > zipf_gen(upper,zipf_factor);
    std::mt19937 gen_;
    std::uniform_int_distribution<int> disrange(minlen, maxlen);
    std::mt19937_64 gen_random(std::random_device{}());
    for (int i = 0; i< num_searches; i ++){
        unsigned int pos = upper - zipf_gen.operator()(gen_);
        ranges[i] = new T[2];
        ranges[i][0] = array[pos];
        unsigned int pos2 = pos+disrange(gen_random);
//        if (pos2>num_keys-1){
//            std::cout<< "my Lord, i need You!"<< std::endl;
//        }
        ranges[i][1] = array[pos2];
    }
    return ranges;
}

template <class T>
T** get_search_ranges_scrambledzipf(std::vector<T> array, unsigned long int num_keys, unsigned int num_searches,int minlen = 0,int maxlen=100) {

    auto* ranges = new T*[num_searches];
    ScrambledZipfianGenerator zipf_gen(num_keys-maxlen);
    std::mt19937 gen_;
    std::uniform_int_distribution<int> disrange(minlen, maxlen);
    std::mt19937_64 gen_random(std::random_device{}());
    for (unsigned int i = 0; i< num_searches; i ++){
        unsigned int pos = zipf_gen.nextValue();
        ranges[i] = new T[2];
        ranges[i][0] = array[pos];
        ranges[i][1] = array[pos+disrange(gen_random)];
    }
    return ranges;
}


template <class T>
T** get_search_ranges_hotspot(std::vector<T> array, unsigned long int num_keys, unsigned long int num_searches, double hotratio=0.2,double accessratio = 0.9,int minlen = 0,int maxlen=100) {
    unsigned int hotspotlen = num_keys * hotratio;
    unsigned int hotqueryn = num_searches*accessratio;
    unsigned int randomqueryn = num_searches - hotqueryn;
    auto* ranges = new T*[num_searches];
    zipf_distribution<long int, double > zipf_gen(num_keys,0.5);
    std::mt19937_64 gen_random(std::random_device{}());
    std::mt19937 gen_;

    unsigned int hot_start = num_keys - zipf_gen.operator()(gen_);
    while (num_keys-hot_start< hotspotlen )
        hot_start = hot_start-hotspotlen;
    while (hot_start+hotspotlen+maxlen> num_keys )
        hot_start = hot_start-2*maxlen;
    std::uniform_int_distribution<int> dis1(hot_start, hot_start+hotspotlen);
    std::uniform_int_distribution<int> disrange(minlen, maxlen);
    for (unsigned int i = 0; i< hotqueryn; i ++){
        unsigned int pos  = dis1(gen_random);
        if (pos<hot_start )
            std::cout<<"i need You, my lovely Lord"<<std::endl;
        if (pos>hot_start+hotspotlen)
            std::cout<<"i need You, my lovely Lord"<<std::endl;
        ranges[i] = new T[2];
        ranges[i][0] = array[pos];
        ranges[i][1] = array[pos+disrange(gen_random)];

    }

    std::uniform_int_distribution<int> dis2(0, num_keys - maxlen);
    for (unsigned int i = 0;i<randomqueryn;i++){
        unsigned int pos = dis2(gen_random);
        ranges[hotqueryn+i] = new T[2];
        ranges[hotqueryn+i][0] = array[pos];
        ranges[hotqueryn+i][1] = array[pos+disrange(gen_random)];
    }
    return ranges;
}




#endif //DATA_H
