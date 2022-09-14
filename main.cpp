/*
 * F:\我的坚果云\2021-3-6hybridlearnedindex\ExperimentCC12\filminsert
 * path: /home/wamdm/chaohong/clionDir/insertfilm
 * support append-only setting
 * support query
 */
#include <iostream>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>
#include <iomanip>
#include <sstream>

#include <sys/time.h>

#include<ctime>

#include "data.h"
#include "film.h"


#define MAX_INT    (((unsigned int)(-1))>>1)



int test_interleave_insert_query(string filename, double memThreshold, long int datasize, long int numcolumn, int pagesize,unsigned int errbnd,double insert_frac){
    std::cout << "my Lord, i need You, i trust in You! test film interleave insert (append&out-of-order insert)and query" << std::endl;
    struct timeval t1, t2;
    double timeuse;

    unsigned long int numkey = ceil(double(datasize*1024*1024)/numcolumn/8);  //int datanum
    unsigned long actual_numkey = ceil(numkey*1.2);   // num_key and actual_numkey, the former is the at least number to be inserted to the index. the latter is the keys to guarantee the numkey
    auto keys = new key_type[actual_numkey];
    loaddata<key_type>(keys, filename,actual_numkey);
    std::vector<key_type> load_data(keys,keys+actual_numkey);
    double reserveMem = 10;
    memThreshold -= reserveMem;

    double zipf = 0.75;   //zipfian factor
    string workload = "zipf";

    std::cout << "my Lord, i need You, i trust in You!!!!!" << std::endl;
    gettimeofday(&t1, NULL);


    cout<<"the data set is "<<filename<<",  the data size is "<<datasize<<"M"<< "  numkey is "<<numkey<<endl;
    cout<<"the number of keys is "<<load_data.size()<<", the record size is "<<numcolumn<<endl;
    cout<<"the page size is "<<pagesize*8<<",  the available memory is "<< memThreshold<< endl;


    unsigned long init_num_key = ceil(numkey*0.75);

    unsigned long batch_size = 100000;  // ceil(numkey*0.1)
    double out_of_order_fracs[] = {0.5};  //0.1,0.2,0.3,0.4,0.5,    0.01,0.001,0.0001,0.6,0.7,0.8
    unsigned int thresholds[] = {5};
    for (int i = 0; i  < (end(out_of_order_fracs)-begin(out_of_order_fracs));i++){
        double out_of_order_frac = out_of_order_fracs[i];
        for (int k = 0; k < end(thresholds)-begin(thresholds);k++){
            unsigned int threshold = thresholds[k] * 10000;
            filminsert::test_interleave_insert_query(errbnd,numkey,pagesize,filename,memThreshold,
                                                     reserveMem,numcolumn,load_data,actual_numkey,datasize,zipf,init_num_key, insert_frac, out_of_order_frac, batch_size,workload,threshold);
        }

    }


    gettimeofday(&t2, NULL);
    timeuse = (t2.tv_sec - t1.tv_sec) + (double) (t2.tv_usec - t1.tv_usec) / 1000000.0;

    cout << "able o_direct disk access time = " << timeuse << endl;  //输出时间（单位：ｓ）
    return 0;
}


int test_out_of_order_insertion(string filename, double memThreshold, long int datasize, long int numcolumn, int pagesize,unsigned short int errbnd){
    double insert_frac[] = {0.5};  //0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9

    for (int i = 0; i < (end(insert_frac)-begin(insert_frac));i++) {
        double insert_ratio = insert_frac[i];   //// the total fraction of append and out-of-order inserts
        test_interleave_insert_query(filename, memThreshold, datasize, numcolumn, pagesize,errbnd,insert_ratio);
    }
    return 0;
}

int test_query_baseline(string filename, double memThreshold, long int datasize, long int numcolumn, int pagesize,unsigned short int errbnd){
    std::cout << "my Lord, i need You, i trust in You! test the baselines" << std::endl;
    struct timeval t1, t2;
    double timeuse;

    unsigned long int numkey = ceil(double(static_cast<unsigned long >(datasize)*1024*1024)/static_cast<unsigned long >(numcolumn)/8);
    auto keys = new key_type[numkey];
    loaddata<key_type>(keys, filename,numkey);
    std::vector<key_type> load_data(keys,keys+numkey);
    delete []keys;
    double reserveMem = 5;
    memThreshold = memThreshold-reserveMem;
    unsigned int queryn = 100000;  // 100000
    string workload = "zipfrandom";  //"random"
    double zipf = 0.75;   //zipfian factor
    vector<key_type> load_pointquery ;
    vector<vector<key_type>> load_rangequery;


    std::cout << "my Lord, i need You, i trust in You!!!!!" << std::endl;
    gettimeofday(&t1, NULL);
    double zipfs[] = {0.75}; //,0.25,0.5,0.75,1.0,1.25,1.5
    string workloads[] = {"zipf"};  ///"zipf","random","zipfrandom","hotspot"
    gettimeofday(&t1, NULL);
    for (int i = 0; i < (end(workloads)-begin(workloads));i++){
        gettimeofday(&t1, NULL);
        cout<<"the data set is "<<filename<<",  the data size is "<<datasize<<"M"<< "  numkey is "<<numkey<<endl;
        cout<<"the number of keys is "<<numkey<<", the record size is "<<numcolumn<<endl;
        cout<<"the page size is "<<pagesize*8<<",  the available memory is "<< memThreshold<< endl;

        workload = workloads[i];   //zipfian factor

//        load_rangequery = loadrquery<vector<key_type>,key_type>(filename,queryn,numkey,zipf,workload);
        load_pointquery = loadpquery<key_type>(filename,queryn,numkey,zipf,workload);

        filminsert::test_filmadaquery(errbnd,numkey,pagesize,filename,memThreshold,reserveMem,numcolumn,load_data,load_pointquery,load_rangequery,queryn,numkey,datasize,zipf,workload);

        gettimeofday(&t2, NULL);
        timeuse = (t2.tv_sec - t1.tv_sec) + (double) (t2.tv_usec - t1.tv_usec) / 1000000.0;
        cout << "able o_direct disk access time = " << timeuse << endl;  //输出时间（单位：ｓ）
    }

    vector <key_type>().swap(load_data);
    vector <key_type>().swap(load_pointquery);
    vector <vector<key_type>>().swap(load_rangequery);




    return 0;
}


int test_interleave_baselines(string filename, double memThreshold, long int datasize, long int numcolumn, int pagesize,unsigned short int errbnd){
    std::cout << "my Lord, i need You, i trust in You! test the interleave insert&query on baselines" << std::endl;
    struct timeval t1, t2;
    double timeuse;

    unsigned long int numkey = ceil(double(datasize*1024*1024)/numcolumn/8);  //int datanum,
    unsigned long actual_numkey = ceil(numkey*1.2);   // num_key and actual_numkey, the former is the at least number to be inserted to the index. the latter is the keys to guarantee the numkey
    auto keys = new key_type[actual_numkey];
    loaddata<key_type>(keys, filename,actual_numkey);
    std::vector<key_type> load_data(keys,keys+actual_numkey);
    delete []keys;
    double reserveMem = 10;
    memThreshold = memThreshold-reserveMem;
//    string workload = "random";  //"random"
    string workloads[] = {"zipf","random","zipfrandom","hotspot"};  ///"zipf","random","zipfrandom","hotspot"
    double zipf = 0.75;   //zipfian factor

    std::cout << "my Lord, i need You, i trust in You!!!!!" << std::endl;
    gettimeofday(&t1, NULL);
    unsigned long init_num_key = ceil(numkey*0.75);

    unsigned int batch_size = 100000;   //100000
//    double insertfracs[] = {0.5};
    gettimeofday(&t1, NULL);
    for (int i = 0; i  < (end(workloads)-begin(workloads));i++){
//        double insertfrac = insertfracs[i];
        double insertfrac = 0.5;
        string workload = workloads[i];
        filminsert::test_interleave_insert_query_baseline(errbnd,numkey,pagesize,filename,memThreshold,
                                                 reserveMem,numcolumn,load_data,actual_numkey,datasize,zipf,init_num_key, insertfrac, 0.0, batch_size,workload);

        gettimeofday(&t1, NULL);
        cout<<"the data set is "<<filename<<",  the data size is "<<datasize<<"M"<< "  numkey is "<<numkey<<endl;
        cout<<"the number of keys is "<<numkey<<", the record size is "<<numcolumn<<endl;
        cout<<"the page size is "<<pagesize*8<<",  the available memory is "<< memThreshold<< endl;



        gettimeofday(&t2, NULL);
        timeuse = (t2.tv_sec - t1.tv_sec) + (double) (t2.tv_usec - t1.tv_usec) / 1000000.0;
        cout << "able o_direct disk access time = " << timeuse << endl;  //输出时间（单位：ｓ）

    }

    vector <key_type>().swap(load_data);

    return 0;
}



int main() {

    std::cout << "my Lord, i need You, i trust in You! start the experiment " << std::endl;
    unsigned int errbnd = 16;
    double memThreshold = 512;
    int pagesize = 1024*64/8;
    string filename = "wiki_ts";  //queryname
    long int numcolumn = 16;
    long int datasize = 4096;
    double zipf = 0.5;
//    test_lru_overhead(datasize,memThreshold);

    datasize = 2048*2;
    memThreshold = 1024*2;
    filename = "books";
    //    filename = "wiki_ts";
//    test_appending(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);
//    filename = "astro_ra";


    errbnd = 16 ;

    // out-of-order insertion 的代码片段
    datasize = 1024*4;
    memThreshold = 2048;
    filename = "astro_ra";

//    memThreshold = 1024+512;
//    filename = "wiki_ts";
//    test_query_baseline(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);


//    test_out_of_order_insertion(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);
    test_interleave_baselines(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);
//
//    memThreshold = 1024;
//    filename = "wiki_ts";
//    test_interleave_baselines(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);

//    memThreshold = 768;
//    test_interleave_baselines(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);
//
//    memThreshold = 1536;
//    test_interleave_baselines(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);

//    memThreshold = 1280;
//    test_interleave_baselines(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);

//    memThreshold = 2048+512;
//    test_interleave_baselines(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);

//    memThreshold = 1792;
//    test_interleave_baselines(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);


//    test_out_of_order_insertion(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);


//    datasize = 1024*4;
//    memThreshold = 1024*2;
//    filename = "wiki_ts";
//    test_interleave_baselines(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);
//    filename = "books";
//    test_interleave_baselines(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);
//    filename = "wiki_ts";
//    test_interleave_baselines(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);
//    test_interleave_append_query(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);
    // append, out-of-order insert, query

//    test_appending(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);

//    filename = "astro_ra";

//    filename = "random0.5";
//    test_appending(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);
//    test_appending(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);
    return 0;
}



/*
 *
 */