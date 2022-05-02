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
#include <regex>

#include "data.h"
#include "film.h"


#define MAX_INT    (((unsigned int)(-1))>>1)


typedef long int key_type;  // if (filename== "books" || filename == "wiki_ts")
//typedef double key_type;   //if (filename== "astro_ra" || filename == "random0.5")



int test_appending(string filename, double memThreshold, long int datasize, long int numcolumn, int pagesize,unsigned int errbnd){
    std::cout << "my Lord, i need You, i trust in You! test the film" << std::endl;
    struct timeval t1, t2;
    double timeuse;

    unsigned long int numkey = ceil(double(datasize*1024*1024)/numcolumn/8);  //int datanum,
    vector<key_type> load_data = loaddata<key_type>(filename,numkey);
    double reserveMem = 5;
    memThreshold -= reserveMem;
    int queryn = 100000;
    double zipf = 1.25;   //zipfian factor
    vector<key_type> load_pointquery ;

    gettimeofday(&t1, NULL);
    unsigned int errs[] = {16};  //4,16,64,256,1024,4096,16384,65536,262144

    gettimeofday(&t1, NULL);
    for (int i = 0; i < (end(errs)-begin(errs));i++){
        gettimeofday(&t1, NULL);
        cout<<"the data set is "<<filename<<",  the data size is "<<datasize<<"M"<< "  numkey is "<<numkey<<endl;
        cout<<"the number of keys is "<<numkey<<", the record size is "<<numcolumn<<endl;
        cout<<"the page size is "<<pagesize*8<<",  the available memory is "<< memThreshold<< endl;

        errbnd = errs[i];   //zipfian factor
        load_pointquery = loadpquery<key_type>(filename,queryn,numkey,zipf,"point");
        filminsert::test_filmappending(errbnd,numkey,pagesize,filename,memThreshold,reserveMem,numcolumn,load_data,load_pointquery,queryn,numkey,datasize,zipf);


        gettimeofday(&t2, NULL);
        timeuse = (t2.tv_sec - t1.tv_sec) + (double) (t2.tv_usec - t1.tv_usec) / 1000000.0;

        cout << "able o_direct disk access time = " << timeuse << endl;  //输出时间（单位：ｓ）
    }
    return 0;

}


vector<string> split(string text) {
    regex ws_re("\\s+");
    vector<string> vector(sregex_token_iterator(text.begin(), text.end(), ws_re, -1),sregex_token_iterator());
    return vector;
}


vector<key_type> random_fun(vector<key_type> &row_vector,unsigned int data_num,double update_ratio){

//    srand(5);
    rand();
    unsigned int update_num = data_num * update_ratio;
    vector<key_type> update_vec;
    update_vec.reserve(update_num);
    std::ostringstream ossb,ossc;
    ossb << update_ratio;
    ossc << update_num;
    string filepath = "/home/wamdm/chaohong/clionDir/Opt_FILMupdate/update_data/update_books" + ossb.str() + "_" + ossc.str() + ".txt";
    string rowfilepath = "/home/wamdm/chaohong/clionDir/Opt_FILMupdate/update_data/row_books" + ossb.str() + "_" + ossc.str() + ".txt";
/*
    ofstream f(filepath, ios::app);
    ofstream rowf(rowfilepath,ios::app);
    for(unsigned int i = 0; i < update_num; i ++){
        unsigned int sample_i = rand() % ((data_num-update_num-100)) +1;
//        unsigned int sample_i= (rand()%(data_num-update_num-10));
        update_vec.emplace_back(row_vector[sample_i]);
        std::vector<key_type>::iterator it = row_vector.begin()+sample_i;
        f<<row_vector[sample_i]<<" ";
        row_vector.erase(it);

    }
    for (unsigned long int i = 0; i < row_vector.size();i ++){
        rowf<<row_vector[i]<<" ";
    }
*/


    ifstream fin(filepath);
    ifstream rowfin(rowfilepath);

    std::string s;
    std::string rows;
    getline(fin, s);
    getline(rowfin,rows);
    key_type tmp;
    key_type rowtmp;
    vector<string> v = split(s);
    vector<string> rowv = split(rows);
    for (auto s1 :v) {
        //getline(fin, s);
        //istringstream istr1(s1);

        std::istringstream iss (s1);
        iss >> tmp;
        update_vec.emplace_back(tmp);
    }
    row_vector.clear();
    for (auto s1 :rowv) {
        //getline(fin, s);
        //istringstream istr1(s1);
        std::istringstream iss (s1);
        iss >> rowtmp;
        row_vector.emplace_back(rowtmp);
    }

    //size_t insert_n = row_vector.size();
    return update_vec;
}


int test_updating(string filename, double memThreshold, long int datasize, long int numcolumn, int pagesize,unsigned int errbnd,double update_ratio,double zipf){
    //std::cout << "my Lord, i need You, i trust in You!  test_updating" << std::endl;
    struct timeval t1, t2;
    double timeuse;

    unsigned long int numkey = ceil(double(datasize*1024*1024)/numcolumn/8);  //int datanum
    vector<key_type> load_data = loaddata<key_type>(filename,numkey);
    double reserveMem = 5;
    memThreshold -= reserveMem;
    int queryn = 100000;

    vector<key_type> load_pointquery ;
    vector<key_type> update_vec;

    gettimeofday(&t1, NULL);

    cout<<"the data set is "<<filename<<",  the data size is "<<datasize<<"M"<< "  numkey is "<<numkey<<endl;
    cout<<"the number of keys is "<<load_data.size()<<", the record size is "<<numcolumn<<endl;
    cout<<"the page size is "<<pagesize*8<<",  the available memory is "<< memThreshold << endl;
    load_pointquery = loadpquery<key_type>(filename,queryn,numkey,zipf,"point");
    update_vec = random_fun(load_data,numkey,update_ratio);

    filminsert::test_filmupdating(errbnd,numkey,pagesize,filename,memThreshold,reserveMem,numcolumn,load_data,update_vec,load_pointquery,queryn,numkey,datasize,zipf,update_ratio);

    gettimeofday(&t2, NULL);
    timeuse = (t2.tv_sec - t1.tv_sec) + (double) (t2.tv_usec - t1.tv_usec) / 1000000.0;

    cout << "able o_direct disk access time = " << timeuse << endl;
    return 0;
}



int main() {

    //std::cout << "my Lord, i need You, i trust in You! start the experiment " << std::endl;
    unsigned int errbnd = 8;
    double memThreshold = 512;
    int pagesize = 1024*64/8;
    string filename = "wiki_ts";  //queryname
    long int numcolumn = 16;
    long int datasize = 4096;
    double zipfs[] = {0.25};
//    test_lru_overhead(datasize,memThreshold);

    datasize = 4096;
    memThreshold = 2048;
    filename = "books";
    //    filename = "wiki_ts";
    //    filename = "random0.5";
    //    filename = "astro_ra";
    filename = "wiki_ts";
    test_appending(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);
    errbnd = 16 ;
    /*
    double update_ratios[] = {0.0001,0.0001};  //0.00001,0.0001,0.001,0.00001,0.0001,0.001  0.0001,0.0001,0.0001,0.001,0.001,0.001,0.01,0.01,0.01

    for (int i = 0; i < (end(update_ratios)-begin(update_ratios));i++) {
        double update_ratio = update_ratios[i];   //zipfian factor
        for (int j = 0; j <(end(zipfs)-begin(zipfs));j++ ){
            double zipf = zipfs[j];
            test_updating(filename, memThreshold, datasize, numcolumn, pagesize, errbnd,update_ratio,zipf);
        }
    }
    */

    test_appending(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);

    return 0;
}

