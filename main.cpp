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

typedef long int key_type;  // if (filename== "books" || filename == "wiki_ts")
//typedef double key_type;   //if (filename== "astro_ra" || filename == "random0.5")


int test_appending(string filename, double memThreshold, long int datasize, long int numcolumn, int pagesize,unsigned int errbnd){
    std::cout << "my Lord, i need You, i trust in You! test the film" << std::endl;
    struct timeval t1, t2;
    double timeuse;

//    typedef long int key_type;  // if (filename== "books" || filename == "wiki_ts")
////    typedef double datatype;   if (filename== "astro_ra" || filename == "random0.5")

    unsigned long int numkey = ceil(double(datasize*1024*1024)/numcolumn/8);  //int datanum,
    vector<key_type> load_data = loaddata<key_type>(filename,numkey);

    int queryn = 100000;
    double zipf = 0.0;   //zipfian factor
//    vector<key_type> load_pointquery1 ;
    vector<key_type> load_pointquery ;

    std::cout << "my Lord, i need You, i trust in You!!!!!" << std::endl;
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
        filminsert::test_filmappending(errbnd,numkey,pagesize,filename,memThreshold,numcolumn,load_data,load_pointquery,queryn,numkey,datasize,zipf);


        gettimeofday(&t2, NULL);
        timeuse = (t2.tv_sec - t1.tv_sec) + (double) (t2.tv_usec - t1.tv_usec) / 1000000.0;

        cout << "able o_direct disk access time = " << timeuse << endl;  //输出时间（单位：ｓ）
    }
    return 0;

}




int main() {

    std::cout << "my Lord, i need You, i trust in You! start the experiment " << std::endl;

//    typedef long int key_type;  // if (filename== "books" || filename == "wiki_ts")
////    typedef double datatype;   if (filename== "astro_ra" || filename == "random0.5")
    unsigned int errbnd = 8;
    double memThreshold = 512;
    int pagesize = 1024*64/8;
    string filename = "wiki_ts";  //queryname
    long int numcolumn = 16;
    long int datasize = 4096;
    double zipf = 0.5;
//    test_lru_overhead(datasize,memThreshold);

    datasize = 4096;
    memThreshold = 2048;

    test_appending(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);
    return 0;
}

