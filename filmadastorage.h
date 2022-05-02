

//
// Created by chaochao on 2021/12/23.
//

#ifndef EXPERIMENTCC12_FILMADASTORAGE_H
#define EXPERIMENTCC12_FILMADASTORAGE_H

#include <fcntl.h>
#include <unistd.h>
#include <memory.h>
//#define _GNU_SOURCE
#include <vector>
#include "filmadalru.h"
#include "film.h"

using namespace std;

namespace filmstorage {
    typedef std::map<std::string,double> infomap;
    template <class key_type,class data_type,class index_type, class lru_type,class disk_type>
    class filmmemory{
    public:
        int totalnum;
        unsigned int inpageid = 0;
        double threshold;
        double reserveMem;
        index_type *index;  //
        lru_type *lru;
        vector<pair<unsigned int,unsigned int>> evicttable;  // pageid , offset in page
        vector<key_type> evictkey;   // the evicted key, used in debug
        vector<pair<unsigned int,unsigned int>*> evictPoss;


        filmmemory( unsigned int numkey,double Threshold, index_type *filmtree , lru_type *LRU){
            totalnum = numkey;
            threshold = Threshold;
            evicttable.reserve(numkey);
            index = filmtree;
            lru = LRU;

        }


        void insert(const std::vector<key_type> keys,const std::vector<key_type> payload,const int error,const int error_recursize ){
            index->build(keys,payload,error,error_recursize,lru);    // batch update
        }

        void append(const std::vector<key_type> keys,const std::vector<key_type> payload,const int error,const int error_recursize ){
            index->update_append(keys,payload,error,error_recursize,lru);    // append-only insertion
            totalnum += 1;
        }

        void cache_append(const std::vector<key_type> keys,const std::vector<key_type> payload,const int error,const int error_recursize ){
            index->batch_append(keys,payload,error,error_recursize,lru);    // append-only insertion optimized by cache
            totalnum += 1;
        }

        void update(const std::vector<key_type> keys,const std::vector<key_type> payload){
            index->update_random(keys,payload,lru);    // out-of-order insertion
        }

        struct memoryusage
        {
            /// the total used memory
            double  totalusemem;
            /// the usage of each component
            infomap meminfo;
            /// Constructor
            inline memoryusage(double total,infomap eachmem)
                    : totalusemem(total), meminfo(eachmem)
            { }
            inline memoryusage()
                    : totalusemem(), meminfo()
            { }

        };

        struct inmempage  // a buffer page in memory
        {
            unsigned int pageid;
            vector<key_type> inmemdata;
            int pagesize;
            int freespace;
            int recordnum;
            inmempage(unsigned int idinpage,int sizepage,int numrecord){
                pageid = idinpage;
                pagesize = sizepage;
                freespace = pagesize;
                recordnum = numrecord;
            }
            inmempage(){
            }


            bool isfull(){  // judge whether the page is full
                if (freespace == 0)
                    return true;
                else
                    return false;
            }

        };

        // define the current page in memory, when the page is full, write to disk
        inmempage *inpage = NULL;
        vector<inmempage*> inmempages ;

        void createinmempage(int sizepage,int numrecord){

            inmempage *curinpage = new inmempage(inpageid,sizepage,numrecord);  //create a page in memory;
            inpageid += 1;
            inpage = curinpage;
        }

        memoryusage simu_computeMemUsage(int transnum){

            infomap indexdatausage = index->show_verify();
            double dataV = ( index->valuesize+1)*8;
            double datausage = double((index->inkeynum-transnum)*(dataV))/1024/1024;
            double addusage = double((index->exkeynum+index->inkeynum)*(1+8)+(index->exkeynum+transnum)*(8+4)+(index->exkeynum+transnum)*sizeof(key_type) )/1024/1024;
            double indexusage = indexdatausage["indexusage"];
            double hashlrusize = sizeof(key_type) + 8+ 16 + 16;
            double no_hashsize = 16+4;  //prev,next, slot-short int
            double lruusage = double(no_hashsize*(index->inkeynum-transnum) + hashlrusize * index->vstats.leaves)/1024/1024;
            double totalmem = lruusage + datausage + indexusage + addusage;
            double leaves = indexdatausage["leaves"];
            double levels = indexdatausage["levels"];
            double innernodes = indexdatausage["inners"];
            double exkeynum = index->exkeynum;

            infomap devied_mem;
            devied_mem.insert(pair<std::string,double>("datausage",datausage));
            devied_mem.insert(pair<std::string,double>("indexusage",indexusage));
            devied_mem.insert(pair<std::string,double>("lruusage",lruusage));
            devied_mem.insert(pair<std::string,double>("addusage",addusage));  //key,flag,pageID,offset
            devied_mem.insert(pair<std::string,double>("leaves",leaves));
            devied_mem.insert(pair<std::string,double>("levels",levels));
            devied_mem.insert(pair<std::string,double>("innernodes",innernodes));
            devied_mem.insert(pair<std::string,double>("totalusage",totalmem));  //key,flag,pageID,offset
            devied_mem.insert(pair<std::string,double>("exkeynum",exkeynum));
            memoryusage res_memusage(totalmem,devied_mem);
            return res_memusage;

        }

        memoryusage runtimecomputeMemUsage(){

            infomap indexdatausage = index->runtimeshow_verify();
            double lruusage = indexdatausage["lruusage"] ;
            double indexusage = indexdatausage["indexusage"];
            double datausage = indexdatausage["datausage"];
            double addusage = indexdatausage["addusage"];
            double totalmem = lruusage + datausage + indexusage + addusage;
            double leaves = indexdatausage["leaves"];
            double levels = indexdatausage["levels"];
            double innernodes = indexdatausage["inners"];
            double exkeynum = index->exkeynum;

            infomap devied_mem;
            devied_mem.insert(pair<std::string,double>("datausage",datausage));
            devied_mem.insert(pair<std::string,double>("indexusage",indexusage));
            devied_mem.insert(pair<std::string,double>("lruusage",lruusage));
            devied_mem.insert(pair<std::string,double>("addusage",addusage));  //key,flag,pageID,offset
            devied_mem.insert(pair<std::string,double>("leaves",leaves));
            devied_mem.insert(pair<std::string,double>("levels",levels));
            devied_mem.insert(pair<std::string,double>("innernodes",innernodes));
            devied_mem.insert(pair<std::string,double>("totalusage",totalmem));  //key,flag,pageID,offset
            devied_mem.insert(pair<std::string,double>("exkeynum",exkeynum));
            memoryusage res_memusage(totalmem,devied_mem);
            return res_memusage;

        }


        memoryusage computeMemUsage(){

            infomap indexdatausage = index->show_verify();
            double lruusage = indexdatausage["lruusage"] ;
            double indexusage = indexdatausage["indexusage"];
            double datausage = indexdatausage["datausage"];
            double addusage = indexdatausage["addusage"];
            double totalmem = lruusage + datausage + indexusage + addusage;
            double leaves = indexdatausage["leaves"];
            double levels = indexdatausage["levels"];
            double innernodes = indexdatausage["inners"];
            double exkeynum = index->exkeynum;

            infomap devied_mem;
            devied_mem.insert(pair<std::string,double>("datausage",datausage));
            devied_mem.insert(pair<std::string,double>("indexusage",indexusage));
            devied_mem.insert(pair<std::string,double>("lruusage",lruusage));
            devied_mem.insert(pair<std::string,double>("addusage",addusage));  //key,flag,pageID,offset
            devied_mem.insert(pair<std::string,double>("leaves",leaves));
            devied_mem.insert(pair<std::string,double>("levels",levels));
            devied_mem.insert(pair<std::string,double>("innernodes",innernodes));
            devied_mem.insert(pair<std::string,double>("totalusage",totalmem));
            devied_mem.insert(pair<std::string,double>("exkeynum",exkeynum));
            memoryusage res_memusage(totalmem,devied_mem);
            return res_memusage;

        }

        pair<unsigned int, unsigned int>* writeevicttable(pair<int, short int> pospage,pair<unsigned int, unsigned int>* evictpos){
            pair<unsigned int,unsigned int>* oldpospage = evictpos;
            evictpos->first = pospage.first;
            evictpos->second = pospage.second;
            return oldpospage;

        }

        pair<unsigned int,unsigned int>*  writeevicttable(pair<unsigned int, unsigned int> pospage){
            int evictpos = evicttable.size();
//            evictkey.push_back(key);
            evicttable.push_back(pospage);
            auto eptr = &evicttable[evictpos];
            return eptr;
        }

        unsigned int evictpagestodisk(disk_type *diskpage){
            unsigned int num = inmempages.size();
            if (num > 0) {

                key_type *buf = new key_type[diskpage->pagesize*num];
                int buf_size = diskpage->pagesize * sizeof(key_type) * num;
                long int fixed_buf_size = diskpage->pagesize * sizeof(key_type);  // 磁盘页固定的大小
                int ret = posix_memalign((void **) &buf, 512, buf_size);

                int fd = open(diskpage->pagefile, O_RDWR | O_DIRECT, 0755);
                auto seekoff = diskpage->nextpageid * fixed_buf_size;
                lseek(fd, diskpage->nextpageid * fixed_buf_size, SEEK_SET);

                if (seekoff <0){
                    cout<< "seekoff <0" << endl;
                }
                int offset = 0;
                for (auto &v: inmempages) {
//                    size_t len = sizeof(*buf);
                    memcpy(buf + offset * diskpage->pagesize, &v->inmemdata[0], v->inmemdata.size() * sizeof(key_type));
                    ++offset;

                }
                ret = write(fd, buf, buf_size);
                if (ret <= 0){
                    cout << "ret <= 0" << endl;
                }
                diskpage->nextpageid += num;
                free(buf);
                close(fd);
                for (int i = 0;i < inmempages.size();i++){
                    delete inmempages[i];
                    inmempages[i] = NULL;
                }
                inmempages.resize(0);
            }
            return num;
        }

        unsigned short int runtimeevictpagestodisk(disk_type *diskpage){
            unsigned short int pnum = inmempages.size();
            unsigned short int num = inmempages.size();
            int fd = open(diskpage->pagefile, O_RDWR | O_DIRECT, 0755);
            unsigned int fixed_buf_size = diskpage->pagesize * sizeof(key_type);
            lseek(fd, diskpage->nextpageid * fixed_buf_size, SEEK_SET);
            int i = 0;
            while (pnum > diskpage->blocknum) {
                pnum -= diskpage->blocknum;
                key_type *buf = new key_type[diskpage->pagesize* diskpage->blocknum];

                int buf_size = diskpage->pagesize * sizeof(key_type) * diskpage->blocknum;
                int ret = posix_memalign((void **) &buf, 512, buf_size);

                int offset = 0;
                for(int k = 0;k<diskpage->blocknum;k++){
                    memcpy(buf + offset * diskpage->pagesize, &inmempages[i*diskpage->blocknum+k]->inmemdata[0], inmempages[i*diskpage->blocknum+k]->inmemdata.size() * sizeof(key_type));
                    ++offset;
                }
                i += 1;
                ret = write(fd, buf, buf_size);
                free(buf);
            }

            key_type *buf = new key_type[diskpage->pagesize*pnum];
            int buf_size = diskpage->pagesize * sizeof(key_type) * pnum;
            int ret = posix_memalign((void **) &buf, 512, buf_size);
            int offset = 0;
            for(int k = 0;k<pnum;k++){
                memcpy(buf + offset * diskpage->pagesize, &inmempages[i*diskpage->blocknum+k]->inmemdata[0], inmempages[i*diskpage->blocknum+k]->inmemdata.size() * sizeof(key_type));
                ++offset;
            }
            ret = write(fd, buf, buf_size);
            free(buf);
            diskpage->nextpageid += num;
            close(fd);

            for (int mi = 0;mi < inmempages.size();mi++){
                delete inmempages[mi];
                inmempages[mi] = NULL;
            }
            inmempages.resize(0);
            return num;
        }

        void evictkeytoinpage(key_type ekey,data_type edata, disk_type *diskpage,pair<unsigned int,unsigned int>* evictpos){

            if (!(inpage->recordnum --))
            {
                inmempages.emplace_back(inpage);
//                diskpage->odirectenterpage(inpage->inmemdata);// write to disk
                createinmempage(diskpage->pagesize,diskpage->recordnum);
                inpage->recordnum -=1;

            }
            unsigned int offset = inpage->inmemdata.size();
            inpage->inmemdata.push_back(ekey);
            inpage->freespace -= 1;
            inpage->inmemdata.insert( inpage->inmemdata.begin()+(inpage->pagesize - inpage->freespace),edata.begin(),edata.end());
            inpage->freespace -= (index->valuesize);

            writeevicttable(pair<unsigned int,unsigned int>(inpage->pageid, offset),evictpos);

        }

        int rangeevictkeytoinpage(key_type ekey,data_type edata, disk_type *diskpage,pair<unsigned int,unsigned int>* evictpos){
            int num =0;
            if (!(inpage->recordnum --))
            {

                diskpage->odirectenterpage(inpage->inmemdata);// write to disk
                createinmempage(diskpage->pagesize,diskpage->recordnum);
                inpage->recordnum -=1;
                num += 1;
            }
            short int offset = inpage->inmemdata.size();
            inpage->inmemdata.push_back(ekey);
            inpage->freespace -= 1;
            inpage->inmemdata.insert( inpage->inmemdata.begin()+(inpage->pagesize - inpage->freespace),edata.begin(),edata.end());
            inpage->freespace -= (index->valuesize);
//            short int offset2 = (inpage->pagesize - index->valuesize-1 - inpage->freespace);
            writeevicttable(pair<unsigned int,unsigned short int>(inpage->pageid, offset),evictpos);
            return num;
        }

        pair<unsigned int,unsigned short int>* evictkeytoinpage(key_type ekey,data_type edata, disk_type *diskpage){

            if (!(inpage->recordnum --))
            {
                diskpage->odirectenterpage(inpage->inmemdata);// write to disk
                createinmempage(diskpage->pagesize,index->recordsize);
                inpage->recordnum -=1;
            }
            inpage->inmemdata.push_back(ekey);
            inpage->freespace -= 1;
            inpage->inmemdata.insert( inpage->inmemdata.begin()+(inpage->pagesize - inpage->freespace),edata.begin(),edata.end());
            inpage->freespace -= (index->recordsize-1);

            auto posevict = writeevicttable(pair<unsigned int,unsigned short int>(inpage->pageid, (inpage->pagesize - index->recordsize - inpage->freespace)));
            return posevict;
        }

        void filmtransfer(int transleaves,disk_type *diskpage){   //perform the transfer procedure,
            //

            if (inpage == NULL){
                createinmempage(diskpage->pagesize,diskpage->recordnum);    //create a page in memory;
            }
            if (index->m_transleaf == NULL){
                index->m_transleaf = index->leaflevel.leafpieces[0];
            }
            if  (transleaves == 1 && this->index->Error > 256)
            {
                auto midtransflag = simu_computeMemUsage(index->m_transleaf->slotkey.size());
                if (midtransflag.totalusemem > threshold ){
                    for (int k = 0; k < index->m_transleaf->slotkey.size();k++)
                    {
                        if (!(inpage->recordnum --))
                        {
                            diskpage->odirectenterpage(inpage->inmemdata);
                            delete inpage;
                            inpage = NULL;
                            createinmempage(diskpage->pagesize,diskpage->recordnum);
                            inpage->recordnum -=1;
                        }

                        inpage->inmemdata.push_back(index->m_transleaf->slotkey[k]);
                        inpage->freespace -= 1;
                        auto transdata = (adalru::Node<unsigned int, std::vector<key_type> >*) index->m_transleaf->slotdata[k].second;
                        inpage->inmemdata.insert( inpage->inmemdata.begin()+(inpage->pagesize - inpage->freespace),transdata->value.begin(),transdata->value.end());
                        inpage->freespace -= (index->valuesize);
                        index->m_transleaf->intrachain.remove_node(transdata);//pop from intrachain
                        delete transdata;
                        transdata = NULL;
                        index->m_transleaf->slotdata[k].first = false;
                        auto eptr = writeevicttable(pair<int,unsigned short int>(inpage->pageid, (inpage->pagesize - index->valuesize-1 - inpage->freespace)));
                        index->m_transleaf->slotdata[k].second  = eptr ;  //
                    }
                    lru->remove(index->m_transleaf->startkey);// pop from interchain
                    index->inkeynum -= index->m_transleaf->slotkey.size();
                    index->exkeynum += index->m_transleaf->slotkey.size();
                    index->m_transleaf = index->leaflevel.leafpieces[index->leafsplit+1];//evict data from the first leaf
                    if (index->m_transleaf->startkey == 0){
                        cout<< "index->m_transleaf->startkey == 0 "<< endl;
                    }
                    index->leafsplit += transleaves;
                }
                else{
                    auto midtrans = computeMemUsage();
                    int split = 0;
                    while (midtrans.totalusemem > threshold){
                        double ratio = (midtrans.totalusemem-threshold)/(midtrans.totalusemem-midtransflag.totalusemem ) ;
                        int trannum = ceil(index->m_transleaf->slotkey.size() * ratio+100);
                        for (int k = 0; k < trannum;k++)
                        {
                            if (!(inpage->recordnum --))
                            {
                                diskpage->odirectenterpage(inpage->inmemdata);

                                delete inpage;
                                inpage = NULL;
                                createinmempage(diskpage->pagesize,diskpage->recordnum);
                                inpage->recordnum -=1;
                            }

                            inpage->inmemdata.push_back(index->m_transleaf->slotkey[split+k]);
                            inpage->freespace -= 1;
                            auto transdata = (adalru::Node<unsigned int, std::vector<key_type> >*) index->m_transleaf->slotdata[split+k].second;
                            inpage->inmemdata.insert( inpage->inmemdata.begin()+(inpage->pagesize - inpage->freespace),transdata->value.begin(),transdata->value.end());
                            inpage->freespace -= (index->valuesize);
                            index->m_transleaf->intrachain.remove_node(transdata);//pop from intrachain
                            delete transdata;
                            transdata = NULL;
                            index->m_transleaf->slotdata[split+k].first = false;
                            auto eptr = writeevicttable(pair<int,unsigned short int>(inpage->pageid, (inpage->pagesize - index->valuesize-1 - inpage->freespace)));
                            index->m_transleaf->slotdata[split+k].second  = eptr ;  //
                        }
                        split += trannum;
                        index->inkeynum -= trannum;
                        index->exkeynum += trannum;
                        midtrans = computeMemUsage();
                    }

                }

            }
            else{
                for (int i = 0; i< transleaves; i++){
                    for (int k = 0; k < index->m_transleaf->slotkey.size();k++)
                    {
                        if (!(inpage->recordnum --))
                        {
                            diskpage->odirectenterpage(inpage->inmemdata);
                            delete inpage;
                            inpage = NULL;
                            createinmempage(diskpage->pagesize,diskpage->recordnum);
                            inpage->recordnum -=1;
                        }

                        inpage->inmemdata.push_back(index->m_transleaf->slotkey[k]);
                        inpage->freespace -= 1;
                        auto transdata = (adalru::Node<unsigned int, std::vector<key_type> >*) index->m_transleaf->slotdata[k].second;
                        inpage->inmemdata.insert( inpage->inmemdata.begin()+(inpage->pagesize - inpage->freespace),transdata->value.begin(),transdata->value.end());
                        inpage->freespace -= (index->valuesize);
                        index->m_transleaf->intrachain.remove_node(transdata);//pop from intrachain
                        delete transdata;
                        transdata = NULL;
                        index->m_transleaf->slotdata[k].first = false;
                        auto eptr = writeevicttable(pair<int,unsigned short int>(inpage->pageid, (inpage->pagesize - index->valuesize-1 - inpage->freespace)));
                        index->m_transleaf->slotdata[k].second  = eptr ;  //
                    }
                    lru->remove(index->m_transleaf->startkey);// pop from interchain
                    index->inkeynum -= index->m_transleaf->slotkey.size();
                    index->exkeynum += index->m_transleaf->slotkey.size();
                    index->m_transleaf = index->leaflevel.leafpieces[index->leafsplit+i+1];//evict data from the first leaf
                    if (index->m_transleaf->startkey == 0){
                        cout<< " index->m_transleaf->startkey == 0"<< endl;
                    }
                }
                index->leafsplit += transleaves;
            }


        }


        template<typename stat_type>
        bool runtimejudgetrans(stat_type range_stats){
            struct timeval ct1, ct2;
            double ctimeuse;
            gettimeofday(&ct1, NULL);
            memoryusage res_memusage = runtimecomputeMemUsage();
            gettimeofday(&ct2, NULL);
            ctimeuse = (ct2.tv_sec - ct1.tv_sec) + (double) (ct2.tv_usec - ct1.tv_usec) / 1000000.0;
            range_stats.computetimeuse += ctimeuse;
            if (res_memusage.totalusemem > threshold) //perform transfer process
                return true;
            else
                return  false;
        }

        pair<bool,memoryusage> judgetransfer(){
            memoryusage res_memusage;
            res_memusage = computeMemUsage();

            map<string,double>::iterator iter;

            for(iter = res_memusage.meminfo.begin(); iter != res_memusage.meminfo.end(); iter++)
                cout<<iter->first<<" "<<iter->second<<" *** ";
            cout<<endl;
            if (res_memusage.totalusemem > threshold){//perform transfer process
                return pair<bool,memoryusage>(true,res_memusage);
            }
            return  pair<bool,memoryusage>(false,res_memusage);
        }
    };

    typedef std::map<std::string,double> infomap;
    template <class key_type>
    class filmdisk{
    public:
        const char *pagefile;
        int pagesize;
        int nextpageid;
        int recordnum;
        int recordsize;
        int blocknum;
        filmdisk( const char* diskfile,int sizepage,int numrecord,int sizerecord) {
            pagefile = diskfile;
            pagesize = sizepage;
            recordnum = numrecord;
            blocknum = 1024*1024/(sizepage*8);
            nextpageid = 0;
            recordsize = sizerecord;
        }


        vector<key_type> readfromdisk(pair<unsigned int,unsigned short int> diskpos ,int sizerecord){   // use the buffer of memory
            FILE *fdisk;// 读取磁盘文件
            fdisk = fopen(pagefile,"rb+");
            fseek(fdisk,diskpos.first*pagesize*8,SEEK_SET);
            vector<key_type> pagedata;
            key_type tmp;

            for (int i = 0; i < pagesize; i ++){
                fread(&tmp, sizeof(long int), 1, fdisk);
                pagedata.push_back(tmp);
//                cout<< tmp << " ";
            }

            vector<key_type> res;

            for (int i = 0; i < sizerecord; i ++){
                res.push_back(pagedata[diskpos.second+i]);
            }
//            cout<<"Jesus, You are my refuge! "<<endl;
            fclose(fdisk);
            return res;

        }

        // use the system o_direct mode to read file
        vector<key_type> odirectreadfromdisk(pair<unsigned int,unsigned short int> diskpos ,int sizerecord){
            int fd;
            key_type *buf;
            vector<key_type> res;
            unsigned long int buf_size = pagesize*8;
            int ret = posix_memalign((void **)&buf, 512, buf_size);
            memset(buf, 'c', buf_size);

            fd = open(pagefile, O_RDWR | O_DIRECT , 0755);
            /*
            if (fd < 0){
                perror("open ./direct_io.data failed");
                exit(1);
            }
             */

            ret = pread(fd, buf, buf_size,diskpos.first*buf_size);
            if (ret <= 0){
                cout << "Jesus, i need You!" << endl;
            }
            for (int i = 0; i < sizerecord; i ++){
                res.push_back(buf[diskpos.second+i]);
            }
//            cout<<"Jesus, You are my refuge! "<<endl;
            free(buf);
            close(fd);
            return res;
        }

        // use the system o_direct mode to read file
        vector<key_type> odirectreadfromdisk(pair<unsigned int,unsigned int>* diskpos){
            int fd;
            key_type *buf;
            vector<key_type> res;
            unsigned long int buf_size = pagesize*8;
            int ret = posix_memalign((void **)&buf, 512, buf_size);
            memset(buf, 'c', buf_size);

            fd = open(pagefile, O_RDWR | O_DIRECT , 0755);

            ret = pread(fd, buf, buf_size,diskpos->first*buf_size);
            if (ret <= 0){
                cout << "Jesus, i need You!" << endl;
            }
            for (int i = 0; i < recordsize; i ++){
                res.push_back(buf[diskpos->second+i]);
            }
//            cout<<"Jesus, You are my refuge! "<<endl;
            free(buf);
            close(fd);
            return res;
        }

        int enterpage(vector<key_type> datas)
        {

            FILE *fdisk;
            fdisk = fopen(pagefile,"rb+");
            fseek(fdisk,nextpageid * (pagesize)*8,SEEK_SET);
            //fseek(fdisk, 0, SEEK_END);
            nextpageid += 1;
//            cout<<datas.size()<<endl;
            for(int i = 0 ; i < datas.size() ; i++) {
                fwrite(&datas[i], sizeof(key_type), 1, fdisk);
//                cout << datas[i] <<" ";
            }

            fflush(fdisk);
            fclose(fdisk);

        }

        void odirectenterpage(vector<key_type> datas)
        {

            int fd;
            key_type *buf = new key_type[pagesize];

            unsigned long int buf_size = pagesize*sizeof(key_type);
            int ret = posix_memalign((void **)&buf, 512, buf_size);
//            memcpy(buf, &datas[0], datas.size()*sizeof(key_type));
//            int ret = posix_memalign((void **)&buf, 512, buf_size);
//            memset(buf, 'c', 4096);
            memcpy(buf, &datas[0], datas.size()*sizeof(key_type));

            fd = open(pagefile, O_RDWR | O_DIRECT , 0755);
            /*
            if (fd < 0){
                perror("open ./direct_io.data failed");
                exit(1);
            }
             */

            ret = pwrite(fd, buf, buf_size,nextpageid * buf_size);
            if (ret <= 0){
                cout << "ret <= 0" << endl;
            }
            nextpageid += 1;

//            cout<<"Jesus, You are my refuge! odirect write "<<endl;

            free(buf);
            close(fd);

//            cout<< "*********************** Jesus, finished one transfer*************************"<<endl;
        }

    };
}


#endif //EXPERIMENTCC12_FILMADASTORAGE_H
