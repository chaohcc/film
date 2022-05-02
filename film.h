//
// Created by CCMa on 2022/1/28.
//

#ifndef FILMINSERT_FILM_H
#define FILMINSERT_FILM_H

#include <stdlib.h>
#include <queue>
#include <sys/types.h>
// *** Required Headers

#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>
#include <memory>
#include <cstddef>
#include <cassert>
#include <map>
#include <string>
#include <vector>
#include <cstdio>
#include <stdio.h>



#include "pwlf.h"
#include "filmadastorage.h"
#include "filmadalru.h"

using namespace std;

namespace filminsert{
#define PGM_SUB_EPS(x, epsilon) ((x) <= (epsilon) ? 0 : ((x) - (epsilon)))
#define PGM_SUB_EPS2(x, epsilon,size) ((x) <= (epsilon) || (x) - (epsilon)>=(size) ? 0 :((x) - (epsilon)) )
#define PGM_ADD_EPS(x, epsilon, size) ((x) + (epsilon) + 2 >= (size) ? (size) : (x) + (epsilon) + 2)


    // record the access_stats, the number of access memory, the number of access disk and the number of cross memory and disk;
    struct access_stats {
        int memnum = 0;
        int disknum =0;
        int crossnum =0;
        int diskpagenum =0;
        int crosspagenum = 0;
        int writenum = 0;
        double querytimeuse =0.0;
        double computetimeuse =0.0;
        double lrutime = 0.0;
        double disktime = 0.0;
        double rdisktime = 0.0;
        double wdisktime = 0.0;
        double memtime = 0.0;
        int pagecross = 0;  // 磁盘页 seek 的跨度
        double zipffactor = 0.0;
        string qtype = "point";

        void print_stats(){
            struct timeval t1;
            gettimeofday(&t1, NULL);
            cout<<"ID " <<double(t1.tv_sec) << " qtype "<< qtype << " zipffactor "<< zipffactor  << " memnum " << memnum << " disknum " << disknum << " crossnum "<< crossnum << " diskpagenum ";
            cout<< diskpagenum<<" crosspagenum " << crosspagenum <<" pagecross " << pagecross<<" writenum " << writenum<< " querytimeuse " <<querytimeuse <<" lrutimeuse " <<lrutime <<" disktimeuse " << disktime <<" readdisktimeuse " <<rdisktime <<" writedisktimeuse " <<wdisktime <<"\n";
            ofstream savefile;
            savefile.open("/home/wamdm/chaohong/clionDir/Opt_FILMupdate/result/adalru_film_performance.txt",ios::app);

            savefile<<"ID " <<double(t1.tv_sec) << " qtype "<< qtype  << " zipffactor "<< zipffactor  << " memnum " << memnum << " disknum " << disknum << " crossnum "<< crossnum << " diskpagenum ";
            savefile << diskpagenum<<" crosspagenum " << crosspagenum <<" pagecross " << pagecross<<" writenum " << writenum<<" querytimeuse " <<querytimeuse  <<" lrutimeuse " <<lrutime <<" disktimeuse " <<disktime << " readdisktimeuse " <<rdisktime <<" writedisktimeuse " <<wdisktime <<"\n";
            savefile << flush;
            savefile.close();

        }
    };


    /**
 * A struct that stores the result of a query to a @ref PGMIndex, that is, a range [@ref lo, @ref hi)
 * centered around an approximate position @ref pos of the sought key.
 */

    inline int cmp(const void *a, const void *b) {
        uint64_t aa = *(uint64_t *) a;
        uint64_t bb = *(uint64_t *) b;
        if (aa < bb)
            return -1;
        if (aa > bb)
            return 1;
        return 0;
    }
    struct ApproxPos {
        int pos; ///< The predicted position of the key.
        int lo;  ///< The lower bound of the range.
        int hi;  ///< The upper bound of the range.
    };

    struct film_stats {
        /// Number of items in the film
        int leaves;
        /// Number of inners in the B+ tree
        vector<int> inners;
        int innernum;
        /// Number of levels
        int numlevel;

        /// Zero initialized
        inline film_stats()
                : leaves(0), numlevel(0) {
        }

    };

    template <typename key_type,typename value_type>
    class FILMinsert{
    public:

        struct Leafpiece;  /

        struct sortPiece:Leafpiece{

            vector<key_type> slotkey;

            vector< pair<bool, void* > > slotdata;   //adalru::Node<key_type, std::vector<key_type> >*sss
            adalru::localLRU <unsigned int,std::vector<key_type> >  intrachain;    //unsigned short int

            sortPiece() = default;

            friend inline bool operator<(const sortPiece &p, const key_type &key) { return p.key < key; }
            friend inline bool operator<(const key_type &key, const sortPiece &p) { return key < p.key; }
            friend inline bool operator<(const sortPiece &p1, const sortPiece &p2) { return p1.key < p2.key; }


            /**
         * Returns the approximate position of the specified key.
         * @param k the key whose position must be approximated
         * @return the approximate position of the specified key
         */
            inline int operator()(const key_type key)  const{
                int begin = 0, end = slotkey.size()-1,mid;  // find last Equal or small, 找到最后一个大于或等于 key 的
                while(begin <= end){
                    mid = (end + begin) / 2;
                    if(slotkey[mid] >= key) {
                        end = mid-1;
                    } else
                        begin = mid +1;
                }
                return end;
            }



            inline bool insert(key_type key, vector<key_type> payload) {
                if (key == 1296099)
                {
                    cout << "what's wrong, my lovely Lord!"<< endl;
                }
                if (this->slotkey.empty()){
                    this->slotkey.emplace_back(key);
                    this->slotdata.emplace_back( intrachain.put(0,payload));
                    return true;
                }
                else{
                    int begin = 0, end = slotkey.size()-1,mid;
                    while(begin <= end){
                        mid = (end + begin) / 2;
                        if(this->slotkey[mid] >= key) {
                            end = mid-1;
                        } else
                            begin = mid +1;
                    }

                    intrachain.modify(end);
                    this->slotdata.insert( slotdata.begin()+end+1,intrachain.put(end+1,payload));
                    this->slotkey.insert(slotkey.begin()+end+1,key);
                    return true;
                }

            }

            inline pair<bool,unsigned long> is_in_leaf(const key_type &key,const int error){
                auto pos = (*this)(key);
                pair<bool,unsigned long> is_in(false,0);
                int slotlo = PGM_SUB_EPS(pos, error + 1) ;  //,(--itlevel)->size()
                int slothi = PGM_ADD_EPS(pos, error,(this->slotkey.size()-1));
                for (; slotlo<=slothi;slotlo++){
                    if (slotlo < 0)
                    {
                        //cout<< "slotlo < 0" << endl;
                    }
                    if (key == this->slotkey[slotlo])
                    {   is_in.first = true;
                        is_in.second = slotlo;
                        break;}
                }

                return is_in;

            }

        };
        struct Innerpiece;
        struct Innerlevel;
        struct Leaflevel;
        unsigned int totalnum;
        unsigned int inkeynum=0;
        unsigned int exkeynum=0;
        int valuesize = 0;
        int Error ;
        int ErrorRecursive ;
        int leafsplit = 0;
        Leafpiece *m_transleaf = NULL;
        Leafpiece *m_tailleaf = NULL;
        Innerpiece *i_innerpiece = NULL;
        queue< pair<key_type,unsigned int> > buffqueue;
        film_stats vstats;

        Innerpiece* root;
        sortPiece* sort_list = NULL;
//        typedef std::pair< std::vector<Leafpiece*>,vector<internal::insertPWLF<key_type, int>*> > leafleveltype;
        Leaflevel leaflevel;   //  the leafpieces composing the leaflevel
//        typedef std::pair< std::vector<Innerpiece>, vector<internal::insertPWLF<key_type, int>*> > innerleveltype;
        vector<Innerlevel*> innerlevels;


        /**
         * constructs an empty film
         */
//        FILMinsert() = default;

        /**
         * consttructs index on the given keys vector
         * @param keys the vector of keys to be indexed (sorted). @param payload the value of each key
         */

        FILMinsert(unsigned int n ,int err1,int err2)
        {
            totalnum =n,Error = err1, ErrorRecursive =err2;
//            internal::insertPWLF<key_type, int> opt(err1);
//            leaflevel.second.push_back(&opt) ;

        }

        /*
         * the method to construct film, that build FILM
         */

        template<class lru_type>
        inline void update_append(std::vector<key_type> keys,std::vector<key_type> payload,unsigned int error,unsigned int error_recursize,lru_type &interchain ){

            if (innerlevels.size()>1)
                root = &innerlevels.back()->innerpieces[0];
            this->verify_film(&vstats);
            cout<<"film has been built successfully" << endl;

        }

        template<class lru_type>    // append-only fashion (batch-by-batch) the batch is set to use the optimization of cache
        inline void batch_append(std::vector<key_type> keys,std::vector<key_type> payload,unsigned int error,unsigned int error_recursize,lru_type &interchain ){
            // build leaf level
            // size_t error,std::vector<key_type> keys,std::vector<key_type> payload,
            //filmtype *filmada, unsigned int &inkeynum,leaf_type* m_taillea
            std::pair<size_t,std::vector<key_type>>  n_parts = internal::batch_segmentation(error,keys,payload,this,inkeynum,this->m_tailleaf,interchain);   //<key_type,Leafpiece>
//            std::cout<< "You are my refuge, from ever to ever"<<endl;
            // on the basis of leaf level, build the upper level until just has one piece, that the make_segmentation just return 1.
            if (innerlevels.size()>=1)
                root = &innerlevels.back()->innerpieces[0];
            this->verify_film(&vstats);
            cout<<"film has been built successfully" << endl;

        }

        template<class lru_type>
        inline void update_random(std::vector<key_type> keys,std::vector<key_type> payload,lru_type &interchain){

            if (sort_list == NULL)
            {
                sort_list = new sortPiece();
            }
            for (unsigned int i = 0; i < keys.size(); i ++){
                this->update_one(keys[i],payload);
            }
            this->inkeynum += keys.size();
            // add sort_list into interchain

            interchain->put((*sort_list).slotkey[0],sort_list);
        }


        /// film recursive find the key has much information which is needs to be return.
        struct result_find
        {
            /// find result flags
            bool find;

            /// this flag indicate whether the queried data is in memory or on disk
            bool  flags;

            /// The leaf node the findkey belong to
            Leafpiece *findleaf;

            /// the location of key in leaf node
            int slot;

            /// Constructor of a result with a specific flag, this can also be used
            /// as for implicit conversion.
            inline result_find(bool f = true)
                    : find (true),flags(f), findleaf(),slot()
            { }

            /// Constructor with a lastkey value.
            inline result_find(bool find_f,bool in_or_exf ,  Leafpiece * find_leaf,int loc_slot)
                    : find(find_f), flags(in_or_exf),  findleaf(find_leaf),slot(loc_slot)
            { }

            /// Test if this result object has a given flag set.
            bool has(bool f) const
            {
                return (flags & f) != 0;
            }

            bool in_or_out(){
                pair<bool,void*> find_data = findleaf->slotdata[slot];
                if (findleaf->slotdata[slot].first == true)
                    return true;
                else
                    return false;
            }
        };

        /** search method
         * the methods used in query, @param key is the queried key, return the leafpiece responsible for the given key
         */
        result_find search_one(const key_type key)  const{
            //
            if (vstats.innernum > 1){
                auto itlevel = innerlevels.end() -1;

                auto pos = (*root)(key);
                auto lo = PGM_SUB_EPS(pos, ErrorRecursive + 1) + (*(--itlevel))->innerpieces.begin();  //,(--itlevel)->size()
                auto hi = PGM_ADD_EPS(pos, ErrorRecursive,(*itlevel)->innerpieces.size()-1) + (*itlevel)->innerpieces.begin();
                auto curit = (*(itlevel--))->innerpieces.begin();

                static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece)
                if (key >= hi->startkey){curit = hi;}
                else{
                    if constexpr (4 <= linear_search_threshold ){
                        for (; std::next(lo)->startkey <= key; ++lo)
                            continue;
                        curit = lo;
                    }
                    else{ curit = std::lower_bound(lo, hi, key);
                    }
                }

                for (; itlevel >= innerlevels.begin();--itlevel  ){
                    pos = (*curit)(key);
                    lo = PGM_SUB_EPS(pos, ErrorRecursive + 1) + (*itlevel)->innerpieces.begin();
                    hi = PGM_ADD_EPS(pos, ErrorRecursive,(*itlevel)->innerpieces.size()-1) + (*itlevel)->innerpieces.begin();

                    static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece)
                    if (key >= hi->startkey){curit = hi;}
                    else{
                        if (key <= (*lo).startkey){
                            pos = std::min<size_t>((*curit)(key), (*curit).intercept);
                            lo = PGM_SUB_EPS(pos, ErrorRecursive + 1) + (*itlevel)->innerpieces.begin();
                            hi = PGM_ADD_EPS(pos, ErrorRecursive,(*itlevel)->innerpieces.size()-1) + (*itlevel)->innerpieces.begin();
                        }
                        if constexpr (4 <= linear_search_threshold ){   //ErrorRecursive
                            for (; std::next(lo)->startkey <= key; ++lo)
                                continue;
                            curit = lo;
                        }
                        else{ curit = std::lower_bound(lo, hi, key);
                        }
                    }
                }

                // leaf level
                auto curleaf = leaflevel.leafpieces.begin();
                pos = (*curit)(key);
//                auto leaflo = PGM_SUB_EPS(pos, Error + 2) + leaflevel.begin();
                auto leaflo = PGM_SUB_EPS2(pos, Error+1 ,leaflevel.leafpieces.size()) + leaflevel.leafpieces.begin();
                auto leafhi = PGM_ADD_EPS(pos, Error,leaflevel.leafpieces.size()-1) + leaflevel.leafpieces.begin();
//                if (key <= (*leaflo)->startkey) {
//                    pos = std::min<size_t>((*curit)(key), (*curit).intercept);
//                    leaflo = PGM_SUB_EPS(pos, Error + 1) +  leaflevel.leafpieces.begin();
//                    leafhi = PGM_ADD_EPS(pos, Error, leaflevel.leafpieces.size() - 1) +  leaflevel.leafpieces.begin();
//                }

                Leafpiece* x = *leaflo;
                Leafpiece* y = *leafhi;
                if (key >= (*leafhi)->startkey){curleaf = leafhi;}
                else{
                    if constexpr (4 <= linear_search_threshold ){   //ErrorRecursive
                        for (; (*leaflo)->endkey < key; leaflo++)
                            continue;
                        curleaf = leaflo;
                    }
                    else{ curleaf = std::lower_bound(leaflo, leafhi, key);
                    }
                }
                Leafpiece* a = *curleaf;
                pos = (*a)(key);

                int slotlo = PGM_SUB_EPS(pos, ErrorRecursive + 1) ;  //,(--itlevel)->size()
                int slothi = PGM_ADD_EPS(pos, ErrorRecursive,(a->slotkey.size()-1));
                int resslot = slotlo;
                for (; resslot<slothi;resslot++){
                    if (key == a->slotkey[resslot])
                    {
                        break;
                    }
                }
                if (key != a->slotkey[resslot])
                {
                    if (sort_list==NULL) {
                        cout << "sort_list==NULL" << endl;
                    }
                    else {  // find in sort_list
                        auto sortpos = (*sort_list)(key);
                        if (sort_list->slotkey[sortpos+1] != key){
                            auto sss = sort_list->slotkey[sortpos+1];

                        }
                        else{
                            result_find index_res = result_find(false,sort_list->slotdata[sortpos].first,sort_list,sortpos);   // sort_list data
                            return index_res;
                        }
                    }
                }

                result_find index_res = result_find(true,a->slotdata[resslot].first,a,resslot);   // regular data
                return index_res;
            }
            else{
                static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece)

                auto curleaf = leaflevel.leafpieces.begin();
                auto pos = (*root)(key);
//                auto l = PGM_SUB_EPS2(pos, Error + 1,leaflevel.size());
                auto leaflo = PGM_SUB_EPS2(pos, Error + 1,leaflevel.leafpieces.size()) + leaflevel.leafpieces.begin();
                auto leafhi = PGM_ADD_EPS(pos, Error,leaflevel.leafpieces.size()-1) + leaflevel.leafpieces.begin();
                auto h = PGM_ADD_EPS(pos, Error,leaflevel.leafpieces.size()-1);
//
//                Leafpiece*  y= *leafhi;
//                Leafpiece*  z= *leaflo;
                if (key >= (*leafhi)->startkey){curleaf = leafhi;}
                else{
                    if constexpr (4 <= linear_search_threshold ){   //ErrorRecursive
                        for (; (*leaflo)->endkey < key; leaflo++){
//                            Leafpiece*  x= *leaflo;
                            continue;
                        }

                        curleaf = leaflo;
                    }
                    else{ curleaf = std::lower_bound(leaflo, leafhi, key);
                    }
                }
                Leafpiece* a = *curleaf;
                pos = (*a)(key);

                int slotlo = PGM_SUB_EPS(pos, ErrorRecursive + 1) ;  //,(--itlevel)->size()
                int slothi = PGM_ADD_EPS(pos, ErrorRecursive,(a->slotkey.size()-1));
                int resslot = slotlo;
                for (; resslot<slothi;resslot++){
                    if (key == a->slotkey[resslot])
                    {
                        break;
                    }
                }
                result_find index_res = result_find(true,a->slotdata[resslot].first,a,resslot);   // regular data
                return index_res;
            }
        }

        bool update_one(const key_type key, const vector<key_type> payload)  {

            if (key < this->m_tailleaf->endkey){

                this->sort_list->insert(key,payload);
//                cout << "thank You, my Lord! " << endl ;
            }
        }

        /**
          * statistics information
          **/

        inline size_t leafpieces_count() const{
            if (leaflevel.leafpieces.empty())
                return 0;
            else{
                size_t leaves = leaflevel.leafpieces.size();
                return leaves;
            }
        }

        inline vector<int > each_inner_count() const{
            vector<int > inners(innerlevels.size());
            for (int i = innerlevels.size()-1; i >=0; --i)
                inners[i] = innerlevels[i]->innerpieces.size();
            return inners;
        }

        inline int inner_count() const{
            int innersnum = 0;
            for (int i = innerlevels.size()-1; i >=0; --i)
                innersnum += innerlevels[i]->innerpieces.size();
            return innersnum;
        }

        size_t height() const {
            return innerlevels.size()+1;
        }

        inline void verify_film(film_stats *fstats){
            fstats->leaves = leafpieces_count();
            fstats->innernum = inner_count();
            auto a = each_inner_count();
            fstats->inners.assign( a.begin(), a.end());
            fstats->numlevel = height();
        }


        const map<string, double> runtimeshow_verify()
        {
            std::map<std::string,double> treeinfo;
//            this->verify_film(&vstats);
            double innodesize = 24;  //
            double leafnodesize = 32 ; //
            double hashlrusize = 48;
            double no_hashsize = 20;  //prev,next, slot-short int
            double dataV = (valuesize+1)*8;
            double datausage = double(inkeynum*(dataV))/1048576;
            double addusage = double((exkeynum+inkeynum)*(1+8)+exkeynum*(4+4)+exkeynum*sizeof(key_type) )/1048576;  //
            // addusage: (exkeynum+inkeynum)*(1+8) —— bitmap + pointer;  exkeynum*8 —— pageid(ungined int),offset(unsigned int); exkeynum*sizeof(key_type) ——exkey;
            double indexusge = double(vstats.leaves*leafnodesize+vstats.innernum*innodesize + 8 )/1048576;  //  leaf nodes, inner nodes, root
            double lruusage = double(no_hashsize*inkeynum + hashlrusize * vstats.leaves)/1048576;

            treeinfo.insert(pair<std::string,int>("leaves",vstats.leaves));
            treeinfo.insert(pair<std::string,int>("inners",vstats.innernum));
            treeinfo.insert(pair<std::string,int>("levels",vstats.numlevel));
            treeinfo.insert(pair<std::string,double>("datausage",datausage));
            treeinfo.insert(pair<std::string,double>("addusage",addusage));  //key,flag,pageID,offset
            treeinfo.insert(pair<std::string,double>("indexusage",indexusge));
            treeinfo.insert(pair<std::string,double>("lruusage",lruusage));
            return treeinfo;
        }


        const map<string, double> show_verify()
        {
            std::map<std::string,double> treeinfo;
            this->verify_film(&vstats);
            double innodesize = sizeof(key_type) + 8+8;  // startkey, double(slope,intercept)
            double leafnodesize = 2 * sizeof(key_type) + 8 + 8 ; //
            double hashlrusize = sizeof(key_type) + 8+ 16 + 16;
            double no_hashsize = 16+4;  //prev,next, slot-short int
            double dataV = (valuesize+1)*8;
            double datausage = double(inkeynum*(dataV))/1024/1024;
            double addusage = double((exkeynum+inkeynum)*(1+8)+exkeynum*(4+4)+exkeynum*sizeof(key_type) )/double(1048576);  // exkey, flag, pageID,offset   (pageid-int-8 offset short int), flag for data in memory

            double indexusge = (vstats.leaves*leafnodesize+vstats.innernum*innodesize + 8 )/double(1048576);  //  leaf nodes, inner nodes, root
            double lruusage = double(no_hashsize*inkeynum + hashlrusize * vstats.leaves)/1024/1024;

            treeinfo.insert(pair<std::string,int>("leaves",vstats.leaves));
            treeinfo.insert(pair<std::string,int>("inners",vstats.innernum));
            treeinfo.insert(pair<std::string,int>("levels",vstats.numlevel));
            treeinfo.insert(pair<std::string,double>("datausage",datausage));
            treeinfo.insert(pair<std::string,double>("addusage",addusage));  //key,flag,pageID,offset
            treeinfo.insert(pair<std::string,double>("indexusage",indexusge));
            treeinfo.insert(pair<std::string,double>("lruusage",lruusage));
            return treeinfo;
        }

        vector<int> traverse_leaves(){  // traverse the leaf level, find the min leaf and the max leaf
            leaflevel;
            int minslot = leaflevel.leafpieces[0]->slotkey.size();
            int maxslot = leaflevel.leafpieces[0]->slotkey.size();
//            for(auto& v : leaflevel){
//                if (v->slotkey.size() < minslot)
//                    minslot = v->slotkey.size();
//                if (v->slotkey.size() > maxslot)
//                    maxslot = v->slotkey.size();
//            }
            for (int i = 1; i<leaflevel.leafpieces.size()-1;i++){
                if (leaflevel.leafpieces[i]->slotkey.size() < minslot)
                    minslot = leaflevel.leafpieces[i]->slotkey.size();
                if (leaflevel.leafpieces[i]->slotkey.size() > maxslot)
                    maxslot = leaflevel.leafpieces[i]->slotkey.size();
            }
            int lastslot = leaflevel.leafpieces[leaflevel.leafpieces.size()-1]->slotkey.size();
            cout << " the last slot: " << lastslot << endl;
            vector<int> leafinfo;
            leafinfo.push_back(minslot);
            leafinfo.push_back(maxslot);
            leafinfo.push_back(lastslot);
            return leafinfo;
        }


    };

#pragma pack(push, 1)
    /**
     * struct in FILM
     */

    template<typename key_type,typename value_type>
    struct FILMinsert<key_type,value_type>::Leafpiece{

        key_type startkey;
        key_type endkey;
        vector<key_type> slotkey;
//        vector< pair<bool, adalru::Node<key_type, std::vector<key_type> >*  > >slotdata;
        vector< pair<bool, void* > > slotdata;   //adalru::Node<key_type, std::vector<key_type> >*
        adalru::localLRU <unsigned int,std::vector<key_type> >  intrachain;    //unsigned short int
        double slope;
        double intercept;

        Leafpiece() = default;

        Leafpiece(key_type startkey,key_type endkey,double slope, double intercept):
                startkey(startkey),endkey(endkey),slope(slope),intercept(intercept){};

        explicit Leafpiece(double slope, double intercept):
                startkey(),endkey(),slope(slope),intercept(intercept){};

        explicit Leafpiece(key_type endkey,double slope, double intercept):
                startkey(),endkey(endkey),slope(slope),intercept(intercept){};


        inline explicit Leafpiece(const typename internal::insertPWLF<key_type, int>::CanonicalSegment &cs):
                startkey(cs.get_first_x()),endkey(cs.get_last_x()){
            auto[cs_slope, cs_intercept] = cs.get_floating_point_segment(startkey);
            if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
                throw std::overflow_error("Change the type of Segment::intercept to int64");
            slope = cs_slope;
            intercept = cs_intercept;
            slotkey.assign(cs.slotkey.begin(),cs.slotkey.end());
        };

        friend inline bool operator<(const Leafpiece &p, const key_type &key) { return p.key < key; }
        friend inline bool operator<(const key_type &key, const Leafpiece &p) { return key < p.key; }
        friend inline bool operator<(const Leafpiece &p1, const Leafpiece &p2) { return p1.key < p2.key; }

        operator key_type() { return startkey; };

        /**
     * Returns the approximate position of the specified key.
     * @param k the key whose position must be approximated
     * @return the approximate position of the specified key
     */
        inline size_t operator()(const key_type &key) const {
            auto pos = int64_t(slope * (key - startkey)) + intercept;
            return pos > 0 ? size_t(pos) : 0ull;
        }

        inline void update(typename internal::insertPWLF<key_type, int>::CanonicalSegment &cs){
            startkey = cs.first;
            endkey = cs.last;
            auto[cs_slope, cs_intercept] = cs.get_floating_point_segment(startkey);
            if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
                throw std::overflow_error("Change the type of Segment::intercept to int64");
            slope = cs_slope;
            intercept = cs_intercept;
        }


        inline pair<bool,unsigned long> is_in_leaf(const key_type &key,const int error){
            auto pos = (*this)(key);
            pair<bool,unsigned long> is_in(false,0);
            int slotlo = PGM_SUB_EPS(pos, error + 1) ;  //,(--itlevel)->size()
            int slothi = PGM_ADD_EPS(pos, error,(this->slotkey.size()-1));
            for (; slotlo<=slothi;slotlo++){
                if (slotlo < 0)
                {
                    //cout<< "slotlo < 0 " << endl;
                }
                if (key == this->slotkey[slotlo])
                {   is_in.first = true;
                    is_in.second = slotlo;
                    break;}
            }

            return is_in;
        }

    };



    template<typename key_type,typename value_type>
    struct FILMinsert<key_type,value_type>::Innerpiece{

        key_type startkey;
        double slope;
        double intercept;

        inline Innerpiece() = default;

        inline Innerpiece(key_type startkey,double slope, double intercept):
                startkey(startkey),slope(slope),intercept(intercept){};

        inline explicit Innerpiece(double slope, double intercept):
                startkey(),slope(slope),intercept(intercept){};


        inline explicit Innerpiece(const typename internal::insertPWLF<key_type, int>::CanonicalSegment &cs):
                startkey(cs.get_first_x()){
            auto[cs_slope, cs_intercept] = cs.get_floating_point_segment(startkey);
            if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
                throw std::overflow_error("Change the type of Segment::intercept to int64");
            slope = cs_slope;
            intercept = cs_intercept;
        };

        friend inline bool operator<(const Innerpiece &p, const key_type &key) { return p.startkey < key; }
        friend inline bool operator<(const key_type &key, const Innerpiece &p) { return key < p.startkey; }
        friend inline bool operator<(const Innerpiece &p1, const Innerpiece &p2) { return p1.startkey < p2.startkey; }

        operator key_type() { return startkey; };

        /**
     * Returns the approximate position of the specified key.
     * @param k the key whose position must be approximated
     * @return the approximate position of the specified key
     */
        inline size_t operator()(const key_type &key) const {
            auto pos = int64_t(slope * (key - startkey)) + intercept;
            return pos > 0 ? size_t(pos) : 0ull;
        }

        inline void update(typename internal::insertPWLF<key_type, int>::CanonicalSegment &cs){
            startkey = cs.first;
            auto[cs_slope, cs_intercept] = cs.get_floating_point_segment(startkey);
            if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
                throw std::overflow_error("Change the type of Segment::intercept to int64");
            slope = cs_slope;
            intercept = cs_intercept;
        }
    };

    template<typename key_type,typename value_type>
    struct FILMinsert<key_type,value_type>::Innerlevel{
        std::vector<FILMinsert<key_type,value_type>::Innerpiece > innerpieces;
        vector<internal::insertPWLF<key_type, int>*> opt;
        unsigned int pos = 0;
        unsigned int nextpos = 0;
    };

    template<typename key_type,typename value_type>
    struct FILMinsert<key_type,value_type>::Leaflevel{
        std::vector<FILMinsert<key_type,value_type>::Leafpiece* > leafpieces;
        vector<internal::insertPWLF<key_type, int>*> opt;
        unsigned int pos = 0;
    };

#pragma pack(pop)
    typedef long int key_type;    // books or wiki_ts
//    typedef double key_type;
    typedef filminsert::FILMinsert< key_type, vector<key_type> > filmadatype;

    typedef filminsert::FILMinsert< key_type, vector<key_type> > filmadatype;
    typedef pair<filmadatype::Leafpiece*,unsigned short int>  filmadalrupair;
    typedef  adalru::hashLRU <key_type,filmadatype::Leafpiece* , adalru::Node<key_type ,filmadatype::Leafpiece*>* > filmadalrutype;
    typedef adalru::localLRU <unsigned short,std::vector<key_type> > locallrutype;
    typedef filmstorage::filmdisk<key_type> filmadadisk;
    typedef filmstorage::filmmemory<key_type,vector<key_type>,filmadatype, filmadalrutype,filmadadisk> filmadamemory;
    typedef filmadamemory::memoryusage memoryusage;


    void runtimeevictkeytopage(filmadamemory *filmmem,filmadatype *filmindex,filmadadisk *diskpage , access_stats *r_stats){
        struct timeval lt1, lt2,dt1,dt2,rdt1,rdt2,wdt1,wdt2;
        double ltimeuse,rdtimeuse,wdtimeuse;
        double dtimeuse;
        if (filmmem->evictPoss.size() != 0){
            for (int i = 0; i< filmmem->evictPoss.size(); i++){
                gettimeofday(&lt1, NULL);
                auto writeevict =  filmmem->evictPoss[i];

                auto evictleaf = filmmem->lru->get_tail();
                auto evictslotV = evictleaf->intrachain.poptail();
                gettimeofday(&lt2, NULL);
                ltimeuse = (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                r_stats->lrutime += ltimeuse;

                filmmem->evictkeytoinpage(evictleaf->slotkey[evictslotV->key], evictslotV->value,
                                          diskpage, writeevict);
                evictleaf->slotdata[evictslotV->key].first = false;
                evictleaf->slotdata[evictslotV->key].second = writeevict;
                delete evictslotV;
                evictslotV = NULL;

            }
            filmindex->inkeynum -= filmmem->evictPoss.size();
            filmindex->exkeynum += filmmem->evictPoss.size();

            gettimeofday(&wdt1, NULL);
            int pagenum = filmmem->runtimeevictpagestodisk(diskpage);
            r_stats->writenum += pagenum;
            gettimeofday(&wdt2, NULL);
            wdtimeuse = (wdt2.tv_sec - wdt1.tv_sec) + (double) (wdt2.tv_usec - wdt1.tv_usec) / 1000000.0;
            r_stats->wdisktime += wdtimeuse;
        }
    }

    // excute point query
    pair<access_stats , vector<vector<key_type>> > filmadapointquery(filmadatype *filmindex,vector<key_type> pqueries, int queryn,filmadadisk *pagedisk , filmadamemory *filmmem){
        struct timeval qt1, qt2;
        double timeuse;
        vector<vector<key_type>> totalres;
        access_stats point_stats;
        bool transflag = filmmem->runtimejudgetrans(point_stats);
        int periodV = filmmem->reserveMem*1024*1024/((filmindex->valuesize+1)*sizeof(key_type)*10);

        gettimeofday(&qt1, NULL);
        for (unsigned int i = 0; i < pqueries.size(); i++) {
            if (point_stats.disknum >0 && i % periodV == 0  ){
                transflag = filmmem->runtimejudgetrans(point_stats);
                while (transflag)
                {
                    runtimeevictkeytopage(filmmem,filmindex,pagedisk,&point_stats);
                    transflag = filmmem->runtimejudgetrans(point_stats);
                }

            }

            vector<key_type> res;
            key_type querykey = pqueries[i];

            auto index_res = filmindex->search_one(querykey);

            if ( index_res.find == false){
                if (index_res.flags){
                    point_stats.memnum += 1 ;
                    res.emplace_back(querykey);
                    auto finddata = (adalru::Node<unsigned int, std::vector<key_type> >*) filmindex->sort_list->slotdata[index_res.slot].second;
                    res.insert(res.begin()+1,finddata->value.begin(),finddata->value.end());
                    filmindex->sort_list->intrachain.moveTohead(finddata);// update intrachain
                    //                totalres.emplace_back(res);

                    filmmem->lru->put(filmindex->sort_list->slotkey[0],filmindex->sort_list);
                }
                else{
                    point_stats.disknum += 1;
                    point_stats.diskpagenum += 1;
                    auto writeevict = (pair<unsigned int,unsigned int>*)index_res.findleaf->slotdata[index_res.slot].second;  //writeevict 指向的是 pageid and offset

                    res = pagedisk->odirectreadfromdisk(writeevict);    // if readfromdisk indicating doesn't use o_direct;
                    //                totalres.push_back(res);
                    res.erase(res.begin());
                    index_res.findleaf->slotdata[index_res.slot] = index_res.findleaf->intrachain.put(index_res.slot,res);

                    filmmem->lru->put(index_res.findleaf->startkey,index_res.findleaf);
                    filmmem->evictPoss.emplace_back(writeevict);
                    filmindex->inkeynum++;
                    filmindex->exkeynum--;

                }

            }
            else{
                if (index_res.flags){
                    point_stats.memnum += 1;
                    auto finddata = (adalru::Node<unsigned int, std::vector<key_type> >*) index_res.findleaf->slotdata[index_res.slot].second;
                    res.emplace_back(querykey);
                    res.insert(res.begin()+1,finddata->value.begin(),finddata->value.end());
                    index_res.findleaf->intrachain.moveTohead(finddata);// update intrachain
                    //                totalres.emplace_back(res);

                    filmmem->lru->put(index_res.findleaf->startkey,index_res.findleaf);
                }
                else{
                    point_stats.disknum += 1;
                    point_stats.diskpagenum += 1;
                    auto writeevict = (pair<unsigned int,unsigned int>*)index_res.findleaf->slotdata[index_res.slot].second;

                    res = pagedisk->odirectreadfromdisk(writeevict);    // if readfromdisk indicating doesn't use o_direct;
                    //                totalres.push_back(res);
                    res.erase(res.begin());
                    index_res.findleaf->slotdata[index_res.slot] = index_res.findleaf->intrachain.put(index_res.slot,res);

                    filmmem->lru->put(index_res.findleaf->startkey,index_res.findleaf);
                    filmmem->evictPoss.emplace_back(writeevict);
                    filmindex->inkeynum++;
                    filmindex->exkeynum--;
                }
            }

        }

        gettimeofday(&qt2, NULL);
        //cout<<"my Lord, You are my refuge~~~~forever! "<<endl;
        timeuse = (qt2.tv_sec - qt1.tv_sec) + (double) (qt2.tv_usec - qt1.tv_usec) / 1000000.0;
        point_stats.querytimeuse += timeuse;
        point_stats.querytimeuse -= point_stats.computetimeuse;
        return pair<access_stats, vector<vector<key_type>>>(point_stats,totalres);
    }



    int test_filmappending(unsigned int errbnd, size_t datanum, int pagesize, string dataset,double mem_threshold,double reserveMem,int recordSize,vector<key_type> keys,vector<key_type> queries, int queryn, size_t numkey,int datasize,double zipf){
        filmadatype filminsert(datanum,errbnd,errbnd);
        filmadalrutype interchain(numkey);
        filminsert.valuesize = recordSize-1;
        std::ostringstream osse,ossr,ossd,ossm,ossp;
        osse << errbnd;   ossp << pagesize;  ossm << mem_threshold;    ossr << (recordSize);  ossd << datanum;
        string diskpath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/diskpath/";
        string file_str = diskpath+ dataset +"_filminsertpages"+ossd.str()+"_"+ ossm.str()+"_"+ossp.str()+"_"+ossr.str()+"_"+osse.str();
        const char* diskfile = file_str.c_str() ;
        int numrecord = pagesize/(recordSize);
        filmstorage::filmmemory<key_type,vector<key_type>,filmadatype, filmadalrutype,filmadadisk> memoryfilm(numkey,mem_threshold,&filminsert,&interchain) ;
        filmstorage::filmdisk<key_type> diskfilm(diskfile,pagesize,numrecord,recordSize);
        memoryfilm.reserveMem = reserveMem;
        fstream fs;
        fs.open(diskfile,ios::in);
        if (fs){
            remove(diskfile);
        }

        int fd = -1;
        int ret = -1;
        uint64_t file_size = 2*datasize*1024*1024ULL;

        fd = open(diskfile, O_CREAT|O_RDWR, 0666);
        if(fd < 0){
            printf("fd < 0");
            return -1;
        }

        //ret = fallocate(fd, 0, 0, file_size);
        ret = posix_fallocate(fd, 0, file_size);
        if(ret < 0 ){
            printf("ret = %d, errno = %d,  %s\n", ret, errno, strerror(errno));
            return -1;
        }

        printf("fallocate create %.2fG file\n", file_size/1024/1024/1024.0);
        close(fd);


        vector<key_type> payload(filminsert.valuesize);

        struct timeval bt1, bt2;
        double buildtimeuse;
        gettimeofday(&bt1, NULL);
        memoryfilm.append(keys,payload,errbnd,errbnd);
        gettimeofday(&bt2, NULL);
        buildtimeuse = (bt2.tv_sec - bt1.tv_sec) + (double) (bt2.tv_usec - bt1.tv_usec) / 1000000.0;
        //cout << "insert time use = " << buildtimeuse << " , thank You, my Lord " << endl;
        ofstream savefile;
        savefile.open("/home/wamdm/chaohong/clionDir/Opt_FILMupdate/result/adalru_film_performance.txt",ios::app);
        savefile<< "method " << "film_ada_lru "<<"available_memory " << mem_threshold << " qtype "<< "point " << "error " << errbnd  << " pagesize "<< (pagesize*8/1024) << "k recordsize ";
        savefile<< recordSize << " build_time " << buildtimeuse << " dataset " << dataset <<" datasize " << datasize << " keynum "<<numkey << " ";
        savefile << "\n";
        savefile << flush;
        savefile.close();


//      show the statistic information of film
        auto filminfo = filminsert.show_verify();
        map<string, double>::reverse_iterator   iter2;
        for(iter2 = filminfo.rbegin(); iter2 != filminfo.rend(); iter2++){
            cout<<iter2->first<<" "<<iter2->second<<" *** ";
        }
        cout<< endl;
        auto leafinfo = filminsert.traverse_leaves();
        cout<< dataset << "  " ;
        for (int i = 0; i< leafinfo.size(); i ++){
            cout<< leafinfo[i]<< "  ";
        }
        cout << endl;

        memoryfilm.lru->capacity = memoryfilm.lru->size;

        pair<bool,memoryusage> transflag = memoryfilm.judgetransfer();
        int transleaves;
        double ratio;
        int leaves = transflag.second.meminfo["leaves"];
        while (transflag.first)
        {
            if (filminsert.leafsplit == 0){
                ratio = memoryfilm.threshold/transflag.second.totalusemem;
                transleaves = (leaves - filminsert.leafsplit) * (1-ratio) ;  // the transfer number of leafs
            }
            else{
                ratio = memoryfilm.threshold/transflag.second.totalusemem;
                transleaves = (leaves - filminsert.leafsplit) * (1-ratio);
                if (transleaves ==0){
                    if (errbnd>=32)
                        transleaves = 1;
                    else
                        transleaves = 3;                }
            }
            memoryfilm.filmtransfer(transleaves,&diskfilm);
            transflag = memoryfilm.judgetransfer();
        }
        cout<< "finish data transfer"<< endl;
        savefile.open("/home/wamdm/chaohong/clionDir/Opt_FILMupdate/result/adalru_film_performance.txt",ios::app);
        map<string,double>::iterator iter;
        savefile<< "method " << "film_ada_lru ";
        for(iter = transflag.second.meminfo.begin(); iter != transflag.second.meminfo.end(); iter++)
            savefile<<iter->first<<" "<<iter->second<<" ";
        savefile <<"\n";
        savefile << flush;
        savefile.close();

        pair<access_stats, vector<vector<key_type>>> point_performance = filmadapointquery(&filminsert, queries, queryn,&diskfilm,&memoryfilm);
        point_performance.first.zipffactor = zipf;
        cout<< " finished the point query test with film-adaptivelru, the performance is:"<<endl;
        point_performance.first.print_stats();
        point_performance.second.clear();
        for (int i = 0;i<memoryfilm.index->leaflevel.leafpieces.size();i++){
            delete memoryfilm.index->leaflevel.leafpieces[i];
            memoryfilm.index->leaflevel.leafpieces[i] = NULL;
        }

        fs.open(diskfile,ios::in);
        if (fs){
            remove(diskfile);
        }
        return 0;
    }


    int test_filmupdating(unsigned int errbnd, size_t datanum, int pagesize, string dataset,double mem_threshold,double reserveMem,int recordSize,vector<key_type> keys,vector<key_type> update_keys,vector<key_type> queries, int queryn, size_t numkey,int datasize,double zipf,double up_ratio){
        filmadatype filmupdate(datanum,errbnd,errbnd);
        filmadalrutype interchain(numkey);
        filmupdate.valuesize = recordSize-1;
        std::ostringstream osse,ossr,ossd,ossm,ossp;
        osse << errbnd;   ossp << pagesize;  ossm << mem_threshold;    ossr << (recordSize);  ossd << datanum;
        string diskpath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/diskpath/";
        string file_str = diskpath+ dataset +"_filminsertpages"+ossd.str()+"_"+ ossm.str()+"_"+ossp.str()+"_"+ossr.str()+"_"+osse.str();
        const char* diskfile = file_str.c_str() ;
        int numrecord = pagesize/(recordSize);
        filmstorage::filmmemory<key_type,vector<key_type>,filmadatype, filmadalrutype,filmadadisk> memoryfilm(numkey,mem_threshold,&filmupdate,&interchain) ;
        filmstorage::filmdisk<key_type> diskfilm(diskfile,pagesize,numrecord,recordSize);
        memoryfilm.reserveMem = reserveMem;
        fstream fs;
        fs.open(diskfile,ios::in);
        if (fs){
            remove(diskfile);
        }

        int fd = -1;
        int ret = -1;
        uint64_t file_size = 2*datasize*1024*1024ULL;

        fd = open(diskfile, O_CREAT|O_RDWR, 0666);
        if(fd < 0){
            printf("fd < 0");
            return -1;
        }

        //ret = fallocate(fd, 0, 0, file_size);
        ret = posix_fallocate(fd, 0, file_size);
        if(ret < 0 ){
            printf("ret = %d, errno = %d,  %s\n", ret, errno, strerror(errno));
            return -1;
        }

        printf("fallocate create %.2fG file\n", file_size/1024/1024/1024.0);
        close(fd);

        vector<key_type> payload(filmupdate.valuesize);

        struct timeval bt1, bt2;
        double buildtimeuse;

        memoryfilm.append(keys,payload,errbnd,errbnd);
        gettimeofday(&bt1, NULL);
        memoryfilm.update(update_keys,payload);
        gettimeofday(&bt2, NULL);
        buildtimeuse = (bt2.tv_sec - bt1.tv_sec) + (double) (bt2.tv_usec - bt1.tv_usec) / 1000000.0;
        cout << "insert time use = " << buildtimeuse << " , thank You, my Lord " << endl;
        ofstream savefile;
        savefile.open("/home/wamdm/chaohong/clionDir/Opt_FILMupdate/result/adalru_film_performance.txt",ios::app);
        savefile<< "method " << "film_ada_lru "<<"available_memory " << mem_threshold << " qtype "<< "point " << "error " << errbnd  << " pagesize "<< (pagesize*8/1024) << "k recordsize ";
        savefile<< recordSize << " build_time " << buildtimeuse << " dataset " << dataset <<" datasize " << datasize << " keynum "<<numkey << " ";
        savefile << "\n";
        savefile << flush;
        savefile.close();


//      show the statistic information of film
        auto filminfo = filmupdate.show_verify();
        map<string, double>::reverse_iterator   iter2;
        for(iter2 = filminfo.rbegin(); iter2 != filminfo.rend(); iter2++){
            cout<<iter2->first<<" "<<iter2->second<<" *** ";
        }
        cout<< endl;
        auto leafinfo = filmupdate.traverse_leaves();
        cout<< dataset << "  " ;
        for (int i = 0; i< leafinfo.size(); i ++){
            cout<< leafinfo[i]<< "  ";
        }
        cout << endl;


        memoryfilm.lru->capacity = memoryfilm.lru->size;

        pair<bool,memoryusage> transflag = memoryfilm.judgetransfer();
        int transleaves;
        double ratio;
        int leaves = transflag.second.meminfo["leaves"];
        while (transflag.first)
        {
            if (filmupdate.leafsplit == 0){
                ratio = memoryfilm.threshold/transflag.second.totalusemem;
                transleaves = (leaves - filmupdate.leafsplit) * (1-ratio) ;  // the transfer number of leafs
            }
            else{
                ratio = memoryfilm.threshold/transflag.second.totalusemem;
                transleaves = (leaves - filmupdate.leafsplit) * (1-ratio);
                if (transleaves ==0){
                    if (errbnd>=32)
                        transleaves = 1;
                    else
                        transleaves = 3;                }
                //cout<< "Jesus, please help me!"<< endl;
            }
            memoryfilm.filmtransfer(transleaves,&diskfilm);
            transflag = memoryfilm.judgetransfer();
        }
        cout<< "finish data transfer"<< endl;
        savefile.open("/home/wamdm/chaohong/clionDir/Opt_FILMupdate/result/adalru_film_performance.txt",ios::app);
        map<string,double>::iterator iter;
        savefile<< "method " << "film_ada_lru ";
        for(iter = transflag.second.meminfo.begin(); iter != transflag.second.meminfo.end(); iter++)
            savefile<<iter->first<<" "<<iter->second<<" ";
        savefile <<"\n";
        savefile << flush;
        savefile.close();

        pair<access_stats, vector<vector<key_type>>> point_performance = filmadapointquery(&filmupdate, queries, queryn,&diskfilm,&memoryfilm);
        point_performance.first.zipffactor = zipf;
        cout<< " finished the point query test with film-adaptivelru, the performance is:"<<endl;
        point_performance.first.print_stats();
        point_performance.second.clear();
        for (int i = 0;i<memoryfilm.index->leaflevel.leafpieces.size();i++){
            delete memoryfilm.index->leaflevel.leafpieces[i];
            memoryfilm.index->leaflevel.leafpieces[i] = NULL;
        }

        fs.open(diskfile,ios::in);
        if (fs){
            remove(diskfile);
        }
        return 0;
    }

}


#endif //FILMINSERT_FILM_H
