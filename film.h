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


#include "pwlf.h"
#include "filmadastorage.h"
#include "filmadalru.h"
#include "zipf.h"

//typedef long int key_type;  // if (filename== "books" || filename == "wiki_ts", "synthetic", "YCSB")
typedef double key_type;   //if (filename== "astro_ra")

#define forceinline inline __attribute__((__always_inline__))
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
        double xlookuptime = 0.0;
        double lrutime = 0.0;
        double disktime = 0.0;
        double rdisktime = 0.0;
        double wdisktime = 0.0;
        double sort_piecetime = 0.0;
        unsigned int pagecross = 0;  // 磁盘页 seek 的跨度
        double zipffactor = 0.0;
        string qtype = "point";
        string workload = "hotspot";  // "zipf", "random", "hotspot","randomzipf"
        double insert_frac = 0.0;
        double out_of_order_frac = 0.0;

        void print_stats(){
            struct timeval t1;
            gettimeofday(&t1, NULL);
            cout<<"ID " <<double(t1.tv_sec) << " qtype "<< qtype << " zipffactor "<< zipffactor << " insert_frac "<< insert_frac << " out_of_order_frac "<< out_of_order_frac << " memnum " << memnum << " disknum " << disknum << " crossnum "<< crossnum << " diskpagenum ";
            cout<< diskpagenum<<" crosspagenum " << crosspagenum <<" pagecross " << pagecross<<" writenum " << writenum<< " querytimeuse " <<querytimeuse <<" lrutimeuse " <<lrutime <<" disktimeuse " << disktime <<" readdisktimeuse " <<rdisktime
            <<" writedisktimeuse " <<wdisktime << " indexloopuptimeuse " << xlookuptime << " sort_piecetime "  << sort_piecetime << " workload " << workload <<"\n";
            ofstream savefile;
            savefile.open("/home/wamdm/chaohong/clionDir/updatefilm/result/adalru_film_performance.txt",ios::app);
//            savefile << "_ ";
            savefile<<"ID " <<double(t1.tv_sec) << " qtype "<< qtype  << " zipffactor "<< zipffactor << " insert_frac "<< insert_frac  << " out_of_order_frac "<< out_of_order_frac  << " memnum " << memnum << " disknum " << disknum << " crossnum "<< crossnum << " diskpagenum ";
            savefile<< diskpagenum<<" crosspagenum " << crosspagenum <<" pagecross " << pagecross<<" writenum " << writenum<<" querytimeuse " <<querytimeuse  <<" lrutimeuse " <<lrutime <<" disktimeuse " <<disktime << " readdisktimeuse " <<rdisktime
            <<" writedisktimeuse " <<wdisktime << " indexloopuptimeuse " << xlookuptime<< " sort_piecetime "  << sort_piecetime << " workload " << workload<<"\n";
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

        struct Leafpiece;  // 声明leafpiece 结构体
        struct sortPiece:Leafpiece{

            vector<key_type> slotkey;
//        vector< pair<bool, adalru::Node<key_type, std::vector<key_type> >*  > >slotdata;
            vector< bool> slotflag;   //adalru::Node<key_type, std::vector<key_type> >*sss
            vector< void* > slotdata;   //adalru::Node<key_type, std::vector<key_type> >*sss
            adalru::localLRU <lruOff_type,key_type* >  intrachain;    //unsigned short int

            sortPiece() = default;

            friend inline bool operator<(const sortPiece &p, const key_type &key) { return p.key < key; }
            friend inline bool operator<(const key_type &key, const sortPiece &p) { return key < p.key; }
            friend inline bool operator<(const sortPiece &p1, const sortPiece &p2) { return p1.key < p2.key; }


            /**
         * Returns the approximate position of the specified key.
         * @param k the key whose position must be approximated
         * @return the approximate position of the specified key
         */
            forceinline int operator()(const key_type key)  const{
                int begin = 0, end = slotkey.size()-1,mid;  // find last Equal or small, 找到最后一个大于或等于 key 的
                while(begin <= end){
                    mid = (end + begin) / 2;
                    if(slotkey[mid] >= key) {
                        end = mid-1;
                    } else
                        begin = mid +1;
                }
                return end+1;
            }


            forceinline bool insert(key_type key, key_type* payload) {

                int begin = 0, end = slotkey.size()-1,mid;  // find last Equal or small, 找到最后一个大于或等于 key 的
                while(begin <= end){
                    mid = (end + begin) / 2;
                    if(this->slotkey[mid] >= key) {
                        end = mid-1;
                    } else
                        begin = mid +1;
                }
                // 遍历 intrachain，所有 大于 end 的 值，都需要加1
                intrachain.modify(end);  // 所有 在 end 之后的，都需要加1
                this->slotdata.insert( slotdata.begin()+end+1,intrachain.put(end+1,payload));
                this->slotkey.insert(slotkey.begin()+end+1,key);
                this->slotflag.insert(slotflag.begin()+end+1,true);
                return true;
            }


            // 判断key 是否在当前leaf 中
            inline pair<bool,unsigned long> is_in_leaf(const key_type &key,const int error){
                pair<bool,unsigned long> is_in(false,0);
                if (key<=this->endkey ){
                    auto pos = (*this)(key);

                    int slotlo = PGM_SUB_EPS(pos, error + 1) ;  //,(--itlevel)->size()
                    int slothi = PGM_ADD_EPS(pos, error,(this->slotkey.size()-1));
                    for (; slotlo<=slothi;slotlo++){

                        if (key == this->slotkey[slotlo])
                        {   is_in.first = true;
                            is_in.second = slotlo;
                            break;}
                    }
//            cout<< " happy new year! Jesus~~~~" << endl;
                    return is_in;
                }
                else
                    return is_in;
            }

        };
        struct Innerpiece;  // 声明leafpiece 结构体
        struct Innerlevel;  // 声明 innerlevel 结构体，每个innerllevel 包含一个 vector<Innerpiece>,一个 opt， 一个 pos；
        struct Leaflevel;
        unsigned int totalnum;
        unsigned int inkeynum=0;
        unsigned int exkeynum=0;
        int valuesize = 0;
        int Error ;
        int ErrorRecursive ;
        int leafsplit = 0;
        Leafpiece *m_transleaf = NULL;
        Leafpiece *m_tailleaf = NULL;   // tailpiece, 用于append increasing keys
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

        ~FILMinsert() {
            // delete leaf piece
        }

        void delete_leaf(Leafpiece* node) {
            if (node == NULL) {
                return;
            } else{
                vector<key_type>().swap(node->slotkey);
                vector<bool>().swap(node->slotflag);
                vector<void*>().swap(node->slotdata);
                // release intrachain

                node->intrachain.deletelru();
//                for (int i = 0; i < node->intrachain.size; i ++){
//                    delete node->intrachain[i];
//                }
//                adalru::localLRU <lruOff_type,key_type* >  intrachain;    //unsigned short int
                delete node;
                node = NULL;
//                data_node_allocator().deallocate(static_cast<data_node_type*>(node), 1);
            }
        }

        void delete_inner(Innerpiece node) {

            delete node;
            node = nullptr;
//                data_node_allocator().deallocate(static_cast<data_node_type*>(node), 1);
        }

        void release(){
            for ( unsigned int leafi = 0; leafi < leaflevel.leafpieces.size();leafi++){
                auto leaf_it = leaflevel.leafpieces[leafi];
                delete_leaf(leaf_it);
            }
            vector<FILMinsert<key_type,value_type>::Leafpiece*>().swap(leaflevel.leafpieces);
            delete leaflevel.opt;
            for (unsigned int leveli = 0; leveli <innerlevels.size();leveli ++ ){
                auto level_it = innerlevels[leveli];
//                for (unsigned int inneri = 0; inneri < level_it->innerpieces.size();inneri ++ ){
//                    auto inner_it = level_it->innerpieces[inneri];
//                    delete_inner(inner_it);
//                }
//                level_it->innerpieces.clear();
                vector<Innerpiece>().swap(level_it->innerpieces);
                delete level_it->opt;
//                level_it->opt.clear();
            }
            innerlevels.clear();
            inkeynum = 0;
            exkeynum = 0;
            m_transleaf = NULL;
            m_tailleaf = NULL;
            if (sort_list != NULL){
                sort_list = NULL;
            }


        }

        /*
         * the method to construct film, that build FILM
         */


        inline std::vector<key_type>  append_one(size_t error,std::vector<key_type> keys,unsigned int k){

            std::vector<key_type> startkeys;

            for (size_t i = 0; i < keys.size(); ++i) {
                pair<key_type,unsigned int> p(keys[i],innerlevels[k]->nextpos++) ;   // i 为 pos
                ++(innerlevels[k]->pos);
                if (!innerlevels[k]->opt->add_point(p.first, p.second)) {  // 如果inner level  不满足error 了，那么再创建一个innerpiece
                    // 当前 innerpiece 不再满足，需要创建new inner piece 并判断该 innerlevel 的上一层level 是否需要更新
                    auto a = innerlevels[k]->opt->get_segment();

                    if (innerlevels[k]->pos > 2)
                        innerlevels[k]->innerpieces.pop_back();
                    innerlevels[k]->innerpieces.emplace_back(a);
//                    auto newinner = innerlevels[k]->innerpieces.back();
//                    cout<< "i need You, my lovely Lord, why there is nan" << endl;

                    // 首先在该层创建一个 new innerpiece， 更新该innerpiece，再递归向向上
//                    typename filmtype::Innerpiece *innerpiece = new typename filmtype::;
                    innerlevels[k]->pos = 0;
                    innerlevels[k]->nextpos -= 2;
                    delete innerlevels[k]->opt;
//                    innerlevels[k]->opt.pop_back();
                    internal::insertPWLF<key_type, int> *inneropt = new internal::insertPWLF<key_type, int>(error);
                    innerlevels[k]->opt = inneropt;
                    if (k==0){
//                        auto aaaaa = leaflevel.leafpieces.size()-2;
                        startkeys.emplace_back(leaflevel.leafpieces[leaflevel.leafpieces.size()-1]->startkey);
                    }

                    else{
//                        auto aaaaa = innerlevels[k-1]->innerpieces.size()-2;
                        startkeys.emplace_back(innerlevels[k-1]->innerpieces[innerlevels[k-1]->innerpieces.size()-2].startkey);
                    }
                    startkeys.emplace_back(p.first);
                    auto rr = append_one(error,startkeys,k);


                    if (innerlevels.back()->innerpieces.size() > 1)
                    {
                        startkeys.clear();
                        startkeys.emplace_back(a.first);
                        startkeys.emplace_back(p.first);
                        Innerpiece innerpiece;// 创建parent piece
                        internal::insertPWLF<key_type, int> *inneropt = new internal::insertPWLF<key_type, int>(error);
//                    std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> > *innerlevel =
//                            new std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> >;
                        Innerlevel *innerlevel = new Innerlevel;
                        innerlevels.emplace_back(innerlevel);
                        innerlevels.back()->opt = inneropt ;
                        auto rr = append_one(error,startkeys,k+1);
//                    cout<< "Jesus, i need You !"  << endl;

                        return startkeys ;
                    }
                    else if (innerlevels.size() > 1 && (k != innerlevels.size()-1) )   // 如上为 由于创建了new innner piece，导致了new innerlevel，如下为，虽然创建了new innerpiece，但只需要更新上层的innner level
                    {
                        // 更新上层的最后一个inner piece
                        startkeys.pop_back();
                        auto rr = append_one(error,startkeys,k+1);
                        startkeys.clear();
//                    cout << "thank You, my Lord! i need You!"<<endl;
//                        startkeys.emplace_back(a.first);
                        return startkeys ;
                    }
                    else if (innerlevels.size() > 1 && (k == innerlevels.size()-1) )   // 如上为 由于创建了new innner piece，导致了new innerlevel，如下为，虽然创建了new innerpiece，但只需要更新上层的innner level
                    {
                        cout << "thank You, my Lord! i need You!"<<endl;
                        startkeys.clear();
                        return startkeys ;
                    }
                    else{
                        startkeys.clear();
                    }
                    a = innerlevels[k]->opt->get_segment();
                    if (innerlevels[k]->pos > 2)
                        innerlevels[k]->innerpieces.pop_back();
                    innerlevels[k]->innerpieces.emplace_back(a);

                    startkeys.emplace_back(a.first);
                    return startkeys;
                }

            }

            auto a = innerlevels[k]->opt->get_segment();
            if (innerlevels[k]->pos > 2)
                innerlevels[k]->innerpieces.pop_back();
            innerlevels[k]->innerpieces.emplace_back(a);

            startkeys.emplace_back(a.first);
            return startkeys ;

        }



        template<class lru_type>   // insert one key on an existing index
        forceinline void append_one(key_type key,key_type* payload,unsigned int error,lru_type interchain){
            std::vector<key_type> startkeys;
            inkeynum++;
            pair<key_type,unsigned int> p(key,leaflevel.opt->points_in_hull) ;

            if (leaflevel.opt->append_point(p.first, p.second)) {
                m_tailleaf->slotkey.emplace_back(p.first);
                m_tailleaf->endkey = p.first;
                m_tailleaf->slotdata.emplace_back( m_tailleaf->intrachain.put(p.second,payload));
                m_tailleaf->slotflag.emplace_back( true);
            }
            else
            {

                auto a = leaflevel.opt->get_segment(m_tailleaf->endkey); // 将生成的new leaf piece插入到leaflevel 中
                m_tailleaf->update(a);
                //leaflevel.leafpieces.emplace_back( m_tailleaf);
                interchain->put(m_tailleaf->startkey,m_tailleaf);
                // 这里是初始化 parent piece 的first key 和 last key
                if (innerlevels.size() == 0){
                    startkeys.emplace_back( m_tailleaf->startkey);
                    startkeys.emplace_back( p.first);
                    Innerpiece innerpiece;// 创建parent piece
                    internal::insertPWLF<key_type, int> *inneropt = new internal::insertPWLF<key_type, int>(error);
//                    std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> > *innerlevel =
//                            new std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> >;
                    Innerlevel *innerlevel = new Innerlevel;
                    innerlevels.emplace_back(innerlevel);
                    innerlevels[0]->opt = inneropt ;
                    cout << " my Lord, Jesus, please have pity on me"<< endl;
                }
                else{
                    // 从innerlevel 的最底层到root 层，判断是否需要更新
                    startkeys.emplace_back( p.first);  //p.first is the break_key, and also the startkey of the next leaf piece
//                    cout << "my lovely Lord, i trust in You!" << endl;
                }
                // use the
                auto rr = append_one(error,startkeys,0);
                startkeys.clear();
//                cout << "Jesus, i need You!!"<< endl;

                m_tailleaf = new Leafpiece;
                m_tailleaf->slotkey.reserve(8192*2);
                leaflevel.opt->append_point(p.first, 0);
                m_tailleaf->slotdata.emplace_back( m_tailleaf->intrachain.put(0,payload));
                m_tailleaf->slotkey.push_back(p.first);
                m_tailleaf->slotflag.emplace_back( true);
                m_tailleaf->startkey = p.first;
                interchain->put(m_tailleaf->startkey,m_tailleaf);
                leaflevel.leafpieces.emplace_back( m_tailleaf);
//                auto newleaf = leaflevel.leafpieces.back();
//                cout<< "i need You, my lovely Lord!" << endl;

            }
        }




        template<class lru_type>  // insert_bulk_load_one_by_one (also one by one), in the beginning, the index must be empty? maybe not  empty is okay
        inline void update_append(std::vector<key_type> keys, key_type* payload,unsigned int error,unsigned int error_recursize,lru_type &interchain ){
            // build leaf level
            // size_t error,std::vector<key_type> keys,std::vector<key_type> payload,
            //filmtype *filmada, unsigned int &inkeynum,leaf_type* m_taillea
            std::pair<size_t,std::vector<key_type>>  n_parts = internal::append_segmentation(error,keys,payload,this,inkeynum,this->m_tailleaf,interchain);   //<key_type,Leafpiece>
//            std::cout<< "You are my refuge, from ever to ever"<<endl;
            // on the basis of leaf level, build the upper level until just has one piece, that the make_segmentation just return 1.
            if (innerlevels.size()>1)
                root = &innerlevels.back()->innerpieces[0];
            this->verify_film(&vstats);
            cout<<"index append test finished,Jesus, my Lord, i need You, praise to You for ever and ever!" << endl;
        }

        /*
        inline void update_random(std::vector<key_type> keys,std::vector<key_type> payload){
            // 首先search , 如果不存在，则插入，如果存在，则 pass
            for (unsigned int i = 0; i < keys.size(); i ++){
                this->update_search_one(keys[i],payload);
            }
            this->inkeynum += keys.size();
        }
        */
        template<class lru_type>
        inline void update_random(std::vector<key_type> keys,key_type* payload,lru_type* interchain){
            // 首先search , 如果不存在，则插入，如果存在，则 pass
            unsigned int i = 0;
            if (sort_list == NULL)
            {
                sort_list = new sortPiece();
                this->sort_list->slotkey.emplace_back(keys[i++]);
                this->sort_list->slotflag.emplace_back(true);
                this->sort_list->slotdata.emplace_back( sort_list->intrachain.put(0,payload));
            }
            for (; i < keys.size(); i ++){
                this->update_one(keys[i],payload,interchain);
            }
            this->inkeynum += keys.size();
            // add sort_list into interchain
            interchain->put((*sort_list).slotkey[0],sort_list);
        }

        //  与 range query 相关的返回结果
        struct range_res_find
        {
            /// find result flags
            bool  flags;

            /// The key to be found
            key_type    firstkey;
            key_type    lastkey;

            /// The leaf nodes that overlap with the queried range
            vector<Leafpiece *> findleafs;

            /// the location of key in leaf node
            unsigned int firstslot;
            unsigned int lastslot;

            /// Constructor of a result with a specific flag, this can also be used
            /// as for implicit conversion.

            forceinline range_res_find(bool f)
                    : flags(f), firstkey(),lastkey(), findleafs(),firstslot(),lastslot()
            { }

            forceinline range_res_find(bool f,key_type firstkey,key_type lastkey)
                    : flags(f), firstkey(firstkey), lastkey(), findleafs(),firstslot(),lastslot()
            { }


            /// Constructor with a lastkey value.
            forceinline range_res_find(bool f ,key_type firstkey,int first_loc,key_type lastkey)
                    : flags(f), firstkey(firstkey), findleafs(),firstslot(first_loc),lastkey(lastkey),lastslot()
            { }

            inline bool in_or_out(){
                return true;
            }

        };

        /// film recursive find the key has much information which is needs to be return.
        struct result_find
        {
            /// find result flags
            bool find;   // 在sorted list 还是在 regular data

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

                if (findleaf->slotflag[slot])
                    return true;
                else
                    return false;
            }
        };


        /** search method
         * the methods used in query, @param key is the queried key, return the leafpiece responsible for the given key
         */

        forceinline result_find search_one(const key_type key,access_stats *query_stats) const {
            // 从 root 开始，即，从root开始的结尾开始
            struct timeval ct1, ct2;
            double ctimeuse;
            if (vstats.innernum > 1){
                auto itlevel = innerlevels.end() -1; // root level 的 ittrater

                auto pos = (*root)(key);
                auto lo = PGM_SUB_EPS(pos, ErrorRecursive + 1) + (*(--itlevel))->innerpieces.begin();  //,(--itlevel)->size()
                auto hi = PGM_ADD_EPS(pos, ErrorRecursive,(*itlevel)->innerpieces.size()-1) + (*itlevel)->innerpieces.begin();
                auto curit = (*(itlevel--))->innerpieces.begin();
                // 在 low 和 high 范围内的cand;idate innerpieces  中 确定所属的 inner piece
                static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece)
                if (key >= hi->startkey){curit = hi;}
                else{
                    if constexpr (4 <= linear_search_threshold ){   //ErrorRecursive
                        for (; std::next(lo)->startkey <= key; ++lo)
                            continue;
                        curit = lo;
                    }
                    else{ curit = std::lower_bound(lo, hi, key);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
                    }
                }
                // 在lo 和 hi 的范围内，定位 key 所属的inner piec, 最后定位到leaf piece
                for (; itlevel >= innerlevels.begin();--itlevel  ){
                    pos = (*curit)(key);
                    lo = PGM_SUB_EPS(pos, ErrorRecursive + 1) + (*itlevel)->innerpieces.begin();
                    hi = PGM_ADD_EPS(pos, ErrorRecursive,(*itlevel)->innerpieces.size()-1) + (*itlevel)->innerpieces.begin();
                    // 在 low 和 high 范围内的candidate innerpieces  中 确定所属的 inner piece
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
                        else{ curit = std::lower_bound(lo, hi, key);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
                        }
                    }
                }

                // leaf level
                auto curleaf = leaflevel.leafpieces.begin();
                pos = (*curit)(key);

                auto leaflo = PGM_SUB_EPS2(pos, Error+1 ,leaflevel.leafpieces.size()) + leaflevel.leafpieces.begin();
                auto leafhi = PGM_ADD_EPS(pos, Error,leaflevel.leafpieces.size()-1) + leaflevel.leafpieces.begin();


                Leafpiece* x = *leaflo;
                Leafpiece* y = *leafhi;
                if (key >= (*leafhi)->startkey){curleaf = leafhi;}
                else{
                    if constexpr (4 <= linear_search_threshold ){   //ErrorRecursive
                        for (; (*leaflo)->endkey < key; leaflo++)
                            continue;
                        curleaf = leaflo;
                    }
                    else{ curleaf = std::lower_bound(leaflo, leafhi, key);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
                    }
                }
                Leafpiece* a = *curleaf;
                pos = (*a)(key);

                int slotlo = PGM_SUB_EPS(pos, ErrorRecursive + 1) ;  //,(--itlevel)->size()
                int slothi = PGM_ADD_EPS(pos, ErrorRecursive,(a->slotkey.size()-1));
                int resslot = slotlo;
                for (; resslot<slothi;resslot++){
                    if (key <= a->slotkey[resslot])
                    {
                        break;
                    }
                }
                gettimeofday(&ct1, NULL);
                if (key != a->slotkey[resslot])
                {

//                    if (sort_list==NULL) {
//                        cout << " i need You, my Lord, please help me!" << endl;
//                    }
 // find in sort_list
                    auto sortpos = (*sort_list)(key);
                    if (sort_list->slotkey[sortpos] == key){
//                            auto sss1 = sort_list->slotkey[sortpos+1];
//                            auto sss2 = sort_list->slotkey[sortpos];
//                            auto sss3 = sort_list->slotkey[sortpos-1];
//                            auto sort_len = sort_list->slotkey.size();
                        result_find index_res = result_find(false,sort_list->slotflag[sortpos],sort_list,sortpos);   // sort_list data
                        return index_res;
                    }
                    else{
                        cout << " i need You, my Lord, please help me! not find the key" << endl;
                    }

                }
                gettimeofday(&ct2, NULL);
                ctimeuse = (ct2.tv_sec - ct1.tv_sec) + (double) (ct2.tv_usec - ct1.tv_usec) / 1000000.0;
                query_stats->sort_piecetime += ctimeuse;
                result_find index_res = result_find(true,a->slotflag[resslot],a,resslot);   // regular data
                return index_res;
            }
            else{

                // 只有一层leaf
                static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece)

                // leaf level
                auto curleaf = leaflevel.leafpieces.begin();
                auto pos = (*root)(key);

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
                    else{ curleaf = std::lower_bound(leaflo, leafhi, key);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
                    }
                }
                Leafpiece* a = *curleaf;
                pos = (*a)(key);

                int slotlo = PGM_SUB_EPS(pos, ErrorRecursive + 1) ;  //,(--itlevel)->size()
                int slothi = PGM_ADD_EPS(pos, ErrorRecursive,(a->slotkey.size()-1));
                int resslot = slotlo;
                for (; resslot<slothi;resslot++){
                    if (key <= a->slotkey[resslot])
                    {
                        break;
                    }
                }
                gettimeofday(&ct1, NULL);
                if (key != a->slotkey[resslot])
                {
//                    if (sort_list==NULL) {
//                        cout << " i need You, my Lord, please help me!" << endl;
//                    }
//                    else {  // find in sort_list
                    auto sortpos = (*sort_list)(key);
                    if (sort_list->slotkey[sortpos] == key){
                        result_find index_res = result_find(false,sort_list->slotflag[sortpos],sort_list,sortpos);   // sort_list data
                        return index_res;
                    }
                    else{
                        cout << " i need You, my Lord, please help me!" << endl;
                    }
//                    }
                }
                gettimeofday(&ct2, NULL);
                ctimeuse = (ct2.tv_sec - ct1.tv_sec) + (double) (ct2.tv_usec - ct1.tv_usec) / 1000000.0;
                query_stats->sort_piecetime += ctimeuse;
                result_find index_res = result_find(true,a->slotflag[resslot],a,resslot);   // regular data
                return index_res;
            }
        }




        forceinline range_res_find search_range(const key_type firstkey, const key_type lastkey,access_stats *query_stats) const {
            struct timeval ct1, ct2;
            double ctimeuse;
            if (vstats.innernum > 1){
                auto itlevel = innerlevels.end() -1; // root level 的 ittrater

                auto pos = (*root)(firstkey);
                auto lo = PGM_SUB_EPS(pos, ErrorRecursive + 1) + (*(--itlevel))->innerpieces.begin();  //,(--itlevel)->size()
                auto hi = PGM_ADD_EPS(pos, ErrorRecursive,(*itlevel)->innerpieces.size()-1) + (*itlevel)->innerpieces.begin();
                auto curit = (*(itlevel--))->innerpieces.begin();
                // 在 low 和 high 范围内的cand;idate innerpieces  中 确定所属的 inner piece
                static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece)
                if (firstkey >= hi->startkey){curit = hi;}
                else{
                    if constexpr (4 <= linear_search_threshold ){   //ErrorRecursive
                        for (; std::next(lo)->startkey <= firstkey; ++lo)
                            continue;
                        curit = lo;
                    }
                    else{ curit = std::lower_bound(lo, hi, firstkey);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
                    }
                }
                // 在lo 和 hi 的范围内，定位 key 所属的inner piec, 最后定位到leaf piece
                for (; itlevel >= innerlevels.begin();--itlevel  ){
                    pos = (*curit)(firstkey);
                    lo = PGM_SUB_EPS(pos, ErrorRecursive + 1) + (*itlevel)->innerpieces.begin();
                    hi = PGM_ADD_EPS(pos, ErrorRecursive,(*itlevel)->innerpieces.size()-1) + (*itlevel)->innerpieces.begin();
                    // 在 low 和 high 范围内的candidate innerpieces  中 确定所属的 inner piece
                    static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece)
                    if (firstkey >= hi->startkey){curit = hi;}
                    else{
                        if (firstkey <= (*lo).startkey){
                            pos = std::min<size_t>((*curit)(firstkey), (*curit).intercept);
                            lo = PGM_SUB_EPS(pos, ErrorRecursive + 1) + (*itlevel)->innerpieces.begin();
                            hi = PGM_ADD_EPS(pos, ErrorRecursive,(*itlevel)->innerpieces.size()-1) + (*itlevel)->innerpieces.begin();
                        }
                        if constexpr (4 <= linear_search_threshold ){   //ErrorRecursive
                            for (; std::next(lo)->startkey <= firstkey; ++lo)
                                continue;
                            curit = lo;
                        }
                        else{ curit = std::lower_bound(lo, hi, firstkey);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
                        }
                    }
                }

                // leaf level
                auto curleaf = leaflevel.leafpieces.begin();
                pos = (*curit)(firstkey);

                auto leaflo = PGM_SUB_EPS2(pos, Error+1 ,leaflevel.leafpieces.size()) + leaflevel.leafpieces.begin();
                auto leafhi = PGM_ADD_EPS(pos, Error,leaflevel.leafpieces.size()-1) + leaflevel.leafpieces.begin();


//                Leafpiece* x = *leaflo;
//                Leafpiece* y = *leafhi;
                if (firstkey >= (*leafhi)->startkey){curleaf = leafhi;}
                else{
                    if constexpr (4 <= linear_search_threshold ){   //ErrorRecursive
                        for (; (*leaflo)->endkey < firstkey; leaflo++)
                            continue;
                        curleaf = leaflo;
                    }
                    else{ curleaf = std::lower_bound(leaflo, leafhi, firstkey);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
                    }
                }
                Leafpiece* a = *curleaf;
                pos = (*a)(firstkey);

                int slotlo = PGM_SUB_EPS(pos, ErrorRecursive + 1) ;  //,(--itlevel)->size()
                int slothi = PGM_ADD_EPS(pos, ErrorRecursive,(a->slotkey.size()-1));
                int resslot = slotlo;
                for (; resslot<slothi;resslot++){
                    if (firstkey <= a->slotkey[resslot])
                    {
                        break;
                    }
                }
                /*
                gettimeofday(&ct1, NULL);
                if (firstkey != a->slotkey[resslot])
                {

                    if (sort_list==NULL) {
                        cout << " i need You, my Lord, please help me!" << endl;
                    }
                    // find in sort_list
                    auto sortpos = (*sort_list)(firstkey);
                    if (sort_list->slotkey[sortpos] == firstkey){
                        result_find index_res = result_find(false,sort_list->slotflag[sortpos],sort_list,sortpos);   // sort_list data
                        return index_res;
                    }
                    else{
                        cout << " i need You, my Lord, please help me! not find the key" << endl;
                    }

                }
                gettimeofday(&ct2, NULL);
                ctimeuse = (ct2.tv_sec - ct1.tv_sec) + (double) (ct2.tv_usec - ct1.tv_usec) / 1000000.0;
                query_stats->sort_piecetime += ctimeuse;
                */

                range_res_find rangeresult = range_res_find(true, firstkey, resslot, lastkey);

                // 在接下来的leaf piece 中 定位lastkey
                auto cur_res = (*curleaf)->is_in_leaf(lastkey, Error);
                while (cur_res.first == false) {
                    rangeresult.findleafs.push_back((*curleaf));
                    ++curleaf;
                    cur_res = (*curleaf)->is_in_leaf(lastkey, Error);
                }
                rangeresult.findleafs.push_back(*curleaf);
                rangeresult.lastslot = cur_res.second;

                return rangeresult;
            }
            else{

                // 只有一层leaf
                static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece)

                // leaf level
                auto curleaf = leaflevel.leafpieces.begin();
                auto pos = (*root)(firstkey);

                auto leaflo = PGM_SUB_EPS2(pos, Error + 1,leaflevel.leafpieces.size()) + leaflevel.leafpieces.begin();
                auto leafhi = PGM_ADD_EPS(pos, Error,leaflevel.leafpieces.size()-1) + leaflevel.leafpieces.begin();
                auto h = PGM_ADD_EPS(pos, Error,leaflevel.leafpieces.size()-1);
//
//                Leafpiece*  y= *leafhi;
//                Leafpiece*  z= *leaflo;
                if (firstkey >= (*leafhi)->startkey){curleaf = leafhi;}
                else{
                    if constexpr (4 <= linear_search_threshold ){   //ErrorRecursive
                        for (; (*leaflo)->endkey < firstkey; leaflo++){
//                            Leafpiece*  x= *leaflo;
                            continue;
                        }

                        curleaf = leaflo;
                    }
                    else{ curleaf = std::lower_bound(leaflo, leafhi, firstkey);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
                    }
                }
                Leafpiece* a = *curleaf;
                pos = (*a)(firstkey);

                int slotlo = PGM_SUB_EPS(pos, ErrorRecursive + 1) ;  //,(--itlevel)->size()
                int slothi = PGM_ADD_EPS(pos, ErrorRecursive,(a->slotkey.size()-1));
                int resslot = slotlo;
                for (; resslot<slothi;resslot++){
                    if (firstkey <= a->slotkey[resslot])
                    {
                        break;
                    }
                }
                /*
                gettimeofday(&ct1, NULL);
                if (firstkey != a->slotkey[resslot])
                {
//                    if (sort_list==NULL) {
//                        cout << " i need You, my Lord, please help me!" << endl;
//                    }
//                    else {  // find in sort_list
                    auto sortpos = (*sort_list)(firstkey);
                    if (sort_list->slotkey[sortpos] == firstkey){
                        result_find index_res = result_find(false,sort_list->slotflag[sortpos],sort_list,sortpos);   // sort_list data
                        return index_res;
                    }
                    else{
                        cout << " i need You, my Lord, please help me!" << endl;
                    }
                }
                gettimeofday(&ct2, NULL);
                ctimeuse = (ct2.tv_sec - ct1.tv_sec) + (double) (ct2.tv_usec - ct1.tv_usec) / 1000000.0;
                query_stats->sort_piecetime += ctimeuse;
                */

                range_res_find rangeresult = range_res_find(true, firstkey, resslot, lastkey);

                // 在接下来的leaf piece 中 定位lastkey
                auto cur_res = (*curleaf)->is_in_leaf(lastkey, Error);
                while (cur_res.first == false) {
                    rangeresult.findleafs.push_back((*curleaf));
                    ++curleaf;
                    cur_res = (*curleaf)->is_in_leaf(lastkey, Error);
                }
                rangeresult.findleafs.push_back(*curleaf);
                rangeresult.lastslot = cur_res.second;

                return rangeresult;
            }




        }


        template<class lru_type>
        forceinline bool update_one(const key_type key, key_type* payload,lru_type interchain)  {
            // 判断一下，如果 小于 last leaf 的end key, 则插入到sortlist 中
            if (key < this->m_tailleaf->endkey){
                // 插入到 sort_list 中
                this->sort_list->insert(key,payload);
//                cout << "thank You, my Lord! " << endl ;
            }
            else{
                this->append_one(key, payload,Error, interchain);
            }
        }

/*
        bool update_search_one(const key_type key, const vector<key_type> payload) const {
            // 从 root 开始，即，从root开始的结尾开始
            if (vstats.innernum > 1){
                auto itlevel = innerlevels.end() -1; // root level 的 ittrater
                auto pos = (*root)(key);
                auto lo = PGM_SUB_EPS(pos, ErrorRecursive + 1) + (*(--itlevel))->innerpieces.begin();  //,(--itlevel)->size()
                auto hi = PGM_ADD_EPS(pos, ErrorRecursive,(*itlevel)->innerpieces.size()-1) + (*itlevel)->innerpieces.begin();
                auto curit = (*(itlevel--))->innerpieces.begin();
                // 在 low 和 high 范围内的cand;idate innerpieces  中 确定所属的 inner piece
                static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece)
                if (key >= hi->startkey){curit = hi;}
                else{
                    if constexpr (4 <= linear_search_threshold ){   //ErrorRecursive
                        for (; std::next(lo)->startkey <= key; ++lo)
                            continue;
                        curit = lo;
                    }
                    else{ curit = std::lower_bound(lo, hi, key);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
                    }
                }
                // 在lo 和 hi 的范围内，定位 key 所属的inner piec, 最后定位到leaf piece
                for (; itlevel >= innerlevels.begin();--itlevel  ){
                    pos = (*curit)(key);
                    lo = PGM_SUB_EPS(pos, ErrorRecursive + 1) + (*itlevel)->innerpieces.begin();
                    hi = PGM_ADD_EPS(pos, ErrorRecursive,(*itlevel)->innerpieces.size()-1) + (*itlevel)->innerpieces.begin();
                    // 在 low 和 high 范围内的candidate innerpieces  中 确定所属的 inner piece
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
                        else{ curit = std::lower_bound(lo, hi, key);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
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
                    else{ curleaf = std::lower_bound(leaflo, leafhi, key);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
                    }
                }
                Leafpiece* a = *curleaf;
                pos = (*a)(key);

                int slotlo = PGM_SUB_EPS(pos, ErrorRecursive + 1) ;  //,(--itlevel)->size()
                int slothi = PGM_ADD_EPS(pos, ErrorRecursive,(a->slotkey.size()-1));
                int resslot = slotlo;
                for (; resslot<slothi;resslot++){
                    if (key >= a->slotkey[resslot])
                    {
                        break;
                    }
                }
                if (key != a->slotkey[resslot])
                {
                    // binary serch find in sort_list;
                    if (a->sort_list.empty()){
                        // directly insert in sort_list
                        a->sort_list.emplace_back(pair<key_type,vector<key_type>>(key,payload));
//                        cout<< " i need You, my Lord, please help me!"<< endl;
                        return true;
                    }
                    else{
                        int begin = 0, end = a->sort_list.size()-1,mid;  // find last Equal or small, 找到最后一个大于或等于 key 的
                        while(begin <= end){
                            mid = (end + begin) / 2;
                            if(a->sort_list[mid].first >=key) {
                                end = mid-1;
                            } else
                                begin = mid +1;
                        }
                        if(a->sort_list[end+1].first!=key){
                            a->sort_list.insert(a->sort_list.begin()+begin,pair<key_type,vector<key_type> >(key,payload));
                            return true;
                        }
                        return false;
                    }
                }
                return false;
            }
            else{

                // 只有一层leaf
                static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Innerpiece);   //sizeof(Innerpiece)
//                if (key == 80728301)
//                    cout<<key<< "  my Lord, i need You!"<<endl;
                // leaf level
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
                    else{ curleaf = std::lower_bound(leaflo, leafhi, key);   // upper_bound 返回在前闭后开区间查找key 的上届，返回大于 key 的第一个元素位置； lower_bound 返回第一个大于或等于 key 的元素的位置
                    }
                }
                Leafpiece* a = *curleaf;
                pos = (*a)(key);

                int slotlo = PGM_SUB_EPS(pos, ErrorRecursive + 1) ;  //,(--itlevel)->size()
                int slothi = PGM_ADD_EPS(pos, ErrorRecursive,(a->slotkey.size()-1));
                int resslot = slotlo;
                for (; resslot<slothi;resslot++){
                    if (key >= a->slotkey[resslot])
                    {
                        break;
                    }
                }
                if (key != a->slotkey[resslot])
                {
                    // binary serch find in sort_list;
                    if (a->sort_list.empty()){
                        // directly insert in sort_list
                        a->sort_list.emplace_back(pair<key_type,vector<key_type>>(key,payload));
//                        cout<< " i need You, my Lord, please help me!"<< endl;
                        return true;
                    }
                    else{
                        int begin = 0, end = a->sort_list.size()-1,mid;  // find last Equal or small, 找到最后一个大于或等于 key 的
                        while(begin <= end){
                            mid = (end + begin) / 2;
                            if(a->sort_list[mid].first >=key) {
                                end = mid-1;
                            } else
                                begin = mid +1;
                        }
                        if(a->sort_list[end+1].first!=key){
                            a->sort_list.insert(a->sort_list.begin()+begin,pair<key_type,vector<key_type> >(key,payload));
                            return true;
                        }
                        return false;
                    }
                }
                return false;

            }

        }
*/


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
            double innodesize = 24;  // startkey, double(slope,intercept)
            double leafnodesize = 32 ; // startkey, endkey, double(slope,intercept),  slotkey(算入了data usage 中，因为在introchain 中 只存储了value)
            double hashlrusize = 48;
            double no_hashsize = 20;  //prev,next, slot-short int
            double dataV = (valuesize+1)*8;
            double datausage = double(inkeynum*(dataV))/1048576;
            double addusage = double((exkeynum+inkeynum)*(1+8)+exkeynum*(sizeof(pageid_type)+sizeof(pageOff_type))+exkeynum*sizeof(key_type) )/1048576;  //
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
            double leafnodesize = 2 * sizeof(key_type) + 8 + 8 ; // startkey, endkey, double(slope,intercept),  slotkey(算入了data usage 中，因为在introchain 中 只存储了value)
            double hashlrusize = sizeof(key_type) + 8+ 16 + 16;
            double no_hashsize = 16+4;  //prev,next, slot-short int
            double dataV = (valuesize+1)*8;
            double datausage = double(inkeynum*(dataV))/1024/1024;
            double addusage = double((exkeynum+inkeynum)*(1+8)+exkeynum*(sizeof(pageid_type)+sizeof(pageOff_type))+exkeynum*sizeof(key_type) )/double(1048576);  // exkey, flag, pageID,offset   (pageid-int-8 offset short int), flag for data in memory
            // addusage: (exkeynum+inkeynum)*(1+8) —— bitmap + pointer;  exkeynum*8 —— pageid(ungined int),offset(unsigned int); exkeynum*sizeof(key_type) ——exkey;
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
            for (unsigned int i = 1; i<leaflevel.leafpieces.size()-1;i++){
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
     * 开始定义 FILM 中用的的struct
     */

    template<typename key_type,typename value_type>
    struct FILMinsert<key_type,value_type>::Leafpiece{

        key_type startkey;
        key_type endkey;
        vector<key_type> slotkey;

        vector< bool > slotflag;   //adalru::Node<key_type, std::vector<key_type> >*
        vector< void*  > slotdata;   //adalru::Node<key_type, std::vector<key_type> >*
        adalru::localLRU <lruOff_type,key_type* >  intrachain;    //unsigned short int
//        vector< pair<key_type,vector<key_type>> > sort_list;    // 用于 存储 任意位置插入的数据
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
        forceinline size_t operator()(const key_type &key) const {
            auto pos = int64_t(slope * (key - startkey)) + intercept;
            return pos > 0 ? size_t(pos) : 0ull;
        }

        forceinline void update(typename internal::insertPWLF<key_type, int>::CanonicalSegment &cs){
            startkey = cs.first;
            endkey = cs.last;
            auto[cs_slope, cs_intercept] = cs.get_floating_point_segment(startkey);
            if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
                throw std::overflow_error("Change the type of Segment::intercept to int64");
            slope = cs_slope;
            intercept = cs_intercept;
        }

        // 判断key 是否在当前leaf 中
        inline pair<bool,unsigned long> is_in_leaf(const key_type &key,const int error){
            pair<bool,unsigned long> is_in(false,0);
            if (key<=this->endkey ){
                auto pos = (*this)(key);

                int slotlo = PGM_SUB_EPS(pos, error + 1) ;  //,(--itlevel)->size()
                int slothi = PGM_ADD_EPS(pos, error,(this->slotkey.size()-1));
                for (; slotlo<=slothi;slotlo++){
//                    if (slotlo < 0)
//                    {
//                        cout<< "my Lord, please help, i trust in You!" << endl;
//                    }
                    if (key < this->slotkey[slotlo] || key >=this->startkey)
                    {
                        is_in.first = true;
                        is_in.second = slotlo;
                        break;
                    }
                    else
                    {
                        if (key == this->endkey && key<this->startkey){
                            is_in.first = true;
                            is_in.second = slotlo;
                            break;
                        }
//                        auto aaa = this->slotkey[slotlo];
//                        cout<< "my Lord, please help, i trust in You!" << endl;
                    }
//                    else{
//                        // check in sort piece
//                        cout << "i need You, my Lord" << endl;
////                        auto sortpos = (*sort_list)(key);
////                        if (sort_list->slotkey[sortpos] == key){
////                            is_in.first = true;
////                            is_in.second = slotlo;
////                            break;
////                        }
//                    }
                }

                return is_in;
            }
            else
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
        forceinline size_t operator()(const key_type &key) const {
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
//        vector<internal::insertPWLF<key_type, int>*> opt;
        internal::insertPWLF<key_type, int>* opt;
        unsigned int pos = 0;
        unsigned int nextpos = 0;
    };

    template<typename key_type,typename value_type>
    struct FILMinsert<key_type,value_type>::Leaflevel{
        std::vector<FILMinsert<key_type,value_type>::Leafpiece* > leafpieces;
//        vector<internal::insertPWLF<key_type, int>*> opt;
        internal::insertPWLF<key_type, int>* opt;
        unsigned int pos = 0;
    };

#pragma pack(pop)

    typedef filminsert::FILMinsert< key_type, key_type* > filmadatype;

    typedef pair<filmadatype::Leafpiece*,unsigned short int>  filmadalrupair;  //lru 中value 的type ，是一个pair，first is 所属的leaf，second 是slot in leaf
    typedef  adalru::hashLRU <key_type,filmadatype::Leafpiece* , adalru::Node<key_type ,filmadatype::Leafpiece*>* > filmadalrutype;
    typedef adalru::localLRU <unsigned short,key_type* > locallrutype;
    typedef filmstorage::filmdisk<key_type> filmadadisk;
    typedef filmstorage::filmmemory<key_type,key_type*,filmadatype, filmadalrutype,filmadadisk> filmadamemory;
    typedef filmadamemory::memoryusage memoryusage;
    //    adalru::hashLRU <key_type,Leafpiece* , adalru::Node<key_type ,Leafpiece*>* >

    void runtimeevictkeytopage(filmadamemory *filmmem,filmadatype *filmindex,filmadadisk *diskpage , access_stats *r_stats){
        struct timeval lt1, lt2,dt1,dt2,rdt1,rdt2,wdt1,wdt2;
        double ltimeuse,rdtimeuse,wdtimeuse;  // 读磁盘的时间，写磁盘的时间
        double dtimeuse;
        if (filmmem->evictPoss.size() != 0){
            for (int i = 0; i< filmmem->evictPoss.size(); i++){
                gettimeofday(&lt1, NULL);
                auto writeevict =  filmmem->evictPoss[i];

                // evict data, get the tail node, remove the tail of intrachain
                auto evictleaf = filmmem->lru->get_tail();
                auto evictslotV = evictleaf->intrachain.poptail();   // the tail of the accessed leaf
                gettimeofday(&lt2, NULL);
                ltimeuse = (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                r_stats->lrutime += ltimeuse;

                filmmem->evictkeytoinpage(evictleaf->slotkey[evictslotV->key], evictslotV->value,
                                          diskpage, writeevict);
                evictleaf->slotflag[evictslotV->key] = false;
                evictleaf->slotdata[evictslotV->key] = writeevict;
                delete evictslotV;
                evictslotV = NULL;

            }
            filmindex->inkeynum -= filmmem->evictPoss.size();
            filmindex->exkeynum += filmmem->evictPoss.size();
            filmmem->evictPoss.clear();
            // 批量将inmempages写出磁盘
            gettimeofday(&wdt1, NULL);
            int pagenum = filmmem->runtimeevictpagestodisk(diskpage);
            r_stats->writenum += pagenum;
            gettimeofday(&wdt2, NULL);
            wdtimeuse = (wdt2.tv_sec - wdt1.tv_sec) + (double) (wdt2.tv_usec - wdt1.tv_usec) / 1000000.0;
            r_stats->wdisktime += wdtimeuse;
        }

    }

    double runtimeevictkeytopage2(filmadamemory *filmmem,double totalusemem, filmadatype *filmindex,filmadadisk *diskpage , access_stats *r_stats){
        struct timeval wdt1,wdt2,lt1,lt2;
        double wdtimeuse =0.0;  // 写磁盘的时间
        double ltimeuse = 0.0;

        if (filmmem->evictPoss.size() != 0){
            gettimeofday(&wdt1, NULL);
            for (unsigned int i = 0; i< filmmem->evictPoss.size(); i++){
                auto writeevict =  filmmem->evictPoss[i];
                // evict data, get the tail node, remove the tail of intrachain
                gettimeofday(&lt1, NULL);
                auto evictleaf = filmmem->lru->get_tail();
                auto evictslotV = evictleaf->intrachain.poptail();   // the tail of the accessed leaf
                gettimeofday(&lt2, NULL);
                ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                filmmem->evictkeytoinpage(evictleaf->slotkey[evictslotV->key], evictslotV->value,
                                          diskpage, writeevict);
                evictleaf->slotflag[evictslotV->key] = false;
                evictleaf->slotdata[evictslotV->key] = writeevict;

                delete evictslotV;
                evictslotV = NULL;

            }
            filmindex->inkeynum -= filmmem->evictPoss.size();
            filmindex->exkeynum += filmmem->evictPoss.size();
            filmmem->evictPoss.clear();
            // 批量将inmempages写出磁盘
            int pagenum = filmmem->runtimeevictpagestodisk(diskpage);
            r_stats->writenum += pagenum;
            gettimeofday(&wdt2, NULL);
            wdtimeuse += (wdt2.tv_sec - wdt1.tv_sec) + (double) (wdt2.tv_usec - wdt1.tv_usec) / 1000000.0;
            wdtimeuse -= ltimeuse;
        }
        else{  // just evict from the original table
            double ratio = (totalusemem-filmmem->threshold)/filmmem->threshold  ;
            unsigned int batch_evict = ceil(filmindex->inkeynum*ratio*0.6) + 300 ;
            if (ratio>0.2){
                batch_evict = ceil(filmindex->inkeynum*ratio*0.5);
            }
            gettimeofday(&wdt1, NULL);
            for (unsigned int i = 0; i< batch_evict; i++)   // the number of records to be evicted in batches
            {
                // evict data, get the tail node, remove the tail of intrachain
                gettimeofday(&lt1, NULL);
                auto evictleaf = filmmem->lru->get_tail();
                auto evictslotV = evictleaf->intrachain.poptail();   // the tail of the accessed leaf
                gettimeofday(&lt2, NULL);
                ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                auto writeevict = filmmem->evictkeytoinpage(evictleaf->slotkey[evictslotV->key], evictslotV->value,
                                          diskpage);
                evictleaf->slotflag[evictslotV->key] = false;
                evictleaf->slotdata[evictslotV->key] = writeevict;
                delete evictslotV;
                evictslotV = NULL;
            }
            filmindex->inkeynum -= batch_evict;
            filmindex->exkeynum += batch_evict;

            int pagenum = filmmem->runtimeevictpagestodisk(diskpage);
            r_stats->writenum += pagenum;
            gettimeofday(&wdt2, NULL);
            wdtimeuse += (wdt2.tv_sec - wdt1.tv_sec) + (double) (wdt2.tv_usec - wdt1.tv_usec) / 1000000.0;
            wdtimeuse -= ltimeuse;
        }
        r_stats->wdisktime += wdtimeuse;
        r_stats->lrutime += ltimeuse;
        return wdtimeuse;
    }



    // 采用unordered_map 实现 dict 的功能
    template<class v_type>   // int bplustree:  pair< oribttype::leaf_node*,pair<short int,int>
    class prepass_dict{
    public:
        int pagenum = 0;
//        std::unordered_map<int,vector<v_type>> prepass;
        std::map< pageid_type,vector<v_type>> prepass;
        struct dict_find {
            /// find result flags
            bool flags;
            vector<v_type> *findinfo;

            inline dict_find()
                    : flags() {}
        };

        void insert_dict(pageid_type pageid,v_type keyinfo){  //keyinfo(first: offset in page, evicted posi pair<int,short>*)
//            typename std::unordered_map<int,vector<v_type>>::const_iterator got = prepass.find(pageid);
            typename std::map<pageid_type,vector<v_type>>::const_iterator got = prepass.find(pageid);
            if (got == prepass.end())
            {
                prepass[pageid].push_back(keyinfo);
                pagenum +=1;
            }
            else
                prepass[pageid].push_back(keyinfo);
        }

        void insert_dict(pageid_type pageid,pageOff_type off, pair<pageid_type ,pageOff_type> * evictposi){  //keyinfo(first: offset in page, evicted posi pair<int,short>*)
//            typename std::unordered_map<int,vector<v_type>>::const_iterator got = prepass.find(pageid);
            typename std::map<pageid_type,vector<v_type>>::const_iterator got = prepass.find(pageid);
            if (got == prepass.end())
            {
                prepass[pageid].push_back(off);
                pagenum +=1;
            }
            else
                prepass[pageid].push_back(off);
        }


        dict_find find_in_dict(pageid_type pageid){
            dict_find res;
//            typename std::unordered_map<int,vector<v_type>>::const_iterator got = prepass.find(pageid);
            typename std::map<pageid_type,vector<v_type>>::const_iterator got = prepass.find(pageid);
            if (got == prepass.end()){
                std::cout<<"not find the key in prepass, my lovely Lord, waiting for you! "<<endl;
                res.flags = false;
                return res;
            }
            else
            {
                res.flags = true;
                res.findinfo = got;
                return res;
            }
        }
    };

    typedef prepass_dict< pair< filmadatype::Leafpiece* ,pair< lruOff_type, pair<pageid_type, pageOff_type>* > > > pre_dict ;

// in the range query, judge whether the keys in each overlapped leaf in memory or in disk
    void range_is_in_or_out(int start, int end, filmadatype::Leafpiece* curleaf ,vector<pair<key_type,key_type*>> *totalres , pre_dict* dict,int sizevalue ,filmadamemory *filmmem,access_stats *r_stats){
        struct timeval lt1, lt2;
        double ltimeuse;
        for (int k = start;k < end;k++)
        {
            bool flag = curleaf->slotflag[k];
            if (flag){
                pair<key_type,key_type*> resitem;
                resitem.first = curleaf->slotkey[k];
                auto finddata = (adalru::Node<lruOff_type, key_type* >*) curleaf->slotdata[k];
                resitem.second = finddata->value;
                totalres->emplace_back(resitem);
                gettimeofday(&lt1, NULL);

                curleaf->intrachain.moveTohead(finddata);// update intrachain
                filmmem->lru->put(curleaf->startkey,curleaf);                 // 更新global lru
                gettimeofday(&lt2, NULL);
                ltimeuse = (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                r_stats->lrutime += ltimeuse;
            }
            else{  // prepass into dict, that access disk
                auto evictpos = (pair<pageid_type ,pageOff_type>*)curleaf->slotdata[k];
                pair<pageid_type ,pageOff_type>* diskpos = evictpos;
                pageid_type pageid = diskpos->first;
                auto off = diskpos->second;
                //dict 中插入的是 pageid，所属leafpiece，（在leaf 中的slot 位置，保存evictposi的地址（ pointer ）——用于后期存储被移出key 的外存地址）
                dict->insert_dict(pageid, pair< filmadatype::Leafpiece*,pair<lruOff_type,pair<pageid_type ,pageOff_type>* > >(curleaf,pair<lruOff_type,pair<pageid_type ,pageOff_type>*>(k,evictpos)));

//                cout<<"my Lord, my heart is designed for You!"<<endl;
            }
        }
    }

    int range_prepass(filmadamemory *filmmem,filmadatype::range_res_find *index_result,vector<pair<key_type,key_type*>> *result,pre_dict* dict,filmadatype *filmindex,access_stats *r_stats){  // 将range query 由btree 得到的信息，进行怕热pass，遍历find leaf中的每一个leaf，如果slotdata.first = true, 则直接放入reslut 中，否则将其放入prepassdict 中

        int leafnum = index_result->findleafs.size();
        if (leafnum == 1){  //all the keys in range belonging to a certain leaf
            range_is_in_or_out(index_result->firstslot,index_result->lastslot+1,index_result->findleafs[0],result,dict,filmindex->valuesize,filmmem,r_stats);
        }
        else {  // the first leaf and last leaf need to be deal with seperately
            // deal with firstleaf
            //index_result->findleafs[0];
            range_is_in_or_out(index_result->firstslot,index_result->findleafs[0]->slotkey.size(),index_result->findleafs[0],result,dict,filmindex->valuesize,filmmem,r_stats);
            // deal with the middle leafs
            for (int i=1; i < leafnum-1; i++){
                range_is_in_or_out(0,index_result->findleafs[i]->slotkey.size(),index_result->findleafs[i],result,dict,filmindex->valuesize,filmmem,r_stats);
            }
            //deal with lastleaf
            //index_result->findleafs[leafnum-1];
            range_is_in_or_out(0,index_result->lastslot+1,index_result->findleafs[leafnum-1],result,dict,filmindex->valuesize,filmmem,r_stats);
        }
        return 0;
    }


// 一次至少读取1024k 的数据
    int range_read_from_disk(vector<pair<key_type,key_type*>> *result,pre_dict* dict,filmadamemory *filmmem,filmadatype *filmindex,filmadadisk *pagedisk, access_stats *r_stats) {
        // 打开磁盘文件
        struct timeval lt1, lt2, dt1, dt2, rdt1, rdt2;
        double ltimeuse = 0.0;  // 读磁盘的时间，写磁盘的时间
        double dtimeuse, rdtimeuse;
        int cross = dict->prepass.rbegin()->first - dict->prepass.begin()->first + 1;
//        cout<< cross << " ";
        r_stats->pagecross += dict->pagenum;
        int fd;
        key_type *buf;

        pair<key_type, key_type *> res;

        unsigned long int fix_buf_size = pagedisk->pagesize * 8;  // 每个磁盘页的大小


        fd = open(pagedisk->pagefile, O_RDWR | O_DIRECT, 0755);
        if (cross * pagedisk->pagesize <= 131072)// 262144 2048k 判断一下，如果最大page 和最小page，相差超过1024k 131072，则采用seek的形式进行读取，
        {
            unsigned long int buf_size = fix_buf_size * cross;  // 一次性读入
            int ret = posix_memalign((void **) &buf, 512, buf_size);
            memset(buf, 'c', buf_size);

            gettimeofday(&rdt1, NULL);
            pageid_type firstp = dict->prepass.begin()->first;
            unsigned long aoffset = firstp * fix_buf_size;
            ret = pread(fd, buf, buf_size, static_cast<unsigned long> (aoffset)); // 从 prepass 的第一个页开始，读取所有跨着的页
            gettimeofday(&rdt2, NULL);
            rdtimeuse = (rdt2.tv_sec - rdt1.tv_sec) + (double) (rdt2.tv_usec - rdt1.tv_usec) / 1000000.0;
            r_stats->rdisktime += rdtimeuse;
            gettimeofday(&dt1, NULL);
            for (auto &v: dict->prepass) {
                // buf 中为所有要读取的数据
                int diff = v.first - dict->prepass.begin()->first;
                unsigned long int offset = diff * (pagedisk->recordsize * pagedisk->recordnum);// 把当前页的数据 单独出来
                for (int k = 0; k <
                                v.second.size(); k++) {   //v.second.size()  是要从当前页中读取多少个data， v.second[i] 是每个data 的相关信息是一个pair，v.second[i].first 是所属的leaf， v.second[i].second。first是slot in leaf，second 是 evictposi指向page id 和off
//                    auto slot1 = v.second[k].second.first;
                    auto slot = v.second[k].second.first;
                    pageOff_type off = v.second[k].second.second->second;   // offset in page
                    auto writeevict = v.second[k].second.second;
                    filmadatype::Leafpiece *curpiece = v.second[k].first;
                    if (writeevict->first == filmmem->inpage->pageid) {
                        //vector<key_type> inmemdata = filmmem->inpage->inmemdata;
                        key_type *inmemdata = filmmem->inpage->inmemdata;
                        res.first = inmemdata[writeevict->second];
                        res.second = new key_type[filmindex->valuesize];
                        for (int i = 0; i < filmindex->valuesize; i++) {
                            res.second[i] = inmemdata[writeevict->second + 1 + i];
                        }
                    } else {
                        res.first = buf[offset + off];
                        res.second = new key_type[filmindex->valuesize];
                        for (int i = 0; i < filmindex->valuesize; i++) {
                            res.second[i] = buf[offset + off + 1 + i];
                            // 这里需要算出来 当前的页 距离第一个页的相对位置
                        }
                    }

                    filmmem->evictPoss.emplace_back(writeevict);
                    filmindex->inkeynum++;
                    filmindex->exkeynum--;
                    result->push_back(res);

                    if (res.first != curpiece->slotkey[slot]) {
                        cout << "i need You, my Lord, Thank You for all the thing! " << endl;
                    }


                    // 将res 插入到curpiece 的intrachain，// 将被访问的record 的数据放入内存中，即其所属的leaf，并修改其slotdata 的flag 为 true
                    gettimeofday(&lt1, NULL);
                    curpiece->slotdata[slot] = curpiece->intrachain.put(slot, res.second);
                    curpiece->slotflag[slot] = true;

                    // 将curpiece 插入到interchain 中，即mem->lru中
                    filmmem->lru->put(curpiece->startkey, curpiece);// 更新 interchain
                    gettimeofday(&lt2, NULL);
                    ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
                }
            }
            gettimeofday(&dt2, NULL);
            dtimeuse = (dt2.tv_sec - dt1.tv_sec) + (double) (dt2.tv_usec - dt1.tv_usec) / 1000000.0;
            r_stats->disktime += dtimeuse;
            r_stats->disktime -= ltimeuse;
            r_stats->lrutime += ltimeuse;
        } else {
            int ret = posix_memalign((void **) &buf, 512, fix_buf_size);
            memset(buf, 'c', fix_buf_size);
            // for loop to read the data of each page needed
            pageid_type abslotepageid = 0;
            gettimeofday(&dt1, NULL);
            for (auto &v: dict->prepass) {
                pageid_type pageid = v.first;  // v.first 是 pageid， v.second 是 这个page 要读取的数据的信息，
                if (pageid == filmmem->inpage->pageid) {
                    cout << " Jesus, Son of David, please have pity on me, read from in-memory page" << endl;
                    key_type *inmemdata = filmmem->inpage->inmemdata;

                    for (int vi = 0; vi <
                                     v.second.size(); vi++) {   //v.second.size()  是要从当前页中读取多少个data， v.second[i] 是每个data 的相关信息是一个pair，v.second[i].first 是所属的leaf， v.second[i].second。first是slot in leaf，second 是 evictposi指向page id 和off
                        auto slot = v.second[vi].second.first;   // v.first 是 page id，v.second 是要从该页中读进来的每个record 的信息，是一个pair，
                        // 该pair 的first 是所属的leaf，该pair 的second 是一个pair，
                        // 包括：first 是在该leaf 中的slot，second 是page 信息，包括page id 和在该page 中的offset
                        short int off = v.second[vi].second.second->second;  //offset in page
                        auto writeevict = v.second[vi].second.second;
                        filmadatype::Leafpiece *curpiece = v.second[vi].first;
                        res.first = inmemdata[off];
                        res.second = new key_type[filmindex->valuesize];
                        for (int ki = 0; ki < filmindex->valuesize; ki++) {
                            res.second[ki] = inmemdata[off + 1 + ki];
                        }
                        filmmem->evictPoss.emplace_back(writeevict);
                        filmindex->inkeynum++;
                        filmindex->exkeynum--;
                        result->push_back(res);

                        if (res.first != curpiece->slotkey[slot]) {
                            cout << "i need You, my Lord, Thank You for all the thing! " << endl;
                        }
                        // 将res 插入到curpiece 的intrachain，// 将被访问的record 的数据放入内存中，即其所属的leaf，并修改其slotdata 的flag 为 true
                        gettimeofday(&lt1, NULL);
                        curpiece->slotdata[slot] = curpiece->intrachain.put(slot, res.second);
                        curpiece->slotflag[slot] = true;
                        // 将curpiece 插入到interchain 中，即mem->lru中
                        filmmem->lru->put(curpiece->startkey, curpiece);// 更新 interchain
                        gettimeofday(&lt2, NULL);
                        ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
//                      cout<< " i need You, my Lord!!! we all need You! "<< endl;
                    }
                }
                else {

                    unsigned long seekdis = fix_buf_size * (pageid - abslotepageid);
                    lseek(fd, seekdis, SEEK_CUR);
                    ret = read(fd, buf, fix_buf_size);
                    abslotepageid = pageid + 1;

                    //buf 中为所有属于这个页的数据，根据v.second 包含的offset in page 的信息，将对应的值读出来，并放入leaf 中
                    for (int vi = 0; vi <
                                     v.second.size(); vi++) {   //v.second.size()  是要从当前页中读取多少个data， v.second[i] 是每个data 的相关信息是一个pair，v.second[i].first 是所属的leaf， v.second[i].second。first是slot in leaf，second 是 evictposi指向page id 和off
                        auto slot = v.second[vi].second.first;
                        short int off = v.second[vi].second.second->second;  //offset in page
                        auto writeevict = v.second[vi].second.second;
                        filmadatype::Leafpiece *curpiece = v.second[vi].first;
                        res.first = buf[off];
                        res.second = new key_type[filmindex->valuesize];
                        for (int ki = 0; ki < filmindex->valuesize; ki++) {
                            res.second[ki] = buf[off + 1 + ki];
                        }
                        filmmem->evictPoss.emplace_back(writeevict);
                        filmindex->inkeynum++;
                        filmindex->exkeynum--;
                        result->push_back(res);

                        if (res.first != curpiece->slotkey[slot]) {
                            cout << "i need You, my Lord, Thank You for all the thing! " << endl;
                        }

                        // 将res 插入到curpiece 的intrachain，// 将被访问的record 的数据放入内存中，即其所属的leaf，并修改其slotdata 的flag 为 true
                        gettimeofday(&lt1, NULL);
                        curpiece->slotdata[slot] = curpiece->intrachain.put(slot, res.second);
                        curpiece->slotflag[slot] = true;
                        // 将curpiece 插入到interchain 中，即mem->lru中
                        filmmem->lru->put(curpiece->startkey, curpiece);// 更新 interchain
                        gettimeofday(&lt2, NULL);
                        ltimeuse += (lt2.tv_sec - lt1.tv_sec) + (double) (lt2.tv_usec - lt1.tv_usec) / 1000000.0;
//                cout<< " i need You, my Lord!!! we all need You! "<< endl;
                    }
                }
            }
            gettimeofday(&dt2, NULL);
            dtimeuse = (dt2.tv_sec - dt1.tv_sec) + (double) (dt2.tv_usec - dt1.tv_usec) / 1000000.0;
            r_stats->disktime += dtimeuse;
            r_stats->disktime -= ltimeuse;
            r_stats->lrutime += ltimeuse;
        }
        //        cout<<"Jesus, You are my refuge! "<<endl;
        delete[] buf;
        close(fd);

        return 0;
    }

    // excute point query
    pair<access_stats , vector<vector<key_type>> > filmadapointquery(filmadatype *filmindex,vector<key_type> pqueries, int queryn,filmadadisk *pagedisk , filmadamemory *filmmem){
        struct timeval qt1, qt2;
        double timeuse;
        vector<vector<key_type>> totalres;
        access_stats point_stats;
        pair<bool,double >  transflag = filmmem->runtimejudgetrans(&point_stats);
        int periodV = filmmem->reserveMem*1024*1024/((filmindex->valuesize+1)*sizeof(key_type)*10);

        gettimeofday(&qt1, NULL);
        for (unsigned int i = 0; i < pqueries.size(); i++) {
            if (point_stats.disknum >0 && i % periodV == 0  ){
                transflag = filmmem->runtimejudgetrans(&point_stats);
                while (transflag.first)
                {
                    runtimeevictkeytopage(filmmem,filmindex,pagedisk,&point_stats);
                    transflag = filmmem->runtimejudgetrans(&point_stats);
                }
                //                cout<< "my Lord, i need You!" <<endl;
            }

            pair<key_type,key_type*> res;
            key_type querykey = pqueries[i];

            auto index_res = filmindex->search_one(querykey,&point_stats);

            if ( index_res.find == false){
                point_stats.memnum += 1;
                res.first = querykey;
//                res.insert(res.begin()+1,index_res.findleaf->sort_list[index_res.slot].second.begin(),index_res.findleaf->sort_list[index_res.slot].second.end());
            }
            else{
                if (index_res.flags){// 从内存中读数据
                    point_stats.memnum += 1;
                    auto finddata = (adalru::Node<lruOff_type, key_type* >*) index_res.findleaf->slotdata[index_res.slot];
                    res.first = querykey;
                    res.second = finddata->value;
                    //                totalres.emplace_back(res);
                    // 更新 interchain
                    filmmem->lru->put(index_res.findleaf->startkey,index_res.findleaf);

                }
                else{ //从磁盘读数据
                    point_stats.disknum += 1;
                    point_stats.diskpagenum += 1;
                    auto writeevict = (pair<pageid_type ,pageOff_type>*)index_res.findleaf->slotdata[index_res.slot];  //writeevict 指向的是 pageid and offset
                    if (writeevict->first == filmmem->inpage->pageid){
                        res.second = new key_type[filmindex->valuesize];
                        key_type* inmemdata = filmmem->inpage->inmemdata;
                        res.first = inmemdata[writeevict->second];
                        for (int i = 0; i < filmindex->valuesize; i ++){
                            res.second[i]= inmemdata[writeevict->second+1+i];
                        }
                    }
                    else{
                        res = pagedisk->odirectreadfromdisk(writeevict);    // if readfromdisk indicating doesn't use o_direct;
                    }

                    index_res.findleaf->slotdata[index_res.slot] = index_res.findleaf->intrachain.put(index_res.slot,res.second);
                    index_res.findleaf->slotflag[index_res.slot] = true;

                    filmmem->lru->put(index_res.findleaf->startkey,index_res.findleaf);// 更新 interchain
                    filmmem->evictPoss.emplace_back(writeevict);
                    filmindex->inkeynum++;
                    filmindex->exkeynum--;
                }
            }

        }

        gettimeofday(&qt2, NULL);
        cout<<"my Lord, You are my refuge~~~~forever! "<<endl;
        timeuse = (qt2.tv_sec - qt1.tv_sec) + (double) (qt2.tv_usec - qt1.tv_usec) / 1000000.0;
        point_stats.querytimeuse += timeuse;
        point_stats.querytimeuse -= point_stats.computetimeuse;
        return pair<access_stats, vector<vector<key_type>>>(point_stats,totalres);
    }

    // interleave point and range queries
    access_stats filmadaquery(filmadatype *filmindex, vector<key_type> pointqueries,vector<vector<key_type>> rangequeries, unsigned int queryn,filmadadisk *pagedisk , filmadamemory *filmmem){
        struct timeval tqt1, tqt2,xt1,xt2,pt1, pt2,plt1,plt2,rdt1,rdt2;
        double pltimeuse,prdtimeuse,pwdtimeuse,timeuse,xtimeuse;  // 读磁盘的时间，写磁盘的时间

        access_stats query_stats;

        auto transflag = filmmem->runtimejudgetrans(&query_stats);
        int periodV = filmmem->reserveMem*1024*1024/((filmindex->valuesize+1)*sizeof(key_type)*10);
        unsigned int querynum = pointqueries.size();
        gettimeofday(&tqt1, NULL);
        for (unsigned int i = 0; i < querynum; i++) {
            if (query_stats.disknum >0 && i % periodV == 0  ){
                transflag = filmmem->runtimejudgetrans(&query_stats);
                while (transflag.first)
                {
                    runtimeevictkeytopage(filmmem,filmindex,pagedisk,&query_stats);
                    transflag = filmmem->runtimejudgetrans(&query_stats);

                }
            }

            // 执行一次point query
            pair<key_type,key_type*> res;
            key_type querykey = pointqueries[i];
//            if (querykey == 3176006677686915072 || querykey == 3004159460022892544 ||querykey ==  2919356911931049984 ||querykey == 2831290469306449920){
//                cout << "Jesus, i need You!" << endl;
//            }
            //test_interleave_baselines(filename, memThreshold, datasize, numcolumn, pagesize,errbnd);
            gettimeofday(&xt1, NULL);
            auto index_res = filmindex->search_one(querykey,&query_stats);
//            if (querykey ==  50970007)
//                cout << "i need You, my Lord" << endl;
            gettimeofday(&xt2, NULL);
            xtimeuse = (xt2.tv_sec - xt1.tv_sec) + (double) (xt2.tv_usec - xt1.tv_usec) / 1000000.0;
            query_stats.xlookuptime += xtimeuse;
            if (index_res.find == false) {   // 从 sort_piece 中 读取数据
                if (index_res.flags) {
                    query_stats.memnum += 1;
                    auto finddata = (adalru::Node<lruOff_type, key_type *> *) filmindex->sort_list->slotdata[index_res.slot];
                    res.first = querykey;
                    res.second = finddata->value;

                    gettimeofday(&plt1, NULL);
                    filmindex->sort_list->intrachain.moveTohead(finddata);// update intrachain
                    //                totalres.emplace_back(res);
                    // 更新 interchain
                    filmmem->lru->put(filmindex->sort_list->slotkey[0], filmindex->sort_list);
                    gettimeofday(&plt2, NULL);
                    pltimeuse =
                            (plt2.tv_sec - plt1.tv_sec) + (double) (plt2.tv_usec - plt1.tv_usec) / 1000000.0;
                    query_stats.lrutime += pltimeuse;
                } else {
                    query_stats.disknum += 1;
                    query_stats.diskpagenum += 1;
                    auto writeevict = (pair<pageid_type, pageOff_type> *) index_res.findleaf->slotdata[index_res.slot];  //writeevict 指向的是 pageid and offset
                    res = pagedisk->odirectreadfromdisk(
                            writeevict);    // if readfromdisk indicating doesn't use o_direct;
                    //                totalres.push_back(res);
                    index_res.findleaf->slotdata[index_res.slot] = index_res.findleaf->intrachain.put(
                            index_res.slot, res.second);
                    index_res.findleaf->slotflag[index_res.slot] = true;

                    filmmem->lru->put(index_res.findleaf->startkey, index_res.findleaf);// 更新 interchain
                    filmmem->evictPoss.emplace_back(writeevict);
                    filmindex->inkeynum++;
                    filmindex->exkeynum--;
                }
            } else {
                if (index_res.flags) {// 从内存中读数据
                    query_stats.memnum += 1;
                    auto finddata = (adalru::Node<lruOff_type, key_type *> *) index_res.findleaf->slotdata[index_res.slot];
                    res.first = querykey;
                    res.second = finddata->value;

                    gettimeofday(&plt1, NULL);
                    index_res.findleaf->intrachain.moveTohead(finddata);// update intrachain

                    // 更新 interchain
                    filmmem->lru->put(index_res.findleaf->startkey, index_res.findleaf);
                    gettimeofday(&plt2, NULL);
                    pltimeuse =
                            (plt2.tv_sec - plt1.tv_sec) + (double) (plt2.tv_usec - plt1.tv_usec) / 1000000.0;
                    query_stats.lrutime += pltimeuse;
                } else { //从磁盘读数据
                    query_stats.disknum += 1;
                    query_stats.diskpagenum += 1;
                    gettimeofday(&rdt1, NULL);
                    auto writeevict = (pair<pageid_type, pageOff_type> *) index_res.findleaf->slotdata[index_res.slot];  //writeevict 指向的是 pageid and offset

                    if (writeevict->first == filmmem->inpage->pageid) {
                        res.second = new key_type[filmindex->valuesize];
                        key_type *inmemdata = filmmem->inpage->inmemdata;
                        res.first = inmemdata[writeevict->second];
                        for (int ki = 0; ki < filmindex->valuesize; ki++) {
                            res.second[ki] = inmemdata[writeevict->second + 1 + ki];
                        }
                    } else {
                        res = pagedisk->odirectreadfromdisk(
                                writeevict);    // if readfromdisk indicating doesn't use o_direct;
                    }
                    gettimeofday(&rdt2, NULL);
                    prdtimeuse = (rdt2.tv_sec - rdt1.tv_sec) +
                                 (double) (rdt2.tv_usec - rdt1.tv_usec) / 1000000.0;
                    query_stats.rdisktime += prdtimeuse;
//                totalres.push_back(res);
                    if (res.first != querykey) {
                        cout << "not find, what's wrong? my Lord, i need You~~~~" << querykey << endl;
                    }
                    gettimeofday(&plt1, NULL);
                    index_res.findleaf->slotdata[index_res.slot] = index_res.findleaf->intrachain.put(
                            index_res.slot, res.second);
                    index_res.findleaf->slotflag[index_res.slot] = true;
                    filmmem->lru->put(index_res.findleaf->startkey, index_res.findleaf);// 更新 interchain
                    gettimeofday(&plt2, NULL);
                    pltimeuse =
                            (plt2.tv_sec - plt1.tv_sec) + (double) (plt2.tv_usec - plt1.tv_usec) / 1000000.0;
                    query_stats.lrutime += pltimeuse;
                    filmmem->evictPoss.emplace_back(writeevict);
                    filmindex->inkeynum++;
                    filmindex->exkeynum--;
                }
            }

/*
            // 执行一次range query
            vector<pair<key_type,key_type*>> rres;
            key_type firstkey = rangequeries[i][0];
            key_type lastkey = rangequeries[i][1];
            gettimeofday(&xt1, NULL);
            auto index_result = filmindex->search_range(firstkey,lastkey);
            gettimeofday(&xt2, NULL);
            xtimeuse = (xt2.tv_sec - xt1.tv_sec) + (double) (xt2.tv_usec - xt1.tv_usec) / 1000000.0;
            query_stats.xlookuptime += xtimeuse;
            pre_dict *dict = new pre_dict();
            range_prepass( filmmem,&index_result,&rres,dict,filmindex,&query_stats);
// 判断 是否需要访问磁盘
            if (dict->pagenum == 0){  // 数据都在内存
//                cout<< "Jesus, happy birthday!" << endl;
                query_stats.memnum += 1;
                // 遍历rres，释放数组
            }
            else if(rres.empty()) { //数据都在磁盘
                range_read_from_disk(&rres, dict, filmmem, filmindex, pagedisk,&query_stats);// 根据prepass 的信息，从磁盘中读数据

                query_stats.disknum += 1;
                query_stats.diskpagenum += dict->pagenum;
                //cout << "i want in Your heart, my Lord!" << endl;
            }
            else{   // 数据 在内存和磁盘中都有
                query_stats.crossnum += 1;
                query_stats.crosspagenum += dict->pagenum;
                //                cout<< "Jesus, sister nana needs You!"<<endl;
                range_read_from_disk(&rres,dict,filmmem,filmindex,pagedisk,&query_stats);
            }
            delete dict;

            // the end of range query
            */
        }


        gettimeofday(&tqt2, NULL);
        cout<<"my Lord, You are my refuge~~~~forever! "<<endl;
        timeuse = (tqt2.tv_sec - tqt1.tv_sec) + (double) (tqt2.tv_usec - tqt1.tv_usec) / 1000000.0;
        query_stats.querytimeuse += timeuse;
        query_stats.querytimeuse -= query_stats.computetimeuse;

        return query_stats;

    }


    int test_filmadaquery(unsigned int errbnd,size_t datanum,int pagesize,string dataset,double mem_threshold,double reserveMem,int recordSize,vector<key_type >keys,vector<key_type > pointqueries,vector<vector<key_type >> rangequeries, unsigned int queryn, unsigned long numkey,int datasize,double zipf,std::string workload){
        filmadatype filmada(datanum,errbnd,errbnd);
        filmadalrutype interchain(numkey);
        filmada.valuesize = recordSize-1;
        std::ostringstream osse,ossr,ossd,ossm,ossp;
        osse << errbnd;   ossp << pagesize;  ossm << mem_threshold;    ossr << (recordSize);  ossd << datanum;
        string diskpath = "/home/wamdm/chaohong/clionDir/FeasFearCPP/diskpath/";
        string file_str = diskpath+ dataset +"_filmadapages"+ossd.str()+"_"+ ossm.str()+"_"+ossp.str()+"_"+ossr.str()+"_"+osse.str();
        const char* diskfile = file_str.c_str() ;
        int numrecord = pagesize/(recordSize);
        filmstorage::filmmemory<key_type,key_type*,filmadatype, filmadalrutype,filmadadisk> memoryfilm(numkey,mem_threshold,&filmada,&interchain) ;
        memoryfilm.reserveMem = reserveMem;
        filmstorage::filmdisk<key_type> diskfilm(diskfile,pagesize,numrecord,recordSize);

        fstream fs;
        fs.open(diskfile,ios::in);
        if (fs){
            remove(diskfile);
        }
        fs.close();

        int fd = -1;
        int ret = -1;
        uint64_t file_size = 3*datasize*1024*1024ULL;
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

        key_type payload[filmada.valuesize]{};

        // 执行数据插入        // 记录数据插入的时间
        struct timeval bt1, bt2;
        double buildtimeuse;
        gettimeofday(&bt1, NULL);
//        memoryfilm.insert(keys,payload,errbnd,errbnd,&diskfilm);    // 执行数据插入
        memoryfilm.append(keys, payload, errbnd, errbnd);    // 执行数据append ，逐个更新
        gettimeofday(&bt2, NULL);
        buildtimeuse = (bt2.tv_sec - bt1.tv_sec) + (double) (bt2.tv_usec - bt1.tv_usec) / 1000000.0;

        ofstream savefile;
        savefile.open("/home/wamdm/chaohong/clionDir/updatefilm/result/adalru_film_performance.txt",ios::app);
        savefile<< "method " << "film_ada_lru "<<"available_memory " << mem_threshold+reserveMem << " qtype "<< workload << " error " << errbnd  << " pagesize "<< (pagesize*8/1024) << "k recordsize ";
        savefile<< recordSize << " build_time " << buildtimeuse << " dataset " << dataset <<" datasize " << datasize << " keynum "<<numkey << " ";
        savefile << flush;
        savefile.close();

        memoryfilm.lru->capacity = memoryfilm.lru->size;
        // 判断是否需要transfer，并执行transfer
        // 判断bulk load 之后，是否需要transfer，并执行transfer
        pair<bool, memoryusage> transflag = memoryfilm.judgetransfer();
        unsigned int transleaves;
        double ratio;
        unsigned int leaves = transflag.second.meminfo["leaves"];
        if (memoryfilm.inpage == NULL) {
            memoryfilm.createinmempage(diskfilm.pagesize, diskfilm.recordnum);    //create a page in memory;
        }

        // evitct according to lru
        while (transflag.first) {
            memoryfilm.filmtransfer(transflag.second.meminfo["totalusage"], &diskfilm);
            transflag = memoryfilm.judgetransfer();
        }


        cout<< "Jesus, thank You! finish data transfer"<< endl;
        savefile.open("/home/wamdm/chaohong/clionDir/updatefilm/result/adalru_film_performance.txt",ios::app);
        map<string,double>::iterator iter;
        savefile<<"init_write_time " << diskfilm.initwtime / 1000000.0<< " ";
        savefile<<"pages_init_write " << diskfilm.nextpageid << " ";
        savefile<< "method " << "film_ada_lru ";
        for (iter = transflag.second.meminfo.begin(); iter != transflag.second.meminfo.end(); iter++)
            savefile << iter->first << " " << iter->second << " ";

        savefile << flush;
        savefile.close();


        access_stats query_performance = filmadaquery(&filmada, pointqueries,rangequeries, queryn,&diskfilm,&memoryfilm);

        query_performance.qtype = "query";
        query_performance.zipffactor = zipf;
        query_performance.workload = workload;
        cout<< " Jesus, finished the range+point query test with film-adaptivelru, the performance is:"<<endl;
        query_performance.print_stats();

//        for (int i = 0;i<memoryfilm.index->leaflevel.size();i++){
//            delete memoryfilm.index->leaflevel[i];
//            memoryfilm.index->leaflevel[i] = NULL;
//        }
        interchain.deletelru();


        fs.open(diskfile,ios::in);
        if (fs){
            remove(diskfile);
        }
        malloc_trim(0);


        return 0;


    }


    int test_filmappending(unsigned int errbnd, size_t datanum, int pagesize, string dataset,double mem_threshold,double reserveMem,int recordSize,vector<key_type> keys,vector<key_type> queries, int queryn, size_t numkey,int datasize,double zipf){
        filmadatype filminsert(datanum,errbnd,errbnd);
        filmadalrutype interchain(numkey);
        filminsert.valuesize = recordSize-1;
        std::ostringstream osse,ossr,ossd,ossm,ossp;
        osse << errbnd;   ossp << pagesize;  ossm << mem_threshold;    ossr << (recordSize);  ossd << datanum;
        string diskpath = "/home/wamdm/chaohong/clionDir/updatefilm/diskpath/";
        string file_str = diskpath+ dataset +"_filminsertpages"+ossd.str()+"_"+ ossm.str()+"_"+ossp.str()+"_"+ossr.str()+"_"+osse.str();
        const char* diskfile = file_str.c_str() ;
        int numrecord = pagesize/(recordSize);
        filmstorage::filmmemory<key_type,key_type*,filmadatype, filmadalrutype,filmadadisk> memoryfilm(numkey,mem_threshold,&filminsert,&interchain) ;
        filmstorage::filmdisk<key_type> diskfilm(diskfile,pagesize,numrecord,recordSize);
        memoryfilm.reserveMem = reserveMem;
        fstream fs;
        fs.open(diskfile,ios::in);
        if (fs){
            fs.close();
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


//        vector<key_type> payload(filminsert.valuesize);
        key_type payload[filminsert.valuesize]{};

        // 执行数据插入        // 记录数据插入的时间
        struct timeval bt1, bt2;
        double buildtimeuse;
        gettimeofday(&bt1, NULL);
        memoryfilm.append(keys,payload,errbnd,errbnd);    // 执行数据append ，逐个更新
        gettimeofday(&bt2, NULL);
        buildtimeuse = (bt2.tv_sec - bt1.tv_sec) + (double) (bt2.tv_usec - bt1.tv_usec) / 1000000.0;
        cout << "insert time use = " << buildtimeuse << " , thank You, my Lord " << endl;
        ofstream savefile;

        savefile.open("/home/wamdm/chaohong/clionDir/updatefilm/result/adalru_film_performance.txt",ios::app);
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
        // 判断是否需要transfer，并执行transfer
        pair<bool,memoryusage> transflag = memoryfilm.judgetransfer();
        unsigned int transleaves;
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
                cout<< "Jesus, please help me!"<< endl;
            }
            memoryfilm.filmtransfer(transleaves,&diskfilm);
            transflag = memoryfilm.judgetransfer();
        }
        cout<< "Jesus, thank You! finish data transfer"<< endl;

        savefile.open("/home/wamdm/chaohong/clionDir/updatefilm/result/adalru_film_performance.txt",ios::app);
        map<string,double>::iterator iter;
        savefile<< "method " << "film_ada_lru ";
        for(iter = transflag.second.meminfo.begin(); iter != transflag.second.meminfo.end(); iter++)
            savefile<<iter->first<<" "<<iter->second<<" ";
        savefile<<"\n";
        savefile << flush;
        savefile.close();

        pair<access_stats, vector<vector<key_type>>> point_performance = filmadapointquery(&filminsert, queries, queryn,&diskfilm,&memoryfilm);
        point_performance.first.zipffactor = zipf;
        cout<< " Jesus, finished the point query test with film-adaptivelru, the performance is:"<<endl;
        point_performance.first.print_stats();
        point_performance.second.clear();
        for (int i = 0;i<memoryfilm.index->leaflevel.leafpieces.size();i++){
            delete memoryfilm.index->leaflevel.leafpieces[i];
            memoryfilm.index->leaflevel.leafpieces[i] = NULL;
        }


        fs.open(diskfile,ios::in);
        if (fs){
            fs.close();
            remove(diskfile);
        }

        return 0;



    }



    int test_interleave_insert_query_baseline(unsigned int errbnd,size_t numkey,int pagesize,string dataset,double mem_threshold,double reserveMem,
                                     int recordSize,vector<key_type>keys,unsigned long actual_numkey,int datasize,double zipf, unsigned long init_num_keys, double insert_frac,double out_of_order_frac, unsigned long batch_size,string workload) {

        filmadatype filmada(0, errbnd, errbnd);
        filmadalrutype *interchain = new filmadalrutype(actual_numkey);
        filmada.valuesize = recordSize - 1;
//        unsigned int num_out_inserts = 80000;   // 阈值：out_of_order insertion  160000/2, 160000, 80000*3, 160000*2
        std::ostringstream osse, ossr, ossd, ossm, ossp;
        osse << errbnd;
        ossp << pagesize;
        ossm << mem_threshold;
        ossr << (recordSize);
        ossd << actual_numkey;
        string diskpath = "/home/wamdm/chaohong/clionDir/updatefilm/diskpath/";
        string file_str =
                diskpath + dataset + "_filmadapages" + ossd.str() + "_" + ossm.str() + "_" + ossp.str() + "_" +
                ossr.str() + "_" + osse.str();
        const char *diskfile = file_str.c_str();
        int numrecord = pagesize / (recordSize);
        filmstorage::filmmemory<key_type, key_type *, filmadatype, filmadalrutype, filmadadisk> memoryfilm(numkey,
                                                                                                           mem_threshold,
                                                                                                           &filmada,
                                                                                                           interchain);
        memoryfilm.reserveMem = reserveMem;
        filmstorage::filmdisk<key_type> diskfilm(diskfile, pagesize, numrecord, recordSize);

        fstream fs;
        fs.open(diskfile, ios::in);
        if (fs) {
            fs.close();
            remove(diskfile);
        }

        int fd = -1;
        int ret = -1;
        uint64_t file_size = 2 * datasize * 1024 * 1024ULL;
        fd = open(diskfile, O_CREAT | O_RDWR, 0666);
        if (fd < 0) {
            printf("fd < 0");
            return -1;
        }

        //ret = fallocate(fd, 0, 0, file_size);
        ret = posix_fallocate(fd, 0, file_size);
        if (ret < 0) {
            printf("ret = %d, errno = %d,  %s\n", ret, errno, strerror(errno));
            return -1;
        }

        printf("fallocate create %.2fG file\n", file_size / 1024 / 1024 / 1024.0);
        close(fd);

        key_type payload[filmada.valuesize]{};

        //step1 bulk load the init_num_keys
        vector<key_type> init_keys;
        init_keys.assign(keys.begin(), keys.begin() + init_num_keys);

        struct timeval oribt1, oribt2;
        double buildtimeuse;
        double initwritetime = 0.0;
        gettimeofday(&oribt1, NULL);

        memoryfilm.append(init_keys, payload, errbnd, errbnd);    // 执行数据append ，逐个更新
        gettimeofday(&oribt2, NULL);
        buildtimeuse = (oribt2.tv_sec - oribt1.tv_sec) + (double) (oribt2.tv_usec - oribt1.tv_usec) / 1000000.0;
        vector<key_type>().swap(init_keys);

        // 判断bulk load 之后，是否需要transfer，并执行transfer
        pair<bool, memoryusage> transflag = memoryfilm.judgetransfer();
        unsigned int transleaves;
        double ratio;
        unsigned int leaves = transflag.second.meminfo["leaves"];
        if (memoryfilm.inpage == NULL) {
            memoryfilm.createinmempage(diskfilm.pagesize, diskfilm.recordnum);    //create a page in memory;
        }

        // evitct according to lru
        while (transflag.first) {
            initwritetime += memoryfilm.filmtransfer(transflag.second.meminfo["totalusage"], &diskfilm);
            transflag = memoryfilm.judgetransfer();
        }

        ofstream savefile;
        savefile.open("/home/wamdm/chaohong/clionDir/updatefilm/result/adalru_film_performance.txt", ios::app);
        savefile << "method " << "film_ada_lru " << "available_memory " << mem_threshold + reserveMem << " qtype "
                 << workload << " error " << errbnd << " pagesize " << (pagesize * 8 / 1024) << "k recordsize ";
        savefile << recordSize << " build_time " << buildtimeuse << " dataset " << dataset << " datasize " << datasize
                 << " keynum " << numkey << " ";
        savefile << endl;
        savefile << flush;
        savefile.close();

        savefile.open("/home/wamdm/chaohong/clionDir/updatefilm/result/adalru_film_performance.txt", ios::app);
        map<string, double>::iterator iter;
        savefile << "init_state " << "film ";
        savefile << "init_write_time " << diskfilm.initwtime / 1000000.0  << " ";
        savefile << "pages_init_write " << diskfilm.nextpageid << " ";
        savefile << "init_num_keys " << init_num_keys << " ";
        savefile << "method " << "film_ada_lru ";
        for (iter = transflag.second.meminfo.begin(); iter != transflag.second.meminfo.end(); iter++)
            savefile << iter->first << " " << iter->second << " ";
        savefile << "\n";
        savefile << flush;
        savefile.close();

        struct timeval tqt1, tqt2, xt1, xt2, plt1, plt2, prdt1, prdt2;
        double pltimeuse, prdtimeuse, pwdtimeuse, timeuse, xtimeuse;  // 读磁盘的时间，写磁盘的时间


        // Run workload
        unsigned int i = init_num_keys;
        long long cumulative_inserts = 0, cumulative_updates = 0, retrain_updates = 0, cumulative_appends = 0;
        long long cumulative_lookups = 0;
        unsigned long int num_inserts_per_batch = static_cast<unsigned long >(batch_size * insert_frac);
        unsigned long int num_lookups_per_batch =
                (batch_size - num_inserts_per_batch) ;  //
        double cumulative_insert_time = 0;
        double cumulative_lookup_time = 0;
        double cumulative_append_time = 0;
        string lookup_distribution = workload;
        struct timeval workload_start_time, workload_run_time, insert_start_time, insert_end_time,
                lookup_start_time, lookup_end_time;
        double workload_elapsed_time = 0.0, lookup_elapsed_time = 0.0, insert_elapsed_time = 0.0,  append_elapsed_time = 0.0;
        double time_limit = 3600;
//        double querycumulative_writetime = 0.0;
        double cumulative_query_writetime = 0.0;
        double cumulative_insert_writetime = 0; // 记录插入process总共导致了多少 写磁盘的时间
        gettimeofday(&workload_start_time, NULL);
        unsigned int batch_no = 0;


        access_stats query_stats;
        bool print_batch_stats = true;
        query_stats.qtype = "query";
        query_stats.zipffactor = zipf;
        query_stats.workload = workload;
        query_stats.insert_frac = insert_frac;
        query_stats.out_of_order_frac = out_of_order_frac ;

        pair<bool, double> runtime_transflag = memoryfilm.runtimejudgetrans(&query_stats);
        int periodV = memoryfilm.reserveMem * 1024 * 1024 / ((filmada.valuesize + 1) * sizeof(key_type) * 10);

        while (true) {
            batch_no++;
            // generate the search key within the inserted key space
            if (i > 0) {



//                unsigned int lookupnum = num_lookups_per_batch/2;  // divide by 2, due to one point and one range query
                unsigned int lookupnum = num_lookups_per_batch;  // divide by 2, due to one point and one range query

                key_type *lookup_keys = nullptr;  // generate lookup keys
                if (lookup_distribution == "random") {
                    lookup_keys = get_search_keys(keys, i, lookupnum);
                } else if (lookup_distribution == "zipfrandom") {
                    lookup_keys = get_search_keys_scrambledzipf(keys, i, lookupnum);
                } else if (lookup_distribution == "zipf") {
                    lookup_keys = get_search_keys_zipf(keys, i,lookupnum, zipf);
                } else if (lookup_distribution == "hotspot") {
                    lookup_keys = get_search_keys_hotspot(keys, i, lookupnum, zipf);
                } else {
                    std::cerr << "--lookup_distribution must be either 'hotspot', 'randomzipf','random' or 'zipf'"
                              << std::endl;
                    return 1;
                }




                // generate lookup ranges
                /*
                key_type **lookup_ranges = nullptr;  // generate lookup keys
                if (lookup_distribution == "random") {
                    lookup_ranges = get_search_ranges(keys, i, lookupnum);
                } else if (lookup_distribution == "zipfrandom") {
                    lookup_ranges = get_search_ranges_scrambledzipf(keys, i, lookupnum);
                } else if (lookup_distribution == "zipf") {
                    lookup_ranges = get_search_ranges_zipf(keys, i, lookupnum, zipf);
                } else if (lookup_distribution == "hotspot") {
                    lookup_ranges = get_search_ranges_hotspot(keys, i,lookupnum, zipf);
                } else {
                    std::cerr << "--lookup_distribution must be either 'hotspot', 'randomzipf','random' or 'zipf'"
                              << std::endl;
                    return 1;
                }
                 */


                // step1: do lookup;
                query_stats.computetimeuse = 0.0;
                double batch_query_writetime = 0.0;
                gettimeofday(&lookup_start_time, NULL);
                // lookup
                for (unsigned int qi = 0; qi < lookupnum; qi++) {
                    if (query_stats.disknum > 0 && qi % periodV == 0) {
                        runtime_transflag = memoryfilm.runtimejudgetrans(&query_stats);
                        while (runtime_transflag.first) {
                            batch_query_writetime += runtimeevictkeytopage2(&memoryfilm, runtime_transflag.second,
                                                                            &filmada, &diskfilm, &query_stats);
                            runtime_transflag = memoryfilm.runtimejudgetrans(&query_stats);
                        }
                    }


                    // point query

                    pair<key_type, key_type *> res;
                    key_type querykey = lookup_keys[qi];

                    gettimeofday(&xt1, NULL);
                    auto index_res = filmada.search_one(querykey, &query_stats);
                    gettimeofday(&xt2, NULL);
                    xtimeuse = (xt2.tv_sec - xt1.tv_sec) + (double) (xt2.tv_usec - xt1.tv_usec) / 1000000.0;
                    query_stats.xlookuptime += xtimeuse;
                    if (index_res.find == false) {   // 从 sort_piece 中 读取数据
                        if (index_res.flags) {
                            query_stats.memnum += 1;
                            auto finddata = (adalru::Node<lruOff_type, key_type *> *) filmada.sort_list->slotdata[index_res.slot];
                            res.first = querykey;
                            res.second = finddata->value;

                            gettimeofday(&plt1, NULL);
                            filmada.sort_list->intrachain.moveTohead(finddata);// update intrachain
                            //                totalres.emplace_back(res);
                            // 更新 interchain
                            memoryfilm.lru->put(filmada.sort_list->slotkey[0], filmada.sort_list);
                            gettimeofday(&plt2, NULL);
                            pltimeuse =
                                    (plt2.tv_sec - plt1.tv_sec) + (double) (plt2.tv_usec - plt1.tv_usec) / 1000000.0;
                            query_stats.lrutime += pltimeuse;
                        } else {
                            query_stats.disknum += 1;
                            query_stats.diskpagenum += 1;
                            auto writeevict = (pair<pageid_type, pageOff_type> *) index_res.findleaf->slotdata[index_res.slot];  //writeevict 指向的是 pageid and offset
                            res = diskfilm.odirectreadfromdisk(
                                    writeevict);    // if readfromdisk indicating doesn't use o_direct;
                            //                totalres.push_back(res);
                            index_res.findleaf->slotdata[index_res.slot] = index_res.findleaf->intrachain.put(
                                    index_res.slot, res.second);
                            index_res.findleaf->slotflag[index_res.slot] = true;

                            memoryfilm.lru->put(index_res.findleaf->startkey, index_res.findleaf);// 更新 interchain
                            memoryfilm.evictPoss.emplace_back(writeevict);
                            filmada.inkeynum++;
                            filmada.exkeynum--;
                        }
                    } else {
                        if (index_res.flags) {// 从内存中读数据
                            query_stats.memnum += 1;
                            auto finddata = (adalru::Node<lruOff_type, key_type *> *) index_res.findleaf->slotdata[index_res.slot];
                            res.first = querykey;
                            res.second = finddata->value;

                            gettimeofday(&plt1, NULL);
                            index_res.findleaf->intrachain.moveTohead(finddata);// update intrachain

                            // 更新 interchain
                            memoryfilm.lru->put(index_res.findleaf->startkey, index_res.findleaf);
                            gettimeofday(&plt2, NULL);
                            pltimeuse =
                                    (plt2.tv_sec - plt1.tv_sec) + (double) (plt2.tv_usec - plt1.tv_usec) / 1000000.0;
                            query_stats.lrutime += pltimeuse;
                        } else { //从磁盘读数据
                            query_stats.disknum += 1;
                            query_stats.diskpagenum += 1;
                            gettimeofday(&prdt1, NULL);
                            auto writeevict = (pair<pageid_type, pageOff_type> *) index_res.findleaf->slotdata[index_res.slot];  //writeevict 指向的是 pageid and offset

                            if (writeevict->first == memoryfilm.inpage->pageid) {
                                res.second = new key_type[filmada.valuesize];
                                key_type *inmemdata = memoryfilm.inpage->inmemdata;
                                res.first = inmemdata[writeevict->second];
                                for (int ki = 0; ki < filmada.valuesize; ki++) {
                                    res.second[ki] = inmemdata[writeevict->second + 1 + ki];
                                }
                            } else {
                                res = diskfilm.odirectreadfromdisk(
                                        writeevict);    // if readfromdisk indicating doesn't use o_direct;
                            }
                            gettimeofday(&prdt2, NULL);
                            prdtimeuse = (prdt2.tv_sec - prdt1.tv_sec) +
                                         (double) (prdt2.tv_usec - prdt1.tv_usec) / 1000000.0;
                            query_stats.rdisktime += prdtimeuse;
//                totalres.push_back(res);
                            if (res.first != querykey) {
                                cout << "not find, what's wrong? my Lord, i need You~~~~" << querykey << endl;
                            }
                            gettimeofday(&plt1, NULL);
                            index_res.findleaf->slotdata[index_res.slot] = index_res.findleaf->intrachain.put(
                                    index_res.slot, res.second);
                            index_res.findleaf->slotflag[index_res.slot] = true;
                            memoryfilm.lru->put(index_res.findleaf->startkey, index_res.findleaf);// 更新 interchain
                            gettimeofday(&plt2, NULL);
                            pltimeuse =
                                    (plt2.tv_sec - plt1.tv_sec) + (double) (plt2.tv_usec - plt1.tv_usec) / 1000000.0;
                            query_stats.lrutime += pltimeuse;
                            memoryfilm.evictPoss.emplace_back(writeevict);
                            filmada.inkeynum++;
                            filmada.exkeynum--;
                        }
                    }



                    // 执行一次range query
                    /*
                    vector<pair<key_type, key_type *>> rres;
                    key_type firstkey = lookup_ranges[qi][0];
                    key_type lastkey = lookup_ranges[qi][1];
//                    if (lastkey == 19.271010892963488 || lastkey == 19.271038535787397)
//                        cout<< "Jesus, please have pity on me"<< endl;

                    gettimeofday(&xt1, NULL);
                    auto index_result = filmada.search_range(firstkey, lastkey, &query_stats);
                    gettimeofday(&xt2, NULL);

                    xtimeuse = (xt2.tv_sec - xt1.tv_sec) + (double) (xt2.tv_usec - xt1.tv_usec) / 1000000.0;
                    query_stats.xlookuptime += xtimeuse;
                    pre_dict *dict = new pre_dict();
                    range_prepass(&memoryfilm, &index_result, &rres, dict, &filmada, &query_stats);
                    // 判断 是否需要访问磁盘
                    if (dict->pagenum == 0) {  // 数据都在内存
                        //          cout<< "Jesus, happy birthday!" << endl;
                        query_stats.memnum += 1;
                        // 遍历rres，释放数组
                    } else if (rres.empty()) { //数据都在磁盘
                        range_read_from_disk(&rres, dict, &memoryfilm, &filmada, &diskfilm,
                                             &query_stats);// 根据prepass 的信息，从磁盘中读数据

                        query_stats.disknum += 1;
                        query_stats.diskpagenum += dict->pagenum;
                        //cout << "i want in Your heart, my Lord!" << endl;
                    } else {   // 数据 在内存和磁盘中都有
                        query_stats.crossnum += 1;
                        query_stats.crosspagenum += dict->pagenum;
                        //                cout<< "Jesus, sister nana needs You!"<<endl;
                        range_read_from_disk(&rres, dict, &memoryfilm, &filmada, &diskfilm, &query_stats);
                    }
                    delete dict;
                    // the end of range query
                     */


                }

                // 记录lookup 的时间
                gettimeofday(&lookup_end_time, NULL);
                lookup_elapsed_time = (lookup_end_time.tv_sec - lookup_start_time.tv_sec) +
                                      (double) (lookup_end_time.tv_usec - lookup_start_time.tv_usec) / 1000000.0;
                lookup_elapsed_time -= query_stats.computetimeuse;
                cumulative_lookup_time += lookup_elapsed_time;
                cumulative_lookups += num_lookups_per_batch;
                cumulative_query_writetime += batch_query_writetime;

                delete[] lookup_keys;
//                for (unsigned int ri = 0; ri < lookupnum; ri++)
//                    delete[]lookup_ranges[ri];
//                delete[] lookup_ranges;


                // step2: Do inserts
                unsigned long int num_actual_inserts =
                        std::min(num_inserts_per_batch, (actual_numkey - i));
                unsigned long int num_keys_after_batch = i + num_actual_inserts;


                // randomly choose some keys to be inserted out-of-order
                unsigned long int num_actual_updates = 0;


                gettimeofday(&insert_start_time, NULL);
                for (; i < num_keys_after_batch; i++) {
                    filmada.append_one(keys[i], payload, errbnd, interchain);
                }
                filmada.root = &filmada.innerlevels.back()->innerpieces[0];
                auto a = filmada.leaflevel.opt->get_segment(filmada.m_tailleaf->endkey);
                filmada.m_tailleaf->update(a);

                gettimeofday(&insert_end_time, NULL);
                append_elapsed_time = (insert_end_time.tv_sec - insert_start_time.tv_sec) +
                                      (double) (insert_end_time.tv_usec - insert_start_time.tv_usec) / 1000000.0;
                cumulative_append_time += append_elapsed_time;


                // step3: judge whether to evict data, after doing inserts
                transflag = memoryfilm.judgetransfer(&query_stats);
                query_stats.computetimeuse = 0.0;
                double batch_insert_writetime = 0;
                while (transflag.first) {
//                memoryfilm.runtimefilmtransfer(&diskfilm);
                    batch_insert_writetime += runtimeevictkeytopage2(&memoryfilm, transflag.second.totalusemem,
                                                                     &filmada, &diskfilm, &query_stats);
                    transflag = memoryfilm.judgetransfer(&query_stats);
                }

                // 将写出磁盘的时间算到插入时间里，更新插入的时间
                insert_elapsed_time = append_elapsed_time;
                insert_elapsed_time += batch_insert_writetime;
                cumulative_insert_time += insert_elapsed_time;
                cumulative_inserts += num_actual_inserts;
                cumulative_insert_writetime += batch_insert_writetime;
                // the total time
                gettimeofday(&workload_run_time, NULL);
                workload_elapsed_time = (workload_run_time.tv_sec - workload_start_time.tv_sec) +
                                        (double) (workload_run_time.tv_usec - workload_start_time.tv_usec) / 1000000.0;
//            if (workload_elapsed_time > time_limit && (filmada.inkeynum + filmada.exkeynum) > actual_numkey ) {
//                break;
//            }

                savefile.open("/home/wamdm/chaohong/clionDir/updatefilm/result/adalru_film_performance.txt", ios::app);
                map<string, double>::iterator iter;
//                savefile << "_ ";
                savefile << "init_write_time " << diskfilm.initwtime / 1000000.0  << " ";
                savefile << "pages_init_write " << diskfilm.nextpageid << " ";
                savefile << "method " << "film_ada_lru ";
                for (iter = transflag.second.meminfo.begin(); iter != transflag.second.meminfo.end(); iter++)
                    savefile << iter->first << " " << iter->second << " ";
//                savefile << "\n";
                savefile << flush;
                savefile.close();
                query_stats.print_stats();
                // step4: print
                if (print_batch_stats) {
                    int num_batch_operations = num_lookups_per_batch + num_actual_inserts;
                    double batch_time = lookup_elapsed_time + insert_elapsed_time;
                    long long cumulative_operations = cumulative_lookups + cumulative_inserts;
                    double cumulative_time = cumulative_lookup_time + cumulative_insert_time;

                    std::cout << "Batch " << batch_no
                              << ", cumulative_ops: " << cumulative_operations
                              << "\n\tbatch_throughput:\t"<< num_lookups_per_batch / lookup_elapsed_time
                              << " lookups/sec,\t"<< num_actual_inserts / (insert_elapsed_time )
                              << " inserts/sec,\t"<< (num_actual_inserts ) / (append_elapsed_time)
                              << " appends/sec,\t"
                              << num_batch_operations / batch_time<< " ops/sec"
                              << num_batch_operations / (batch_time)<< " ops_with_evict/sec"
                              << "\n\tcumulative_throughput:\t"<< cumulative_lookups / cumulative_lookup_time
                              << " lookups/sec,\t"
                              << cumulative_inserts / cumulative_insert_time<< " inserts/sec,\t"
                              << cumulative_operations / cumulative_time << " ops/sec "
                              << cumulative_query_writetime << " cumulative_query_writetime "
                              << cumulative_insert_writetime << " cumulative_insert_writetime ";

                    std::cout << "\n\tcumulative_execution_1:\t"
                              << (filmada.inkeynum + filmada.exkeynum) * (filmada.valuesize + 1) * 8 /
                                 double(1024 * 1024)
                              << " G datasize, \t"
                              << filmada.inkeynum + filmada.exkeynum << " #cumulative_keys,\t"
                              << filmada.inkeynum << " #inmem_keys,\t"
                              << filmada.exkeynum << " #exmem_keys"
                              << "\n\tcumulative_execution_1:\t"
                              << cumulative_updates << " #cumulative_updates,\t"
                              << cumulative_inserts << " #cumulative_inserts,\t" << cumulative_insert_time
                              << " cumulative_insert_time,\t"
                              << cumulative_lookups << " #cumulative_lookups,\t" << cumulative_lookup_time
                              << " cumulative_lookup_time,\t"
                              << (diskfilm.initwtime / 1000000.0) + query_stats.wdisktime << " cumulative_writetime "
                              << cumulative_insert_writetime << " cumulative_insert_writetime ";

                    std::cout << std::endl;

                    savefile.open("/home/wamdm/chaohong/clionDir/updatefilm/result/adalru_film_performance.txt",
                                  ios::app);

                    savefile << batch_no << " Batch " <<
                             cumulative_operations << " cumulative_ops "
                             << "_" << " batch_throughput "
                             << num_lookups_per_batch / lookup_elapsed_time << " lookups/sec "
                             << num_actual_inserts / insert_elapsed_time << " inserts/sec "
                             << (num_actual_inserts - num_actual_updates) / append_elapsed_time << " appends/sec "
                             << num_batch_operations / batch_time << " ops/sec "
                             << batch_query_writetime << " batch_query_writetime "
                             << batch_insert_writetime << " batch_insert_writetime "
                             << "_" << " cumulative_throughput "
                             << cumulative_lookups / cumulative_lookup_time << " cum_lookups/sec "
                             << cumulative_inserts / cumulative_insert_time << " cum_inserts/sec "
                             << cumulative_operations / cumulative_time << " cum_ops/sec "
                             << cumulative_query_writetime << " cumulative_query_writetime "
                             << cumulative_insert_writetime << " cumulative_insert_writetime ";

                    savefile << (filmada.inkeynum + filmada.exkeynum) * (filmada.valuesize + 1) * 8 / double(1024 * 1024)
                             << "G datasize "
                             << filmada.inkeynum + filmada.exkeynum << " #cumulative_keys "
                             << filmada.inkeynum << " #inmem_keys "
                             << filmada.exkeynum << " #exmem_keys "
                             << cumulative_updates << " #cumulative_updates "
                             << cumulative_inserts << " #cumulative_inserts "
                             << cumulative_insert_time << " cumulative_insert_time "
                             << cumulative_lookups << " #cumulative_lookups " << cumulative_lookup_time
                             << " cumulative_lookup_time "
                             << (diskfilm.initwtime / 1000000.0) + query_stats.wdisktime << " cumulative_writetime ";

                    savefile << "\n";
                    savefile << flush;
                    savefile.close();

                }


                // Check for workload end conditions
                if ((cumulative_lookup_time + cumulative_insert_time) > time_limit || (filmada.inkeynum + filmada.exkeynum) >= numkey) {
                    break;
                }

                if (num_actual_inserts < num_inserts_per_batch) {
                    // End if we have inserted all keys in a workload with inserts
                    break;
                }
            }
            cout << "my lovely Lord, finished the interleave inserts and queries of batch " << batch_no << endl;

        }
        fs.open(diskfile,ios::in);
        if (fs){
            fs.close();
            remove(diskfile);
        }
        malloc_trim(0);
        return 0;
    }



    int test_interleave_insert_query(unsigned int errbnd,size_t numkey,int pagesize,string dataset,double mem_threshold,double reserveMem,
                                     int recordSize,vector<key_type>keys,unsigned long actual_numkey,int datasize,double zipf, unsigned long init_num_keys, double insert_frac,double out_of_order_frac, unsigned long batch_size,string workload,unsigned int threshold) {

        filmadatype filmada(0, errbnd, errbnd);
        filmadalrutype *interchain = new filmadalrutype(actual_numkey);
        filmada.valuesize = recordSize - 1;
        unsigned int num_out_inserts = threshold;   // 阈值：out_of_order insertion  160000/2, 160000, 80000*3, 160000*2
        std::ostringstream osse, ossr, ossd, ossm, ossp;
        osse << errbnd;
        ossp << pagesize;
        ossm << mem_threshold;
        ossr << (recordSize);
        ossd << actual_numkey;
        string diskpath = "/home/wamdm/chaohong/clionDir/updatefilm/diskpath/";
        string file_str =
                diskpath + dataset + "_filmadapages" + ossd.str() + "_" + ossm.str() + "_" + ossp.str() + "_" +
                ossr.str() + "_" + osse.str();
        const char *diskfile = file_str.c_str();
        int numrecord = pagesize / (recordSize);
        filmstorage::filmmemory<key_type, key_type *, filmadatype, filmadalrutype, filmadadisk> memoryfilm(numkey,
                                                                                                           mem_threshold,
                                                                                                           &filmada,
                                                                                                           interchain);
        memoryfilm.reserveMem = reserveMem;
        filmstorage::filmdisk<key_type> diskfilm(diskfile, pagesize, numrecord, recordSize);

        fstream fs;
        fs.open(diskfile, ios::in);
        if (fs) {
            fs.close();
            remove(diskfile);
        }

        int fd = -1;
        int ret = -1;
        uint64_t file_size = 2 * datasize * 1024 * 1024ULL;
        fd = open(diskfile, O_CREAT | O_RDWR, 0666);
        if (fd < 0) {
            printf("fd < 0");
            return -1;
        }

        //ret = fallocate(fd, 0, 0, file_size);
        ret = posix_fallocate(fd, 0, file_size);
        if (ret < 0) {
            printf("ret = %d, errno = %d,  %s\n", ret, errno, strerror(errno));
            return -1;
        }

        printf("fallocate create %.2fG file\n", file_size / 1024 / 1024 / 1024.0);
        close(fd);

        key_type payload[filmada.valuesize]{};

        //step1 bulk load the init_num_keys
        vector<key_type> init_keys;
        init_keys.assign(keys.begin(), keys.begin() + init_num_keys);

        struct timeval oribt1, oribt2;
        double buildtimeuse;
        double initwritetime = 0.0;
        gettimeofday(&oribt1, NULL);

        memoryfilm.append(init_keys, payload, errbnd, errbnd);    // 执行数据append ，逐个更新
        gettimeofday(&oribt2, NULL);
        buildtimeuse = (oribt2.tv_sec - oribt1.tv_sec) + (double) (oribt2.tv_usec - oribt1.tv_usec) / 1000000.0;
        vector<key_type>().swap(init_keys);

        // 判断bulk load 之后，是否需要transfer，并执行transfer
        pair<bool, memoryusage> transflag = memoryfilm.judgetransfer();
        unsigned int transleaves;
        double ratio;
        unsigned int leaves = transflag.second.meminfo["leaves"];
        if (memoryfilm.inpage == NULL) {
            memoryfilm.createinmempage(diskfilm.pagesize, diskfilm.recordnum);    //create a page in memory;
        }
        /*
         // transfer by leaves
        while (transflag.first) {
            if (filmada.leafsplit == 0) {
                ratio = memoryfilm.threshold / transflag.second.totalusemem;
                transleaves = (leaves - filmada.leafsplit) * (1 - ratio);  // the transfer number of leafs
            } else {
                ratio = memoryfilm.threshold / transflag.second.totalusemem;
                transleaves = (leaves - filmada.leafsplit) * (1 - ratio);
                if (transleaves == 0) {
                    if (errbnd >= 64) transleaves = 1;
                    else transleaves = 3;
                }
                cout << "Jesus, please help me!" << endl;
            }
            initwritetime += memoryfilm.filmtransfer(transleaves, &diskfilm);
            transflag = memoryfilm.judgetransfer();
        }
*/
        // evitct according to lru
        while (transflag.first) {
            initwritetime += memoryfilm.filmtransfer(transflag.second.meminfo["totalusage"], &diskfilm);
            transflag = memoryfilm.judgetransfer();
        }

        ofstream savefile;
        savefile.open("/home/wamdm/chaohong/clionDir/updatefilm/result/adalru_film_performance.txt", ios::app);
        savefile << "method " << "film_ada_lru " << "available_memory " << mem_threshold + reserveMem << " qtype "
                 << workload << " error " << errbnd << " pagesize " << (pagesize * 8 / 1024) << "k recordsize ";
        savefile << recordSize << " build_time " << buildtimeuse << " dataset " << dataset << " datasize " << datasize
                 << " keynum " << numkey << " ";
        savefile << endl;
        savefile << flush;
        savefile.close();

        savefile.open("/home/wamdm/chaohong/clionDir/updatefilm/result/adalru_film_performance.txt", ios::app);
        map<string, double>::iterator iter;
        savefile << "init_state " << "film ";
        savefile << "init_write_time " << diskfilm.initwtime / 1000000.0  << " ";
        savefile << "pages_init_write " << diskfilm.nextpageid << " ";
        savefile << "init_num_keys " << init_num_keys << " ";
        savefile << "method " << "film_ada_lru ";
        savefile << "threshold "<< num_out_inserts << " ";
        savefile << "out_of_order_frac "<< out_of_order_frac << " ";
        for (iter = transflag.second.meminfo.begin(); iter != transflag.second.meminfo.end(); iter++)
            savefile << iter->first << " " << iter->second << " ";
        savefile << "\n";
        savefile << flush;
        savefile.close();

        struct timeval tqt1, tqt2, xt1, xt2, plt1, plt2, prdt1, prdt2;
        double pltimeuse, prdtimeuse, pwdtimeuse, timeuse, xtimeuse;  // 读磁盘的时间，写磁盘的时间


        // Run workload
        unsigned int i = init_num_keys;
        long long cumulative_inserts = 0, cumulative_updates = 0, retrain_updates = 0, cumulative_appends = 0;
        long long cumulative_lookups = 0;
        unsigned long int num_inserts_per_batch = static_cast<unsigned long >(batch_size * insert_frac);
        unsigned long int num_lookups_per_batch =
                (batch_size - num_inserts_per_batch) ;  //
        double cumulative_insert_time = 0;
        double cumulative_lookup_time = 0;
        double cumulative_update_time = 0, cumulative_append_time = 0;
        string lookup_distribution = workload;
        struct timeval workload_start_time, workload_run_time, insert_start_time, insert_end_time,
                lookup_start_time, lookup_end_time;
        double workload_elapsed_time = 0.0, lookup_elapsed_time = 0.0, insert_elapsed_time = 0.0, update_elapsed_time = 0.0, append_elapsed_time = 0.0,rebuild_elapsed_time = 0.0;
        double time_limit = 3600;
//        double querycumulative_writetime = 0.0;
        double cumulative_query_writetime = 0.0;
        double cumulative_insert_writetime = 0; // 记录插入process总共导致了多少 写磁盘的时间
        gettimeofday(&workload_start_time, NULL);
        unsigned int batch_no = 0;


        access_stats query_stats;
        bool print_batch_stats = true;
        query_stats.qtype = "query";
        query_stats.zipffactor = zipf;
        query_stats.workload = workload;
        query_stats.insert_frac = insert_frac;
        query_stats.out_of_order_frac = out_of_order_frac ;

        pair<bool, double> runtime_transflag = memoryfilm.runtimejudgetrans(&query_stats);
        int periodV = memoryfilm.reserveMem * 1024 * 1024 / ((filmada.valuesize + 1) * sizeof(key_type) * 10);

        while (true) {
            batch_no++;
            // generate the search key within the inserted key space
            if (i > 0) {

//                key_type maxinkey = init_keys[init_num_keys-1];// the maxkey  19.270965511316433
                key_type *lookup_keys = nullptr;  // generate lookup keys
                unsigned int lookupnum = num_lookups_per_batch/2;  // divide by 2, due to one point and one range query
                if (lookup_distribution == "random") {
                    lookup_keys = get_search_keys(keys, i, lookupnum);
                } else if (lookup_distribution == "zipfrandom") {
                    lookup_keys = get_search_keys_scrambledzipf(keys, i, lookupnum);
                } else if (lookup_distribution == "zipf") {
                    lookup_keys = get_search_keys_zipf(keys, i,lookupnum, zipf);
                } else if (lookup_distribution == "hotspot") {
                    lookup_keys = get_search_keys_hotspot(keys, i, lookupnum, zipf);
                } else {
                    std::cerr << "--lookup_distribution must be either 'hotspot', 'randomzipf','random' or 'zipf'"
                              << std::endl;
                    return 1;
                }

                // generate lookup ranges
                key_type **lookup_ranges = nullptr;  // generate lookup keys
                if (lookup_distribution == "random") {
                    lookup_ranges = get_search_ranges(keys, i, lookupnum);
                } else if (lookup_distribution == "zipfrandom") {
                    lookup_ranges = get_search_ranges_scrambledzipf(keys, i, lookupnum);
                } else if (lookup_distribution == "zipf") {
                    lookup_ranges = get_search_ranges_zipf(keys, i, lookupnum, zipf);
                } else if (lookup_distribution == "hotspot") {
                    lookup_ranges = get_search_ranges_hotspot(keys, i,lookupnum, zipf);
                } else {
                    std::cerr << "--lookup_distribution must be either 'hotspot', 'zipfrandom','random' or 'zipf'"
                              << std::endl;
                    return 1;
                }

                // step1: do lookup;
                query_stats.computetimeuse = 0.0;
                double batch_query_writetime = 0.0;
                gettimeofday(&lookup_start_time, NULL);
                // lookup
                for (unsigned int qi = 0; qi < lookupnum; qi++) {
                    if (query_stats.disknum > 0 && qi % periodV == 0) {
                        runtime_transflag = memoryfilm.runtimejudgetrans(&query_stats);
                        while (runtime_transflag.first) {
                            batch_query_writetime += runtimeevictkeytopage2(&memoryfilm, runtime_transflag.second,
                                                                            &filmada, &diskfilm, &query_stats);
                            runtime_transflag = memoryfilm.runtimejudgetrans(&query_stats);
                        }
                    }

                    // point query
                    pair<key_type, key_type *> res;
                    key_type querykey = lookup_keys[qi];

                    gettimeofday(&xt1, NULL);
                    auto index_res = filmada.search_one(querykey, &query_stats);
                    gettimeofday(&xt2, NULL);
                    xtimeuse = (xt2.tv_sec - xt1.tv_sec) + (double) (xt2.tv_usec - xt1.tv_usec) / 1000000.0;
                    query_stats.xlookuptime += xtimeuse;
                    if (index_res.find == false) {   // 从 sort_piece 中 读取数据
                        if (index_res.flags) {
                            query_stats.memnum += 1;
                            auto finddata = (adalru::Node<lruOff_type, key_type *> *) filmada.sort_list->slotdata[index_res.slot];
                            res.first = querykey;
                            res.second = finddata->value;

                            gettimeofday(&plt1, NULL);
                            filmada.sort_list->intrachain.moveTohead(finddata);// update intrachain
                            //                totalres.emplace_back(res);
                            // 更新 interchain
                            memoryfilm.lru->put(filmada.sort_list->slotkey[0], filmada.sort_list);
                            gettimeofday(&plt2, NULL);
                            pltimeuse =
                                    (plt2.tv_sec - plt1.tv_sec) + (double) (plt2.tv_usec - plt1.tv_usec) / 1000000.0;
                            query_stats.lrutime += pltimeuse;
                        } else {
                            query_stats.disknum += 1;
                            query_stats.diskpagenum += 1;
                            auto writeevict = (pair<pageid_type, pageOff_type> *) index_res.findleaf->slotdata[index_res.slot];  //writeevict 指向的是 pageid and offset
                            res = diskfilm.odirectreadfromdisk(
                                    writeevict);    // if readfromdisk indicating doesn't use o_direct;
                            //                totalres.push_back(res);
                            index_res.findleaf->slotdata[index_res.slot] = index_res.findleaf->intrachain.put(
                                    index_res.slot, res.second);
                            index_res.findleaf->slotflag[index_res.slot] = true;

                            memoryfilm.lru->put(index_res.findleaf->startkey, index_res.findleaf);// 更新 interchain
                            memoryfilm.evictPoss.emplace_back(writeevict);
                            filmada.inkeynum++;
                            filmada.exkeynum--;
                        }
                    } else {
                        if (index_res.flags) {// 从内存中读数据
                            query_stats.memnum += 1;
                            auto finddata = (adalru::Node<lruOff_type, key_type *> *) index_res.findleaf->slotdata[index_res.slot];
                            res.first = querykey;
                            res.second = finddata->value;

                            gettimeofday(&plt1, NULL);
                            index_res.findleaf->intrachain.moveTohead(finddata);// update intrachain

                            // 更新 interchain
                            memoryfilm.lru->put(index_res.findleaf->startkey, index_res.findleaf);
                            gettimeofday(&plt2, NULL);
                            pltimeuse =
                                    (plt2.tv_sec - plt1.tv_sec) + (double) (plt2.tv_usec - plt1.tv_usec) / 1000000.0;
                            query_stats.lrutime += pltimeuse;
                        } else { //从磁盘读数据
                            query_stats.disknum += 1;
                            query_stats.diskpagenum += 1;
                            gettimeofday(&prdt1, NULL);
                            auto writeevict = (pair<pageid_type, pageOff_type> *) index_res.findleaf->slotdata[index_res.slot];  //writeevict 指向的是 pageid and offset

                            if (writeevict->first == memoryfilm.inpage->pageid) {
                                res.second = new key_type[filmada.valuesize];
                                key_type *inmemdata = memoryfilm.inpage->inmemdata;
                                res.first = inmemdata[writeevict->second];
                                for (int ki = 0; ki < filmada.valuesize; ki++) {
                                    res.second[ki] = inmemdata[writeevict->second + 1 + ki];
                                }
                            } else {
                                res = diskfilm.odirectreadfromdisk(
                                        writeevict);    // if readfromdisk indicating doesn't use o_direct;
                            }
                            gettimeofday(&prdt2, NULL);
                            prdtimeuse = (prdt2.tv_sec - prdt1.tv_sec) +
                                         (double) (prdt2.tv_usec - prdt1.tv_usec) / 1000000.0;
                            query_stats.rdisktime += prdtimeuse;
//                totalres.push_back(res);
                            if (res.first != querykey) {
                                cout << "not find, what's wrong? my Lord, i need You~~~~" << querykey << endl;
                            }
                            gettimeofday(&plt1, NULL);
                            index_res.findleaf->slotdata[index_res.slot] = index_res.findleaf->intrachain.put(
                                    index_res.slot, res.second);
                            index_res.findleaf->slotflag[index_res.slot] = true;
                            memoryfilm.lru->put(index_res.findleaf->startkey, index_res.findleaf);// 更新 interchain
                            gettimeofday(&plt2, NULL);
                            pltimeuse =
                                    (plt2.tv_sec - plt1.tv_sec) + (double) (plt2.tv_usec - plt1.tv_usec) / 1000000.0;
                            query_stats.lrutime += pltimeuse;
                            memoryfilm.evictPoss.emplace_back(writeevict);
                            filmada.inkeynum++;
                            filmada.exkeynum--;
                        }
                    }

                    // 执行一次range query
                    vector<pair<key_type, key_type *>> rres;
                    key_type firstkey = lookup_ranges[qi][0];
                    key_type lastkey = lookup_ranges[qi][1];
//                    if (lastkey == 19.271010892963488 || lastkey == 19.271038535787397)
//                        cout<< "Jesus, please have pity on me"<< endl;

                    gettimeofday(&xt1, NULL);
                    auto index_result = filmada.search_range(firstkey, lastkey, &query_stats);
                    gettimeofday(&xt2, NULL);

                    xtimeuse = (xt2.tv_sec - xt1.tv_sec) + (double) (xt2.tv_usec - xt1.tv_usec) / 1000000.0;
                    query_stats.xlookuptime += xtimeuse;
                    pre_dict *dict = new pre_dict();
                    range_prepass(&memoryfilm, &index_result, &rres, dict, &filmada, &query_stats);
                    // determine whether to access disk
                    if (dict->pagenum == 0) {  // the request data  数据都在内存
                        //          cout<< "Jesus, happy birthday!" << endl;
                        query_stats.memnum += 1;
                        // 遍历rres，释放数组
                    } else if (rres.empty()) { //数据都在磁盘
                        range_read_from_disk(&rres, dict, &memoryfilm, &filmada, &diskfilm,
                                             &query_stats);// read data from disk according to 根据prepass 的信息，从磁盘中读数据

                        query_stats.disknum += 1;
                        query_stats.diskpagenum += dict->pagenum;
                        //cout << "i want in Your heart, my Lord!" << endl;
                    } else {   //the request data are in both memory and disk
                        query_stats.crossnum += 1;
                        query_stats.crosspagenum += dict->pagenum;
                        //                cout<< "Jesus, sister nana needs You!"<<endl;
                        range_read_from_disk(&rres, dict, &memoryfilm, &filmada, &diskfilm, &query_stats);
                    }
                    delete dict;
                    // the end of range query
                }

                // record the start time of lookup
                gettimeofday(&lookup_end_time, NULL);
                lookup_elapsed_time = (lookup_end_time.tv_sec - lookup_start_time.tv_sec) +
                                      (double) (lookup_end_time.tv_usec - lookup_start_time.tv_usec) / 1000000.0;
                lookup_elapsed_time -= query_stats.computetimeuse;
                cumulative_lookup_time += lookup_elapsed_time;
                cumulative_lookups += num_lookups_per_batch;
                cumulative_query_writetime += batch_query_writetime;

                delete[] lookup_keys;
                for (unsigned int ri = 0; ri < lookupnum; ri++)
                    delete[]lookup_ranges[ri];
                delete[] lookup_ranges;


                // step2: Do inserts
                unsigned long int num_actual_inserts =
                        std::min(num_inserts_per_batch, (actual_numkey - i));
                unsigned long int num_keys_after_batch = i + num_actual_inserts;

                // randomly choose some keys to be inserted out-of-order
                unsigned long int num_actual_updates = num_actual_inserts * out_of_order_frac;
                vector<key_type> update_vec, append_vec;
                update_vec.reserve(num_actual_updates + 200);
                append_vec.reserve(num_actual_inserts);
                std::random_device rd; // obtain a random number from hardware
                std::mt19937_64 gen(rd()); // seed the generator
                auto upper = 1 / out_of_order_frac + 1;
                std::uniform_real_distribution<> dist(1, upper); // define the range  [a, b)
                key_type maxinkey = keys[num_keys_after_batch - 1];
                if (out_of_order_frac == 0.0){
                    append_vec.assign(keys.begin()+i, keys.begin() + num_keys_after_batch);
                    i += num_actual_inserts;
                }
                else{
                    for (; i < num_keys_after_batch; i++) {
                        if (1 == static_cast<int>(dist(gen)))
                            update_vec.push_back(keys[i]);
                        else
                            append_vec.emplace_back(keys[i]);
                    }
                }

                auto ulen = update_vec.size();
                auto alen = append_vec.size();

                gettimeofday(&insert_start_time, NULL);
                for (unsigned int ai = 0; ai < append_vec.size(); ai++) {
                    filmada.append_one(append_vec[ai], payload, errbnd, interchain);
                }
                filmada.root = &filmada.innerlevels.back()->innerpieces[0];
                auto a = filmada.leaflevel.opt->get_segment(filmada.m_tailleaf->endkey);
                filmada.m_tailleaf->update(a);

                gettimeofday(&insert_end_time, NULL);
                append_elapsed_time = (insert_end_time.tv_sec - insert_start_time.tv_sec) +
                                      (double) (insert_end_time.tv_usec - insert_start_time.tv_usec) / 1000000.0;
                cumulative_append_time += append_elapsed_time;
                cumulative_updates += ulen;
                retrain_updates += ulen;
                /*
                // if out_of_order insertion, judge whether arrive at the threshold to rebuild
                 */

                if (out_of_order_frac) {
                    gettimeofday(&insert_start_time, NULL);   // num_out_inserts
                    if (filmada.sort_list == NULL && update_vec.size() > num_out_inserts ||
                            retrain_updates > num_out_inserts) {// rebuild FILM
                        retrain_updates = 0;
                        // record the start time of rebuild
                        gettimeofday(&lookup_start_time, NULL);
                        filmada.release();
                        interchain->deletelru();
                        filmadalrutype *newinterchain = new filmadalrutype(actual_numkey);
                        interchain = newinterchain;
                        memoryfilm.lru = interchain;
                        init_keys.assign(keys.begin(), keys.begin() + i);
                        memoryfilm.append(init_keys, payload, errbnd, errbnd);    // 执行数据append
                        vector<key_type>().swap(init_keys);
                        if (memoryfilm.inpage != NULL) {
                            memoryfilm.inpageid = 0;
                            delete memoryfilm.inpage;
                            memoryfilm.createinmempage(diskfilm.pagesize, diskfilm.recordnum);
                            vector<pair<pageid_type, pageOff_type>>().swap(memoryfilm.evicttable);
                            vector<pair<pageid_type, pageOff_type> *>().swap(memoryfilm.evictPoss);
                            memoryfilm.evicttable.reserve(numkey);
                            diskfilm.nextpageid = 0;
                        }
                        // record rebuild end time
                        gettimeofday(&lookup_end_time, NULL);
                        rebuild_elapsed_time = (lookup_end_time.tv_sec - lookup_start_time.tv_sec) +
                                              (double) (lookup_end_time.tv_usec - lookup_start_time.tv_usec) / 1000000.0;
                        cout << "thank you, my Lovely Lord, retrain the model" << endl;
                    } else
                        filmada.update_random(update_vec, payload, interchain);
                    filmada.root = &filmada.innerlevels.back()->innerpieces[0];
                    a = filmada.leaflevel.opt->get_segment(filmada.m_tailleaf->endkey);
                    filmada.m_tailleaf->update(a);
                    gettimeofday(&insert_end_time, NULL);
                    update_elapsed_time = (insert_end_time.tv_sec - insert_start_time.tv_sec) +
                                          (double) (insert_end_time.tv_usec - insert_start_time.tv_usec) / 1000000.0;
                    cumulative_update_time += update_elapsed_time;
                }



                vector<key_type>().swap(update_vec);
                vector<key_type>().swap(append_vec);

                // step3: judge whether to evict data, after doing inserts
                transflag = memoryfilm.judgetransfer(&query_stats);
                query_stats.computetimeuse = 0.0;
                double batch_insert_writetime = 0;
                while (transflag.first) {
//                memoryfilm.runtimefilmtransfer(&diskfilm);
                    batch_insert_writetime += runtimeevictkeytopage2(&memoryfilm, transflag.second.totalusemem,
                                                                     &filmada, &diskfilm, &query_stats);
                    transflag = memoryfilm.judgetransfer(&query_stats);
                }

                // 将写出磁盘的时间算到插入时间里，更新插入的时间
                insert_elapsed_time = update_elapsed_time + append_elapsed_time;
                insert_elapsed_time += batch_insert_writetime;
//                insert_elapsed_time -= query_stats.computetimeuse;
                cumulative_insert_time += insert_elapsed_time;
                cumulative_inserts += num_actual_inserts;
                cumulative_appends += alen;
                cumulative_insert_writetime += batch_insert_writetime;
                // the total time
                gettimeofday(&workload_run_time, NULL);
                workload_elapsed_time = (workload_run_time.tv_sec - workload_start_time.tv_sec) +
                                        (double) (workload_run_time.tv_usec - workload_start_time.tv_usec) / 1000000.0;
//            if (workload_elapsed_time > time_limit && (filmada.inkeynum + filmada.exkeynum) > actual_numkey ) {
//                break;
//            }

                savefile.open("/home/wamdm/chaohong/clionDir/updatefilm/result/adalru_film_performance.txt", ios::app);
                map<string, double>::iterator iter;
//                savefile << "_ ";
                savefile << "init_write_time " << diskfilm.initwtime / 1000000.0  << " ";
                savefile << "pages_init_write " << diskfilm.nextpageid << " ";
                savefile << "method " << "film_ada_lru ";
                for (iter = transflag.second.meminfo.begin(); iter != transflag.second.meminfo.end(); iter++)
                    savefile << iter->first << " " << iter->second << " ";
//                savefile << "\n";
                savefile << flush;
                savefile.close();
                query_stats.print_stats();
                // step4: print 输出
                if (print_batch_stats) {
                    int num_batch_operations = num_lookups_per_batch + num_actual_inserts;
                    double batch_time = lookup_elapsed_time + insert_elapsed_time;
                    long long cumulative_operations = cumulative_lookups + cumulative_inserts;
                    double cumulative_time = cumulative_lookup_time + cumulative_insert_time;

                    std::cout << "Batch " << batch_no
                              << ", cumulative_ops: " << cumulative_operations
                              << "\n\tbatch_throughput:\t"
                              << num_lookups_per_batch / lookup_elapsed_time
                              << " lookups/sec,\t"
                              << num_actual_inserts / (insert_elapsed_time )
                              << " inserts/sec,\t"
                              << (num_actual_inserts - num_actual_updates) / (append_elapsed_time)
                              << " appends/sec,\t"
                              << num_actual_updates / update_elapsed_time
                              << " updates/sec,\t" << num_batch_operations / batch_time
                              << " ops/sec" << num_batch_operations / (batch_time)
                              << " ops_with_evict/sec"
                              << "\n\tcumulative_throughput:\t"
                              << cumulative_lookups / cumulative_lookup_time
                              << " lookups/sec,\t"
                              << cumulative_inserts / cumulative_insert_time
                              << " inserts/sec,\t"
                              << cumulative_updates / cumulative_update_time
                              << " updates/sec,\t"
                              << cumulative_operations / cumulative_time << " ops/sec "
                              << cumulative_query_writetime << " cumulative_query_writetime "
                              << cumulative_insert_writetime << " cumulative_insert_writetime "
                              << update_elapsed_time << " update_elapsed_time "
                              << cumulative_update_time << " cumulative_update_time "
                              << rebuild_elapsed_time << " rebuild_elapsed_time ";

                    std::cout << "\n\tcumulative_execution_1:\t"
                              << (filmada.inkeynum + filmada.exkeynum) * (filmada.valuesize + 1) * 8 /
                                 double(1024 * 1024)
                              << " G datasize, \t"
                              << filmada.inkeynum + filmada.exkeynum << " #cumulative_keys,\t"
                              << filmada.inkeynum << " #inmem_keys,\t"
                              << filmada.exkeynum << " #exmem_keys"
                              << "\n\tcumulative_execution_1:\t"
                              << cumulative_updates << " #cumulative_updates,\t" << cumulative_update_time
                              << " cumulative_update_time,\t"
                              << cumulative_inserts << " #cumulative_inserts,\t" << cumulative_insert_time
                              << " cumulative_insert_time,\t"
                              << cumulative_lookups << " #cumulative_lookups,\t" << cumulative_lookup_time
                              << " cumulative_lookup_time,\t"
                              << (diskfilm.initwtime / 1000000.0) + query_stats.wdisktime << " cumulative_writetime "
                              << cumulative_insert_writetime << " cumulative_insert_writetime "
                              << update_elapsed_time << " update_elapsed_time "
                              << cumulative_update_time << " cumulative_update_time "
                              << rebuild_elapsed_time << " rebuild_elapsed_time ";

                    std::cout << std::endl;

                    savefile.open("/home/wamdm/chaohong/clionDir/updatefilm/result/adalru_film_performance.txt",ios::app);

                    savefile << batch_no << " Batch " <<
                             cumulative_operations << " cumulative_ops "
                             << "_" << " batch_throughput "
                             << num_lookups_per_batch / lookup_elapsed_time << " lookups/sec "
                             << num_actual_inserts / insert_elapsed_time << " inserts/sec "
                             << (num_actual_inserts - num_actual_updates) / append_elapsed_time << " appends/sec "
                             << num_actual_updates / update_elapsed_time << " updates/sec "
                             << num_batch_operations / batch_time << " ops/sec "
                             << batch_query_writetime << " batch_query_writetime "
                             << batch_insert_writetime << " batch_insert_writetime "
                             << "_" << " cumulative_throughput "
                             << cumulative_lookups / cumulative_lookup_time << " cum_lookups/sec "
                             << cumulative_inserts / cumulative_insert_time << " cum_inserts/sec "
                             << cumulative_updates / cumulative_update_time << " cum_updates/sec "
                             << cumulative_operations / cumulative_time << " cum_ops/sec "
                             << cumulative_query_writetime << " cumulative_query_writetime "
                             << cumulative_insert_writetime << " cumulative_insert_writetime ";

                    savefile << (filmada.inkeynum + filmada.exkeynum) * (filmada.valuesize + 1) * 8 / double(1024 * 1024)
                            << "G datasize "
                            << filmada.inkeynum + filmada.exkeynum << " #cumulative_keys "
                            << filmada.inkeynum << " #inmem_keys "
                            << filmada.exkeynum << " #exmem_keys "
                            << cumulative_updates << " #cumulative_updates "
                            << retrain_updates << " #retrain_updates "
                            << cumulative_update_time<< " cumulative_update_time "
                            << cumulative_inserts << " #cumulative_inserts "
                            << cumulative_insert_time << " cumulative_insert_time "
                            << cumulative_lookups << " #cumulative_lookups " << cumulative_lookup_time
                            << " cumulative_lookup_time "
                            << (diskfilm.initwtime / 1000000.0) + query_stats.wdisktime << " cumulative_writetime ";

                    savefile << "\n";
                    savefile << flush;
                    savefile.close();

                }


                // Check for workload end conditions
                if (workload_elapsed_time > time_limit || (filmada.inkeynum + filmada.exkeynum) >= numkey) {
                    break;
                }

                if (num_actual_inserts < num_inserts_per_batch) {
                    // End if we have inserted all keys in a workload with inserts
                    break;
                }
            }
            cout << "my lovely Lord, finished the interleave inserts and queries of batch " << batch_no << endl;

        }
        fs.open(diskfile,ios::in);
        if (fs){
            fs.close();
            remove(diskfile);
        }
        malloc_trim(0);
        return 0;
    }


}


#endif //FILMINSERT_FILM_H
