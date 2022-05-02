//
// Created by CCMa on 2022/1/28.
//

#ifndef FILMINSERT_PWLF_H
#define FILMINSERT_PWLF_H



#include <cmath>
#include <limits>
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <climits>
//#include <algorithm>
#include <omp.h>
#include <iostream>
#include "film.h"



using namespace std;
namespace filminsert::internal{
//    std::conditional_t< 条件std::is_floating_point_v<T>, 条件为真 则 long double,条件为假 则 std::conditional_t<(sizeof(T) < 8), int64_t, __int128>>;
    template<typename T>    // std::is_floating_point_v<T>  检查T 是否为浮点类型，是则为true， else false
    using LargeSigned = typename std::conditional_t<std::is_floating_point_v<T>,
            long double,
            std::conditional_t<(sizeof(T) <8), int64_t, __int128>>;  // chaochao modify  "< " to <=

    template<typename X, typename Y>
    class insertPWLF{
        using SX = LargeSigned<X>;
        using SY = LargeSigned<Y>;
        struct Slope {
            SX dx{};
            SY dy{};

            bool operator<(const Slope &p) const { return dy * p.dx < dx * p.dy; }  // operator 运算符重载
            bool operator>(const Slope &p) const { return dy * p.dx > dx * p.dy; }
            bool operator==(const Slope &p) const { return dy * p.dx == dx * p.dy; }
            bool operator!=(const Slope &p) const { return dy * p.dx != dx * p.dy; }
            explicit operator long double() const { return dy / (long double) dx; }
        };

        struct Point {
            X x{};
            Y y{};

            Slope operator-(const Point &p) const { return {SX(x) - p.x, SY(y) - p.y}; }
        };

    public:
        const Y epsilon;
        std::vector<Point> lower;
        std::vector<Point> upper;
        X first_x = 0;
        X last_x = 0;
        size_t lower_start = 0;
        size_t upper_start = 0;
        size_t points_in_hull = 0;
        Point rectangle[4];
        queue< Point > buffqueue;   // buffqueue.size() 返回队列的个数，buffqueue.empty() 判断队列是否为空



        auto cross(const Point &O, const Point &A, const Point &B) const {
            auto OA = A - O;
            auto OB = B - O;

//        long int tmp = OA.dx * OB.dy - OA.dy * OB.dx;
//        std::cout<< tmp <<std::endl;
            return OA.dx * OB.dy - OA.dy * OB.dx;
        }

    public:

        class CanonicalSegment;

        explicit insertPWLF(Y epsilon) : epsilon(epsilon), lower(), upper() {
            if (epsilon < 0)
                throw std::invalid_argument("epsilon cannot be negative");

            upper.reserve(1u << 16);  // reserve, 告知容器应该准备保存多少元素，不改变容器中元素的数量，仅影响容器预先分配多少的内存空间
            lower.reserve(1u << 16);

        }
        explicit insertPWLF() : epsilon(), lower(), upper() {

            upper.reserve(1u << 16);  // reserve, 告知容器应该准备保存多少元素，不改变容器中元素的数量，仅影响容器预先分配多少的内存空间
            lower.reserve(1u << 16);
        }

        bool add_point(const X &x, const Y &y) {
//        if (points_in_hull > 0 && x <= last_x)
//            throw std::logic_error("Points must be increasing by x.");

            last_x = x;
            constexpr auto max_y = std::numeric_limits<Y>::max();
            constexpr auto min_y = std::numeric_limits<Y>::lowest();
            Point p1{x, y >= max_y - epsilon ? max_y : y + epsilon};
            Point p2{x, y <= min_y + epsilon ? min_y : y - epsilon};

            if (points_in_hull == 0) {
                first_x = x;
                rectangle[0] = p1;
                rectangle[1] = p2;
                upper.clear();
                lower.clear();
                upper.push_back(p1);
                lower.push_back(p2);
                upper_start = lower_start = 0;
                ++points_in_hull;
                return true;
            }

            if (points_in_hull == 1) {
                rectangle[2] = p2;
                rectangle[3] = p1;
                upper.push_back(p1);
                lower.push_back(p2);
                ++points_in_hull;
                return true;
            }

            auto slope1 = rectangle[2] - rectangle[0];
            auto slope2 = rectangle[3] - rectangle[1];
            bool outside_line1 = p1 - rectangle[2] < slope1;
            bool outside_line2 = p2 - rectangle[3] > slope2;

            if (outside_line1 || outside_line2) {
                points_in_hull = 0;
                return false;
            }

            if (p1 - rectangle[1] < slope2) {
                // Find extreme slope
                auto min = lower[lower_start] - p1;
                auto min_i = lower_start;
                for (auto i = lower_start + 1; i < lower.size(); i++) {
                    auto val = lower[i] - p1;
                    if (val > min)
                        break;
                    min = val;
                    min_i = i;
                }

                rectangle[1] = lower[min_i];
                rectangle[3] = p1;
                lower_start = min_i;

                // Hull update
                auto end = upper.size();
                for (; end >= upper_start + 2 && cross(upper[end - 2], upper[end - 1], p1) <= 0; --end)
                    continue;
                upper.resize(end);
                upper.push_back(p1);
            }

            if (p2 - rectangle[0] > slope1) {
                // Find extreme slope
                auto max = upper[upper_start] - p2;
                auto max_i = upper_start;
                for (auto i = upper_start + 1; i < upper.size(); i++) {
                    auto val = upper[i] - p2;
                    if (val < max)
                        break;
                    max = val;
                    max_i = i;
                }

                rectangle[0] = upper[max_i];
                rectangle[2] = p2;
                upper_start = max_i;

                // Hull update
                auto end = lower.size();
                for (; end >= lower_start + 2 && cross(lower[end - 2], lower[end - 1], p2) >= 0; --end)
                    continue;
                lower.resize(end);
                lower.push_back(p2);
            }

            ++points_in_hull;
            return true;
        }

        bool append_point(const X &x, const Y &y) {
//        if (points_in_hull > 0 && x <= last_x)
//            throw std::logic_error("Points must be increasing by x.");

            last_x = x;
            constexpr auto max_y = std::numeric_limits<Y>::max();
            constexpr auto min_y = std::numeric_limits<Y>::lowest();
            Point p1{x, y >= max_y - epsilon ? max_y : y + epsilon};
            Point p2{x, y <= min_y + epsilon ? min_y : y - epsilon};

            if (points_in_hull == 0) {
                buffqueue.push(p1);
                buffqueue.push(p2);
                if (buffqueue.size() == 4){
                    upper_start = lower_start = 0;
                    upper.clear();
                    lower.clear();
                    rectangle[0] =  buffqueue.front();
                    first_x = rectangle[0].x;
                    buffqueue.pop();
                    rectangle[1] = buffqueue.front();
                    buffqueue.pop();
                    rectangle[3] =  buffqueue.front();
                    buffqueue.pop();
                    rectangle[2] =  buffqueue.front();
                    buffqueue.pop();
                    upper.push_back(rectangle[0]);
                    lower.push_back(rectangle[1]);
                    upper.push_back(rectangle[3]);
                    lower.push_back(rectangle[2] );
                    points_in_hull += 2;
                }
//                else{
//                    cout<< "my Lord, i need You! please have pity on me!!"  << endl;
//                }

                return true;
            }


            auto slope1 = rectangle[2] - rectangle[0];
            auto slope2 = rectangle[3] - rectangle[1];
            bool outside_line1 = p1 - rectangle[2] < slope1;
            bool outside_line2 = p2 - rectangle[3] > slope2;

            if (outside_line1 || outside_line2) {
                points_in_hull = 0;
                return false;
            }

            if (p1 - rectangle[1] < slope2) {
                // Find extreme slope
                auto min = lower[lower_start] - p1;
                auto min_i = lower_start;
                for (auto i = lower_start + 1; i < lower.size(); i++) {
                    auto val = lower[i] - p1;
                    if (val > min)
                        break;
                    min = val;
                    min_i = i;
                }

                rectangle[1] = lower[min_i];
                rectangle[3] = p1;
                lower_start = min_i;

                // Hull update
                auto end = upper.size();
                for (; end >= upper_start + 2 && cross(upper[end - 2], upper[end - 1], p1) <= 0; --end)
                    continue;
                upper.resize(end);
                upper.push_back(p1);
            }

            if (p2 - rectangle[0] > slope1) {
                // Find extreme slope
                auto max = upper[upper_start] - p2;
                auto max_i = upper_start;
                for (auto i = upper_start + 1; i < upper.size(); i++) {
                    auto val = upper[i] - p2;
                    if (val < max)
                        break;
                    max = val;
                    max_i = i;
                }

                rectangle[0] = upper[max_i];
                rectangle[2] = p2;
                upper_start = max_i;

                // Hull update
                auto end = lower.size();
                for (; end >= lower_start + 2 && cross(lower[end - 2], lower[end - 1], p2) >= 0; --end)
                    continue;
                lower.resize(end);
                lower.push_back(p2);
            }

//            ++points_in_hull;
            return true;
        }


        CanonicalSegment get_segment() {
            if (points_in_hull == 1)
                return CanonicalSegment(rectangle[0], rectangle[1], first_x);
            return CanonicalSegment(rectangle, first_x);
        }

        CanonicalSegment get_segment( X break_x) {
            if (points_in_hull == 1)
                return CanonicalSegment(rectangle[0], rectangle[1], first_x);
            return CanonicalSegment(rectangle, first_x,break_x);
        }

        CanonicalSegment get_segment( X break_x,std::vector<X> keys_vec) {
            if (points_in_hull == 1)
                return CanonicalSegment(rectangle[0], rectangle[1], first_x);
            return CanonicalSegment(rectangle, first_x,break_x,keys_vec);
        }


        void reset() {
            points_in_hull = 0;
            lower.clear();
            upper.clear();
        }
    };

    template<typename X, typename Y>
    class insertPWLF<X, Y>::CanonicalSegment {
        friend class insertPWLF;


    public:
        std::vector<X> slotkey;
        Point rectangle[4];
        X first;
        X last;

        CanonicalSegment( X first) : first(first) {};
        CanonicalSegment(const Point &p0, const Point &p1, X first) : rectangle{p0, p1, p0, p1}, first(first) {};

        CanonicalSegment(const Point &p0, const Point &p1, X first, X last) : rectangle{p0, p1, p0, p1}, first(first) ,last(last) {};

        CanonicalSegment(const Point (&rectangle)[4], X first)
                : rectangle{rectangle[0], rectangle[1], rectangle[2], rectangle[3]}, first(first) {};

        CanonicalSegment(const Point (&rectangle)[4], X first, X last)
                : rectangle{rectangle[0], rectangle[1], rectangle[2], rectangle[3]}, first(first), last(last) {};



        bool one_point() const {
            return rectangle[0].x == rectangle[2].x && rectangle[0].y == rectangle[2].y
                   && rectangle[1].x == rectangle[3].x && rectangle[1].y == rectangle[3].y;
        }

    public:

        CanonicalSegment() = default;

        X get_first_x() const { return first; }
        X get_last_x() const { return last; }

        std::pair<long double, long double> get_intersection() const {
            auto &p0 = rectangle[0];
            auto &p1 = rectangle[1];
            auto &p2 = rectangle[2];
            auto &p3 = rectangle[3];
            auto slope1 = p2 - p0;
            auto slope2 = p3 - p1;

            if (one_point() || slope1 == slope2)
                return {p0.x, p0.y};

            auto p0p1 = p1 - p0;
            auto a = slope1.dx * slope2.dy - slope1.dy * slope2.dx;
            auto b = (p0p1.dx * slope2.dy - p0p1.dy * slope2.dx) / static_cast<long double>(a);
            auto i_x = p0.x + b * slope1.dx;
            auto i_y = p0.y + b * slope1.dy;
            return {i_x, i_y};
        }

        std::pair<long double, SY> get_floating_point_segment(const X &origin) const {
            if (one_point())
                return {0, (rectangle[0].y + rectangle[1].y) / 2};

            if constexpr (std::is_integral_v<X> && std::is_integral_v<Y>) {
                auto slope = rectangle[3] - rectangle[1];
                auto intercept_n = slope.dy * (SX(origin) - rectangle[1].x);
                auto intercept_d = slope.dx;
                auto rounding_term = ((intercept_n < 0) ^ (intercept_d < 0) ? -1 : +1) * intercept_d / 2;
                auto intercept = (intercept_n + rounding_term) / intercept_d + rectangle[1].y;
                return {static_cast<long double>(slope), intercept};
            }

            auto[i_x, i_y] = get_intersection();
            auto[min_slope, max_slope] = get_slope_range();
            auto slope = (min_slope + max_slope) / 2.;
            auto intercept = i_y - (i_x - origin) * slope;
            return {slope, intercept};
        }

        std::pair<long double, long double> get_slope_range() const {
            if (one_point())
                return {0, 1};

            auto min_slope = static_cast<long double>(rectangle[2] - rectangle[0]);
            auto max_slope = static_cast<long double>(rectangle[3] - rectangle[1]);
            return {min_slope, max_slope};
        }
    };


    template<typename key_type,typename filmtype>
    std::pair<size_t,std::vector<key_type> > append_segmentation(size_t error,std::vector<key_type> keys,
                                                                 filmtype *filmada,unsigned int k){

        size_t c = 0;
        size_t start = 0;
        std::vector<key_type> startkeys;
//        if (k>0){
//            cout << " Jesus, i need You!" << endl;
//        }
//        if (k>1){
//            cout << " Jesus, i need You!" << endl;
//        }


        for (size_t i = 0; i < keys.size(); ++i) {
            pair<key_type,int> p(keys[i],filmada->innerlevels[k]->nextpos++) ;   // i 为 pos
            ++(filmada->innerlevels[k]->pos);
            if (!filmada->innerlevels[k]->opt[0]->add_point(p.first, p.second)) {  // 如果inner level  不满足error 了，那么再创建一个innerpiece
//                if (k>0){
////                    cout << " Jesus, i need You!" << endl;
//                }
//                if (k>1){
////                    cout << " Jesus, i need You!" << endl;
//                }
                // 当前 innerpiece 不再满足，需要创建new inner piece 并判断该 innerlevel 的上一层level 是否需要更新
                auto a = filmada->innerlevels[k]->opt[0]->get_segment();

                if (filmada->innerlevels[k]->pos > 2)
                    filmada->innerlevels[k]->innerpieces.pop_back();
                filmada->innerlevels[k]->innerpieces.emplace_back(a);
                // 首先在该层创建一个 new innerpiece， 更新该innerpiece，再递归向向上
//                    typename filmtype::Innerpiece *innerpiece = new typename filmtype::;
                filmada->innerlevels[k]->pos = 0;
                filmada->innerlevels[k]->nextpos -= 2;
                delete filmada->innerlevels[k]->opt[0];
                filmada->innerlevels[k]->opt.pop_back();
                insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
                filmada->innerlevels[k]->opt.emplace_back(inneropt);
                if (k==0){
                    auto aaaaa = filmada->leaflevel.leafpieces.size()-2;
                    startkeys.emplace_back(filmada->leaflevel.leafpieces[filmada->leaflevel.leafpieces.size()-1]->startkey);
                }

                else{
                    auto aaaaa = filmada->innerlevels[k-1]->innerpieces.size()-2;
                    startkeys.emplace_back(filmada->innerlevels[k-1]->innerpieces[filmada->innerlevels[k-1]->innerpieces.size()-2].startkey);
                }
                startkeys.emplace_back(p.first);
                auto rr = append_segmentation(error,startkeys,filmada,k);

                ++c;
                if (filmada->innerlevels.back()->innerpieces.size() > 1)
                {
                    startkeys.clear();
                    startkeys.emplace_back(a.first);
                    startkeys.emplace_back(p.first);
                    typename filmtype:: Innerpiece innerpiece;// 创建parent piece
                    insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
//                    std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> > *innerlevel =
//                            new std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> >;
                    typename filmtype::Innerlevel *innerlevel = new typename filmtype::Innerlevel;
                    filmada->innerlevels.emplace_back(innerlevel);
                    filmada->innerlevels.back()->opt.push_back(inneropt) ;
                    auto rr = append_segmentation(error,startkeys,filmada,k+1);
//                    cout<< "Jesus, i need You !"  << endl;
                    startkeys.emplace_back(a.first);
                    return std::pair<size_t,std::vector<key_type> >(++c,startkeys) ;
                }
                else if (filmada->innerlevels.size() > 1 & k != filmada->innerlevels.size()-1 )   // 如上为 由于创建了new innner piece，导致了new innerlevel，如下为，虽然创建了new innerpiece，但只需要更新上层的innner level
                {
                    // 更新上层的最后一个inner piece
                    startkeys.pop_back();
                    auto rr = append_segmentation(error,startkeys,filmada,k+1);
                    startkeys.clear();
//                    cout << "thank You, my Lord! i need You!"<<endl;
                    startkeys.emplace_back(a.first);
                    return std::pair<size_t,std::vector<key_type> >(++c,startkeys) ;
                }
                else if (filmada->innerlevels.size() > 1 & k == filmada->innerlevels.size()-1 )   // 如上为 由于创建了new innner piece，导致了new innerlevel，如下为，虽然创建了new innerpiece，但只需要更新上层的innner level
                {
                    cout << "thank You, my Lord! i need You!"<<endl;
                    startkeys.clear();
                    return std::pair<size_t,std::vector<key_type> >(++c,startkeys) ;
                }
                else{
                    startkeys.clear();
                }
                a = filmada->innerlevels[k]->opt[0]->get_segment();
                if (filmada->innerlevels[k]->pos > 2)
                    filmada->innerlevels[k]->innerpieces.pop_back();
                filmada->innerlevels[k]->innerpieces.emplace_back(a);

                startkeys.emplace_back(a.first);
                return std::pair<size_t,std::vector<key_type> >(++c,startkeys) ;
            }

        }

        auto a = filmada->innerlevels[k]->opt[0]->get_segment();
        if (filmada->innerlevels[k]->pos > 2)
            filmada->innerlevels[k]->innerpieces.pop_back();
        filmada->innerlevels[k]->innerpieces.emplace_back(a);

        startkeys.emplace_back(a.first);
        return std::pair<size_t,std::vector<key_type> >(++c,startkeys) ;

    }



    template<typename key_type, typename leaf_type,typename filmtype,typename filmadalrutype>
    std::pair<size_t,std::vector<key_type> > append_segmentation(size_t error,std::vector<key_type> keys,std::vector<key_type> payload,
                                                                 filmtype *filmada, unsigned int &inkeynum,leaf_type* m_tailleaf,filmadalrutype *interchain){
        size_t c = 0;
        size_t start = 0;
        std::vector<key_type> startkeys;
        unsigned int pos = 0;
        pair<key_type,unsigned int> p(keys[0],pos++) ;

        inkeynum +=1;

        if (filmada->m_tailleaf == NULL)
        {
            leaf_type *cur_leaf = new leaf_type;
            cur_leaf->slotkey.reserve(8192*40);

            filmada->m_tailleaf = cur_leaf;
            m_tailleaf = cur_leaf;
            insertPWLF<key_type, int> *leafopt = new insertPWLF<key_type, int> (error);

            filmada->leaflevel.opt.push_back(leafopt) ;

        }
//        auto xx = filmada->leaflevel.second[0]->buffqueue.size();

        filmada->leaflevel.opt[0]->append_point(p.first, p.second);
        filmada->m_tailleaf->slotdata.emplace_back( filmada->m_tailleaf->intrachain.put(p.second,payload));
        filmada->m_tailleaf->slotkey.push_back(p.first);

        for (inkeynum; inkeynum < keys.size(); ++inkeynum) {
            pair<key_type,unsigned int> next_p(keys[inkeynum],pos++) ;
            if (inkeynum != start && next_p.first == p.first)
                continue;
            p = next_p;
            if (filmada ->leaflevel.opt[0]->append_point(p.first, p.second)) {
                filmada->m_tailleaf->slotkey.emplace_back(p.first);
                filmada->m_tailleaf->slotdata.emplace_back( filmada->m_tailleaf->intrachain.put(p.second,payload));
            }
            else
            {
                start = inkeynum;
//                if (start == 85741){
//                    cout<< " i need You, my Lord!"<< endl;
//                }
                auto a = filmada ->leaflevel.opt[0]->get_segment(keys[--inkeynum]); // 将生成的new leaf piece插入到leaflevel 中
                filmada->m_tailleaf->update(a);
                filmada ->leaflevel.leafpieces.emplace_back( filmada->m_tailleaf);
                interchain->put(filmada->m_tailleaf->startkey,filmada->m_tailleaf);
                // 这里是初始化 parent piece 的first key 和 last key
                if (filmada->innerlevels.size() == 0){
                    startkeys.emplace_back( filmada->m_tailleaf->startkey);
                    startkeys.emplace_back( p.first);
                    typename filmtype:: Innerpiece innerpiece;// 创建parent piece
                    insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
//                    std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> > *innerlevel =
//                            new std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> >;
                    typename filmtype::Innerlevel *innerlevel = new typename filmtype::Innerlevel;
                    filmada->innerlevels.emplace_back(innerlevel);
                    filmada->innerlevels[0]->opt.push_back(inneropt) ;
                    cout << " my Lord, Jesus, please have pity on me"<< endl;
                }
                else{
                    // 从innerlevel 的最底层到root 层，判断是否需要更新
                    startkeys.emplace_back( p.first);
//                    cout << "my lovely Lord, i trust in You!" << endl;
                }
                auto rr = append_segmentation(error,startkeys,filmada,0);
                startkeys.clear();
//                cout << "Jesus, i need You!!"<< endl;
                pos = 0;

                filmada->m_tailleaf = new leaf_type;
                filmada->m_tailleaf->slotkey.reserve(8192*40);
                ++c;
            }
        }

        auto a = filmada ->leaflevel.opt[0]->get_segment( keys.back());
        filmada->m_tailleaf->update(a);
        filmada->leaflevel.leafpieces.emplace_back( filmada->m_tailleaf);
        interchain->put(filmada->m_tailleaf->startkey,filmada->m_tailleaf);
        return std::pair<size_t,vector<key_type>> (++c,startkeys);

    }

    template<typename key_type,typename filmtype>
    std::pair<size_t,std::vector<key_type> > batch_segmentation(size_t error,std::vector<key_type> keys,
                                                                 filmtype *filmada,unsigned int k){

        size_t c = 0;
        size_t start = 0;
        std::vector<key_type> startkeys;
//        if (k>0){
//            cout << " Jesus, i need You!" << endl;
//        }
//        if (k>1){
//            cout << " Jesus, i need You!" << endl;
//        }


        for (size_t i = 0; i < keys.size(); ++i) {
            pair<key_type,int> p(keys[i],filmada->innerlevels[k]->nextpos++) ;   // i 为 pos
            ++(filmada->innerlevels[k]->pos);
            if (!filmada->innerlevels[k]->opt[0]->add_point(p.first, p.second)) {  // 如果inner level  不满足error 了，那么再创建一个innerpiece
//                if (k>0){
////                    cout << " Jesus, i need You!" << endl;
//                }
//                if (k>1){
////                    cout << " Jesus, i need You!" << endl;
//                }
                // 当前 innerpiece 不再满足，需要创建new inner piece 并判断该 innerlevel 的上一层level 是否需要更新
                auto a = filmada->innerlevels[k]->opt[0]->get_segment();

                if (filmada->innerlevels[k]->pos > 2)
                    filmada->innerlevels[k]->innerpieces.pop_back();
                filmada->innerlevels[k]->innerpieces.emplace_back(a);
                // 首先在该层创建一个 new innerpiece， 更新该innerpiece，再递归向向上
//                    typename filmtype::Innerpiece *innerpiece = new typename filmtype::;
                filmada->innerlevels[k]->pos = 0;
                filmada->innerlevels[k]->nextpos -= 2;
                delete filmada->innerlevels[k]->opt[0];
                filmada->innerlevels[k]->opt.pop_back();
                insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
                filmada->innerlevels[k]->opt.emplace_back(inneropt);
                if (k==0){
                    auto aaaaa = filmada->leaflevel.leafpieces.size()-2;
                    startkeys.emplace_back(filmada->leaflevel.leafpieces[filmada->leaflevel.leafpieces.size()-1]->startkey);
                }

                else{
                    auto aaaaa = filmada->innerlevels[k-1]->innerpieces.size()-2;
                    startkeys.emplace_back(filmada->innerlevels[k-1]->innerpieces[filmada->innerlevels[k-1]->innerpieces.size()-2].startkey);
                }
                startkeys.emplace_back(p.first);
                auto rr = append_segmentation(error,startkeys,filmada,k);

                ++c;
                if (filmada->innerlevels.back()->innerpieces.size() > 1)
                {
                    startkeys.clear();
                    startkeys.emplace_back(a.first);
                    startkeys.emplace_back(p.first);
                    typename filmtype:: Innerpiece innerpiece;// 创建parent piece
                    insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
//                    std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> > *innerlevel =
//                            new std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> >;
                    typename filmtype::Innerlevel *innerlevel = new typename filmtype::Innerlevel;
                    filmada->innerlevels.emplace_back(innerlevel);
                    filmada->innerlevels.back()->opt.push_back(inneropt) ;
                    auto rr = append_segmentation(error,startkeys,filmada,k+1);
//                    cout<< "Jesus, i need You !"  << endl;
                    startkeys.emplace_back(a.first);
                    return std::pair<size_t,std::vector<key_type> >(++c,startkeys) ;
                }
                else if (filmada->innerlevels.size() > 1 & k != filmada->innerlevels.size()-1 )   // 如上为 由于创建了new innner piece，导致了new innerlevel，如下为，虽然创建了new innerpiece，但只需要更新上层的innner level
                {
                    // 更新上层的最后一个inner piece
                    startkeys.pop_back();
                    auto rr = append_segmentation(error,startkeys,filmada,k+1);
                    startkeys.clear();
//                    cout << "thank You, my Lord! i need You!"<<endl;
                    startkeys.emplace_back(a.first);
                    return std::pair<size_t,std::vector<key_type> >(++c,startkeys) ;
                }
                else if (filmada->innerlevels.size() > 1 & k == filmada->innerlevels.size()-1 )   // 如上为 由于创建了new innner piece，导致了new innerlevel，如下为，虽然创建了new innerpiece，但只需要更新上层的innner level
                {
                    cout << "thank You, my Lord! i need You!"<<endl;
                    startkeys.clear();
                    return std::pair<size_t,std::vector<key_type> >(++c,startkeys) ;
                }
                else{
                    startkeys.clear();
                }
                a = filmada->innerlevels[k]->opt[0]->get_segment();
                if (filmada->innerlevels[k]->pos > 2)
                    filmada->innerlevels[k]->innerpieces.pop_back();
                filmada->innerlevels[k]->innerpieces.emplace_back(a);

                startkeys.emplace_back(a.first);
                return std::pair<size_t,std::vector<key_type> >(++c,startkeys) ;
            }

        }

        auto a = filmada->innerlevels[k]->opt[0]->get_segment();
        if (filmada->innerlevels[k]->pos > 2)
            filmada->innerlevels[k]->innerpieces.pop_back();
        filmada->innerlevels[k]->innerpieces.emplace_back(a);

        startkeys.emplace_back(a.first);
        return std::pair<size_t,std::vector<key_type> >(++c,startkeys) ;

    }

    template<typename key_type, typename leaf_type,typename filmtype,typename filmadalrutype>
    std::pair<size_t,std::vector<key_type> > batch_segmentation_use_p_struct(size_t error,std::vector<key_type> keys,std::vector<key_type> payload,
                                                                 filmtype *filmada, unsigned int &inkeynum,leaf_type* m_tailleaf,filmadalrutype *interchain){
        size_t c = 0;
        size_t start = 0;
        std::vector<key_type> startkeys;
        unsigned int pos = 0;
//        pair<key_type,unsigned int> p(keys[0],pos++) ;

//        inkeynum +=1;

        if (filmada->m_tailleaf == NULL)
        {
            leaf_type *cur_leaf = new leaf_type;
            cur_leaf->slotkey.reserve(8192*40);

            filmada->m_tailleaf = cur_leaf;
            m_tailleaf = cur_leaf;
            insertPWLF<key_type, int> *leafopt = new insertPWLF<key_type, int> (error);

            filmada->leaflevel.opt.push_back(leafopt) ;

        }
//        auto xx = filmada->leaflevel.second[0]->buffqueue.size();

//        filmada->leaflevel.opt[0]->append_point(p.first, p.second);
//        filmada->m_tailleaf->slotdata.emplace_back( filmada->m_tailleaf->intrachain.put(p.second,payload));
//        filmada->m_tailleaf->slotkey.push_back(p.first);
        unsigned cachen = 4;
        pair<key_type,unsigned int> cache[cachen] ;
        for (inkeynum; inkeynum < keys.size(); ) {
            // 需要在这里添加cache 的优化
            //pair<key_type,unsigned int> next_p(keys[inkeynum],pos++) ;
            for (int cachei = 0; cachei < cachen; cachei ++)
            {
                pair<key_type,unsigned int> next_p(keys[inkeynum+cachei],pos++) ;
                cache[cachei] = next_p;
            }

            //if (inkeynum != start && next_p.first == p.first)
                //continue;
            //if (inkeynum != start )
                //continue;
            for (int cachei = 0; cachei < cachen; cachei ++){
//                if (cache[cachei].first == 80919816)
//                    cout << "Jesus, i need You!!" <<  endl;
                if (filmada ->leaflevel.opt[0]->append_point(cache[cachei].first, cache[cachei].second)) {
                    filmada->m_tailleaf->slotkey.emplace_back(cache[cachei].first);
                    filmada->m_tailleaf->slotdata.emplace_back( filmada->m_tailleaf->intrachain.put(cache[cachei].second,payload));
                    inkeynum ++;
                }
                else
                { // 插入不成功，即需要create new piece
                    start = inkeynum;
//                if (start == 85741){
//                    cout<< " i need You, my Lord!"<< endl;
//                }
                    auto a = filmada ->leaflevel.opt[0]->get_segment(keys[inkeynum]); // 将生成的new leaf piece插入到leaflevel 中， 这个是last key
                    filmada->m_tailleaf->update(a);
                    filmada ->leaflevel.leafpieces.emplace_back( filmada->m_tailleaf);
                    interchain->put(filmada->m_tailleaf->startkey,filmada->m_tailleaf);
                    // 这里是初始化 parent piece 的first key 和 last key
                    if (filmada->innerlevels.size() == 0){
                        startkeys.emplace_back( filmada->m_tailleaf->startkey);
                        startkeys.emplace_back( cache[cachei].first);
                        typename filmtype:: Innerpiece innerpiece;// 创建parent piece
                        insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);
//                    std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> > *innerlevel =
//                            new std::pair< std::vector<typename filmtype:: Innerpiece>, std::vector<internal::insertPWLF<key_type, int>*> >;
                        typename filmtype::Innerlevel *innerlevel = new typename filmtype::Innerlevel;
                        filmada->innerlevels.emplace_back(innerlevel);
                        filmada->innerlevels[0]->opt.push_back(inneropt) ;
                        cout << " my Lord, Jesus, please have pity on me"<< endl;
                    }
                    else{
                        // 从innerlevel 的最底层到root 层，判断是否需要更新
                        startkeys.emplace_back( cache[cachei].first);
//                    cout << "my lovely Lord, i trust in You!" << endl;
                    }
                    auto rr = append_segmentation(error,startkeys,filmada,0);
                    startkeys.clear();
//                cout << "Jesus, i need You!!"<< endl;
                    pos = 0;

                    filmada->m_tailleaf = new leaf_type;
                    filmada->m_tailleaf->slotkey.reserve(8192*40);
                    ++c;
                    cachei = cachen;
                }
            }

        }
        if (inkeynum != start ){
            auto a = filmada ->leaflevel.opt[0]->get_segment( keys.back());
            filmada->m_tailleaf->update(a);
            filmada->leaflevel.leafpieces.emplace_back( filmada->m_tailleaf);
            interchain->put(filmada->m_tailleaf->startkey,filmada->m_tailleaf);
        }

        return std::pair<size_t,vector<key_type>> (++c,startkeys);

    }

    template<typename key_type, typename leaf_type,typename filmtype,typename filmadalrutype>
    std::pair<size_t,std::vector<key_type> > batch_segmentation(size_t error,std::vector<key_type> keys,std::vector<key_type> payload,
                                                                filmtype *filmada, unsigned int &inkeynum,leaf_type* m_tailleaf,filmadalrutype *interchain){
        size_t c = 0;
        size_t start = 0;
        std::vector<key_type> startkeys;
        unsigned int pos = 0;
//        pair<key_type,unsigned int> p(keys[0],pos++) ;

        if (filmada->m_tailleaf == NULL)
        {
            leaf_type *cur_leaf = new leaf_type;
            cur_leaf->slotkey.reserve(8192*20);

            filmada->m_tailleaf = cur_leaf;
            m_tailleaf = cur_leaf;
            insertPWLF<key_type, int> *leafopt = new insertPWLF<key_type, int> (error);

            filmada->leaflevel.opt.push_back(leafopt) ;

        }


//        filmada->leaflevel.opt[0]->append_point(keys[0], pos);
//        filmada->m_tailleaf->slotdata.emplace_back( filmada->m_tailleaf->intrachain.put(pos++,payload));
//        filmada->m_tailleaf->slotkey.push_back(keys[0]);
//        inkeynum +=1;
//        auto iter = filmada->m_tailleaf->slotkey.begin()+1;

        auto iter = filmada
                ->m_tailleaf->slotkey.begin();
        unsigned cachen = 256;
        unsigned int cachei = 0 ;
        key_type cachekey[cachen] ;
//        unsigned int cachepos[cachen];
        for (inkeynum; inkeynum < keys.size(); ) {
            // 需要在这里添加cache 的优化
            //pair<key_type,unsigned int> next_p(keys[inkeynum],pos++) ;
            for (cachei = 0; (inkeynum + cachei) <keys.size() and cachei < cachen; cachei ++)
            {
//                pair<key_type,unsigned int> next_p(keys[inkeynum+cachei],pos++) ;
                cachekey[cachei] = keys[inkeynum+cachei];
            }
            for (cachei = 0; cachei < cachen; cachei ++){
               // if (cachekey[cachei] == 1336030)
                    //cout << "Jesus, i need You!!" <<  endl;
                if (filmada ->leaflevel.opt[0]->append_point(cachekey[cachei], pos)) {
//                    filmada->m_tailleaf->slotkey.emplace_back(cachekey[cachei]);
                    filmada->m_tailleaf->slotdata.emplace_back( filmada->m_tailleaf->intrachain.put(pos++,payload));
                    inkeynum++;
                }
                else
                { // 插入不成功，即需要create new piece
                    start = inkeynum ;
                    filmada->m_tailleaf->slotkey.insert(iter,cachekey,cachekey+cachei);

                    auto a = filmada ->leaflevel.opt[0]->get_segment(keys[inkeynum-1]); // 将生成的new leaf piece插入到leaflevel 中， 这个是last key
                    filmada->m_tailleaf->update(a);
                    filmada ->leaflevel.leafpieces.emplace_back( filmada->m_tailleaf);
                    interchain->put(filmada->m_tailleaf->startkey,filmada->m_tailleaf);
                    // 这里是初始化 parent piece 的first key 和 last key
                    if (filmada->innerlevels.size() == 0){
                        startkeys.emplace_back( filmada->m_tailleaf->startkey);
                        startkeys.emplace_back( cachekey[cachei]);
                        typename filmtype:: Innerpiece innerpiece;// 创建parent piece
                        insertPWLF<key_type, int> *inneropt = new insertPWLF<key_type, int>(error);

                        typename filmtype::Innerlevel *innerlevel = new typename filmtype::Innerlevel;
                        filmada->innerlevels.emplace_back(innerlevel);
                        filmada->innerlevels[0]->opt.push_back(inneropt) ;
                        cout << " my Lord, Jesus, please have pity on me"<< endl;
                    }
                    else{
                        // 从innerlevel 的最底层到root 层，判断是否需要更新
                        startkeys.emplace_back( cachekey[cachei]);
//                    cout << "my lovely Lord, i trust in You!" << endl;
                    }
                    auto rr = append_segmentation(error,startkeys,filmada,0);
                    startkeys.clear();
//                cout << "Jesus, i need You!!"<< endl;
                    pos = 0;
                    filmada->m_tailleaf = new leaf_type;
                    filmada->m_tailleaf->slotkey.reserve(8192*20);
                    iter = filmada->m_tailleaf->slotkey.begin();
                    ++c;
                    cachei = cachen;
                }

            }
            if (pos !=0){
                //auto xxxxx = *iter;
                filmada->m_tailleaf->slotkey.insert(iter,cachekey,cachekey+cachei);
                iter = filmada->m_tailleaf->slotkey.end();
            }

        }
        if (inkeynum != start ){
            auto a = filmada ->leaflevel.opt[0]->get_segment( keys.back());
            filmada->m_tailleaf->update(a);
            filmada->leaflevel.leafpieces.emplace_back( filmada->m_tailleaf);
            interchain->put(filmada->m_tailleaf->startkey,filmada->m_tailleaf);
        }

        return std::pair<size_t,vector<key_type>> (++c,startkeys);

    }

}




#endif //FILMINSERT_PWLF_H
