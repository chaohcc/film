//
// Created by chaochao on 2021/12/23.
//

#ifndef EXPERIMENTCC12_FILMADALRU_H
#define EXPERIMENTCC12_FILMADALRU_H

#include <unordered_map>
#define MAX_INT    (((unsigned int)(-1))>>1)
// define the node in doubly-linked list
using namespace std;
namespace adalru {
    template<class Key, class Value>
    struct Node {
        Key key;
        Value value;
        Node *prev;
        Node *next;

        Node(Key k, Value v) : key(k), value(v), prev(nullptr), next(nullptr) {};

        Node() : key(), value(), prev(nullptr), next(nullptr) {};
//    Node(): prev(nullptr), next(nullptr){ };
    };

    template<class Key, class Value, class mapvalue>   //mapvalue 是指向在 lru 中的node 的指针  在hash map 中使用
    class hashLRU {

    public:
        int size = 0;
        int capacity;
        std::unordered_map<Key, mapvalue> map;
        Node<Key, Value> *head;
        Node<Key, Value> *tail;

        hashLRU(int def_capcity) {
            capacity = def_capcity;
            head = new Node<Key, Value>;
            tail = new Node<Key, Value>;
            head->next = tail;
            tail->prev = head;
        }

        hashLRU() {
            capacity = MAX_INT;
            head = new Node<Key, Value>;
            tail = new Node<Key, Value>;
            head->next = tail;
            tail->prev = head;
        }


        // get the k node in LRU and move it to the head of LRU
        Node<Key, Value> *get(Key k) {
            Node<Key, Value> *node = new Node<Key, Value>;
            if (map.find(k) != map.end()) {
                node = map[k];
                moveTohead(node);
                return node;
            }// k is existed in map
            else {
                return node;
            }

        }

        // 将node 插入到LRU 的head
        Node<Key, Value> *appendhead(Node<Key, Value> *node) {
            node->prev = head;
            node->next = head->next;
            head->next->prev = node;
            head->next = node;

        }

        // 将 找到的node 移动到head
        void moveTohead(Node<Key, Value> *node) {
            if (node->prev == head) return;
            node->prev->next = node->next;
            node->next->prev = node->prev;
            appendhead(node);
        }

        // put the k node into the head of LRU
        void put(Key k, Value v) {
            Node<Key, Value> *node = new Node<Key, Value>(k, v);
            if (map.find(k) == map.end()) // k is not existed in map
            {
                map[k] = node;
                // 判断 size, 如果size = capacity,说明LRU 满了
                if (size == capacity) {
                    poptail();
                }
                size += 1;
                appendhead(node);
            } else {
//                map[k]->value = v;
                moveTohead(map[k]);
            }
        }

        //remove the k node from LRU
        void remove(Key k) {
            // 首先找到k 所属的node
            if (map.find(k) == map.end()) return; // 说明要删除的node 不存在
            Node<Key, Value> *node = new Node<Key, Value>;
            node = map[k];
            node->prev->next = node->next;
            node->next->prev = node->prev;
            map.erase(k);
            size -= 1;
        }

        //remove the k node from LRU
        void removenode(Node<Key, Value> *node) {
            // 首先找到k 所属的node
            node->prev->next = node->next;
            node->next->prev = node->prev;
            map.erase(node->key);
            size -= 1;
        }

        // pop the tail of the LRU, that the least recent used item
        Node<Key, Value> *poptail() {
//            map.erase(tail->prev->key);
            map.erase(tail->prev->key);
            Node<Key, Value> *node;
            node = tail->prev;
            tail->prev->prev->next = tail;
            tail->prev = tail->prev->prev;
            size -= 1;
            return node;
        }

        //get the tail node from local LRU that from this leaf evict key
        Value get_tail() {
            auto tailnode = tail->prev;
            if (tailnode->value->intrachain.size==1)
                removenode(tailnode);
            return tailnode->value;
        }

    };

    template<class Key, class Value>
    class localLRU {
    public:
        int size = 0;

        Node<Key, Value> *head;
        Node<Key, Value> *tail;

        localLRU() {
            head = new Node<Key, Value>;
            tail = new Node<Key, Value>;
            head->next = tail;
            tail->prev = head;
        }

        // 将node 插入到LRU 的head
        Node<Key, Value> *appendhead(Node<Key, Value> *node) {
            node->prev = head;
            node->next = head->next;
            head->next->prev = node;
            head->next = node;
        }

        // 将 找到的node 移动到head
        void moveTohead(Node<Key, Value> *node) {
            if (node->prev == head) return;
            node->prev->next = node->next;
            node->next->prev = node->prev;
            appendhead(node);
        }

        // put the k node into the head of LRU
        pair<bool,Node<Key, Value>*>  put(Key k, Value v) {
            Node<Key, Value> *node = new Node<Key, Value>(k, v);
            size += 1;
            appendhead(node);
//            if (k > 4294946590)
//                cout<< "my Lord ,i need You!"<< endl;
            return pair<bool,Node<Key, Value>*> (true,node);
        }

//        pair<bool,Node<Key, Value>*>  put(Value data) {
//
//            Node<Key, Value> *node = new Node<Key, Value>(data[0], data);
//            size += 1;
//            appendhead(node);
//            return pair<bool,Node<Key, Value>*> (true,node);
//        }


        //remove the k node from LRU
       void remove_node(Node<Key, Value> *node) {
            node->prev->next = node->next;
            node->next->prev = node->prev;
            size -= 1;
        }



        // pop the tail of the local LRU, that the least recent used item
        Node<Key, Value> *poptail() {
            Node<Key, Value> *node;
            node = tail->prev;
            tail->prev->prev->next = tail;
            tail->prev = tail->prev->prev;
            size -= 1;
            return node;
        }

    };

}
#endif //EXPERIMENTCC12_FILMADALRU_H
