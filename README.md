# film
a fully learned index for larger-than-memory databases


#

### Getting started
FILM can be used as a header-only library
You will need to compile the program with at least the C++17 standard (e.g., '-std = c++17') in the [CMakeLists.txt]

In this repository, you can compile and run film with [main.cpp]

There are some examples in [main.cpp]
    -- test_interleave_insert_query: test the workload that interleave between inserts and queries,
    -- test_out_of_order_insertion: test out-of-order insertion, control by parameter out_of_order_frac
    -- test_query_baseline: test the query workload that with a fixed number of point and/or range queries



### Datasets and data type
You can also run this benchmark on your own dataset.
change the filepath in [data.h]
the format of your dataset need to be in either binary format or text format (one key per line).
the datasets we used in paper are in text format, binary format are also supported in [data.sh] with the 'load_binary_data'.

the currently support data type of keys are long int or double.
need to modify (typedef double key_type;) or (typedef long int key_type;)  in [film.h]

the header file about datasets is data.h

-- 'loaddata'
-- get the searched keys at runtime
    ----'get_search_keys_zipf'
    ----'get_search_keys'
    ----'get_search_keys_scrambledzipf'
    ----'get_search_keys_hotspot'

-- get the searched ranges at runtime
    ----'get_search_ranges_zipf'
    ----'get_search_ranges'
    ----'get_search_ranges_scrambledzipf'
    ----'get_search_ranges_hotspot'

-- 'loadpquery'  (pre-generated searched keys)
-- 'loadrquery' (pre-generated searched keys)

film.h:
    --- test_interleave_insert_query
    --- search_one
    --- search_range
    --- append_one

filmadalru.h: the header file about adaptive LRU
    --- globalchain: hashLRU
    --- localchian: localLRU

filmadastorage.h:
the header file about data transfer, disk access, cold data eviction




## the preformance of record size
the model size comparison of differnet methods in terms of record size.

<img src="https://user-images.githubusercontent.com/51820918/155710404-2abb0e9a-7a74-4718-b47c-8bceddbb463c.png" width="50%" height="50%">


![recordSizewiki_ts_add](https://user-images.githubusercontent.com/51820918/155705150-5a7aa409-503d-4ef0-9e06-ef00f2fc7db8.png)

## the range query performance with different amount of available memory 

<img src="https://user-images.githubusercontent.com/51820918/155710245-68bd16c0-8e0c-487d-9b74-51d34bd0871b.png" width="50%" height="50%">


## dataset
the books and wiki_ts are come from SOSD. ref: https://github.com/learnedsystems/SOSD

the optimal solution of generating piece-wise-linear functions has well studied by computional geometry [ref: Joseph O’Rourke. 1981. An on-line algorithm for fitting straight lines between data ranges. Commun. ACM 24, 9 (1981), 574–578.], and the PGM-Index has implemented it in C++ implementation[ref: https://github.com/gvinciguerra/PGM-index], the learned model of FILM is on the basis of the segmentation from PGM-index.

