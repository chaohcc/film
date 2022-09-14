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

the header file about datasets and workloads are data.h, zipf.h

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


the core source file of film are in film.h, filmadalru.h, filmadastorage, pwlf.h

        -- film.h: the core source file of film
                --- test_interleave_insert_query
                --- search_one
                --- search_range
                --- append_one

        -- filmadalru.h: the header file about adaptive LRU
                --- globalchain: hashLRU
                --- localchian: localLRU

        -- filmadastorage.h: the header file about data transfer, disk access, cold data eviction

        -- pwlf.h: header file about subrange partion, build learned model
                --- piece design
                --- append segmentation