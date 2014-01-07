// edouard severing 2010

#ifndef DR_INTERACE_H_INCLUDED
#define DR_INTERACE_H_INCLUDED

#include <iostream>
#include <vector>
#include <string>
#include <list>
#include <algorithm>
#include "base_algorithm.h"

using namespace std;

struct run_item{
    int root_ind;
    consistency_count_vector count_vector;
    consistency_vector cons_vector;
    vector<string> sequence_block;
    vector< pos_vector > positions_block;
};

typedef vector< run_item > run_item_container;

typedef list< int > int_list;

struct run_query{
    int ind;
    string name;
};

typedef list< run_query > query_list;

class dr_inter{
    public:
        // constructor just resets everything
        dr_inter() { do_reset(); }
        // tries to load
        void try_loading( const char* rFileName, const char* aFileName );
        // check if files have been loaded
        bool are_files_loaded() { return fLoaded; }

        int max_root() { return base_object.sequence_entries.size() - 1; }


        bool can_run();

        void set_run_mode( int rM ) { run_mode = rM; }
        int get_run_mode() { return run_mode; }

        bool is_root_set() { return ( queries.size() > 0 ); }

        void clear_queries() { queries.clear(); }

        void add_query_index( int ind );
        void add_all_query();
        int number_of_queries() { return queries.size(); }

        void execute();

        void add_sequence_block( run_item& leItem );
        void add_position_block( run_item& leItem );

        string get_name( int ind ) { return base_object.sequence_entries[ ind ].header; }
        int get_lenght( int ind ) { return base_object.sequence_entries[ ind ].size_p1 - 1; }

        void start_result_iteration() { rIterator = ritems.begin(); }
        bool next_result( run_item& rItem );

        bool add_query_by_name( string leName );

    private:

        void do_reset();
        alignment_base base_object;

        string rFileName; // keeps the
        string aFileName;

        bool fLoaded;

        int current_root;

        int run_mode;

        int_list queries;
        run_item_container ritems;
        run_item_container::iterator rIterator;

};



#endif // DR_INTERACE_H_INCLUDED
