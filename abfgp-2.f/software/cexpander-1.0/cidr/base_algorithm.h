// edouard severing 2010

#ifndef BASE_ALGORITHM_H_INCLUDED
#define BASE_ALGORITHM_H_INCLUDED

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctype.h>

using namespace std;

// simple function for line-start
bool line_st( string& line, string st );

// sequence string size item

struct sequence_entry{
    string header;
    unsigned int size_p1;
    string sequence;
};

// container for sequence entries

typedef vector<sequence_entry> sequence_entry_vector;


// object for alignment-reading

struct alignment_entry{
    string sequence1;
    string sequence2;
    string a1;
    string a2;
};

// checker function
bool alignment_entry_ok( alignment_entry& leEntry );

// empty alignment_entry

alignment_entry empty_alignment_entry();

// read next_alignment
alignment_entry get_next_alignment( ifstream& infile );

// link_vector
struct link_vector{
    vector< vector< int > > link_pair;
};

// initiator of link_vector
void initiate_link_vector( link_vector& l );


// link_vector container
typedef vector< link_vector > link_vector_container;


// converter function
void alignment_to_link_vector( string& a1, string& a2,
    link_vector& links, int l1, int l2 );

typedef vector<double> consistency_count_vector;

typedef vector<bool> consistency_vector;

typedef vector<int> pos_vector;

class alignment_base{
    public:
        // report file read

        int read_reportfile( const char* Fname );
        // clear sequence entries

        void clear_entries() { sequence_entries.clear(); }

        int header_index( string sname );

        void a_to_sequence( string& a, sequence_entry& leEntry );

        // alignment input (NEXT)
        int read_alignmentfile( const char* Fname );

        // full_read

        bool full_read( const char* rFile, const char* aFile );


        void set_count_gaps_as_true( bool value ) { count_gaps_as_true = value;}

        bool get_count_gaps_as_true() { return count_gaps_as_true; }

    //private:

        // sets_up the link_vector ( triangle )
        void setup_vectors();
        // calculate position of pair-link
        int link_pos( int index1, int index2 );
        void link_pos_focus( int& rPos, int& Focus_1, int index1, int index2 );

        bool is_tri_consistent( int a, int b, int c, int posOnA );
        sequence_entry_vector sequence_entries;

        pos_vector retrieve_positions( consistency_vector& leVector, int root_nr,  int subject_nr );

        string retrieve_sequence( consistency_vector& leVector, int root_nr,  int subject_nr );

        // consistency triangle_count
        float n1_triangles();
        consistency_count_vector full_consistency_count( int nr, bool dofrac );
        consistency_vector full_consistency( int nr );
        // contains the link vector
        link_vector_container link_vectors;
        bool count_gaps_as_true;
};


#endif // BASE_ALGORITHM_H_INCLUDED
