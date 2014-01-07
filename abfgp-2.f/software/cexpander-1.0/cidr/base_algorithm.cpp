#include "base_algorithm.h"

bool line_st( string& line, string st )
{
    if ( line.substr( 0, st.size() ) == st )
        return true;
    else
        return false;
}

bool alignment_entry_ok( alignment_entry& leEntry )
{
    if ( leEntry.sequence1.size() < 1 )
        return false;

    if ( leEntry.sequence2.size() < 1 )
        return false;

    if ( leEntry.a1.size() < 1 )
        return false;

    if ( leEntry.a2.size() < 1 )
        return false;

    // consistency

    if ( leEntry.a2.size() != leEntry.a1.size() )
        return false;

    // if we are here we passed the test

    return true;
}

// empty alignment_entry

alignment_entry empty_alignment_entry()
{
    alignment_entry tmp;
    tmp.a1 = "";
    tmp.a2 = "";
    tmp.sequence1 = "";
    tmp.sequence2 = "";

    return tmp;
}


// read next_alignment
alignment_entry get_next_alignment( ifstream& infile )
{
    // fetch pre-emptied struct
    alignment_entry tmp = empty_alignment_entry();

    string istring; // string to be read

    while ( ! infile.eof() )
    {
        getline( infile, istring );

        if ( line_st( istring, "Sequence1: " ) )
            tmp.sequence1 = istring.substr( 11 );

        else if ( line_st( istring, "Sequence2: " ) )
                tmp.sequence2 = istring.substr( 11 );
        else if ( line_st( istring, "a1:" ) )
                tmp.a1 = istring.substr( 3 );
        else if ( line_st( istring, "a2:" ) )
        {
            tmp.a2 = istring.substr( 3 );
            break; // here we hope that the alignent entry has been filled
        }
    }

    return tmp;
}

void initiate_link_vector( link_vector& l )
{
    l.link_pair.resize( 2 );
}
// converter function
void alignment_to_link_vector( string& a1, string& a2,
    link_vector& links, int l1, int l2 )
{
    // setup vectors by length
    initiate_link_vector( links );

    links.link_pair[0].resize( l1 );
    links.link_pair[1].resize( l2 );

    //
    int p1 = 0; // residue 1
    int p2 = 0; // residue 2


    for ( int i=0; i != a1.size(); i++ )
    {
        if ( isalpha( a1[i] ) && isalpha( a2[i] ) )
        {
            p1++;
            p2++;

            links.link_pair[ 0 ][ p1 ] = p2;
            links.link_pair[ 1 ][ p2 ] = p1;
        }
        else if ( isalpha( a1[i] ) && !isalpha( a2[i] ) )
        {
            p1++;
            links.link_pair[ 0 ][ p1 ] = 0;
        }
        else
        {
            p2++;
            links.link_pair[ 1 ][ p2 ] = 0;
        }
    }
}


//
// here start the methods of the alignment-base class
//
int alignment_base::read_reportfile( const char* fname )
{
    // do entrie clear
    clear_entries();

    ifstream myfile;
    string istring;
    //sstream mystringstream;

    myfile.open( fname );

    if ( !myfile.is_open() )
    {
        cerr << "can't open: " << fname << endl;
        return -1;
    }


    while (!myfile.eof())
    {
        // WARNING THERE IS YET NO ADDITIONAL ERROR CHECKING IN THIS SECTION
        // prepare new-entry
        sequence_entry new_entry;

        //  try to read header
        getline( myfile, new_entry.header, '\t' );

        if ( new_entry.header.size() != 0 )
        {
            // get next part
            getline( myfile, istring );

            if ( istring.size() != 0 )
            {
                // get size int and add to buffer
                new_entry.size_p1 = atoi( istring.c_str() ) + 1;
                // reset sequence
                new_entry.sequence = ""; // highly class specific_approach

                sequence_entries.push_back( new_entry );
            }
        }
    }

    myfile.close();

    return sequence_entries.size();

}


int alignment_base::header_index( string sname )
{
    if ( sequence_entries.size() < 1 )
        return -1;

    int r = sequence_entries.size() -1;
    int l = 0;

    if ( sname < sequence_entries[0].header )
        return -1;

    if ( sname > sequence_entries[r].header )
        return -1;

    int fnd = -1;

    while ( l <= r )
    {
        int mid = ( l + r ) / 2;

        if ( sname > sequence_entries[mid].header )
            l = mid + 1;
        else if ( sname < sequence_entries[ mid ].header )
            r = mid - 1;
        else
        {
            fnd = mid;
            break; // found
        }
    }
    return fnd;
}

int alignment_base::read_alignmentfile( const char* Fname )
{
    ifstream myfile( Fname );

    if ( !myfile.is_open() )
    {
        cerr << "can't open: " << Fname << endl;
        return -1;
    }

    int total = 0;
    while (!myfile.eof())
    {
        alignment_entry next_entry = get_next_alignment( myfile );

        // check alignment
        if ( !alignment_entry_ok( next_entry ) )
            break;

        total++;

        // here the conversion code comes in
        int F1 = header_index( next_entry.sequence1 );
        int F2 = header_index( next_entry.sequence2 );

        // here get sequence if neccessary

        a_to_sequence( next_entry.a1, sequence_entries[ F1 ] );
        a_to_sequence( next_entry.a2, sequence_entries[ F2 ] );

        int pos = link_pos( F1, F2 );

        if ( F1 < F2 )
            alignment_to_link_vector( next_entry.a1,
                                      next_entry.a2,
                                      link_vectors[ pos ],
                                      sequence_entries[ F1 ].size_p1,
                                      sequence_entries[ F2 ].size_p1
                                      );
        else
            alignment_to_link_vector( next_entry.a2,
                                      next_entry.a1,
                                      link_vectors[ pos ],
                                      sequence_entries[ F2 ].size_p1,
                                      sequence_entries[ F1 ].size_p1
                                      );


    }


    myfile.close();
    return total;
}

bool alignment_base::full_read( const char* rFile, const char* aFile )
{
    if ( read_reportfile( rFile ) < 1 )
        return false;

    setup_vectors();

    if ( read_alignmentfile( aFile ) < 1 )
        return false;

    // if we reach this place (something has been read)
    return true;
}

void alignment_base::setup_vectors()
{
     int N = sequence_entries.size();

     link_vectors.resize( ( ( N * N ) - N ) / 2 );
}

int alignment_base::link_pos( int index1, int index2 )
{
    int mn = min( index1, index2 );
    int mx = max( index1, index2 );

    return mn + ( ( mx * mx ) - mx ) / 2;
}

void alignment_base::link_pos_focus( int& rPos, int& Focus_1, int index1, int index2 )
{
    rPos = link_pos( index1, index2 );

    if ( index1 > index2 )
        Focus_1 = 1;
    else
        Focus_1 = 0;
}


bool alignment_base::is_tri_consistent( int a, int b, int c, int posOnA )
{
    // long but to prevent errors ( core algorithm will be splitted if neccesasary)
    int index_ab;
    int vab;

    link_pos_focus( index_ab, vab, a, b );

    int index_ac;
    int vac;

    link_pos_focus( index_ac, vac, a, c );

    int index_bc;
    int vbc;

    link_pos_focus( index_bc, vbc, b, c );

    int plink_ab = link_vectors[ index_ab ].link_pair[ vab ][ posOnA ];
    int plink_ac = link_vectors[ index_ac ].link_pair[ vac ][ posOnA ];

    if ( plink_ab == 0 || plink_ac == 0 )
        return false; // reverted no gaps-consideration for now; // gapped

    if ( link_vectors[ index_bc ].link_pair[ vbc ][ plink_ab ] == plink_ac )
        return true;
    else
        return false;
}

float alignment_base::n1_triangles()
{
    float l = sequence_entries.size() - 1;
    return ( ( l * l ) - l ) / 2;
}

consistency_count_vector alignment_base::full_consistency_count( int nr, bool dofrac )
{
    consistency_count_vector revect;
    revect.resize( sequence_entries[ nr ].size_p1 );
    revect[0] = 0;

    float d;

    if ( dofrac )
        d = n1_triangles();
    else
        d = 1;

    for ( int p = 1; p != sequence_entries[ nr ].size_p1; p++ )
    {
        revect[ p ] = 0; // initiate
        for ( int i = 0; i != sequence_entries.size(); i++ )
        {
            if ( i!= nr )
            {
                for ( int j = i + 1; j != sequence_entries.size(); j++ )
                {
                    if ( j!= nr )
                    {
                        if ( is_tri_consistent( nr, i, j, p ) )
                            revect[ p ] = revect[ p ] + 1;
                    }
                }
            }

        }

        revect[ p ] = revect[ p ] / d;
    }
    return revect;
}

consistency_vector alignment_base::full_consistency( int nr )
{
    consistency_vector revect;

    revect.resize( sequence_entries[ nr ].size_p1 );
    revect[0] = false; // initiate

    for ( int p = 1; p != sequence_entries[ nr ].size_p1; p++ )
    {
        revect[ p ] = true;
        for ( int i = 0; i != sequence_entries.size(); i++ )
        {
            if ( i != nr )
            {
                for ( int j = i + 1; j != sequence_entries.size(); j++ )
                {
                    if ( j != nr )
                    {

                        if ( ! is_tri_consistent( nr, i, j, p ) )
                        {
                            revect[ p ] = false;
                            //cout << "specific break" << endl;
                            break; // because it is not consistent
                        }
                    }
                }
            }
            if ( !revect[ p ] ) break; // break out of loop ( because false )
        }
    }
    return revect;
}

void alignment_base::a_to_sequence( string& a, sequence_entry& leEntry )
{
    if ( leEntry.sequence == "" )
    {
        // pre-allocate space
        leEntry.sequence.resize( leEntry.size_p1 - 1 );
        int p = -1;
        for ( int i = 0; i != a.size(); i++ )
        {
            if ( isalpha( a[i] ) )
            {
                p++;
                leEntry.sequence[ p ] = a[i];
            }
        }
    }
}

pos_vector alignment_base::retrieve_positions( consistency_vector& leVector, int root_nr, int subject_nr )
{
    pos_vector tmp;
    tmp.clear();

    for ( int i = 1; i != leVector.size(); i++ )
    {
        if ( leVector[ i ] == 1 )
        {
            if ( root_nr != subject_nr )
            {
                int k, focus;

                link_pos_focus( k, focus, root_nr, subject_nr );

                tmp.push_back( link_vectors[ k ].link_pair[ focus ][ i ] );
            }
            else tmp.push_back( i );
        }
    }
    return tmp;
}

string alignment_base::retrieve_sequence( consistency_vector& leVector, int root_nr, int subject_nr )
{
    string tmp = "";

    for ( int i = 1; i != leVector.size(); i++ )
    {
        if ( leVector[ i ] == 1 )
        {
            if ( root_nr != subject_nr )
            {
                int k, focus;

                link_pos_focus( k, focus, root_nr, subject_nr );

                tmp+= sequence_entries[ subject_nr ].sequence[ link_vectors[ k ].link_pair[ focus ][ i ] - 1 ];
            }
            else tmp+= sequence_entries[ root_nr ].sequence[ i - 1 ];
        }
    }
    return tmp;
}
