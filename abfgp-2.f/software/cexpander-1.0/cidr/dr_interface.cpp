#ifndef DR_INTERFACE_CPP_INCLUDED
#define DR_INTERFACE_CPP_INCLUDED

#include "dr_interace.h"

void dr_inter::try_loading( const char* rFile, const char* aFile )
{
     // full reset
     do_reset();

    if ( base_object.full_read( rFile, aFile ) )
    {
        rFileName = rFile;
        aFileName = aFile;
        fLoaded = true; // just indicate that file has been loaded
    }
    else fLoaded = false; // just indicate thate file has been loaded

}



bool dr_inter::can_run()
{
    if ( are_files_loaded()  && is_root_set() ) return true;

    return false;
}

void dr_inter::add_query_index( int ind )
{
    if ( ind == -1 )
        add_all_query();
    else
        // add if not in-list
        if ( find( queries.begin(), queries.end(), ind ) ==  queries.end() )
            queries.push_back( ind );
}

void dr_inter::add_all_query()
{
    // clear first
    clear_queries();

    // do all
    for ( int i = 0; i!= base_object.sequence_entries.size(); i++ )
        queries.push_back( i );

}

void dr_inter::execute()
{
    // clear run_items
    ritems.clear();

    // calcs conistency_vector count
    // calcs conistency_vector
    // generate sequence block
    // generate positions block
    int_list::iterator it;

    for (it = queries.begin(); it != queries.end(); it++ )
    {
        // make new object
        run_item tmp;
        tmp.root_ind = *it;
        cerr << "running "  << *it << endl;
        // calc consistency_vector
        tmp.cons_vector = base_object.full_consistency( *it );
        // calc consistency vector count (hard set fraction true)
        tmp.count_vector = base_object.full_consistency_count( *it, true );

        add_sequence_block( tmp );
        add_position_block( tmp );

        ritems.push_back( tmp );

    }
}

void dr_inter::add_sequence_block( run_item& leItem )
{
    leItem.sequence_block.clear();

    for (int i=0; i != base_object.sequence_entries.size(); i++ )
        leItem.sequence_block.push_back(
            base_object.retrieve_sequence( leItem.cons_vector,
                                           leItem.root_ind,
                                           i )
                                        );
}

void dr_inter::add_position_block( run_item& leItem )
{
    leItem.positions_block.clear();

    for (int i=0; i != base_object.sequence_entries.size(); i++ )
        leItem.positions_block.push_back(
            base_object.retrieve_positions( leItem.cons_vector,
                                            leItem.root_ind,
                                            i )
                                        );
}

bool dr_inter::next_result( run_item& rItem )
{
    if ( rIterator == ritems.end() )
        return false;

    rItem = *rIterator;

    rIterator++;
    return true;
}

bool dr_inter::add_query_by_name( string leName )
{
    int tmp = base_object.header_index( leName );

    if ( tmp < 0 )
        return false;

    add_query_index( tmp );
    return true;
}

void dr_inter::do_reset()
{
    // clear names
    rFileName = "";
    aFileName = "";
    // no file loaded
    fLoaded = false;
    // clear the queries
    clear_queries();
    // clear ritems
    ritems.clear();
}

#endif // DR_INTERFACE_CPP_INCLUDED
