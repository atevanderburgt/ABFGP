#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include "dr_interace.h"

using namespace std;

void main_title()
{
    string titlename = "Cexpander49-STC-dr (SaTan's Cave)\nEdouard Severing\n2010";

    cout << titlename << endl;
    cout << endl;
}

void handle_question( string question )
{
    cerr << question << endl << "-> ";
}

template < class T >
bool input_in_range_cc( T leftB, T rightB, T& rvalue )
{
    //string tmp;

    //getline( cin, tmp );
    string tmp;
    getline( cin, tmp );
    stringstream s(tmp);

    if ( !( s >> rvalue ) )
    {
        cerr << "(error) value convertion error !!!"  << endl;
        return false;
    }


    else if ( rvalue < leftB || rvalue > rightB)
    {
        cerr << "(error) value not in range !!!" << endl;
        return false;
    }
    else return true;
}

int first_vis( int pos, int end_pos, int i_pos, string& test_string )
{

    while ( pos != end_pos )
    {
        if ( test_string[ pos ] > 32 )
            return pos;

        pos = pos + i_pos;
    }

    return -1; // nothing
}


string strip_string( string test_string )
{
    int lp = first_vis( 0, test_string.size(), 1, test_string );
    int rp = first_vis( test_string.size() - 1, -1, -1, test_string );

    if ( lp >= rp )
        return "";

    return test_string.substr( lp, 1 + rp - lp );
}

void cexp_load_files( dr_inter& leInterface )
{
    string f1;
    string f2;
    handle_question( "Name of report file:" );
    getline( cin, f1 );
    handle_question( "Name of alignment file:" );
    getline( cin, f2 );

    // try loading
    leInterface.try_loading( f1.c_str(), f2.c_str() );

    if (!leInterface.are_files_loaded())
        cerr << "(error) files were not loaded !" << endl;

    else cerr << "files loaded:" << endl;

}

/* DEPRICATED
void cexp_set_run_mode( dr_inter& leInterface )
{
    int tmp;
    handle_question( "enter run mode (1 consitency only) (2 triangle fractions only) ( 3 both)" );

    if ( input_in_range_cc( 1, 3, tmp ) )
    {
        leInterface.set_run_mode( tmp );
        cerr << "run mode set:" << endl;
    }
}
*/

void cexp_add_query_index( dr_inter& leInterface )
{
    if ( !leInterface.are_files_loaded() )
    {
        cerr << "(error) files are not loaded" << endl;
        return; //exit from function
    }

    int tmp;
    handle_question( "enter query index-nr (-1 adds all)" );

    if (  input_in_range_cc( -1, leInterface.max_root(), tmp )  )
    {
        //
        leInterface.add_query_index( tmp );
        cerr << "new query list size:" << leInterface.number_of_queries() << endl;
    }
}

void cexp_add_query_name( dr_inter& leInterface )
{
    if ( !leInterface.are_files_loaded() )
    {
        cerr << "(error) files are not loaded" << endl;
        return;
    }

    string leName;

    handle_question( "enter query name" );

    getline( cin, leName );

    leName = strip_string( leName );

    if ( leName.size() > 0 )
    {
        if (!leInterface.add_query_by_name( leName ) )
            cerr << "(error) not found: " << leName << endl;
        else
            cerr << "new query list size:" << leInterface.number_of_queries() << endl;
    }
}


void cexp_execute( dr_inter& leInterface )
{
    if (!leInterface.can_run() )
    {
        cerr << "(error) not ready to run !!!" << endl;
        return;
    }

    cerr << "launching ...."  << endl;
    leInterface.execute();
    cerr << "Done "  << endl;

}
bool  bi_comp( string s1, string s2 )
{
    if ( s1 == s2 || s1 == ( "$" + s2 )  )
        return true;
    else
        return false;
}

template < class T >
    void cout_vector_linear( vector< T > leVector )
    {
        for( int i = 1; i != leVector.size(); i++ )
            cout << i << "\t" << leVector[i] << endl;
    }

template < class T >
    void cout_vector_block( vector< T > leVector )
    {
        for ( int i = 1; i != leVector.size(); i++ )
        {
            if ( ( i % 4 ) == 1 )
                cout << i << "\t";

            cout << leVector[i];

            if ( ( i % 4 ) == 0 )
                cout << endl;
            else
                cout << "\t";

        }

        if ( ( ( leVector.size() - 1 ) % 4 ) != 0 )
            cout << endl;
    }


void dump_consistency_vector_count( dr_inter& leInterface, bool linear )
{
    leInterface.start_result_iteration();

    run_item r;
    //cerr << "k" << endl;

    while ( leInterface.next_result( r ) )
    {
        cout << ">CVC_" << r.root_ind << endl;

        if ( linear )
            cout_vector_linear( r.count_vector );
        else
            cout_vector_block( r.count_vector );

        cout << endl;
    }
    cout << endl;
}

void dump_consistency_vector( dr_inter& leInterface, bool linear )
{
    leInterface.start_result_iteration();

    run_item r;

    while ( leInterface.next_result( r ) )
    {
        cout << ">CV_" << r.root_ind << endl;

        if ( linear )
            cout_vector_linear( r.cons_vector );
        else
            cout_vector_block( r.cons_vector );

        cout << endl;
    }
    cout << endl;
}

void dump_sequence_block( dr_inter& leInterface )
{
    leInterface.start_result_iteration();

    run_item r;

    while ( leInterface.next_result( r ) )
    {
        for ( int i = 0; i != r.sequence_block.size(); i++ )
        {
            cout << ">SEQ_" << r.root_ind << "_"  << i << endl;
            cout << r.sequence_block[ i ] << endl << endl;
        }
    }

    cout << endl;
}



void dump_names( dr_inter& leInterface, bool to_stdout )
{
    for ( int i = 0; i <= leInterface.max_root(); i++ )
    {
        if ( to_stdout )
            cout << "#\t" << i << "\t" << leInterface.get_name( i ) << endl;
        else
            cerr << "#\t" << i << "\t" << leInterface.get_name( i ) << endl;
    }

    if ( to_stdout )
        cout << endl;

    else
        cerr << endl;
}

void dump_position_block( dr_inter& leInterface, bool linear )
{
    leInterface.start_result_iteration();

    run_item r;

    while ( leInterface.next_result( r ) )
    {
        for ( int i = 0; i != r.positions_block.size(); i++ )
        {
            cout << ">PB_" << r.root_ind << "_" << i << endl;

            if ( linear )
                cout_vector_linear( r.positions_block[i] );
            else
                cout_vector_block( r.positions_block[i] );
            cout << endl;
        }

    }
    cerr << endl;
}

void dump_all( dr_inter& leInterface, bool v_lin, bool c_lin, bool p_lin )
{
    cout << "$names_start" << endl;
    dump_names( leInterface, true );
    cout << "$names_end" << endl << endl;

    cout << "$cvc_start" << endl;
    dump_consistency_vector_count( leInterface, c_lin );
    cout << "$cvc_end" << endl << endl;

    cout << "$cv_start" << endl << endl;
    dump_consistency_vector( leInterface, v_lin );
    cout << "$cv_end" << endl << endl;

    cout << "$PB_start" << endl << endl;
    dump_position_block( leInterface, p_lin );
    cout << "$PB_end" << endl << endl;

    cout << "$seq_start" << endl << endl;
    dump_sequence_block( leInterface );
    cout << "$seq_end" << endl;
}



int main()
{
    // main object_interface
    dr_inter my_interface;

    string in_put = "#";

    main_title();

    bool vlinear = false;
    bool clinear = false;
    bool plinear = false;

    while (true) // endless loop
    {
        //cerr << "<<" << in_put[0] << ">>" << endl;

        //if ( in_put.size() != 0 )
        //{
        cerr << "cex> ";
            //in_put = "";
        //}
        getline( cin, in_put );


        // strip

        in_put = strip_string( in_put );

        // to upper

        for (int i = 0; i != in_put.size(); i++ )
            in_put[ i ] = toupper( in_put[i] );



        // only exit command
        if ( in_put.size() > 0 )
        {
            cerr << endl;
            if ( bi_comp( in_put, "EXIT" ) )
                break;
            // file loading
            else if ( bi_comp( in_put, "LOAD" ) )
                cexp_load_files( my_interface );

            // set root


            else if ( bi_comp( in_put, "CLEARQUERIES" ) )
            {
                my_interface.clear_queries();
                cerr << "queries cleared" << endl;
            }
            else if ( bi_comp( in_put, "ADDQUERY"  ) )
            {
                cexp_add_query_index( my_interface );
                //cin.ignore( MAX_INT, '\n' ) // dirty way for flushing stream;
            }
            else if ( bi_comp( in_put, "ADDQUERYNAME" ) )
            {
                cexp_add_query_name( my_interface );
            }
            else if ( bi_comp( in_put, "RUN"  ) )
                cexp_execute( my_interface );

            else if ( bi_comp( in_put, "CVC_LINEAR" ) )
            {
                vlinear = true;
                cerr << "linear-count vector" << endl;
            }
            else if ( bi_comp( in_put, "CVC_BLOCK" ) )
            {
                vlinear = false;
                cerr << "block-count vector" << endl;
            }
            else if ( bi_comp( in_put, "DUMPCVC" ) )
            {
                dump_consistency_vector_count( my_interface, vlinear );
                cerr << "cvc-dumped" << endl;
            }

            else if ( bi_comp( in_put, "CV_LINEAR" ) )
            {
                clinear = true;
                cerr << "linear consistency-vector" << endl;
            }
            else if ( bi_comp( in_put, "CV_BLOCK" ) )
            {
                clinear = false;
                cerr << "block consistency-vector" << endl;
            }

            else if ( bi_comp( in_put, "DUMPCV" ) )
            {
                dump_consistency_vector( my_interface, clinear );
                cerr << "cv-dumped" << endl;
            }
            else if ( bi_comp( in_put, "DUMPSB" ) )
            {
                dump_sequence_block( my_interface );
                cerr << "sb-dumped" << endl;
            }
            else if ( bi_comp( in_put, "ENTRIES" ) )
            {
                dump_names( my_interface, false );
                cerr << "names-dumped" << endl;
            }
            else if ( bi_comp( in_put, "DUMPENTRIES" ) )
                dump_names( my_interface, true );
            else if ( bi_comp( in_put, "DUMPALL" ) )
            {
                dump_all( my_interface, vlinear, clinear, plinear );
                cerr << "all-dumped" << endl;
            }
            else if ( bi_comp( in_put, "DUMPPB" ) )
            {
                dump_position_block( my_interface, plinear );
                cerr << "pb-dumped" << endl;
            }
            else if ( bi_comp( in_put, "PB_BLOCK" ) )
            {
                plinear = false;
                cerr << "position block block" << endl;
            }
            else if ( bi_comp( in_put, "PB_LINEAR" ) )
            {
                plinear = true;
                cerr << "position block linear" << endl;
            }
            else
                cerr << "(error) Unknown command: " << in_put << endl;
        }
        else cerr << endl;

    }
    cerr << endl;
    return 0;

}
