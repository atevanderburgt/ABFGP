/***************************************************************************
 *   Copyright (C) 2006 by edouard severing   *
 *   edouard@localhost   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <cstdlib>
#include "alignucs.h"
#include "ebasics.h"
#include "alignprots.h"
#include "defragM.h"

using namespace std;

bool nucleotides = false;

bool localOnly = false;
bool chargeOverHang = false;

std::string inputfile = "";
std::string pMatrix = "blosum62";

int nt_gapOpen = 12;
int nt_gapExtention = 2;

int p_gapOpen = 10;
int p_gapExtention = 1;

int nt_match = 5;
int nt_missmatch = -4;
int nt_illchar = -4; 
std::string nt_DNAmatrix = "";

int ColW = 50;

bool showskip = true;

generalSequence Query1;
generalSequence Query2;

bool calcRoutes = false;
bool plotRoutes = false;
int maxRoutes = 100000;

// new 13-10-2006

bool computerOutput = false; 
std::string aTitle = "celL blaU alignment";

alignnucleotides AN;
alignproteins AP;

std::string curateProt(std::string S)
{
	for ( int j = 0; j != S.length(); j++ )
	{
		if ( S[j] < 65 || S[j] > 90 )
			S[j] = 'X';
	}
	return S;
}

void display_help()
{
	cout << "USAGE cbalignp parameters\n";
	cout << "\n";
	cout << "(Edouard Severing 2006 )\n";
	cout << "\n";
	cout << "Available options\n";
	cout << "\n";
	cout << "-i [string] input file required\n"; 
	cout << "-l          local alignment [sw] (def off )\n";
	cout << "-o [int]    affine gap open (+) (def n:12; p: 10 )\n";
	cout << "-e [int]    affine gap extention (+) (def n: 2; p: 1)\n";
	cout << "-O          full global (def off )\n";
	cout << "-s          score nucleotide match (def 5 )\n";
	cout << "-x          score nucleotide missmatch ( def -4 )\n";
	cout << "-I          score illchar ( def -4 )\n";
	cout << "-n          align nucleotides (def protein )\n";
	cout << "-D          load DNA matrix \n";
	cout << "-m          protein matrix (def blosum62)\n";
	cout << "-C          column width alignment (def 50)\n"; 
	cout << "-S          Don't Show skip ( def false )\n";
	cout << "-R          Calc routes ( def false ) \n";
	cout << "-P          Plot routes ( def false ) \n";
	cout << "-N          Max number of routes to extend ( def 100000 ) \n";
	cout << "-y          Computer Output ( not human friendly ) ( def off ) \n";
	cout << "-t          alignment title ( default celL blaU alignment )\n";
	cout << "\n";
	cout << "To save to alignment to a file use > fileout at the end of the\n";
	cout << "command line\n";
	cout << "\n"; 
}

void doAlignment()
{
	
	if ( nucleotides )
	{
		// new title 13-10-2006
		AN.setTitle( aTitle );
		
		AN.setmodus( localOnly, chargeOverHang,showskip );
		if ( nt_DNAmatrix == "" )
			AN.setScores( nt_match, nt_missmatch, nt_illchar );
		else
			AN.loadDNAMatrix( nt_DNAmatrix );
			
		AN.setPenalty( nt_gapOpen, nt_gapExtention ); 
		AN.setSequences( Query1, Query2 );
		AN.doalignment();
		
		// block added for COMPUTEROUTPUT 13-10-2006
		if (! computerOutput )
			AN.dumpAlignment(ColW);
		else
			AN.dumpAlignment_COMPUTERS();
		
		
		if ( calcRoutes )
		{
			TmatrixPaths Mpaths;
			TtraceMatrix tracematrix;
			
			traceM_clear( tracematrix ); // clear 
			traceM_build( tracematrix, AN.nCols(), AN.nRows() );
			traceM_init( tracematrix );
			traceM_construct( Mpaths, AN, maxRoutes );
			
			cout << "number of routes: " << Mpaths.size() << endl; 
			if ( Mpaths.size() > 0 )
				cout << Mpaths[0].path << endl;
			cout << AN.getT( AN.nCols() - 1, AN.nRows() - 1) << endl;
			
			if ( plotRoutes )
			{
				traceM_fill( tracematrix, Mpaths );
				traceM_Draw( tracematrix, AN );
			}
		}
	
	}
	
	else
	{
		// new title 13-10-2006
		AP.setTitle( aTitle );
		AP.setmodus( localOnly, chargeOverHang, showskip );
		AP.setPenalty( p_gapOpen, p_gapExtention );
		AP.loadInHouseMatrix( pMatrix );
		Query1.iSequence = curateProt( Query1.iSequence );
		Query2.iSequence = curateProt( Query2.iSequence ); 
		
		AP.setSequences( Query1, Query2 );
		AP.doalignment();
		
		// block changed 13-10-2006 for computer DUMPING 
		if ( ! computerOutput )
			AP.dumpAlignment(ColW);
		else
			AP.dumpAlignment_COMPUTERS();
		
		if ( calcRoutes )
		{
			TmatrixPaths Mpaths;
			TtraceMatrix tracematrix;
			
			traceM_clear( tracematrix ); // clear 
			traceM_build( tracematrix, AP.nCols(), AP.nRows() );
			traceM_init( tracematrix );
			traceM_construct( Mpaths, AP, maxRoutes );
			
			cout << "number of routes: " << Mpaths.size() << endl; 
			if ( Mpaths.size() > 0 )
				cout << Mpaths[0].path << endl;
			cout << AP.getT( AP.nCols() - 1, AP.nRows() - 1) << endl;
			
			if ( plotRoutes )
			{
				traceM_fill( tracematrix, Mpaths );
				traceM_Draw( tracematrix, AP );
			}
		}
	}
}
	

int main(int argc, char *argv[])
{

	if ( argc < 2 )
	{
		display_help();
		return EXIT_SUCCESS;
	}
	
	clineItemList commands;
  
	commands.push_back( clI( "-i","",0,false,1) ); // inputfile
	commands.push_back( clI( "-n","",0,false,0) ); // nucleotide alignment ? 
	
	commands.push_back( clI( "-o","",0,false,2) ); // gap open
	commands.push_back( clI( "-e","",0,false,2) ); // gap extention
	
	commands.push_back( clI( "-s","",0,false,2) ); // match -score ( nucleotides )
	commands.push_back( clI( "-x","",0,false,2) ); // miss-match score ( nucleotides )
	
	commands.push_back( clI( "-m","",0,false,1) ); // matrix loading 
	
	commands.push_back( clI( "-l","",0,false,0) ); // local only
	commands.push_back( clI( "-O","",0,false,0) ); // charge overhang
	
	commands.push_back( clI( "-C","",0,false,2) ); // column width 
	commands.push_back( clI( "-S","",0,false,0) );
	
	commands.push_back( clI( "-R","",0,false,0) );
	commands.push_back( clI( "-P","",0,false,0) );
	commands.push_back( clI( "-N","",0,false,2) );
	commands.push_back( clI( "-I","",0,false,2) ); 
	commands.push_back( clI( "-D","",0,false,1) ); 
	commands.push_back( clI( "-y","",0,false,0) );
	commands.push_back( clI( "-t","",0,false,1) );
	
	
	std::string cRes = handleCommandLine( commands, argc, argv );
  
	if ( cRes != "" )
	{
		cout << "err!!! " << cRes << endl;
		return EXIT_SUCCESS; // program has no error
	}
	
	// new 13-10-2006 -y and -t 
	
	if ( commands[ gItemNr( commands, "-t" ) ].switched )
		aTitle = commands[ gItemNr( commands, "-t" ) ].Istring;
	
	if ( commands[ gItemNr( commands, "-y" ) ].switched )
		computerOutput = true;
	if ( commands[ gItemNr(commands, "-S" ) ].switched )
		showskip = false;
	// parameter fill in 
	if ( commands[ gItemNr(commands, "-i" ) ].switched )
		inputfile = commands[ gItemNr( commands, "-i" ) ].Istring;
	
	nucleotides = commands[ gItemNr( commands, "-n" ) ].switched; // we want nucleotide alignment ? 
	
	if ( nucleotides )
	{
		if ( commands[ gItemNr( commands, "-o" ) ].switched )
			nt_gapOpen = commands[ gItemNr( commands, "-o" ) ].Ivalue;
		
		if ( commands[ gItemNr( commands, "-e" ) ].switched )
			nt_gapExtention = commands[ gItemNr( commands, "-e" )].Ivalue;
		
		if ( commands[ gItemNr( commands, "-s" ) ].switched )
			nt_match = commands[ gItemNr( commands,"-s" )].Ivalue; 
			
		if ( commands[ gItemNr( commands, "-x" ) ].switched )
			nt_missmatch = commands[ gItemNr( commands, "-x" ) ].Ivalue;
		
		if ( commands[ gItemNr( commands, "-I" ) ].switched )
			nt_illchar = commands[ gItemNr( commands, "-I" ) ].Ivalue;
			
		if ( commands[ gItemNr( commands, "-D" ) ].switched )
			nt_DNAmatrix = commands[gItemNr( commands, "-D" ) ].Istring; 
	}
	else
	{
		if ( commands[ gItemNr( commands, "-o" ) ].switched )
			p_gapOpen = commands[ gItemNr( commands, "-o" ) ].Ivalue;
		
		if ( commands[ gItemNr( commands, "-e" ) ].switched )
			p_gapExtention = commands[ gItemNr( commands, "-e" )].Ivalue;
	
		if ( commands[ gItemNr( commands, "-m" ) ].switched )
			pMatrix = commands[gItemNr( commands, "-m" ) ].Istring;
	}	
	
	// local and overhang charge
	
	chargeOverHang = commands[ gItemNr( commands, "-O" ) ].switched;
	localOnly = commands[gItemNr( commands, "-l" ) ].switched;
	
	if ( commands[gItemNr( commands, "-C" ) ].switched )
		ColW = commands[gItemNr( commands, "-C" ) ].Ivalue; 
	
	calcRoutes = commands[gItemNr( commands, "-R" )].switched;
	plotRoutes = commands[gItemNr( commands, "-P" )].switched;
	
	if ( commands[gItemNr( commands, "-N" ) ].switched )
		maxRoutes = commands[gItemNr( commands, "-N" ) ].Ivalue;
	
	char *Fname = new char[ inputfile.length() + 1 ];
	inputfile.copy(Fname, inputfile.length() );
	Fname[ inputfile.length() ] = 0;
	
	/*vector< generalSequence> S = readFastaFile( Fname );
	
	delete [] Fname;
	 
	if ( S.size() < 2 )
	{
		if ( S.size() == 0 )
			cout << "err!!! : no sequences loaded " << endl;
		else
			cout << "err!!! : not enough sequences " << endl;
			
		return EXIT_SUCCESS;
	}
	
	Query1 = S[0];
	Query2 = S[1]; 
	
	doAlignment();
	*/
	
	ifstream myFile;
	myFile.open( Fname );
	int launch = 0;
	while (!myFile.eof())
	{
	  Query1 = nextFasta( myFile );
	  Query2 = nextFasta( myFile ); 
	  
	  if (Query1.iHeader != "" && Query2.iHeader != "" )
	  {
	    if (Query2.iHeader != "" && Query2.iHeader != "" ) 
	    {
	      launch+=1;
	      //cerr << "launch: " << launch << endl;
	      cout << "<<<<<" << endl;
	      //cout << "launch: " << launch << endl;
	      doAlignment();
	      //cout << "done: " << launch << endl;
	      cout << ">>>>>" << endl;
	    }
	  }
	}
	myFile.close();
	/*alignproteins N;	
  
	N.setSequences( S[0], S[1] );
	N.loadInHouseMatrix("blosum62");
	N.setPenalty( 10,1);
	N.setmodus( true, false );
	N.doalignment();
  
	mAlignment W = N.getalignment();
  
  //cout << W.a1 << endl;
 // cout << W.mL << endl;
  //cout << W.a2 << endl;
	N.dumpAlignment( 50 ); */
	return EXIT_SUCCESS;
}
