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

#ifndef __ebasicsH
#define __ebasicsH

#include <string>
#include <iostream>
#include <vector>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>

using namespace std;

vector< string > split( string& S, char c, bool keepEmpty );

vector< string > split_low( string& S, bool keepEmpty );

int scanLow( string& S, int pos, int step );

int scanFor( string& S, int pos, int step );

int scanNot( string& S, int pos, char c, int step ); 

string stripL_low( string S ); 

string stripR_low( string S );

string strip_low( string S); 

string stripL( string S, char c );

string stripR( string S, char c ); 

string strip( string S, char c ); 

string revString( string S );

int cCount(string& S, char ch);

int cCount_NNot(string& S, char ch);

string stringOf(int N, char ch);

string iSpace( int N );

int Idigs( int N );

struct generalSequence 
{
	string iHeader;
	string iSequence;
};

vector<generalSequence> readFastaFile(char* Fname );

struct clineItem{
string com;
string Istring;
bool switched;
int Ivalue; 
short T;
};

generalSequence nextFasta( ifstream& leStream );

typedef vector< clineItem >  clineItemList;

clineItem clI( string com, string Idefault, bool switched , int IvDefault, short T );

int gItemNr( clineItemList& L, string com );

clineItemList handleCommandLine( clineItemList Ilist );

clineItem clI( string com, string Idefault, bool switched ,int  IvDefault,short T );

int gItemNr( clineItemList& L, string com );

string handleCommandLine( clineItemList& Ilist, int argc, char* argv[] );

struct TPos{
	int x;
	int y;
};

TPos cPos( int x, int y );


#endif 


