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

#include "alignucs.h"


int alignnucleotides::nucToNr( char N )
{
	int r = 4;
	if ( N == 'u' || N == 'U' || N == 'T' || N == 't' )
		r = 0;
	else
		if ( N == 'A' || N == 'a' )
			r = 1;
		else
			if ( N == 'C' || N == 'c' )
				r = 2;
			else
				if ( N == 'G' || N == 'g' )
					r = 3; 
	return r;
}
	
int alignnucleotides::score( int i, int j )
{
	// old version save
	/*
	if ( Sequence1.iSequence[ i - 1 ] == Sequence2.iSequence[ j - 1] )
		return Match;
	else
		return MissMatch;*/
		
	return DNAmatrix[ nucToNr( Sequence1.iSequence[i - 1] ) ][ nucToNr( Sequence2.iSequence[j-1] ) ]; // new version 
}

char alignnucleotides::mType( int i, int j )
{
	if ( Sequence1.iSequence[ i - 1] == Sequence2.iSequence[ j - 1] )
		return '|';
	else
		return '.';
}

void alignnucleotides::setScores( int sMatch, int sMissMatch, int ichar )
{
	Match = sMatch;
	MissMatch = sMissMatch;
	illchar = ichar;
	setUniformMode(); // set uniform model 
}

void alignnucleotides::calcScores()
{
	alignment.positives = cCount( alignment.mL, '|' );
	alignment.identities = alignment.positives;
	alignment.gaps = cCount( alignment.mL, ' ' );
}

void alignnucleotides::display_specific()
{
	if ( uniformMode )
	{
		cout << "Model: Uniform\n";
		cout << "match:  " << Match << endl;
		cout << "missmatch:   " << MissMatch << endl;
		cout << "illchar:     " << illchar << endl;
	}
	else
	dumpDNAMatrix();	
} 

void alignnucleotides::setUniformMode()
{
	for ( int x = 0; x != 5; x++ )
		for ( int y = 0; y != 5; y++ )
		{
			if ( x == 4 || y == 4 )
				DNAmatrix[x][y] = illchar;
			else
				if ( x == y )
					DNAmatrix[x][y] = Match; 
				else
					DNAmatrix[x][y] = MissMatch; 
		}
	uniformMode = true;
} 

void alignnucleotides::forceIll( int ill )
{
	for ( int k = 0; k != 5; k++ )
	{
		DNAmatrix[k][4] = ill;
		DNAmatrix[4][k] = ill;
	}
} 

void alignnucleotides::loadDNAMatrix( std::string Fname )
{
	char *filename = new char[ Fname.length() + 1 ];
	Fname.copy( filename, Fname.length() );
	filename[Fname.length()] = 0;
	
	ifstream myfile;
	myfile.open( filename );
	std::string Istring;
	int N1, N2, sc;
	vector <std::string> sp;
	
	if ( myfile.is_open() )
	{
		while (! myfile.eof() )
		{
			getline( myfile, Istring );
			if ( Istring.length() > 0 )
			{
				sp = split( Istring, '>', false );
				if ( sp.size() == 3 )
				{
					N1 = nucToNr( sp[0][0] );
					N2 = nucToNr( sp[1][0] );
					sc = atoi( sp[2].c_str() );
						
					if ( N1 == 4 || N2 == 4 )
						forceIll( sc );
					else
						DNAmatrix[N1][N2] = sc;
				}
			}
		}
		matrixFname = Fname;
		uniformMode = false;
		myfile.close();
	}
	delete [] filename; 
}

void alignnucleotides::dumpDNAMatrix()
{
	cout << "Matrix-file: " << matrixFname << endl; 
	for ( int y = 0; y != 5; y++ )
	{
		for ( int x = 0; x != 5; x++ )
		{
			if ( x > 0 )
				cout << "\t";
			cout << DNAmatrix[x][y];
		}
		cout << endl;
	}
} 
			