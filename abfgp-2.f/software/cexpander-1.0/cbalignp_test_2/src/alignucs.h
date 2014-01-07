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

#ifndef __alignucsH
#define __alignucsH

#include "gmatrix.h"

typedef int TDNAmatrix[5][5];

class alignnucleotides: public generalmatrix
{
	public:
	
		alignnucleotides() {}
		~alignnucleotides() {}
	
		int nucToNr( char N );
		int score( int i, int j );
		char mType( int i, int j );
		void setScores( int sMatch, int sMissMatch, int ichar );
		
		
		void calcScores(); 
		void display_specific();

		void setUniformMode();
		void forceIll( int ill );
		
		void loadDNAMatrix( std::string Fname );
		void dumpDNAMatrix();
		
	protected:
	
	int Match;
	int MissMatch;
	int illchar;
	TDNAmatrix DNAmatrix;
	std::string matrixFname;
	bool uniformMode;
	
};
	
#endif 
